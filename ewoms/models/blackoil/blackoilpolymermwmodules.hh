// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief Contains the classes required to extend the black-oil model by polymer.
 */
#ifndef EWOMS_BLACK_OIL_POLYMER_MW_MODULE_HH
#define EWOMS_BLACK_OIL_POLYMER_MW_MODULE_HH

#include "blackoilproperties.hh"
#include <ewoms/io/vtkblackoilpolymermwmodule.hh>
#include <ewoms/models/common/quantitycallbacks.hh>

#include <opm/material/common/Tabulated1DFunction.hpp>

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyadsTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlymaxTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyrockTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyviscTable.hpp>
#endif

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>

#include <string>

namespace Ewoms {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by polymer, and related functionality related to polymer molecular
 *        weight transport and degradation.
 */
template <class TypeTag, bool enablePolymerMWV = GET_PROP_VALUE(TypeTag, EnablePolymerMW)>
class BlackOilPolymerMWModule
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef Opm::MathToolbox<Evaluation> Toolbox;

    typedef typename Opm::Tabulated1DFunction<Scalar> TabulatedFunction;

    static constexpr unsigned polymerConcentrationIdx = Indices::polymerConcentrationIdx;
    static constexpr unsigned contiPolymerEqIdx = Indices::contiPolymerEqIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;


    static constexpr unsigned enablePolymerMW = enablePolymerMWV;
    static constexpr unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr unsigned numPhases = FluidSystem::numPhases;

public:
    enum AdsorptionBehaviour { Desorption = 1, NoDesorption = 2 };

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize all internal data structures needed by the polymer module
     */
    static void initFromDeck(const Opm::Deck& deck, const Opm::EclipseState& eclState)
    {
        if (!deck.hasKeyword("POLYMW"))
            return;

        // TODO: some other sanity check will also be done here.

        // some sanity checks: if polymers are enabled, the POLYMER keyword must be
        // present, if polymers are disabled the keyword must not be present.
        if (enablePolymerMW && !deck.hasKeyword("POLYMER")) {
            throw std::runtime_error("Non-trivial polymer treatment requested at compile time, but "
                                     "the deck does not contain the POLYMER keyword");
        }
        else if (!enablePolymerMW && deck.hasKeyword("POLYMW")) {
            throw std::runtime_error("Polymer treatment disabled at compile time, but the deck "
                                     "contains the POLYMER keyword");
        }

        const auto& tableManager = eclState.getTableManager();

        unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
        setNumSatRegions(numSatRegions);

        // initialize the objects which deal with the PLYROCK keyword
        const auto& plyrockTables = tableManager.getPlyrockTables();
        if (!plyrockTables.empty()) {
            assert(numSatRegions == plyrockTables.size());
            for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++ satRegionIdx) {
                const auto& plyrockTable = plyrockTables.template getTable<Opm::PlyrockTable>(satRegionIdx);
                setPlyrock(satRegionIdx,
                           plyrockTable.getDeadPoreVolumeColumn()[satRegionIdx],
                           plyrockTable.getResidualResistanceFactorColumn()[satRegionIdx],
                           plyrockTable.getRockDensityFactorColumn()[satRegionIdx],
                           static_cast<AdsorptionBehaviour>(plyrockTable.getAdsorbtionIndexColumn()[satRegionIdx]),
                           plyrockTable.getMaxAdsorbtionColumn()[satRegionIdx]);
            }
        } else {
            throw std::runtime_error("PLYROCK must be specified in POLYMER runs\n");
        }

        // initialize the objects which deal with the PLYADS keyword
        const auto& plyadsTables = tableManager.getPlyadsTables();
        if (!plyadsTables.empty()) {
            assert(numSatRegions == plyadsTables.size());
            for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++ satRegionIdx) {
                const auto& plyadsTable = plyadsTables.template getTable<Opm::PlyadsTable>(satRegionIdx);
                // Copy data
                const auto& c = plyadsTable.getPolymerConcentrationColumn();
                const auto& ads = plyadsTable.getAdsorbedPolymerColumn();
                plyadsAdsorbedPolymer_[satRegionIdx].setXYContainers(c, ads);
            }
        } else {
            throw std::runtime_error("PLYADS must be specified in POLYMER runs\n");
        }


        unsigned numPvtRegions = tableManager.getTabdims().getNumPVTTables();
        setNumPvtRegions(numPvtRegions);

        // initialize the objects which deal with the PLYVISC keyword
        const auto& plyviscTables = tableManager.getPlyviscTables();
        if (!plyviscTables.empty()) {

            assert(numPvtRegions == plyviscTables.size());
            for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++ pvtRegionIdx) {
                const auto& plyadsTable = plyviscTables.template getTable<Opm::PlyviscTable>(pvtRegionIdx);
                // Copy data
                const auto& c = plyadsTable.getPolymerConcentrationColumn();
                const auto& visc = plyadsTable.getViscosityMultiplierColumn();
                plyviscViscosityMultiplierTable_[pvtRegionIdx].setXYContainers(c, visc);
            }
        } else {
            throw std::runtime_error("PLYVISC must be specified in POLYMER runs\n");
        }

        // initialize the objects which deal with the PLYMAX keyword
        const auto& plymaxTables = tableManager.getPlymaxTables();
        unsigned numMixRegions = plymaxTables.size();
        setNumMixRegions(numMixRegions);
        if (!plymaxTables.empty()) {
            for (unsigned mixRegionIdx = 0; mixRegionIdx < numMixRegions; ++ mixRegionIdx) {
                const auto& plymaxTable = plymaxTables.template getTable<Opm::PlymaxTable>(mixRegionIdx);
                setPlymax(mixRegionIdx, plymaxTable.getPolymerConcentrationColumn()[mixRegionIdx]);
            }
        } else {
            throw std::runtime_error("PLYMAX must be specified in POLYMER runs\n");
        }

        if (deck.hasKeyword("PLMIXPAR")) {
            // initialize the objects which deal with the PLMIXPAR keyword
            for (unsigned mixRegionIdx = 0; mixRegionIdx < numMixRegions; ++ mixRegionIdx) {
                const auto& plmixparRecord = deck.getKeyword("PLMIXPAR").getRecord(mixRegionIdx);
                setPlmixpar(mixRegionIdx, plmixparRecord.getItem("TODD_LONGSTAFF").getSIDouble(0));
            }
        } else {
            throw std::runtime_error("PLMIXPAR must be specified in POLYMER runs\n");
        }
    }
#endif

    /*!
     * \brief Specify the number of satuation regions.
     *
     * This must be called before setting the PLYROCK and PLYADS of any region.
     */
    static void setNumSatRegions(unsigned numRegions)
    {
        plyrockDeadPoreVolume_.resize(numRegions);
        plyrockResidualResistanceFactor_.resize(numRegions);
        plyrockRockDensityFactor_.resize(numRegions);
        plyrockAdsorbtionIndex_.resize(numRegions);
        plyrockMaxAdsorbtion_.resize(numRegions);
        plyadsAdsorbedPolymer_.resize(numRegions);
    }

    /*!
     * \brief Specify the polymer rock properties a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    static void setPlyrock(unsigned satRegionIdx,
                           const Scalar& plyrockDeadPoreVolume,
                           const Scalar& plyrockResidualResistanceFactor,
                           const Scalar& plyrockRockDensityFactor,
                           const Scalar& plyrockAdsorbtionIndex,
                           const Scalar& plyrockMaxAdsorbtion)
    {
        plyrockDeadPoreVolume_[satRegionIdx] = plyrockDeadPoreVolume;
        plyrockResidualResistanceFactor_[satRegionIdx] = plyrockResidualResistanceFactor;
        plyrockRockDensityFactor_[satRegionIdx] = plyrockRockDensityFactor;
        plyrockAdsorbtionIndex_[satRegionIdx] = plyrockAdsorbtionIndex;
        plyrockMaxAdsorbtion_[satRegionIdx] = plyrockMaxAdsorbtion;
    }

    /*!
     * \brief Specify the polymer rock properties a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    static void setAdsrock(unsigned satRegionIdx,
                           const TabulatedFunction& plyadsAdsorbedPolymer)
    {
        plyadsAdsorbedPolymer_[satRegionIdx] = plyadsAdsorbedPolymer;
    }

    /*!
     * \brief Specify the number of pvt regions.
     *
     * This must be called before setting the PLYVISC of any region.
     */
    static void setNumPvtRegions(unsigned numRegions)
    {
        plyviscViscosityMultiplierTable_.resize(numRegions);
    }

    /*!
     * \brief Specify the polymer viscosity a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    static void setPlyvisc(unsigned satRegionIdx,
                           const TabulatedFunction& plyviscViscosityMultiplierTable)
    {
        plyviscViscosityMultiplierTable_[satRegionIdx] = plyviscViscosityMultiplierTable;
    }

    /*!
     * \brief Specify the number of mix regions.
     *
     * This must be called before setting the PLYMAC and PLMIXPAR of any region.
     */
    static void setNumMixRegions(unsigned numRegions)
    {
        plymaxMaxConcentration_.resize(numRegions);
        plymixparToddLongstaff_.resize(numRegions);
    }

    /*!
     * \brief Specify the maximum polymer concentration a single region.
     *
     * The index of specified here must be in range [0, numMixRegionIdx)
     */
    static void setPlymax(unsigned mixRegionIdx,
                          const Scalar& plymaxMaxConcentration)
    {
        plymaxMaxConcentration_[mixRegionIdx] = plymaxMaxConcentration;
    }

    /*!
     * \brief Specify the maximum polymer concentration a single region.
     *
     * The index of specified here must be in range [0, numMixRegionIdx)
     */
    static void setPlmixpar(unsigned mixRegionIdx,
                            const Scalar& plymixparToddLongstaff)
    {
        plymixparToddLongstaff_[mixRegionIdx] = plymixparToddLongstaff;
    }



    /*!
     * \brief Register all run-time parameters for the black-oil polymer module.
     */
    static void registerParameters()
    {
        if (!enablePolymerMW)
            // polymers have been disabled at compile time
            return;

        // TODO: it needs a new VTK module
        // Ewoms::VtkBlackOilPolymerMWModule<TypeTag>::registerParameters();
        Ewoms::VtkBlackOilPolymerModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all polymer specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if (!enablePolymerMW)
            // polymers have been disabled at compile time
            return;

        // TODO: it needs a new VTK module
        // model.addOutputModule(new Ewoms::VtkBlackOilPolymerMWModule<TypeTag>(simulator));
        model.addOutputModule(new Ewoms::VtkBlackOilPolymerModule<TypeTag>(simulator));
    }

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if (!enablePolymerMW)
            // polymers have been disabled at compile time
            return false;

        return pvIdx == polymerConcentrationIdx;
    }

    static std::string primaryVarName(unsigned pvIdx OPM_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        return "polymer_waterconcentration";
    }

    static Scalar primaryVarWeight(unsigned pvIdx OPM_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enablePolymerMW)
            return false;

        return eqIdx == contiPolymerEqIdx;
    }

    static std::string eqName(unsigned eqIdx OPM_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        return "conti^polymer";
    }

    static Scalar eqWeight(unsigned eqIdx OPM_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    // must be called after water storage is computed
    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if (!enablePolymerMW)
            return;

        const auto& fs = intQuants.fluidState();

        LhsEval surfaceVolumeWater =
                Toolbox::template decay<LhsEval>(fs.saturation(waterPhaseIdx))
                * Toolbox::template decay<LhsEval>(fs.invB(waterPhaseIdx))
                * Toolbox::template decay<LhsEval>(intQuants.porosity());

        // avoid singular matrix if no water is present.
        surfaceVolumeWater = Opm::max(surfaceVolumeWater, 1e-10);

        // polymer in water phase
        storage[contiPolymerEqIdx] += surfaceVolumeWater
                * Toolbox::template decay<LhsEval>(intQuants.polymerConcentration())
                * (1.0 - Toolbox::template decay<LhsEval>(intQuants.polymerDeadPoreVolume()));

        // polymer in solid phase
        storage[contiPolymerEqIdx] +=
                Toolbox::template decay<LhsEval>(1.0 - intQuants.porosity())
                * Toolbox::template decay<LhsEval>(intQuants.polymerRockDensity())
                * Toolbox::template decay<LhsEval>(intQuants.polymerAdsorption());

    }

    static void computeFlux(RateVector& flux,
                            const ElementContext& elemCtx,
                            unsigned scvfIdx,
                            unsigned timeIdx)

    {
        if (!enablePolymerMW)
            return;

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        unsigned upIdx = extQuants.upstreamIndex(FluidSystem::waterPhaseIdx);
        unsigned inIdx = extQuants.interiorIndex();
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const unsigned contiWaterEqIdx = Indices::conti0EqIdx + Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);


        if (upIdx == inIdx) {
            flux[contiPolymerEqIdx] =
                    extQuants.volumeFlux(waterPhaseIdx)
                    *up.fluidState().invB(waterPhaseIdx)
                    *up.polymerViscosityCorrection()
                    /extQuants.polymerShearFactor()
                    *up.polymerConcentration();

            // modify water
            flux[contiWaterEqIdx] /=
                    extQuants.waterShearFactor();
        } else {
            flux[contiPolymerEqIdx] =
                    extQuants.volumeFlux(waterPhaseIdx)
                    *Opm::decay<Scalar>(up.fluidState().invB(waterPhaseIdx))
                    *Opm::decay<Scalar>(up.polymerViscosityCorrection())
                    /Opm::decay<Scalar>(extQuants.polymerShearFactor())
                    *Opm::decay<Scalar>(up.polymerConcentration());

            // modify water
            flux[contiWaterEqIdx] /=
                    Opm::decay<Scalar>(extQuants.waterShearFactor());
        }

    }

    /*!
     * \brief Assign the polymer specific primary variables to a PrimaryVariables object
     */
    static void assignPrimaryVars(PrimaryVariables& priVars,
                                  Scalar polymerConcentration)
    {
        if (!enablePolymerMW)
            return;

        priVars[polymerConcentrationIdx] = polymerConcentration;
    }

    /*!
     * \brief Do a Newton-Raphson update the primary variables of the polymers.
     */
    static void updatePrimaryVars(PrimaryVariables& newPv,
                                  const PrimaryVariables& oldPv,
                                  const EqVector& delta)
    {
        if (!enablePolymerMW)
            return;

        // do a plain unchopped Newton update
        newPv[polymerConcentrationIdx] = oldPv[polymerConcentrationIdx] - delta[polymerConcentrationIdx];
    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables& oldPv OPM_UNUSED,
                                     const EqVector& delta OPM_UNUSED)
    {
        // do not consider consider the cange of polymer primary variables for
        // convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    /*!
     * \brief Return how much a residual is considered an error
     */
    static Scalar computeResidualError(const EqVector& resid)
    {
        // do not weight the residual of polymers when it comes to convergence
        return std::abs(Toolbox::scalarValue(resid[contiPolymerEqIdx]));
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if (!enablePolymerMW)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);
        const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
        outstream << priVars[polymerConcentrationIdx];
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if (!enablePolymerMW)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);
        PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
        PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

        instream >> priVars0[polymerConcentrationIdx];

        // set the primary variables for the beginning of the current time step.
        priVars1 = priVars0[polymerConcentrationIdx];
    }

    static const Scalar plyrockDeadPoreVolume(const ElementContext& elemCtx,
                                              unsigned scvIdx,
                                              unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockDeadPoreVolume_[satnumRegionIdx];
    }

    static const Scalar plyrockResidualResistanceFactor(const ElementContext& elemCtx,
                                                        unsigned scvIdx,
                                                        unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockResidualResistanceFactor_[satnumRegionIdx];
    }

    static const Scalar plyrockRockDensityFactor(const ElementContext& elemCtx,
                                                 unsigned scvIdx,
                                                 unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockRockDensityFactor_[satnumRegionIdx];
    }

    static const Scalar plyrockAdsorbtionIndex(const ElementContext& elemCtx,
                                               unsigned scvIdx,
                                               unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockAdsorbtionIndex_[satnumRegionIdx];
    }

    static const Scalar plyrockMaxAdsorbtion(const ElementContext& elemCtx,
                                             unsigned scvIdx,
                                             unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockMaxAdsorbtion_[satnumRegionIdx];
    }

    static const TabulatedFunction& plyadsAdsorbedPolymer(const ElementContext& elemCtx,
                                                          unsigned scvIdx,
                                                          unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyadsAdsorbedPolymer_[satnumRegionIdx];
    }

    static const TabulatedFunction& plyviscViscosityMultiplierTable(const ElementContext& elemCtx,
                                                                    unsigned scvIdx,
                                                                    unsigned timeIdx)
    {
        unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyviscViscosityMultiplierTable_[pvtnumRegionIdx];
    }

    static const TabulatedFunction& plyviscViscosityMultiplierTable(unsigned pvtnumRegionIdx)
    {
        return plyviscViscosityMultiplierTable_[pvtnumRegionIdx];
    }

    static const Scalar plymaxMaxConcentration(const ElementContext& elemCtx,
                                               unsigned scvIdx,
                                               unsigned timeIdx)
    {
        unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plymaxMaxConcentration_[polymerMixRegionIdx];
    }

    static const Scalar plymixparToddLongstaff(const ElementContext& elemCtx,
                                               unsigned scvIdx,
                                               unsigned timeIdx)
    {
        unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plymixparToddLongstaff_[polymerMixRegionIdx];
    }

private:
    static std::vector<Scalar> plyrockDeadPoreVolume_;
    static std::vector<Scalar> plyrockResidualResistanceFactor_;
    static std::vector<Scalar> plyrockRockDensityFactor_;
    static std::vector<Scalar> plyrockAdsorbtionIndex_;
    static std::vector<Scalar> plyrockMaxAdsorbtion_;
    static std::vector<TabulatedFunction> plyadsAdsorbedPolymer_;
    static std::vector<TabulatedFunction> plyviscViscosityMultiplierTable_;
    static std::vector<Scalar> plymaxMaxConcentration_;
    static std::vector<Scalar> plymixparToddLongstaff_;
};



template <class TypeTag, bool enablePolymerMWV>
std::vector<typename BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::Scalar>
BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::plyrockDeadPoreVolume_;

template <class TypeTag, bool enablePolymerMWV>
std::vector<typename BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::Scalar>
BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::plyrockResidualResistanceFactor_;

template <class TypeTag, bool enablePolymerMWV>
std::vector<typename BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::Scalar>
BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::plyrockRockDensityFactor_;

template <class TypeTag, bool enablePolymerMWV>
std::vector<typename BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::Scalar>
BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::plyrockAdsorbtionIndex_;

template <class TypeTag, bool enablePolymerMWV>
std::vector<typename BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::Scalar>
BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::plyrockMaxAdsorbtion_;

template <class TypeTag, bool enablePolymerMWV>
std::vector<typename BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::TabulatedFunction>
BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::plyadsAdsorbedPolymer_;

template <class TypeTag, bool enablePolymerMWV>
std::vector<typename BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::TabulatedFunction>
BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::plyviscViscosityMultiplierTable_;

template <class TypeTag, bool enablePolymerMWV>
std::vector<typename BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::Scalar>
BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::plymaxMaxConcentration_;

template <class TypeTag, bool enablePolymerMWV>
std::vector<typename BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::Scalar>
BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::plymixparToddLongstaff_;

/*!
 * \ingroup BlackOil
 * \class Ewoms::BlackOilPolymerMWIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        polymers extension of the black-oil model.
 */
template <class TypeTag, bool enablePolymerMWV = GET_PROP_VALUE(TypeTag, EnablePolymerMW)>
class BlackOilPolymerMWIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef BlackOilPolymerMWModule<TypeTag> PolymerMWModule;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    static constexpr int polymerConcentrationIdx = Indices::polymerConcentrationIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;


public:

    /*!
     * \brief Update the intensive properties needed to handle polymers from the
     *        primary variables
     *
     */
    /* void polymerPropertiesUpdate_(const ElementContext& elemCtx,
                                  unsigned dofIdx,
                                  unsigned timeIdx)
    {
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        polymerConcentration_ = priVars.makeEvaluation(polymerConcentrationIdx, timeIdx);
        const Scalar cmax = PolymerMWModule::plymaxMaxConcentration(elemCtx, dofIdx, timeIdx);

        // permeability reduction due to polymer
        const Scalar& maxAdsorbtion = PolymerMWModule::plyrockMaxAdsorbtion(elemCtx, dofIdx, timeIdx);
        const auto& plyadsAdsorbedPolymer = PolymerMWModule::plyadsAdsorbedPolymer(elemCtx, dofIdx, timeIdx);
        polymerAdsorption_ = plyadsAdsorbedPolymer.eval(polymerConcentration_, *//*extrapolate=*//* true);
        if (PolymerMWModule::plyrockAdsorbtionIndex(elemCtx, dofIdx, timeIdx) == PolymerMWModule::NoDesorption ) {
            const Scalar& maxPolymerAdsorption = elemCtx.problem().maxPolymerAdsorption(elemCtx, dofIdx, timeIdx);
            polymerAdsorption_ = std::max(Evaluation(maxPolymerAdsorption) , polymerAdsorption_);
        }

        // compute resitanceFactor
        const Scalar& residualResistanceFactor = PolymerMWModule::plyrockResidualResistanceFactor(elemCtx, dofIdx, timeIdx);
        const Evaluation resistanceFactor = 1.0 + (residualResistanceFactor - 1.0) * polymerAdsorption_ / maxAdsorbtion;

        // compute effective viscosities
        auto& fs = asImp_().fluidState_;
        const Evaluation& muWater = fs.viscosity(waterPhaseIdx);
        const auto& viscosityMultiplier = PolymerMWModule::plyviscViscosityMultiplierTable(elemCtx, dofIdx, timeIdx);
        Evaluation viscosityMixture = viscosityMultiplier.eval(polymerConcentration_, */ /*extrapolate=*/ /* true) * muWater;

        // Do the Todd-Longstaff mixing
        const Scalar plymixparToddLongstaff = PolymerMWModule::plymixparToddLongstaff(elemCtx, dofIdx, timeIdx);
        Evaluation viscosityPolymer = viscosityMultiplier.eval(cmax, */ /*extrapolate=*/ /* true) * muWater;
        Evaluation viscosityPolymerEffective = pow(viscosityMixture, plymixparToddLongstaff) * pow(viscosityPolymer, 1.0 - plymixparToddLongstaff);
        Evaluation viscosityWaterEffective = pow(viscosityMixture, plymixparToddLongstaff) * pow(muWater, 1.0 - plymixparToddLongstaff);

        Evaluation cbar = polymerConcentration_ / cmax;
        // waterViscosity / effectiveWaterViscosity
        waterViscosityCorrection_ = muWater * ( (1.0 - cbar) / viscosityWaterEffective + cbar / viscosityPolymerEffective );

        // adjust water mobility
        asImp_().mobility_[waterPhaseIdx] *= waterViscosityCorrection_ / resistanceFactor;

        // effectiveWaterViscosity / effectivePolymerViscosity
        polymerViscosityCorrection_ =  (muWater / waterViscosityCorrection_) / viscosityPolymerEffective;

        // update rock properties
        polymerDeadPoreVolume_ = PolymerMWModule::plyrockDeadPoreVolume(elemCtx, dofIdx, timeIdx);
        polymerRockDensity_ = PolymerMWModule::plyrockRockDensityFactor(elemCtx, dofIdx, timeIdx);
    } */

    /* const Evaluation& polymerConcentration() const
    { return polymerConcentration_; }

    const Scalar& polymerDeadPoreVolume() const
    { return polymerDeadPoreVolume_; }

    const Evaluation& polymerAdsorption() const
    { return polymerAdsorption_; }

    const Scalar& polymerRockDensity() const
    { return polymerRockDensity_; }

    // effectiveWaterViscosity / effectivePolymerViscosity
    const Evaluation& polymerViscosityCorrection() const
    { return polymerViscosityCorrection_; }

    // waterViscosity / effectiveWaterViscosity
    const Evaluation& waterViscosityCorrection() const
    { return waterViscosityCorrection_; } */


protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation polymerConcentration_;
    Scalar polymerDeadPoreVolume_;
    Scalar polymerRockDensity_;
    Evaluation polymerAdsorption_;
    Evaluation polymerViscosityCorrection_;
    Evaluation waterViscosityCorrection_;


};

// TODO: not sure if we need this
/* template <class TypeTag>
class BlackOilPolymerMWIntensiveQuantities<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    void polymerPropertiesUpdate_(const ElementContext& elemCtx OPM_UNUSED,
                                  unsigned scvIdx OPM_UNUSED,
                                  unsigned timeIdx OPM_UNUSED)
    { }

    const Evaluation& polymerConcentration() const
    { throw std::runtime_error("polymerConcentration() called but polymers are disabled"); }

    const Evaluation& polymerDeadPoreVolume() const
    { throw std::runtime_error("polymerDeadPoreVolume() called but polymers are disabled"); }

    const Evaluation& polymerAdsorption() const
    { throw std::runtime_error("polymerAdsorption() called but polymers are disabled"); }

    const Evaluation& polymerRockDensity() const
    { throw std::runtime_error("polymerRockDensity() called but polymers are disabled"); }

    const Evaluation& polymerViscosityCorrection() const
    { throw std::runtime_error("polymerViscosityCorrection() called but polymers are disabled"); }

    const Evaluation& waterViscosityCorrection() const
    { throw std::runtime_error("waterViscosityCorrection() called but polymers are disabled"); }
}; */


/*!
 * \ingroup BlackOil
 * \class Ewoms::BlackOilPolymerExtensiveQuantities
 *
 * \brief Provides the polymer specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enablePolymerMWV = GET_PROP_VALUE(TypeTag, EnablePolymerMW)>
class BlackOilPolymerMWExtensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef Opm::MathToolbox<Evaluation> Toolbox;

    typedef BlackOilPolymerMWModule<TypeTag> PolymerMWModule;


    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr unsigned waterPhaseIdx =  FluidSystem::waterPhaseIdx;

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldVector<Evaluation, dimWorld> DimEvalVector;

public:
    /*!
     * \brief Method which calculates the shear factor based on flow velocity
     *
     * This is the variant of the method which assumes that the problem is specified
     * using transmissibilities, i.e., *not* via permeabilities.
     */
private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }
};

/* template <class TypeTag>
class BlackOilPolymerMWExtensiveQuantities<TypeTag, false>
{
public:
}; */


} // namespace Ewoms

#endif
