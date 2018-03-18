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
    // TODO: make here we should define a typedef for BlackOilPolymerMWIntensiveQuantities for later use
    // TODO: and also BlackOilPolymerMWExtensiveQuantities for later use
    // TODO: probably also the similar things in BlackOilPolymerModule
    // TODO: unless I want to change all the name of the functions to be different from the old polymer implementation
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
    // TODO: there will be an index for polymer molecular weight
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;

    // TODO: why not bool
    static constexpr unsigned enablePolymerMW = enablePolymerMWV;
    static constexpr unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr unsigned numPhases = FluidSystem::numPhases;

public:
    enum AdsorptionBehaviour { Desorption = 1, NoDesorption = 2 };

    // the struct contains
    struct PlyvmhCoefficients {
        Scalar k_mh;
        Scalar a_mh;
        Scalar gamma;
        Scalar kappa;
    };


#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize all internal data structures needed by the polymer module
     */
    static void initFromDeck(const Opm::Deck& deck, const Opm::EclipseState& eclState)
    {
        // To use this module, you need to have keyword POLYMW for the polymer molecular weight
        // As a rule, keyword POLYMER must be specified together.
        // Since the viscosity is calculated based on the Mark-Houwink equation and Huggins equation, we need to
        // have keyword PLYVMH for the coefficients of the formulation.
        // For the regions, PLYVMH shares the same region with the mixing region specified with PLMIXNUM.

        // TODO: as the first part of the implementation, we use constant polymer weight, without tracking the
        // TODO: transport of the polymer molecular weight.

        // TODO: for some keywords work with the oridinary POLYMER keyword while not work with polymer molecular weight,
        // TODO: we should give a warning to tell these are not used and will take no effects.

        // TODO: the order of the checking needs some design here

        // sanity check
        if (enablePolymerMW && !deck.hasKeyword("POLYMW")) {
            throw std::runtime_error("Polymer molecular weight related functionality is requested at compile time, but "
                                     "the deck does not contain the POLYMW keyword");
        } else if (!enablePolymerMW && deck.hasKeyword("POLYMW")) {
            throw std::runtime_error("Polymer molecular weight related functionality is disabled at compile time, but "
                                     " the deck contains the POLYMW keyword");
        }

        if ( deck.hasKeyword("POLYMW") && !deck.hasKeyword("POLYMER") ) {
            throw std::logic_error("POLYMER keyword must be specified along with POLYMW keyword to activate polymer "
                                   "related functionality ");
        }

        const Opm::UnitSystem& unitSystem = deck.getActiveUnitSystem();

        if (unitSystem.getName() != "Metric") {
            throw std::logic_error("Only METRIC unit system is supported for the polymer molecular weight "
                                   "related simulation for now");
        }

        if (!deck.hasKeyword("POLYMW"))
            return;

        const auto& tableManager = eclState.getTableManager();

        unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
        setNumSatRegions(numSatRegions);

        // initialize the objects which deal with the PLYROCK keyword
        // TODO: the adsorb related typo needs to be fixed from the parser side
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
            throw std::runtime_error("PLYROCK must be specified in POLYMW runs\n");
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
                plyadsAdsorpedPolymer_[satRegionIdx].setXYContainers(c, ads);
            }
        } else {
            throw std::runtime_error("PLYADS must be specified in POLYMW runs\n");
        }

        // for the several keywords that are used for standard polymer simulation, we should not use them.
        // For the moment, we throw to let the user to remove them
        if (deck.hasKeyword("PLYVISC") || deck.hasKeyword("PLMIXPAR")) {
            Opm::OpmLog::warning("PLYVISC and PLMIXPAR should not be used in POLYMW runs, "
                    "they will not take effects. Different viscosity model based on PLYVMH is used \n");
        }

        if (deck.hasKeyword("PLYSHLOG") || deck.hasKeyword("PLYSHEAR") || deck.hasKeyword("SHRATE")) {
            Opm::OpmLog::warning("Shear calculation based on PLYSHLOG, PLYSHEAR or SHRATE is not supported "
                                 "with POLYMW runs, they will not take any effect. \n");
        }

        // initialize the objects which deal with the PLYMAX keyword
        const auto& plymaxTables = tableManager.getPlymaxTables();
        const unsigned numMixRegions = plymaxTables.size();
        setNumMixRegions(numMixRegions);
        if (!plymaxTables.empty()) {
            for (unsigned mixRegionIdx = 0; mixRegionIdx < numMixRegions; ++ mixRegionIdx) {
                const auto& plymaxTable = plymaxTables.template getTable<Opm::PlymaxTable>(mixRegionIdx);
                setPlymax(mixRegionIdx, plymaxTable.getPolymerConcentrationColumn()[mixRegionIdx]);
            }
        } else {
            throw std::runtime_error("PLYMAX must be specified in POLYMW runs\n");
        }

        // TODO: for now, we use the mixing region first for keyword PLYVMH
        // initialize the PLYVMH related
        const Opm::DeckKeyword& plyvmhKeyword = deck.getKeyword("PLYVMH");
        assert(plyvmhKeyword.size() == numMixRegions);
        if (plyvmhKeyword.size() > 0) {
            for (size_t region_idx = 0; region_idx < plyvmhKeyword.size(); ++region_idx) {
                const Opm::DeckRecord& record = plyvmhKeyword.getRecord(region_idx);
                plyvmhCoefficients_[region_idx].k_mh = record.getItem("K_MH").getSIDouble(0);
                plyvmhCoefficients_[region_idx].a_mh = record.getItem("A_MH").getSIDouble(0);
                plyvmhCoefficients_[region_idx].gamma = record.getItem("GAMMA").getSIDouble(0);
                plyvmhCoefficients_[region_idx].kappa = record.getItem("KAPPA").getSIDouble(0);
            }
        } else {
            throw std::runtime_error("PLYVMH keyword must be specified in POLYMW rus \n");
        }
    }
#endif

    /*!
     * \brief Specify the number of saturation regions.
     *
     * This must be called before setting the PLYROCK and PLYADS of any region.
     */
    static void setNumSatRegions(unsigned numRegions)
    {
        plyrockDeadPoreVolume_.resize(numRegions);
        plyrockResidualResistanceFactor_.resize(numRegions);
        plyrockRockDensityFactor_.resize(numRegions);
        plyrockAdsorptionIndex_.resize(numRegions);
        plyrockMaxAdsorption_.resize(numRegions);
        plyadsAdsorpedPolymer_.resize(numRegions);
    }

    /*!
     * \brief Specify the polymer rock properties a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    // TODO: the name of the function needs to be changed, since it is related to region
    // TODO: the function name should show it
    static void setPlyrock(unsigned satRegionIdx,
                           const Scalar& plyrockDeadPoreVolume,
                           const Scalar& plyrockResidualResistanceFactor,
                           const Scalar& plyrockRockDensityFactor,
                           const Scalar& plyrockAdsorptionIndex,
                           const Scalar& plyrockMaxAdsorption)
    {
        plyrockDeadPoreVolume_[satRegionIdx] = plyrockDeadPoreVolume;
        plyrockResidualResistanceFactor_[satRegionIdx] = plyrockResidualResistanceFactor;
        plyrockRockDensityFactor_[satRegionIdx] = plyrockRockDensityFactor;
        plyrockAdsorptionIndex_[satRegionIdx] = plyrockAdsorptionIndex;
        plyrockMaxAdsorption_[satRegionIdx] = plyrockMaxAdsorption;
    }

    /*!
     * \brief Specify the number of mix regions.
     *
     * This must be called before setting the PLYMAC and PLMIXPAR of any region.
     */
    static void setNumMixRegions(unsigned numRegions)
    {
        plymaxMaxConcentration_.resize(numRegions);
        plyvmhCoefficients_.resize(numRegions);
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
     * \brief Register all run-time parameters for the black-oil polymer module.
     */
    static void registerParameters()
    {
        if (!enablePolymerMW)
            // polymers have been disabled at compile time
            return;

        // TODO: it needs a new VTK module
        Ewoms::VtkBlackOilPolymerMWModule<TypeTag>::registerParameters();
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
        model.addOutputModule(new Ewoms::VtkBlackOilPolymerMWModule<TypeTag>(simulator));
    }

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if (!enablePolymerMW)
            // polymers have been disabled at compile time
            return false;

        // TODO: we need to check two possible indices
        return pvIdx == polymerConcentrationIdx;
    }

    static std::string primaryVarName(unsigned pvIdx OPM_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: basically, it will two possiblities
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

        // TODO: two indices to check here
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

        // TODO: the intensive quantities are defined based on multiple inheritance
        // TODO: either we need to use different naming for the functions or we need to find a way
        // TODO: to specify where the function it is from
        // TODO: how about intQuants.BlackOilPolymerMWIntensiveQuantities<TypeTag>::polymerAdsorption


        const auto& fs = intQuants.fluidState();

        LhsEval surfaceVolumeWater =
                Toolbox::template decay<LhsEval>(fs.saturation(waterPhaseIdx))
                * Toolbox::template decay<LhsEval>(fs.invB(waterPhaseIdx))
                * Toolbox::template decay<LhsEval>(intQuants.porosity());

        // avoid singular matrix if no water is present.
        surfaceVolumeWater = Opm::max(surfaceVolumeWater, 1e-10);

        // polymer in water phase
        // TODO: I have to specify which polymerConcentration I am using, polymer module or polymerMW module
        storage[contiPolymerEqIdx] += surfaceVolumeWater
                * Toolbox::template decay<LhsEval>(intQuants.polymerConcentrationMW())
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
                    * up.fluidState().invB(waterPhaseIdx)
                    * up.polymerConcentrationMW()
                    / up.waterPolymerViscosityCorrectionMW();
        } else {
            flux[contiPolymerEqIdx] =
                    extQuants.volumeFlux(waterPhaseIdx)
                    * Opm::decay<Scalar>(up.fluidState().invB(waterPhaseIdx))
                    * Opm::decay<Scalar>(up.polymerViscosityCorrection())
                    * Opm::decay<Scalar>(up.polymerConcentrationMW())
                    / Opm::decay<Scalar>(up.waterPolymerViscosityCorrectionMW());
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
        // TODO: there will be two primary variables needed to be updated
    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables& oldPv OPM_UNUSED,
                                     const EqVector& delta OPM_UNUSED)
    {
        // do not consider consider the change of polymer primary variables for
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
        // TODO: it will be two equations and two residuals
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
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockDeadPoreVolume_[satnumRegionIdx];
    }

    static const Scalar plyrockResidualResistanceFactor(const ElementContext& elemCtx,
                                                        unsigned scvIdx,
                                                        unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockResidualResistanceFactor_[satnumRegionIdx];
    }

    static const Scalar plyrockRockDensityFactor(const ElementContext& elemCtx,
                                                 unsigned scvIdx,
                                                 unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockRockDensityFactor_[satnumRegionIdx];
    }

    static const Scalar plyrockAdsorptionIndex(const ElementContext& elemCtx,
                                               unsigned scvIdx,
                                               unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockAdsorptionIndex_[satnumRegionIdx];
    }

    static const Scalar plyrockMaxAdsorption(const ElementContext& elemCtx,
                                             unsigned scvIdx,
                                             unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockMaxAdsorption_[satnumRegionIdx];
    }

    static const TabulatedFunction& plyadsAdsorpedPolymer(const ElementContext& elemCtx,
                                                          unsigned scvIdx,
                                                          unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyadsAdsorpedPolymer_[satnumRegionIdx];
    }

    static const Scalar plymaxMaxConcentration(const ElementContext& elemCtx,
                                               unsigned scvIdx,
                                               unsigned timeIdx)
    // TODO: mixing regions are used for polymer maximum concentration
    {
        unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plymaxMaxConcentration_[polymerMixRegionIdx];
    }


    static const PlyvmhCoefficients& plyvmhCoefficients(const ElementContext& elemCtx,
                                                        const unsigned scvIdx,
                                                        const unsigned timeIdx)
    {
        const unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyvmhCoefficients_[polymerMixRegionIdx];
    }

private:
    static std::vector<Scalar> plyrockDeadPoreVolume_;
    static std::vector<Scalar> plyrockResidualResistanceFactor_;
    static std::vector<Scalar> plyrockRockDensityFactor_;
    static std::vector<Scalar> plyrockAdsorptionIndex_;
    static std::vector<Scalar> plyrockMaxAdsorption_;
    static std::vector<TabulatedFunction> plyadsAdsorpedPolymer_;
    static std::vector<Scalar> plymaxMaxConcentration_;

    static std::vector<PlyvmhCoefficients> plyvmhCoefficients_;
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
BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::plyrockAdsorptionIndex_;

template <class TypeTag, bool enablePolymerMWV>
std::vector<typename BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::Scalar>
BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::plyrockMaxAdsorption_;

template <class TypeTag, bool enablePolymerMWV>
std::vector<typename BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::TabulatedFunction>
BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::plyadsAdsorpedPolymer_;

template <class TypeTag, bool enablePolymerMWV>
std::vector<typename BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::Scalar>
BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::plymaxMaxConcentration_;

template <class TypeTag, bool enablePolymerMWV>
std::vector<typename BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::PlyvmhCoefficients>
BlackOilPolymerMWModule<TypeTag, enablePolymerMWV>::plyvmhCoefficients_;

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
    // TODO: there is where the viscosity calculation will happen
    void polymerPropertiesUpdateMW_(const ElementContext& elemCtx,
                                  unsigned dofIdx,
                                  unsigned timeIdx)
    {
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        polymerConcentration_ = priVars.makeEvaluation(polymerConcentrationIdx, timeIdx);
        const Scalar cmax = PolymerMWModule::plymaxMaxConcentration(elemCtx, dofIdx, timeIdx);

        // permeability reduction due to polymer
        const Scalar& maxAdsorption = PolymerMWModule::plyrockMaxAdsorption(elemCtx, dofIdx, timeIdx);
        const auto& plyadsAdsorpedPolymer = PolymerMWModule::plyadsAdsorpedPolymer(elemCtx, dofIdx, timeIdx);
        polymerAdsorption_ = plyadsAdsorpedPolymer.eval(polymerConcentration_, /*extrapolate=*/ true);
        if (PolymerMWModule::plyrockAdsorptionIndex(elemCtx, dofIdx, timeIdx) == PolymerMWModule::NoDesorption ) {
            const Scalar& maxPolymerAdsorption = elemCtx.problem().maxPolymerAdsorption(elemCtx, dofIdx, timeIdx);
            polymerAdsorption_ = std::max(Evaluation(maxPolymerAdsorption) , polymerAdsorption_);
        }

        // compute resitanceFactor
        const Scalar& residualResistanceFactor = PolymerMWModule::plyrockResidualResistanceFactor(elemCtx, dofIdx, timeIdx);
        const Evaluation resistanceFactor = 1.0 + (residualResistanceFactor - 1.0) * polymerAdsorption_ / maxAdsorption;

        // update rock properties
        polymerDeadPoreVolume_ = PolymerMWModule::plyrockDeadPoreVolume(elemCtx, dofIdx, timeIdx);
        polymerRockDensity_ = PolymerMWModule::plyrockRockDensityFactor(elemCtx, dofIdx, timeIdx);

        // compute effective viscosities
        // TODO: here to compute the viscosity to update the waterViscosityCorrection_ and polyerViscosityCorrection_
        // TODO: the following part should be put into a function
        const auto& plyvmhCoefficients = PolymerMWModule::plyvmhCoefficients(elemCtx, dofIdx, timeIdx);
        const Scalar k_mh = plyvmhCoefficients.k_mh;
        const Scalar a_mh = plyvmhCoefficients.a_mh;
        const Scalar gamma = plyvmhCoefficients.gamma;
        const Scalar kappa = plyvmhCoefficients.kappa;

        // TODO: later, the polymer molecular weight should be changed to Evaluation
        const Scalar polymerMolecularWeight = 40.;
        // TODO: the following pow function should have a version for Evaluation type
        const Scalar intrinsicViscosity = k_mh * std::pow(polymerMolecularWeight * 1000., a_mh) * 1000.;
        // the meaning of this variable is not very clear, while it is x in the formulation employed
        const Evaluation x = 1.e-6 * polymerConcentration_ * intrinsicViscosity;
        // TODO: this is for multiplication, while the original one looks like for division correction,
        // TODO: it is possible they are used for mobility, then it make sense to use division correction
        // TODO: then we can just invert it
        waterPolymerViscosityCorrection_ = 1.0 + gamma * (x + kappa * x * x);

        // adjust water mobility
        asImp_().mobility_[waterPhaseIdx] *= waterPolymerViscosityCorrection_ / resistanceFactor;
    }

    const Evaluation& polymerConcentrationMW() const
    { return polymerConcentration_; }

    const Scalar& polymerDeadPoreVolumeMW() const
    { return polymerDeadPoreVolume_; }

    const Evaluation& polymerAdsorptionMW() const
    { return polymerAdsorption_; }

    const Scalar& polymerRockDensityMW() const
    { return polymerRockDensity_; }

    // viscosity correction for the water polymer mixture
    // since in this model, we are not considering the mixing parameter, so the polymer-water mixture will
    // share one viscosity
    const Evaluation& waterPolymerViscosityCorrectionMW() const
    { return waterPolymerViscosityCorrection_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation polymerConcentration_;
    Scalar polymerDeadPoreVolume_;
    Scalar polymerRockDensity_;
    Evaluation polymerAdsorption_;
    // TODO: with the new implementation, I think the following two will be one
    // TODO: basically, it means the mixing is always full, mixing parameter is 1.0 here,
    // TODO: although we do not use mixing parameter in the input and simulation
    Evaluation waterPolymerViscosityCorrection_;
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

    const Evaluation& polymerConcentrationMW() const
    { throw std::runtime_error("polymerConcentrationMW() called but polymers are disabled"); }

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
/* template <class TypeTag, bool enablePolymerMWV = GET_PROP_VALUE(TypeTag, EnablePolymerMW)>
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
private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }
}; */

/* template <class TypeTag>
class BlackOilPolymerMWExtensiveQuantities<TypeTag, false>
{
public:
}; */


} // namespace Ewoms

#endif
