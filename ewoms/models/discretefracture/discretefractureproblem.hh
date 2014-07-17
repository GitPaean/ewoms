/*
  Copyright (C) 2014 by Andreas Lauser

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
*/
/*!
 * \file
 *
 * \copydoc Ewoms::DiscreteFractureProblem
 */
#ifndef EWOMS_DISCRETE_FRACTURE_PROBLEM_HH
#define EWOMS_DISCRETE_FRACTURE_PROBLEM_HH

#include "discretefractureproperties.hh"

#include <ewoms/models/common/multiphasebaseproblem.hh>

#include <opm/core/utility/Average.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(HeatConductionLawParams);
NEW_PROP_TAG(EnableGravity);
NEW_PROP_TAG(VelocityModule);
}}

namespace Ewoms {
/*!
 * \ingroup Discretization
 *
 * \brief The base class for the problems of ECFV discretizations which deal
 *        with a multi-phase flow through a porous medium.
 */
template<class TypeTag>
class DiscreteFractureProblem
    : public MultiPhaseBaseProblem<TypeTag>
{
    typedef Ewoms::MultiPhaseBaseProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

    enum { dimWorld = GridView::dimensionworld };
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Problem::FvBaseProblem(Simulator &)
     */
    DiscreteFractureProblem(Simulator &simulator)
        : ParentType(simulator)
    {}

    /*!
     * \brief Returns the intrinsic permeability of a face due to a fracture.
     *
     * This method is specific to the finite volume discretizations. If left unspecified,
     * it calls the intrinsicPermeability() methods for the face's interior and exterior
     * finite volume and averages them harmonically. Note that if this function is
     * defined, the intrinsicPermeability() method does not need to be defined by the
     * problem (if a finite-volume discretization is used).
     */
    template <class Context>
    void fractureFaceIntrinsicPermeability(DimMatrix &result,
                                           const Context &context,
                                           int localFaceIdx, int timeIdx) const
    {
        const auto &scvf = context.stencil(timeIdx).interiorFace(localFaceIdx);
        int interiorElemIdx = scvf.interiorIndex();
        int exteriorElemIdx = scvf.exteriorIndex();
        const DimMatrix &K1 = asImp_().fractureIntrinsicPermeability(context, interiorElemIdx, timeIdx);
        const DimMatrix &K2 = asImp_().fractureIntrinsicPermeability(context, exteriorElemIdx, timeIdx);

        // entry-wise harmonic mean. this is almost certainly wrong if
        // you have off-main diagonal entries in your permeabilities!
        for (int i = 0; i < dimWorld; ++i)
            for (int j = 0; j < dimWorld; ++j)
                result[i][j] = Opm::utils::harmonicAverage(K1[i][j], K2[i][j]);
    }
    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$ at a given position due to a fracture
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    const DimMatrix &fractureIntrinsicPermeability(const Context &context,
                                                   int spaceIdx, int timeIdx) const
    {
        OPM_THROW(std::logic_error,
                   "Not implemented: Problem::fractureIntrinsicPermeability()");
    }

    /*!
     * \brief Returns the porosity [] inside fractures for a given control volume.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    Scalar fracturePorosity(const Context &context,
                            int spaceIdx, int timeIdx) const
    {
        OPM_THROW(std::logic_error,
                   "Not implemented: Problem::fracturePorosity()");
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }
    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // namespace Ewoms

#endif