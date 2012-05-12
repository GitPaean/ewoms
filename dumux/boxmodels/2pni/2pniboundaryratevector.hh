// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Implements a boundary vector for the fully implicit non-isothermal two-phase model.
 */
#ifndef DUMUX_BOX_2PNI_BOUNDARY_RATE_VECTOR_HH
#define DUMUX_BOX_2PNI_BOUNDARY_RATE_VECTOR_HH

#include "2pniproperties.hh"

#include <dumux/boxmodels/2p/2pboundaryratevector.hh>
#include <dumux/common/valgrind.hh>

#include <dune/common/fvector.hh>

namespace Dumux
{
/*!
 * \ingroup 2P2CNIModel
 *
 * \brief Implements a boundary vector for the fully implicit non-isothermal two-phase model.
 */
template <class TypeTag>
class TwoPNIBoundaryRateVector
    : public TwoPBoundaryRateVector<TypeTag>
{
    typedef TwoPBoundaryRateVector<TypeTag> ParentType;
    friend class TwoPBoundaryRateVector<TypeTag>;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { energyEqIdx = Indices::energyEqIdx };

public:
    /*!
     * \brief Default constructor
     */
    TwoPNIBoundaryRateVector()
        : ParentType()
    { };

    /*!
     * \brief Constructor with assignment from scalar
     */
    TwoPNIBoundaryRateVector(Scalar value)
        : ParentType(value)
    { };

    /*!
     * \brief Copy constructor
     */
    TwoPNIBoundaryRateVector(const TwoPNIBoundaryRateVector &value)
        : ParentType(value)
    { };


protected:
    template <class FluidState>
    void enthalpyFlux_(const FluxVariables &fluxVars,
                       const VolumeVariables &insideVolVars,
                       const FluidState &fs,
                       const typename FluidSystem::ParameterCache &paramCache,
                       int phaseIdx,
                       Scalar density)
    {
        Scalar enthalpy =
            (fs.pressure(phaseIdx) > insideVolVars.fluidState().pressure(phaseIdx))
            ? FluidSystem::enthalpy(fs, paramCache, phaseIdx)
            : insideVolVars.fluidState().enthalpy(phaseIdx);
        
        (*this)[energyEqIdx] += 
            fluxVars.filterVelocityNormal(phaseIdx)
            * enthalpy
            * density;
        Valgrind::CheckDefined((*this)[energyEqIdx]);
    }
};

} // end namepace

#endif