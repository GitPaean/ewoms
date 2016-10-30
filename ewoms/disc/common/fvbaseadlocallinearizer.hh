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
 * \copydoc Ewoms::FvBaseAdLocalLinearizer
 */
#ifndef EWOMS_FV_BASE_AD_LOCAL_LINEARIZER_HH
#define EWOMS_FV_BASE_AD_LOCAL_LINEARIZER_HH

#include <opm/material/densead/Math.hpp>

#include "fvbaseproperties.hh"

#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Ewoms {
// forward declaration
template<class TypeTag>
class FvBaseAdLocalLinearizer;

namespace Properties {
// declare the property tags required for the finite differences local linearizer
NEW_TYPE_TAG(AutoDiffLocalLinearizer);

NEW_PROP_TAG(LocalLinearizer);
NEW_PROP_TAG(Evaluation);

NEW_PROP_TAG(LocalResidual);
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(Model);
NEW_PROP_TAG(PrimaryVariables);
NEW_PROP_TAG(ElementContext);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Evaluation);
NEW_PROP_TAG(GridView);

// set the properties to be spliced in
SET_TYPE_PROP(AutoDiffLocalLinearizer, LocalLinearizer,
              Ewoms::FvBaseAdLocalLinearizer<TypeTag>);

//! Set the function evaluation w.r.t. the primary variables
SET_PROP(AutoDiffLocalLinearizer, Evaluation)
{
private:
    static const int numEq = GET_PROP_VALUE(TypeTag, NumEq);

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::DenseAd::Evaluation<Scalar, numEq> type;
};

} // namespace Properties

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief Calculates the local residual and its Jacobian for a single element of the grid.
 *
 * This class uses automatic differentiation to calculate the partial derivatives (the
 * alternative is finite differences).
 */
template<class TypeTag>
class FvBaseAdLocalLinearizer
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, LocalLinearizer) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) LocalResidual;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef Dune::FieldVector<Scalar, numEq> ScalarVectorBlock;
    typedef Dune::FieldMatrix<Scalar, numEq, numEq> ScalarMatrixBlock;

    typedef Dune::BlockVector<ScalarVectorBlock> ScalarLocalBlockVector;
    typedef Dune::Matrix<ScalarMatrixBlock> ScalarLocalBlockMatrix;

#if __GNUC__ == 4 && __GNUC_MINOR__ <= 6
public:
    // make older GCCs happy by providing a public copy constructor (this is necessary
    // for their implementation of std::vector, although the method is never called...)
    FvBaseAdLocalLinearizer(const FvBaseAdLocalLinearizer&)
        : internalElemContext_(0)
    {}

#else
    // copying local residual objects around is a very bad idea, so we explicitly prevent
    // it...
    FvBaseAdLocalLinearizer(const FvBaseAdLocalLinearizer&) = delete;
#endif
public:
    FvBaseAdLocalLinearizer()
        : internalElemContext_(0)
    { }

    ~FvBaseAdLocalLinearizer()
    { delete internalElemContext_; }

    /*!
     * \brief Register all run-time parameters for the local jacobian.
     */
    static void registerParameters()
    { }

    /*!
     * \brief Initialize the local Jacobian object.
     *
     * At this point we can assume that everything has been allocated,
     * although some objects may not yet be completely initialized.
     *
     * \param simulator The simulator object of the simulation.
     */
    void init(Simulator &simulator)
    {
        simulatorPtr_ = &simulator;
        delete internalElemContext_;
        internalElemContext_ = new ElementContext(simulator);
    }

    /*!
     * \brief Compute an element's local Jacobian matrix and evaluate its residual.
     *
     * The local Jacobian for a given context is defined as the derivatives of the
     * residuals of all degrees of freedom featured by the stencil with regard to the
     * primary variables of the stencil's "primary" degrees of freedom. Adding the local
     * Jacobians for all elements in the grid will give the global Jacobian 'grad f(x)'.
     *
     * \param element The grid element for which the local residual and its local
     *                Jacobian should be calculated.
     */
    void linearize(const Element &element)
    {
        internalElemContext_->updateAll(element);

        linearize(*internalElemContext_);
    }

    /*!
     * \brief Compute an element's local Jacobian matrix and evaluate its residual.
     *
     * The local Jacobian for a given context is defined as the derivatives of the
     * residuals of all degrees of freedom featured by the stencil with regard to the
     * primary variables of the stencil's "primary" degrees of freedom. Adding the local
     * Jacobians for all elements in the grid will give the global Jacobian 'grad f(x)'.
     *
     * After calling this method the ElementContext is in an undefined state, so do not
     * use it anymore!
     *
     * \param elemCtx The element execution context for which the local residual and its
     *                local Jacobian should be calculated.
     */
    void linearize(ElementContext &elemCtx)
    {
        // update the weights of the primary variables for the context
        model_().updatePVWeights(elemCtx);

        resize_(elemCtx);
        reset_(elemCtx);

        // calculate the local residual
        localResidual_.eval(elemCtx);

        // calculate the local jacobian matrix
        int numPrimaryDof = elemCtx.numPrimaryDof(/*timeIdx=*/0);
        for (int dofIdx = 0; dofIdx < numPrimaryDof; dofIdx++) {
            // convert the local Jacobian matrix and the right hand side from the data
            // structures used by the automatic differentiation code to the conventional
            // ones used by the linear solver.
            updateLocalLinearization_(elemCtx, dofIdx);
        }

    }

    /*!
     * \brief Return reference to the local residual.
     */
    LocalResidual &localResidual()
    { return localResidual_; }

    /*!
     * \brief Return reference to the local residual.
     */
    const LocalResidual &localResidual() const
    { return localResidual_; }

    /*!
     * \brief Returns the local Jacobian matrix of the residual of a sub-control volume.
     *
     * \param domainScvIdx The local index of the sub control volume to which the primary
     *                     variables are associated with
     * \param rangeScvIdx The local index of the sub control volume which contains the
     *                    local residual
     */
    const ScalarMatrixBlock &jacobian(int domainScvIdx, int rangeScvIdx) const
    { return jacobian_[domainScvIdx][rangeScvIdx]; }

    /*!
     * \brief Returns the local residual of a sub-control volume.
     *
     * \param dofIdx The local index of the sub control volume
     */
    const ScalarVectorBlock &residual(int dofIdx) const
    { return residual_[dofIdx]; }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    const Simulator &simulator_() const
    { return *simulatorPtr_; }
    const Problem &problem_() const
    { return simulatorPtr_->problem(); }
    const Model &model_() const
    { return simulatorPtr_->model(); }

    /*!
     * \brief Returns the numeric difference method which is applied.
     */
    static int numericDifferenceMethod_()
    { return EWOMS_GET_PARAM(TypeTag, int, NumericDifferenceMethod); }

    /*!
     * \brief Resize all internal attributes to the size of the
     *        element.
     */
    void resize_(const ElementContext &elemCtx)
    {
        int numDof = elemCtx.numDof(/*timeIdx=*/0);
        int numPrimaryDof = elemCtx.numPrimaryDof(/*timeIdx=*/0);

        residual_.resize(numDof);
        jacobian_.setSize(numDof, numPrimaryDof);
    }

    /*!
     * \brief Reset the all relevant internal attributes to 0
     */
    void reset_(const ElementContext &elemCtx)
    {
        residual_ = 0.0;
        jacobian_ = 0.0;
    }

    /*!
     * \brief Updates the current local Jacobian matrix with the partial derivatives of
     *        all equations in regard to the primary variable 'pvIdx' at vertex 'dofIdx'.
     */
    void updateLocalLinearization_(const ElementContext &elemCtx,
                                   int primaryDofIdx)
    {
        const auto& resid = localResidual_.residual();

        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
            residual_[primaryDofIdx][eqIdx] = resid[primaryDofIdx][eqIdx].value();

        int numDof = elemCtx.numDof(/*timeIdx=*/0);
        for (int dofIdx = 0; dofIdx < numDof; dofIdx++) {
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++) {
                for (int pvIdx = 0; pvIdx < numEq; pvIdx++) {
                    // A[dofIdx][primaryDofIdx][eqIdx][pvIdx] is the partial derivative of
                    // the residual function 'eqIdx' for the degree of freedom 'dofIdx' with
                    // regard to the primary variable 'pvIdx' of the degree of freedom
                    // 'primaryDofIdx'
                    jacobian_[dofIdx][primaryDofIdx][eqIdx][pvIdx] =
                        resid[dofIdx][eqIdx].derivative(pvIdx);
                    Valgrind::CheckDefined(jacobian_[dofIdx][primaryDofIdx][eqIdx][pvIdx]);
                }
            }
        }
    }

    Simulator *simulatorPtr_;
    Model *modelPtr_;

    ElementContext *internalElemContext_;

    LocalResidual localResidual_;

    ScalarLocalBlockVector residual_;
    ScalarLocalBlockMatrix jacobian_;
};

} // namespace Ewoms

#endif
