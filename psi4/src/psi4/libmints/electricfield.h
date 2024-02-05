/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmints_electricfield_h_
#define _psi_src_lib_libmints_electricfield_h_

#include <vector>
#include "typedefs.h"
#include "psi4/libmints/vector3.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/integral.h"

namespace psi {
class Molecule;

/*! \ingroup MINTS
 *  \class ElectricFieldInt
 *  \brief Computes electric field integrals.
 *
 *  Use an IntegralFactory to create this object.
 */
class ElectricFieldInt : public OneBodyAOInt {

   public:
    //! Constructor. Do not call directly use an IntegralFactory.
    ElectricFieldInt(std::vector<SphericalTransform> &, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                     int deriv = 0);
    //! Virtual destructor
    ~ElectricFieldInt() override;

    //! Does the method provide first derivatives?
    bool has_deriv1() override { return false; }

    static Vector3 nuclear_contribution(const Vector3 &origin, std::shared_ptr<Molecule> mol);

    /// Computes all integrals and stores them in result by default this method throws
    void compute(std::vector<SharedMatrix>& result) override;

    /** Compute field integrals at coords with a functor to obtain
    a) the expectation value of the electric field at all coords (ContractOverDensityFieldFunctor)
    b) the induction operator matrix by contraction with dipoles (ContractOverDipolesFunctor)
    */
    template <typename ContractionFunctor>
    void compute_with_functor(ContractionFunctor functor, SharedMatrix coords);

    void set_origin(const Vector3& _origin) override;
};

class ContractOverDipolesFunctor {
    /**
    Contracts the electric field integrals with a matrix of dipoles,
    result is added to the matrix F
    */
   protected:
    /// Pointer to the resulting operator matrix
    double **pF_;
    /// The array of dipoles
    double **dipoles_;

   public:
    ContractOverDipolesFunctor(SharedMatrix dipoles, SharedMatrix F) : pF_(F->pointer()), dipoles_(dipoles->pointer()) {
        if (F->rowdim() != F->coldim()) throw PSIEXCEPTION("Invalid Fock matrix in ContractOverCharges");
        if (dipoles->coldim() != 3) throw PSIEXCEPTION("Dipole matrix must have 3 columns.");
    }

    void operator()(int bf1, int bf2, int center, double integralx, double integraly, double integralz) {
        pF_[bf1][bf2] +=
            integralx * dipoles_[center][0] + integraly * dipoles_[center][1] + integralz * dipoles_[center][2];
    }
};

class ContractOverDensityFieldFunctor {
    /**
    Contracts the electric field integrals with a density matrix D
    and writes the result (expectation value) to field
    */
   protected:
    /// The array of the density matrix
    double **pD_;
    /// The array of electric fields
    double **field_;

   public:
    ContractOverDensityFieldFunctor(SharedMatrix field, SharedMatrix D) : pD_(D->pointer()), field_(field->pointer()) {
        if (D->rowdim() != D->coldim()) throw PSIEXCEPTION("Invalid density matrix in ContractOverDensityFieldFunctor");
        if (field->coldim() != 3) throw PSIEXCEPTION("Field matrix must have 3 columns.");
    }
    void operator()(int bf1, int bf2, int center, double integralx, double integraly, double integralz) {
        field_[center][0] += pD_[bf1][bf2] * integralx;
        field_[center][1] += pD_[bf1][bf2] * integraly;
        field_[center][2] += pD_[bf1][bf2] * integralz;
    }
};

}  // namespace psi

#endif
