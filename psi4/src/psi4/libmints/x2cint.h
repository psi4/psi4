/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmints_x2cint_h_
#define _psi_src_lib_libmints_x2cint_h_

#include "psi4/libmints/matrix.h"

#define X2CDEBUG 0

namespace psi {

/*! \ingroup MINTS
 *  \class X2CInt
 *  \brief Computes the 1e-X2C kinetic and potential integrals.
 */
class X2CInt
{
public:
    X2CInt();
    ~X2CInt();

    /*! @{
     * Computes the X2C kinetic and potential integrals
     * @param S Shared matrix object that will hold the X2C overlap integrals.
     * @param T Shared matrix object that will hold the X2C kinetic energy integrals.
     * @param V Shared matrix object that will hold the X2C potential energy integrals.
     * @param options an Options object used to read basis set information.
     */
    void compute(std::shared_ptr<BasisSet> basis,
                 std::shared_ptr<BasisSet> x2c_basis, SharedMatrix S,
                 SharedMatrix T, SharedMatrix V);
    /*! @} */

private:
    /// The name of the basis set
    std::string basis_;
    /// The name of the basis set
    std::string x2c_basis_;
    /// Do basis set projection?
    bool do_project_;

    /// Integral factory
    std::shared_ptr<IntegralFactory> integral_;

    /// The basis used to solve the Dirac equation
    std::shared_ptr<BasisSet> aoBasis_;
    /// The basis onto which we project the final FW (Foldy-Wouthuysen) Hamiltonian
    /// This is used only if we project (do_project)
    std::shared_ptr<BasisSet> aoBasis_contracted_;

    /// Matrix factory for matrices of dimension 2 nbf x 2 nbf
    std::shared_ptr<MatrixFactory> ssFactory_;
    /// Matrix factory for matrices of dimension nbf x nbf
    std::shared_ptr<MatrixFactory> soFactory_;
    /// Dimension of the orbital basis
    Dimension nsopi_;
    /// Dimension of the constracted orbital basis
    Dimension nsopi_contracted_;

    // Matrices in the orbital basis
    /// The overlap matrix in the orbital basis
    SharedMatrix sMat;
    /// The kinetic energy matrix in the orbital basis
    SharedMatrix tMat;
    /// The potential energy matrix in the orbital basis
    SharedMatrix vMat;
    /// The spin-free relativistic potential (W = pVp) matrix in the orbital basis
    SharedMatrix wMat;
    /// The X matrix
    SharedMatrix xMat;
    /// The R matrix
    SharedMatrix rMat;
    /// The XR matrix
    SharedMatrix xrMat;
    /// The X2C overlap matrix
    SharedMatrix S_x2c_;
    /// The X2C kinetic energy matrix
    SharedMatrix T_x2c_;
    /// The X2C potential energy matrix
    SharedMatrix V_x2c_;

    // Matrices in the orbital basis doubled
    /// The four-component Hamiltonian of the modified Dirac equation
    SharedMatrix dMat;
    /// The four-component overlap matrix of the modified Dirac equation
    SharedMatrix SXMat;
    /// Eigenvectors of the modified Dirac equation
    SharedMatrix C_LS_Mat;
    /// Eigenvalues of the modified Dirac equation
    SharedVector E_LS_Mat;

    /// Setup the basis objects, integral factories, etc.
    void setup(std::shared_ptr<BasisSet> basis,
               std::shared_ptr<BasisSet> x2c_basis);
    /// Compute the S, T, V, and W integrals
    void compute_integrals();
    /// Compute the Hamiltonian and overlap matrices of the modified Dirac equation
    void form_dirac_h();
    /// Diagonalize the modified Dirac equation
    void diagonalize_dirac_h();
    /// Form the matrix X
    void form_X();
    /// Form the matrices R and XR
    void form_R();
    /// Form the FW Hamiltonian for positive energy states
    void form_h_FW_plus();
    /// Write the FW Hamiltonian for positive energy states
    void write_integrals_to_disk();
    /// Test the FW Hamiltonian
    void test_h_FW_plus();
    /// Basis set projection
    void project();
};

}

#endif // _psi_src_lib_libmints_x2cint_h_
