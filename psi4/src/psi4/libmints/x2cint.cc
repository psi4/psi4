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

#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"
#include "psi4/physconst.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/libmints/rel_potential.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/x2cint.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/libmints/basisset.h"

namespace psi {

X2CInt::X2CInt()
{
}

X2CInt::~X2CInt()
{
}

void X2CInt::compute(std::shared_ptr<BasisSet> basis,
                     std::shared_ptr<BasisSet> x2c_basis, SharedMatrix S,
                     SharedMatrix T, SharedMatrix V) {
    // tstart();
    setup(basis, x2c_basis);
    compute_integrals();
    form_dirac_h();
    diagonalize_dirac_h();
    form_X();
    form_R();
    form_h_FW_plus();

    if(do_project_){
        project();
    }

    test_h_FW_plus();

    S->copy(S_x2c_);
    T->copy(T_x2c_);
    V->copy(V_x2c_);
    // tstop();
}

void X2CInt::setup(std::shared_ptr<BasisSet> basis,
                   std::shared_ptr<BasisSet> x2c_basis) {
    outfile->Printf(
        "         "
        "------------------------------------------------------------");
    outfile->Printf("\n         Spin-Free X2C Integrals at the One-Electron Level (SFX2C-1e)");
    outfile->Printf("\n                 by Prakash Verma and Francesco A. Evangelista");
    outfile->Printf("\n         ------------------------------------------------------------\n");

    // basis set constructed from BASIS.
    basis_ = basis->name();
    aoBasis_contracted_ = basis;

    // basis set constructed from BASIS_RELATIVISTIC.
    x2c_basis_ = x2c_basis->name();
    aoBasis_ = x2c_basis;
    do_project_ = true;

    // Print X2C options
    outfile->Printf("\n  ==> X2C Options <==\n");
    outfile->Printf("\n    Computational Basis: %s",basis_.c_str());
    outfile->Printf("\n    X2C Basis: %s",x2c_basis_.c_str());
    outfile->Printf("\n    The X2C Hamiltonian will be computed in the X2C Basis\n");

    // The integral factory oversees the creation of integral objects
    integral_ = std::shared_ptr<IntegralFactory>(new IntegralFactory(aoBasis_, aoBasis_, aoBasis_, aoBasis_));

    // Create an SO basis...we need the point group for this part.
    // SOBasisSet object for the computational basis
    std::shared_ptr<SOBasisSet> soBasis(new SOBasisSet(aoBasis_, integral_));

    // Obtain the dimension object to initialize the factory.
    nsopi_ = soBasis->dimension();
    nsopi_contracted_ = nsopi_;

    // Create a Dimension object for the spinors
    Dimension nsspi = nsopi_ + nsopi_;

    // Matrix factory for matrices of dimension nbf x nbf
    soFactory_ = std::shared_ptr<MatrixFactory>(new MatrixFactory);
    soFactory_->init_with(nsopi_,nsopi_);

    // Matrix factory for matrices of dimension 2 nbf x 2 nbf
    ssFactory_ = std::shared_ptr<MatrixFactory>(new MatrixFactory);
    ssFactory_->init_with(nsspi,nsspi);
}

void X2CInt::compute_integrals()
{
    // Create the integral objects
    std::shared_ptr<OneBodySOInt> sOBI(integral_->so_overlap());
    std::shared_ptr<OneBodySOInt> tOBI(integral_->so_kinetic());
    std::shared_ptr<OneBodySOInt> vOBI(integral_->so_potential());
    std::shared_ptr<OneBodySOInt> wOBI(integral_->so_rel_potential());

    // Form the one-electron integral matrices from the matrix factory
    sMat = SharedMatrix(soFactory_->create_matrix("Overlap"));
    tMat = SharedMatrix(soFactory_->create_matrix("Kinetic"));
    vMat = SharedMatrix(soFactory_->create_matrix("Potential"));
    wMat = SharedMatrix(soFactory_->create_matrix("Relativistic Potential"));

    // Compute the one electron integrals, telling each object where to store the result
    sOBI->compute(sMat);
    tOBI->compute(tMat);
    vOBI->compute(vMat);
    wOBI->compute(wMat);

#if X2CDEBUG
    sMat->print();
    tMat->print();
    vMat->print();
    wMat->print();
//    sMat_cont->print();
#endif
}

void X2CInt::form_dirac_h()
{
    // Form the Dirac Hamiltonian
    //      | V       T       |
    //  d = |                 |
    //      | T  1/4c^2 W - T |
    //
    // and the overlap matrix
    //         | S        0 |
    // SXMat=  |            |
    //         | 0  T/2c**2 |
    //

    dMat = SharedMatrix(ssFactory_->create_matrix("Dirac Hamiltonian"));
    SXMat = SharedMatrix(ssFactory_->create_matrix("SX Hamiltonian"));

    for (int h = 0; h < dMat->nirrep(); ++h){
        int nsopi_h = dMat->rowdim(h) / 2;
        for (int p = 0; p < nsopi_h; ++p){
            for (int q = 0; q < nsopi_h; ++q){
                double Vpq = vMat->get(h,p,q);
                double Tpq = tMat->get(h,p,q);
                double Wpq = wMat->get(h,p,q);
                double Spq = sMat->get(h,p,q);


                // Set the large-large component block
                SXMat->set(h,p,q,Spq);
                // Set the small-small component block
                SXMat->set(h,p + nsopi_h,q + nsopi_h,0.5 *Tpq / (pc_c_au * pc_c_au));

                // Set the large-large component block
                dMat->set(h,p,q,Vpq);
                // Set the small-large component block
                dMat->set(h,p + nsopi_h,q,Tpq);
                // Set the large-small component block
                dMat->set(h,p,q + nsopi_h,Tpq);
                // Set the small-small component block
                dMat->set(h,p + nsopi_h,q + nsopi_h,0.25 * Wpq/ (pc_c_au * pc_c_au)  - Tpq);
            }
        }
    }

#if X2CDEBUG
    dMat->print();
#endif
}

void X2CInt::diagonalize_dirac_h()
{
    /* Prakash
     * C_LS_Mat   Eigenvector that contains both large and small component
     * E_LS_Mat   EigenValues  if Dirac Hamiltonian
     */

    C_LS_Mat = SharedMatrix(ssFactory_->create_matrix("Dirac EigenVectors"));
    E_LS_Mat = SharedVector(new Vector("Dirac EigenValues",C_LS_Mat->rowspi()));
    SharedMatrix dtmpMat(ssFactory_->create_matrix("Dirac tmp Hamiltonian"));


    /* Prakash
      * Step 3.
      * Diagonalize the Dirac Hamiltonian with X2C 1e Hamiltonian
      */

    SXMat->power(-1.0/2.0);
    dMat->transform(SXMat);
    dMat->diagonalize(dtmpMat,E_LS_Mat);                  //diagonalize Dfock
    C_LS_Mat->gemm(false, false, 1.0, SXMat, dtmpMat, 0.0 );  // C = X C'

#if X2CDEBUG
    C_LS_Mat->print();
    E_LS_Mat->print();
#endif
}

void X2CInt::form_X()
{
    /* Prakash
     * Step 4.
     * Divide positive energy part C_LS_Mat into two parts C_Large and C_Small
     * X can be constructed as X * C_Large = C_Small
     * more appropriately as   X = C_Small *  C_Large {^-1}
     */
    SharedMatrix clMat(soFactory_->create_matrix("Large EigenVectors"));
    SharedMatrix csMat(soFactory_->create_matrix("Small EigenVectors"));
    xMat = SharedMatrix(soFactory_->create_matrix("X matrix"));

    /*
     * collect the correct Matrix element from Eivenvectos
     */
    for (int h = 0; h < clMat->nirrep(); ++h){
        int nsopi_h = clMat->rowdim(h);
        for (int p = 0; p < nsopi_h; ++p){
            for (int q = 0; q < nsopi_h; ++q){
                double Lpq = C_LS_Mat->get(h,p,q + nsopi_h);
                double Spq = C_LS_Mat->get(h,p + nsopi_h,q + nsopi_h);
                // Set the large-large component block
                clMat->set(h,p,q,Lpq);
                // Set the small-small component block
                csMat->set(h,p,q,Spq);
            }
        }
    }

    /*
     * Take the inverse of large or upper component
     */
    clMat->general_invert();           // Find C_Large inverse

    /*
     * FORM X = C_small * (C_large)^{-1}
     */
    xMat->gemm(false, false, 1.0, csMat, clMat, 0.0 );

#if X2CDEBUG
    xMat->print();
#endif
}

void X2CInt::form_R()
{
    /* Prakash
     * step 5
     * S_{tilda} = S + 1/2c**2  X^{dagger} T X
     * since we need X^{dagger} T later so we keep that around for little longer
     */


    /*
     * FORM X^ T X
     * where X = C_small * ( C_large ) ^{-1}
     * and  T is the kinetic energy
     */
    SharedMatrix XTX(soFactory_->create_matrix("XTX matrix"));
    XTX->transform(xMat,tMat,xMat);

    /*
     * FORM  S_{tilda} = S + (X^ T X) / 2c**2
     * S is the overlap matrix
     */
    SharedMatrix S_tilde(soFactory_->create_matrix("S tilde matrix"));
    XTX->scale(1.0/(2.0*pc_c_au * pc_c_au));             // scale by 1/2c**2
    S_tilde->copy(sMat);                                 // S_tilda = S + X^ T X
    S_tilde->add(XTX);

#if X2CDEBUG
    S_tilde->print();
#endif

    /*Prakash
      * Step 6
      * R = S^{-1/2} (S ^{-1/2} S_tilda S^{-1/2})^{-1/2} S^{1/2}
      */

    SharedMatrix  S_inv_half(soFactory_->create_matrix("Eigenvector S matrix"));
    SharedMatrix  sTmp1(soFactory_->create_matrix("S tmp1 matrix"));
    SharedMatrix  sTmp2(soFactory_->create_matrix("S tmp2 matrix"));

    /*
     * FORM  S^{-1/2}
     * Matrix = U diag(Matrix) U^
     */

    S_inv_half->copy(sMat);
    S_inv_half->power(-1.0/2.0);
#if X2CDEBUG
    S_inv_half->print();
#endif

    /*
     * FORM  (S ^{-1/2} S_tilda S^{-1/2})^{-1/2}
     */
//    clMat->gemm(false, false, 1.0, S_inv_half, S_tilde, 0.0 );  // S^{-1/2} S_tilda
//    sTmp->gemm(false, false, 1.0, clMat, S_inv_half,0.0);      // S^{-1/2} S_tilda S^{-1/2}
    // TOCHECK -> this line below should take care of the two above.
    sTmp1->transform(S_tilde,S_inv_half);
    sTmp1->power(-1.0/2.0);


    /*
     * S^{-1/2} (S ^{-1/2} S_tilda S^{-1/2})^{-1/2} S^{1/2}
     */

    sTmp2->gemm(false, false,  1.0, S_inv_half, sTmp1, 0.0);    // S^{-1/2} * (sTmp1)
    S_inv_half->general_invert();                              // S^{1/2}
    rMat = SharedMatrix(soFactory_->create_matrix("R matrix"));
    rMat->gemm(false, false, 1.0, sTmp2, S_inv_half, 0.0 );   // R=S^{-1/2} (S_inv S_tilda S_inv)^{-1/2}S^{1/2}
#if X2CDEBUG
    rMat->print();
#endif

    /*
     * FORM XR matrix
     */

    xrMat = SharedMatrix(soFactory_->create_matrix("XR matrix"));
    xrMat->gemm(false, false, 1.0, xMat, rMat, 0.0 );     // XR = X R matrix
#if X2CDEBUG
    xrMat->print();
#endif
}

void X2CInt::form_h_FW_plus()
{
    // Check if the matrices are allocated and have the correct size
    S_x2c_ = SharedMatrix(soFactory_->create_matrix(PSIF_SO_S));
    T_x2c_ = SharedMatrix(soFactory_->create_matrix(PSIF_SO_T));
    V_x2c_ = SharedMatrix(soFactory_->create_matrix(PSIF_SO_V));

    /*
     * S^{FW}_{+} = S
     */

    S_x2c_->copy(sMat);

    /*Prakash
     * step 7
     * construct h^{FW}_{+}
     *                   =      R^ T  XR
     *                     + (XR)^ T  R
     *                     - (XR)^ T  XR
     *                     +    R^ V  R
     *                     + (XR)^ W' XR
     *
     * where W' is the scaled version of W i.e W/4c**2
     */

    SharedMatrix  Tmp1(soFactory_->create_matrix("Temporary matrix"));


    /*
     *   R^ T XR + (XR)^ T R - (XR)^ T (XR)
     */

    Tmp1->transform(rMat,tMat,xrMat);
    T_x2c_->copy(Tmp1);

    Tmp1->transpose_this();
    T_x2c_->add(Tmp1);

    Tmp1->zero();
    Tmp1->transform(tMat,xrMat);
    T_x2c_->subtract(Tmp1);

    /*
     *   R^ V R
     */
    Tmp1->zero();
    Tmp1->transform(vMat,rMat);
    V_x2c_->copy(Tmp1);

    /*
     *   (XR)^ W' XR
     */
    Tmp1->zero();
    Tmp1->transform(wMat,xrMat);
    Tmp1->scale(1.0/(4.0*pc_c_au * pc_c_au));
    V_x2c_->add(Tmp1);

#if X2CDEBUG
    S_x2c_->print();
    T_x2c_->print();
    V_x2c_->print();
#endif
}

void X2CInt::write_integrals_to_disk()
{
    /*
     *  Write T and V to disk
     */
    S_x2c_->save(_default_psio_lib_,PSIF_OEI);
    T_x2c_->save(_default_psio_lib_,PSIF_OEI);
    V_x2c_->save(_default_psio_lib_,PSIF_OEI);
}

void X2CInt::project()
{
    // Integral factory for the BASIS/X2C_BASIS mixed basis
    std::shared_ptr<IntegralFactory> integral_contracted(new IntegralFactory(aoBasis_contracted_, aoBasis_, aoBasis_, aoBasis_));

    std::shared_ptr<SOBasisSet> soBasis_contracted(new SOBasisSet(aoBasis_contracted_, integral_contracted));

    nsopi_contracted_ = soBasis_contracted->dimension();

    std::shared_ptr<MatrixFactory> soFactory_contracted(new MatrixFactory);
    soFactory_contracted->init_with(nsopi_contracted_,nsopi_);

    // Form the overlap matrix in the BASIS/X2C_BASIS basis
    std::shared_ptr<OneBodySOInt> sOBI_cu(integral_contracted->so_overlap());
    SharedMatrix S_cu(soFactory_contracted->create_matrix("Overlap"));
    sOBI_cu->compute(S_cu);

    SharedMatrix S_inv = sMat->clone();
    S_inv->general_invert();

    SharedMatrix D(new Matrix("D",nsopi_,nsopi_contracted_));
    //  Form D = S_uu^{-1} S_uc.  Notice that we transpose S_cu
    D->gemm(false,true,1.0,S_inv,S_cu,0.0);

    S_x2c_->transform(D);
    T_x2c_->transform(D);
    V_x2c_->transform(D);
}

void X2CInt::test_h_FW_plus()
{
    /*
     * Diagonalize the Hamiltonian
     */
    SharedMatrix Evec_x2c(S_x2c_->clone());
    SharedVector Eval_x2c(new Vector("Eigenvalues of h_FW^{+}",sMat->rowspi()));
    SharedMatrix S_inv_half(S_x2c_->clone());
    SharedMatrix H_x2c = T_x2c_->clone();
    H_x2c->add(V_x2c_);
    S_inv_half->power(-0.5);
    H_x2c->transform(S_inv_half);
    H_x2c->diagonalize(Evec_x2c,Eval_x2c);

    double sum = 0.0;
    for (int h = 0; h < dMat->nirrep(); ++h){
        int nsopi_h = dMat->rowdim(h) / 2;
        int maxp = nsopi_contracted_[h];
        if (maxp != nsopi_h){
            outfile->Printf("\n    Comparing only %d out of %d elements of H_Dirac\n",maxp,nsopi_h);
        }
        for (int p = 0; p < maxp; ++p){
            double eval_dirac = E_LS_Mat->get(h,nsopi_h + p);
            double eval_x2c = Eval_x2c->get(h,p);
            sum += std::fabs(eval_dirac-eval_x2c);
//            outfile->Printf("\n %d,%2d = %20.12f %20.12f %20.12f",h,p,eval_dirac,eval_x2c,eval_dirac-eval_x2c);
        }
    }
    outfile->Printf("\n    The 1-norm of |H_X2C - H_Dirac| is: %.12f\n",sum);
    if (sum > 1.0e-6){
        outfile->Printf("\n    WARNING: The X2C and Dirac Hamiltonians have substatially different eigenvalues!\n");
        if (do_project_){
            outfile->Printf("             This is probably caused by the recontraction of the basis set.\n\n");
        }else{
            outfile->Printf("             There is something wrong with the X2C module.\n\n");
        }
        outfile->Flush();
    }
}

}
