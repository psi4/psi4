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

#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/molecule.h"
#include "psi4/physconst.h"
#include "psi4/libmints/pointgrp.h"
#include "occwave.h"


using namespace std;


namespace psi{ namespace occwave{

void OCCWave::ekt_ip()
{

//outfile->Printf("\n ekt_ip is starting... \n");
//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
     // Memory allocation
     SharedMatrix GFock_primeA = std::shared_ptr<Matrix>(new Matrix("Alpha GF prime", nirrep_, nmopi_, nmopi_));
     SharedMatrix GFock_copyA = std::shared_ptr<Matrix>(new Matrix("Alpha GF copy", nirrep_, nmopi_, nmopi_));
     SharedMatrix g1symm_copyA = std::shared_ptr<Matrix>(new Matrix("Alpha OPDM copy", nirrep_, nmopi_, nmopi_));
     SharedMatrix g1HalfA = std::shared_ptr<Matrix>(new Matrix("g^-1/2", nirrep_, nmopi_, nmopi_));
     SharedMatrix UvecA = std::shared_ptr<Matrix>(new Matrix("UvecA", nirrep_, nmopi_, nmopi_));
     SharedMatrix Uvec_primeA = std::shared_ptr<Matrix>(new Matrix("Uvec_primeA", nirrep_, nmopi_, nmopi_));
     SharedMatrix PSA = std::shared_ptr<Matrix>(new Matrix("Alpha pole strength matrix", nirrep_, nmopi_, nmopi_));
     SharedMatrix gc_transA = std::shared_ptr<Matrix>(new Matrix("Alpha C'*g", nirrep_, nmopi_, nmopi_));
     SharedMatrix tempA = std::shared_ptr<Matrix>(new Matrix("Alpha temp", nirrep_, nmopi_, nmopi_));
     SharedVector Diag_g1A = std::shared_ptr<Vector>(new Vector("Diag OO-block OPDM", nirrep_, nmopi_));
     SharedVector ps_vecA = std::shared_ptr<Vector>(new Vector("alpha pole strength vector", nirrep_, nmopi_));
     SharedVector eorbA = std::shared_ptr<Vector>(new Vector("eorbA", nirrep_, nmopi_));

     // For Non-OO methods
     if (orb_opt_ == "FALSE" && reference_ == "RESTRICTED") GFock->scale(0.5);
     else if (orb_opt_ == "FALSE" && reference_ == "UNRESTRICTED") {
              GFockA->scale(0.5);
              GFockB->scale(0.5);
     }

     // Make sure GFM is symmetric
     if (sym_gfm_ == "TRUE" && reference_ == "RESTRICTED") {
         SharedMatrix temp(GFock->transpose());
         GFock_copyA->copy(GFock);
         GFock_copyA->add(temp);
         GFock_copyA->scale(0.5);
         GFock->copy(GFock_copyA);
         if (print_ >=2 ) GFock->print();

         // Symm OPDM
         SharedMatrix temp2(g1symm->transpose());
         g1symm_copyA->copy(g1symm);
         g1symm_copyA->add(temp2);
         g1symm_copyA->scale(0.5);
         g1symm->copy(g1symm_copyA);
         if (print_ >=2 ) g1symm->print();
     }


     // Diagonalize OPDM
     UvecA->zero();
     Diag_g1A->zero();
     if (reference_ == "RESTRICTED") g1symm->diagonalize(UvecA, Diag_g1A);
     else if (reference_ == "UNRESTRICTED") g1symmA->diagonalize(UvecA, Diag_g1A);

     /*
     // Print g^(-1/2)
     for (int h = 0; h < nirrep_; ++h) {
          for (int i = 0; i < nmopi_[h]; ++i) {
               outfile->Printf("\t h, i, Diag_g1A: %3d %3d %20.14f \n", h, i, Diag_g1A->get(h, i));

          }
     }
     */

     // Make sure all eigenvalues are positive
     for (int h = 0; h < nirrep_; ++h) {
          for (int i = 0; i < nmopi_[h]; ++i) {
               if (Diag_g1A->get(h, i) < 0.0) Diag_g1A->set(h, i, -1.0*Diag_g1A->get(h, i));
          }
     }

     /*
     // Again Print g^(-1/2)
     for (int h = 0; h < nirrep_; ++h) {
          for (int i = 0; i < nmopi_[h]; ++i) {
               outfile->Printf("\t h, i, Diag_g1A: %3d %3d %20.14f \n", h, i, Diag_g1A->get(h, i));

          }
     }
     */

     // Form g^(-1/2)
     for (int h = 0; h < nirrep_; ++h) {
          for (int i = 0; i < nmopi_[h]; ++i) {
               Diag_g1A->set(h, i, 1/sqrt(Diag_g1A->get(h, i)));
          }
     }

     g1HalfA->zero();
     for (int h = 0; h < nirrep_; ++h) {
          for (int i = 0; i < nmopi_[h]; ++i) {
               g1HalfA->set(h, i, i, Diag_g1A->get(h, i));
          }
     }

     tempA->zero();
     tempA->gemm(false, true, 1.0, g1HalfA, UvecA, 0.0);
     g1HalfA->gemm(false, false, 1.0, UvecA, tempA, 0.0);

     // Build GFock prime matrix
     GFock_primeA->zero();
     if (reference_ == "RESTRICTED") tempA->gemm(true, false, 1.0, g1HalfA, GFock, 0.0);
     else if (reference_ == "UNRESTRICTED") tempA->gemm(true, false, 1.0, g1HalfA, GFockA, 0.0);
     GFock_primeA->gemm(false, false, 1.0, tempA, g1HalfA, 0.0);

     // Diagonalize GFock to get orbital energies
     eorbA->zero();
     Uvec_primeA->zero();
     GFock_primeA->diagonalize(Uvec_primeA, eorbA);
     UvecA->gemm(false, false, 1.0, g1HalfA, Uvec_primeA, 0.0);

     // Pole strength
     PSA->zero();
     gc_transA->zero();
     if (reference_ == "RESTRICTED") tempA->gemm(false, false, 1.0, g1symm, UvecA, 0.0);
     else if (reference_ == "UNRESTRICTED") tempA->gemm(false, false, 1.0, g1symmA, UvecA, 0.0);
     gc_transA = tempA->transpose();
     PSA->gemm(false, false, 1.0, gc_transA, tempA, 0.0);
     ps_vecA->zero();
     for (int h = 0; h < nirrep_; ++h) {
          for (int i = 0; i < nmopi_[h]; ++i) {
               ps_vecA->set(h, i, PSA->get(h, i, i));
          }
     }
     if (reference_ == "RESTRICTED") ps_vecA->scale(0.5);

    // Sort pole strength
    Array1d *evals_A = new Array1d("Alpha ORB C1", nmo_);
    Array1d *ps_vec2A = new Array1d("Sorted Pole strength", nmo_);
    Array1i *irrep_A = new Array1i("IrrepA", nmo_);
    evals_A->zero();
    ps_vec2A->zero();
    irrep_A->zero();


    // Copy ps vec
    int count = 0;
    for (int h = 0; h < nirrep_; ++h){
	 for (int i = 0; i < nmopi_[h]; ++i){
              evals_A->set(count, eorbA->get(h, i));
              ps_vec2A->set(count, ps_vecA->get(h, i));
              irrep_A->set(count, h);
              count++;
          }
    }

    // Sort to descending order
    for (int i = 0; i < nmo_; ++i) {
         for(int j = nmo_-1; j > i; --j) {
             if (ps_vec2A->get(j-1) < ps_vec2A->get(j)) {
                 double dum = evals_A->get(j-1);
                 evals_A->set(j-1, evals_A->get(j));
                 evals_A->set(j, dum);

                 int dum2 = irrep_A->get(j-1);
                 irrep_A->set(j-1, irrep_A->get(j));
                 irrep_A->set(j, dum2);

                 double dum3 = ps_vec2A->get(j-1);
                 ps_vec2A->set(j-1, ps_vec2A->get(j));
                 ps_vec2A->set(j, dum3);
             }
         }
    }

    // Re-Sort occupied orbitals to energy order
    Array1d *eoccA = new Array1d("Alpha occupied orb", nooA);
    Array1d *ps_occA = new Array1d("occupied Pole strength", nooA);
    Array1i *irrep_occA = new Array1i("occupied IrrepA", nooA);
    eoccA->zero();
    ps_occA->zero();
    irrep_occA->zero();


    // Copy
    for (int i = 0; i < nooA; ++i) {
         eoccA->set(i, evals_A->get(i));
         ps_occA->set(i, ps_vec2A->get(i));
         irrep_occA->set(i, irrep_A->get(i));
    }

    // Sort to ascending order
    for (int i = 0; i < nooA; ++i) {
         for(int j = nooA-1; j > i; --j) {
             if (eoccA->get(j-1) > eoccA->get(j)) {
                 double dum = eoccA->get(j-1);
                 eoccA->set(j-1, eoccA->get(j));
                 eoccA->set(j, dum);

                 int dum2 = irrep_occA->get(j-1);
                 irrep_occA->set(j-1, irrep_occA->get(j));
                 irrep_occA->set(j, dum2);

                 double dum3 = ps_occA->get(j-1);
                 ps_occA->set(j-1, ps_occA->get(j));
                 ps_occA->set(j, dum3);
             }
         }
    }

    // Print IPs
    outfile->Printf("\n\tEKT-OCC Ionization Potentials (Alpha Spin Case) \n");
    outfile->Printf("\t------------------------------------------------------------------- \n");


    Molecule& mol = *reference_wavefunction_->molecule().get();
    CharacterTable ct = mol.point_group()->char_table();
    string pgroup = mol.point_group()->symbol();

 // print alpha IPs
 if (print_ < 2) {
    outfile->Printf( "\tState    Symmetry   -IP (a.u.)       IP (eV)        Pole Strength \n");
    outfile->Printf("\t------------------------------------------------------------------- \n");

    for (int i = 0; i < nooA; ++i){
         int h = irrep_occA->get(i);
	 outfile->Printf("\t%3d %10s %15.6f %15.6f %15.6f \n", i+1, ct.gamma(h).symbol(),
                          eoccA->get(i), -eoccA->get(i)*pc_hartree2ev, ps_occA->get(i));

    }
    outfile->Printf("\t------------------------------------------------------------------- \n");

 }// end if

 else if (print_ >= 2) {
    outfile->Printf( "\tState    Symmetry   -IP (a.u.)       IP (eV)        Pole Strength \n");
    outfile->Printf("\t------------------------------------------------------------------- \n");

    for (int i = 0; i < nmo_; ++i){
         int h = irrep_A->get(i);
	 outfile->Printf("\t%3d %10s %15.6f %15.6f %15.6f \n", i+1, ct.gamma(h).symbol(),
                          evals_A->get(i), -evals_A->get(i)*pc_hartree2ev, ps_vec2A->get(i));

    }
    outfile->Printf("\t------------------------------------------------------------------- \n");

 }// end else if

//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
if (reference_ == "UNRESTRICTED") {
     // Memory allocation
     SharedMatrix GFock_primeB = std::shared_ptr<Matrix>(new Matrix("Beta OO-block GF prime", nirrep_, nmopi_, nmopi_));
     SharedMatrix GFock_copyB = std::shared_ptr<Matrix>(new Matrix("Beta GF copy", nirrep_, nmopi_, nmopi_));
     SharedMatrix g1symm_copyB = std::shared_ptr<Matrix>(new Matrix("Alpha OPDM copy", nirrep_, nmopi_, nmopi_));
     SharedMatrix g1HalfB = std::shared_ptr<Matrix>(new Matrix("g^-1/2", nirrep_, nmopi_, nmopi_));
     SharedMatrix UvecB = std::shared_ptr<Matrix>(new Matrix("UvecB", nirrep_, nmopi_, nmopi_));
     SharedMatrix Uvec_primeB = std::shared_ptr<Matrix>(new Matrix("Uvec_primeB", nirrep_, nmopi_, nmopi_));
     SharedMatrix PSB = std::shared_ptr<Matrix>(new Matrix("Beta pole strength matrix", nirrep_, nmopi_, nmopi_));
     SharedMatrix gc_transB = std::shared_ptr<Matrix>(new Matrix("Beta C'*g", nirrep_, nmopi_, nmopi_));
     SharedMatrix tempB = std::shared_ptr<Matrix>(new Matrix("Beta temp", nirrep_, nmopi_, nmopi_));
     SharedVector Diag_g1B = std::shared_ptr<Vector>(new Vector("DiagA OO-block OPDM", nirrep_, nmopi_));
     SharedVector ps_vecB = std::shared_ptr<Vector>(new Vector("Beta pole strength vector", nirrep_, nmopi_));
     SharedVector eorbB = std::shared_ptr<Vector>(new Vector("eorbB", nirrep_, nmopi_));

     // Make sure GFM is symmetric
     if (sym_gfm_ == "TRUE") {
         // alpha spin
         SharedMatrix temp_A(GFockA->transpose());
         GFock_copyA->copy(GFockA);
         GFock_copyA->add(temp_A);
         GFock_copyA->scale(0.5);
         GFockA->copy(GFock_copyA);

         // beta spin
         SharedMatrix temp_B(GFockB->transpose());
         GFock_copyB->copy(GFockB);
         GFock_copyB->add(temp_B);
         GFock_copyB->scale(0.5);
         GFockB->copy(GFock_copyB);

         // Symm OPDM alpha
         SharedMatrix temp2A(g1symmA->transpose());
         g1symm_copyA->copy(g1symmA);
         g1symm_copyA->add(temp2A);
         g1symm_copyA->scale(0.5);
         g1symmA->copy(g1symm_copyA);
         if (print_ >=2 ) g1symmA->print();

         // Symm OPDM alpha
         SharedMatrix temp2B(g1symmB->transpose());
         g1symm_copyB->copy(g1symmB);
         g1symm_copyB->add(temp2B);
         g1symm_copyB->scale(0.5);
         g1symmB->copy(g1symm_copyB);
         if (print_ >=2 ) g1symmB->print();
     }

     // Diagonalize OPDM
     UvecB->zero();
     Diag_g1B->zero();
     g1symmB->diagonalize(UvecB, Diag_g1B);

     // Make sure all eigenvalues are positive
     for (int h = 0; h < nirrep_; ++h) {
          for (int i = 0; i < nmopi_[h]; ++i) {
               if (Diag_g1B->get(h, i) < 0.0) Diag_g1B->set(h, i, -1.0*Diag_g1B->get(h, i));
          }
     }

     // Form g^(-1/2)
     for (int h = 0; h < nirrep_; ++h) {
          for (int i = 0; i < nmopi_[h]; ++i) {
               Diag_g1B->set(h, i, 1/sqrt(Diag_g1B->get(h, i)));
          }
     }

     g1HalfB->zero();
     for (int h = 0; h < nirrep_; ++h) {
          for (int i = 0; i < nmopi_[h]; ++i) {
               g1HalfB->set(h, i, i, Diag_g1B->get(h, i));
          }
     }

     tempB->zero();
     tempB->gemm(false, true, 1.0, g1HalfB, UvecB, 0.0);
     g1HalfB->gemm(false, false, 1.0, UvecB, tempB, 0.0);

     // Build GFock prime matrix
     GFock_primeB->zero();
     tempB->gemm(true, false, 1.0, g1HalfB, GFockB, 0.0);
     GFock_primeB->gemm(false, false, 1.0, tempB, g1HalfB, 0.0);

     // Diagonalize GFock to get orbital energies
     eorbB->zero();
     Uvec_primeB->zero();
     GFock_primeB->diagonalize(Uvec_primeB, eorbB);
     UvecB->gemm(false, false, 1.0, g1HalfB, Uvec_primeB, 0.0);

     // Pole strength
     PSB->zero();
     gc_transB->zero();
     tempB->gemm(false, false, 1.0, g1symmB, UvecB, 0.0);
     gc_transB = tempB->transpose();
     PSB->gemm(false, false, 1.0, gc_transB, tempB, 0.0);
     ps_vecB->zero();
     for (int h = 0; h < nirrep_; ++h) {
          for (int i = 0; i < nmopi_[h]; ++i) {
               ps_vecB->set(h, i, PSB->get(h, i, i));
          }
     }

    // Sort pole strength
    Array1d *evals_B = new Array1d("Alpha ORB C1", nmo_);
    Array1d *ps_vec2B = new Array1d("Sorted Pole strength", nmo_);
    Array1i *irrep_B = new Array1i("IrrepB", nmo_);
    evals_B->zero();
    ps_vec2B->zero();
    irrep_B->zero();

    // Copy ps vec
    int count = 0;
    for (int h = 0; h < nirrep_; ++h){
	 for (int i = 0; i < nmopi_[h]; ++i){
              evals_B->set(count, eorbB->get(h, i));
              ps_vec2B->set(count, ps_vecB->get(h, i));
              irrep_B->set(count, h);
              count++;
          }
    }

    // Sort to descending order
    for (int i = 0; i < nmo_; ++i) {
         for(int j = nmo_-1; j > i; --j) {
             if (ps_vec2B->get(j-1) < ps_vec2B->get(j)) {
                 double dum = evals_B->get(j-1);
                 evals_B->set(j-1, evals_B->get(j));
                 evals_B->set(j, dum);

                 int dum2 = irrep_B->get(j-1);
                 irrep_B->set(j-1, irrep_B->get(j));
                 irrep_B->set(j, dum2);

                 double dum3 = ps_vec2B->get(j-1);
                 ps_vec2B->set(j-1, ps_vec2B->get(j));
                 ps_vec2B->set(j, dum3);
             }
         }
    }

    // Re-Sort occupied orbitals to energy order
    Array1d *eoccB = new Array1d("Beta occupied orb", nooB);
    Array1d *ps_occB = new Array1d("Beta occupied Pole strength", nooB);
    Array1i *irrep_occB = new Array1i("occupied IrrepB", nooB);
    eoccB->zero();
    ps_occB->zero();
    irrep_occB->zero();


    // Copy
    for (int i = 0; i < nooB; ++i) {
         eoccB->set(i, evals_B->get(i));
         ps_occB->set(i, ps_vec2B->get(i));
         irrep_occB->set(i, irrep_B->get(i));
    }

    // Sort to ascending order
    for (int i = 0; i < nooB; ++i) {
         for(int j = nooB-1; j > i; --j) {
             if (eoccB->get(j-1) > eoccB->get(j)) {
                 double dum = eoccB->get(j-1);
                 eoccB->set(j-1, eoccB->get(j));
                 eoccB->set(j, dum);

                 int dum2 = irrep_occB->get(j-1);
                 irrep_occB->set(j-1, irrep_occB->get(j));
                 irrep_occB->set(j, dum2);

                 double dum3 = ps_occB->get(j-1);
                 ps_occB->set(j-1, ps_occB->get(j));
                 ps_occB->set(j, dum3);
             }
         }
    }

    // Print IPs
    outfile->Printf("\n\tEKT-OCC Ionization Potentials (Beta Spin Case) \n");
    outfile->Printf("\t------------------------------------------------------------------- \n");


 // print alpha IPs
 if (print_ < 2) {
    outfile->Printf( "\tState    Symmetry   -IP (a.u.)       IP (eV)        Pole Strength \n");
    outfile->Printf("\t------------------------------------------------------------------- \n");

    for (int i = 0; i < nooB; ++i){
         int h = irrep_occB->get(i);
	 outfile->Printf("\t%3d %10s %15.6f %15.6f %15.6f \n", i+1, ct.gamma(h).symbol(),
                          eoccB->get(i), -eoccB->get(i)*pc_hartree2ev, ps_occB->get(i));

    }
    outfile->Printf("\t------------------------------------------------------------------- \n");

 }// end if

 else if (print_ >= 2) {
    outfile->Printf( "\tState    Symmetry   -IP (a.u.)       IP (eV)        Pole Strength \n");
    outfile->Printf("\t------------------------------------------------------------------- \n");

    for (int i = 0; i < nmo_; ++i){
         int h = irrep_B->get(i);
	 outfile->Printf("\t%3d %10s %15.6f %15.6f %15.6f \n", i+1, ct.gamma(h).symbol(),
                          evals_B->get(i), -evals_B->get(i)*pc_hartree2ev, ps_vec2B->get(i));

    }
    outfile->Printf("\t------------------------------------------------------------------- \n");

 }// end else if

       GFock_primeB.reset();
       g1HalfB.reset();
       UvecB.reset();
       Uvec_primeB.reset();
       PSB.reset();
       gc_transB.reset();
       tempB.reset();
       Diag_g1B.reset();
       ps_vecB.reset();
       eorbB.reset();

       delete irrep_B;
       delete ps_vec2B;
       delete evals_B;
       delete irrep_occB;
       delete ps_occB;
       delete eoccB;

}// if (reference_ == "UNRESTRICTED")

     // For Non-OO methods
     if (orb_opt_ == "FALSE" && reference_ == "RESTRICTED") GFock->scale(2.0);
     else if (orb_opt_ == "FALSE" && reference_ == "UNRESTRICTED") {
              GFockA->scale(2.0);
              GFockB->scale(2.0);
     }

       GFock_primeA.reset();
       g1HalfA.reset();
       UvecA.reset();
       Uvec_primeA.reset();
       PSA.reset();
       gc_transA.reset();
       tempA.reset();
       Diag_g1A.reset();
       ps_vecA.reset();
       eorbA.reset();

       delete irrep_A;
       delete ps_vec2A;
       delete evals_A;
       delete irrep_occA;
       delete ps_occA;
       delete eoccA;

//outfile->Printf("\n ekt_ip is done. \n");

} // end ekt_ip
}} // End Namespaces
