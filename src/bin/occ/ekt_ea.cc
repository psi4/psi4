/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <libqt/qt.h>
#include <libtrans/integraltransform.h>

#include "physconst.h"
#include "occwave.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace occwave{
  
void OCCWave::ekt_ea()
{   

//outfile->Printf("\n ekt_ea is starting... \n"); 

     SharedMatrix GFock_copyA = boost::shared_ptr<Matrix>(new Matrix("Alpha GF copy", nirrep_, nmopi_, nmopi_));
     SharedMatrix g1symm_copyA = boost::shared_ptr<Matrix>(new Matrix("Alpha OPDM copy", nirrep_, nmopi_, nmopi_));
     SharedMatrix GFock_copyB = boost::shared_ptr<Matrix>(new Matrix("Beta GF copy", nirrep_, nmopi_, nmopi_));
     SharedMatrix g1symm_copyB = boost::shared_ptr<Matrix>(new Matrix("Alpha OPDM copy", nirrep_, nmopi_, nmopi_));
     G1tilde = boost::shared_ptr<Matrix>(new Matrix("MO-basis 1-G1", nirrep_, nmopi_, nmopi_));
     G1tildeA = boost::shared_ptr<Matrix>(new Matrix("MO-basis Alpha 1-G1", nirrep_, nmopi_, nmopi_));
     G1tildeB = boost::shared_ptr<Matrix>(new Matrix("MO-basis Beta 1-G1", nirrep_, nmopi_, nmopi_));

     // Make sure  GFM is symmetric
     if (ekt_ip_ == "FALSE") {
      if (sym_gfm_ == "TRUE") {
       if (reference_ == "RESTRICTED") {
         SharedMatrix temp(GFock->transpose());
         GFock_copyA->copy(GFock);
         GFock_copyA->add(temp);
         GFock_copyA->scale(0.5);
         GFock->copy(GFock_copyA);

         // Symm OPDM
         SharedMatrix temp2(g1symm->transpose());
         g1symm_copyA->copy(g1symm);
         g1symm_copyA->add(temp2);
         g1symm_copyA->scale(0.5);
         g1symm->copy(g1symm_copyA);
       } 

       else if (reference_ == "UNRESTRICTED") {
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

         // Symm OPDM alpha
         SharedMatrix temp2B(g1symmB->transpose());
         g1symm_copyB->copy(g1symmB);
         g1symm_copyB->add(temp2B);
         g1symm_copyB->scale(0.5);
         g1symmB->copy(g1symm_copyB);
       } 
      } 
     }

     // For Non-OO methods
     if (orb_opt_ == "FALSE" && reference_ == "RESTRICTED") GFock->scale(0.5);  
     else if (orb_opt_ == "FALSE" && reference_ == "UNRESTRICTED") {
              GFockA->scale(0.5);  
              GFockB->scale(0.5);  
     }

     // Form virtual space F matrix
     gfock_ea();

     // Form 1-gamma
     if (reference_ == "RESTRICTED") {
         G1tilde->identity();
         G1tilde->scale(2.0);// Note that it is necessary since G1 = G1A + G1B
         G1tilde->subtract(g1symm);
     }
 
     else if (reference_ == "UNRESTRICTED") {
         G1tildeA->identity();
         G1tildeB->identity();
         G1tildeA->subtract(g1symmA);
         G1tildeB->subtract(g1symmB);
     }

//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
     // Memory allocation
     SharedMatrix GFock_primeA = boost::shared_ptr<Matrix>(new Matrix("Alpha GF prime", nirrep_, nmopi_, nmopi_));
     SharedMatrix g1HalfA = boost::shared_ptr<Matrix>(new Matrix("g^-1/2", nirrep_, nmopi_, nmopi_));
     SharedMatrix UvecA = boost::shared_ptr<Matrix>(new Matrix("UvecA", nirrep_, nmopi_, nmopi_));
     SharedMatrix Uvec_primeA = boost::shared_ptr<Matrix>(new Matrix("Uvec_primeA", nirrep_, nmopi_, nmopi_));
     SharedMatrix PSA = boost::shared_ptr<Matrix>(new Matrix("Alpha pole strength matrix", nirrep_, nmopi_, nmopi_));
     SharedMatrix gc_transA = boost::shared_ptr<Matrix>(new Matrix("Alpha C'*g", nirrep_, nmopi_, nmopi_));
     SharedMatrix tempA = boost::shared_ptr<Matrix>(new Matrix("Alpha temp", nirrep_, nmopi_, nmopi_));
     SharedVector Diag_g1A = boost::shared_ptr<Vector>(new Vector("Diag OO-block OPDM", nirrep_, nmopi_));
     SharedVector ps_vecA = boost::shared_ptr<Vector>(new Vector("alpha pole strength vector", nirrep_, nmopi_));
     SharedVector eorbA = boost::shared_ptr<Vector>(new Vector("eorbA", nirrep_, nmopi_));
      
     // Diagonalize OPDM
     UvecA->zero();
     Diag_g1A->zero();
     if (reference_ == "RESTRICTED") G1tilde->diagonalize(UvecA, Diag_g1A);
     else if (reference_ == "UNRESTRICTED") G1tildeA->diagonalize(UvecA, Diag_g1A);


     // Make sure all eigenvalues are positive
     for (int h = 0; h < nirrep_; ++h) {
          for (int i = 0; i < nmopi_[h]; ++i) {
               if (Diag_g1A->get(h, i) < 0.0) Diag_g1A->set(h, i, -1.0*Diag_g1A->get(h, i));
          }
     }


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
     if (reference_ == "RESTRICTED") tempA->gemm(true, false, 1.0, g1HalfA, Ftilde, 0.0);
     else if (reference_ == "UNRESTRICTED") tempA->gemm(true, false, 1.0, g1HalfA, FtildeA, 0.0);
     GFock_primeA->gemm(false, false, 1.0, tempA, g1HalfA, 0.0);

     // Diagonalize GFock to get orbital energies
     eorbA->zero();
     Uvec_primeA->zero();
     GFock_primeA->diagonalize(Uvec_primeA, eorbA);
     UvecA->gemm(false, false, 1.0, g1HalfA, Uvec_primeA, 0.0);
     
     // Pole strength
     PSA->zero();
     gc_transA->zero();
     if (reference_ == "RESTRICTED") tempA->gemm(false, false, 1.0, G1tilde, UvecA, 0.0); 
     else if (reference_ == "UNRESTRICTED") tempA->gemm(false, false, 1.0, G1tildeA, UvecA, 0.0); 
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

    // Re-Sort virtual orbitals to energy order
    Array1d *evirA = new Array1d("Alpha virtual orb", nvoA);
    Array1d *ps_virA = new Array1d("occupied Pole strength", nvoA);
    Array1i *irrep_virA = new Array1i("virtual IrrepA", nvoA);
    evirA->zero();
    ps_virA->zero();
    irrep_virA->zero();
 

    // Copy 
    for (int i = 0; i < nvoA; ++i) {
         evirA->set(i, evals_A->get(i + nooA));
         ps_virA->set(i, ps_vec2A->get(i + nooA));
         irrep_virA->set(i, irrep_A->get(i + nooA));
    }

    // Sort to ascending order
    for (int i = 0; i < nvoA; ++i) {
         for(int j = nvoA-1; j > i; --j) {
             if (evirA->get(j-1) > evirA->get(j)) {
                 double dum = evirA->get(j-1);
                 evirA->set(j-1, evirA->get(j));
                 evirA->set(j, dum);

                 int dum2 = irrep_virA->get(j-1);
                 irrep_virA->set(j-1, irrep_virA->get(j));
                 irrep_virA->set(j, dum2);

                 double dum3 = ps_virA->get(j-1);
                 ps_virA->set(j-1, ps_virA->get(j));
                 ps_virA->set(j, dum3);
             }
         }
    }
 
    // Print EAs
    outfile->Printf("\n\tEKT-OCC Electron Affinities (Alpha Spin Case) \n"); 
    outfile->Printf("\t------------------------------------------------------------------- \n"); 
    
	  
    Molecule& mol = *reference_wavefunction_->molecule().get();
    CharacterTable ct = mol.point_group()->char_table();
    string pgroup = mol.point_group()->symbol();

 // print alpha EAs
 if (print_ < 2) {
    outfile->Printf( "\tState    Symmetry   -EA (a.u.)       EA (eV)        Pole Strength \n");
    outfile->Printf("\t------------------------------------------------------------------- \n"); 
       
    for (int i = 0; i < nvoA; ++i){
         int h = irrep_virA->get(i);
	 outfile->Printf("\t%3d %10s %15.6f %15.6f %15.6f \n", i+1, ct.gamma(h).symbol(), 
                          evirA->get(i), -evirA->get(i)*pc_hartree2ev, ps_virA->get(i));
	    
    }
    outfile->Printf("\t------------------------------------------------------------------- \n"); 
       
 }// end if

 else if (print_ >= 2) {
    outfile->Printf( "\tState    Symmetry   -EA (a.u.)       EA (eV)        Pole Strength \n");
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
     SharedMatrix GFock_primeB = boost::shared_ptr<Matrix>(new Matrix("Beta OO-block GF prime", nirrep_, nmopi_, nmopi_));
     SharedMatrix g1HalfB = boost::shared_ptr<Matrix>(new Matrix("g^-1/2", nirrep_, nmopi_, nmopi_));
     SharedMatrix UvecB = boost::shared_ptr<Matrix>(new Matrix("UvecB", nirrep_, nmopi_, nmopi_));
     SharedMatrix Uvec_primeB = boost::shared_ptr<Matrix>(new Matrix("Uvec_primeB", nirrep_, nmopi_, nmopi_));
     SharedMatrix PSB = boost::shared_ptr<Matrix>(new Matrix("Beta pole strength matrix", nirrep_, nmopi_, nmopi_));
     SharedMatrix gc_transB = boost::shared_ptr<Matrix>(new Matrix("Beta C'*g", nirrep_, nmopi_, nmopi_));
     SharedMatrix tempB = boost::shared_ptr<Matrix>(new Matrix("Beta temp", nirrep_, nmopi_, nmopi_));
     SharedVector Diag_g1B = boost::shared_ptr<Vector>(new Vector("DiagA OO-block OPDM", nirrep_, nmopi_));
     SharedVector ps_vecB = boost::shared_ptr<Vector>(new Vector("Beta pole strength vector", nirrep_, nmopi_));
     SharedVector eorbB = boost::shared_ptr<Vector>(new Vector("eorbB", nirrep_, nmopi_));

     // Diagonalize OPDM
     UvecB->zero();
     Diag_g1B->zero();
     G1tildeB->diagonalize(UvecB, Diag_g1B);

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
     tempB->gemm(true, false, 1.0, g1HalfB, FtildeB, 0.0);
     GFock_primeB->gemm(false, false, 1.0, tempB, g1HalfB, 0.0);

     // Diagonalize GFock to get orbital energies
     eorbB->zero();
     Uvec_primeB->zero();
     GFock_primeB->diagonalize(Uvec_primeB, eorbB);
     UvecB->gemm(false, false, 1.0, g1HalfB, Uvec_primeB, 0.0);
     
     // Pole strength
     PSB->zero();
     gc_transB->zero();
     tempB->gemm(false, false, 1.0, G1tildeB, UvecB, 0.0); 
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

    // Re-Sort virtual orbitals to energy order
    Array1d *evirB = new Array1d("Alpha virtual orb", nvoB);
    Array1d *ps_virB = new Array1d("occupied Pole strength", nvoB);
    Array1i *irrep_virB = new Array1i("virtual IrrepA", nvoB);
    evirB->zero();
    ps_virB->zero();
    irrep_virB->zero();
 

    // Copy 
    for (int i = 0; i < nvoB; ++i) {
         evirB->set(i, evals_B->get(i + nooB));
         ps_virB->set(i, ps_vec2B->get(i + nooB));
         irrep_virB->set(i, irrep_B->get(i + nooB));
    }

    // Sort to ascending order
    for (int i = 0; i < nvoB; ++i) {
         for(int j = nvoB-1; j > i; --j) {
             if (evirB->get(j-1) > evirB->get(j)) {
                 double dum = evirB->get(j-1);
                 evirB->set(j-1, evirB->get(j));
                 evirB->set(j, dum);

                 int dum2 = irrep_virB->get(j-1);
                 irrep_virB->set(j-1, irrep_virB->get(j));
                 irrep_virB->set(j, dum2);

                 double dum3 = ps_virB->get(j-1);
                 ps_virB->set(j-1, ps_virB->get(j));
                 ps_virB->set(j, dum3);
             }
         }
    }

    // Print EAs
    outfile->Printf("\n\tEKT-OCC Electron Affinities (Beta Spin Case) \n"); 
    outfile->Printf("\t------------------------------------------------------------------- \n"); 
    

// print alpha EAs
 if (print_ < 2) {
    outfile->Printf( "\tState    Symmetry   -EA (a.u.)       EA (eV)        Pole Strength \n");
    outfile->Printf("\t------------------------------------------------------------------- \n"); 
       
    for (int i = 0; i < nvoB; ++i){
         int h = irrep_virB->get(i);
	 outfile->Printf("\t%3d %10s %15.6f %15.6f %15.6f \n", i+1, ct.gamma(h).symbol(), 
                          evirB->get(i), -evirB->get(i)*pc_hartree2ev, ps_virB->get(i));
	    
    }
    outfile->Printf("\t------------------------------------------------------------------- \n"); 
       
 }// end if

 else if (print_ >= 2) {
    outfile->Printf( "\tState    Symmetry   -EA (a.u.)       EA (eV)        Pole Strength \n");
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
       delete irrep_virB;
       delete ps_virB;
       delete evirB;

}// if (reference_ == "UNRESTRICTED") 

     // For Non-OO methods
     /*
     if (orb_opt_ == "FALSE" && reference_ == "RESTRICTED" && dertype == "FIRST") GFock->scale(2.0);  
     else if (orb_opt_ == "FALSE" && reference_ == "UNRESTRICTED" && dertype == "FIRST") {
              GFockA->scale(2.0);  
              GFockB->scale(2.0);  
     }
     */

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
       delete irrep_virA;
       delete ps_virA;
       delete evirA;

//outfile->Printf("\n ekt_ea is done. \n"); 

} // end ekt_ip
}} // End Namespaces
