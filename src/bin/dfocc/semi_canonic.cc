/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

// Semicanonicalizing RHF Fock matrix by diagonalizing active-occupied (AOCC-AOCC) and active-virtual (AVIR-AVIR) blocks
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{
 
void DFOCC::semi_canonic()
{
        // tell oter functions tat orbitals are already semi canonical.
        orbs_already_sc = 1;

	SharedTensor2d UooA = boost::shared_ptr<Tensor2d>(new Tensor2d("UooA", naoccA, naoccA));
	SharedTensor2d UvvA = boost::shared_ptr<Tensor2d>(new Tensor2d("UvvA", navirA, navirA));
	SharedTensor2d FockooA = boost::shared_ptr<Tensor2d>(new Tensor2d("Fock <I|J>", naoccA, naoccA));
	SharedTensor2d FockvvA = boost::shared_ptr<Tensor2d>(new Tensor2d("Fock <A|B>", navirA, navirA));
	SharedTensor1d eigooA = boost::shared_ptr<Tensor1d>(new Tensor1d("epsilon <I|J>", naoccA));
	SharedTensor1d eigvvA = boost::shared_ptr<Tensor1d>(new Tensor1d("epsilon <A|B>", navirA));
     	
        // Fockoo alpha spin case
        #pragma omp parallel for
	for(int i = 0 ; i < naoccA; ++i){
            for(int j = 0 ; j < naoccA; ++j){
                FockooA->set(i, j, FockA->get(i + nfrzc, j + nfrzc));
            }
	}
	
	// Fockvv alpha spin case
	#pragma omp parallel for
	for(int a = 0 ; a < navirA; ++a){
            for(int b = 0 ; b < navirA; ++b){
                int aa = a + noccA;
                int bb = b + noccA;
                FockvvA->set(a, b, FockA->get(aa, bb));
            }
	}

	// Diagonalize Fock  
	FockooA->diagonalize(UooA, eigooA, cutoff);
	FockvvA->diagonalize(UvvA, eigvvA, cutoff);

        // Print orbital energies
	if (occ_orb_energy == "TRUE" && mo_optimized == 1) {
	  psi::fprintf(outfile,"\n\n\tOCC Alpha Orbital Energies (a.u.) \n"); 
	  psi::fprintf(outfile,"\t  ---------------------------------- \n"); 
	  fflush(outfile);
	  
	  // print occ orb energy
	 psi::fprintf(outfile, "\tAlpha occupied orbitals\n");
	 for (int i = 0; i < naoccA; i++){
	      psi::fprintf(outfile,"\t%2d %20.10f \n",i,eigooA->get(i));
	      fflush(outfile);   
	 }// end loop over naocc
	  
	  // print vir orb energy
	  psi::fprintf(outfile, "\n\tAlpha virtual orbitals\n");
	  for (int i = 0; i < navirA; i++){
	      psi::fprintf(outfile,"\t%2d %20.10f \n",i + noccA,eigvvA->get(i));
	      fflush(outfile);
	  }// end loop over naocc
	  
	}// end main if

        // Build U	
	UorbA->zero();
	
	//set to identity: it is necessary if we have frozen core or frozen virtual orbitals.
	UorbA->identity();
	
	// Uoo contribution alpha spin case
        #pragma omp parallel for
	for(int i = 0 ; i < naoccA; ++i){
            for(int j = 0 ; j < naoccA; ++j){
                UorbA->set(i + nfrzc, j + nfrzc, UooA->get(i,j));
            }
	}
	
	// Uvv contribution alpha spin case
	#pragma omp parallel for
	for(int a = 0 ; a < navirA; ++a){
            for(int b = 0 ; b < navirA; ++b){
                int aa = a + noccA;
                int bb = b + noccA;
                UorbA->set(aa, bb, UvvA->get(a,b));
            }
	}

        // Get new MOs
        SharedTensor2d Ca_new = boost::shared_ptr<Tensor2d>(new Tensor2d("New alpha MO coefficients", nso_, nmo_));
	Ca_new->gemm(false, false, CmoA, UorbA, 1.0, 0.0); 
	CmoA->copy(Ca_new);
	Ca_new.reset();

	if (print_ > 2) {
	  UorbA->print();
	  CmoA->print();
	}

        UooA.reset();
	UvvA.reset();
	FockooA.reset();
	FockvvA.reset();
	eigooA.reset();
	eigvvA.reset();

        // Build Cocc
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < noccA; i++) {
                 CoccA->set(mu, i, CmoA->get(mu, i));
             }
        }
        if (print_ > 2) CoccA->print();

        // Build Cvir
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < nvirA; a++) {
                 CvirA->set(mu, a, CmoA->get(mu, a + noccA));
             }
        }
        if (print_ > 2) CvirA->print();
 
        // Build active Caocc
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < naoccA; i++) {
                 CaoccA->set(mu, i, CmoA->get(mu, i + nfrzc));
             }
        }
        if (print_ > 2) CaoccA->print();

        // Build active Cvir
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < navirA; a++) {
                 CavirA->set(mu, a, CmoA->get(mu, a + noccA));
             }
        }
        if (print_ > 2) CavirA->print();

//==========================================================================================
//========================= UHF REFERENCE ==================================================
//==========================================================================================
     if (reference_ == "UNRESTRICTED") {
       	SharedTensor2d UooB = boost::shared_ptr<Tensor2d>(new Tensor2d("UooB", naoccB, naoccB));
	SharedTensor2d UvvB = boost::shared_ptr<Tensor2d>(new Tensor2d("UvvB", navirB, navirB));
	SharedTensor2d FockooB = boost::shared_ptr<Tensor2d>(new Tensor2d("Fock <i|j>", naoccB, naoccB));
	SharedTensor2d FockvvB = boost::shared_ptr<Tensor2d>(new Tensor2d("Fock <a|b>", navirB, navirB));
	SharedTensor1d eigooB = boost::shared_ptr<Tensor1d>(new Tensor1d("epsilon <i|j>", naoccB));
	SharedTensor1d eigvvB = boost::shared_ptr<Tensor1d>(new Tensor1d("epsilon <a|b>", navirB));
     	
        // Fockoo beta spin case
        #pragma omp parallel for
	for(int i = 0 ; i < naoccB; ++i){
            for(int j = 0 ; j < naoccB; ++j){
                FockooB->set(i, j, FockB->get(i + nfrzc, j + nfrzc));
            }
	}
	
	// Fockvv beta spin case
	#pragma omp parallel for
	for(int a = 0 ; a < navirB; ++a){
            for(int b = 0 ; b < navirB; ++b){
                int aa = a + noccB;
                int bb = b + noccB;
                FockvvB->set(a, b, FockB->get(aa, bb));
            }
	}

	// Diagonalize Fock  
	FockooB->diagonalize(UooB, eigooB, cutoff);
	FockvvB->diagonalize(UvvB, eigvvB, cutoff);

        // Print orbital energies
	if (occ_orb_energy == "TRUE" && mo_optimized == 1) {
	  psi::fprintf(outfile,"\n\n\tOCC Beta Orbital Energies (a.u.) \n"); 
	  psi::fprintf(outfile,"\t  ---------------------------------- \n"); 
	  fflush(outfile);
	  
	  // print occ orb energy
	 psi::fprintf(outfile, "\tBeta occupied orbitals\n");
	 for (int i = 0; i < naoccB; i++){
	      psi::fprintf(outfile,"\t%2d %20.10f \n",i,eigooB->get(i));
	      fflush(outfile);   
	 }// end loop over naocc
	  
	  // print vir orb energy
	  psi::fprintf(outfile, "\n\tBeta virtual orbitals\n");
	  for (int i = 0; i < navirB; i++){
	      psi::fprintf(outfile,"\t%2d %20.10f \n",i + noccB,eigvvB->get(i));
	      fflush(outfile);
	  }// end loop over naocc
	  
	}// end main if

        // Build U	
	UorbB->zero();
	
	//set to identity: it is necessary if we have frozen core or frozen virtual orbitals.
	UorbB->identity();
	
	// Uoo contribution beta spin case
        #pragma omp parallel for
	for(int i = 0 ; i < naoccB; ++i){
            for(int j = 0 ; j < naoccB; ++j){
                UorbB->set(i + nfrzc, j + nfrzc, UooB->get(i,j));
            }
	}
	
	// Uvv contribution beta spin case
	#pragma omp parallel for
	for(int a = 0 ; a < navirB; ++a){
            for(int b = 0 ; b < navirB; ++b){
                int aa = a + noccB;
                int bb = b + noccB;
                UorbB->set(aa, bb, UvvB->get(a,b));
            }
	}

        // Get new MOs
        SharedTensor2d Cb_new = boost::shared_ptr<Tensor2d>(new Tensor2d("New beta MO coefficients", nso_, nmo_));
	Cb_new->gemm(false, false, CmoB, UorbB, 1.0, 0.0); 
	CmoB->copy(Cb_new);
	Cb_new.reset();

	if (print_ > 2) {
	  UorbB->print();
	  CmoB->print();
	}

        UooB.reset();
	UvvB.reset();
	FockooB.reset();
	FockvvB.reset();
	eigooB.reset();
	eigvvB.reset();

        // Build Cocc
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < noccB; i++) {
                 CoccB->set(mu, i, CmoB->get(mu, i));
             }
        }
        if (print_ > 2) CoccB->print();

        // Build Cvir
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < nvirB; a++) {
                 CvirB->set(mu, a, CmoB->get(mu, a + noccB));
             }
        }
        if (print_ > 2) CvirB->print();
 
        // Build active Caocc
        for (int mu = 0; mu < nso_; mu++) {
             for (int i = 0; i < naoccB; i++) {
                 CaoccB->set(mu, i, CmoB->get(mu, i + nfrzc));
             }
        }
        if (print_ > 2) CaoccB->print();

        // Build active Cvir
        for (int mu = 0; mu < nso_; mu++) {
             for (int a = 0; a < navirB; a++) {
                 CavirB->set(mu, a, CmoB->get(mu, a + noccB));
             }
        }
        if (print_ > 2) CavirB->print();
     }// end uhf	
}
}} // End Namespaces


