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

// Semicanonicalizing RHF Fock matrix by diagonalizing active-occupied (AOCC-AOCC) and active-virtual (AVIR-AVIR) blocks
#include "occwave.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/pointgrp.h"

using namespace std;


namespace psi{ namespace occwave{

void OCCWave::semi_canonic()
{
        // tell other functions that orbitals are already semi canonical.
        orbs_already_sc = 1;

	SharedMatrix UooA = std::shared_ptr<Matrix>(new Matrix(nirrep_, occpiA, occpiA));
	SharedMatrix UvvA = std::shared_ptr<Matrix>(new Matrix(nirrep_, virtpiA, virtpiA));
	SharedMatrix FockooA = std::shared_ptr<Matrix>(new Matrix(nirrep_, occpiA, occpiA));
	SharedMatrix FockvvA = std::shared_ptr<Matrix>(new Matrix(nirrep_, virtpiA, virtpiA));
	SharedVector eigooA = std::shared_ptr<Vector>(new Vector(nirrep_, occpiA));
	SharedVector eigvvA = std::shared_ptr<Vector>(new Vector(nirrep_, virtpiA));

	UooA->zero();
	UvvA->zero();
	FockooA->zero();
	FockvvA->zero();

	// OCC-OCC
	for(int h = 0; h < nirrep_; h++){
	  for(int i = 0; i < occpiA[h]; i++){
	    eigooA->set(h,i,0.0);
	  }
	}

	// VIR-VIR
	for(int h = 0; h < nirrep_; h++){
	  for(int i = 0; i < virtpiA[h]; i++){
	    eigvvA->set(h,i,0.0);
	  }
	}

       // Fockoo alpha spin case
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
                FockooA->set(h, i, j, FockA->get(h, i, j));
            }
	  }
	}

	// Fockvv alpha spin case
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiA[h]; ++a){
            for(int b = 0 ; b < virtpiA[h]; ++b){
                int aa = a + occpiA[h];
                int bb = b + occpiA[h];
                FockvvA->set(h, a, b, FockA->get(h, aa, bb));
            }
	  }
	}

	// Diagonalize Fock
	FockooA->diagonalize(UooA, eigooA);
	FockvvA->diagonalize(UvvA, eigvvA);

        // Print orbital energies
	if (occ_orb_energy == "TRUE" && mo_optimized == 1) {
	  outfile->Printf("\n\n\t  OMP2 Alpha Orbital Energies (a.u.) \n");
	  outfile->Printf("\t  ---------------------------------- \n");


	  Molecule& mol = *reference_wavefunction_->molecule().get();
	  CharacterTable ct = mol.point_group()->char_table();
          string pgroup = mol.point_group()->symbol();

	  // print occ orb energy
	  outfile->Printf( "\t  Alpha occupied orbitals\n");
	  for (int h=0; h<nirrep_; h++){
	    int count=1;
	    for (int i = 0; i < occpiA[h]; i++){
	      outfile->Printf("\t %2d%-3s %20.10f \n",count,ct.gamma(h).symbol(),eigooA->get(h,i));

	      count++;
	    }// end loop over occpi
	  }// end loop over h


	  // print vir orb energy
	  outfile->Printf( "\n\t  Alpha virtual orbitals\n");
	  for (int h=0; h<nirrep_; h++){
	    int count=1;
	    for (int i = 0; i < virtpiA[h]; i++){
	      outfile->Printf("\t %2d%-3s %20.10f \n",count,ct.gamma(h).symbol(),eigvvA->get(h,i));

	      count++;
	    }// end loop over occpi
	  }// end loop over h

	}// end main if

        // Build U
	UorbA->zero();

	//set to identity: it is necessary if we have frozen core or frozen virtual orbitals.
	//UorbA->identity();

	// Uoo contribution alpha spin case
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
                UorbA->set(h, i, j, UooA->get(h, i, j));
            }
	  }
	}

	// Uvv contribution alpha spin case
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiA[h]; ++a){
            for(int b = 0 ; b < virtpiA[h]; ++b){
                int aa = a + occpiA[h];
                int bb = b + occpiA[h];
                UorbA->set(h, aa, bb, UvvA->get(h, a, b));
            }
	  }
	}

        // Get new MOs
        Ca_new = std::shared_ptr<Matrix>(new Matrix("New alpha MO coefficients", nirrep_, nsopi_, nmopi_));
	Ca_new->zero();
	Ca_new->gemm(false, false, 1.0, Ca_, UorbA, 0.0);
	Ca_->zero();
	Ca_->copy(Ca_new);
	Ca_new.reset();

	if (print_ > 2) {
	  UorbA->print();
	  Ca_->print();
	}

        UooA.reset();
	UvvA.reset();
	FockooA.reset();
	FockvvA.reset();
	eigooA.reset();
	eigvvA.reset();

     // UHF REFERENCE
     if (reference_ == "UNRESTRICTED") {
	SharedMatrix UooB = std::shared_ptr<Matrix>(new Matrix(nirrep_, occpiB, occpiB));
	SharedMatrix UvvB = std::shared_ptr<Matrix>(new Matrix(nirrep_, virtpiB, virtpiB));
	SharedMatrix FockooB = std::shared_ptr<Matrix>(new Matrix(nirrep_, occpiB, occpiB));
	SharedMatrix FockvvB = std::shared_ptr<Matrix>(new Matrix(nirrep_, virtpiB, virtpiB));
	SharedVector eigooB = std::shared_ptr<Vector>(new Vector(nirrep_, occpiB));
	SharedVector eigvvB = std::shared_ptr<Vector>(new Vector(nirrep_, virtpiB));

	UooB->zero();
	UvvB->zero();
	FockooB->zero();
	FockvvB->zero();

	// occ-occ
	for(int h = 0; h < nirrep_; h++){
	  for(int i = 0; i < occpiB[h]; i++){
	    eigooB->set(h,i,0.0);
	  }
	}

	// vir-vir
	for(int h = 0; h < nirrep_; h++){
	  for(int i = 0; i < virtpiB[h]; i++){
	    eigvvB->set(h,i,0.0);
	  }
	}

	// Fockoo beta spin case
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < occpiB[h]; ++j){
                FockooB->set(h, i, j, FockB->get(h, i, j));
            }
	  }
	}

	// Fockvv beta spin case
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiB[h]; ++a){
            for(int b = 0 ; b < virtpiB[h]; ++b){
                int aa = a + occpiB[h];
                int bb = b + occpiB[h];
                FockvvB->set(h, a, b, FockB->get(h, aa, bb));
            }
	  }
	}

	// Diagonalize Fock
	FockooB->diagonalize(UooB, eigooB);
	FockvvB->diagonalize(UvvB, eigvvB);

        // print orbital energies
	if (occ_orb_energy == "TRUE" && mo_optimized == 1 && reference_ == "UNRESTRICTED") {
	  outfile->Printf("\n\n\t  OMP2 Beta Orbital Energies (a.u.) \n");
	  outfile->Printf("\t  --------------------------------- \n");



	  Molecule& mol = *reference_wavefunction_->molecule().get();
	  CharacterTable ct = mol.point_group()->char_table();
          string pgroup = mol.point_group()->symbol();

	  // print occ orb energy
	  outfile->Printf( "\t  Beta occupied orbitals\n");
	  for (int h=0; h<nirrep_; h++){
	    int count=1;
	    for (int i = 0; i < occpiB[h]; i++){
	      outfile->Printf("\t %2d%-3s %20.10f \n",count,ct.gamma(h).symbol(),eigooB->get(h,i));

	      count++;
	    }// end loop over occpi
	  }// end loop over h


	  // print vir orb energy
	  outfile->Printf( "\n\t  Beta virtual orbitals\n");
	  for (int h=0; h<nirrep_; h++){
	    int count=1;
	    for (int i = 0; i < virtpiB[h]; i++){
	      outfile->Printf("\t %2d%-3s %20.10f \n",count,ct.gamma(h).symbol(),eigvvB->get(h,i));

	      count++;
	    }// end loop over occpi
	  }// end loop over h

	}// end main if

        // Build U
	UorbB->zero();
	UorbB->identity();

	// Uoo contribution beta spin case
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < occpiB[h]; ++j){
                UorbB->set(h, i, j, UooB->get(h, i, j));
            }
	  }
	}

	// Uvv contribution beta spin case
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiB[h]; ++a){
            for(int b = 0 ; b < virtpiB[h]; ++b){
                int aa = a + occpiB[h];
                int bb = b + occpiB[h];
                UorbB->set(h, aa, bb, UvvB->get(h, a, b));
            }
	  }
	}

        // Get new MOs
	Cb_new = std::shared_ptr<Matrix>(new Matrix("New beta MO coefficients", nirrep_, nsopi_, nmopi_));
	Cb_new->zero();
	Cb_new->gemm(false, false, 1.0, Cb_, UorbB, 0.0);
	Cb_->zero();
	Cb_->copy(Cb_new);
	Cb_new.reset();


	if (print_ > 2) {
          UorbB->print();
	  Cb_->print();
	}

	UooB.reset();
	UvvB.reset();
	FockooB.reset();
	FockvvB.reset();
	eigooB.reset();
	eigvvB.reset();
     }// end uhf
}
}} // End Namespaces
