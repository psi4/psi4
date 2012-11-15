// Semicanonicalizing RHF Fock matrix by diagonalizing active-occupied (AOCC-AOCC) and active-virtual (AVIR-AVIR) blocks

/** Standard library includes */
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip> 
#include <vector> 

/** Required PSI4 includes */ 
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>


/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>

#include "ocepawave.h"
#include "defines.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace ocepawave{
 
void OCEPAWave::semi_canonic()
{

	SharedMatrix UooA = boost::shared_ptr<Matrix>(new Matrix(nirrep_, aoccpiA, aoccpiA));
	SharedMatrix UvvA = boost::shared_ptr<Matrix>(new Matrix(nirrep_, avirtpiA, avirtpiA));
	SharedMatrix FockooA = boost::shared_ptr<Matrix>(new Matrix(nirrep_, aoccpiA, aoccpiA));
	SharedMatrix FockvvA = boost::shared_ptr<Matrix>(new Matrix(nirrep_, avirtpiA, avirtpiA));
	SharedVector eigooA = boost::shared_ptr<Vector>(new Vector(nirrep_, aoccpiA));
	SharedVector eigvvA = boost::shared_ptr<Vector>(new Vector(nirrep_, avirtpiA));

	UooA->zero();
	UvvA->zero();
	FockooA->zero();
	FockvvA->zero();
     	
	// OCC-OCC 
	for(int h = 0; h < nirrep_; h++){
	  for(int i = 0; i < aoccpiA[h]; i++){
	    eigooA->set(h,i,0.0);
	  }
	}
	
	// VIR-VIR
	for(int h = 0; h < nirrep_; h++){
	  for(int i = 0; i < avirtpiA[h]; i++){
	    eigvvA->set(h,i,0.0);
	  }
	}

       // Fockoo alpha spin case
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
                FockooA->set(h, i, j, FockA->get(h, i, j));
            }
	  }
	}
	
	// Fockvv alpha spin case
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < avirtpiA[h]; ++a){
            for(int b = 0 ; b < avirtpiA[h]; ++b){
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
	if (ocepa_orb_energy == "TRUE" && mo_optimized == 1) {
	  fprintf(outfile,"\n\n\t  OCEPA Alpha Orbital Energies (a.u.) \n"); 
	  fprintf(outfile,"\t  ---------------------------------- \n"); 
	  fflush(outfile);
	  
	  Molecule& mol = *reference_wavefunction_->molecule().get();
	  CharacterTable ct = mol.point_group()->char_table();
          string pgroup = mol.point_group()->symbol();
	  
	  // print occ orb energy
	  fprintf(outfile, "\t  Alpha occupied orbitals\n");
	  for (int h=0; h<nirrep_; h++){
	    int count=1;
	    for (int i = 0; i < aoccpiA[h]; i++){
	      fprintf(outfile,"\t %2d%-3s %20.10f \n",count,ct.gamma(h).symbol(),eigooA->get(h,i));
	      fflush(outfile);   
	      count++;
	    }// end loop over aoccpi
	  }// end loop over h
	  
	  
	  // print vir orb energy
	  fprintf(outfile, "\n\t  Alpha virtual orbitals\n");
	  for (int h=0; h<nirrep_; h++){
	    int count=1;
	    for (int i = 0; i < avirtpiA[h]; i++){
	      fprintf(outfile,"\t %2d%-3s %20.10f \n",count,ct.gamma(h).symbol(),eigvvA->get(h,i));
	      fflush(outfile);
	      count++;
	    }// end loop over aoccpi
	  }// end loop over h
	  
	}// end main if

        // Build U	
	UorbA->zero();
	
	//set to identity
	UorbA->identity();
	
	// Uoo contribution alpha spin case
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
                UorbA->set(h, i, j, UooA->get(h, i, j));
            }
	  }
	}
	
	// Uvv contribution alpha spin case
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < avirtpiA[h]; ++a){
            for(int b = 0 ; b < avirtpiA[h]; ++b){
                int aa = a + occpiA[h];
                int bb = b + occpiA[h];
                UorbA->set(h, aa, bb, UvvA->get(h, a, b));
            }
	  }
	}

        // Get new MOs
        Ca_new = boost::shared_ptr<Matrix>(new Matrix("New alpha MO coefficients", nirrep_, nsopi_, nmopi_));
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
	SharedMatrix UooB = boost::shared_ptr<Matrix>(new Matrix(nirrep_, aoccpiB, aoccpiB));
	SharedMatrix UvvB = boost::shared_ptr<Matrix>(new Matrix(nirrep_, avirtpiB, avirtpiB));
	SharedMatrix FockooB = boost::shared_ptr<Matrix>(new Matrix(nirrep_, aoccpiB, aoccpiB));
	SharedMatrix FockvvB = boost::shared_ptr<Matrix>(new Matrix(nirrep_, avirtpiB, avirtpiB));
	SharedVector eigooB = boost::shared_ptr<Vector>(new Vector(nirrep_, aoccpiB));
	SharedVector eigvvB = boost::shared_ptr<Vector>(new Vector(nirrep_, avirtpiB));

	UooB->zero();
	UvvB->zero();
	FockooB->zero();
	FockvvB->zero();

	// occ-occ 
	for(int h = 0; h < nirrep_; h++){
	  for(int i = 0; i < aoccpiB[h]; i++){
	    eigooB->set(h,i,0.0);
	  }
	}
	
	// vir-vir
	for(int h = 0; h < nirrep_; h++){
	  for(int i = 0; i < avirtpiB[h]; i++){
	    eigvvB->set(h,i,0.0);
	  }
	}

	// Fockoo beta spin case
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int j = 0 ; j < aoccpiB[h]; ++j){
                FockooB->set(h, i, j, FockB->get(h, i, j));
            }
	  }
	}

	// Fockvv beta spin case
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < avirtpiB[h]; ++a){
            for(int b = 0 ; b < avirtpiB[h]; ++b){
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
	if (ocepa_orb_energy == "TRUE" && mo_optimized == 1 && reference_ == "UNRESTRICTED") {
	  fprintf(outfile,"\n\n\t  OCEPA Beta Orbital Energies (a.u.) \n"); 
	  fprintf(outfile,"\t  --------------------------------- \n"); 
	  fflush(outfile);
	  
	  
	  Molecule& mol = *reference_wavefunction_->molecule().get();
	  CharacterTable ct = mol.point_group()->char_table();
          string pgroup = mol.point_group()->symbol();
	  
	  // print occ orb energy
	  fprintf(outfile, "\t  Beta occupied orbitals\n");
	  for (int h=0; h<nirrep_; h++){
	    int count=1;
	    for (int i = 0; i < aoccpiB[h]; i++){
	      fprintf(outfile,"\t %2d%-3s %20.10f \n",count,ct.gamma(h).symbol(),eigooB->get(h,i));
	      fflush(outfile);
	      count++;
	    }// end loop over aoccpi
	  }// end loop over h
	  
	  
	  // print vir orb energy
	  fprintf(outfile, "\n\t  Beta virtual orbitals\n");
	  for (int h=0; h<nirrep_; h++){
	    int count=1;
	    for (int i = 0; i < avirtpiB[h]; i++){
	      fprintf(outfile,"\t %2d%-3s %20.10f \n",count,ct.gamma(h).symbol(),eigvvB->get(h,i));
	      fflush(outfile);
	      count++;
	    }// end loop over aoccpi
	  }// end loop over h
	  
	}// end main if
     	
        // Build U
	UorbB->zero();
	UorbB->identity();
	
	// Uoo contribution beta spin case
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int j = 0 ; j < aoccpiB[h]; ++j){
                UorbB->set(h, i, j, UooB->get(h, i, j));
            }
	  }
	}

	// Uvv contribution beta spin case
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < avirtpiB[h]; ++a){
            for(int b = 0 ; b < avirtpiB[h]; ++b){
                int aa = a + occpiB[h];
                int bb = b + occpiB[h];
                UorbB->set(h, aa, bb, UvvB->get(h, a, b));
            }
	  }
	}

        // Get new MOs
	Cb_new = boost::shared_ptr<Matrix>(new Matrix("New beta MO coefficients", nirrep_, nsopi_, nmopi_));
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

        

/********************************************************************************************/
/************************** Save MO coefficients to Chkpt file ******************************/
/********************************************************************************************/
        /*
	C_pitzerA = Ca_->to_block_matrix();    
	C_pitzerB = Cb_->to_block_matrix();    
	chkpt_->wt_alpha_scf(C_pitzerA);
	chkpt_->wt_beta_scf(C_pitzerB);
	free_block(C_pitzerA);
	free_block(C_pitzerB);
	*/

}
}} // End Namespaces

