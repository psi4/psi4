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

#include "omp2wave.h"
#include "defines.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace omp2wave{
 
void OMP2Wave::semi_canonic()
{

/********************************************************************************************/
/************************** memalloc ********************************************************/
/********************************************************************************************/ 
	SharedMatrix UooA = shared_ptr<Matrix>(new Matrix(nirreps, aoccpiA, aoccpiA));
	SharedMatrix UooB = shared_ptr<Matrix>(new Matrix(nirreps, aoccpiB, aoccpiB));
	SharedMatrix UvvA = shared_ptr<Matrix>(new Matrix(nirreps, avirtpiA, avirtpiA));
	SharedMatrix UvvB = shared_ptr<Matrix>(new Matrix(nirreps, avirtpiB, avirtpiB));
	SharedMatrix FockooA = shared_ptr<Matrix>(new Matrix(nirreps, aoccpiA, aoccpiA));
	SharedMatrix FockooB = shared_ptr<Matrix>(new Matrix(nirreps, aoccpiB, aoccpiB));
	SharedMatrix FockvvA = shared_ptr<Matrix>(new Matrix(nirreps, avirtpiA, avirtpiA));
	SharedMatrix FockvvB = shared_ptr<Matrix>(new Matrix(nirreps, avirtpiB, avirtpiB));
	SharedVector eigooA = shared_ptr<Vector>(new Vector(nirreps, aoccpiA));
	SharedVector eigooB = shared_ptr<Vector>(new Vector(nirreps, aoccpiB));
	SharedVector eigvvA = shared_ptr<Vector>(new Vector(nirreps, avirtpiA));
	SharedVector eigvvB = shared_ptr<Vector>(new Vector(nirreps, avirtpiB));
      
/********************************************************************************************/
/************************** Initialize ******************************************************/
/********************************************************************************************/
	UooA->zero();
	UooB->zero();
	UvvA->zero();
	UvvB->zero();
	FockooA->zero();
	FockooB->zero();
	FockvvA->zero();
	FockvvB->zero();
	
	// OCC-OCC 
	for(int h = 0; h < nirreps; h++){
	  for(int i = 0; i < aoccpiA[h]; i++){
	    eigooA->set(h,i,0.0);
	  }
	}
	
	// occ-occ 
	for(int h = 0; h < nirreps; h++){
	  for(int i = 0; i < aoccpiB[h]; i++){
	    eigooB->set(h,i,0.0);
	  }
	}
	
	// VIR-VIR
	for(int h = 0; h < nirreps; h++){
	  for(int i = 0; i < avirtpiA[h]; i++){
	    eigvvA->set(h,i,0.0);
	  }
	}
	
	// vir-vir
	for(int h = 0; h < nirreps; h++){
	  for(int i = 0; i < avirtpiB[h]; i++){
	    eigvvB->set(h,i,0.0);
	  }
	}
	
/********************************************************************************************/
/************************** Fockoo & Fockvv *************************************************/
/********************************************************************************************/
	// Fockoo alpha spin case
	for(int i = nfrzc; i < nooA; i++){
	  for(int j = nfrzc; j < nooA; j++){
	    int i2 = c1topitzerA[i];
	    int j2 = c1topitzerA[j];
	    int i3 = pitzer2symblk[i2];
	    int j3 = pitzer2symblk[j2];
	    int hi=mosym[i2];
	    int hj=mosym[j2];
	    
	    if (hi == hj) {
	      int i4 = i3 - frzcpi[hi];
	      int j4 = j3 - frzcpi[hj];
	      FockooA->set(hi,i4,j4,FockA->get(hi,i3,j3)); 
	    }     
	  }
	}
	
	// Fockoo beta spin case
	for(int i = nfrzc; i < nooB; i++){
	  for(int j = nfrzc; j < nooB; j++){
	    int i2 = c1topitzerB[i];
	    int j2 = c1topitzerB[j];
	    int i3 = pitzer2symblk[i2];
	    int j3 = pitzer2symblk[j2];
	    int hi=mosym[i2];
	    int hj=mosym[j2];
	    
	    if (hi == hj) {
	      int i4 = i3 - frzcpi[hi];
	      int j4 = j3 - frzcpi[hj];
	      FockooB->set(hi,i4,j4,FockB->get(hi,i3,j3)); 
	    }     
	  }
	}
	
	// Fockvv alpha spin case
	for(int a = nooA; a < npop; a++){
	  for(int b = nooA; b < npop; b++){
	    int a2 = c1topitzerA[a];
	    int b2 = c1topitzerA[b];
	    int a3 = pitzer2symblk[a2];
	    int b3 = pitzer2symblk[b2];
	    int ha=mosym[a2];
	    int hb=mosym[b2];
	    
	    if (ha == hb) {
	      int a4 = a3 - occpiA[ha];
	      int b4 = b3 - occpiA[hb];
	      FockvvA->set(ha,a4,b4,FockA->get(ha,a3,b3)); 
	    }
	  }
	}
	
	// Fockvv beta spin case
	for(int a = nooB; a < npop; a++){
	  for(int b = nooB; b < npop; b++){
	    int a2 = c1topitzerB[a];
	    int b2 = c1topitzerB[b];
	    int a3 = pitzer2symblk[a2];
	    int b3 = pitzer2symblk[b2];
	    int ha=mosym[a2];
	    int hb=mosym[b2];
	    
	    if (ha == hb) {
	      int a4 = a3 - occpiB[ha];
	      int b4 = b3 - occpiB[hb];
	      FockvvB->set(ha,a4,b4,FockB->get(ha,a3,b3)); 
	    }
	  }
	}
	
/********************************************************************************************/
/************************** Diagonalize Fockoo & Fockvv *************************************/
/********************************************************************************************/
	// Diagonalize Fockoo  
	FockooA->diagonalize(UooA, eigooA);
	FockooB->diagonalize(UooB, eigooB);
	
	// Diagonalize Fockvv  
	FockvvA->diagonalize(UvvA, eigvvA);
	FockvvB->diagonalize(UvvB, eigvvB);

/********************************************************************************************/
/************************** OMP2 Alpha Orbital Energies *************************************/
/********************************************************************************************/
	if (omp2_orb_energy == "TRUE" && mo_optimized == 1) {
	  fprintf(outfile,"\n\n\t  OMP2 Alpha Orbital Energies (a.u.) \n"); 
	  fprintf(outfile,"\t  ---------------------------------- \n"); 
	  fflush(outfile);
	  
	  Molecule& mol = *reference_wavefunction_->molecule().get();
	  CharacterTable ct = mol.point_group()->char_table();
          string pgroup = mol.point_group()->symbol();
	  
	  // print occ orb energy
	  fprintf(outfile, "\t  Alpha occupied orbitals\n");
	  for (int h=0; h<nirreps; h++){
	    int count=1;
	    for (int i = 0; i < aoccpiA[h]; i++){
	      fprintf(outfile,"\t %2d%-3s %20.10f \n",count,ct.gamma(h).symbol(),eigooA->get(h,i));
	      fflush(outfile);   
	      count++;
	    }// end loop over aoccpi
	  }// end loop over h
	  
	  
	  // print vir orb energy
	  fprintf(outfile, "\n\t  Alpha virtual orbitals\n");
	  for (int h=0; h<nirreps; h++){
	    int count=1;
	    for (int i = 0; i < avirtpiA[h]; i++){
	      fprintf(outfile,"\t %2d%-3s %20.10f \n",count,ct.gamma(h).symbol(),eigvvA->get(h,i));
	      fflush(outfile);
	      count++;
	    }// end loop over aoccpi
	  }// end loop over h
	  
	}// end main if
	
	
/********************************************************************************************/
/************************** OMP2 Beta Orbital Energies **************************************/
/********************************************************************************************/
	if (omp2_orb_energy == "TRUE" && mo_optimized == 1) {
	  fprintf(outfile,"\n\n\t  OMP2 Beta Orbital Energies (a.u.) \n"); 
	  fprintf(outfile,"\t  --------------------------------- \n"); 
	  fflush(outfile);
	  
	  
	  Molecule& mol = *reference_wavefunction_->molecule().get();
	  CharacterTable ct = mol.point_group()->char_table();
          string pgroup = mol.point_group()->symbol();
	  
	  // print occ orb energy
	  fprintf(outfile, "\t  Beta occupied orbitals\n");
	  for (int h=0; h<nirreps; h++){
	    int count=1;
	    for (int i = 0; i < aoccpiB[h]; i++){
	      fprintf(outfile,"\t %2d%-3s %20.10f \n",count,ct.gamma(h).symbol(),eigooB->get(h,i));
	      fflush(outfile);
	      count++;
	    }// end loop over aoccpi
	  }// end loop over h
	  
	  
	  // print vir orb energy
	  fprintf(outfile, "\n\t  Beta virtual orbitals\n");
	  for (int h=0; h<nirreps; h++){
	    int count=1;
	    for (int i = 0; i < avirtpiB[h]; i++){
	      fprintf(outfile,"\t %2d%-3s %20.10f \n",count,ct.gamma(h).symbol(),eigvvB->get(h,i));
	      fflush(outfile);
	      count++;
	    }// end loop over aoccpi
	  }// end loop over h
	  
	}// end main if

/********************************************************************************************/
/************************** Uorbrot *********************************************************/
/********************************************************************************************/
	UorbA->zero();
	UorbB->zero();
	
	//set to identity
	UorbA->identity();
	UorbB->identity();
	
	// Fockoo contribution alpha spin case
	for(int i = nfrzc; i < nooA; i++){
	  for(int j = nfrzc; j < nooA; j++){
	    int i2 = c1topitzerA[i];
	    int j2 = c1topitzerA[j];
	    int i3 = pitzer2symblk[i2];
	    int j3 = pitzer2symblk[j2];
	    int hi=mosym[i2];
	    int hj=mosym[j2];
	    
	    if (hi == hj) {
	      int i4 = i3 - frzcpi[hi];
	      int j4 = j3 - frzcpi[hj];
	      UorbA->set(hi,i3,j3,UooA->get(hi,i4,j4)); 
	    }
	  }
	}
	
	// Fockoo contribution beta spin case
	for(int i = nfrzc; i < nooB; i++){
	  for(int j = nfrzc; j < nooB; j++){
	    int i2 = c1topitzerB[i];
	    int j2 = c1topitzerB[j];
	    int i3 = pitzer2symblk[i2];
	    int j3 = pitzer2symblk[j2];
	    int hi=mosym[i2];
	    int hj=mosym[j2];
	    
	    if (hi == hj) {
	      int i4 = i3 - frzcpi[hi];
	      int j4 = j3 - frzcpi[hj];
	      UorbB->set(hi,i3,j3,UooB->get(hi,i4,j4)); 
	    }
	  }
	}
	
	// Fockvv contribution alpha spin case
	for(int a = nooA; a < npop; a++){
	  for(int b = nooA; b < npop; b++){
	    int a2 = c1topitzerA[a];
	    int b2 = c1topitzerA[b];
	    int a3 = pitzer2symblk[a2];
	    int b3 = pitzer2symblk[b2];
	    int ha=mosym[a2];
	    int hb=mosym[b2];
	    
	    if (ha == hb) {
	      int a4 = a3 - occpiA[ha];
	      int b4 = b3 - occpiA[hb];
	      UorbA->set(ha,a3,b3,UvvA->get(ha,a4,b4)); 
	    }
	  }
	}
	
	// Fockvv contribution beta spin case
	for(int a = nooB; a < npop; a++){
	  for(int b = nooB; b < npop; b++){
	    int a2 = c1topitzerB[a];
	    int b2 = c1topitzerB[b];
	    int a3 = pitzer2symblk[a2];
	    int b3 = pitzer2symblk[b2];
	    int ha=mosym[a2];
	    int hb=mosym[b2];
	    
	    if (ha == hb) {
	      int a4 = a3 - occpiB[ha];
	      int b4 = b3 - occpiB[hb];
	      UorbB->set(ha,a3,b3,UvvB->get(ha,a4,b4)); 
	    }
	  }
	}

/********************************************************************************************/
/************************** Build new MO coeff. *********************************************/
/********************************************************************************************/
        Ca_new = shared_ptr<Matrix>(new Matrix("New alpha MO coefficients", nirreps, sopi, mopi));
	Cb_new = shared_ptr<Matrix>(new Matrix("New beta MO coefficients", nirreps, sopi, mopi));
	Ca_new->zero();
	Cb_new->zero();
	Ca_new->gemm(false, false, 1.0, Ca_, UorbA, 0.0); 
	Cb_new->gemm(false, false, 1.0, Cb_, UorbB, 0.0); 
	Ca_->zero();
	Cb_->zero();
	Ca_->copy(Ca_new);
	Cb_->copy(Cb_new);
	Ca_new.reset();
	Cb_new.reset();

	if (print_ > 1) {
	  UorbA->print();
          UorbB->print();
	  Ca_->print();
	  Cb_->print();
	}

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
	
/********************************************************************************************/
/************************** free array ******************************************************/
/********************************************************************************************/
        UooA.reset();
	UooB.reset();
	UvvA.reset();
	UvvB.reset();
	FockooA.reset();
	FockooB.reset();
	FockvvA.reset();
	FockvvB.reset();
	eigooA.reset();
	eigooB.reset();
	eigvvA.reset();
	eigvvB.reset();

/********************************************************************************************/ 
/********************************************************************************************/ 	

}
}} // End Namespaces

