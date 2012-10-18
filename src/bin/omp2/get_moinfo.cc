// p_so(pitzer) = p_symblk + PitzerOffset[h]; where h=mosym[p_symblk]
// p_symblk = pitzer2symblk[p_so(pitzer)];
// c1topitzer 
// pitzer2c1 
// c1toqt 
// qt2c1 

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


/** Required PSI3 includes */ 
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
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

void OMP2Wave::get_moinfo()
{      
  //fprintf(outfile,"\n get_moinfo is starting... \n"); fflush(outfile);
//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference == "RHF") {


/********************************************************************************************/
/************************** MO info *********************************************************/
/********************************************************************************************/
	// Read in mo info
	/*
	nso_ = chkpt_->rd_nso_();
	nmo_ = chkpt_->rd_nmo_();
	nao = chkpt_->rd_nao();
	nfrzc = chkpt_->rd_nfzc();
	nfrzv = chkpt_->rd_nfzv();
	nirrep_ = chkpt_->rd_nirrep_();

	irreplabels = chkpt_->rd_irr_labs();
	nsopi_ = chkpt_->rd_nsopi_();  
	nmopi_ = chkpt_->rd_orbspi(); 
	doccpi_ = chkpt_->rd_clsdpi();
	soccpi_ = chkpt_->rd_openpi();
	frzcpi_ = chkpt_->rd_frzcpi_();
	frzvpi_ = chkpt_->rd_frzvpi_();
        */

        nso_     = reference_wavefunction_->nso();
        nirrep_ = reference_wavefunction_->nirrep();
        nmo_     = reference_wavefunction_->nmo();
        nmopi_    = reference_wavefunction_->nmopi();
        nsopi_    = reference_wavefunction_->nsopi();
        doccpi_  = reference_wavefunction_->doccpi();
        soccpi_  = reference_wavefunction_->soccpi();
        frzcpi_  = reference_wavefunction_->frzcpi();
        frzvpi_  = reference_wavefunction_->frzvpi();


        // get nfrzc and nfrzv
        nfrzc = 0;
        nfrzv = 0;
        for(int h=0; h<nirrep_; h++) {
	  nfrzc += frzcpi_[h];
	  nfrzv += frzvpi_[h];
	}
	
	// form occpi and virtpi
	occpiA = init_int_array(nirrep_);
	virtpiA = init_int_array(nirrep_);
	memset(occpiA,0, sizeof(int)*nirrep_);
	memset(virtpiA,0, sizeof(int)*nirrep_);
	for(int h=0; h<nirrep_; h++) {
	  virtpiA[h] = nmopi_[h] - doccpi_[h];
	  occpiA[h] = doccpi_[h];
	}
	
	//active occ and virt
	adoccpi = init_int_array(nirrep_);
	aoccpiA = init_int_array(nirrep_);
	avirtpiA = init_int_array(nirrep_);
	memset(adoccpi,0, sizeof(int)*nirrep_);
	memset(aoccpiA,0, sizeof(int)*nirrep_);
	memset(avirtpiA,0, sizeof(int)*nirrep_);
	for(int h=0; h<nirrep_; h++) {
	  adoccpi[h] = doccpi_[h] - frzcpi_[h];
	  avirtpiA[h] = virtpiA[h] - frzvpi_[h];
	  aoccpiA[h] = doccpi_[h] - frzcpi_[h];
	}
	
	// Read in nuclear repulsion energy
	//Enuc = chkpt_->rd_enuc();
	Enuc = Process::environment.molecule()->nuclear_repulsion_energy();
	
	// Read SCF energy
	//Escf=chkpt_->rd_escf();
        Escf=reference_wavefunction_->reference_energy();
	Eref=Escf;
	Eelec=Escf-Enuc;
	
	/* Build mosym arrays */
	mosym = new int [nmo_];
	memset(mosym,0,sizeof(int)*nmo_);
	for(int h=0, q=0; h < nirrep_; h++){
	  for(int p=0; p < nmopi_[h]; p++){
	    mosym[q++] = h;
	  }
	}
	
	/* Build sosym arrays */
	sosym = new int [nso_];
	memset(sosym,0,sizeof(int)*nmo_);
	for(int h=0, q=0; h < nirrep_; h++){
	  for(int p=0; p < nsopi_[h]; p++){
	    sosym[q++] = h;
	  }
	}

	// find nooA
	nooA=0;
	for(int h=0; h < nirrep_; h++){
	  for(int p=0; p < doccpi_[h]; p++){
	    nooA++;
	  }
	}

	// PitzerOffset 
	PitzerOffset = new int[nirrep_];
	memset(PitzerOffset,0,sizeof(int)*nirrep_);
	for(int h=1; h < nirrep_; h++){
	  PitzerOffset[h] = PitzerOffset[h-1] + nmopi_[h-1];
	}
	
	nvoA=nmo_-nooA;   	// Number of virtual orbitals
	nacooA=nooA-nfrzc; 	// Number of active occupied orbitals
	nacso=nmo_-nfrzc-nfrzv; 	// Number of active  orbitals
	nacvoA=nvoA-nfrzv; 	// Number of active virtual orbitals
	npop=nmo_-nfrzv;         // Number of populated orbitals

	ntri_so = 0.5*nso_*(nso_+1);
        ntri = 0.5*nmo_*(nmo_+1);
	dimtei = 0.5*ntri*(ntri+1);

/********************************************************************************************/
/************************** pitzer2symblk ***************************************************/
/********************************************************************************************/
      pitzer2symirrep = new int[nmo_];
      pitzer2symblk = new int[nmo_];
      occ2symblkA = new int[nooA];
      virt2symblkA = new int[nvoA];
      memset(pitzer2symirrep,0,sizeof(int)*nmo_);
      memset(pitzer2symblk,0,sizeof(int)*nmo_);
      memset(occ2symblkA,0,sizeof(int)*nooA);
      memset(virt2symblkA,0,sizeof(int)*nvoA);

      // pitzer2symblk
      int ij,myoffset;
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<nmopi_[h]; ++i) {
            pitzer2symirrep[ij] = h;
            pitzer2symblk[ij] = ij-myoffset;
            ij++;
        }
        myoffset += nmopi_[h];
      }
      
      // occ2symblkA
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<occpiA[h]; ++i) {
            occ2symblkA[ij] = ij-myoffset;
            ij++;
        }
        myoffset += occpiA[h];
      }
      
      
      // vir2symblkA
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<virtpiA[h]; ++i) {
            virt2symblkA[ij] = ij-myoffset;
            ij++;
        }
        myoffset += virtpiA[h];
      }
    
/********************************************************************************************/
/************************** occ_off & vir_off ***********************************************/
/********************************************************************************************/ 
    occ_offA = new int[nirrep_];
    vir_offA = new int[nirrep_];
    memset(occ_offA, 0, sizeof(int)*nirrep_);
    memset(vir_offA, 0, sizeof(int)*nirrep_);
    int ocountA = occpiA[0]; 
    int vcountA = virtpiA[0];
    for(int h=1; h < nirrep_; h++) {
      occ_offA[h] = ocountA;
      ocountA += occpiA[h];
      vir_offA[h] = vcountA;
      vcountA += virtpiA[h];
    }
      
/********************************************************************************************/
/************************** Read orbital coefficients ***************************************/
/********************************************************************************************/
        // read orbital coefficients from chkpt
	Ca_ = SharedMatrix(reference_wavefunction_->Ca());
	Ca_ref = boost::shared_ptr<Matrix>(new Matrix("Ref alpha MO coefficients", nirrep_, nsopi_, nmopi_));
	
	// read orbital coefficients from external files
	if (read_mo_coeff == "TRUE"){
	  fprintf(outfile,"\n\tReading MO coefficients in pitzer order from external files CmoA.psi...\n");  
	  fflush(outfile);
	  double **C_pitzerA = block_matrix(nso_,nmo_);
	  memset(C_pitzerA[0], 0, sizeof(double)*nso_*nmo_);
	
	  // read binary data
	  ifstream InFile1;
	  InFile1.open("CmoA.psi", ios::in | ios::binary);
	  InFile1.read( (char*)C_pitzerA[0], sizeof(double)*nso_*nmo_);
	  InFile1.close();
	  
	  //set C_scf
	  Ca_->set(C_pitzerA);
	  free_block(C_pitzerA);
        }
        
        // Build Reference MOs
        Ca_ref->copy(Ca_);
	
	if(print_ > 1) {
	  Ca_->print();
	}

/********************************************************************************************/
/************************** Create all required matrice *************************************/
/********************************************************************************************/
        // Build Hso
	Hso = boost::shared_ptr<Matrix>(new Matrix("SO-basis One-electron Ints", nirrep_, nsopi_, nsopi_));
	Tso = boost::shared_ptr<Matrix>(new Matrix("SO-basis Kinetic Energy Ints", nirrep_, nsopi_, nsopi_));
	Vso = boost::shared_ptr<Matrix>(new Matrix("SO-basis Potential Energy Ints", nirrep_, nsopi_, nsopi_));
	Hso->zero();
	Tso->zero();
	Vso->zero();
	
	// Read SO-basis one-electron integrals
	double *so_ints = init_array(ntri_so);
        IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_T, so_ints, ntri_so, 0, 0, outfile);
        Tso->set(so_ints);
        IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_V, so_ints, ntri_so, 0, 0, outfile);
        Vso->set(so_ints);
        free(so_ints);
	Hso->copy(Tso); 
	Hso->add(Vso);
	
}// end if (reference == "RHF") 

  
//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference == "UHF") {


/********************************************************************************************/
/************************** MO info *********************************************************/
/********************************************************************************************/
	// Read in mo info
	/*
	nso_ = chkpt_->rd_nso_();
	nmo_ = chkpt_->rd_nmo_();
	nao = chkpt_->rd_nao();
	nfrzc = chkpt_->rd_nfzc();
	nfrzv = chkpt_->rd_nfzv();
	nirrep_ = chkpt_->rd_nirrep_();

	irreplabels = chkpt_->rd_irr_labs();
	nsopi_ = chkpt_->rd_nsopi_();  
	nmopi_ = chkpt_->rd_orbspi(); 
	doccpi_ = chkpt_->rd_clsdpi();
	soccpi_ = chkpt_->rd_openpi();
	frzcpi_ = chkpt_->rd_frzcpi_();
	frzvpi_ = chkpt_->rd_frzvpi_();
        */

        nirrep_ = reference_wavefunction_->nirrep();
        nso_     = reference_wavefunction_->nso();
        nmo_     = reference_wavefunction_->nmo();
        nmopi_    = reference_wavefunction_->nmopi();
        nsopi_    = reference_wavefunction_->nsopi();
        doccpi_  = reference_wavefunction_->doccpi();
        soccpi_  = reference_wavefunction_->soccpi();
        frzcpi_  = reference_wavefunction_->frzcpi();
        frzvpi_  = reference_wavefunction_->frzvpi();

        // get nfrzc and nfrzv
        nfrzc = 0;
        nfrzv = 0;
        for(int h=0; h<nirrep_; h++) {
	  nfrzc += frzcpi_[h];
	  nfrzv += frzvpi_[h];
	}
	
	// form occpi and virtpi
	occpiA = init_int_array(nirrep_);
	occpiB = init_int_array(nirrep_);
	virtpiA = init_int_array(nirrep_);
	virtpiB = init_int_array(nirrep_);
	memset(occpiA,0, sizeof(int)*nirrep_);
	memset(occpiB,0, sizeof(int)*nirrep_);
	memset(virtpiA,0, sizeof(int)*nirrep_);
	memset(virtpiB,0, sizeof(int)*nirrep_);
	for(int h=0; h<nirrep_; h++) {
	  virtpiA[h] = nmopi_[h] - soccpi_[h] - doccpi_[h];
	  virtpiB[h] = nmopi_[h] - doccpi_[h];
	  occpiB[h] = doccpi_[h];
	  occpiA[h] = doccpi_[h] + soccpi_[h];
	}
	
	//active occ and virt
	adoccpi = init_int_array(nirrep_);
	aoccpiA = init_int_array(nirrep_);
	aoccpiB = init_int_array(nirrep_);
	avirtpiA = init_int_array(nirrep_);
	avirtpiB = init_int_array(nirrep_);
	memset(adoccpi,0, sizeof(int)*nirrep_);
	memset(aoccpiA,0, sizeof(int)*nirrep_);
	memset(aoccpiB,0, sizeof(int)*nirrep_);
	memset(avirtpiA,0, sizeof(int)*nirrep_);
	memset(avirtpiB,0, sizeof(int)*nirrep_);
	for(int h=0; h<nirrep_; h++) {
	  adoccpi[h] = doccpi_[h] - frzcpi_[h];
	  avirtpiA[h] = virtpiA[h] - frzvpi_[h];
	  avirtpiB[h] = virtpiB[h] - frzvpi_[h];
	  aoccpiB[h] = doccpi_[h] - frzcpi_[h];
	  aoccpiA[h] = doccpi_[h] + soccpi_[h] - frzcpi_[h];
	}


	// Read in nuclear repulsion energy
	//Enuc = chkpt_->rd_enuc();
	Enuc = Process::environment.molecule()->nuclear_repulsion_energy();
	
	// Read SCF energy
	//Escf=chkpt_->rd_escf();
        Escf=reference_wavefunction_->reference_energy();
	Eref=Escf;
	Eelec=Escf-Enuc;
	
	/* Build mosym arrays */
	mosym = new int [nmo_];
	memset(mosym,0,sizeof(int)*nmo_);
	for(int h=0, q=0; h < nirrep_; h++){
	  for(int p=0; p < nmopi_[h]; p++){
	    mosym[q++] = h;
	  }
	}
	
	/* Build sosym arrays */
	sosym = new int [nso_];
	memset(sosym,0,sizeof(int)*nmo_);
	for(int h=0, q=0; h < nirrep_; h++){
	  for(int p=0; p < nsopi_[h]; p++){
	    sosym[q++] = h;
	  }
	}

	// find nooB
	nooB=0;
	for(int h=0; h < nirrep_; h++){
	  for(int p=0; p < doccpi_[h]; p++){
	    nooB++;
	  }
	}
	
	// find nooA 
	nooA=nooB;
	for(int h=0; h < nirrep_; h++){
	  for(int p=0; p < soccpi_[h]; p++){
	    nooA++;
	  }
	}
	
	
	// PitzerOffset 
	PitzerOffset = new int[nirrep_];
	memset(PitzerOffset,0,sizeof(int)*nirrep_);
	for(int h=1; h < nirrep_; h++){
	  PitzerOffset[h] = PitzerOffset[h-1] + nmopi_[h-1];
	}
	
	nvoA=nmo_-nooA;   	// Number of virtual orbitals
	nvoB=nmo_-nooB;   	// Number of virtual orbitals
	nacooA=nooA-nfrzc; 	// Number of active occupied orbitals
	nacooB=nooB-nfrzc; 	// Number of active occupied orbitals
	nacso=nmo_-nfrzc-nfrzv; 	// Number of active  orbitals
	nacvoA=nvoA-nfrzv; 	// Number of active virtual orbitals
	nacvoB=nvoB-nfrzv; 	// Number of active virtual orbitals
	npop=nmo_-nfrzv;         // Number of populated orbitals

	ntri_so = 0.5*nso_*(nso_+1);
        ntri = 0.5*nmo_*(nmo_+1);
	dimtei = 0.5*ntri*(ntri+1);

/********************************************************************************************/
/************************** pitzer2symblk ***************************************************/
/********************************************************************************************/
      pitzer2symirrep = new int[nmo_];
      pitzer2symblk = new int[nmo_];
      occ2symblkA = new int[nooA];
      occ2symblkB = new int[nooB];
      virt2symblkA = new int[nvoA];
      virt2symblkB = new int[nvoB];
      memset(pitzer2symirrep,0,sizeof(int)*nmo_);
      memset(pitzer2symblk,0,sizeof(int)*nmo_);
      memset(occ2symblkA,0,sizeof(int)*nooA);
      memset(occ2symblkB,0,sizeof(int)*nooB);
      memset(virt2symblkA,0,sizeof(int)*nvoA);
      memset(virt2symblkB,0,sizeof(int)*nvoB);

      // pitzer2symblk
      int ij,myoffset;
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<nmopi_[h]; ++i) {
            pitzer2symirrep[ij] = h;
            pitzer2symblk[ij] = ij-myoffset;
            ij++;
        }
        myoffset += nmopi_[h];
      }
      
      // occ2symblkA
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<occpiA[h]; ++i) {
            occ2symblkA[ij] = ij-myoffset;
            ij++;
        }
        myoffset += occpiA[h];
      }
      
      // occ2symblkB
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<occpiB[h]; ++i) {
            occ2symblkB[ij] = ij-myoffset;
            ij++;
        }
        myoffset += occpiB[h];
      }
      
      // vir2symblkA
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<virtpiA[h]; ++i) {
            virt2symblkA[ij] = ij-myoffset;
            ij++;
        }
        myoffset += virtpiA[h];
      }
      
      // vir2symblkB
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<virtpiB[h]; ++i) {
            virt2symblkB[ij] = ij-myoffset;
            ij++;
        }
        myoffset += virtpiB[h];
      }
    
/********************************************************************************************/
/************************** occ_off & vir_off ***********************************************/
/********************************************************************************************/ 
    occ_offA = new int[nirrep_];
    occ_offB = new int[nirrep_];
    vir_offA = new int[nirrep_];
    vir_offB = new int[nirrep_];
    memset(occ_offA, 0, sizeof(int)*nirrep_);
    memset(occ_offB, 0, sizeof(int)*nirrep_);
    memset(vir_offA, 0, sizeof(int)*nirrep_);
    memset(vir_offB, 0, sizeof(int)*nirrep_);
    int ocountA = occpiA[0]; 
    int ocountB = occpiB[0]; 
    int vcountA = virtpiA[0];
    int vcountB = virtpiB[0];
    for(int h=1; h < nirrep_; h++) {
      occ_offA[h] = ocountA;
      occ_offB[h] = ocountB;
      ocountA += occpiA[h];
      ocountB += occpiB[h];
      
      vir_offA[h] = vcountA;
      vir_offB[h] = vcountB;
      vcountA += virtpiA[h];
      vcountB += virtpiB[h];
    }
      
/********************************************************************************************/
/************************** Read orbital coefficients ***************************************/
/********************************************************************************************/
        // read orbital coefficients from chkpt
	Ca_ = SharedMatrix(reference_wavefunction_->Ca());
        Cb_ = SharedMatrix(reference_wavefunction_->Cb());
	Ca_ref = boost::shared_ptr<Matrix>(new Matrix("Ref alpha MO coefficients", nirrep_, nsopi_, nmopi_));
	Cb_ref = boost::shared_ptr<Matrix>(new Matrix("Ref beta MO coefficients", nirrep_, nsopi_, nmopi_));
	
	// read orbital coefficients from external files
	if (read_mo_coeff == "TRUE"){
	  fprintf(outfile,"\n\tReading MO coefficients in pitzer order from external files CmoA.psi and CmoB.psi...\n");  
	  fflush(outfile);
	  double **C_pitzerA = block_matrix(nso_,nmo_);
	  double **C_pitzerB = block_matrix(nso_,nmo_);
	  memset(C_pitzerA[0], 0, sizeof(double)*nso_*nmo_);
	  memset(C_pitzerB[0], 0, sizeof(double)*nso_*nmo_);
	
	  // read binary data
	  ifstream InFile1;
	  InFile1.open("CmoA.psi", ios::in | ios::binary);
	  InFile1.read( (char*)C_pitzerA[0], sizeof(double)*nso_*nmo_);
	  InFile1.close();
	  
	  // read binary data
	  ifstream InFile2;
	  InFile2.open("CmoB.psi", ios::in | ios::binary);
	  InFile2.read( (char*)C_pitzerB[0], sizeof(double)*nso_*nmo_);
	  InFile2.close();
	
	  //set C_scf
	  Ca_->set(C_pitzerA);
	  Cb_->set(C_pitzerB);
	  
	  free_block(C_pitzerA);
	  free_block(C_pitzerB);
        }
        
        // Build Reference MOs
        Ca_ref->copy(Ca_);
	Cb_ref->copy(Cb_);
	
	if(print_ > 1) {
	  Ca_->print();
	  Cb_->print();
	}

/********************************************************************************************/
/************************** Create all required matrice *************************************/
/********************************************************************************************/
        // Build Hso
	Hso = boost::shared_ptr<Matrix>(new Matrix("SO-basis One-electron Ints", nirrep_, nsopi_, nsopi_));
	Tso = boost::shared_ptr<Matrix>(new Matrix("SO-basis Kinetic Energy Ints", nirrep_, nsopi_, nsopi_));
	Vso = boost::shared_ptr<Matrix>(new Matrix("SO-basis Potential Energy Ints", nirrep_, nsopi_, nsopi_));
	Hso->zero();
	Tso->zero();
	Vso->zero();
	
	// Read SO-basis one-electron integrals
	double *so_ints = init_array(ntri_so);
        IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_T, so_ints, ntri_so, 0, 0, outfile);
        Tso->set(so_ints);
        IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_V, so_ints, ntri_so, 0, 0, outfile);
        Vso->set(so_ints);
        free(so_ints);
	Hso->copy(Tso); 
	Hso->add(Vso);
	
}// end if (reference == "UHF") 

	//fprintf(outfile,"\n get_moinfo is done. \n"); fflush(outfile);

/********************************************************************************************/
/********************************************************************************************/

}
}} // End Namespaces

