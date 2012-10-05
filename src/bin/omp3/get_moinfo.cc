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

#include "ccfiles.h"

#include "omp3wave.h"
#include "defines.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace omp3wave{

void OMP3Wave::get_moinfo()
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
	nso = chkpt_->rd_nso();
	nmo = chkpt_->rd_nmo();
	nao = chkpt_->rd_nao();
	nfrzc = chkpt_->rd_nfzc();
	nfrzv = chkpt_->rd_nfzv();
	nirreps = chkpt_->rd_nirreps();

	irreplabels = chkpt_->rd_irr_labs();
	sopi = chkpt_->rd_sopi();  
	mopi = chkpt_->rd_orbspi(); 
	doccpi = chkpt_->rd_clsdpi();
	soccpi = chkpt_->rd_openpi();
	frzcpi = chkpt_->rd_frzcpi();
	frzvpi = chkpt_->rd_frzvpi();
        */

        nso     = reference_wavefunction_->nso();
        nirreps = reference_wavefunction_->nirrep();
        nmo     = reference_wavefunction_->nmo();
        mopi    = reference_wavefunction_->nmopi();
        sopi    = reference_wavefunction_->nsopi();
        doccpi  = reference_wavefunction_->doccpi();
        soccpi  = reference_wavefunction_->soccpi();
        frzcpi  = reference_wavefunction_->frzcpi();
        frzvpi  = reference_wavefunction_->frzvpi();

        // get nfrzc and nfrzv
        nfrzc = 0;
        nfrzv = 0;
        for(int h=0; h<nirreps; h++) {
	  nfrzc += frzcpi[h];
	  nfrzv += frzvpi[h];
	}
	
	// form occpi and virtpi
	occpiA = init_int_array(nirreps);
	virtpiA = init_int_array(nirreps);
	memset(occpiA,0, sizeof(int)*nirreps);
	memset(virtpiA,0, sizeof(int)*nirreps);
	for(int h=0; h<nirreps; h++) {
	  virtpiA[h] = mopi[h] - doccpi[h];
	  occpiA[h] = doccpi[h];
	}
	
	//active occ and virt
	adoccpi = init_int_array(nirreps);
	aoccpiA = init_int_array(nirreps);
	avirtpiA = init_int_array(nirreps);
	memset(adoccpi,0, sizeof(int)*nirreps);
	memset(aoccpiA,0, sizeof(int)*nirreps);
	memset(avirtpiA,0, sizeof(int)*nirreps);
	for(int h=0; h<nirreps; h++) {
	  adoccpi[h] = doccpi[h] - frzcpi[h];
	  avirtpiA[h] = virtpiA[h] - frzvpi[h];
	  aoccpiA[h] = doccpi[h] - frzcpi[h];
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
	mosym = new int [nmo];
	memset(mosym,0,sizeof(int)*nmo);
	for(int h=0, q=0; h < nirreps; h++){
	  for(int p=0; p < mopi[h]; p++){
	    mosym[q++] = h;
	  }
	}
	
	/* Build sosym arrays */
	sosym = new int [nso];
	memset(sosym,0,sizeof(int)*nmo);
	for(int h=0, q=0; h < nirreps; h++){
	  for(int p=0; p < sopi[h]; p++){
	    sosym[q++] = h;
	  }
	}

	// find nooA
	nooA=0;
	for(int h=0; h < nirreps; h++){
	  for(int p=0; p < doccpi[h]; p++){
	    nooA++;
	  }
	}

	// PitzerOffset 
	PitzerOffset = new int[nirreps];
	memset(PitzerOffset,0,sizeof(int)*nirreps);
	for(int h=1; h < nirreps; h++){
	  PitzerOffset[h] = PitzerOffset[h-1] + mopi[h-1];
	}
	
	nvoA=nmo-nooA;   	// Number of virtual orbitals
	nacooA=nooA-nfrzc; 	// Number of active occupied orbitals
	nacso=nmo-nfrzc-nfrzv; 	// Number of active  orbitals
	nacvoA=nvoA-nfrzv; 	// Number of active virtual orbitals
	npop=nmo-nfrzv;         // Number of populated orbitals

	ntri_so = 0.5*nso*(nso+1);
        ntri = 0.5*nmo*(nmo+1);
	dimtei = 0.5*ntri*(ntri+1);

/********************************************************************************************/
/************************** pitzer2symblk ***************************************************/
/********************************************************************************************/
      pitzer2symirrep = new int[nmo];
      pitzer2symblk = new int[nmo];
      occ2symblkA = new int[nooA];
      virt2symblkA = new int[nvoA];
      memset(pitzer2symirrep,0,sizeof(int)*nmo);
      memset(pitzer2symblk,0,sizeof(int)*nmo);
      memset(occ2symblkA,0,sizeof(int)*nooA);
      memset(virt2symblkA,0,sizeof(int)*nvoA);

      // pitzer2symblk
      int ij,myoffset;
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirreps; ++h) {
        for (int i=0; i<mopi[h]; ++i) {
            pitzer2symirrep[ij] = h;
            pitzer2symblk[ij] = ij-myoffset;
            ij++;
        }
        myoffset += mopi[h];
      }
      
      // occ2symblkA
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirreps; ++h) {
        for (int i=0; i<occpiA[h]; ++i) {
            occ2symblkA[ij] = ij-myoffset;
            ij++;
        }
        myoffset += occpiA[h];
      }
      
      
      // vir2symblkA
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirreps; ++h) {
        for (int i=0; i<virtpiA[h]; ++i) {
            virt2symblkA[ij] = ij-myoffset;
            ij++;
        }
        myoffset += virtpiA[h];
      }
    
/********************************************************************************************/
/************************** occ_off & vir_off ***********************************************/
/********************************************************************************************/ 
    occ_offA = new int[nirreps];
    vir_offA = new int[nirreps];
    memset(occ_offA, 0, sizeof(int)*nirreps);
    memset(vir_offA, 0, sizeof(int)*nirreps);
    int ocountA = occpiA[0]; 
    int vcountA = virtpiA[0];
    for(int h=1; h < nirreps; h++) {
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
	Ca_ref = boost::shared_ptr<Matrix>(new Matrix("Ref alpha MO coefficients", nirreps, sopi, mopi));
	
	// read orbital coefficients from external files
	if (read_mo_coeff == "TRUE"){
	  fprintf(outfile,"\n\tReading MO coefficients in pitzer order from external files CmoA.psi...\n");  
	  fflush(outfile);
	  double **C_pitzerA = block_matrix(nso,nmo);
	  memset(C_pitzerA[0], 0, sizeof(double)*nso*nmo);
	
	  // read binary data
	  ifstream InFile1;
	  InFile1.open("CmoA.psi", ios::in | ios::binary);
	  InFile1.read( (char*)C_pitzerA[0], sizeof(double)*nso*nmo);
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
	Hso = boost::shared_ptr<Matrix>(new Matrix("SO-basis One-electron Ints", nirreps, sopi, sopi));
	Tso = boost::shared_ptr<Matrix>(new Matrix("SO-basis Kinetic Energy Ints", nirreps, sopi, sopi));
	Vso = boost::shared_ptr<Matrix>(new Matrix("SO-basis Potential Energy Ints", nirreps, sopi, sopi));
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
	nso = chkpt_->rd_nso();
	nmo = chkpt_->rd_nmo();
	nao = chkpt_->rd_nao();
	nfrzc = chkpt_->rd_nfzc();
	nfrzv = chkpt_->rd_nfzv();
	nirreps = chkpt_->rd_nirreps();

	irreplabels = chkpt_->rd_irr_labs();
	sopi = chkpt_->rd_sopi();  
	mopi = chkpt_->rd_orbspi(); 
	doccpi = chkpt_->rd_clsdpi();
	soccpi = chkpt_->rd_openpi();
	frzcpi = chkpt_->rd_frzcpi();
	frzvpi = chkpt_->rd_frzvpi();
        */

        nso     = reference_wavefunction_->nso();
        nirreps = reference_wavefunction_->nirrep();
        nmo     = reference_wavefunction_->nmo();
        mopi    = reference_wavefunction_->nmopi();
        sopi    = reference_wavefunction_->nsopi();
        doccpi  = reference_wavefunction_->doccpi();
        soccpi  = reference_wavefunction_->soccpi();
        frzcpi  = reference_wavefunction_->frzcpi();
        frzvpi  = reference_wavefunction_->frzvpi();

        // get nfrzc and nfrzv
        nfrzc = 0;
        nfrzv = 0;
        for(int h=0; h<nirreps; h++) {
	  nfrzc += frzcpi[h];
	  nfrzv += frzvpi[h];
	}
	
	// form occpi and virtpi
	occpiA = init_int_array(nirreps);
	occpiB = init_int_array(nirreps);
	virtpiA = init_int_array(nirreps);
	virtpiB = init_int_array(nirreps);
	memset(occpiA,0, sizeof(int)*nirreps);
	memset(occpiB,0, sizeof(int)*nirreps);
	memset(virtpiA,0, sizeof(int)*nirreps);
	memset(virtpiB,0, sizeof(int)*nirreps);
	for(int h=0; h<nirreps; h++) {
	  virtpiA[h] = mopi[h] - soccpi[h] - doccpi[h];
	  virtpiB[h] = mopi[h] - doccpi[h];
	  occpiB[h] = doccpi[h];
	  occpiA[h] = doccpi[h] + soccpi[h];
	}
	
	//active occ and virt
	adoccpi = init_int_array(nirreps);
	aoccpiA = init_int_array(nirreps);
	aoccpiB = init_int_array(nirreps);
	avirtpiA = init_int_array(nirreps);
	avirtpiB = init_int_array(nirreps);
	memset(adoccpi,0, sizeof(int)*nirreps);
	memset(aoccpiA,0, sizeof(int)*nirreps);
	memset(aoccpiB,0, sizeof(int)*nirreps);
	memset(avirtpiA,0, sizeof(int)*nirreps);
	memset(avirtpiB,0, sizeof(int)*nirreps);
	for(int h=0; h<nirreps; h++) {
	  adoccpi[h] = doccpi[h] - frzcpi[h];
	  avirtpiA[h] = virtpiA[h] - frzvpi[h];
	  avirtpiB[h] = virtpiB[h] - frzvpi[h];
	  aoccpiB[h] = doccpi[h] - frzcpi[h];
	  aoccpiA[h] = doccpi[h] + soccpi[h] - frzcpi[h];
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
	mosym = new int [nmo];
	memset(mosym,0,sizeof(int)*nmo);
	for(int h=0, q=0; h < nirreps; h++){
	  for(int p=0; p < mopi[h]; p++){
	    mosym[q++] = h;
	  }
	}
	
	/* Build sosym arrays */
	sosym = new int [nso];
	memset(sosym,0,sizeof(int)*nmo);
	for(int h=0, q=0; h < nirreps; h++){
	  for(int p=0; p < sopi[h]; p++){
	    sosym[q++] = h;
	  }
	}

	// find nooB
	nooB=0;
	for(int h=0; h < nirreps; h++){
	  for(int p=0; p < doccpi[h]; p++){
	    nooB++;
	  }
	}
	
	// find nooA 
	nooA=nooB;
	for(int h=0; h < nirreps; h++){
	  for(int p=0; p < soccpi[h]; p++){
	    nooA++;
	  }
	}
	
	
	// PitzerOffset 
	PitzerOffset = new int[nirreps];
	memset(PitzerOffset,0,sizeof(int)*nirreps);
	for(int h=1; h < nirreps; h++){
	  PitzerOffset[h] = PitzerOffset[h-1] + mopi[h-1];
	}
	
	nvoA=nmo-nooA;   	// Number of virtual orbitals
	nvoB=nmo-nooB;   	// Number of virtual orbitals
	nacooA=nooA-nfrzc; 	// Number of active occupied orbitals
	nacooB=nooB-nfrzc; 	// Number of active occupied orbitals
	nacso=nmo-nfrzc-nfrzv; 	// Number of active  orbitals
	nacvoA=nvoA-nfrzv; 	// Number of active virtual orbitals
	nacvoB=nvoB-nfrzv; 	// Number of active virtual orbitals
	npop=nmo-nfrzv;         // Number of populated orbitals

	ntri_so = 0.5*nso*(nso+1);
        ntri = 0.5*nmo*(nmo+1);
	dimtei = 0.5*ntri*(ntri+1);

/********************************************************************************************/
/************************** pitzer2symblk ***************************************************/
/********************************************************************************************/
      pitzer2symirrep = new int[nmo];
      pitzer2symblk = new int[nmo];
      occ2symblkA = new int[nooA];
      occ2symblkB = new int[nooB];
      virt2symblkA = new int[nvoA];
      virt2symblkB = new int[nvoB];
      memset(pitzer2symirrep,0,sizeof(int)*nmo);
      memset(pitzer2symblk,0,sizeof(int)*nmo);
      memset(occ2symblkA,0,sizeof(int)*nooA);
      memset(occ2symblkB,0,sizeof(int)*nooB);
      memset(virt2symblkA,0,sizeof(int)*nvoA);
      memset(virt2symblkB,0,sizeof(int)*nvoB);

      // pitzer2symblk
      int ij,myoffset;
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirreps; ++h) {
        for (int i=0; i<mopi[h]; ++i) {
            pitzer2symirrep[ij] = h;
            pitzer2symblk[ij] = ij-myoffset;
            ij++;
        }
        myoffset += mopi[h];
      }
      
      // occ2symblkA
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirreps; ++h) {
        for (int i=0; i<occpiA[h]; ++i) {
            occ2symblkA[ij] = ij-myoffset;
            ij++;
        }
        myoffset += occpiA[h];
      }
      
      // occ2symblkB
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirreps; ++h) {
        for (int i=0; i<occpiB[h]; ++i) {
            occ2symblkB[ij] = ij-myoffset;
            ij++;
        }
        myoffset += occpiB[h];
      }
      
      // vir2symblkA
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirreps; ++h) {
        for (int i=0; i<virtpiA[h]; ++i) {
            virt2symblkA[ij] = ij-myoffset;
            ij++;
        }
        myoffset += virtpiA[h];
      }
      
      // vir2symblkB
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirreps; ++h) {
        for (int i=0; i<virtpiB[h]; ++i) {
            virt2symblkB[ij] = ij-myoffset;
            ij++;
        }
        myoffset += virtpiB[h];
      }
    
/********************************************************************************************/
/************************** occ_off & vir_off ***********************************************/
/********************************************************************************************/ 
    occ_offA = new int[nirreps];
    occ_offB = new int[nirreps];
    vir_offA = new int[nirreps];
    vir_offB = new int[nirreps];
    memset(occ_offA, 0, sizeof(int)*nirreps);
    memset(occ_offB, 0, sizeof(int)*nirreps);
    memset(vir_offA, 0, sizeof(int)*nirreps);
    memset(vir_offB, 0, sizeof(int)*nirreps);
    int ocountA = occpiA[0]; 
    int ocountB = occpiB[0]; 
    int vcountA = virtpiA[0];
    int vcountB = virtpiB[0];
    for(int h=1; h < nirreps; h++) {
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
	Ca_ref = boost::shared_ptr<Matrix>(new Matrix("Ref alpha MO coefficients", nirreps, sopi, mopi));
	Cb_ref = boost::shared_ptr<Matrix>(new Matrix("Ref beta MO coefficients", nirreps, sopi, mopi));
	
	// read orbital coefficients from external files
	if (read_mo_coeff == "TRUE"){
	  fprintf(outfile,"\n\tReading MO coefficients in pitzer order from external files CmoA.psi and CmoB.psi...\n");  
	  fflush(outfile);
	  double **C_pitzerA = block_matrix(nso,nmo);
	  double **C_pitzerB = block_matrix(nso,nmo);
	  memset(C_pitzerA[0], 0, sizeof(double)*nso*nmo);
	  memset(C_pitzerB[0], 0, sizeof(double)*nso*nmo);
	
	  // read binary data
	  ifstream InFile1;
	  InFile1.open("CmoA.psi", ios::in | ios::binary);
	  InFile1.read( (char*)C_pitzerA[0], sizeof(double)*nso*nmo);
	  InFile1.close();
	  
	  // read binary data
	  ifstream InFile2;
	  InFile2.open("CmoB.psi", ios::in | ios::binary);
	  InFile2.read( (char*)C_pitzerB[0], sizeof(double)*nso*nmo);
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
	Hso = boost::shared_ptr<Matrix>(new Matrix("SO-basis One-electron Ints", nirreps, sopi, sopi));
	Tso = boost::shared_ptr<Matrix>(new Matrix("SO-basis Kinetic Energy Ints", nirreps, sopi, sopi));
	Vso = boost::shared_ptr<Matrix>(new Matrix("SO-basis Potential Energy Ints", nirreps, sopi, sopi));
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

