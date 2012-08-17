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
  
/********************************************************************************************/
/************************** MO info *********************************************************/
/********************************************************************************************/
	// Read in mo info
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
	Enuc = chkpt_->rd_enuc();
	
	// Read SCF energy
	Escf=chkpt_->rd_escf();
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
/************************** Copy all required matrice from chk file *************************/
/********************************************************************************************/	
	//read orbital energies
	evalsA = chkpt_->rd_alpha_evals(); 
	evalsB = chkpt_->rd_beta_evals(); 
	
      // remove degeneracies
      for (int p=1; p<nmo; p++){
	for (int q=0; q<p; q++){
	  if (evalsA[p] == evalsA[q]) {
	    evalsA[p]+=1e-4;
	    evalsA[q]-=1e-4;
	  }
	}
      }  
      
      // remove degeneracies
      for (int p=1; p<nmo; p++){
	for (int q=0; q<p; q++){
	  if (evalsB[p] == evalsB[q]) {
	    evalsB[p]+=1e-4;
	    evalsB[q]-=1e-4;
	  }
	}
      }  
	
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
/************************** c1topitzer ******************************************************/
/********************************************************************************************/
      c1topitzerA = new int[nmo];
      c1topitzerB = new int[nmo];
      pitzer2c1A = new int[nmo];
      pitzer2c1B = new int[nmo];
      memset(c1topitzerA,0,sizeof(int)*nmo);
      memset(c1topitzerB,0,sizeof(int)*nmo);
      memset(pitzer2c1A,0,sizeof(int)*nmo);
      memset(pitzer2c1B,0,sizeof(int)*nmo);
     
      // form evals_c1
      evals_c1A = new double[nmo];
      evals_c1B = new double[nmo];
      memset(evals_c1A,0,sizeof(double)*nmo);
      memset(evals_c1B,0,sizeof(double)*nmo);
      for(int p=0; p<nmo; p++) {
	evals_c1A[p] = evalsA[p];
	evals_c1B[p] = evalsB[p];
      }
      
      // sort evals A
      for(int p=0; p<nmo; p++) {
	for(int q=nmo-1; q>p; q--) {
	  if (evals_c1A[q-1] > evals_c1A[q]) {
	    double dum = evals_c1A[q-1];    
	    evals_c1A[q-1] = evals_c1A[q];    
	    evals_c1A[q] = dum;
	  }
	}
      }
     
      //set to reg order  
      for (int p=0; p<nmo; p++){
	for (int q=0; q<nmo; q++){
	  if (evalsA[p] == evals_c1A[q]) {
	    pitzer2c1A[p]=q;
	    c1topitzerA[q]=p;
	  }
	}
      }  
      
      // sort evals B
      for(int p=0; p<nmo; p++) {
	for(int q=nmo-1; q>p; q--) {
	  if (evals_c1B[q-1] > evals_c1B[q]) {
	    double dum = evals_c1B[q-1];    
	    evals_c1B[q-1] = evals_c1B[q];    
	    evals_c1B[q] = dum;
	  }
	}
      }
     
      //set to reg order  
      for (int p=0; p<nmo; p++){
	for (int q=0; q<nmo; q++){
	  if (evalsB[p] == evals_c1B[q]) {
	    pitzer2c1B[p]=q;
	    c1topitzerB[q]=p;
	  }
	}
      }  
      
  
/********************************************************************************************/
/************************** c1toqt **********************************************************/
/********************************************************************************************/ 
      /*
      // c1 to qt
      c1toqtA = new int[nmo];
      c1toqtB = new int[nmo];
      qt2c1A = new int[nmo];
      qt2c1B = new int[nmo];
      memset(c1toqtA,0,sizeof(int)*nmo);
      memset(c1toqtB,0,sizeof(int)*nmo);
      memset(qt2c1A,0,sizeof(int)*nmo);
      memset(qt2c1B,0,sizeof(int)*nmo);
      
      // Alpha spin case
      
      //frzc
      int itemp=0;
      if (nfrzc != 0){
      for(int h=0; h<nirreps; h++){ 
	if (h > 0) itemp+=frzcpi[h-1];
	for (int p=0; p<frzcpi[h]; p++){
	  if (frzcpi[h] != 0){
	    int pp = p + PitzerOffset[h];  // convert sym-block order to pitzer order
	    int p2 = pitzer2c1A[pp];
	    c1toqtA[p2]=p+itemp;
	    qt2c1A[p+itemp]=p2;
	  }
	}
      }  
      }      
      
      //adocc
      itemp=nfrzc;
      for(int h=0; h<nirreps; h++){ 
	if (h > 0) itemp+=adoccpiA[h-1];
	for (int p=0; p<adoccpiA[h]; p++){
	  if (aoccpiA[h] != 0){
	    int pp = p + PitzerOffset[h] + frzcpi[h];  // convert sym-block order to pitzer order
	    int p2 = pitzer2c1A[pp];
	    c1toqtA[p2]=p+itemp;
	    qt2c1A[p+itemp]=p2;
	  }
	}
      }         
      
      //avirt
      itemp=nooA;
      for(int h=0; h<nirreps; h++){ 
	if (h > 0) itemp+=avirtpiA[h-1];
	for (int p=0; p<avirtpiA[h]; p++){
	   if (avirtpiA[h] != 0){
	      int pp = p + PitzerOffset[h] + doccpi[h];  // convert sym-block order to pitzer order
	      int p2 = pitzer2c1A[pp];
	      c1toqtA[p2]=p+itemp;
	      qt2c1A[p+itemp]=p2;
	   }
	}
      }  
      
      
      //frzv
      if (nfrzv != 0){
      itemp=npop;
      for(int h=0; h<nirreps; h++){ 
	if (h > 0) itemp+=frzvpi[h-1];
	for (int p=0; p<frzvpi[h]; p++){
	   if (frzvpi[h] != 0){
	      int pp = p + PitzerOffset[h] + avirtpi[h] + doccpi[h];  // convert sym-block order to pitzer order
	      int p2 = pitzer2c1[pp];
	      c1toqt[p2]=p+itemp;
	      qt2c1[p+itemp]=p2;
	   }
	}
      } 
      }
      */
      
/********************************************************************************************/
/************************** mosym_c1 ********************************************************/
/********************************************************************************************/ 
	/* Build mosym arrays */
	mosym_c1 = new int [nmo];
	memset(mosym_c1,0,sizeof(int)*nmo);
	for(int p=0; p<nmo; p++) {
	  int p2 = c1topitzerA[p];
	  mosym_c1[p] = mosym[p2];
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
	
	//fprintf(outfile,"\n get_moinfo is done. \n"); fflush(outfile);

/********************************************************************************************/
/********************************************************************************************/

}
}} // End Namespaces

