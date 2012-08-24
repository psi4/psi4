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
#include <libtrans/mospace.h>
#include <libtrans/integraltransform.h>

/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>

#include "omp3wave.h"
#include "defines.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace omp3wave{

void OMP3Wave::response_pdms()
{   
     //fprintf(outfile,"\n response_pdms is starting... \n"); fflush(outfile);
/********************************************************************************************/
/************************** Initialize ******************************************************/
/********************************************************************************************/ 
	gamma1corrA->zero();
	gamma1corrB->zero();
	g1symmA->zero();
	g1symmB->zero();

/********************************************************************************************/
/************************** Build G intermediates *******************************************/
/********************************************************************************************/  
        timer_on("G int");
	G_int(); 
        timer_off("G int");
   
/********************************************************************************************/
/************************** (OO) Block whole ************************************************/
/********************************************************************************************/ 
        timer_on("OPDM");
	// alpha contrb.
	for (int i=nfrzc; i<nooA;i++) {
	  for (int j=nfrzc; j<nooA;j++) { 
	    int i2 = c1topitzerA[i];
	    int j2 = c1topitzerA[j];
	    int i_sym = pitzer2symblk[i2];
	    int j_sym = pitzer2symblk[j2];
	    int h1=mosym[i2];
	    int h2=mosym[j2];
	    if (h1 == h2){
	      int i3 = i_sym - frzcpi[h1];
	      int j3 = j_sym - frzcpi[h1];
	      g1symmA->set(h1,i_sym,j_sym,GooA->get(h1,i3,j3) + GooA->get(h1,j3,i3));
	    }
	  }
	}
	
	// beta contrb.
	for (int i=nfrzc; i<nooB;i++) {
	  for (int j=nfrzc; j<nooB;j++) { 
	    int i2 = c1topitzerB[i];
	    int j2 = c1topitzerB[j];
	    int i_sym = pitzer2symblk[i2];
	    int j_sym = pitzer2symblk[j2];
	    int h1=mosym[i2];
	    int h2=mosym[j2];
	    if (h1 == h2){
	      int i3 = i_sym - frzcpi[h1];
	      int j3 = j_sym - frzcpi[h1];
	      g1symmB->set(h1,i_sym,j_sym,GooB->get(h1,i3,j3) + GooB->get(h1,j3,i3));
	    }
	  }
	}
	
/********************************************************************************************/
/************************** (VV) Block whole ************************************************/
/********************************************************************************************/ 
	 // alpha contrb.
	 for (int a=nooA; a<npop;a++) {
	  for (int b=nooA; b<npop;b++) { 
	    int a2 = c1topitzerA[a];
	    int b2 = c1topitzerA[b];
	    int a_sym = pitzer2symblk[a2];
	    int b_sym = pitzer2symblk[b2];
	    int h1=mosym[a2];
	    int h2=mosym[b2];
	    if (h1 == h2){
	      int a3 = a_sym - occpiA[h1];
	      int b3 = b_sym - occpiA[h1];
	      g1symmA->set(h1,a_sym,b_sym,GvvA->get(h1,a3,b3) + GvvA->get(h1,b3,a3));
	    }
	  }
	 }
	 
	 // beta contrb.
	 for (int a=nooB; a<npop;a++) {
	  for (int b=nooB; b<npop;b++) { 
	    int a2 = c1topitzerB[a];
	    int b2 = c1topitzerB[b];
	    int a_sym = pitzer2symblk[a2];
	    int b_sym = pitzer2symblk[b2];
	    int h1=mosym[a2];
	    int h2=mosym[b2];
	    if (h1 == h2){
	      int a3 = a_sym - occpiB[h1];
	      int b3 = b_sym - occpiB[h1];
	      g1symmB->set(h1,a_sym,b_sym,GvvB->get(h1,a3,b3) + GvvB->get(h1,b3,a3));
	    }
	  }
	 }
	 
	g1symmA->scale(-0.5);
	g1symmB->scale(-0.5);
	gamma1corrA->copy(g1symmA); // correlation opdm
	gamma1corrB->copy(g1symmB); // correlation opdm
   
/********************************************************************************************/
/************************** HF Contribution *************************************************/
/********************************************************************************************/
	// alpha contrb.
	for(int h=0; h<nirreps; h++) {
	  if (occpiA[h] != 0) {
	    for (int i=0; i<occpiA[h];i++) {
	      g1symmA->add(h,i,i,1.0);
	    }
	  }
	}
	
	// beta contrb.
	for(int h=0; h<nirreps; h++) {
	  if (occpiB[h] != 0) {
	    for (int i=0; i<occpiB[h];i++) {
	      g1symmB->add(h,i,i,1.0);
	    }
	  }
	}
        timer_off("OPDM");

/********************************************************************************************/
/************************** TPDMs ***********************************************************/
/********************************************************************************************/ 
        timer_on("V int");
        V_2nd_order(); 
        timer_off("V int");
        timer_on("TPDM OOVV");
	twopdm_oovv();
        timer_off("TPDM OOVV");
        timer_on("TPDM OOOO");
	twopdm_oooo();
        timer_off("TPDM OOOO");

        if (twopdm_abcd_type == "COMPUTE") {
           timer_on("TPDM VVVV");
           twopdm_vvvv();
           timer_off("TPDM VVVV");
        }

        timer_on("TPDM OVOV");
        twopdm_ovov();
        timer_off("TPDM OVOV");
        timer_on("TPDM VOVO");
        twopdm_vovo();
        timer_off("TPDM VOVO");
        timer_on("TPDM OVVO");
        twopdm_ovvo();
        timer_off("TPDM OVVO");
        timer_on("TPDM REF");
	twopdm_ref(); 
        timer_off("TPDM REF");
        timer_on("TPDM CORR OPDM");
	twopdm_corr_opdm();
        timer_off("TPDM CORR OPDM");

/********************************************************************************************/
/************************** Print ***********************************************************/
/********************************************************************************************/  
      if (print_ > 1) {
	g1symmA->print();
	g1symmB->print();
      }
      //fprintf(outfile,"\n response_pdms done... \n"); fflush(outfile);

} // end of response_pdms


/*=======================*/
/*  twopdm_oovv()        */
/*=======================*/
void OMP3Wave::twopdm_oovv()
{      
    dpdbuf4 G, T;
    
    psio_->open(PSIF_OMP3_DPD, PSIO_OPEN_OLD);  
    psio_->open(PSIF_OMP3_DENSITY, PSIO_OPEN_OLD);
    
    // Alpha-Alpha spin case
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    dpd_buf4_copy(&T, PSIF_OMP3_DENSITY, "TPDM <OO|VV>");
    dpd_buf4_close(&T);
    
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    dpd_buf4_scm(&G, 0.25);
    dpd_buf4_close(&G);
    
    // Beta-Beta spin case
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
    dpd_buf4_copy(&T, PSIF_OMP3_DENSITY, "TPDM <oo|vv>");
    dpd_buf4_close(&T);
    
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
    dpd_buf4_scm(&G, 0.25);
    dpd_buf4_close(&G);
    
    // Alpha-Beta spin case
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    dpd_buf4_copy(&T, PSIF_OMP3_DENSITY, "TPDM <Oo|Vv>");
    dpd_buf4_close(&T);
    
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
    dpd_buf4_scm(&G, 0.25);
    dpd_buf4_close(&G);
    
    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
      dpd_buf4_print(&G, outfile, 1);
      dpd_buf4_close(&G);
      
      dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
      dpd_buf4_print(&G, outfile, 1);
      dpd_buf4_close(&G);
      
      dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
      dpd_buf4_print(&G, outfile, 1);
      dpd_buf4_close(&G);
    }    
    
    
    psio_->close(PSIF_OMP3_DPD, 1);  
    psio_->close(PSIF_OMP3_DENSITY, 1);

} // end of twopdm_oovv



/*=======================*/
/*  twopdm_oooo()        */
/*=======================*/
void OMP3Wave::twopdm_oooo()
{      
    dpdbuf4 G, T, V;
     
    psio_->open(PSIF_OMP3_DENSITY, PSIO_OPEN_OLD);
    
    // Alpha-Alpha spin case
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "V_2 <OO|OO>");
    dpd_buf4_copy(&V, PSIF_OMP3_DENSITY, "TPDM <OO|OO>");
    dpd_buf4_close(&V);
    
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
    dpd_buf4_scm(&G, 0.25);
    dpd_buf4_close(&G);
    
    // Beta-Beta spin case
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "V_2 <oo|oo>");
    dpd_buf4_copy(&V, PSIF_OMP3_DENSITY, "TPDM <oo|oo>");
    dpd_buf4_close(&V);
    
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");
    dpd_buf4_scm(&G, 0.25);
    dpd_buf4_close(&G);
    
    // Alpha-Beta spin case
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "V_2 <Oo|Oo>");
    dpd_buf4_copy(&V, PSIF_OMP3_DENSITY, "TPDM <Oo|Oo>");
    dpd_buf4_close(&V);
    
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
    dpd_buf4_scm(&G, 0.25);
    dpd_buf4_close(&G);
    
    
    psio_->close(PSIF_OMP3_DENSITY, 1);

} // end of twopdm_oooo


/*=======================*/
/*  twopdm_vvvv()        */
/*=======================*/
void OMP3Wave::twopdm_vvvv()
{      
    // NOTE: contract444 can handle only TN and NT type contractions, which means (0,0) and (1,1) type target indices,
    //  with out-of-core algorithm!!!!     
    dpdbuf4  T, L, G, V;
     
    psio_->open(PSIF_OMP3_DPD, PSIO_OPEN_OLD); 
    psio_->open(PSIF_OMP3_DENSITY, PSIO_OPEN_OLD);
   
    // Alpha-Alpha spin case
    // G_ABCD(2) = 1/8 \sum_{M,N} T_MN^CD(1) L_AB^MN(1) = 1/8 \sum_{M,N} T_MN^AB(1) T_MN^CD(1)
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_init(&L, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");    
    dpd_contract444(&L, &T, &G, 1, 1, 0.125, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_close(&G);
    
    // Beta-Beta spin case
    // G_abcd(2) = 1/8 \sum_{m,n} T_mn^cd(1) L_ab^mn(1) = 1/8 \sum_{m,n} T_mn^ab(1) T_mn^cd(1)
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    dpd_buf4_init(&L, PSIF_OMP3_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
    dpd_contract444(&L, &T, &G, 1, 1, 0.125, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_close(&G);
    
    
    // Alpha-Beta spin case
    // G_AbCd(2) = 1/4 \sum_{M,n} T_Mn^Cd(1) L_Ab^Mn(1) = 1/4 \sum_{M,n} T_Mn^Ab(1) T_Mn^Cd(1)
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_init(&L, PSIF_OMP3_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
    dpd_contract444(&L, &T, &G, 1, 1, 0.25, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_close(&G);
    

    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>"); 
      dpd_buf4_print(&G, outfile, 1);
      dpd_buf4_close(&G);
      
      dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
      dpd_buf4_print(&G, outfile, 1);
      dpd_buf4_close(&G);
      
      dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
      dpd_buf4_print(&G, outfile, 1);
      dpd_buf4_close(&G);
    }    
    
    psio_->close(PSIF_OMP3_DENSITY, 1);
    psio_->close(PSIF_OMP3_DPD, 1);

} // end of twopdm_vvvv


/*=======================*/
/*  twopdm_ovov()        */
/*=======================*/
void OMP3Wave::twopdm_ovov()
{      
    dpdbuf4 G, T, V;
     
    psio_->open(PSIF_OMP3_DENSITY, PSIO_OPEN_OLD);
    
    // Build G_IAJB
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "V_2 <OV|OV>");
    dpd_buf4_copy(&V, PSIF_OMP3_DENSITY, "TPDM <OV|OV>");
    dpd_buf4_close(&V);
    
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    dpd_buf4_scm(&G, -0.5);
    dpd_buf4_close(&G);
    
    // Build G_iajb
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "V_2 <ov|ov>");
    dpd_buf4_copy(&V, PSIF_OMP3_DENSITY, "TPDM <ov|ov>");
    dpd_buf4_close(&V);
    
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
    dpd_buf4_scm(&G, -0.5);
    dpd_buf4_close(&G);
    
    // Build G_IaJb
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "V_2 <Ov|Ov>");
    dpd_buf4_copy(&V, PSIF_OMP3_DENSITY, "TPDM <Ov|Ov>");
    dpd_buf4_close(&V);
    
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
    dpd_buf4_scm(&G, -0.5);
    dpd_buf4_close(&G);
    
    psio_->close(PSIF_OMP3_DENSITY, 1);

} // end of twopdm_ovov


/*=======================*/
/*  twopdm_vovo()        */
/*=======================*/
void OMP3Wave::twopdm_vovo()
{      
    dpdbuf4 G, T, V;
     
    psio_->open(PSIF_OMP3_DENSITY, PSIO_OPEN_OLD);
    
    // G_AIBJ = G_IAJB
    // G_aibj = G_iajb
    
    // G_AiBj = -1/2 V_iAjB
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "V_2 <oV|oV>");
    dpd_buf4_sort(&V, PSIF_OMP3_DENSITY , qpsr, ID("[V,o]"), ID("[V,o]"), "TPDM <Vo|Vo>");
    dpd_buf4_close(&V);
    
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");
    dpd_buf4_scm(&G, -0.5);
    dpd_buf4_close(&G);
    
    psio_->close(PSIF_OMP3_DENSITY, 1);

} // end of twopdm_vovo


/*=======================*/
/*  twopdm_ovvo()        */
/*=======================*/
void OMP3Wave::twopdm_ovvo()
{      
    dpdbuf4 G, T, V;
     
    psio_->open(PSIF_OMP3_DENSITY, PSIO_OPEN_OLD);
    
    // G_IABJ = -G_IAJB
    // G_iabj = -G_iajb
    
    // G_IaBj = 1/2 V_IajB
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "V_2 <Ov|oV>");
    dpd_buf4_sort(&V, PSIF_OMP3_DENSITY , pqsr, ID("[O,v]"), ID("[V,o]"), "TPDM <Ov|Vo>");
    dpd_buf4_close(&V);
    
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
    dpd_buf4_scm(&G, 0.5);
    dpd_buf4_close(&G);
   
    // VoOv block is here! 
    // G_AiJb = G_JbAi so I do not need to VoOv block, however for proper contraction in the GFock.cc I need it. 
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
    dpd_buf4_sort(&G, PSIF_OMP3_DENSITY , rspq, ID("[V,o]"), ID("[O,v]"), "TPDM <Vo|Ov>");
    dpd_buf4_close(&G);

        
    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
      dpd_buf4_print(&G, outfile, 1);
      dpd_buf4_close(&G);
    }
    
    psio_->close(PSIF_OMP3_DENSITY, 1);

} // end of twopdm_ovvo



/*=======================*/
/*  twopdm_ref()         */
/*=======================*/
void OMP3Wave::twopdm_ref()
{
    dpdbuf4 G;
    
    psio_->open(PSIF_OMP3_DENSITY, PSIO_OPEN_OLD);
    
    // Alpha-Alpha spin case
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
	dpd_buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
		if (i == k && j == l) G.matrix[h][ij][kl] += 0.25;
		if (i == l && j == k) G.matrix[h][ij][kl] -= 0.25;
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }    
    dpd_buf4_close(&G);
    
    
    // Beta-Beta spin case
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");    
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
	dpd_buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
		if (i == k && j == l) G.matrix[h][ij][kl] += 0.25;
		if (i == l && j == k) G.matrix[h][ij][kl] -= 0.25;
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }    
    dpd_buf4_close(&G);
    
    
    // Alpha-Beta spin case
     dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                 ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>"); 
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
	dpd_buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
		if (i == k && j == l) G.matrix[h][ij][kl] += 0.25;
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }    
    dpd_buf4_close(&G);
    
    psio_->close(PSIF_OMP3_DENSITY, 1);
    
} // end of twopdm_ref


}} // End Namespaces




