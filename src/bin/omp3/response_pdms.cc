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

        // Build G intermediates 
        timer_on("G int");
	G_int(); 
        timer_off("G int");

 if (reference == "RHF") {
        // Initialize
	gamma1corr->zero();
	g1symm->zero();

        // OPDM
        timer_on("OPDM");
	// OO-block alpha contrb.
	#pragma omp parallel for
	for(int h = 0; h < nirreps; ++h){
	  for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
                g1symm->set(h, i, j, GooA->get(h, i, j) + GooA->get(h, j, i));
            }
	  }
	}

	// VV-block alpha contrb.
        #pragma omp parallel for
	for(int h = 0; h < nirreps; ++h){
	  for(int a = 0 ; a < avirtpiA[h]; ++a){
            for(int b = 0 ; b < avirtpiA[h]; ++b){
                int aa = a + occpiA[h];
                int bb = b + occpiA[h];
                g1symm->set(h, aa, bb, GvvA->get(h, a, b) + GvvA->get(h, b, a));
            }
	  }
	}

	g1symm->scale(-1.0);
	gamma1corr->copy(g1symm); // correlation opdm
  
        // REF contribution 
	// alpha contrb.
        #pragma omp parallel for
	for(int h=0; h<nirreps; h++) {
	  if (occpiA[h] != 0) {
	    for (int i=0; i<occpiA[h];i++) {
	      g1symm->add(h,i,i,2.0);
	    }
	  }
	}
        timer_off("OPDM");

        //print
        if (print_ > 1) {
	  g1symm->print();
        }

        // TPDM
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
        timer_on("TPDM REF");
	twopdm_ref(); 
        timer_off("TPDM REF");
        timer_on("TPDM CORR OPDM");
	twopdm_corr_opdm();
        timer_off("TPDM CORR OPDM");
 }// end if (reference == "RHF") 

 else if (reference == "UHF") {
        // Initialize
	gamma1corrA->zero();
	gamma1corrB->zero();
	g1symmA->zero();
	g1symmB->zero();

        // OPDM
        timer_on("OPDM");
	// OO-block alpha contrb.
	#pragma omp parallel for
	for(int h = 0; h < nirreps; ++h){
	  for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
                g1symmA->set(h, i, j, GooA->get(h, i, j) + GooA->get(h, j, i));
            }
	  }
	}

	// OO-block beta contrb.
	#pragma omp parallel for
	for(int h = 0; h < nirreps; ++h){
	  for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int j = 0 ; j < aoccpiB[h]; ++j){
                g1symmB->set(h, i, j, GooB->get(h, i, j) + GooB->get(h, j, i));
            }
	  }
	}

	// VV-block alpha contrb.
        #pragma omp parallel for
	for(int h = 0; h < nirreps; ++h){
	  for(int a = 0 ; a < avirtpiA[h]; ++a){
            for(int b = 0 ; b < avirtpiA[h]; ++b){
                int aa = a + occpiA[h];
                int bb = b + occpiA[h];
                g1symmA->set(h, aa, bb, GvvA->get(h, a, b) + GvvA->get(h, b, a));
            }
	  }
	}

        // VV-block beta contrb.
        #pragma omp parallel for
	for(int h = 0; h < nirreps; ++h){
	  for(int a = 0 ; a < avirtpiB[h]; ++a){
            for(int b = 0 ; b < avirtpiB[h]; ++b){
                int aa = a + occpiB[h];
                int bb = b + occpiB[h];
                g1symmB->set(h, aa, bb, GvvB->get(h, a, b) + GvvB->get(h, b, a));
            }
	  }
	}

	g1symmA->scale(-0.5);
	g1symmB->scale(-0.5);
	gamma1corrA->copy(g1symmA); // correlation opdm
	gamma1corrB->copy(g1symmB); // correlation opdm
  
        // REF contribution 
	// alpha contrb.
        #pragma omp parallel for
	for(int h=0; h<nirreps; h++) {
	  if (occpiA[h] != 0) {
	    for (int i=0; i<occpiA[h];i++) {
	      g1symmA->add(h,i,i,1.0);
	    }
	  }
	}
	
	// beta contrb.
        #pragma omp parallel for
	for(int h=0; h<nirreps; h++) {
	  if (occpiB[h] != 0) {
	    for (int i=0; i<occpiB[h];i++) {
	      g1symmB->add(h,i,i,1.0);
	    }
	  }
	}
        timer_off("OPDM");

        //print
        if (print_ > 1) {
	  g1symmA->print();
	  g1symmB->print();
        }

        // TPDM
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



 }// end if (reference == "UHF") 

  //fprintf(outfile,"\n response_pdms done... \n"); fflush(outfile);

} // end of response_pdms


/*=======================*/
/*  twopdm_oovv()        */
/*=======================*/
void OMP3Wave::twopdm_oovv()
{      

    dpdbuf4 G, T, V, Tau;
    
    psio_->open(PSIF_OMP3_DPD, PSIO_OPEN_OLD);  
    psio_->open(PSIF_OMP3_DENSITY, PSIO_OPEN_OLD);

 if (reference == "RHF") {
    // G (IJ,AB) = 1/4 V(IB,AJ)  
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,V]"), ID("[V,O]"),
                  ID("[O,V]"), ID("[V,O]"), 0, "V_2 <OV|VO>");
    dpd_buf4_sort(&V, PSIF_OMP3_DENSITY , psrq, ID("[O,O]"), ID("[V,V]"), "TPDM <OO|VV>");
    dpd_buf4_close(&V);

    // G (IJ,AB) += 1/4 (2T_IJ^AB - T_JI^AB)
    dpd_buf4_init(&Tau, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau <OO|VV>");
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    dpd_buf4_axpy(&Tau, &G, 1.0); // 1.0*Tau + G -> G
    dpd_buf4_close(&Tau);
    dpd_buf4_scm(&G, 0.25);
    dpd_buf4_close(&G);
 }// end if (reference == "RHF") 

 else if (reference == "UHF") {
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
    
 }// end if (reference == "UHF") 

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
    
 if (reference == "RHF") {
    // G(IJ,KL) = 1/8 (V_IJKL +  V_ILKJ)
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "V_2 <OO|OO>");
    dpd_buf4_sort(&V, PSIF_OMP3_DENSITY , psrq, ID("[O,O]"), ID("[O,O]"), "TPDM <OO|OO>");
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
    dpd_buf4_axpy(&V, &G, 1.0); // 1.0*V + G -> G
    dpd_buf4_close(&V);
    dpd_buf4_scm(&G, 0.125);
    dpd_buf4_close(&G);
 }// end if (reference == "RHF") 

 else if (reference == "UHF") {
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
    
 }// end if (reference == "UHF") 
    
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
   
 if (reference == "RHF") {
    // NOTE: A VVVV-sort takes too long time, hence I will not symmetrize the 
    // TPDM VVVV-block. However, in future for analytical gradients 
    // I need to symmetrize it.
    
    // G_ABCD(2) = 1/2\sum_{M,N} T_MN^CD(1) (2T_MN^AB(1) - T_MN^BA(1))
    dpd_buf4_init(&T, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_init(&L, PSIF_OMP3_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau_1 <OO|VV>");
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");    
    dpd_contract444(&L, &T, &G, 1, 1, 0.5, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_close(&G);

    //Print 
    if (print_ > 3) {
      dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>"); 
      dpd_buf4_print(&G, outfile, 1);
      dpd_buf4_close(&G);
    }  
 }// end if (reference == "RHF") 

 else if (reference == "UHF") {
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
 }// end if (reference == "UHF") 
    
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
    
 if (reference == "RHF") {
    // Fully-symmetric
    // G(IA,JB) = -1/4 (V_IAJB +  V_IBJA)
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "V_2 <OV|OV>");
    dpd_buf4_sort(&V, PSIF_OMP3_DENSITY , psrq, ID("[O,V]"), ID("[O,V]"), "TPDM <OV|OV>");
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    dpd_buf4_axpy(&V, &G, 1.0); // 1.0*V + G -> G
    dpd_buf4_close(&V);
    dpd_buf4_scm(&G, -0.25);
    dpd_buf4_close(&G);

    /*
    // Partially-symmetric
    // G(IA,JB) = -1/2 V_IAJB 
    dpd_buf4_init(&V, PSIF_OMP3_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "V_2 <OV|OV>");
    dpd_buf4_copy(&V, PSIF_OMP3_DENSITY , "TPDM <OV|OV>");
    dpd_buf4_close(&V);
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    dpd_buf4_scm(&G, -0.5);
    dpd_buf4_close(&G);
    */
 }// end if (reference == "RHF") 

 else if (reference == "UHF") {
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
    
 }// end if (reference == "UHF") 
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

 if (reference == "RHF") {
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
	dpd_buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
		if (i == k && j == l) G.matrix[h][ij][kl] += 1.0;
		if (i == l && j == k) G.matrix[h][ij][kl] -= 0.25;
		if (i == j && k == l) G.matrix[h][ij][kl] -= 0.25;
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    
    dpd_buf4_close(&G);
    
 }// end RHF 

 else if (reference == "UHF") {
    // Alpha-Alpha spin case
    dpd_buf4_init(&G, PSIF_OMP3_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
	dpd_buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
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
        #pragma omp parallel for
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
        #pragma omp parallel for
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
    
 }// end UHF 

    psio_->close(PSIF_OMP3_DENSITY, 1);
    
} // end of twopdm_ref

}} // End Namespaces




