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
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libtrans/mospace.h>
#include <libtrans/integraltransform.h>

/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>

#include "defines.h"
#include "occwave.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace occwave{
  
//=======================//
// ref_energy            //
//=======================//
void OCCWave::ref_energy()
{
     double Ehf;     
     Ehf=0.0;

 if (reference_ == "RESTRICTED") {
    for (int h=0; h<nirrep_; h++){
      for (int i=0; i<occpiA[h];i++) {
	Ehf+=HmoA->get(h,i,i) + FockA->get(h,i,i);
      }
    }         
    Eref = Ehf + Enuc;
 }// end rhf
 
 else if (reference_ == "UNRESTRICTED") { 
     
     // alpha contribution
     for (int h=0; h<nirrep_; h++){
      for (int i=0; i<occpiA[h];i++) {
	Ehf+=HmoA->get(h,i,i) + FockA->get(h,i,i);
      }
    }  
    
    // beta contribution
     for (int h=0; h<nirrep_; h++){
      for (int i=0; i<occpiB[h];i++) {
	Ehf+=HmoB->get(h,i,i) + FockB->get(h,i,i);
      }
    }  
    
    Eref = (0.5 * Ehf) + Enuc; 
 }// end uhf
    
} // end of ref_energy


//=======================//
// omp2_mp2_energy       //
//=======================//
void OCCWave::omp2_mp2_energy()
{
     dpdbuf4 K, T, Tau, Tss;

     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);

     Ecorr = 0.0;

     Escsmp2AA = 0.0;
     Escsmp2AB = 0.0;
     Escsmp2BB = 0.0;
     Escsmp2 = 0.0;

     Esosmp2AB = 0.0;
     Esosmp2 = 0.0;

     Escsnmp2AA = 0.0;
     Escsnmp2BB = 0.0;
     Escsnmp2 = 0.0;
     
     Escsmimp2AA = 0.0;
     Escsmimp2AB = 0.0;
     Escsmimp2BB = 0.0;
     Escsmimp2 = 0.0;

     Escsmp2vdwAA = 0.0;
     Escsmp2vdwAB = 0.0;
     Escsmp2vdwBB = 0.0;
     Escsmp2vdw = 0.0;

     Esospimp2AB = 0.0;
     Esospimp2 = 0.0;

 if (reference_ == "RESTRICTED") {
     // Same-spin contribution
     dpd_buf4_init(&Tss, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TAA <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     Emp2AA = 0.5 * dpd_buf4_dot(&Tss, &K);     
     dpd_buf4_close(&Tss);

     Escsmp2AA = ss_scale * Emp2AA; 
     Escsnmp2AA = 1.76 * Emp2AA; 
     Escsmimp2AA = 1.29 * Emp2AA; 
     Escsmp2vdwAA = 0.5 * Emp2AA; 

     Emp2BB = Emp2AA;    
     Escsmp2BB = ss_scale * Emp2BB;  
     Escsnmp2BB = 1.76 * Emp2BB; 
     Escsmimp2BB = 1.29 * Emp2BB; 
     Escsmp2vdwBB = 0.50 * Emp2BB; 
   
     // Opposite-spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T <OO|VV>");
     Emp2AB = dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);

     Escsmp2AB = os_scale * Emp2AB;  
     if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB; 
     else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;  
     Escsmimp2AB = 0.40 * Emp2AB; 
     Escsmp2vdwAB = 1.28 * Emp2AB; 
     Esospimp2AB = 1.40 * Emp2AB; 
     
 }// end rhf


 else if (reference_ == "UNRESTRICTED") { 

     // Compute Energy
     // Alpha-Alpha spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     Emp2AA = 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Escsmp2AA = ss_scale * Emp2AA; 
     Escsnmp2AA = 1.76 * Emp2AA; 
     Escsmimp2AA = 1.29 * Emp2AA; 
     Escsmp2vdwAA = 0.50 * Emp2AA; 
     
     
     // Alpha-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
     Emp2AB = dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Escsmp2AB = os_scale * Emp2AB;  
     if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB; 
     else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;  
     Escsmimp2AB = 0.40 * Emp2AB; 
     Escsmp2vdwAB = 1.28 * Emp2AB; 
     Esospimp2AB = 1.40 * Emp2AB; 
     
     
     // Beta-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     Emp2BB = 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Escsmp2BB = ss_scale * Emp2BB;  
     Escsnmp2BB = 1.76 * Emp2BB; 
     Escsmimp2BB = 1.29 * Emp2BB; 
     Escsmp2vdwBB = 0.50 * Emp2BB; 
     
 }// end uhf

     Ecorr = Emp2AA + Emp2AB + Emp2BB;
     Emp2 = Eref + Ecorr;
     Escsmp2 = Eref + Escsmp2AA + Escsmp2AB + Escsmp2BB;
     Esosmp2 = Eref + Esosmp2AB;     
     Escsnmp2 = Eref + Escsnmp2AA + Escsnmp2BB;
     Escsmimp2 = Eref + Escsmimp2AA + Escsmimp2AB + Escsmimp2BB;
     Escsmp2vdw = Eref + Escsmp2vdwAA + Escsmp2vdwAB + Escsmp2vdwBB;
     Esospimp2 = Eref + Esospimp2AB;     

     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OCC_DPD, 1);    
} // end of omp2_mp2_energy


//=======================//
// omp3_mp2_energy       //
//=======================//
void OCCWave::omp3_mp2_energy()
{
     dpdbuf4 K, T;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
     
     Ecorr = 0.0;

     Escsmp2AA = 0.0;
     Escsmp2AB = 0.0;
     Escsmp2BB = 0.0;
     Escsmp2 = 0.0;

     Esosmp2AB = 0.0;
     Esosmp2 = 0.0;

     Escsnmp2AA = 0.0;
     Escsnmp2BB = 0.0;
     Escsnmp2 = 0.0;
     
     Escsmimp2AA = 0.0;
     Escsmimp2AB = 0.0;
     Escsmimp2BB = 0.0;
     Escsmimp2 = 0.0;

     Escsmp2vdwAA = 0.0;
     Escsmp2vdwAB = 0.0;
     Escsmp2vdwBB = 0.0;
     Escsmp2vdw = 0.0;

     Esospimp2AB = 0.0;
     Esospimp2 = 0.0;
     
 if (reference_ == "RESTRICTED") {
     // Compute Energy
     // Alpha-Alpha spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1AA <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     Emp2AA = 0.5 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     Emp2BB = Emp2AA;  
     
     Escsmp2AA = ss_scale * Emp2AA; 
     Escsnmp2AA = 1.76 * Emp2AA; 
     Escsmimp2AA = 1.29 * Emp2AA; 
     Escsmp2vdwAA = 0.50 * Emp2AA; 
     
     Escsmp2BB = ss_scale * Emp2BB;  
     Escsnmp2BB = 1.76 * Emp2BB; 
     Escsmimp2BB = 1.29 * Emp2BB; 
     Escsmp2vdwBB = 0.50 * Emp2BB; 
     
     // Alpha-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
     Emp2AB = dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Escsmp2AB = os_scale * Emp2AB;  
     Esosmp2AB = sos_scale * Emp2AB; 
     //if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB; 
     //else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;  
     Escsmimp2AB = 0.40 * Emp2AB; 
     Escsmp2vdwAB = 1.28 * Emp2AB; 
     Esospimp2AB = 1.40 * Emp2AB; 

     Ecorr = Emp2AA + Emp2BB + Emp2AB;
 
 }// end rhf
 
 else if (reference_ == "UNRESTRICTED") { 
     // Compute Energy
     // Alpha-Alpha spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     Ecorr += 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Emp2AA = Ecorr;    
     Escsmp2AA = ss_scale * Emp2AA; 
     Escsnmp2AA = 1.76 * Emp2AA; 
     Escsmimp2AA = 1.29 * Emp2AA; 
     Escsmp2vdwAA = 0.50 * Emp2AA; 
     
     
     // Alpha-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
     Ecorr += dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Emp2AB = Ecorr - Emp2AA;
     Escsmp2AB = os_scale * Emp2AB;  
     Esosmp2AB = sos_scale * Emp2AB; 
     //if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB; 
     //else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;  
     Escsmimp2AB = 0.40 * Emp2AB; 
     Escsmp2vdwAB = 1.28 * Emp2AB; 
     Esospimp2AB = 1.40 * Emp2AB; 
     
     
     // Beta-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     Ecorr += 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Emp2BB = Ecorr - Emp2AA - Emp2AB;  
     Escsmp2BB = ss_scale * Emp2BB;  
     Escsnmp2BB = 1.76 * Emp2BB; 
     Escsmimp2BB = 1.29 * Emp2BB; 
     Escsmp2vdwBB = 0.50 * Emp2BB; 
     
 }// end uhf

     Emp2 = Eref + Ecorr;
     Escsmp2 = Eref + Escsmp2AA + Escsmp2AB + Escsmp2BB;
     Esosmp2 = Eref + Esosmp2AB;     
     Escsnmp2 = Eref + Escsnmp2AA + Escsnmp2BB;
     Escsmimp2 = Eref + Escsmimp2AA + Escsmimp2AB + Escsmimp2BB;
     Escsmp2vdw = Eref + Escsmp2vdwAA + Escsmp2vdwAB + Escsmp2vdwBB;
     Esospimp2 = Eref + Esospimp2AB;     
     
     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OCC_DPD, 1);    
     
} // end of omp3_mp2_energy


/*=======================*/
/*  mp3_energy()         */
/*=======================*/
void OCCWave::mp3_energy()
{
     dpdbuf4 K, T;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
     
     Ecorr = 0.0;
 
 if (reference_ == "RESTRICTED") {
     // Compute Energy
     // Alpha-Alpha spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2AA <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     Emp3AA = 0.5 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     Emp3BB = Emp3AA;    
     
     
     // Alpha-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     Emp3AB = dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);     

 }// end rhf

 else if (reference_ == "UNRESTRICTED") {
     // Compute Energy
     // Alpha-Alpha spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     Emp3AA = 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);     
     
     
     // Alpha-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
     Emp3AB = dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);     
 
     // Beta-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     Emp3BB = 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);     

 }// end uhf
     
     Ecorr = Emp3AA + Emp3BB + Emp3AB;
     Emp3 = Eref + Ecorr;
     Escsmp3 = Escsmp2 + (e3_scale * (Emp3 - Emp2) );
     Esosmp3 = Esosmp2 + (e3_scale * (Emp3 - Emp2) );
     Escsnmp3 = Escsnmp2 + (e3_scale * (Emp3 - Emp2) );
     Escsmimp3 = Escsmimp2 + (e3_scale * (Emp3 - Emp2) );
     Escsmp3vdw = Escsmp2vdw + (e3_scale * (Emp3 - Emp2) );
     Esospimp3 = Esospimp2 + (e3_scale * (Emp3 - Emp2) );
     
     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OCC_DPD, 1);    
     
} // end of mp3_energy


//=======================//
// ocepa_mp2_energy      //
//=======================//
void OCCWave::ocepa_mp2_energy()
{
     dpdbuf4 K, T;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
     
     Ecorr = 0.0;

     Escsmp2AA = 0.0;
     Escsmp2AB = 0.0;
     Escsmp2BB = 0.0;
     Escsmp2 = 0.0;

     Esosmp2AB = 0.0;
     Esosmp2 = 0.0;

     Escsnmp2AA = 0.0;
     Escsnmp2BB = 0.0;
     Escsnmp2 = 0.0;
     
     Escsmimp2AA = 0.0;
     Escsmimp2AB = 0.0;
     Escsmimp2BB = 0.0;
     Escsmimp2 = 0.0;

     Escsmp2vdwAA = 0.0;
     Escsmp2vdwAB = 0.0;
     Escsmp2vdwBB = 0.0;
     Escsmp2vdw = 0.0;

     Esospimp2AB = 0.0;
     Esospimp2 = 0.0;
     
 if (reference_ == "RESTRICTED") {
     // Compute Energy
     // Alpha-Alpha spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2AA <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     Emp2AA = 0.5 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     Emp2BB = Emp2AA;  
     
     Escsmp2AA = ss_scale * Emp2AA; 
     Escsnmp2AA = 1.76 * Emp2AA; 
     Escsmimp2AA = 1.29 * Emp2AA; 
     Escsmp2vdwAA = 0.50 * Emp2AA; 
     
     Escsmp2BB = ss_scale * Emp2BB;  
     Escsnmp2BB = 1.76 * Emp2BB; 
     Escsmimp2BB = 1.29 * Emp2BB; 
     Escsmp2vdwBB = 0.50 * Emp2BB; 
     
     // Alpha-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     Emp2AB = dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Escsmp2AB = os_scale * Emp2AB;  
     Esosmp2AB = sos_scale * Emp2AB; 
     Escsmimp2AB = 0.40 * Emp2AB; 
     Escsmp2vdwAB = 1.28 * Emp2AB; 
     Esospimp2AB = 1.40 * Emp2AB; 
     Ecorr = Emp2AA + Emp2BB + Emp2AB;
 
 }// end rhf
 
 else if (reference_ == "UNRESTRICTED") { 
     // Compute Energy
     // Alpha-Alpha spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     Ecorr += 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Emp2AA = Ecorr;    
     Escsmp2AA = ss_scale * Emp2AA; 
     Escsnmp2AA = 1.76 * Emp2AA; 
     Escsmimp2AA = 1.29 * Emp2AA; 
     Escsmp2vdwAA = 0.50 * Emp2AA; 
     
     
     // Alpha-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
     Ecorr += dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Emp2AB = Ecorr - Emp2AA;
     Escsmp2AB = os_scale * Emp2AB;  
     Esosmp2AB = sos_scale * Emp2AB; 
     Escsmimp2AB = 0.40 * Emp2AB; 
     Escsmp2vdwAB = 1.28 * Emp2AB; 
     Esospimp2AB = 1.40 * Emp2AB; 
     
     
     // Beta-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     Ecorr += 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);
     
     Emp2BB = Ecorr - Emp2AA - Emp2AB;  
     Escsmp2BB = ss_scale * Emp2BB;  
     Escsnmp2BB = 1.76 * Emp2BB; 
     Escsmimp2BB = 1.29 * Emp2BB; 
     Escsmp2vdwBB = 0.50 * Emp2BB; 
     
 }// end uhf

     Emp2 = Eref + Ecorr;
     Escsmp2 = Eref + Escsmp2AA + Escsmp2AB + Escsmp2BB;
     Esosmp2 = Eref + Esosmp2AB;     
     Escsnmp2 = Eref + Escsnmp2AA + Escsnmp2BB;
     Escsmimp2 = Eref + Escsmimp2AA + Escsmimp2AB + Escsmimp2BB;
     Escsmp2vdw = Eref + Escsmp2vdwAA + Escsmp2vdwAB + Escsmp2vdwBB;
     Esospimp2 = Eref + Esospimp2AB;     
     
     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OCC_DPD, 1);    
     
} // end of mp2_energy


/*=======================*/
/*  cepa_energy()         */
/*=======================*/
void OCCWave::cepa_energy()
{
     dpdbuf4 K, T;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
     
     Ecorr = 0.0;
 
 if (reference_ == "RESTRICTED") {
     // Compute Energy
     // Alpha-Alpha spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2AA <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     EcepaAA = 0.5 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     EcepaBB = EcepaAA;    
     
     // Alpha-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     EcepaAB = dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);     

 }// end rhf

 else if (reference_ == "UNRESTRICTED") {
     // Compute Energy
     // Alpha-Alpha spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     EcepaAA = 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);     
     
     
     // Alpha-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
     EcepaAB = dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);     
 
     // Beta-Beta spin contribution
     dpd_buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     EcepaBB = 0.25 * dpd_buf4_dot(&T, &K);     
     dpd_buf4_close(&T);
     dpd_buf4_close(&K);     

 }// end uhf
     
     Ecorr = EcepaAA + EcepaBB + EcepaAB;
     Ecepa = Eref + Ecorr;
     Escscepa = Eref + ((cepa_ss_scale_ * (EcepaAA + EcepaBB)) + (cepa_os_scale_ * EcepaAB));
     Esoscepa = Eref + (cepa_sos_scale_ * EcepaAB);
     
     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OCC_DPD, 1);    
     
} // end of cepa_energy


}} // End Namespaces


