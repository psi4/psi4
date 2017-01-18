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

#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/matrix.h"
#include "defines.h"
#include "occwave.h"


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
     Emp2_t1 = 0.0;

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
     global_dpd_->buf4_init(&Tss, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TAA <OO|VV>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     Emp2AA = 0.5 * global_dpd_->buf4_dot(&Tss, &K);
     global_dpd_->buf4_close(&Tss);

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
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T <OO|VV>");
     Emp2AB = global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

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
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     Emp2AA = 0.25 * global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

     Escsmp2AA = ss_scale * Emp2AA;
     Escsnmp2AA = 1.76 * Emp2AA;
     Escsmimp2AA = 1.29 * Emp2AA;
     Escsmp2vdwAA = 0.50 * Emp2AA;


     // Alpha-Beta spin contribution
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
     Emp2AB = global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

     Escsmp2AB = os_scale * Emp2AB;
     if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB;
     else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;
     Escsmimp2AB = 0.40 * Emp2AB;
     Escsmp2vdwAB = 1.28 * Emp2AB;
     Esospimp2AB = 1.40 * Emp2AB;


     // Beta-Beta spin contribution
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     Emp2BB = 0.25 * global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

     Escsmp2BB = ss_scale * Emp2BB;
     Escsnmp2BB = 1.76 * Emp2BB;
     Escsmimp2BB = 1.29 * Emp2BB;
     Escsmp2vdwBB = 0.50 * Emp2BB;

 if (reference == "ROHF" && orb_opt_ == "FALSE" && wfn_type_ == "OMP2") {
    // Singles-contribution
    // Alpha
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int a = 0 ; a < avirtpiA[h]; ++a){
                Emp2_t1 += t1A->get(h, i, a) * FockA->get(h, a + occpiA[h], i + frzcpi_[h]);
            }
        }
    }

    // beta
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int a = 0 ; a < avirtpiB[h]; ++a){
                Emp2_t1 += t1B->get(h, i, a) * FockB->get(h, a + occpiB[h], i + frzcpi_[h]);
            }
        }
    }
 }// end if (reference == "ROHF")


 }// end uhf

     Ecorr = Emp2AA + Emp2AB + Emp2BB + Emp2_t1;
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
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1AA <OO|VV>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     Emp2AA = 0.5 * global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
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
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
     Emp2AB = global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

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
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     Ecorr += 0.25 * global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

     Emp2AA = Ecorr;
     Escsmp2AA = ss_scale * Emp2AA;
     Escsnmp2AA = 1.76 * Emp2AA;
     Escsmimp2AA = 1.29 * Emp2AA;
     Escsmp2vdwAA = 0.50 * Emp2AA;


     // Alpha-Beta spin contribution
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
     Ecorr += global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

     Emp2AB = Ecorr - Emp2AA;
     Escsmp2AB = os_scale * Emp2AB;
     Esosmp2AB = sos_scale * Emp2AB;
     //if (mo_optimized == 0) Esosmp2AB = sos_scale * Emp2AB;
     //else if (mo_optimized == 1) Esosmp2AB = sos_scale2 * Emp2AB;
     Escsmimp2AB = 0.40 * Emp2AB;
     Escsmp2vdwAB = 1.28 * Emp2AB;
     Esospimp2AB = 1.40 * Emp2AB;


     // Beta-Beta spin contribution
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     Ecorr += 0.25 * global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

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
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2AA <OO|VV>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     Emp3AA = 0.5 * global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     Emp3BB = Emp3AA;


     // Alpha-Beta spin contribution
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     Emp3AB = global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

 }// end rhf

 else if (reference_ == "UNRESTRICTED") {
     // Compute Energy
     // Alpha-Alpha spin contribution
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     Emp3AA = 0.25 * global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);


     // Alpha-Beta spin contribution
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
     Emp3AB = global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

     // Beta-Beta spin contribution
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     Emp3BB = 0.25 * global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

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
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2AA <OO|VV>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     Emp2AA = 0.5 * global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
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
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     Emp2AB = global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

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
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     Ecorr += 0.25 * global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

     Emp2AA = Ecorr;
     Escsmp2AA = ss_scale * Emp2AA;
     Escsnmp2AA = 1.76 * Emp2AA;
     Escsmimp2AA = 1.29 * Emp2AA;
     Escsmp2vdwAA = 0.50 * Emp2AA;


     // Alpha-Beta spin contribution
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
     Ecorr += global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

     Emp2AB = Ecorr - Emp2AA;
     Escsmp2AB = os_scale * Emp2AB;
     Esosmp2AB = sos_scale * Emp2AB;
     Escsmimp2AB = 0.40 * Emp2AB;
     Escsmp2vdwAB = 1.28 * Emp2AB;
     Esospimp2AB = 1.40 * Emp2AB;


     // Beta-Beta spin contribution
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     Ecorr += 0.25 * global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

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
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2AA <OO|VV>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     EcepaAA = 0.5 * global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     EcepaBB = EcepaAA;

     // Alpha-Beta spin contribution
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     EcepaAB = global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

 }// end rhf

 else if (reference_ == "UNRESTRICTED") {
     // Compute Energy
     // Alpha-Alpha spin contribution
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     EcepaAA = 0.25 * global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);


     // Alpha-Beta spin contribution
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
     EcepaAB = global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

     // Beta-Beta spin contribution
     global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     EcepaBB = 0.25 * global_dpd_->buf4_dot(&T, &K);
     global_dpd_->buf4_close(&T);
     global_dpd_->buf4_close(&K);

 }// end uhf

     Ecorr = EcepaAA + EcepaBB + EcepaAB;
     Ecepa = Eref + Ecorr;
     Escscepa = Eref + ((cepa_ss_scale_ * (EcepaAA + EcepaBB)) + (cepa_os_scale_ * EcepaAB));
     Esoscepa = Eref + (cepa_sos_scale_ * EcepaAB);

     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OCC_DPD, 1);

} // end of cepa_energy


}} // End Namespaces
