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
#include "defines.h"
#include "occwave.h"


using namespace std;


namespace psi{ namespace occwave{

void OCCWave::ccl_energy()
{
    //outfile->Printf("\n ccl_energy is starting... \n");
    // Two-electron contribution
    dpdbuf4 G, K;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

 if (reference_ == "RESTRICTED") {
    // OOOO-Block contribution
    // E += 2*G_ijkl <ij|kl>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO|OO>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
    Ecc_rdm += 2.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);

  if (wfn_type_ != "OMP2") {
    // VVVV-Block contribution
    // E += 2*G_ABCD <AB|CD>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
    Ecc_rdm += 2.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);
  }// end if (wfn_type_ != "OMP2") {


    // OOVV-Block contribution
    // E += 8*G_IJAB <IJ|AB>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    Ecc_rdm += 8.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);


    // OVOV-Block contribution
    // E += 4*G_IAJB <IA|JB>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    Ecc_rdm += 4.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);

 }// end rhf

 else if (reference_ == "UNRESTRICTED") {
    // OOOO-Block contribution
    // E += G_IJKL <IJ||KL>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO||OO>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
    Ecc_rdm += global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);

    // E += G_ijkl <ij||kl>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo||oo>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");
    Ecc_rdm += global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);

    // E += 4*G_IjKl <Ij||Kl>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
    Ecc_rdm += 4.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);


  if (wfn_type_ != "OMP2") {
    // VVVV-Block contribution
    /*
    // E += G_ABCD <AB||CD>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV||VV>");
    dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
    Ecc_rdm += dpd_buf4_dot(&G, &K);
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    */

     // E += 2*G_ABCD <AB|CD>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
    Ecc_rdm += 2.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);

    /*
    // E += G_abcd <ab||cd>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv||vv>");
    dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
    Ecc_rdm += dpd_buf4_dot(&G, &K);
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    */

     // E += 2*G_abcd <ab|cd>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv|vv>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
    Ecc_rdm += 2.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);

    // E += 4*G_AbCd <Ab||Cd>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
    Ecc_rdm += 4.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);
  }// end if (wfn_type_ != "OMP2") {


    // OOVV-Block contribution
    // E += 2*G_IJAB <IJ||AB>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    Ecc_rdm += 2.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);

    // E += 2*G_ijab <ij|ab>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
    Ecc_rdm += 2.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);

    // E += 8*G_IjAb <Ij||Ab>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
    Ecc_rdm += 8.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);


    // OVOV-Block contribution
    // E += 4*G_IAJB <IA||JB>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV||OV>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    Ecc_rdm += 4.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);

    // E += 4*G_iajb <ia||jb>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov||ov>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
    Ecc_rdm += 4.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);

    // E += 4*G_IaJb <Ia||Jb>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
    Ecc_rdm += 4.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);

    // VOVO-Block contribution
    // E += 4*G_AiBj <Ai||Bj>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0,ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "MO Ints <Vo|Vo>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");
    Ecc_rdm += 4.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);

    // OVVO-Block contribution
    // E += 8*G_IaBj <Ia||Bj>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0,ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "MO Ints <Ov|Vo>");
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
    Ecc_rdm += 8.0 * global_dpd_->buf4_dot(&G, &K);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);

 }// end uhf

    psio_->close(PSIF_LIBTRANS_DPD, 1);
    psio_->close(PSIF_OCC_DENSITY, 1);

  if (wfn_type_ == "OMP2") {
    EcorrL=Ecc_rdm-Escf;
    Emp2L=Ecc_rdm;
    DE = Emp2L - Emp2L_old;
    Emp2L_old = Emp2L;
  }

  else if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
    EcorrL=Ecc_rdm-Escf;
    Emp3L=Ecc_rdm;
    DE = Emp3L - Emp3L_old;
    Emp3L_old = Emp3L;
  }

  else if (wfn_type_ == "OCEPA") {
    EcorrL = Ecc_rdm - Escf;
    EcepaL = Ecc_rdm;
    DE = EcepaL - EcepaL_old;
    EcepaL_old = EcepaL;
  }

    Etpdm = Ecc_rdm - Eopdm;

    /*
    outfile->Printf("\tOPDM energy (a.u.)          : %12.14f\n", Eopdm);
    outfile->Printf("\tTPDM energy (a.u.)          : %12.14f\n", Etpdm);
    outfile->Printf("\tTotal PDM energy (a.u.)     : %12.14f\n", Ecc_rdm);

    */
    //outfile->Printf("\n ccl_energy is done... \n");

} // end of ccl_energy
}} // End Namespaces
