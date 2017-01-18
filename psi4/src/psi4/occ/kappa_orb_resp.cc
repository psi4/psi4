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
#include "psi4/libmints/matrix.h"
#include "psi4/libpsio/psio.hpp"
#include "occwave.h"
#include "defines.h"


using namespace std;


namespace psi{ namespace occwave{

void OCCWave::kappa_orb_resp()
{
//outfile->Printf("\n kappa_orb_resp is starting... \n");

if (reference_ == "RESTRICTED") {
    // Set the kappa to -negative of the mo grad
    kappaA->copy(wogA);
    kappaA->scale(-1.0);

    // Open dpd files
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    dpdbuf4 K;

    // Sort some integrals
    // (OV|OV) -> (VO|VO)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                 ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , qpsr, ID("[V,O]"), ID("[V,O]"), "MO Ints (VO|VO)");
    global_dpd_->buf4_close(&K);

    // (ai|bj) -> (aj|bi)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[V,O]"),
                  ID("[V,O]"), ID("[V,O]"), 0, "MO Ints (VO|VO)");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , psrq, ID("[V,O]"), ID("[V,O]"), "MO Ints (aj|bi)");
    global_dpd_->buf4_close(&K);

    // <OV|OV> -> <VO|VO>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                 ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , qpsr, ID("[V,O]"), ID("[V,O]"), "MO Ints <VO|VO>");
    global_dpd_->buf4_close(&K);


    // Build the MO Hessian
    Aorb = new Array2d("MO Hessian Matrix", nidpA, nidpA);
    Aorb->zero();
    // A(ai,bj) = 8*(ai|bj)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[V,O]"),
                  ID("[V,O]"), ID("[V,O]"), 0, "MO Ints (VO|VO)");
        int h =0;
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        //#pragma omp parallel for
        for(int ai = 0; ai < K.params->rowtot[h]; ++ai){
            int a = K.params->roworb[h][ai][0];
            int i = K.params->roworb[h][ai][1];
            for(int bj = 0; bj < K.params->coltot[h]; ++bj){
                int b = K.params->colorb[h][bj][0];
                int j = K.params->colorb[h][bj][1];
                Aorb->set(ai, bj, 8.0 * K.matrix[h][ai][bj]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    global_dpd_->buf4_close(&K);

    // A(ai,bj) -= 2*<ai|bj>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[V,O]"),
                  ID("[V,O]"), ID("[V,O]"), 0, "MO Ints <VO|VO>");
        h = 0;
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        //#pragma omp parallel for
        for(int ai = 0; ai < K.params->rowtot[h]; ++ai){
            int a = K.params->roworb[h][ai][0];
            int i = K.params->roworb[h][ai][1];
            for(int bj = 0; bj < K.params->coltot[h]; ++bj){
                int b = K.params->colorb[h][bj][0];
                int j = K.params->colorb[h][bj][1];
                Aorb->add(ai, bj, -2.0 * K.matrix[h][ai][bj]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    global_dpd_->buf4_close(&K);

    // A(ai,bj) -= 2*(aj|bi)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[V,O]"),
                  ID("[V,O]"), ID("[V,O]"), 0, "MO Ints (aj|bi)");
        h = 0;
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        //#pragma omp parallel for
        for(int ai = 0; ai < K.params->rowtot[h]; ++ai){
            int a = K.params->roworb[h][ai][0];
            int i = K.params->roworb[h][ai][1];
            for(int bj = 0; bj < K.params->coltot[h]; ++bj){
                int b = K.params->colorb[h][bj][0];
                int j = K.params->colorb[h][bj][1];
                Aorb->add(ai, bj, -2.0 * K.matrix[h][ai][bj]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    global_dpd_->buf4_close(&K);

    // Close dpd files
    psio_->close(PSIF_LIBTRANS_DPD, 1);

    // Add Fock contribution
    for(int x = 0; x < nidpA; x++) {
	int a = idprowA[x];
	int i = idpcolA[x];
	int h = idpirrA[x];
        double value = FockA->get(h, a + occpiA[h], a + occpiA[h]) - FockA->get(h, i, i);
	Aorb->add(x, x, 2.0 * value);
    }
    if (print_ > 2) Aorb->print();

    // Solve the orb-resp equations
    pcg_conver = 0;// here 0 means successfull
    if (lineq == "CDGESV") Aorb->cdgesv(kappaA, pcg_conver);
    else if (lineq == "FLIN") {
         double det = 0.0;
         Aorb->lineq_flin(kappaA, &det);
         if (fabs(det) < DIIS_MIN_DET) {
             outfile->Printf( "Warning!!! MO Hessian matrix is near-singular\n");
             outfile->Printf( "Determinant is %6.3E\n", det);

             pcg_conver = 1;// here 1 means unsuccessful
         }
    }
    else if (lineq == "POPLE") Aorb->lineq_pople(kappaA, 6, cutoff);
    delete Aorb;

    // If LINEQ FAILED!
    if (pcg_conver != 0) {
       // Build kappa again
       for(int x = 0; x < nidpA; x++) {
	  int a = idprowA[x];
	  int i = idpcolA[x];
	  int h = idpirrA[x];
	  double value = FockA->get(h, a + occpiA[h], a + occpiA[h]) - FockA->get(h, i, i);
	  kappaA->set(x, -wogA->get(x) / (2.0*value));
       }

       outfile->Printf("\tWarning!!! MO Hessian matrix is near-singular, switching to MSD. \n");

    } // end if pcg_conver = 0

        // find biggest_kappa
	biggest_kappaA=0;
	for (int i=0; i<nidpA;i++) {
	    if (fabs(kappaA->get(i)) > biggest_kappaA) biggest_kappaA=fabs(kappaA->get(i));
	}

        // Scale
	if (biggest_kappaA > step_max) {
	    for (int i=0; i<nidpA;i++) kappaA->set(i, kappaA->get(i) *(step_max/biggest_kappaA));
	}

        // find biggest_kappa again
	if (biggest_kappaA > step_max)
	{
	  biggest_kappaA=0;
	  for (int i=0; i<nidpA;i++)
	  {
	      if (fabs(kappaA->get(i)) > biggest_kappaA)
	      {
		  biggest_kappaA = fabs(kappaA->get(i));
	      }
	  }
	}

        // norm
	rms_kappaA=0;
	rms_kappaA = kappaA->rms();

        // print
        if(print_ > 2) kappaA->print();

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
    // Open dpd files
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    dpdbuf4 K;

    // Sort some integrals
    // (OV|OV) -> (VO|VO)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                 ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , qpsr, ID("[V,O]"), ID("[V,O]"), "MO Ints (VO|VO)");
    global_dpd_->buf4_close(&K);

    // (ov|ov) -> (vo|vo)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                 ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , qpsr, ID("[v,o]"), ID("[v,o]"), "MO Ints (vo|vo)");
    global_dpd_->buf4_close(&K);

    // (AI|BJ) -> (AJ|BI)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[V,O]"),
                  ID("[V,O]"), ID("[V,O]"), 0, "MO Ints (VO|VO)");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , psrq, ID("[V,O]"), ID("[V,O]"), "MO Ints (AJ|BI)");
    global_dpd_->buf4_close(&K);

    // (ai|bj) -> (aj|bi)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[v,o]"),
                  ID("[v,o]"), ID("[v,o]"), 0, "MO Ints (vo|vo)");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , psrq, ID("[v,o]"), ID("[v,o]"), "MO Ints (aj|bi)");
    global_dpd_->buf4_close(&K);

    // <OV|OV> -> <VO|VO>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                 ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , qpsr, ID("[V,O]"), ID("[V,O]"), "MO Ints <VO|VO>");
    global_dpd_->buf4_close(&K);

    // <ov|ov> -> <vo|vo>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                 ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov>");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , qpsr, ID("[v,o]"), ID("[v,o]"), "MO Ints <vo|vo>");
    global_dpd_->buf4_close(&K);

    // (OV|ov) -> (VO|vo)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                 ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
    global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , qpsr, ID("[V,O]"), ID("[v,o]"), "MO Ints (VO|vo)");
    global_dpd_->buf4_close(&K);

    // Build the MO Hessian
    // Alpha-Alpha spin cae
    AorbAA = new Array2d("Alpha-Alpha MO Hessian Matrix", nidpA, nidpA);
    AorbAA->zero();
    // A(AI,BJ) = 4*(AI|BJ)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[V,O]"),
                  ID("[V,O]"), ID("[V,O]"), 0, "MO Ints (VO|VO)");
        int h =0;
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        //#pragma omp parallel for
        for(int ai = 0; ai < K.params->rowtot[h]; ++ai){
            int a = K.params->roworb[h][ai][0];
            int i = K.params->roworb[h][ai][1];
            for(int bj = 0; bj < K.params->coltot[h]; ++bj){
                int b = K.params->colorb[h][bj][0];
                int j = K.params->colorb[h][bj][1];
                AorbAA->set(ai, bj, 4.0 * K.matrix[h][ai][bj]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    global_dpd_->buf4_close(&K);

    // A(AI,BJ) -= 2*<AI|BJ>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[V,O]"),
                  ID("[V,O]"), ID("[V,O]"), 0, "MO Ints <VO|VO>");
        h = 0;
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        //#pragma omp parallel for
        for(int ai = 0; ai < K.params->rowtot[h]; ++ai){
            int a = K.params->roworb[h][ai][0];
            int i = K.params->roworb[h][ai][1];
            for(int bj = 0; bj < K.params->coltot[h]; ++bj){
                int b = K.params->colorb[h][bj][0];
                int j = K.params->colorb[h][bj][1];
                AorbAA->add(ai, bj, -2.0 * K.matrix[h][ai][bj]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    global_dpd_->buf4_close(&K);

    // A(AI,BJ) -= 2*(AJ|BI)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[V,O]"),
                  ID("[V,O]"), ID("[V,O]"), 0, "MO Ints (AJ|BI)");
        h = 0;
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        //#pragma omp parallel for
        for(int ai = 0; ai < K.params->rowtot[h]; ++ai){
            int a = K.params->roworb[h][ai][0];
            int i = K.params->roworb[h][ai][1];
            for(int bj = 0; bj < K.params->coltot[h]; ++bj){
                int b = K.params->colorb[h][bj][0];
                int j = K.params->colorb[h][bj][1];
                AorbAA->add(ai, bj, -2.0 * K.matrix[h][ai][bj]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    global_dpd_->buf4_close(&K);


    // Add Fock contribution
    for(int x = 0; x < nidpA; x++) {
	int a = idprowA[x];
	int i = idpcolA[x];
	int h = idpirrA[x];
        double value = FockA->get(h, a + occpiA[h], a + occpiA[h]) - FockA->get(h, i, i);
	AorbAA->add(x, x, 2.0 * value);
    }
    if (print_ > 2) AorbAA->print();

    // Build the UHF MO Hessian matrix
    Aorb = new Array2d("UHF MO Hessian Matrix", nidp_tot, nidp_tot);
    Aorb->zero();
    // AAAA part
    for (int x=0; x<nidpA;x++) {
      for (int y=0; y<nidpA;y++) {
	Aorb->set(x,y,AorbAA->get(x,y));
      }
    }
    delete AorbAA;

    // Beta-Beta spin cae
    AorbBB = new Array2d("Beta-Beta MO Hessian Matrix", nidpB, nidpB);
    AorbBB->zero();
    // A(ai,bj) = 4*(ai|bj)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[v,o]"),
                  ID("[v,o]"), ID("[v,o]"), 0, "MO Ints (vo|vo)");
        h =0;
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        //#pragma omp parallel for
        for(int ai = 0; ai < K.params->rowtot[h]; ++ai){
            int a = K.params->roworb[h][ai][0];
            int i = K.params->roworb[h][ai][1];
            for(int bj = 0; bj < K.params->coltot[h]; ++bj){
                int b = K.params->colorb[h][bj][0];
                int j = K.params->colorb[h][bj][1];
                AorbBB->set(ai, bj, 4.0 * K.matrix[h][ai][bj]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    global_dpd_->buf4_close(&K);

    // A(ai,bj) -= 2*<ai|bj>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[v,o]"),
                  ID("[v,o]"), ID("[v,o]"), 0, "MO Ints <vo|vo>");
        h = 0;
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        //#pragma omp parallel for
        for(int ai = 0; ai < K.params->rowtot[h]; ++ai){
            int a = K.params->roworb[h][ai][0];
            int i = K.params->roworb[h][ai][1];
            for(int bj = 0; bj < K.params->coltot[h]; ++bj){
                int b = K.params->colorb[h][bj][0];
                int j = K.params->colorb[h][bj][1];
                AorbBB->add(ai, bj, -2.0 * K.matrix[h][ai][bj]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    global_dpd_->buf4_close(&K);

    // A(ai,bj) -= 2*(aj|bi)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[v,o]"),
                  ID("[v,o]"), ID("[v,o]"), 0, "MO Ints (aj|bi)");
        h = 0;
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        //#pragma omp parallel for
        for(int ai = 0; ai < K.params->rowtot[h]; ++ai){
            int a = K.params->roworb[h][ai][0];
            int i = K.params->roworb[h][ai][1];
            for(int bj = 0; bj < K.params->coltot[h]; ++bj){
                int b = K.params->colorb[h][bj][0];
                int j = K.params->colorb[h][bj][1];
                AorbBB->add(ai, bj, -2.0 * K.matrix[h][ai][bj]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    global_dpd_->buf4_close(&K);

    // Add Fock contribution
    for(int x = 0; x < nidpB; x++) {
	int a = idprowB[x];
	int i = idpcolB[x];
	int h = idpirrB[x];
        double value = FockB->get(h, a + occpiB[h], a + occpiB[h]) - FockB->get(h, i, i);
	AorbBB->add(x, x, 2.0 * value);
    }
    if (print_ > 2) AorbBB->print();

    // Build the UHF MO Hessian matrix
    // BBBB part
    for (int x=0; x<nidpB;x++) {
      for (int y=0; y<nidpB;y++) {
	Aorb->set(x+nidpA,y+nidpA,AorbBB->get(x,y));
      }
    }
    delete AorbBB;

    // Alpha-Beta spin cae
    AorbAB = new Array2d("Alpha-Beta MO Hessian Matrix", nidpA, nidpB);
    AorbAB->zero();
    // A(AI,bj) = 4*(AI|bj)
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[v,o]"),
                  ID("[V,O]"), ID("[v,o]"), 0, "MO Ints (VO|vo)");
        h = 0;
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        //#pragma omp parallel for
        for(int ai = 0; ai < K.params->rowtot[h]; ++ai){
            int a = K.params->roworb[h][ai][0];
            int i = K.params->roworb[h][ai][1];
            for(int bj = 0; bj < K.params->coltot[h]; ++bj){
                int b = K.params->colorb[h][bj][0];
                int j = K.params->colorb[h][bj][1];
                AorbAB->set(ai, bj, 4.0 * K.matrix[h][ai][bj]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    global_dpd_->buf4_close(&K);
    if (print_ > 2) AorbAB->print();

    // Close dpd files
    psio_->close(PSIF_LIBTRANS_DPD, 1);


    // Build the UHF MO Hessian matrix
    // AABB part
    for (int x=0; x<nidpA;x++) {
      for (int y=0; y<nidpB;y++) {
	Aorb->set(x,y+nidpA,AorbAB->get(x,y));
      }
    }

    // BBAA part
    for (int x=0; x<nidpB;x++) {
      for (int y=0; y<nidpA;y++) {
	Aorb->set(x+nidpA,y,AorbAB->get(y,x));
      }
    }
    delete AorbAB;

    // Print
    if (print_ > 2) Aorb->print();

    // Build total kappa
    kappa->zero();
    for (int x=0; x<nidpA;x++) kappa->set(x, -wogA->get(x));
    for (int x=0; x<nidpB;x++) kappa->set(x + nidpA, -wogB->get(x));

    // Solve the orb-resp equations
    pcg_conver = 0;// here 0 means successfull
    if (lineq == "CDGESV") Aorb->cdgesv(kappa, pcg_conver);
    else if (lineq == "FLIN") {
         double det = 0.0;
         Aorb->lineq_flin(kappa, &det);
         if (fabs(det) < DIIS_MIN_DET) {
         //if (fabs(det) < 1e-2) {
             outfile->Printf( "Warning!!! MO Hessian matrix is near-singular\n");
             outfile->Printf( "Determinant is %6.3E\n", det);

             pcg_conver = 1;// here 1 means unsuccessful
         }
    }
    else if (lineq == "POPLE") Aorb->lineq_pople(kappa, 6, cutoff);
    delete Aorb;

    // Build kappaA and kappaB
    //kappa->print();
    kappaA->zero();
    kappaB->zero();
    for (int x=0; x<nidpA;x++) kappaA->set(x, kappa->get(x));
    for (int x=0; x<nidpB;x++) kappaB->set(x, kappa->get(x + nidpA));

    // If LINEQ FAILED!
    if (pcg_conver != 0) {
       // Build kappa again
       for(int x = 0; x < nidpA; x++) {
	  int a = idprowA[x];
	  int i = idpcolA[x];
	  int h = idpirrA[x];
	  double value = FockA->get(h, a + occpiA[h], a + occpiA[h]) - FockA->get(h, i, i);
	  kappaA->set(x, -wogA->get(x) / (2.0*value));
       }

	// beta
	for(int x = 0; x < nidpB; x++) {
	  int a = idprowB[x];
	  int i = idpcolB[x];
	  int h = idpirrB[x];
	  double value = FockB->get(h, a + occpiB[h], a + occpiB[h]) - FockB->get(h, i, i);
	  kappaB->set(x, -wogB->get(x) / (2.0*value));
	}
       outfile->Printf("\tWarning!!! MO Hessian matrix is near-singular, switching to MSD. \n");

    } // end if pcg_conver = 0


        // find biggest_kappa
	biggest_kappaA=0;
	for (int i=0; i<nidpA;i++) {
	    if (fabs(kappaA->get(i)) > biggest_kappaA) biggest_kappaA=fabs(kappaA->get(i));
	}

	biggest_kappaB=0;
	for (int i=0; i<nidpB;i++){
	    if (fabs(kappaB->get(i)) > biggest_kappaB) biggest_kappaB=fabs(kappaB->get(i));
	}

        // Scale
	if (biggest_kappaA > step_max) {
	    for (int i=0; i<nidpA;i++) kappaA->set(i, kappaA->get(i) *(step_max/biggest_kappaA));
	}

	if (biggest_kappaB > step_max) {
	    for (int i=0; i<nidpB;i++) kappaB->set(i, kappaB->get(i) *(step_max/biggest_kappaB));
	}

        // find biggest_kappa again
	if (biggest_kappaA > step_max)
	{
	  biggest_kappaA=0;
	  for (int i=0; i<nidpA;i++)
	  {
	      if (fabs(kappaA->get(i)) > biggest_kappaA)
	      {
		  biggest_kappaA = fabs(kappaA->get(i));
	      }
	  }
	}

	if (biggest_kappaB > step_max)
	{
	  biggest_kappaB=0;
	  for (int i=0; i<nidpB;i++)
	  {
	      if (fabs(kappaB->get(i)) > biggest_kappaB)
	      {
		  biggest_kappaB=fabs(kappaB->get(i));
	      }
	  }
	}

        // norm
	rms_kappaA=0;
	rms_kappaB=0;
	rms_kappaA = kappaA->rms();
	rms_kappaB = kappaB->rms();

        // print
        if(print_ > 2){
          kappaA->print();
          kappaB->print();
        }

}// end if (reference_ == "UNRESTRICTED")
 //outfile->Printf("\n kappa_orb_resp done. \n");
}// end kappa_orb_resp
}} // End Namespaces
