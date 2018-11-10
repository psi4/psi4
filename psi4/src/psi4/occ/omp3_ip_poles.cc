/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

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
#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/libtrans/integraltransform.h"

/** Required libmints includes */

#include "psi4/libmints/factory.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/molecule.h"
#include "defines.h"
#include "occwave.h"


namespace psi{ namespace occwave{

void OCCWave::omp3_ip_poles()
{

//outfile->Printf("\n omp3_ip_poles is starting... \n");
//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference_ == "RESTRICTED") {
     // Memory allocation
     auto eOccOrbA = std::make_shared<Vector>("eOccOrbA", nirrep_, occpiA);
     eOccOrbA->zero();

     dpdbuf4 K, T, D;

     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);


    // Build denominators D_jk^ia = E_j + E_k - E_i - E_a
    auto *aOccEvals = new double [nacooA];
    auto *aVirEvals = new double [nacvoA];

    // Pick out the diagonal elements of the Fock matrix, making sure that they are in the order
    // used by the DPD library, i.e. starting from zero for each space and ordering by irrep

    int aOccCount = 0, aVirCount = 0;

    //Diagonal elements of the Fock matrix
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < aoccpiA[h]; ++i) aOccEvals[aOccCount++] = FockA->get(h, i + frzcpi_[h], i + frzcpi_[h]);
        for(int a = 0; a < avirtpiA[h]; ++a) aVirEvals[aVirCount++] = FockA->get(h, occpiA[h] + a, occpiA[h] + a);
    }

    // Build denominators
    // The alpha-alpha spin case: D_IA^JK
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 0, "D <OV|OO>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int a = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int j = D.params->colorb[h][col][0];
                int k = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(aOccEvals[i] + aVirEvals[a] - aOccEvals[j] - aOccEvals[k]);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&D, h);
        global_dpd_->buf4_mat_irrep_close(&D, h);
    }
    global_dpd_->buf4_close(&D);

    delete [] aOccEvals;
    delete [] aVirEvals;


    // NOTE:! The followings are MP2 amplitudes, they should be MP3 amplitudes.
    // Build T_IA^JK
    // T_IA^JK = <IA||JK> / D_IA^JK
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO|OV>");
    global_dpd_->buf4_sort(&K, PSIF_OCC_DPD , rspq, ID("[O,V]"), ID("[O,O]"), "T2 <OV|OO>");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 0, "D <OV|OO>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 0, "T2 <OV|OO>");
    global_dpd_->buf4_dirprd(&D, &T);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&T);

    // Build Beta occ orbital energy
    // e_I = F_II : reference contribution
     for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < occpiA[h]; ++i){
	    eOccOrbA->set(h, i, FockA->get(h, i, i));
         }
     }

    // e_I =  \sum{J,A,B} (2*T_IJ^AB - T_JI^AB)  <IJ|AB>
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                 ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                 ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
	global_dpd_->buf4_mat_irrep_init(&T, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&T, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            int i = K.params->roworb[h][ij][0];
            int j = K.params->roworb[h][ij][1];
            int ji = K.params->rowidx[j][i];
            int hi = K.params->psym[i];
            int ii = i - occ_offA[hi];
            for(int ab = 0; ab < K.params->coltot[h]; ++ab){
                double value = ( (2.0 * T.matrix[h][ij][ab]) - T.matrix[h][ji][ab]) * K.matrix[h][ij][ab];
		eOccOrbA->add(hi, ii, value);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
        global_dpd_->buf4_mat_irrep_close(&T, h);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&T);


    // e_I = \sum{J,K,A} (2*T_IA^JK - T_IA^KJ)  <JK|IA>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO|OV>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 0, "T2 <OV|OO>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
	global_dpd_->buf4_mat_irrep_init(&T, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&T, h);
        for(int jk = 0; jk < K.params->rowtot[h]; ++jk){
            int j = K.params->roworb[h][jk][0];
            int k = K.params->roworb[h][jk][1];
            int kj = K.params->rowidx[k][j];
            for(int ia = 0; ia < K.params->coltot[h]; ++ia){
                int i = K.params->colorb[h][ia][0];
                int hi = K.params->rsym[i];
                int ii = i - occ_offA[hi];
                double value = ( (2.0 * T.matrix[h][ia][jk]) - T.matrix[h][ia][kj]) * K.matrix[h][jk][ia];
		eOccOrbA->add(hi, ii, value);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
        global_dpd_->buf4_mat_irrep_close(&T, h);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&T);

    psio_->close(PSIF_LIBTRANS_DPD, 1);
    psio_->close(PSIF_OCC_DPD, 1);


    // Sort orbitals
    Array1d *evals_A = new Array1d("Alpha OCC ORB C1", nooA);
    Array1i *irrep_A = new Array1i("IrrepA", nooA);
    evals_A->zero();
    irrep_A->zero();


    // Sort Alpha spin-case
    int count = 0;
    for (int h = 0; h < nirrep_; ++h){
	 for (int i = 0; i < occpiA[h]; ++i){
              evals_A->set(count, eOccOrbA->get(h, i));
              irrep_A->set(count, h);
              count++;
          }
    }

    //outfile->Printf("\tI am here. \n");

    for (int i = 0; i < nooA; ++i) {
         for(int j = nooA-1; j > i; --j) {
             if (evals_A->get(j-1) > evals_A->get(j)) {
                 double dum = evals_A->get(j-1);
                 evals_A->set(j-1, evals_A->get(j));
                 evals_A->set(j, dum);
                 int dum2 = irrep_A->get(j-1);
                 irrep_A->set(j-1, irrep_A->get(j));
                 irrep_A->set(j, dum2);
             }
         }
    }


    // Print occupied orbital energies
    if (mo_optimized == 1) outfile->Printf("\n\tOMP3 Occupied Orbital Energies (a.u.) \n");
    else if (mo_optimized == 0) outfile->Printf("\n\tMP3 Occupied Orbital Energies (a.u.) \n");
    outfile->Printf("\t----------------------------------------------- \n");


    Molecule& mol = *reference_wavefunction_->molecule().get();
    CharacterTable ct = mol.point_group()->char_table();
    std::string pgroup = mol.point_group()->symbol();

    // print alpha occ orb energy
    outfile->Printf( "\tAlpha occupied orbitals\n");
    count = 1;
    for (int i = 0; i < nooA; ++i) {
         int h = irrep_A->get(i);
	 outfile->Printf("\t%3d (%-3s) %20.10f \n",count,ct.gamma(h).symbol(),evals_A->get(i));

	 count++;
    }

       eOccOrbA.reset();
       delete evals_A;
       delete irrep_A;

}// end if (reference_ == "RESTRICTED")


//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference_ == "UNRESTRICTED") {

     // Memory allocation
     auto eOccOrbA = std::make_shared<Vector>("eOccOrbA", nirrep_, occpiA);
     auto eOccOrbB = std::make_shared<Vector>("eOccOrbB", nirrep_, occpiB);
     eOccOrbA->zero();
     eOccOrbB->zero();

     dpdbuf4 K, T, D;

     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);


    // Build denominators D_jk^ia = E_j + E_k - E_i - E_a
    auto *aOccEvals = new double [nacooA];
    auto *bOccEvals = new double [nacooB];
    auto *aVirEvals = new double [nacvoA];
    auto *bVirEvals = new double [nacvoB];

    // Pick out the diagonal elements of the Fock matrix, making sure that they are in the order
    // used by the DPD library, i.e. starting from zero for each space and ordering by irrep

    int aOccCount = 0, bOccCount = 0, aVirCount = 0, bVirCount = 0;

    //Diagonal elements of the Fock matrix
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < aoccpiA[h]; ++i) aOccEvals[aOccCount++] = FockA->get(h, i + frzcpi_[h], i + frzcpi_[h]);
	for(int i = 0; i < aoccpiB[h]; ++i) bOccEvals[bOccCount++] = FockB->get(h, i + frzcpi_[h], i + frzcpi_[h]);
        for(int a = 0; a < avirtpiA[h]; ++a) aVirEvals[aVirCount++] = FockA->get(h, occpiA[h] + a, occpiA[h] + a);
	for(int a = 0; a < avirtpiB[h]; ++a) bVirEvals[bVirCount++] = FockB->get(h, occpiB[h] + a, occpiB[h] + a);
    }

    // Build denominators
    // The alpha-alpha spin case: D_IA^JK
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 0, "D <OV|OO>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int a = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int j = D.params->colorb[h][col][0];
                int k = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(aOccEvals[i] + aVirEvals[a] - aOccEvals[j] - aOccEvals[k]);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&D, h);
        global_dpd_->buf4_mat_irrep_close(&D, h);
    }
    global_dpd_->buf4_close(&D);

    // The beta-beta spin case: D_ia^jk
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o,o]"), 0, "D <ov|oo>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int a = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int j = D.params->colorb[h][col][0];
                int k = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(bOccEvals[i] + bVirEvals[a] - bOccEvals[j] - bOccEvals[k]);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&D, h);
        global_dpd_->buf4_mat_irrep_close(&D, h);
    }
    global_dpd_->buf4_close(&D);

    // The alpha-beta spin case: D_Ia^Jk
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "D <Ov|Oo>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int a = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int j = D.params->colorb[h][col][0];
                int k = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(aOccEvals[i] + bVirEvals[a] - aOccEvals[j] - bOccEvals[k]);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&D, h);
        global_dpd_->buf4_mat_irrep_close(&D, h);
    }
    global_dpd_->buf4_close(&D);


    // The beta-alpha spin case: D_Ai^Jk
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "D <Vo|Oo>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int a = D.params->roworb[h][row][0];
            int i = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int j = D.params->colorb[h][col][0];
                int k = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(aVirEvals[a] + bOccEvals[i] - aOccEvals[j] - bOccEvals[k]);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&D, h);
        global_dpd_->buf4_mat_irrep_close(&D, h);
    }
    global_dpd_->buf4_close(&D);


    delete [] aOccEvals;
    delete [] bOccEvals;
    delete [] aVirEvals;
    delete [] bVirEvals;


    // Build T_IA^JK
    // T_IA^JK = <IA||JK> / D_IA^JK
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO||OV>");
    global_dpd_->buf4_sort(&K, PSIF_OCC_DPD , rspq, ID("[O,V]"), ID("[O,O]"), "T2 <OV|OO>");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 0, "D <OV|OO>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 0, "T2 <OV|OO>");
    global_dpd_->buf4_dirprd(&D, &T);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&T);


    // T_ia^jk = <ia||jk> / D_ia^jk
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo||ov>");
    global_dpd_->buf4_sort(&K, PSIF_OCC_DPD , rspq, ID("[o,v]"), ID("[o,o]"), "T2 <ov|oo>");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o,o]"), 0, "D <ov|oo>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o,o]"), 0, "T2 <ov|oo>");
    global_dpd_->buf4_dirprd(&D, &T);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&T);


    // T_Ia^Jk = <Ia|Jk> / D_Ia^Jk
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "MO Ints <Oo|Ov>");
    global_dpd_->buf4_sort(&K, PSIF_OCC_DPD , rspq, ID("[O,v]"), ID("[O,o]"), "T2 <Ov|Oo>");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "D <Ov|Oo>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "T2 <Ov|Oo>");
    global_dpd_->buf4_dirprd(&D, &T);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&T);


    // T_Ai^Jk = <Ai|Jk> /D_Ai^Jk
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "MO Ints <Oo|Vo>");
    global_dpd_->buf4_sort(&K, PSIF_OCC_DPD , rspq, ID("[V,o]"), ID("[O,o]"), "T2 <Vo|Oo>");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "D <Vo|Oo>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "T2 <Vo|Oo>");
    global_dpd_->buf4_dirprd(&D, &T);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&T);


    // Build Beta occ orbital energy
    // e_I = F_II : reference contribution
     for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < occpiA[h]; ++i){
	    eOccOrbA->set(h, i, FockA->get(h, i, i));
         }
     }

    // e_I = 1/2 * \sum{J,A,B} T_IJ^AB  <IJ||AB>
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                 ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                 ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
	global_dpd_->buf4_mat_irrep_init(&T, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&T, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            int i = K.params->roworb[h][ij][0];
            int hi = K.params->psym[i];
            int ii = i - occ_offA[hi];
            for(int ab = 0; ab < K.params->coltot[h]; ++ab){
		eOccOrbA->add(hi, ii, 0.5 * K.matrix[h][ij][ab] * T.matrix[h][ij][ab]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
        global_dpd_->buf4_mat_irrep_close(&T, h);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&T);

    // e_I = \sum{j,A,b} T_Ij^Ab  <Ij|Ab>
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
	global_dpd_->buf4_mat_irrep_init(&T, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&T, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            int i = K.params->roworb[h][ij][0];
            int hi = K.params->psym[i];
            int ii = i - occ_offA[hi];
            for(int ab = 0; ab < K.params->coltot[h]; ++ab){
		eOccOrbA->add(hi, ii, K.matrix[h][ij][ab] * T.matrix[h][ij][ab]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
        global_dpd_->buf4_mat_irrep_close(&T, h);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&T);

    // e_I = 1/2 * \sum{J,K,A} T_IA^JK  <JK||IA>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO||OV>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 0, "T2 <OV|OO>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
	global_dpd_->buf4_mat_irrep_init(&T, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&T, h);
        for(int jk = 0; jk < K.params->rowtot[h]; ++jk){
            for(int ia = 0; ia < K.params->coltot[h]; ++ia){
                int i = K.params->colorb[h][ia][0];
                int hi = K.params->rsym[i];
                int ii = i - occ_offA[hi];
		eOccOrbA->add(hi, ii, 0.5 * K.matrix[h][jk][ia] * T.matrix[h][ia][jk]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
        global_dpd_->buf4_mat_irrep_close(&T, h);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&T);

    // e_I = \sum{J,k,a} T_Ia^Jk  <Jk|Ia>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "MO Ints <Oo|Ov>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "T2 <Ov|Oo>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
	global_dpd_->buf4_mat_irrep_init(&T, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&T, h);
        for(int jk = 0; jk < K.params->rowtot[h]; ++jk){
            for(int ia = 0; ia < K.params->coltot[h]; ++ia){
                int i = K.params->colorb[h][ia][0];
                int hi = K.params->rsym[i];
                int ii = i - occ_offA[hi];
		eOccOrbA->add(hi, ii, K.matrix[h][jk][ia] * T.matrix[h][ia][jk]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
        global_dpd_->buf4_mat_irrep_close(&T, h);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&T);


    // Build Beta occ orbital energy
    // e_i = F_ii : reference contribution
     for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < occpiB[h]; ++i){
	    eOccOrbB->set(h, i, FockB->get(h, i , i));
         }
     }

    // e_i = 1/2 * \sum{j,a,b} T_ij^ab  <ij||ab>
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                 ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                 ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
	global_dpd_->buf4_mat_irrep_init(&T, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&T, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            int i = K.params->roworb[h][ij][0];
            int hi = K.params->psym[i];
            int ii = i - occ_offB[hi];
            for(int ab = 0; ab < K.params->coltot[h]; ++ab){
		eOccOrbB->add(hi, ii, 0.5 * K.matrix[h][ij][ab] * T.matrix[h][ij][ab]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
        global_dpd_->buf4_mat_irrep_close(&T, h);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&T);

    // e_i = \sum{J,A,b} T_Ji^Ab  <Ji|Ab>
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
	global_dpd_->buf4_mat_irrep_init(&T, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&T, h);
        for(int ji = 0; ji < K.params->rowtot[h]; ++ji){
            int i = K.params->roworb[h][ji][1];
            int hi = K.params->qsym[i];
            int ii = i - occ_offB[hi];
            for(int ab = 0; ab < K.params->coltot[h]; ++ab){
		eOccOrbB->add(hi, ii, K.matrix[h][ji][ab] * T.matrix[h][ji][ab]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
        global_dpd_->buf4_mat_irrep_close(&T, h);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&T);

    // e_i = 1/2 * \sum{j,k,a} T_ia^jk  <jk||ia>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo||ov>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o,o]"), 0, "T2 <ov|oo>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
	global_dpd_->buf4_mat_irrep_init(&T, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&T, h);
        for(int jk = 0; jk < K.params->rowtot[h]; ++jk){
            for(int ia = 0; ia < K.params->coltot[h]; ++ia){
                int i = K.params->colorb[h][ia][0];
                int hi = K.params->rsym[i];
                int ii = i - occ_offB[hi];
		eOccOrbB->add(hi, ii, 0.5 * K.matrix[h][jk][ia] * T.matrix[h][ia][jk]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
        global_dpd_->buf4_mat_irrep_close(&T, h);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&T);

    // e_i = \sum{J,k,A} T_Ai^Jk  <Jk|Ai>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "MO Ints <Oo|Vo>");
    global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "T2 <Vo|Oo>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
	global_dpd_->buf4_mat_irrep_init(&T, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&T, h);
        for(int jk = 0; jk < K.params->rowtot[h]; ++jk){
            for(int ai = 0; ai < K.params->coltot[h]; ++ai){
                int i = K.params->colorb[h][ai][1];
                int hi = K.params->ssym[i];
                int ii = i - occ_offB[hi];
		eOccOrbB->add(hi, ii, K.matrix[h][jk][ai] * T.matrix[h][ai][jk]);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
        global_dpd_->buf4_mat_irrep_close(&T, h);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&T);

    psio_->close(PSIF_LIBTRANS_DPD, 1);
    psio_->close(PSIF_OCC_DPD, 1);


    // Sort orbitals
    Array1d *evals_A = new Array1d("Alpha OCC ORB C1", nooA);
    Array1d *evals_B = new Array1d("Beta OCC ORB C1", nooB);
    Array1i *irrep_A = new Array1i("IrrepA", nooA);
    Array1i *irrep_B = new Array1i("IrrepA", nooB);
    evals_A->zero();
    evals_B->zero();
    irrep_A->zero();
    irrep_B->zero();


    // Sort Alpha spin-case
    int count = 0;
    for (int h = 0; h < nirrep_; ++h){
	 for (int i = 0; i < occpiA[h]; ++i){
              evals_A->set(count, eOccOrbA->get(h, i));
              irrep_A->set(count, h);
              count++;
          }
    }

    //outfile->Printf("\tI am here. \n");

    for (int i = 0; i < nooA; ++i) {
         for(int j = nooA-1; j > i; --j) {
             if (evals_A->get(j-1) > evals_A->get(j)) {
                 double dum = evals_A->get(j-1);
                 evals_A->set(j-1, evals_A->get(j));
                 evals_A->set(j, dum);
                 int dum2 = irrep_A->get(j-1);
                 irrep_A->set(j-1, irrep_A->get(j));
                 irrep_A->set(j, dum2);
             }
         }
    }

    // Sort beta spin-case
    count = 0;
    for (int h = 0; h < nirrep_; ++h){
	 for (int i = 0; i < occpiB[h]; ++i){
              evals_B->set(count, eOccOrbB->get(h, i));
              irrep_B->set(count, h);
              count++;
          }
    }
    for (int i = 0; i < nooB; ++i) {
         for(int j = nooB-1; j > i; --j) {
             if (evals_B->get(j-1) > evals_B->get(j)) {
                 double dum = evals_B->get(j-1);
                 evals_B->set(j-1, evals_B->get(j));
                 evals_B->set(j, dum);
                 int dum2 = irrep_B->get(j-1);
                 irrep_B->set(j-1, irrep_B->get(j));
                 irrep_B->set(j, dum2);
             }
         }
    }


    // Print occupied orbital energies
    if (mo_optimized == 1) outfile->Printf("\n\tOMP3 Occupied Orbital Energies (a.u.) \n");
    else if (mo_optimized == 0) outfile->Printf("\n\tMP3 Occupied Orbital Energies (a.u.) \n");
    outfile->Printf("\t----------------------------------------------- \n");


    Molecule& mol = *reference_wavefunction_->molecule().get();
    CharacterTable ct = mol.point_group()->char_table();
    std::string pgroup = mol.point_group()->symbol();

    // print alpha occ orb energy
    outfile->Printf( "\tAlpha occupied orbitals\n");
    count = 1;
    for (int i = 0; i < nooA; ++i) {
         int h = irrep_A->get(i);
	 outfile->Printf("\t%3d (%-3s) %20.10f \n",count,ct.gamma(h).symbol(),evals_A->get(i));

	 count++;
    }

    // print beta occ orb energy
    outfile->Printf( "\n\tBeta occupied orbitals\n");
    count = 1;
    for (int i = 0; i < nooB; ++i) {
         int h = irrep_B->get(i);
	 outfile->Printf("\t%3d (%-3s) %20.10f \n",count,ct.gamma(h).symbol(),evals_B->get(i));

	 count++;
    }

       eOccOrbA.reset();
       eOccOrbB.reset();
       delete evals_A;
       delete evals_B;
       delete irrep_A;
       delete irrep_B;

}// end if (reference_ == "UNRESTRICTED")
//outfile->Printf("\n omp3_ip_poles is done. \n");
} // end omp3_ip_poles
}} // End Namespaces
