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
#include "defines.h"
#include "occwave.h"


using namespace std;


namespace psi{ namespace occwave{

void OCCWave::fock_beta()
{

  // F(pq) = h(pq) + \sum_{m} <pm||qm> in spin-orbital form.

  // F(pq) = h(pq) + \sum_{m} <pm||qm> + \sum_{M} <Mp|Mq> in spin-adapted form.

/************************************************************************************************/
/*********************************** Build Fij **************************************************/
/************************************************************************************************/
  // F(ij) = h(ij) + \sum_{m} <im||jm> + \sum_{M} <Mi|Mj>

  dpdfile2 F;
  dpdbuf4 K;

  /* Prepare the alpha occ-occ Fock matrix files */
  global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "Fock <o|o>");
  global_dpd_->file2_mat_init(&F);

  /* Two-electron contributions */

  /* Prepare the <oo||oo> integral buffers */
   global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo||oo>");

  for(int h=0; h < nirrep_; h++) {

      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirrep_; Gi++) {
	  int Gm=h^Gi;

          /* Loop over the orbitals of the target */
          for(int i=0; i < occpiB[Gi]; i++) {
              int I = occ_offB[Gi] + i;
              for(int j=0; j < occpiB[Gi]; j++) {
                  int J = occ_offB[Gi] + j;
                  for(int m=0; m < occpiB[Gm]; m++) {
                      int M = occ_offB[Gm] + m;

                      int IM = K.params->rowidx[I][M];
                      int JM = K.params->colidx[J][M];

                      F.matrix[Gi][i][j] += K.matrix[h][IM][JM];

                    }
                }
            }
    }
    global_dpd_->buf4_mat_irrep_close(&K, h);
  }
  global_dpd_->buf4_close(&K);


  /* Prepare the <Oo|Oo> integral buffers */
   global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");

  for(int h=0; h < nirrep_; h++) {

      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirrep_; Gi++) {
	  int Gm=h^Gi;

          /* Loop over the orbitals of the target */
          for(int i=0; i < occpiB[Gi]; i++) {
              int I = occ_offB[Gi] + i;
              for(int j=0; j < occpiB[Gi]; j++) {
                  int J = occ_offB[Gi] + j;
                  for(int m=0; m < occpiA[Gm]; m++) {
                      int M = occ_offA[Gm] + m;

                      int MI = K.params->rowidx[M][I];
                      int MJ = K.params->colidx[M][J];

                      F.matrix[Gi][i][j] += K.matrix[h][MI][MJ];

                    }
                }
            }
    }
    global_dpd_->buf4_mat_irrep_close(&K, h);
  }
  global_dpd_->buf4_close(&K);

  /* Close the Integral buffers */
  global_dpd_->file2_mat_wrt(&F);
  global_dpd_->file2_mat_close(&F);
  global_dpd_->file2_close(&F);

/************************************************************************************************/
/*********************************** Build Fab **************************************************/
/************************************************************************************************/
  // F(ab) = h(ab) + \sum_{m} <am|bm> + \sum_{M} <Ma|Mb>

  /* Prepare the alpha vir-vir Fock matrix files */
  global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "Fock <v|v>");
  global_dpd_->file2_mat_init(&F);

  /* Two-electron contributions */

  /* Prepare the <ov||ov> integral buffers */
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov||ov>");

  for(int h=0; h < nirrep_; h++) {

      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Ga=0; Ga < nirrep_; Ga++) {
	  int Gb = Ga;
	  int Gm = h^Ga;

	  /* Loop over orbitals of the target */
	  for(int a=0; a < virtpiB[Ga]; a++) {
	      int A = vir_offB[Ga] + a;
	      for(int b=0; b < virtpiB[Gb]; b++) {
		  int B = vir_offB[Gb] + b;

		  for(int m=0; m < occpiB[Gm]; m++) {
		      int M = occ_offB[Gm] + m;

		      int MA = K.params->rowidx[M][A];
		      int MB = K.params->colidx[M][B];

		      F.matrix[Ga][a][b] += K.matrix[h][MA][MB];

		    }
		}
	    }
    }
    global_dpd_->buf4_mat_irrep_close(&K, h);
  }
  global_dpd_->buf4_close(&K);


   /* Prepare the <Ov|Ov> integral buffers */
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>");

  for(int h=0; h < nirrep_; h++) {

      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Ga=0; Ga < nirrep_; Ga++) {
	  int Gb = Ga;
	  int Gm = h^Ga;

	  /* Loop over orbitals of the target */
	  for(int a=0; a < virtpiB[Ga]; a++) {
	      int A = vir_offB[Ga] + a;
	      for(int b=0; b < virtpiB[Gb]; b++) {
		  int B = vir_offB[Gb] + b;

		  for(int m=0; m < occpiA[Gm]; m++) {
		      int M = occ_offA[Gm] + m;

		      int MA = K.params->rowidx[M][A];
		      int MB = K.params->colidx[M][B];

		      F.matrix[Ga][a][b] += K.matrix[h][MA][MB];

		    }
		}
	    }
    }
    global_dpd_->buf4_mat_irrep_close(&K, h);
  }
  global_dpd_->buf4_close(&K);


  /* Close the buffers */
  global_dpd_->file2_mat_wrt(&F);
  global_dpd_->file2_mat_close(&F);
  global_dpd_->file2_close(&F);

/************************************************************************************************/
/*********************************** Build Fia **************************************************/
/************************************************************************************************/
  // F(ia) = h(ia) + \sum_{m} <im||am> + \sum_{M} <Mi|Ma>

  /* Prepare the alpha and beta occ-vir FockB matrix files */
  global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "Fock <o|v>");
  global_dpd_->file2_mat_init(&F);


  /* Two-electron contributions */

  /* Prepare the <oo||ov> integral buffers */
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo||ov>");

  for(int h=0; h < nirrep_; h++) {

      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirrep_; Gi++) {
	  int Ga = Gi; int Gm = h^Gi;

	  /* Loop over orbitals of the target */
	  for(int i=0; i < occpiB[Gi]; i++) {
	      int I = occ_offB[Gi] + i;
	      for(int a=0; a < virtpiB[Ga]; a++) {
		  int A = vir_offB[Ga] + a;

		  for(int m=0; m < occpiB[Gm]; m++) {
		      int M = occ_offB[Gm] + m;

		      int MI = K.params->rowidx[M][I];
		      int MA = K.params->colidx[M][A];

		      F.matrix[Gi][i][a] += K.matrix[h][MI][MA];

		    }
		}
	    }
    }
    global_dpd_->buf4_mat_irrep_close(&K, h);
  }
  global_dpd_->buf4_close(&K);


   /* Prepare the <Oo|Ov> integral buffers */
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "MO Ints <Oo|Ov>");

  for(int h=0; h < nirrep_; h++) {

      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirrep_; Gi++) {
	  int Ga = Gi; int Gm = h^Gi;

	  /* Loop over orbitals of the target */
	  for(int i=0; i < occpiB[Gi]; i++) {
	      int I = occ_offB[Gi] + i;
	      for(int a=0; a < virtpiB[Ga]; a++) {
		  int A = vir_offB[Ga] + a;

		  for(int m=0; m < occpiA[Gm]; m++) {
		      int M = occ_offA[Gm] + m;

		      int MI = K.params->rowidx[M][I];
		      int MA = K.params->colidx[M][A];

		      F.matrix[Gi][i][a] += K.matrix[h][MI][MA];

		    }
		}
	    }
    }
    global_dpd_->buf4_mat_irrep_close(&K, h);
  }
  global_dpd_->buf4_close(&K);

  /* Close the buffers */
  global_dpd_->file2_mat_wrt(&F);
  global_dpd_->file2_mat_close(&F);
  global_dpd_->file2_close(&F);

/************************************************************************************************/
/*********************************** Set Fock ***************************************************/
/************************************************************************************************/

    // <o|o> block
    global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "Fock <o|o>");
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < occpiB[h]; ++j){
		FockB->set(h, i, j, F.matrix[h][i][j]);
            }
        }
    }
    global_dpd_->file2_close(&F);

    // <v|v> block
    global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "Fock <v|v>");
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < virtpiB[h]; ++i){
            for(int j = 0 ; j < virtpiB[h]; ++j){
               FockB->set(h, i + occpiB[h], j + occpiB[h], F.matrix[h][i][j]);
            }
        }
    }
    global_dpd_->file2_close(&F);

    // <o|v> block
    global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "Fock <o|v>");
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < virtpiB[h]; ++j){
               FockB->set(h, i, j + occpiB[h], F.matrix[h][i][j]);
	       FockB->set(h, j + occpiB[h], i , F.matrix[h][i][j]);
            }
        }
    }
    global_dpd_->file2_close(&F);

     //1e-contr.
    FockB->add(HmoB);


/************************************************************************************************/
/*********************************** Print FockB *************************************************/
/************************************************************************************************/
	if (print_ > 1) FockB->print();

/************************************************************************************************/
/************************************************************************************************/

}
}} // End Namespaces
