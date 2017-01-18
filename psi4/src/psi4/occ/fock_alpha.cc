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

void OCCWave::fock_alpha()
{

//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference_ == "RESTRICTED") {

/************************************************************************************************/
/*********************************** Build Fij **************************************************/
/************************************************************************************************/
  dpdfile2 F;
  dpdbuf4 K;

  /* Prepare the alpha and beta occ-occ Fock matrix files */
  global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "Fock <O|O>");
  global_dpd_->file2_mat_init(&F);

  /* Two-electron contributions */

  /* Prepare the <OO|OO> integral buffers */
   global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,O]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,O]"), 0, "MO Ints <OO|OO>");
  // part-1
  for(int h=0; h < nirrep_; h++) {

      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirrep_; Gi++) {
	  int Gm=h^Gi;

          /* Loop over the orbitals of the target */
          for(int i=0; i < occpiA[Gi]; i++) {
              int I = occ_offA[Gi] + i;
              for(int j=0; j < occpiA[Gi]; j++) {
                  int J = occ_offA[Gi] + j;
                  for(int m=0; m < occpiA[Gm]; m++) {
                      int M = occ_offA[Gm] + m;

                      int IM = K.params->rowidx[I][M];
                      int JM = K.params->colidx[J][M];

                      F.matrix[Gi][i][j] += 2.0 * K.matrix[h][IM][JM];

                    }
                }
            }
    }
    global_dpd_->buf4_mat_irrep_close(&K, h);
  }


      // part-2
      global_dpd_->buf4_mat_irrep_init(&K, 0);
      global_dpd_->buf4_mat_irrep_rd(&K, 0);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirrep_; Gi++) {
          int Gj=Gi;
	  for(int Gm=0; Gm < nirrep_; Gm++) {

          /* Loop over the orbitals of the target */
          for(int i=0; i < occpiA[Gi]; i++) {
              int I = occ_offA[Gi] + i;
              for(int j=0; j < occpiA[Gj]; j++) {
                  int J = occ_offA[Gj] + j;
                  for(int m=0; m < occpiA[Gm]; m++) {
                      int M = occ_offA[Gm] + m;

		      int IJ = K.params->rowidx[I][J];
                      int MM = K.params->colidx[M][M];

                      F.matrix[Gi][i][j] -= K.matrix[0][IJ][MM];

                    }
                }
            }
	  }
    }
    global_dpd_->buf4_mat_irrep_close(&K, 0);

  /* Close the Integral buffers */
  global_dpd_->buf4_close(&K);
  global_dpd_->file2_mat_wrt(&F);
  global_dpd_->file2_mat_close(&F);
  global_dpd_->file2_close(&F);

/************************************************************************************************/
/*********************************** Build Fab **************************************************/
/************************************************************************************************/
  /* Prepare the alpha and beta vir-vir Fock matrix files */
  global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "Fock <V|V>");
  global_dpd_->file2_mat_init(&F);

  /* Two-electron contributions */

  /* Prepare the <OV|OV> integral buffers */
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,V]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,V]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OV|OV>");
  for(int h=0; h < nirrep_; h++) {

      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Ga=0; Ga < nirrep_; Ga++) {
	  int Gb = Ga;
	  int Gm = h^Ga;

	  /* Loop over orbitals of the target */
	  for(int a=0; a < virtpiA[Ga]; a++) {
	      int A = vir_offA[Ga] + a;
	      for(int b=0; b < virtpiA[Gb]; b++) {
		  int B = vir_offA[Gb] + b;

		  for(int m=0; m < occpiA[Gm]; m++) {
		      int M = occ_offA[Gm] + m;

		      int MA = K.params->rowidx[M][A];
		      int MB = K.params->colidx[M][B];

		      F.matrix[Ga][a][b] += 2.0 * K.matrix[h][MA][MB];

		    }
		}
	    }
    }
    global_dpd_->buf4_mat_irrep_close(&K, h);
  }
  global_dpd_->buf4_close(&K);


  /* Prepare the <OO|VV> integral buffers */
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <OO|VV>");
      global_dpd_->buf4_mat_irrep_init(&K, 0);
      global_dpd_->buf4_mat_irrep_rd(&K, 0);

      /* Loop over irreps of the target */
      for(int Ga=0; Ga < nirrep_; Ga++) {
	  int Gb = Ga;
	  for(int Gm=0; Gm < nirrep_; Gm++) {

	  /* Loop over orbitals of the target */
	  for(int a=0; a < virtpiA[Ga]; a++) {
	      int A = vir_offA[Ga] + a;
	      for(int b=0; b < virtpiA[Gb]; b++) {
		  int B = vir_offA[Gb] + b;

		  for(int m=0; m < occpiA[Gm]; m++) {
		      int M = occ_offA[Gm] + m;

		      int MM = K.params->rowidx[M][M];
		      int AB = K.params->colidx[A][B];

		      F.matrix[Ga][a][b] -= K.matrix[0][MM][AB];

		    }
		}
	    }
	  }
    }
    global_dpd_->buf4_mat_irrep_close(&K, 0);

  /* Close the <OV|OV> integral buffers */
  global_dpd_->buf4_close(&K);
  global_dpd_->file2_mat_wrt(&F);
  global_dpd_->file2_mat_close(&F);
  global_dpd_->file2_close(&F);

/************************************************************************************************/
/*********************************** Build Fia **************************************************/
/************************************************************************************************/
  /* Prepare the alpha and beta occ-vir Fock matrix files */
  global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "Fock <O|V>");
  global_dpd_->file2_mat_init(&F);


  /* Two-electron contributions */

  /* Prepare the <OO|OV> integral buffers */
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OO|OV>");
  // part-1
  for(int h=0; h < nirrep_; h++) {

      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirrep_; Gi++) {
	  int Ga = Gi; int Gm = h^Gi;

	  /* Loop over orbitals of the target */
	  for(int i=0; i < occpiA[Gi]; i++) {
	      int I = occ_offA[Gi] + i;
	      for(int a=0; a < virtpiA[Ga]; a++) {
		  int A = vir_offA[Ga] + a;

		  for(int m=0; m < occpiA[Gm]; m++) {
		      int M = occ_offA[Gm] + m;

		      int MI = K.params->rowidx[M][I];
		      int MA = K.params->colidx[M][A];

		      F.matrix[Gi][i][a] += 2.0 * K.matrix[h][MI][MA];

		    }
		}
	    }
    }
    global_dpd_->buf4_mat_irrep_close(&K, h);
  }

      // part-2
      global_dpd_->buf4_mat_irrep_init(&K, 0);
      global_dpd_->buf4_mat_irrep_rd(&K, 0);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirrep_; Gi++) {
	  int Ga = Gi;
	  for(int Gm=0; Gm < nirrep_; Gm++) {

	  /* Loop over orbitals of the target */
	  for(int i=0; i < occpiA[Gi]; i++) {
	      int I = occ_offA[Gi] + i;
	      for(int a=0; a < virtpiA[Ga]; a++) {
		  int A = vir_offA[Ga] + a;

		  for(int m=0; m < occpiA[Gm]; m++) {
		      int M = occ_offA[Gm] + m;

		      int MM = K.params->rowidx[M][M];
		      int IA = K.params->colidx[I][A];

		      F.matrix[Gi][i][a] -= K.matrix[0][MM][IA];

		    }
		}
	    }
	  }
    }
    global_dpd_->buf4_mat_irrep_close(&K, 0);

  /* Close the <OO|OV> integral buffers */
  global_dpd_->buf4_close(&K);
  global_dpd_->file2_mat_wrt(&F);
  global_dpd_->file2_mat_close(&F);
  global_dpd_->file2_close(&F);

/************************************************************************************************/
/*********************************** Set Fock ***************************************************/
/************************************************************************************************/

    // <O|O> block
    global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "Fock <O|O>");
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
		FockA->set(h, i, j, F.matrix[h][i][j]);
            }
        }
    }
    global_dpd_->file2_close(&F);

    // <V|V> block
    global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "Fock <V|V>");
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < virtpiA[h]; ++i){
            for(int j = 0 ; j < virtpiA[h]; ++j){
               FockA->set(h, i + occpiA[h], j + occpiA[h], F.matrix[h][i][j]);
            }
        }
    }
    global_dpd_->file2_close(&F);

    // <O|V> block
    global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "Fock <O|V>");
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < virtpiA[h]; ++j){
               FockA->set(h, i, j + occpiA[h], F.matrix[h][i][j]);
	       FockA->set(h, j + occpiA[h], i , F.matrix[h][i][j]);
            }
        }
    }
    global_dpd_->file2_close(&F);

     //1e-contr.
    FockA->add(HmoA);


/************************************************************************************************/
/*********************************** Print Fock *************************************************/
/************************************************************************************************/
	if (print_ > 2) FockA->print();

}// end if (reference_ == "RESTRICTED")


//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference_ == "UNRESTRICTED") {

  // F(pq) = h(pq) + \sum_{m} <pm||qm> in spin-orbital form.

  // F(PQ) = h(PQ) + \sum_{M} <PM||QM>  + \sum_{m} <Pm|Qm> in spin-adapted form.

/************************************************************************************************/
/*********************************** Build Fij **************************************************/
/************************************************************************************************/
  // F(IJ) = h(IJ) + \sum_{M} <IM||JM> + \sum_{m} <Im|Jm>

  dpdfile2 F;
  dpdbuf4 K;

  /* Prepare the alpha occ-occ Fock matrix files */
  global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "Fock <O|O>");
  global_dpd_->file2_mat_init(&F);

  /* Two-electron contributions */

  /* Prepare the <OO||OO> integral buffers */
   global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO||OO>");

  for(int h=0; h < nirrep_; h++) {

      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirrep_; Gi++) {
	  int Gm=h^Gi;

          /* Loop over the orbitals of the target */
          for(int i=0; i < occpiA[Gi]; i++) {
              int I = occ_offA[Gi] + i;
              for(int j=0; j < occpiA[Gi]; j++) {
                  int J = occ_offA[Gi] + j;
                  for(int m=0; m < occpiA[Gm]; m++) {
                      int M = occ_offA[Gm] + m;

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
          for(int i=0; i < occpiA[Gi]; i++) {
              int I = occ_offA[Gi] + i;
              for(int j=0; j < occpiA[Gi]; j++) {
                  int J = occ_offA[Gi] + j;
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

  /* Close the Integral buffers */
  global_dpd_->file2_mat_wrt(&F);
  global_dpd_->file2_mat_close(&F);
  global_dpd_->file2_close(&F);

/************************************************************************************************/
/*********************************** Build Fab **************************************************/
/************************************************************************************************/
  // F(AB) = h(AB) + \sum_{M} <AM||BM>  + \sum_{m} <Am|Bm>

  /* Prepare the alpha vir-vir Fock matrix files */
  global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "Fock <V|V>");
  global_dpd_->file2_mat_init(&F);

  /* Two-electron contributions */

  /* Prepare the <OV||OV> integral buffers */
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV||OV>");
  for(int h=0; h < nirrep_; h++) {

      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Ga=0; Ga < nirrep_; Ga++) {
	  int Gb = Ga;
	  int Gm = h^Ga;

	  /* Loop over orbitals of the target */
	  for(int a=0; a < virtpiA[Ga]; a++) {
	      int A = vir_offA[Ga] + a;

	      for(int b=0; b < virtpiA[Gb]; b++) {
		  int B = vir_offA[Gb] + b;

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


   /* Prepare the <Vo|Vo> integral buffers */
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "MO Ints <Vo|Vo>");
  for(int h=0; h < nirrep_; h++) {

      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Ga=0; Ga < nirrep_; Ga++) {
	  int Gb = Ga;
	  int Gm = h^Ga;

	  /* Loop over orbitals of the target */
	  for(int a=0; a < virtpiA[Ga]; a++) {
	      int A = vir_offA[Ga] + a;

	      for(int b=0; b < virtpiA[Gb]; b++) {
		  int B = vir_offA[Gb] + b;

		  for(int m=0; m < occpiB[Gm]; m++) {
		      int M = occ_offB[Gm] + m;

		      int AM = K.params->rowidx[A][M];
		      int BM = K.params->colidx[B][M];

		      F.matrix[Ga][a][b] += K.matrix[h][AM][BM];

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
  // F(IA) = h(IA) + \sum_{M} <IM||AM> + \sum_{m} <Im|Am>

  /* Prepare the alpha and beta occ-vir FockA matrix files */
  global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "Fock <O|V>");
  global_dpd_->file2_mat_init(&F);


  /* Two-electron contributions */

  /* Prepare the <OO||OV> integral buffers */
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO||OV>");

  for(int h=0; h < nirrep_; h++) {

      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirrep_; Gi++) {
	  int Ga = Gi; int Gm = h^Gi;

	  /* Loop over orbitals of the target */
	  for(int i=0; i < occpiA[Gi]; i++) {
	      int I = occ_offA[Gi] + i;
	      for(int a=0; a < virtpiA[Ga]; a++) {
		  int A = vir_offA[Ga] + a;

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


   /* Prepare the <Oo|Vo> integral buffers */
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "MO Ints <Oo|Vo>");

  for(int h=0; h < nirrep_; h++) {

      global_dpd_->buf4_mat_irrep_init(&K, h);
      global_dpd_->buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirrep_; Gi++) {
	  int Ga = Gi; int Gm = h^Gi;

	  /* Loop over orbitals of the target */
	  for(int i=0; i < occpiA[Gi]; i++) {
	      int I = occ_offA[Gi] + i;
	      for(int a=0; a < virtpiA[Ga]; a++) {
		  int A = vir_offA[Ga] + a;

		  for(int m=0; m < occpiB[Gm]; m++) {
		      int M = occ_offB[Gm] + m;

		      int IM = K.params->rowidx[I][M];
		      int AM = K.params->colidx[A][M];

		      F.matrix[Gi][i][a] += K.matrix[h][IM][AM];

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

    // <O|O> block
    global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "Fock <O|O>");
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
		FockA->set(h, i, j, F.matrix[h][i][j]);
            }
        }
    }
    global_dpd_->file2_close(&F);

    // <V|V> block
    global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "Fock <V|V>");
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < virtpiA[h]; ++i){
            for(int j = 0 ; j < virtpiA[h]; ++j){
               FockA->set(h, i + occpiA[h], j + occpiA[h], F.matrix[h][i][j]);
            }
        }
    }
    global_dpd_->file2_close(&F);

    // <O|V> block
    global_dpd_->file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "Fock <O|V>");
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < virtpiA[h]; ++j){
               FockA->set(h, i, j + occpiA[h], F.matrix[h][i][j]);
	       FockA->set(h, j + occpiA[h], i , F.matrix[h][i][j]);
            }
        }
    }
    global_dpd_->file2_close(&F);

     //1e-contr.
     FockA->add(HmoA);


/************************************************************************************************/
/*********************************** Print FockA ************************************************/
/************************************************************************************************/
	if (print_ > 1) FockA->print();

}// end if (reference_ == "UNRESTRICTED")


}// end of code
}} // End Namespaces
