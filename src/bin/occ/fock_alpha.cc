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
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <libtrans/mospace.h>
#include <libtrans/integraltransform.h>


/** Required libmints includes */
#include <libmints/factory.h>
#include <libmints/wavefunction.h>
#include <libmints/mints.h>

#include "defines.h"
#include "occwave.h"

using namespace boost;
using namespace psi;
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
  dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "Fock <O|O>");
  dpd_file2_mat_init(&F);

  /* Two-electron contributions */

  /* Prepare the <OO|OO> integral buffers */
   dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,O]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,O]"), 0, "MO Ints <OO|OO>");
  // part-1
  for(int h=0; h < nirrep_; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

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
    dpd_buf4_mat_irrep_close(&K, h);
  }
  
  
      // part-2
      dpd_buf4_mat_irrep_init(&K, 0);
      dpd_buf4_mat_irrep_rd(&K, 0);

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
    dpd_buf4_mat_irrep_close(&K, 0);

  /* Close the Integral buffers */
  dpd_buf4_close(&K);
  dpd_file2_mat_wrt(&F);
  dpd_file2_mat_close(&F);
  dpd_file2_close(&F);
  
/************************************************************************************************/
/*********************************** Build Fab **************************************************/
/************************************************************************************************/  
  /* Prepare the alpha and beta vir-vir Fock matrix files */
  dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "Fock <V|V>");
  dpd_file2_mat_init(&F);

  /* Two-electron contributions */

  /* Prepare the <OV|OV> integral buffers */
  dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,V]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,V]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OV|OV>");
  for(int h=0; h < nirrep_; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

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
    dpd_buf4_mat_irrep_close(&K, h);
  }
  dpd_buf4_close(&K);
  
  
  /* Prepare the <OO|VV> integral buffers */
  dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <OO|VV>");
      dpd_buf4_mat_irrep_init(&K, 0);
      dpd_buf4_mat_irrep_rd(&K, 0);

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
    dpd_buf4_mat_irrep_close(&K, 0);

  /* Close the <OV|OV> integral buffers */
  dpd_buf4_close(&K);
  dpd_file2_mat_wrt(&F);
  dpd_file2_mat_close(&F);
  dpd_file2_close(&F);
  
/************************************************************************************************/
/*********************************** Build Fia **************************************************/
/************************************************************************************************/
  /* Prepare the alpha and beta occ-vir Fock matrix files */
  dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "Fock <O|V>");
  dpd_file2_mat_init(&F);


  /* Two-electron contributions */

  /* Prepare the <OO|OV> integral buffers */
  dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OO|OV>");
  // part-1
  for(int h=0; h < nirrep_; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

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
    dpd_buf4_mat_irrep_close(&K, h);
  }
  
      // part-2
      dpd_buf4_mat_irrep_init(&K, 0);
      dpd_buf4_mat_irrep_rd(&K, 0);

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
    dpd_buf4_mat_irrep_close(&K, 0);

  /* Close the <OO|OV> integral buffers */
  dpd_buf4_close(&K);
  dpd_file2_mat_wrt(&F);
  dpd_file2_mat_close(&F);
  dpd_file2_close(&F);
  
/************************************************************************************************/
/*********************************** Set Fock ***************************************************/
/************************************************************************************************/ 
   
    // <O|O> block
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "Fock <O|O>");
    dpd_file2_mat_init(&F);
    dpd_file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
		FockA->set(h, i, j, F.matrix[h][i][j]);
            }
        }
    }
    dpd_file2_close(&F);
    
    // <V|V> block
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "Fock <V|V>");
    dpd_file2_mat_init(&F);
    dpd_file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < virtpiA[h]; ++i){
            for(int j = 0 ; j < virtpiA[h]; ++j){
               FockA->set(h, i + occpiA[h], j + occpiA[h], F.matrix[h][i][j]);
            }
        }
    }
    dpd_file2_close(&F);
    
    // <O|V> block
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "Fock <O|V>");
    dpd_file2_mat_init(&F);
    dpd_file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < virtpiA[h]; ++j){
               FockA->set(h, i, j + occpiA[h], F.matrix[h][i][j]);
	       FockA->set(h, j + occpiA[h], i , F.matrix[h][i][j]);
            }
        }
    }
    dpd_file2_close(&F);
    
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
  dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "Fock <O|O>");
  dpd_file2_mat_init(&F);

  /* Two-electron contributions */

  /* Prepare the <OO||OO> integral buffers */
   dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO||OO>");
   
  for(int h=0; h < nirrep_; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

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
    dpd_buf4_mat_irrep_close(&K, h);
  }
  dpd_buf4_close(&K);
  
    
  /* Prepare the <Oo|Oo> integral buffers */
   dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");  

  for(int h=0; h < nirrep_; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

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
    dpd_buf4_mat_irrep_close(&K, h);
  }
  dpd_buf4_close(&K);

  /* Close the Integral buffers */
  dpd_file2_mat_wrt(&F);
  dpd_file2_mat_close(&F);
  dpd_file2_close(&F);
  
/************************************************************************************************/
/*********************************** Build Fab **************************************************/
/************************************************************************************************/ 
  // F(AB) = h(AB) + \sum_{M} <AM||BM>  + \sum_{m} <Am|Bm> 
  
  /* Prepare the alpha vir-vir Fock matrix files */
  dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "Fock <V|V>");
  dpd_file2_mat_init(&F);

  /* Two-electron contributions */

  /* Prepare the <OV||OV> integral buffers */
  dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV||OV>");
  for(int h=0; h < nirrep_; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

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
    dpd_buf4_mat_irrep_close(&K, h);
  }
  dpd_buf4_close(&K);
  
    
   /* Prepare the <Vo|Vo> integral buffers */
  dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "MO Ints <Vo|Vo>");
  for(int h=0; h < nirrep_; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

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
    dpd_buf4_mat_irrep_close(&K, h);
  }
  dpd_buf4_close(&K);
    

  /* Close the buffers */  
  dpd_file2_mat_wrt(&F);
  dpd_file2_mat_close(&F);
  dpd_file2_close(&F);
  
/************************************************************************************************/
/*********************************** Build Fia **************************************************/
/************************************************************************************************/
  // F(IA) = h(IA) + \sum_{M} <IM||AM> + \sum_{m} <Im|Am>
  
  /* Prepare the alpha and beta occ-vir FockA matrix files */
  dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "Fock <O|V>");
  dpd_file2_mat_init(&F);


  /* Two-electron contributions */

  /* Prepare the <OO||OV> integral buffers */
  dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO||OV>");
  
  for(int h=0; h < nirrep_; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

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
    dpd_buf4_mat_irrep_close(&K, h);
  }
  dpd_buf4_close(&K);
    
    
   /* Prepare the <Oo|Vo> integral buffers */
  dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "MO Ints <Oo|Vo>");
  
  for(int h=0; h < nirrep_; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

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
    dpd_buf4_mat_irrep_close(&K, h);
  }
  dpd_buf4_close(&K);

  /* Close the buffers */  
  dpd_file2_mat_wrt(&F);
  dpd_file2_mat_close(&F);
  dpd_file2_close(&F);
  
/************************************************************************************************/
/*********************************** Set Fock ***************************************************/
/************************************************************************************************/ 
   
    // <O|O> block
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "Fock <O|O>");
    dpd_file2_mat_init(&F);
    dpd_file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
		FockA->set(h, i, j, F.matrix[h][i][j]);
            }
        }
    }
    dpd_file2_close(&F);
    
    // <V|V> block
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "Fock <V|V>");
    dpd_file2_mat_init(&F);
    dpd_file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < virtpiA[h]; ++i){
            for(int j = 0 ; j < virtpiA[h]; ++j){
               FockA->set(h, i + occpiA[h], j + occpiA[h], F.matrix[h][i][j]);
            }
        }
    }
    dpd_file2_close(&F);
    
    // <O|V> block
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "Fock <O|V>");
    dpd_file2_mat_init(&F);
    dpd_file2_mat_rd(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < virtpiA[h]; ++j){
               FockA->set(h, i, j + occpiA[h], F.matrix[h][i][j]);
	       FockA->set(h, j + occpiA[h], i , F.matrix[h][i][j]);
            }
        }
    }
    dpd_file2_close(&F);
    
     //1e-contr.
     FockA->add(HmoA);    
    
  
/************************************************************************************************/
/*********************************** Print FockA ************************************************/
/************************************************************************************************/
	if (print_ > 1) FockA->print();

}// end if (reference_ == "UNRESTRICTED") 
	

}// end of code
}} // End Namespaces

