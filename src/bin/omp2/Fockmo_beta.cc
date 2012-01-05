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
#include "omp2wave.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace omp2wave{

void OMP2Wave::Fockmo_beta()
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
  dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "Fock <o|o>");
  dpd_file2_mat_init(&F);

  /* Two-electron contributions */

  /* Prepare the <oo||oo> integral buffers */
   dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo||oo>");
   
  for(int h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirreps; Gi++) {
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
    dpd_buf4_mat_irrep_close(&K, h);
  }
  dpd_buf4_close(&K);
    
    
  /* Prepare the <Oo|Oo> integral buffers */
   dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");  

  for(int h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirreps; Gi++) {
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
  // F(ab) = h(ab) + \sum_{m} <am|bm> + \sum_{M} <Ma|Mb> 
  
  /* Prepare the alpha vir-vir Fock matrix files */
  dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "Fock <v|v>");
  dpd_file2_mat_init(&F);

  /* Two-electron contributions */

  /* Prepare the <ov||ov> integral buffers */
  dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov||ov>");
  
  for(int h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Ga=0; Ga < nirreps; Ga++) {
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
    dpd_buf4_mat_irrep_close(&K, h);
  }
  dpd_buf4_close(&K);
  
    
   /* Prepare the <Ov|Ov> integral buffers */
  dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>");
  
  for(int h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Ga=0; Ga < nirreps; Ga++) {
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
  // F(ia) = h(ia) + \sum_{m} <im||am> + \sum_{M} <Mi|Ma>
  
  /* Prepare the alpha and beta occ-vir FockB matrix files */
  dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "Fock <o|v>");
  dpd_file2_mat_init(&F);


  /* Two-electron contributions */

  /* Prepare the <oo||ov> integral buffers */
  dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo||ov>");

  for(int h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirreps; Gi++) {
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
    dpd_buf4_mat_irrep_close(&K, h);
  }
  dpd_buf4_close(&K);
    
    
   /* Prepare the <Oo|Ov> integral buffers */
  dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "MO Ints <Oo|Ov>");

  for(int h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&K, h);
      dpd_buf4_mat_irrep_rd(&K, h);

      /* Loop over irreps of the target */
      for(int Gi=0; Gi < nirreps; Gi++) {
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
   
    // <o|o> block
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "Fock <o|o>");
    dpd_file2_mat_init(&F);
    dpd_file2_mat_rd(&F);
    for(int h = 0; h < nirreps; ++h){
        for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < occpiB[h]; ++j){
		FockB->set(h, i, j, F.matrix[h][i][j]);
            }
        }
    }
    dpd_file2_close(&F);
    
    // <v|v> block
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "Fock <v|v>");
    dpd_file2_mat_init(&F);
    dpd_file2_mat_rd(&F);
    for(int h = 0; h < nirreps; ++h){
        for(int i = 0 ; i < virtpiB[h]; ++i){
            for(int j = 0 ; j < virtpiB[h]; ++j){
               FockB->set(h, i + occpiB[h], j + occpiB[h], F.matrix[h][i][j]);
            }
        }
    }
    dpd_file2_close(&F);
    
    // <o|v> block
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "Fock <o|v>");
    dpd_file2_mat_init(&F);
    dpd_file2_mat_rd(&F);
    for(int h = 0; h < nirreps; ++h){
        for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < virtpiB[h]; ++j){
               FockB->set(h, i, j + occpiB[h], F.matrix[h][i][j]);
	       FockB->set(h, j + occpiB[h], i , F.matrix[h][i][j]);
            }
        }
    }
    dpd_file2_close(&F);
    
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

