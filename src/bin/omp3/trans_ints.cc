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
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <libtrans/mospace.h>
#include <libtrans/integraltransform.h>

/** Required libmints includes */
#include <libmints/factory.h>
#include <libmints/wavefunction.h>
#include <libmints/mints.h>


#include "omp3wave.h"
#include "defines.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace omp3wave{

void OMP3Wave::trans_ints()
{    
    //fprintf(outfile,"\n trans_ints is starting... \n"); fflush(outfile);
/********************************************************************************************/
/************************** Transform 2-electron int. to MO space ***************************/
/********************************************************************************************/  
    ints->update_orbitals();  
    ints->set_print(print_ - 2 >= 0 ? print_ - 2 : 0);
    
    // Trans (OO|OO)
    ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ);
    
    // Trans (OO|OV)
    ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::vir);
    
    // Trans (OV|OO)
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::occ);

    // Trans (OV|OV)
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
    
    // Trans (OO|VV)
    ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir);
    
    // Trans (VV|OO)
    ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::occ);
    
    // Trans (OV|VV)
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::vir, MOSpace::vir);
    
    // Trans (VV|OV)
    ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::vir);
    
    
    // Trans (VV|VV)
    ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::vir, MOSpace::vir);
    
    
/********************************************************************************************/
/************************** sort chem -> phys ***********************************************/
/********************************************************************************************/  
     dpdbuf4 K, G;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     
     // Build MO ints    
     
     // (OO|OO) -> <OO|OO>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[O,O]"), "MO Ints <OO|OO>");
     dpd_buf4_close(&K);
     
     // (oo|oo) -> <oo|oo>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o>=o]+"), ID("[o>=o]+"), 0, "MO Ints (oo|oo)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[o,o]"), ID("[o,o]"), "MO Ints <oo|oo>");
     dpd_buf4_close(&K);
     
     // (OO|oo) -> <Oo|Oo>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,o]"),
                  ID("[O>=O]+"), ID("[o>=o]+"), 0, "MO Ints (OO|oo)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,o]"), ID("[O,o]"), "MO Ints <Oo|Oo>");
     dpd_buf4_close(&K);
     
 
  
     // (OO|OV) -> <OO|OV>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O>=O]+"), ID("[O,V]"), 0, "MO Ints (OO|OV)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[O,V]"), "MO Ints <OO|OV>");
     dpd_buf4_close(&K);
     
     // (oo|ov) -> <oo|ov>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o>=o]+"), ID("[o,v]"), 0, "MO Ints (oo|ov)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[o,o]"), ID("[o,v]"), "MO Ints <oo|ov>");
     dpd_buf4_close(&K);
     
     // (OO|ov) -> <Oo|Ov>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,v]"),
                  ID("[O>=O]+"), ID("[o,v]"), 0, "MO Ints (OO|ov)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,o]"), ID("[O,v]"), "MO Ints <Oo|Ov>");
     dpd_buf4_close(&K);
     
     // (OV|oo) -> <Oo|Vo>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,o]"),
                  ID("[O,V]"), ID("[o>=o]+"), 0, "MO Ints (OV|oo)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,o]"), ID("[V,o]"), "MO Ints <Oo|Vo>");
     dpd_buf4_close(&K);
          
     
     
     
     // (OV|OV) -> <OO|VV>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");
     dpd_buf4_close(&K);
     
     // (ov|ov) -> <oo|vv>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[o,o]"), ID("[v,v]"), "MO Ints <oo|vv>");
     dpd_buf4_close(&K);
     
     // (OV|ov) -> <Oo|Vv>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,o]"), ID("[V,v]"), "MO Ints <Oo|Vv>");
     dpd_buf4_close(&K);
     
     // (OV|ov) -> <Ov|Vo>: <Ia||Bj> = <Ia|Bj> = (IB|ja)
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , psqr, ID("[O,v]"), ID("[V,o]"), "MO Ints <Ov|Vo>");
     dpd_buf4_close(&K);
     
    
     // (OO|VV) -> <OV|OV>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "MO Ints <OV|OV>");
     dpd_buf4_close(&K);
     
     // (oo|vv) -> <ov|ov>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>=o]+"), ID("[v>=v]+"), 0, "MO Ints (oo|vv)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[o,v]"), ID("[o,v]"), "MO Ints <ov|ov>");
     dpd_buf4_close(&K);
     
     // (OO|vv) -> <Ov|Ov>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                  ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,v]"), ID("[O,v]"), "MO Ints <Ov|Ov>");
     dpd_buf4_close(&K);
     
     // (VV|oo) -> <Vo|Vo>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                  ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[V,o]"), ID("[V,o]"), "MO Ints <Vo|Vo>");
     dpd_buf4_close(&K);
     
     
     
     
     // (OV|VV) -> <OV|VV>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,V]"), ID("[V,V]"), "MO Ints <OV|VV>");
     dpd_buf4_close(&K);
     
     // (ov|vv) -> <ov|vv>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v>=v]+"), 0, "MO Ints (ov|vv)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[o,v]"), ID("[v,v]"), "MO Ints <ov|vv>");
     dpd_buf4_close(&K);
     
     // (OV|vv) -> <Ov|Vv>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"),
                  ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,v]"), ID("[V,v]"), "MO Ints <Ov|Vv>");
     dpd_buf4_close(&K);
     
     // (VV|ov) -> <Vo|Vv>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,v]"),
                  ID("[V>=V]+"), ID("[o,v]"), 0, "MO Ints (VV|ov)");
     dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[V,o]"), ID("[V,v]"), "MO Ints <Vo|Vv>");
     dpd_buf4_close(&K);
     
     
       
       // (VV|VV) -> <VV|VV>
      dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V>=V]+"), ID("[V>=V]+"), 0, "MO Ints (VV|VV)");
      dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[V,V]"), ID("[V,V]"), "MO Ints <VV|VV>");
      dpd_buf4_close(&K);
      
      // (vv|vv) -> <vv|vv>
      dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v>=v]+"), ID("[v>=v]+"), 0, "MO Ints (vv|vv)");
      dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[v,v]"), ID("[v,v]"), "MO Ints <vv|vv>");
      dpd_buf4_close(&K);
      
      // (VV|vv) -> <Vv|Vv>
      dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[v,v]"),
                  ID("[V>=V]+"), ID("[v>=v]+"), 0, "MO Ints (VV|vv)");
      dpd_buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[V,v]"), ID("[V,v]"), "MO Ints <Vv|Vv>");
      dpd_buf4_close(&K);
     
/********************************************************************************************/
/************************** Antisymmetrized Ints ********************************************/
/********************************************************************************************/
     // <OO||OO>:  <IJ||KL> =  <IJ|KL> - <IJ|LK> 
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO|OO>");
     dpd_buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <OO||OO>");
     dpd_buf4_close(&K);
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO||OO>");
     dpd_buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO|OO>");
     for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            for(int kl = 0; kl < K.params->coltot[h]; ++kl){
                int k = K.params->colorb[h][kl][0];
                int l = K.params->colorb[h][kl][1];
		int lk = G.params->colidx[l][k];
                K.matrix[h][ij][kl] -= G.matrix[h][ij][lk];
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
        dpd_buf4_mat_irrep_close(&G, h);
     }
     dpd_buf4_close(&K);
     dpd_buf4_close(&G);

 
     // <oo||oo>:  <ij||kl> =  <ij|kl> - <ij|lk> 
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo|oo>");
     dpd_buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <oo||oo>");
     dpd_buf4_close(&K);
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo||oo>");
     dpd_buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo|oo>");
     for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            for(int kl = 0; kl < K.params->coltot[h]; ++kl){
                int k = K.params->colorb[h][kl][0];
                int l = K.params->colorb[h][kl][1];
		int lk = G.params->colidx[l][k];
                K.matrix[h][ij][kl] -= G.matrix[h][ij][lk];
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
        dpd_buf4_mat_irrep_close(&G, h);
     }
     dpd_buf4_close(&K);
     dpd_buf4_close(&G); 
   

     // <OO||OV>:  <IJ||KA> = <IJ|KA> - <JI|KA>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO|OV>");
     dpd_buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <OO||OV>");
     dpd_buf4_close(&K);
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO||OV>");
     dpd_buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO|OV>");
     for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            int i = K.params->roworb[h][ij][0];
            int j = K.params->roworb[h][ij][1];
            int ji = G.params->rowidx[j][i];
            for(int ka = 0; ka < K.params->coltot[h]; ++ka){ 
                K.matrix[h][ij][ka] -= G.matrix[h][ji][ka];
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
        dpd_buf4_mat_irrep_close(&G, h);
     }
     dpd_buf4_close(&K);
     dpd_buf4_close(&G);
    

    
     // <oo||ov>:   <ij||ka> = <ij|ka> - <ji|ka>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo|ov>");
     dpd_buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <oo||ov>");
     dpd_buf4_close(&K);
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo||ov>");
     dpd_buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo|ov>");
     for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            int i = K.params->roworb[h][ij][0];
            int j = K.params->roworb[h][ij][1];
            int ji = G.params->rowidx[j][i];
            for(int ka = 0; ka < K.params->coltot[h]; ++ka){ 
                K.matrix[h][ij][ka] -= G.matrix[h][ji][ka];
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
        dpd_buf4_mat_irrep_close(&G, h);
     }
     dpd_buf4_close(&K);
     dpd_buf4_close(&G);


    
     // <OO||VV>:  <IJ||AB> = <IJ|AB> - <IJ|BA> 
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     dpd_buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <OO||VV>");
     dpd_buf4_close(&K);
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     dpd_buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
	dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            for(int ab = 0; ab < K.params->coltot[h]; ++ab){
                int a = K.params->colorb[h][ab][0];
                int b = K.params->colorb[h][ab][1];
		int ba = G.params->colidx[b][a];
		K.matrix[h][ij][ab] -= G.matrix[h][ij][ba];                
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);

     
     // <oo||vv>:  <ij||ab> = <ij|ab> - <ij|ba> 
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
     dpd_buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <oo||vv>");
     dpd_buf4_close(&K);
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     dpd_buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
     for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
	dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            for(int ab = 0; ab < K.params->coltot[h]; ++ab){
                int a = K.params->colorb[h][ab][0];
                int b = K.params->colorb[h][ab][1];
		int ba = G.params->colidx[b][a];
		K.matrix[h][ij][ab] -= G.matrix[h][ij][ba];                
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);

    
    
     // <OV||OV>:  <IA||JB> = <IA|JB> - (IB|JA)
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
     dpd_buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <OV||OV>");
     dpd_buf4_close(&K);
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV||OV>");
     dpd_buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
     for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        for(int ia = 0; ia < K.params->rowtot[h]; ++ia){
            int i = K.params->roworb[h][ia][0];
            int a = K.params->roworb[h][ia][1];
            int Gi = K.params->psym[i];
            for(int jb = 0; jb < K.params->coltot[h]; ++jb){
                int j = K.params->colorb[h][jb][0];
                int b = K.params->colorb[h][jb][1];
                int Gb = K.params->ssym[b];
                int Gib = Gi^Gb;
	        dpd_buf4_mat_irrep_init(&G, Gib);
                dpd_buf4_mat_irrep_rd(&G, Gib);
		int ib = G.params->rowidx[i][b];
		int ja = G.params->colidx[j][a];
		K.matrix[h][ia][jb] -= G.matrix[Gib][ib][ja];                
                dpd_buf4_mat_irrep_close(&G, Gib);
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
    }
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);


     // <ov||ov>:  <ia||jb> = <ia|jb> - (ib|ja)
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov>");
     dpd_buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <ov||ov>");
     dpd_buf4_close(&K);
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov||ov>");
     dpd_buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
     for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        for(int ia = 0; ia < K.params->rowtot[h]; ++ia){
            int i = K.params->roworb[h][ia][0];
            int a = K.params->roworb[h][ia][1];
            int Gi = K.params->psym[i];
            for(int jb = 0; jb < K.params->coltot[h]; ++jb){
                int j = K.params->colorb[h][jb][0];
                int b = K.params->colorb[h][jb][1];
                int Gb = K.params->ssym[b];
                int Gib = Gi^Gb;
	        dpd_buf4_mat_irrep_init(&G, Gib);
                dpd_buf4_mat_irrep_rd(&G, Gib);
		int ib = G.params->rowidx[i][b];
		int ja = G.params->colidx[j][a];
		K.matrix[h][ia][jb] -= G.matrix[Gib][ib][ja];                
                dpd_buf4_mat_irrep_close(&G, Gib);
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
    }
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
   
 
     // <OV||VV>:  <IA||BC> = <IA|BC> - <IA|CB>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
     dpd_buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <OV||VV>");
     dpd_buf4_close(&K);
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV||VV>");
     dpd_buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
     for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
	dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ia = 0; ia < K.params->rowtot[h]; ++ia){
            for(int bc = 0; bc < K.params->coltot[h]; ++bc){
                int b = K.params->colorb[h][bc][0];
                int c = K.params->colorb[h][bc][1];
		int cb = G.params->colidx[c][b];
		K.matrix[h][ia][bc] -= G.matrix[h][ia][cb];                
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
    

     // <ov||vv>:  <ia||bc> = <ia|bc> - <ia|cb>
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
     dpd_buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <ov||vv>");
     dpd_buf4_close(&K);
     dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov||vv>");
     dpd_buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
     for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
	dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ia = 0; ia < K.params->rowtot[h]; ++ia){
            for(int bc = 0; bc < K.params->coltot[h]; ++bc){
                int b = K.params->colorb[h][bc][0];
                int c = K.params->colorb[h][bc][1];
		int cb = G.params->colidx[c][b];
		K.matrix[h][ia][bc] -= G.matrix[h][ia][cb];                
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&K);
    dpd_buf4_close(&G);
           
     
      // <VV||VV>: <AB||CD> = <AB|CD> - <AB|DC>
      dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
      dpd_buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <VV||VV>");
      dpd_buf4_close(&K);
      dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV||VV>");
      dpd_buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ab = 0; ab < K.params->rowtot[h]; ++ab){
            for(int cd = 0; cd < K.params->coltot[h]; ++cd){
                int c = K.params->colorb[h][cd][0];
                int d = K.params->colorb[h][cd][1];
		int dc = G.params->colidx[d][c];
                K.matrix[h][ab][cd] -= G.matrix[h][ab][dc];
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
        dpd_buf4_mat_irrep_close(&G, h);
     }
     dpd_buf4_close(&K);
     dpd_buf4_close(&G);



      // <vv||vv>: <ab||cd> = <ab|cd> - <ab|dc>
      dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv|vv>");
      dpd_buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <vv||vv>");
      dpd_buf4_close(&K);
      dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv||vv>");  
      dpd_buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv|vv>");
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ab = 0; ab < K.params->rowtot[h]; ++ab){
            for(int cd = 0; cd < K.params->coltot[h]; ++cd){
                int c = K.params->colorb[h][cd][0];
                int d = K.params->colorb[h][cd][1];
		int dc = G.params->colidx[d][c];
                K.matrix[h][ab][cd] -= G.matrix[h][ab][dc];
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
        dpd_buf4_mat_irrep_close(&G, h);
     }
     dpd_buf4_close(&K);
     dpd_buf4_close(&G);


/********************************************************************************************/
/************************** Transform 1-electron int. to MO space ***************************/
/********************************************************************************************/        
      // Trans H matrix
      HmoA->copy(Hso);
      HmoB->copy(Hso);
      HmoA->transform(Ca_);
      HmoB->transform(Cb_);
      
      if (print_ > 1) {
	HmoA->print();
	HmoB->print();
      }
      
      // Trans Fock matrix    
      Fockmo_alpha();      
      Fockmo_beta();   

      denominators();
      psio_->close(PSIF_LIBTRANS_DPD, 1);
      //fprintf(outfile,"\n trans_ints done. \n"); fflush(outfile);
 
}//



void OMP3Wave::denominators()
{
    //fprintf(outfile,"\n denominators is starting... \n"); fflush(outfile);
    dpdbuf4 D;
    dpdfile2 Fo,Fv;
    
    double *aOccEvals = new double [nacooA];
    double *bOccEvals = new double [nacooB];
    double *aVirEvals = new double [nacvoA];
    double *bVirEvals = new double [nacvoB];
    
    // Pick out the diagonal elements of the Fock matrix, making sure that they are in the order
    // used by the DPD library, i.e. starting from zero for each space and ordering by irrep
    
    int aOccCount = 0, bOccCount = 0, aVirCount = 0, bVirCount = 0;
    
    //Diagonal elements of the Fock matrix
    for(int h = 0; h < nirreps; ++h){
        for(int i = 0; i < aoccpiA[h]; ++i) aOccEvals[aOccCount++] = FockA->get(h, i + frzcpi[h], i + frzcpi[h]);
	for(int i = 0; i < aoccpiB[h]; ++i) bOccEvals[bOccCount++] = FockB->get(h, i + frzcpi[h], i + frzcpi[h]);
        for(int a = 0; a < avirtpiA[h]; ++a) aVirEvals[aVirCount++] = FockA->get(h, occpiA[h] + a, occpiA[h] + a); 
	for(int a = 0; a < avirtpiB[h]; ++a) bVirEvals[bVirCount++] = FockB->get(h, occpiB[h] + a, occpiB[h] + a); 
    }
    
    // Build denominators
    // The alpha-alpha spin case
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(aOccEvals[i] + aOccEvals[j] - aVirEvals[a] - aVirEvals[b]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&D, h);
        dpd_buf4_mat_irrep_close(&D, h);
    }
    if (print_ > 2) dpd_buf4_print(&D, outfile, 1);
    dpd_buf4_close(&D);
    

    // The beta-beta spin case 
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(bOccEvals[i] + bOccEvals[j] - bVirEvals[a] - bVirEvals[b]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&D, h);
        dpd_buf4_mat_irrep_close(&D, h);
    }
    if (print_ > 2) dpd_buf4_print(&D, outfile, 1);
    dpd_buf4_close(&D);
    
    
    // The alpha-beta spin case 
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(aOccEvals[i] + bOccEvals[j] - aVirEvals[a] - bVirEvals[b]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&D, h);
        dpd_buf4_mat_irrep_close(&D, h);
    }
    if (print_ > 2) dpd_buf4_print(&D, outfile, 1);
    dpd_buf4_close(&D);
    
    //Print
    if(print_ > 1){
      fprintf(outfile,"\n \n"); fflush(outfile);
      for(int i = 0; i<nacooA; i++) {
	fprintf(outfile,"\taOccEvals[%1d]: %20.14f\n", i, aOccEvals[i]); 
	fflush(outfile);
      }
      
      fprintf(outfile,"\n \n"); fflush(outfile);
      for(int i = 0; i<nacooB; i++) {
	fprintf(outfile,"\tbOccEvals[%1d]: %20.14f\n", i, bOccEvals[i]); 
	fflush(outfile);
      }
      
      fprintf(outfile,"\n \n"); fflush(outfile);
      for(int i = 0; i<nacvoA; i++) {
	fprintf(outfile,"\taVirEvals[%1d]: %20.14f\n", i, aVirEvals[i]); 
	fflush(outfile);
      }
      
      fprintf(outfile,"\n \n"); fflush(outfile);
      for(int i = 0; i<nacvoB; i++) {
	fprintf(outfile,"\tbVirEvals[%1d]: %20.14f\n", i, bVirEvals[i]);
	fflush(outfile);
      }      
    }
    
    delete [] aOccEvals;
    delete [] bOccEvals;
    delete [] aVirEvals;
    delete [] bVirEvals;

 
    // Off-diagonal elements of the Fock matrix    
    // Build Occupied-Occupied block
    // The alpha-alpha spin case 
    dpd_file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    dpd_file2_mat_init(&Fo);
    for(int h = 0; h < nirreps; ++h){
        for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
		if (i != j) Fo.matrix[h][i][j] = FockA->get(h, i + frzcpi[h], j + frzcpi[h]);
		else Fo.matrix[h][i][j] = 0.0;
            }
        }
    }
    dpd_file2_mat_wrt(&Fo);
    dpd_file2_close(&Fo);
    
    if (print_ > 2) {
      dpd_file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
      dpd_file2_mat_init(&Fo);
      dpd_file2_mat_print(&Fo, outfile);
      dpd_file2_close(&Fo);
    }
    
    
    // The beta-beta spin case 
    dpd_file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    dpd_file2_mat_init(&Fo);
    for(int h = 0; h < nirreps; ++h){
        for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int j = 0 ; j < aoccpiB[h]; ++j){
		if (i != j) Fo.matrix[h][i][j] = FockB->get(h, i + frzcpi[h], j + frzcpi[h]);
		else Fo.matrix[h][i][j] = 0.0;
            }
        }
    }
    dpd_file2_mat_wrt(&Fo);
    dpd_file2_close(&Fo);
    
    if (print_ > 2) {
      dpd_file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
      dpd_file2_mat_init(&Fo);
      dpd_file2_mat_print(&Fo, outfile);
      dpd_file2_close(&Fo);
    }
    
    
    // Build Virtual-Virtual block
    // The alpha-alpha spin case 
    dpd_file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    dpd_file2_mat_init(&Fv);
    for(int h = 0; h < nirreps; ++h){
        for(int i = 0 ; i < avirtpiA[h]; ++i){
            for(int j = 0 ; j < avirtpiA[h]; ++j){
                if (i != j) Fv.matrix[h][i][j] = FockA->get(h, i + occpiA[h], j + occpiA[h]);
		else Fv.matrix[h][i][j] = 0.0;
            }
        }
    }
    dpd_file2_mat_wrt(&Fv);
    dpd_file2_close(&Fv);
    
    if (print_ > 2) {
      dpd_file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
      dpd_file2_mat_init(&Fv);
      dpd_file2_mat_print(&Fv, outfile);
      dpd_file2_close(&Fv);
    }
    
    
    // The beta-beta spin case 
    dpd_file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");
    dpd_file2_mat_init(&Fv);
    for(int h = 0; h < nirreps; ++h){
        for(int i = 0 ; i < avirtpiB[h]; ++i){
            for(int j = 0 ; j < avirtpiB[h]; ++j){
                if (i != j) Fv.matrix[h][i][j] = FockB->get(h, i + occpiB[h], j + occpiB[h]);
		else Fv.matrix[h][i][j] = 0.0;
            }
        }
    }
    dpd_file2_mat_wrt(&Fv);
    dpd_file2_close(&Fv);
    
    if (print_ > 2) {
      dpd_file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");
      dpd_file2_mat_init(&Fv);
      dpd_file2_mat_print(&Fv, outfile);
      dpd_file2_close(&Fv);
    }    

//fprintf(outfile,"\n denominators done. \n"); fflush(outfile);
    
}// end denominators


}} // End Namespaces

