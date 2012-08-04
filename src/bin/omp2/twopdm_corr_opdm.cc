/* This code includes correlation opdm contributions. */

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
#include <libqt/qt.h>
#include <libtrans/mospace.h>
#include <libtrans/integraltransform.h>


/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>

#include "omp2wave.h"
#include "defines.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace omp2wave{

void OMP2Wave::twopdm_corr_opdm()
{

    dpdbuf4 G;
    
    psio_->open(PSIF_OMP2_DENSITY, PSIO_OPEN_OLD);
    
    // TPDM OOOO-Block 
    // Alpha-Alpha spin case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");    
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
		
		int hi = G.params->psym[i];
		int hj = G.params->qsym[j];
		int hk = G.params->rsym[k];
		int hl = G.params->ssym[l];
		
		int ii = i - G.params->poff[hi];
		int jj = j - G.params->qoff[hj];
		int kk = k - G.params->roff[hk];
		int ll = l - G.params->soff[hl];
		
		if (i == k && hj == hl) G.matrix[h][ij][kl] += 0.25 * gamma1corrA->get(hj,jj,ll);
		if (j == l && hi == hk) G.matrix[h][ij][kl] += 0.25 * gamma1corrA->get(hi,ii,kk);
		if (i == l && hj == hk) G.matrix[h][ij][kl] -= 0.25 * gamma1corrA->get(hj,jj,kk);
		if (j == k && hi == hl) G.matrix[h][ij][kl] -= 0.25 * gamma1corrA->get(hi,ii,ll);
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    if (print_ > 1) dpd_buf4_print(&G, outfile, 1);
    dpd_buf4_close(&G);
    
  
    // Beta-Beta spin case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");    
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
		
		int hi = G.params->psym[i];
		int hj = G.params->qsym[j];
		int hk = G.params->rsym[k];
		int hl = G.params->ssym[l];
		
		int ii = i - G.params->poff[hi];
		int jj = j - G.params->qoff[hj];
		int kk = k - G.params->roff[hk];
		int ll = l - G.params->soff[hl];
		
		if (i == k && hj == hl) G.matrix[h][ij][kl] += 0.25 * gamma1corrB->get(hj,jj,ll);
		if (j == l && hi == hk) G.matrix[h][ij][kl] += 0.25 * gamma1corrB->get(hi,ii,kk);
		if (i == l && hj == hk) G.matrix[h][ij][kl] -= 0.25 * gamma1corrB->get(hj,jj,kk);
		if (j == k && hi == hl) G.matrix[h][ij][kl] -= 0.25 * gamma1corrB->get(hi,ii,ll);
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    if (print_ > 1) dpd_buf4_print(&G, outfile, 1);
    dpd_buf4_close(&G);
    
    
    // Alpha-Beta spin case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                 ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");  
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
		
		int hi = G.params->psym[i];
		int hj = G.params->qsym[j];
		int hk = G.params->rsym[k];
		int hl = G.params->ssym[l];
		
		int ii = i - G.params->poff[hi];
		int jj = j - G.params->qoff[hj];
		int kk = k - G.params->roff[hk];
		int ll = l - G.params->soff[hl];
		
		if (i == k && hj == hl) G.matrix[h][ij][kl] += 0.25 * gamma1corrB->get(hj,jj,ll);
		if (j == l && hi == hk) G.matrix[h][ij][kl] += 0.25 * gamma1corrA->get(hi,ii,kk);
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    if (print_ > 1) dpd_buf4_print(&G, outfile, 1);
    dpd_buf4_close(&G);
    
 
    // TPDM <OV|OV>
    // Alpha-Alpha spin case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");    
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        //dpd_buf4_mat_irrep_rd(&G, h); // since it is not formed earlier you can not read it.
        for(int ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            for(int jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
		int ha = G.params->qsym[a];
		int hb = G.params->ssym[b];
		int aa = a - G.params->qoff[ha] + occpiA[ha];
		int bb = b - G.params->soff[hb] + occpiA[hb];
		if (i == j && ha == hb) G.matrix[h][ia][jb] = 0.25 * gamma1corrA->get(ha,aa,bb);
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);
    
    
    // Beta-Beta spin case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");    
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        for(int ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            for(int jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
		int ha = G.params->qsym[a];
		int hb = G.params->ssym[b];
		int aa = a - G.params->qoff[ha] + occpiB[ha];
		int bb = b - G.params->soff[hb] + occpiB[hb];
		if (i == j && ha == hb) G.matrix[h][ia][jb] = 0.25 * gamma1corrB->get(ha,aa,bb);
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);
    
    
    // Alpha-Beta spin case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");    
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        for(int ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            for(int jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
		int ha = G.params->qsym[a];
		int hb = G.params->ssym[b];
		int aa = a - G.params->qoff[ha] + occpiB[ha];
		int bb = b - G.params->soff[hb] + occpiB[hb];
		if (i == j && ha == hb) G.matrix[h][ia][jb] = 0.25 * gamma1corrB->get(ha,aa,bb);
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);
    
    
    // TPDM <Vo|Vo>
    // Alpha-Beta spin case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");    
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        for(int ai = 0; ai < G.params->rowtot[h]; ++ai){
            int a = G.params->roworb[h][ai][0];
            int i = G.params->roworb[h][ai][1];
            for(int bj = 0; bj < G.params->coltot[h]; ++bj){
                int b = G.params->colorb[h][bj][0];
                int j = G.params->colorb[h][bj][1];
		int ha = G.params->psym[a];
		int hb = G.params->rsym[b];
		int aa = a - G.params->poff[ha] + occpiA[ha];
		int bb = b - G.params->roff[hb] + occpiA[hb];
		if (i == j && ha == hb) G.matrix[h][ai][bj] = 0.25 * gamma1corrA->get(ha,aa,bb);
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);
    

    // G_IABJ = -GIAJB so I do not need to OVVO block,
    //  however for proper contraction (to avoid construction of <OV||VV>)  in the GFock.cc I need it. 
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    dpd_buf4_sort(&G, PSIF_OMP2_DENSITY , pqsr, ID("[O,V]"), ID("[V,O]"), "TPDM <OV|VO>");
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,V]"), ID("[V,O]"),
                  ID("[O,V]"), ID("[V,O]"), 0, "TPDM <OV|VO>");
    dpd_buf4_scm(&G, -1.0);
    dpd_buf4_close(&G);

    /*
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,V]"), ID("[V,O]"),
                  ID("[O,V]"), ID("[V,O]"), 0, "TPDM <OV|VO>");
    dpd_buf4_print(&G, outfile, 1);
    dpd_buf4_close(&G);
    */

    // G_iabj = -Giajb so I do not need to ovvo block,
    //  however for proper contraction (to avoid construction of <ov||vv>)  in the GFock.cc I need it. 
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
    dpd_buf4_sort(&G, PSIF_OMP2_DENSITY , pqsr, ID("[o,v]"), ID("[v,o]"), "TPDM <ov|vo>");
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,v]"), ID("[v,o]"),
                  ID("[o,v]"), ID("[v,o]"), 0, "TPDM <ov|vo>");
    dpd_buf4_scm(&G, -1.0);
    dpd_buf4_close(&G);
 
    /*
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,v]"), ID("[v,o]"),
                  ID("[o,v]"), ID("[v,o]"), 0, "TPDM <ov|vo>");
    dpd_buf4_print(&G, outfile, 1);
    dpd_buf4_close(&G);
    */

    
    psio_->close(PSIF_OMP2_DENSITY, 1);
   
}
}} // End Namespaces

