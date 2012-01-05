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
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libtrans/mospace.h>
#include <libtrans/integraltransform.h>

/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>

#include "defines.h"
#include "omp2wave.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace omp2wave{
  
void OMP2Wave::t2_1st_sc()
{   
     dpdbuf4 K, T, D;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OMP2_DPD, PSIO_OPEN_OLD);
     
     // Build T2AA
     // T_IJ^AB = <IJ||AB>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
    dpd_buf4_copy(&K, PSIF_OMP2_DPD, "T2_1 <OO|VV>");
    dpd_buf4_close(&K);
    
    
    // T_IJ^AB = T_IJ^AB / D_IJ^AB
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_dirprd(&D, &T);
    dpd_buf4_close(&D);
    if (print_ > 1) dpd_buf4_print(&T, outfile, 1);
    dpd_buf4_close(&T);
    
    
    // Build T2BB
    // T_ij^ab = <ij|ab>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
    dpd_buf4_copy(&K, PSIF_OMP2_DPD, "T2_1 <oo|vv>");
    dpd_buf4_close(&K);
    
    
    // T_ij^ab = T_ij^ab / D_ij^ab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    dpd_buf4_dirprd(&D, &T);
    dpd_buf4_close(&D);
    if (print_ > 1) dpd_buf4_print(&T, outfile, 1);
    dpd_buf4_close(&T);
    
    
    // Build T2AB
    // T_Ij^Ab = <Ij|Ab>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_copy(&K, PSIF_OMP2_DPD, "T2_1 <Oo|Vv>");
    dpd_buf4_close(&K);
    
    
    // T_Ij^Ab = T_Ij^Ab / D_Ij^Ab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_dirprd(&D, &T);
    dpd_buf4_close(&D);
    if (print_ > 1) dpd_buf4_print(&T, outfile, 1);
    dpd_buf4_close(&T);
     
     
     psio_->close(PSIF_LIBTRANS_DPD, 1);
     psio_->close(PSIF_OMP2_DPD, 1);

} // end t2_1st_sc


void OMP2Wave::t2_1st_general()
{   
     //fprintf(outfile,"\n t2_1st_general is starting... \n"); fflush(outfile);
     
     dpdbuf4 K, T, Tnew, D, R;
     dpdfile2 Fo,Fv;
     int nElements;
     
     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     psio_->open(PSIF_OMP2_DPD, PSIO_OPEN_OLD);
     
     
    // Build new T2AA 
    // T_IJ^AB = <IJ||AB>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
    dpd_buf4_copy(&K, PSIF_OMP2_DPD, "T2_1new <OO|VV>");
    dpd_buf4_close(&K);
    
    
    // initalize Tnew and Told
    dpd_buf4_init(&Tnew, PSIF_OMP2_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1new <OO|VV>");
    dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    
    
    // T_IJ^AB = \sum_{E} T_IJ^AE * F_EB + \sum_{E} T_IJ^EB * F_AE
    dpd_file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");  
    dpd_contract424(&T, &Fv, &Tnew, 3, 1, 0, 1.0, 1.0); 
    dpd_contract244(&Fv, &T, &Tnew, 1, 2, 1, 1.0, 1.0); 
    dpd_file2_close(&Fv);
    
    // T_IJ^AB = -\sum_{M} T_IM^AB * F_MJ - \sum_{M} T_MJ^AB * F_MI
    dpd_file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    dpd_contract424(&T, &Fo, &Tnew, 1, 0, 1, -1.0, 1.0);
    dpd_contract244(&Fo, &T, &Tnew, 0, 0, 0, -1.0, 1.0);
    dpd_file2_close(&Fo);
    dpd_buf4_close(&T);
    
  
    // T_IJ^AB = T_IJ^AB / D_IJ^AB
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    dpd_buf4_dirprd(&D, &Tnew);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Tnew);
    
    
    
    
    // Build new T2BB 
    // T_ij^ab = <ij||ab>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
    dpd_buf4_copy(&K, PSIF_OMP2_DPD, "T2_1new <oo|vv>");
    dpd_buf4_close(&K);
    
    
    // initalize Tnew and Told
    dpd_buf4_init(&Tnew, PSIF_OMP2_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1new <oo|vv>");
    dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    
    
    // T_ij^ab = \sum_{e} T_ij^ae * F_eb + \sum_{e} T_ij^eb * F_ae
    dpd_file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");  
    dpd_contract424(&T, &Fv, &Tnew, 3, 1, 0, 1.0, 1.0); 
    dpd_contract244(&Fv, &T, &Tnew, 1, 2, 1, 1.0, 1.0); 
    dpd_file2_close(&Fv);
    
    // T_ij^ab = -\sum_{m} T_im^ab * F_mj - \sum_{m} T_mj^ab * F_mi
    dpd_file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    dpd_contract424(&T, &Fo, &Tnew, 1, 0, 1, -1.0, 1.0);
    dpd_contract244(&Fo, &T, &Tnew, 0, 0, 0, -1.0, 1.0);
    dpd_file2_close(&Fo);
    dpd_buf4_close(&T);
    
  
    // T_ij^ab = T_ij^ab / D_ij^ab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    dpd_buf4_dirprd(&D, &Tnew);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Tnew);
    
    
    
    // Build new T2AB
    // T_Ij^Ab = <Ij||Ab>
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_copy(&K, PSIF_OMP2_DPD, "T2_1new <Oo|Vv>");
    dpd_buf4_close(&K);
    
    
    // initalize Tnew and Told
    dpd_buf4_init(&Tnew, PSIF_OMP2_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1new <Oo|Vv>");
    dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    
    
    // T_Ij^Ab = \sum_{e} T_Ij^Ae * F_be + \sum_{E} T_Ij^Eb * F_AE
    dpd_file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");  
    dpd_contract424(&T, &Fv, &Tnew, 3, 1, 0, 1.0, 1.0); 
    dpd_file2_close(&Fv);
    dpd_file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");  
    dpd_contract244(&Fv, &T, &Tnew, 1, 2, 1, 1.0, 1.0); 
    dpd_file2_close(&Fv);
    
    // T_Ij^Ab = -\sum_{m} T_Im^Ab * F_mj - \sum_{M} T_Mj^Ab * F_MI
    dpd_file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    dpd_contract424(&T, &Fo, &Tnew, 1, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Fo);
    dpd_file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    dpd_contract244(&Fo, &T, &Tnew, 0, 0, 0, -1.0, 1.0);
    dpd_file2_close(&Fo);
    dpd_buf4_close(&T);
    
  
    // T_Ij^Ab = T_Ij^Ab / D_Ij^Ab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    dpd_buf4_dirprd(&D, &Tnew);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Tnew);
    
    
    
    
    // Compute amplitude residual to Check Convergence
    // Alpha-Alpha spin case
    dpd_buf4_init(&Tnew, PSIF_OMP2_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1new <OO|VV>");
    dpd_buf4_copy(&Tnew, PSIF_OMP2_DPD, "RT2_1 <OO|VV>");
    dpd_buf4_init(&R, PSIF_OMP2_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "RT2_1 <OO|VV>");
    dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_axpy(&T, &R, -1.0); // -1.0*T + R -> R
    dpd_buf4_close(&T);
    
    nElements = 0;
    for(int h = 0; h < nirreps; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
    rms_t2AA = 0.0;
    rms_t2AA = dpd_buf4_dot_self(&R);
    dpd_buf4_close(&R);
    rms_t2AA = sqrt(rms_t2AA) / nElements;
    
    // Reset
    dpd_buf4_copy(&Tnew, PSIF_OMP2_DPD, "T2_1 <OO|VV>");
    if (print_ > 1) dpd_buf4_print(&Tnew, outfile, 1);
    dpd_buf4_close(&Tnew);
    
    
    // Beta-Beta spin case
    dpd_buf4_init(&Tnew, PSIF_OMP2_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1new <oo|vv>");
    dpd_buf4_copy(&Tnew, PSIF_OMP2_DPD, "RT2_1 <oo|vv>");
    dpd_buf4_init(&R, PSIF_OMP2_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "RT2_1 <oo|vv>");
    dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    dpd_buf4_axpy(&T, &R, -1.0); // -1.0*T + R -> R
    dpd_buf4_close(&T);
    
    nElements = 0;
    for(int h = 0; h < nirreps; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
    rms_t2BB = 0.0;
    rms_t2BB = dpd_buf4_dot_self(&R);
    dpd_buf4_close(&R);
    rms_t2BB = sqrt(rms_t2BB) / nElements;
    
    // Reset
    dpd_buf4_copy(&Tnew, PSIF_OMP2_DPD, "T2_1 <oo|vv>");
    if (print_ > 1) dpd_buf4_print(&Tnew, outfile, 1);
    dpd_buf4_close(&Tnew);
    
    
    // Alpha-Beta spin case
    dpd_buf4_init(&Tnew, PSIF_OMP2_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1new <Oo|Vv>");
    dpd_buf4_copy(&Tnew, PSIF_OMP2_DPD, "RT2_1 <Oo|Vv>");
    dpd_buf4_init(&R, PSIF_OMP2_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "RT2_1 <Oo|Vv>");
    dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_axpy(&T, &R, -1.0); // -1.0*T + R -> R
    dpd_buf4_close(&T);
    
    nElements = 0;
    for(int h = 0; h < nirreps; h++) nElements += R.params->coltot[h] * R.params->rowtot[h];
    rms_t2AB = 0.0;
    rms_t2AB = dpd_buf4_dot_self(&R);
    dpd_buf4_close(&R);
    rms_t2AB = sqrt(rms_t2AA) / nElements;
    
    // Reset
    dpd_buf4_copy(&Tnew, PSIF_OMP2_DPD, "T2_1 <Oo|Vv>");
    if (print_ > 1) dpd_buf4_print(&Tnew, outfile, 1);
    dpd_buf4_close(&Tnew);
    
    
    
    // close files
    psio_->close(PSIF_LIBTRANS_DPD, 1);
    psio_->close(PSIF_OMP2_DPD, 1);
    
    //fprintf(outfile,"\n t2_1st_general done. \n"); fflush(outfile);

} // end t2_1st_general

}} // End Namespaces

