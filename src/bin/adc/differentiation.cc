#include "psi4-dec.h"
#include <libtrans/integraltransform.h>
#include "adc.h"

namespace psi{ namespace adc{
    
//
//  V     : The converged eigenvector in the S manifold.
//  D     : This gives contribution such that \frac{\partial A(\omega)}{\partial \omega}V(\omega).
//  XOVOV : A 2h-2p intermediate.
//  YOVOV : A 2h-2p intermediate.
//  ZOVOV : A 2h-2p intermediate.
//  BOVOV : A 2h-2p intermediate.
//
    
double 
ADC::rhf_differentiate_omega(int irrep, int root)
{
    char lbl[32];
    double dot;
    dpdfile2 S, D;
    dpdbuf4 A, V, K, Z;

    sprintf(lbl, "V^(%d)_[%d]12", root, irrep);
    dpd_file2_init(&S, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
    sprintf(lbl, "D^(%d)_[%d]12", root, irrep);
    dpd_file2_init(&D, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
    
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    sprintf(lbl, "ZOOVV_[%d]1234", irrep);
    dpd_buf4_init(&Z, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    // ZOVOV_{jiab} <--   \sum_{c} <jc|ab> v_{ic}
    dpd_contract424(&V, &S, &Z, 1, 1, 1,  1, 0);
    dpd_buf4_close(&V);
        
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,O]"), ID("[O,O]"), ID("[V,O]"), 0, "MO Ints <OO|VO>");
    // ZOVOV_{ijab} <-- - \sum_{k} <ij|ak> v_{kb}
    dpd_contract424(&V, &S, &Z, 3, 0, 0, -1, 1);
    dpd_buf4_close(&V);
        
    // dB_{iajb} <-- - (2Z_{ijab}-Z_{ijba}+2Z_{jiab}-Z_{jiba}) / (\omega+e_i-e_a+e_j-e_b)^2
    sprintf(lbl, "BOOVV_[%d]1234", irrep);
    dpd_buf4_scmcopy(&Z, PSIF_ADC_SEM, lbl, -2.0);
    dpd_buf4_sort_axpy(&Z, PSIF_ADC_SEM, pqsr, ID("[O,O]"), ID("[V,V]"), lbl,  1.0);
    dpd_buf4_sort_axpy(&Z, PSIF_ADC_SEM, qprs, ID("[O,O]"), ID("[V,V]"), lbl,  1.0);
    dpd_buf4_sort_axpy(&Z, PSIF_ADC_SEM, qpsr, ID("[O,O]"), ID("[V,V]"), lbl, -2.0);
    dpd_buf4_close(&Z);
    
    dpd_buf4_init(&Z, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    sprintf(lbl, "D_[%d]1234", irrep);
    dpd_buf4_init(&A, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    dpd_buf4_dirprd(&A, &Z);
    dpd_buf4_dirprd(&A, &Z);
    dpd_buf4_close(&A);
    
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    // \dV_{ia} <-- \sum_{jbc} B_{jicb} <ja|cb>
    dpd_contract442(&Z, &V, &D, 1, 1, 1, 0);
    dpd_buf4_close(&V);
    
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,O]"), ID("[O,O]"), ID("[V,O]"), 0, "MO Ints <OO|VO>");
    // \dV_{ia} <-- - \sum_{jkb} <kj|bi> B_{jkab}
    dpd_contract442(&V, &Z, &D, 3, 3, -1, 0); //This is genuine
    dpd_buf4_close(&V);
    dpd_buf4_close(&Z);

    // \frac{\partial \omega^{eigen}}{\partial omega} = V^t\frac{\partial A(\omega)}{\partial \omega}V
    dot = dpd_file2_dot(&S, &D);
    
    dpd_file2_close(&S);
    dpd_file2_close(&D);
        
    return dot;
}

}} // End Namespaces
