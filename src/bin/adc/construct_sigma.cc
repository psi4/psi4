#include "psi4-dec.h"
#include <libtrans/integraltransform.h>
#include "adc.h"

namespace psi{ namespace adc{
    
//
//  DOV   : Tensors that contain the contribution from a 3h-3p diagram that cannot be 
//          be evaluated independently from the trial vector, B.
//  EOV   : An intermediate introduced in order to retain the symmetricity of the response matrix.
//  XOOVV : A 2h-2p intermediate.
//  YOOVV : A 2h-2p intermediate.
//  ZOOVV : A 2h-2p intermediate.
//  BOOVV : A 2h-2p intermediate.
//
    
void 
ADC::rhf_construct_sigma(int irrep, int root)
{
    char lbl[32];
    dpdfile2 B, S, D, E, Bt, C;
    dpdbuf4 A, V, K, X, Y, Z, BT, XT;
            
    sprintf(lbl, "S^(%d)_[%d]12", root, irrep);
    dpd_file2_init(&S, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
    sprintf(lbl, "B^(%d)_[%d]12", root, irrep);
    dpd_file2_init(&B, PSIF_ADC,     irrep, ID('O'), ID('V'), lbl);
    
    // CIS term and the two 3h-3p diagrams are summed into the sigma vector.
    dpd_buf4_init(&A, PSIF_ADC_SEM, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "A3h3p1234");
    dpd_contract422(&A, &B, &S, 0, 0, 1, 0);
    dpd_buf4_close(&A);
    
    // Evaluation of the remaining one 3h-3p diagram
    dpd_buf4_init(&K, PSIF_ADC,          0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "2 K1234 - K1243");
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "MO Ints 2 V1234 - V1243");
    
    sprintf(lbl, "DOV_[%d]12", irrep);
    dpd_file2_init(&D, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
    // D_{ia} <-- \sum_{jb} (2 <ij|ab> - <ij|ba>) b_{jb}
    dpd_dot24(&B, &V, &D, 0, 0, 1, 0);
    // \sigma_{ia} <-- 0.5 \sum_{jb} (2 K_{ijab} - K_{ijba}) D_{jb}
    dpd_dot24(&D, &K, &S, 0, 0, 0.5, 1);
    dpd_file2_close(&D);
    
    sprintf(lbl, "EOV_[%d]12", irrep);
    dpd_file2_init(&E, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
    // E_{ia} <-- \sum_{jb} (2 K_{ijab} - K_{ijba}) b_{jb}
    dpd_dot24(&B, &K, &E, 0, 0, 1, 0);
    // \sigma_{ia} <-- \sum_{jb} (2 <ij|ab> - <ij|ba>) E_{jb}
    dpd_dot24(&E, &V, &S, 0, 0, 0.5, 1);
    dpd_file2_close(&E);
    
    dpd_buf4_close(&K);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    sprintf(lbl, "XOOVV_[%d]1234", irrep);
    dpd_buf4_init(&X, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    // XOVOV_{jiab} <-- \sum_{c} <jc|ab> b_{ic}
    dpd_contract424(&V, &B, &X, 1, 1, 1, 1, 0);
    dpd_buf4_close(&V);
    
    sprintf(lbl, "ZOOVV_[%d]1234", irrep);
    dpd_buf4_copy(&X, PSIF_ADC_SEM, lbl);
    dpd_buf4_close(&X);

    sprintf(lbl, "YOOVV_[%d]1234", irrep);
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,O]"), ID("[O,O]"), ID("[V,O]"), 0, "MO Ints <OO|VO>");
    dpd_buf4_init(&Y, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    // YOVOV_{ijab} <-- \sum_{k} <ij|ak> b_{kb}
    dpd_contract424(&V, &B, &Y, 3, 0, 0, 1, 0);
    dpd_buf4_close(&V);
    
    sprintf(lbl, "ZOOVV_[%d]1234", irrep);
    dpd_buf4_init(&Z, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    // Z_{ijab} <-- X_{jiab} - Y_{ijab}
    dpd_buf4_axpy(&Y, &Z, -1.0);
    dpd_buf4_close(&Y);

    // B_{iajb} <-- (2Z_{ijab}-Z_{ijba}+2Z_{jiab}-Z_{jiba}) / (\omega+e_i-e_a+e_j-e_b)
    sprintf(lbl, "BOOVV_[%d]1234", irrep);
    dpd_buf4_scmcopy(&Z, PSIF_ADC_SEM, lbl, 2.0);
    dpd_buf4_sort_axpy(&Z, PSIF_ADC_SEM, pqsr, ID("[O,O]"), ID("[V,V]"), lbl, -1.0);
    dpd_buf4_sort_axpy(&Z, PSIF_ADC_SEM, qprs, ID("[O,O]"), ID("[V,V]"), lbl, -1.0);
    dpd_buf4_sort_axpy(&Z, PSIF_ADC_SEM, qpsr, ID("[O,O]"), ID("[V,V]"), lbl,  2.0);
    dpd_buf4_close(&Z);
    
    dpd_buf4_init(&Z, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    sprintf(lbl, "D_[%d]1234", irrep);
    dpd_buf4_init(&A, PSIF_ADC_SEM, irrep, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, lbl);
    dpd_buf4_dirprd(&A, &Z);
    dpd_buf4_close(&A);
 
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    // \Sigma_{ia} <-- \sum_{jbc} B_{jicb} <ja|cb>
    dpd_contract442(&Z, &V, &S, 1, 1, 1, 1);
    dpd_buf4_close(&V);
    
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,O]"), ID("[O,O]"), ID("[V,O]"), 0, "MO Ints <OO|VO>");
    // \sigma_{ia} <-- - \sum_{jkb} <kj|bi> B_{jkab}
    dpd_contract442(&V, &Z, &S, 3, 3, -1, 1); //This is genuine
    dpd_buf4_close(&V);
    dpd_buf4_close(&Z);

    dpd_file2_close(&S);
    dpd_file2_close(&B);
}

}} // End Namespaces
