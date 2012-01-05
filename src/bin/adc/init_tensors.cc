#include "psi4-dec.h"
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libmints/mints.h>
#include <libtrans/mospace.h>
#include "ccfiles.h"
#include "adc.h"

namespace psi{ namespace adc {
    
double
ADC::rhf_init_tensors()
{
    bool do_pr;
    double ePR2, sq_norm, energy;
    dpdbuf4 K, V;
    dpdfile2 A;
    // Setting up and initialize the integraltransform object 
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);printf("madeok\n");
    _ints = new IntegralTransform(reference_wavefunction_, spaces, IntegralTransform::Restricted);printf("madeok\n");
    _ints->set_keep_iwl_so_ints(true);
    _ints->set_keep_dpd_so_ints(true);
    dpd_set_default(_ints->get_dpd_id());
    // Make (OV|OV) integrals
    fprintf(outfile, "\n\t==> Transforming (OV|OV) Integrals <==\n");
    _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
    // Make (OO|VV) integrals
    fprintf(outfile, "\n\t==> Transforming (OO|VV) Integrals <==\n");
    _ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir);
    // Make (OO|OV) integrals
    fprintf(outfile, "\n\t==> Transforming (OV|OO) Integrals <==\n");
    _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::occ);
    // Make (OV|VV) integrals
    fprintf(outfile, "\n\t==> Transforming (OV|VV) Integrals <==\n");
    _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::vir, MOSpace::vir);
    fflush(outfile);
    
    // Preparing MP1 amplitudes then calculating MP2 energy 
    // and use of LMO is not considered in this code.
    // In ADC(2) calculation, 2 <ij|ab> - <ij|ba> typed
    // integral list is needed. So making this too.
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    psio_->open(PSIF_ADC, PSIO_OPEN_NEW);
    
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    dpd_buf4_sort(&V, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    dpd_buf4_sort(&V, PSIF_LIBTRANS_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "MO Ints V1243");
    dpd_buf4_scmcopy(&V, PSIF_LIBTRANS_DPD, "MO Ints 2 V1234 - V1243", 2.0);
    dpd_buf4_copy(&V, PSIF_ADC, "K1234");
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "MO Ints 2 V1234 - V1243");
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "MO Ints V1243)");
    dpd_buf4_axpy(&V, &K, -1);
    dpd_buf4_close(&V);
    dpd_buf4_close(&K);

    dpd_buf4_init(&K, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "K1234");
    
    for(int h = 0;h < nirrep_;h++){
        dpd_buf4_mat_irrep_init(&K, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        for(int ij = 0;ij < K.params->rowtot[h];ij++){
            int i = K.params->roworb[h][ij][0];
            int j = K.params->roworb[h][ij][1];
            for(int ab = 0;ab < K.params->coltot[h];ab++){
                int a = K.params->colorb[h][ab][0];
                int b = K.params->colorb[h][ab][1];
                K.matrix[h][ij][ab] /= aocce_[i] + aocce_[j] - avire_[a] - avire_[b];
            }
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
    }
    dpd_buf4_sort(&K, PSIF_ADC, pqsr, ID("[O,O]"), ID("[V,V]"), "K1243");
    dpd_buf4_scmcopy(&K, PSIF_ADC, "2 K1234 - K1243", 2.0);
    dpd_buf4_close(&K);
    dpd_buf4_init(&K, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "2 K1234 - K1243");
    dpd_buf4_init(&V, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "K1243");
    dpd_buf4_axpy(&V, &K, -1);
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    energy = dpd_buf4_dot(&V, &K);
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "K1234");
    sq_norm = 1 + dpd_buf4_dot(&V, &K);
    
    do_pr = options_.get_bool("PR");
    
    fprintf(outfile, "\n\t==> Ground State <==\n");
    if(!do_pr) fprintf(outfile, "->");
    fprintf(outfile, "\tMP2 energy    = %20.14f\n", energy);
    fprintf(outfile, "\t[Squared-norm of MP1 wavefunction    = %10.7f]\n", sq_norm);
    // Partially renormalized MP2 energy and the MP1 wavefunction
    // Reference: IJQC 78 (2000) 226, CPL 443 (2007) 389.
    dpd_file2_init(&A, PSIF_ADC, 0, ID('O'), ID('O'), "RHO_OO");
    dpd_contract442(&V, &K, &A, 0, 0, 1, 0);
    dpd_buf4_close(&K);
    fflush(outfile);
    
    dpd_file2_mat_init(&A);
    dpd_file2_mat_rd(&A);
    dpd_buf4_init(&K, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "tilde K1234");
    
    for(int h = 0;h < nirrep_;h++){
        dpd_buf4_mat_irrep_init(&V, h);
        dpd_buf4_mat_irrep_init(&K, h);
        dpd_buf4_mat_irrep_rd(&V, h);
        for(int ij = 0;ij < V.params->rowtot[h];ij++){
            int i = V.params->roworb[h][ij][0];
            int j = V.params->roworb[h][ij][1];
            int isym = V.params->psym[i];
            int jsym = V.params->qsym[j];
            int I = i - V.params->poff[isym];
            int J = j - V.params->qoff[jsym];
            double Nij = 1 + (A.matrix[isym][I][I] + A.matrix[jsym][J][J]) / 2;
            for(int ab = 0;ab < V.params->coltot[h];ab++)
                K.matrix[h][ij][ab] = V.matrix[h][ij][ab] / Nij;
        }
        dpd_buf4_mat_irrep_wrt(&K, h);
        dpd_buf4_mat_irrep_close(&K, h);
        dpd_buf4_mat_irrep_close(&V, h);
    }
    dpd_buf4_close(&V);
    dpd_file2_mat_close(&A);
    dpd_file2_close(&A);

    dpd_buf4_sort(&K, PSIF_ADC, pqsr, ID("[O,O]"), ID("[V,V]"), "tilde K1243");
    dpd_buf4_scmcopy(&K, PSIF_ADC, "tilde 2 K1234 - K1243", 2.0);
    dpd_buf4_close(&K);
    dpd_buf4_init(&K, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "tilde 2 K1234 - K1243");
    dpd_buf4_init(&V, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "tilde K1243");
    dpd_buf4_axpy(&V, &K, -1.0);
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    ePR2 = dpd_buf4_dot(&V, &K);
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "tilde K1234");
    sq_norm = 1 + dpd_buf4_dot(&V, &K);
    dpd_buf4_close(&K);
    dpd_buf4_close(&V);
    
    if(do_pr) fprintf(outfile, "->");
    fprintf(outfile, "\tPR-MP2 energy = %20.14f\n", ePR2);
    fprintf(outfile, "\t[Squared-norm of PR-MP1 wavefunction = %10.7f]\n\n", sq_norm);
    fflush(outfile);
    
    if(options_.get_bool("PR")) energy = ePR2;
    
    // Reordering each ERIs other than (OO|VV) type from Mulliken to Dirac notation 
    // for convenience of the evaluation of the sigma tensor
    
    // Sort(prqs): <OV|OV> <-- (OO|VV) 
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
    dpd_buf4_sort(&V, PSIF_LIBTRANS_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "MO Ints <OV|OV>");
    dpd_buf4_close(&V);
    
    // Sort(pqrs): <OO|VO> <-- (OV|OO)
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"), ID("[O,V]"), ID("[O>=O]+"), 0, "MO Ints (OV|OO)");
    dpd_buf4_sort(&V, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,O]"), "MO Ints <OO|VO>");
    dpd_buf4_close(&V);

    // Sort(prqs): <OV|VV> <-- (OV|VV)
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
    dpd_buf4_sort(&V, PSIF_LIBTRANS_DPD, prqs, ID("[O,V]"), ID("[V,V]"), "MO Ints <OV|VV>");
    dpd_buf4_close(&V);

    return energy;
} 
    
}} // End Namespaces
