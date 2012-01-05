#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <libciomr/libciomr.h>
#include <libmints/mints.h>
#include "physconst.h"
#include "adc.h"

//
//  A3h3p: The frequency independent terms in the response matrix, composed of the CIS part and 3hole-3particle
//         diagrams which give so-called the diffenential correlation effect in terms of CIS(D) sense, so with  
//         OVOV indices.
//  AOO  : Tensor with OO indices that represents symmetrized contribution from a 3h3p diagram that can be
//         evaluated.
//         independently from the trial vector in the SEM procedure.
//  AVV  : Tensor with OO indices that represents the time reversed contribution from the above one.
//  D    : Diagonal elements packed in DPD fashon of two index tensors, which is used in diagonalization step.
//

namespace psi{ namespace adc{
    
void 
ADC::rhf_prepare_tensors()
{    
    char lbl[32];
    double *omega, **lambda;
    dpdbuf4 Aovov, K, V;
    dpdfile2 Xoo, Xvv, Aoo, Avv, Dov, Cocc, Cvir, B;
    
    fprintf(outfile, "\t==> CIS/ADC(1) Level <==\n\n");
    psio_->open(PSIF_ADC_SEM, PSIO_OPEN_NEW);
    
    // CIS calculation for obtaining the guess energy and vector
    // for the second order calculation. For this purpose, the external
    // exchange operator (EEO) method a.k.a Fock-like contraction
    // is not utilized because the construction of sigma-vector for ADC(2) is
    // far more demanding step.
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    dpd_buf4_scmcopy(&V, PSIF_ADC, "A1234", -1.0);
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    dpd_buf4_sort_axpy(&V, PSIF_ADC, prqs, ID("[O,V]"), ID("[O,V]"), "A1234", 2.0);
    dpd_buf4_close(&V);
    dpd_buf4_init(&Aovov, PSIF_ADC, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "A1234");
    char **irrep_      = Process::environment.molecule()->irrep_labels();
    for(int h = 0;h < nirrep_;h++){
        dpd_buf4_mat_irrep_init(&Aovov, h);
        dpd_buf4_mat_irrep_rd(&Aovov, h);
        for(int ia = 0;ia < Aovov.params->rowtot[h];ia++){
            int i = Aovov.params->roworb[h][ia][0];
            int a = Aovov.params->colorb[h][ia][1];
            for(int jb = 0;jb < Aovov.params->coltot[h];jb++){
                int j = Aovov.params->colorb[h][jb][0];
                int b = Aovov.params->colorb[h][jb][1];
                Aovov.matrix[h][ia][jb] += (avire_[a]-aocce_[i]) * (i==j) * (a==b);
            }
        }
        omega = init_array(Aovov.params->rowtot[h]);
        lambda = block_matrix(Aovov.params->rowtot[h], rpi_[h]);
        if(rpi_[h]) david(Aovov.matrix[h], Aovov.params->coltot[h], rpi_[h], omega, lambda, 1e-14, 0);
        for(int root = 0;root < rpi_[h];root++){
            if(DEBUG_) printf("%d%3s, %10.7f\n", root+1, irrep_[h], omega[root]);
            omega_guess_->set(h,root, omega[root]);
            sprintf(lbl, "B^(%d)_[%d]12", root, h);
            dpd_file2_init(&B, PSIF_ADC, h, ID('O'), ID('V'), lbl);
            dpd_file2_mat_init(&B);
            for(int ia = 0;ia < Aovov.params->rowtot[h];ia++){
                int i = Aovov.params->roworb[h][ia][0];
                int a = Aovov.params->roworb[h][ia][1];
                
                int I = B.params->rowidx[i];
                int A = B.params->colidx[a];
                int Isym = B.params->psym[i];
                B.matrix[Isym][I][A] = lambda[ia][root];
            }
            dpd_file2_mat_wrt(&B);
            dpd_file2_mat_close(&B);
            fprintf(outfile, "\t%d%3s state: %10.7f (a.u.), %10.7f (eV)\n", root+1, irrep_[h], omega[root], omega[root]*_hartree2ev);
            fprintf(outfile, "\t---------------------------------------------\n");
            int nprint;
            if(nxspi_[h] < num_amps_) nprint = nxspi_[h];
            else nprint = num_amps_;
            amps_write(&B, nprint, outfile);
            fprintf(outfile, "\n");
            dpd_file2_close(&B);
        }
        
        free(omega);
        free_block(lambda);
        dpd_buf4_mat_irrep_wrt(&Aovov, h);
        dpd_buf4_mat_irrep_close(&Aovov, h);
    }
    dpd_buf4_close(&Aovov);
    fflush(outfile);

    // Initialize 4-index tensors for this step of calculation.
    dpd_buf4_init(&Aovov, PSIF_ADC, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "A1234");
    dpd_buf4_copy(&Aovov, PSIF_ADC_SEM, "A3h3p1234");
    dpd_buf4_close(&Aovov);
    
    dpd_buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    dpd_buf4_init(&K, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "2 K1234 - K1243");
    
    dpd_file2_init(&Xoo, PSIF_ADC_SEM, 0, ID('O'), ID('O'), "XOO12");
    // XOO_{ij} <-- - \sum_{kab} (2K_{ikab} - K_{ikba}) <jk|ab>
    dpd_contract442(&K, &V, &Xoo, 0, 0, -1, 0);
    dpd_file2_init(&Aoo, PSIF_ADC_SEM, 0, ID('O'), ID('O'), "AOO12");
    // AOO_{ij} <-- ((XOO)_{ij} + (XOO)_{ji}) / 2
    dpd_file2_axpy(&Xoo, &Aoo, 0.5, 0);
    dpd_file2_axpy(&Xoo, &Aoo, 0.5, 1);
    dpd_file2_mat_init(&Aoo);
    dpd_file2_mat_rd(&Aoo);
    dpd_file2_close(&Xoo);
    
    dpd_file2_init(&Xvv, PSIF_ADC_SEM, 0, ID('V'), ID('V'), "XVV12");
    // XVV_{ab} <-- - \sum_{ijc} (2K_{ijac} - K_{ijca}) <ij|bc>
    dpd_contract442(&K, &V, &Xvv, 2, 2, -1, 0);
    dpd_file2_init(&Avv, PSIF_ADC_SEM, 0, ID('V'), ID('V'), "AVV12");
    // AVV_{ab} <-- ((XVV)_{ab} + (XVV)_{ba}) / 2
    dpd_file2_axpy(&Xvv, &Avv, 0.5, 0);
    dpd_file2_axpy(&Xvv, &Avv, 0.5, 1);
    dpd_file2_mat_init(&Avv);
    dpd_file2_mat_rd(&Avv);
    dpd_file2_close(&Xvv);
    
    dpd_buf4_close(&K);
    dpd_buf4_close(&V);

    // A3h3p_{iajb} <-- \delta_{ij}(XVV)_{ab} + \delta_{ab}(XOO)_{ij}
    dpd_buf4_init(&Aovov, PSIF_ADC_SEM, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "A3h3p1234");
    for(int h = 0;h < nirrep_;h++){
        dpd_buf4_mat_irrep_init(&Aovov, h);
        dpd_buf4_mat_irrep_rd(&Aovov, h);
        for(int ia = 0;ia < Aovov.params->rowtot[h];ia++){
            int i = Aovov.params->roworb[h][ia][0];
            int a = Aovov.params->roworb[h][ia][1];
            int I = Aoo.params->rowidx[i];
            int A = Avv.params->rowidx[a];
            int Isym = Aoo.params->psym[i];
            int Asym = Avv.params->psym[a];
            for(int jb = 0;jb < Aovov.params->coltot[h];jb++){
                int j = Aovov.params->colorb[h][jb][0];
                int b = Aovov.params->colorb[h][jb][1];
                int J = Aoo.params->colidx[j];
                int B = Avv.params->colidx[b];
                int Jsym = Aoo.params->qsym[j];
                int Bsym = Avv.params->qsym[b];
                Aovov.matrix[h][ia][jb] += Aoo.matrix[Isym][I][J] * (a==b) * (Isym==Jsym)
                    + Avv.matrix[Asym][A][B] * (i==j) * (Asym==Bsym);
            }
        }
        dpd_buf4_mat_irrep_wrt(&Aovov, h);
        dpd_buf4_mat_irrep_close(&Aovov, h);
    }
    dpd_buf4_close(&Aovov);
    dpd_file2_mat_close(&Aoo);
    dpd_file2_close(&Aoo);
    dpd_file2_mat_close(&Avv);
    dpd_file2_close(&Avv);
    
    // Preparing D tensor for each irrep
    dpd_buf4_init(&Aovov, PSIF_ADC_SEM, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "A3h3p1234");
    for(int h = 0;h < nirrep_;h++){
        dpd_buf4_mat_irrep_init(&Aovov, h);
        dpd_buf4_mat_irrep_rd(&Aovov, h);
        sprintf(lbl, "D_[%d]12", h);
        dpd_file2_init(&Dov, PSIF_ADC_SEM, h, ID('O'), ID('V'), lbl);
        dpd_file2_mat_init(&Dov);
        for(int ia = 0;ia < Aovov.params->rowtot[h];ia++){
            int i = Aovov.params->roworb[h][ia][0];
            int a = Aovov.params->roworb[h][ia][1];
            
            int I = Dov.params->rowidx[i];
            int A = Dov.params->colidx[a];
            int Isym = Dov.params->psym[i];
            Dov.matrix[Isym][I][A] = Aovov.matrix[h][ia][ia];
        }
        dpd_file2_mat_wrt(&Dov);
        dpd_file2_mat_close(&Dov);
        dpd_file2_close(&Dov);
        dpd_buf4_mat_irrep_close(&Aovov, h);
    }
    dpd_buf4_close(&Aovov);

    psio_->close(PSIF_ADC, 1);
    psio_->close(PSIF_ADC_SEM, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

}} // End Namespaces
