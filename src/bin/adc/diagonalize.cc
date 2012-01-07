#include "psi4-dec.h"
#include <liboptions/liboptions.h>
#include <libpsio/psio.h>
#include <libtrans/integraltransform.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <cmath>
#include "adc.h"

namespace psi{ namespace adc{

//
//  This block-Davidson code is based on the in-core version put in lib/libqt/david.cc
//  but written by utilizind DPD algorithm.
//
//  S: Sigma vector for the response matrix with OV indices.
//  F: The correction vectors for the basis of Ritz space.
//  V: Converged eigenvectors.
//
    
void 
ADC::rhf_diagonalize(int irrep, int num_root, bool first, double omega_in, double *eps)
{
    char lbl[32];
    int iter, converged, length, *conv, skip_check, maxdim, *residual_ok;
    double **Alpha, **G, *lambda, *lambda_o, *residual_norm, cutoff;
    dpdfile2 B, S, Bp, Bn, F, L, V;
    dpdfile4 A;
    FILE *iter_adc;
    
    maxdim = 10 * rpi_[irrep];
    iter = 0;
    converged = 0;
    cutoff = pow(10, -(conv_+1));
    length = rpi_[irrep];
    
    residual_ok   = init_int_array(rpi_[irrep]);
    residual_norm = init_array(rpi_[irrep]);
    conv          = init_int_array(rpi_[irrep]);
    
    G        = block_matrix(maxdim, maxdim);
    Alpha    = block_matrix(maxdim, maxdim);
    lambda   = init_array(maxdim);
    lambda_o = init_array(maxdim);
    
    for(int I = 0;I < rpi_[irrep];I++) lambda_o[I] = omega_guess_->get(irrep, I);
    shift_denom4(irrep, omega_in);
    
    ffile(&iter_adc, "iter.dat", 1);
    
    timer_on("SEM");
    while(converged < rpi_[irrep] && iter < sem_max_){
        skip_check = 0;
        fprintf(iter_adc, "\niter = %d, dim = %d\n", iter, length);

        // Evaluating the sigma vectors
        timer_on("Sigma construction");
        for(int I = 0;I < length;I++)
            if(!nopen_) rhf_construct_sigma(irrep, I);
        timer_off("Sigma construction");

        // Making so called Davidson mini-Hamiltonian, or Rayleigh matrix
        for(int I = 0;I < length;I++){
            sprintf(lbl, "S^(%d)_[%d]12", I, irrep);
            dpd_file2_init(&S, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
            for(int J = 0;J <= I;J++){
                sprintf(lbl, "B^(%d)_[%d]12", J, irrep);
                dpd_file2_init(&B, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
                double sum = dpd_file2_dot(&B, &S);
                if(I != J)
                    G[I][J] = G[J][I] = sum;
                else
                    G[I][J] = sum;
                dpd_file2_close(&B);
            }
            dpd_file2_close(&S);
        }
        if(first && !iter)  poles_[irrep][num_root-1].ps_value = G[num_root-1][num_root-1];
        sq_rsp(length, length, G, lambda, 1, Alpha, 1e-12);

        // Constructing the corretion vectors
        for(int k = 0;k < rpi_[irrep];k++){
            sprintf(lbl, "F^(%d)_[%d]12", k, irrep);
            dpd_file2_init(&F, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
            dpd_file2_scm(&F, 0.0);
            for(int I = 0;I < length;I++){
                sprintf(lbl, "B^(%d)_[%d]12", I, irrep);
                dpd_file2_init(&B, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
                dpd_file2_axpy(&B, &F, -Alpha[I][k]*lambda[k], 0);
                dpd_file2_close(&B);
                sprintf(lbl, "S^(%d)_[%d]12", I, irrep);
                dpd_file2_init(&S, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
                dpd_file2_axpy(&S, &F, Alpha[I][k], 0);
                dpd_file2_close(&S);
            }
            shift_denom2(k, irrep, lambda[k]);
            sprintf(lbl, "L^(%d)_[%d]12", k, irrep);
            dpd_file2_init(&L, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
            dpd_file2_dirprd(&L, &F);
            dpd_file2_close(&L);
        
            double norm = dpd_file2_dot_self(&F);
            residual_norm[k] = sqrt(norm);
            if(residual_norm[k] > pow(10, -norm_tol_)) dpd_file2_scm(&F, 1/residual_norm[k]);
            else {
                dpd_file2_scm(&F, 0.0);
                residual_ok[k] = 1;
            }
            dpd_file2_close(&F);
        }

        // Expand the Ritz space by orthogonalizing {F} to {B} according to Gram-Schmidt procedure
        for(int k = 0;k < rpi_[irrep];k++){
            sprintf(lbl, "F^(%d)_[%d]12", k, irrep);
            dpd_file2_init(&F, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
            sprintf(lbl, "Bpp_[%d]12", irrep);
            dpd_file2_copy(&F, PSIF_ADC_SEM, lbl);
            for(int I = 0;I < length;I++){
                sprintf(lbl, "B^(%d)_[%d]12", I, irrep);
                dpd_file2_init(&B, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
                double coeff = - dpd_file2_dot(&F, &B);
                sprintf(lbl, "Bpp_[%d]12", irrep);
                dpd_file2_init(&Bp, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
                dpd_file2_axpy(&B, &Bp, coeff, 0.0);
                dpd_file2_close(&B);
                dpd_file2_close(&Bp);
            }
            dpd_file2_close(&F);
            sprintf(lbl, "Bpp_[%d]12", irrep);
            dpd_file2_init(&B, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
            double norm = dpd_file2_dot_self(&B);
            norm = sqrt(norm);
        
            if(norm > pow(10, -norm_tol_)){
                dpd_file2_scm(&B, 1/norm);
                sprintf(lbl, "B^(%d)_[%d]12", length, irrep);
                dpd_file2_copy(&B, PSIF_ADC, lbl);
                length++;
            }
            dpd_file2_close(&B);
        }
    
        if(maxdim-length < rpi_[irrep] || (nxspi_[irrep]-length) < rpi_[irrep]){
            fprintf(iter_adc, "Subspace too large:maxdim = %d, L = %d\n", maxdim, length);
            fprintf(iter_adc, "Collapsing eigenvectors.\n");

            for(int k = 0;k < rpi_[irrep];k++){
                sprintf(lbl, "Bn^(%d)_[%d]12", k, irrep);
                dpd_file2_init(&Bn, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
                dpd_file2_scm(&Bn, 0.0);
                for(int I = 0;I < length;I++){
                    sprintf(lbl, "B^(%d)_[%d]12", I, irrep);
                    dpd_file2_init(&B, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
                    dpd_file2_axpy(&B, &Bn, Alpha[I][k], 0);
                    dpd_file2_close(&B);
                }
                dpd_file2_close(&Bn);
            }
            for(int k = 0;k < rpi_[irrep];k++){
                sprintf(lbl, "Bn^(%d)_[%d]12", k, irrep);
                dpd_file2_init(&Bn, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
                sprintf(lbl, "B^(%d)_[%d]12", k, irrep);
                dpd_file2_copy(&Bn, PSIF_ADC, lbl);
                dpd_file2_close(&Bn);
            }
            skip_check = 1;
            length = rpi_[irrep];
        }
    
        if(!skip_check){
            zero_int_array(conv, rpi_[irrep]);
            fprintf(iter_adc, "Root          Eigenvalue   Delta     Res_Norm     Conv?\n");
            fprintf(iter_adc, "----     ---------------- -------    --------- ----------\n");
        
            for(int k = 0;k < rpi_[irrep];k++){
                double diff = fabs(lambda[k]-lambda_o[k]);
                if(diff < cutoff && residual_ok[k]){
                    conv[k] = 1;
                    converged++;
                }
                lambda_o[k] = lambda[k];
                fprintf(iter_adc, "%3d  %20.14f %4.3e   %4.3e     %1s\n", k, lambda[k], diff, residual_norm[k], conv[k] == 1 ? "Y" : "N");
                fflush(iter_adc);
            }
        }
        
        int all_conv = 0;
        for(int i = 0;i < num_root;i++) all_conv += conv[i];

        if(all_conv == num_root && converged >= num_root){
            fprintf(iter_adc, "Davidson algorithm converged in %d iterations for %dth root.\n", iter, num_root-1);
            for(int I = 0;I < num_root;I++){
                eps[I] = lambda[I];
                sprintf(lbl, "V^(%d)_[%d]12", I, irrep);
                dpd_file2_init(&V, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
                dpd_file2_scm(&V, 0.0);
                for(int k = 0;k < length;k++){
                    sprintf(lbl, "B^(%d)_[%d]12", k, irrep);
                    dpd_file2_init(&B, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
                    dpd_file2_axpy(&B, &V, Alpha[k][I], 0);
                    dpd_file2_close(&B);
                }
                dpd_file2_close(&V);
            }
            break;
        }
        iter++;
    }
    timer_off("SEM");
    
    fclose(iter_adc);
    free(residual_ok);
    free(residual_norm);
    free(conv);
    free_block(G);
    free_block(Alpha);
    free(lambda);
    free(lambda_o);
    
}

}} // End Namespaces
