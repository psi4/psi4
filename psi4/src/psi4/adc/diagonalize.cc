/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/psi4-dec.h"
#include "psi4/libparallel/ParallelPrinter.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
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
ADCWfn::rhf_diagonalize(int irrep, int num_root, bool first, double omega_in, double *eps)
{
    char lbl[32];
    int iter, converged, prev_length, length, *conv, skip_check, maxdim, *residual_ok;
    double **Alpha, **G, *lambda, *lambda_o, *residual_norm, cutoff;
    dpdfile2 B, S, Bp, Bn, F, L, V;
    dpdfile4 A;

    maxdim = 10 * rpi_[irrep];
    iter = 0;
    converged = 0;
    cutoff = conv_;
    length = rpi_[irrep];
    prev_length = 0;

    residual_ok   = init_int_array(rpi_[irrep]);
    residual_norm = init_array(rpi_[irrep]);
    conv          = init_int_array(rpi_[irrep]);

    G        = block_matrix(maxdim, maxdim);
    Alpha    = block_matrix(maxdim, maxdim);
    lambda   = init_array(maxdim);
    lambda_o = init_array(maxdim);

    for(int I = 0;I < rpi_[irrep];I++) lambda_o[I] = omega_guess_->get(irrep, I);
    shift_denom4(irrep, omega_in);

    std::shared_ptr<OutFile> printer(new OutFile("iter.dat",APPEND));

    timer_on("SEM");
    while(converged < rpi_[irrep] && iter < sem_max_){
        skip_check = 0;
        printer->Printf("\niter = %d, dim = %d\n", iter, length);

        // Evaluating the sigma vectors
        timer_on("Sigma construction");
        for(int I = prev_length;I < length;I++)
            if(!nopen_) rhf_construct_sigma(irrep, I);
        timer_off("Sigma construction");

        // Making so called Davidson mini-Hamiltonian, or Rayleigh matrix
        for(int I = prev_length;I < length;I++){
            sprintf(lbl, "S^(%d)_[%d]12", I, irrep);
            global_dpd_->file2_init(&S, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
            for(int J = 0;J <= I;J++){
                sprintf(lbl, "B^(%d)_[%d]12", J, irrep);
                global_dpd_->file2_init(&B, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
                double sum = global_dpd_->file2_dot(&B, &S);
                if(I != J)
                    G[I][J] = G[J][I] = sum;
                else
                    G[I][J] = sum;
                global_dpd_->file2_close(&B);
            }
            global_dpd_->file2_close(&S);
        }
        if(first && !iter)  poles_[irrep][num_root-1].ps_value = G[num_root-1][num_root-1];
        sq_rsp(length, length, G, lambda, 1, Alpha, 1e-12);

        // Constructing the corretion vectors
        for(int k = 0;k < rpi_[irrep];k++){
            sprintf(lbl, "F^(%d)_[%d]12", k, irrep);
            global_dpd_->file2_init(&F, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
            global_dpd_->file2_scm(&F, 0.0);
            for(int I = 0;I < length;I++){
                sprintf(lbl, "B^(%d)_[%d]12", I, irrep);
                global_dpd_->file2_init(&B, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
                global_dpd_->file2_axpy(&B, &F, -Alpha[I][k]*lambda[k], 0);
                global_dpd_->file2_close(&B);
                sprintf(lbl, "S^(%d)_[%d]12", I, irrep);
                global_dpd_->file2_init(&S, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
                global_dpd_->file2_axpy(&S, &F, Alpha[I][k], 0);
                global_dpd_->file2_close(&S);
            }
            shift_denom2(k, irrep, lambda[k]);
            sprintf(lbl, "L^(%d)_[%d]12", k, irrep);
            global_dpd_->file2_init(&L, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
            global_dpd_->file2_dirprd(&L, &F);
            global_dpd_->file2_close(&L);

            double norm = global_dpd_->file2_dot_self(&F);
            residual_norm[k] = sqrt(norm);
            if(residual_norm[k] > norm_tol_) global_dpd_->file2_scm(&F, 1/residual_norm[k]);
            else {
                global_dpd_->file2_scm(&F, 0.0);
                residual_ok[k] = 1;
            }
            global_dpd_->file2_close(&F);
        }

	prev_length = length;

        // Expand the Ritz space by orthogonalizing {F} to {B} according to Gram-Schmidt procedure
        for(int k = 0;k < rpi_[irrep];k++){
            sprintf(lbl, "F^(%d)_[%d]12", k, irrep);
            global_dpd_->file2_init(&F, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
            sprintf(lbl, "Bpp_[%d]12", irrep);
            global_dpd_->file2_copy(&F, PSIF_ADC_SEM, lbl);
            for(int I = 0;I < length;I++){
                sprintf(lbl, "B^(%d)_[%d]12", I, irrep);
                global_dpd_->file2_init(&B, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
                double coeff = - global_dpd_->file2_dot(&F, &B);
                sprintf(lbl, "Bpp_[%d]12", irrep);
                global_dpd_->file2_init(&Bp, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
                global_dpd_->file2_axpy(&B, &Bp, coeff, 0.0);
                global_dpd_->file2_close(&B);
                global_dpd_->file2_close(&Bp);
            }
            global_dpd_->file2_close(&F);
            sprintf(lbl, "Bpp_[%d]12", irrep);
            global_dpd_->file2_init(&B, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
            double norm = global_dpd_->file2_dot_self(&B);
            norm = sqrt(norm);

            if(norm > norm_tol_){
                global_dpd_->file2_scm(&B, 1/norm);
                sprintf(lbl, "B^(%d)_[%d]12", length, irrep);
                global_dpd_->file2_copy(&B, PSIF_ADC, lbl);
                length++;
            }
            global_dpd_->file2_close(&B);
        }

        if(maxdim-length < rpi_[irrep] || (nxspi_[irrep]-length) < rpi_[irrep]){
            printer->Printf( "Subspace too large:maxdim = %d, L = %d\n", maxdim, length);
            printer->Printf( "Collapsing eigenvectors.\n");

            for(int k = 0;k < rpi_[irrep];k++){
                sprintf(lbl, "Bn^(%d)_[%d]12", k, irrep);
                global_dpd_->file2_init(&Bn, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
                global_dpd_->file2_scm(&Bn, 0.0);
                for(int I = 0;I < prev_length;I++){
                    sprintf(lbl, "B^(%d)_[%d]12", I, irrep);
                    global_dpd_->file2_init(&B, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
                    global_dpd_->file2_axpy(&B, &Bn, Alpha[I][k], 0);
                    global_dpd_->file2_close(&B);
                }
                global_dpd_->file2_close(&Bn);
            }
            for(int k = 0;k < rpi_[irrep];k++){
                sprintf(lbl, "Bn^(%d)_[%d]12", k, irrep);
                global_dpd_->file2_init(&Bn, PSIF_ADC_SEM, irrep, ID('O'), ID('V'), lbl);
                sprintf(lbl, "B^(%d)_[%d]12", k, irrep);
                global_dpd_->file2_copy(&Bn, PSIF_ADC, lbl);
                global_dpd_->file2_close(&Bn);
            }
            skip_check = 1;
            length = rpi_[irrep];
	    prev_length = 0;
        }

        if(!skip_check){
            zero_int_array(conv, rpi_[irrep]);
            printer->Printf("Root          Eigenvalue   Delta     Res_Norm     Conv?\n");
            printer->Printf("----     ---------------- -------    --------- ----------\n");

            for(int k = 0;k < rpi_[irrep];k++){
                double diff = fabs(lambda[k]-lambda_o[k]);
                if(diff < cutoff && residual_ok[k]){
                    conv[k] = 1;
                    converged++;
                }
                lambda_o[k] = lambda[k];
                printer->Printf("%3d  %20.14f %4.3e   %4.3e     %1s\n", k, lambda[k], diff, residual_norm[k], conv[k] == 1 ? "Y" : "N");
            }
        }

        int all_conv = 0;
        for(int i = 0;i < num_root;i++) all_conv += conv[i];

        if(all_conv == num_root && converged >= num_root){
            printer->Printf("Davidson algorithm converged in %d iterations for %dth root.\n", iter, num_root-1);
            for(int I = 0;I < num_root;I++){
                eps[I] = lambda[I];
                sprintf(lbl, "V^(%d)_[%d]12", I, irrep);
                global_dpd_->file2_init(&V, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
                global_dpd_->file2_scm(&V, 0.0);
                for(int k = 0;k < length;k++){
                    sprintf(lbl, "B^(%d)_[%d]12", k, irrep);
                    global_dpd_->file2_init(&B, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
                    global_dpd_->file2_axpy(&B, &V, Alpha[k][I], 0);
                    global_dpd_->file2_close(&B);
                }
                global_dpd_->file2_close(&V);
            }
            break;
        }
        iter++;
    }
    timer_off("SEM");

    free(residual_ok);
    free(residual_norm);
    free(conv);
    free_block(G);
    free_block(Alpha);
    free(lambda);
    free(lambda_o);

}

}} // End Namespaces
