#include "psi4-dec.h"
#include <libmints/mints.h>
#include <lib3index/3index.h>
#include <libqt/qt.h>
#include <liboptions/liboptions.h>
#include <lib3index/3index.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include "physconst.h"
#include "dfadc.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi{ namespace plugin_dfadc {

#define BIG_ 1E50
void 
DFADC::diagonalize(double *&eps, double **&V, int nroot, bool first, bool CISguess)
{       
    eps = init_array((ULI)nroot);
    V   = block_matrix((ULI)nroot, navir_*(ULI)naocc_);
    
    double **Alpha, **G, *lambda_o, *lambda, *Adiag;
    double **B, **S, **F, *Bp, *Bn, *L;
    bool skip_check;
    FILE *iter_adc;
    
    int nthread = 1;
    #ifdef _OPENMP
    nthread = omp_get_max_threads();
    #endif

    // Setting up each parameters
    int iter      = 0;
    int converged = 0;
    int length    = init_ritz_ * nroot;
    int maxdim    = max_ritz_  * nroot; 
    
    bool do_residual  = options_.get_bool("do_residual");
    bool *residual_ok = (bool*)malloc(nroot*sizeof(bool));
    for(int i = 0;i < nroot;i++) residual_ok[i] = !do_residual && CISguess;
    bool *conv        = (bool*)malloc(nroot*sizeof(bool));
    for(int i = 0;i < nroot;i++) conv[i] = false;
    double *residual_norm = init_array(nroot*sizeof(double));
    
    double cutoff    = pow(10, -(conv_+1));
    
    S        = block_matrix((ULI)maxdim, navir_*(ULI)naocc_, false); // Sigma vector
    B        = block_matrix((ULI)maxdim, navir_*(ULI)naocc_, false); // Ritz vector
    F        = block_matrix((ULI)maxdim, navir_*(ULI)naocc_, false); // Correction vector
    G        = block_matrix((ULI)maxdim, (ULI)maxdim, false);        // Rayleigh matrix
    Alpha    = block_matrix((ULI)maxdim, (ULI)maxdim, false);
    lambda   = init_array((ULI)maxdim);
    lambda_o = init_array((ULI)maxdim);
    Adiag    = init_array(navir_*(ULI)naocc_);
    
    // Preparing initial guess vector
    if(CISguess){
        #pragma omp parallel for num_threads(nthread)
        for(int I = 0;I < length;I++)
            C_DCOPY(naocc_*(ULI)navir_, Bcis_[I], 1, B[I], 1);
        C_DCOPY((ULI)nroot, Ecis_, 1, lambda_o, 1);
    }
    else {
        C_DCOPY(navir_*naocc_, diag_, 1, Adiag    , 1);
        for(int I = 0;I < length;I++){
            double minimum = Adiag[0];
            int    pos_min = 0;
            for(int j = 1;j < navir_*naocc_;j++)
                if(minimum > Adiag[j]) {
                    minimum = Adiag[j];
                    pos_min = j;
                }
            B[I][pos_min]  = 1.0;
            lambda_o[I]    = minimum;
            Adiag[pos_min] = BIG_;
        }
        free(Adiag);
    }
        
    ffile(&iter_adc, "iter.dat", 1);
    
    timer_on("SEM");
    while(converged < nroot && iter < sem_max_){
        skip_check = false;
        fprintf(iter_adc, "\nNiter = %d, Dim = %d\n", iter+1, length);
        
        sigma_tensor(B, length, S, CISguess);
        
        // Preparing Davidson mini-Hamiltonian
        double sum;
        for(int I = 0;I < length;I++){
            for(int J = 0;J <= I;J++){
                C_DGEMM('N', 'T', 1, 1, naocc_*navir_, 1.0, B[I], naocc_*navir_, S[J], naocc_*navir_, 0.0, &sum, 1);
                if(I != J) G[I][J] = G[J][I] = sum;
                else       G[I][J] = sum;
            }
        }
        if(first && CISguess) omega_ps_ = G[nroot-1][nroot-1]; // Picking up the CIS(D) like value
        sq_rsp(length, length, G, lambda, 1, Alpha, 1e-12);

        // Forming correction vectors F^k
        for(int k = 0;k < nroot;k++){
            memset((void*)F[k], '\0', navir_*naocc_*sizeof(double));
            for(int I = 0;I < length;I++){
                C_DAXPY(navir_*(ULI)naocc_, -Alpha[I][k]*lambda[k], B[I], 1, F[k], 1);
                C_DAXPY(navir_*(ULI)naocc_,  Alpha[I][k]          , S[I], 1, F[k], 1);
            }
            #pragma omp parallel for schedule(dynamic) num_threads(nthread)
            for(int ia = 0;ia < navir_*naocc_;ia++) {
                double denom = lambda[k] - diag_[ia];
                if(fabs(denom) > pow(10, -norm_tol_))F[k][ia] /= denom;
                else F[k][ia] = 0.0;
            }
            //residual_norm[k] = C_DNRM2(navir_*(ULI)naocc_, F[k], 1);
            C_DGEMM('N', 'T', 1, 1, naocc_*navir_, 1.0, F[k], naocc_*navir_, F[k], naocc_*navir_, 0.0, &residual_norm[k], 1);
            residual_norm[k] = sqrt(residual_norm[k]);

            if(residual_norm[k] > pow(10, -norm_tol_)){
                double reciprocal = 1 / residual_norm[k];
                C_DSCAL(navir_*(ULI)naocc_, reciprocal, F[k], 1);
            }
            else {
                memset((void*)F[k], '\0', navir_*naocc_*sizeof(double));
                residual_ok[k] = true;
            }
        }
        
        // Adding {F} to {B}, to expand the Ritz space
        double coeff, norm;
        for(int k = 0;k < nroot;k++){
            double *Bp = init_array(navir_*(ULI)naocc_);
            C_DCOPY(navir_*(ULI)naocc_, F[k], 1, Bp, 1);
            for(int I = 0;I < length;I++){
                //double coeff = C_DDOT(navir_*(ULI)naocc_, F[k], 1, B[I], 1);
                C_DGEMM('N', 'T', 1, 1, naocc_*navir_, 1.0, F[k], naocc_*navir_, B[I], naocc_*navir_, 0.0, &coeff, 1);
                C_DAXPY(navir_*(ULI)naocc_, -coeff, B[I], 1, Bp, 1);
            }
            //double norm = C_DNRM2(navir_*(ULI)naocc_, Bp, 1);
            C_DGEMM('N', 'T', 1, 1, naocc_*navir_, 1.0, Bp, naocc_*navir_, Bp, naocc_*navir_, 0.0, &norm, 1);
            norm = sqrt(norm);
            if(norm > pow(10, -norm_tol_)){
                double reciprocal = 1 / norm;
                C_DSCAL(navir_*(ULI)naocc_, reciprocal, Bp, 1);
                C_DCOPY(navir_*(ULI)naocc_, Bp, 1, B[length], 1);
                length++;
            }
            free(Bp);
        }
        if(maxdim-length < nroot || naocc_*navir_-length < nroot){
            fprintf(iter_adc, "Subspace too large: maxdim = %3d, length = %3d\n", maxdim, length);
            fprintf(iter_adc, "Collapsing eigenvectors.\n");
            
            double **Bnew = block_matrix((ULI)nroot, navir_*(ULI)naocc_);
            for(int k = 0;k < nroot;k++){
                for(int I = 0;I < length;I++){
                    C_DAXPY(navir_*(ULI)naocc_, Alpha[I][k], B[I], 1, Bnew[k], 1);   
                }
            }
            for(int k = 0;k < nroot;k++) C_DCOPY(navir_*(ULI)naocc_, Bnew[k], 1, B[k], 1);
            skip_check = true;
            length = nroot;
            free_block(Bnew);
        }
        
        if(!skip_check){
            converged = 0;
            for(int i = 0;i < nroot;i++) conv[i] = false;
            fprintf(iter_adc, "Root          Eigenvalue   Delta     Res_Norm     Conv?\n");
            fprintf(iter_adc, "----     ---------------- -------    --------- ----------\n");
            
            for(int k = 0;k < nroot;k++){
                double diff = fabs(lambda[k]-lambda_o[k]);
                if(diff < cutoff && residual_ok[k]){
                    conv[k] = true;
                    converged++;
                }
                lambda_o[k] = lambda[k];
                fprintf(iter_adc, "%3d  %20.14f %4.3e   %4.3e     %1s\n", k, lambda[k], diff, residual_norm[k], conv[k] == true ? "Y" : "N");
                fflush(iter_adc);
            }
        }

        if(converged == nroot){
            memset((void*)V[0], '\0', naocc_*navir_*nroot*sizeof(double));
            fprintf(iter_adc, "Davidson algorithm converged in %d iterations for %dth root.\n", iter+1, nroot);
            for(int I = 0;I < nroot;I++){
                eps[I] = lambda[I];
                for(int k = 0;k < length;k++){
                    C_DAXPY(navir_*(ULI)naocc_, Alpha[k][I], B[k], 1, V[I], 1);
                }
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

void
DFADC::sigma_tensor(double **B, int dim, double **&S, bool do_ADC)
{
    int nthread = 1;
    #ifdef _OPENMP
    nthread = omp_get_max_threads();
    #endif

    double temp1;
    double **temp2 = block_matrix((ULI)naocc_, (ULI)navir_);
    double **temp3 = block_matrix((ULI)navir_, (ULI)navir_);
        
    // \sigma^R_{ia} <-- \sum_{ia} (e_a-e_i) b^R_{ia}
    #pragma omp collapse(3) parallel for schedule(dynamic) num_threads(nthread)
    for(int R = 0;R < dim;R++)
        for(int i = 0;i < naocc_;i++)
            for(int a = 0;a < navir_;a++)
                S[R][i*navir_+a] = B[R][i*navir_+a] * (-occe_[i]+vire_[a]);
        
    for(int R = 0;R < dim;R++){
        
        // \sigma^R_{ia} <-- \sum_{jb} H^{CIS}_{ia,jb} b^R_{jb}
        for(int I = 0;I < ribasis_->nbf();I++){
            
            // J-term, \sigma^R_{ia} <-- 2 Q^I_{ia} \sum_{jb} Q^I_{jb} b^R_{jb}
            C_DGEMM('T', 'N', 1, 1, naocc_*navir_, 1.0, Qpia_[I], 1, B[R], 1, 0.0, &(temp1), 1);
            C_DAXPY(navir_*(ULI)naocc_, 2.0*temp1, Qpia_[I], 1, S[R], 1);
        
            // K-term, \sigma^R_{ia} <-- - \sum_{jb} Q^I_{ij} b^R_{jb} Q^I_{ab}
            C_DGEMM('N', 'T', naocc_, navir_, navir_, -1.0, B[R],    navir_, Qpab_[I],  navir_, 0.0, temp2[0], navir_);
            C_DGEMM('N', 'N', naocc_, navir_, naocc_,  1.0, Qpij_[I], naocc_, temp2[0], navir_, 1.0, S[R],     navir_);
                        
        } // End loop over I
        
        // Compute the contributions to the sigma vectors from all the seocnd order diagrams
        if(do_ADC){
            // \sigma^R_{ia} <-- \sum_{b} b^R_{ib} A_{ab}
            C_DGEMM('N', 'T', naocc_, navir_, navir_, -1.0, B[R], navir_, Aab_[0], navir_, 1.0, S[R], navir_);
            // \sigma^R_{ia} <-- \sum_{j} A_{ij} b^R_{ja}
            C_DGEMM('N', 'N', naocc_, navir_, naocc_, -1.0, Aij_[0], naocc_, B[R], navir_, 1.0, S[R], navir_);
            
            double **Dia = block_matrix((ULI)naocc_, (ULI)navir_);
            for(int I = 0;I < ribasis_->nbf();I++){
                
                // J-term, D_{ia} <-- 2.0 Q^I_{ia} \sum_{jb} Q^I_{jb} b^S_{jb} 
                C_DGEMM('T', 'N', 1, 1, naocc_*navir_, 1.0, Qpia_[I], 1, B[R], 1, 0.0, &(temp1), 1);
                C_DAXPY(navir_*(ULI)naocc_, temp1, Qpia_[I], 1, Dia[0], 1);
                
                // K-term, D_{ia} <-- - \sum_{jb} Q^I_{ib} b^R_{jb} Q^I_{ja}
                C_DGEMM('T', 'N', navir_, navir_, naocc_,  1.0, B[R], navir_, Qpia_[I], navir_, 0.0, temp3[0], navir_);
                C_DGEMM('N', 'N', naocc_, navir_, navir_, -0.5, Qpia_[I], navir_, temp3[0], navir_, 1.0, Dia[0], navir_);
                
            } // End loop over I
            
            // \sigma^R_{ia} <-- 0.5 \sum_{jb} (2 K_{iajb} - K_{ibja}) D_{jb}
            C_DGEMM('N', 'N', naocc_*navir_, 1, naocc_*navir_, 1.0, Kiajb_[0], naocc_*navir_, Dia[0], 1, 1.0, S[R],    1);
            free_block(Dia);

            double **Eia = block_matrix((ULI)naocc_, (ULI)navir_);
            C_DGEMM('N', 'N', naocc_*navir_, 1, naocc_*navir_, 1.0, Kiajb_[0], naocc_*navir_, B[R],   1,  0.0, Eia[0], 1);
            
            for(int I = 0;I < ribasis_->nbf();I++){
                
                // J-term, \sigma^S_{ia} <-- 2.0 Q^I_{ia} \sum_{jb} Q^I_{jb} E_{jb}
                C_DGEMM('T', 'N', 1, 1, naocc_*navir_, 1.0, Qpia_[I], 1, Eia[0], 1, 0.0, &(temp1), 1);
                C_DAXPY(navir_*(ULI)naocc_, temp1, Qpia_[I], 1, S[R], 1);
                
                // K-term, \sigma^S_{ia} <-- \sum_{jb} Q^I_{ib} E_{jb} Q^I_{ja}
                C_DGEMM('T', 'N', navir_, navir_, naocc_,  1.0, Eia[0], navir_, Qpia_[I], navir_, 0.0, temp3[0], navir_);
                C_DGEMM('N', 'N', naocc_, navir_, navir_, -0.5, Qpia_[I], navir_, temp3[0], navir_, 1.0, S[R], navir_);
                
            } // End loop over I
            
            free_block(Eia);
            // Then, all the 3h-3p contributions referred to as the differential correlation terms are accounted for in the sigma vector!
            
            // Next, calculate all the 2h-2p diagrams
            double **temp4 = block_matrix((ULI)ribasis_->nbf(), naocc_*(ULI)navir_);
            
            // Z_{iajb} <-- \sum_{I} B^I_{ja} \sum_{c} b^S_{ic} B^I_{cb} 
            double **Ziajb = block_matrix(naocc_*(ULI)navir_,   naocc_*(ULI)navir_);
            for(int I = 0;I < ribasis_->nbf();I++)
                C_DGEMM('N', 'N', naocc_, navir_, navir_, 1.0, B[R], navir_, Qpab_[I], navir_, 1.0, temp4[I], navir_);
            
            C_DGEMM('T', 'N', naocc_*navir_, naocc_*navir_, ribasis_->nbf(), 1.0, Qpia_[0], naocc_*navir_, temp4[0], naocc_*navir_, 0.0, Ziajb[0], naocc_*navir_);
            
            free_block(temp4);

            temp4 = block_matrix((ULI)ribasis_->nbf(), naocc_*(ULI)navir_);
            
            // Z_{iajb} <-- - \sum_{I} B^I_{ia} \sum_{k} B^I_{jk} b^S_{kb}
            for(int I = 0;I < ribasis_->nbf();I++)
                C_DGEMM('N', 'N', naocc_, navir_, naocc_, 1.0, Qpij_[I], naocc_, B[R], navir_, 1.0, temp4[I], navir_);
            
            C_DGEMM('T', 'N', naocc_*navir_, naocc_*navir_, ribasis_->nbf(), -1.0, Qpia_[0], naocc_*navir_, temp4[0], naocc_*navir_, 1.0, Ziajb[0], naocc_*navir_);
            
            free_block(temp4);
            
            // W_{iajb} <-- (2 Z_{iajb} - Z_{ibja} - Z_{jaib} + 2 Z_{jbia})
            //            / (\omega + e_i + e_j - e_a - e_b)
            // One of the most nastiest steps!!
            double **Wiajb = block_matrix(naocc_*(ULI)navir_, naocc_*(ULI)navir_);
            #pragma omp collapse(4) parallel for schedule(dynamic) num_threads(nthread)
            for(int i = 0;i < naocc_;i++)
                for(int a = 0;a < navir_;a++)
                    for(int j = 0;j < naocc_;j++)
                        for(int b = 0;b < navir_;b++)
                            Wiajb[i*navir_+a][j*navir_+b] = 
                            ( 2*Ziajb[i*navir_+a][j*navir_+b]  -Ziajb[i*navir_+b][j*navir_+a]
                               -Ziajb[j*navir_+a][i*navir_+b]+2*Ziajb[j*navir_+b][i*navir_+a])
                            / (omega_+occe_[i]+occe_[j]-vire_[a]-vire_[b]);
            free_block(Ziajb);
            
            // \sigma^R_{ia} <--   \sum_{Ib} B^I_{ab} \sum_{jc} W_{ibjc} B^I_{jc}
            // \sigma^R_{ia} <-- - \sum_{Ij} B^I_{ji} \sum_{kb} W_{jakb} B^I_{kb}
            for(int I = 0;I < ribasis_->nbf();I++) {
                
                C_DGEMM('N', 'N', naocc_*navir_, 1, naocc_*navir_, 1.0, Wiajb[0], naocc_*navir_, Qpia_[I], 1, 0.0, temp2[0], 1);
                C_DGEMM('N', 'T', naocc_, navir_, navir_,  1.0, temp2[0], navir_, Qpab_[I], navir_, 1.0, S[R], navir_);
                C_DGEMM('T', 'N', naocc_, navir_, naocc_, -1.0, Qpij_[I], naocc_, temp2[0], navir_, 1.0, S[R], navir_);
                
            } // End loop over I

            free_block(Wiajb);
                        
        } // End do_ADC
    } // End loop over R
    
    free_block(temp2);
    free_block(temp3);
}

double
DFADC::accelerate(double *V)
{
    int nthread = 1;
    #ifdef _OPENMP
    nthread = omp_get_max_threads();
    #endif
        
    double partial_omega;
    double **temp1 = block_matrix((ULI)ribasis_->nbf(), naocc_*(ULI)navir_);
    double **temp2 = block_matrix((ULI)naocc_, (ULI)navir_);
                    
    // Calculate V(\omega)^{\dagger}\partial A(\omega)/\partial\omega V^{\dagger}
    // where, A stands for the correlated response matrix
    
    // dZ_{iajb} <-- \sum_{I} B^I_{ja} \sum_{c} V_{ic} B^I_{cb} 
    double **Ziajb = block_matrix(naocc_*(ULI)navir_,   naocc_*(ULI)navir_);
    for(int I = 0;I < ribasis_->nbf();I++)
        C_DGEMM('N', 'N', naocc_, navir_, navir_, 1.0, V, navir_, Qpab_[I], navir_, 1.0, temp1[I], navir_);
                
    C_DGEMM('T', 'N', naocc_*navir_, naocc_*navir_, ribasis_->nbf(), 1.0, Qpia_[0], naocc_*navir_, temp1[0], naocc_*navir_, 0.0, Ziajb[0], naocc_*navir_);
                
    free_block(temp1);
                
    temp1 = block_matrix((ULI)ribasis_->nbf(), naocc_*(ULI)navir_);
                
    // dZ_{iajb} <-- - \sum_{I} B^I_{ia} \sum_{k} B^I_{jk} V_{kb}
    for(int I = 0;I < ribasis_->nbf();I++)
        C_DGEMM('N', 'N', naocc_, navir_, naocc_, 1.0, Qpij_[I], naocc_, V, navir_, 1.0, temp1[I], navir_);
                
    C_DGEMM('T', 'N', naocc_*navir_, naocc_*navir_, ribasis_->nbf(), -1.0, Qpia_[0], naocc_*navir_, temp1[0], naocc_*navir_, 1.0, Ziajb[0], naocc_*navir_);
                
    free_block(temp1);
                
    // dW_{iajb} <-- (2 dZ_{iajb} - dZ_{ibja} - dZ_{jaib} + 2 dZ_{jbia})
    //            / (\omega + e_i + e_j - e_a - e_b)^2
    // One of the most nastiest steps!!
    double **Wiajb = block_matrix(naocc_*(ULI)navir_, naocc_*(ULI)navir_);
    #pragma omp collapse(4) parallel for schedule(dynamic) num_threads(nthread)
    for(int i = 0;i < naocc_;i++)
        for(int a = 0;a < navir_;a++)
            for(int j = 0;j < naocc_;j++)
                for(int b = 0;b < navir_;b++)
                    Wiajb[i*navir_+a][j*navir_+b] = 
                        ( 2*Ziajb[i*navir_+a][j*navir_+b]  -Ziajb[i*navir_+b][j*navir_+a]
                           -Ziajb[j*navir_+a][i*navir_+b]+2*Ziajb[j*navir_+b][i*navir_+a])
                                / (omega_+occe_[i]+occe_[j]-vire_[a]-vire_[b])
                                / (omega_+occe_[i]+occe_[j]-vire_[a]-vire_[b]);
    
    free_block(Ziajb);
                
    double **dV = block_matrix((ULI)naocc_, (ULI)navir_);
    // dV_{ia} <--   \sum_{Ib} B^I_{ab} \sum_{jc} W_{ibjc} B^I_{jc}
    // dV_{ia} <-- - \sum_{Ij} B^I_{ji} \sum_{kb} W_{jakb} B^I_{kb}
    for(int I = 0;I < ribasis_->nbf();I++) {
                    
        C_DGEMM('N', 'N', naocc_*navir_, 1, naocc_*navir_, 1.0, Wiajb[0], naocc_*navir_, Qpia_[I], 1, 0.0, temp2[0], 1);
        C_DGEMM('N', 'T', naocc_, navir_, navir_, -1.0, temp2[0], navir_, Qpab_[I], navir_, 1.0, dV[0], navir_);
        C_DGEMM('T', 'N', naocc_, navir_, naocc_,  1.0, Qpij_[I], naocc_, temp2[0], navir_, 1.0, dV[0], navir_);
                    
    } // End loop over I

    // V(\omega)^{\dagger}\partial A(\omega)/\partial\omega V^{\dagger} = V^{\dagger} dV^{\dagger}
    C_DGEMM('T', 'N', 1, 1, naocc_*navir_, 1.0, V, 1, dV[0], 1, 0.0, &(partial_omega), 1);
    
    free_block(Wiajb);
    free_block(temp2);
    
    return partial_omega;
}

}} // End Namespaces
