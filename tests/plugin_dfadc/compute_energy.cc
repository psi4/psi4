#include "psi4-dec.h"
#include <libmints/mints.h>
#include <lib3index/3index.h>
#include <libciomr/libciomr.h>
#include "physconst.h"
#include "dfadc.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi{ namespace plugin_dfadc {
    
double 
DFADC::compute_energy()
{       
    double **J_mhalf, **Qpia, **Qpij, **Qpab, **Qpai, **Qiap;
    int nQ = ribasis_->nbf();
    
    int nthread = 1;
    #ifdef _OPENMP
    nthread = omp_get_max_threads();
    #endif
    
    fprintf(outfile, "\n\t==> DF-CIS/ADC(1) Level <==\n\n");
    
    diag_ = init_array(naocc_*navir_);
    #pragma omp parallel for num_threads(nthread)
    for(int i = 0;i < naocc_;i++)
        for(int a = 0;a < navir_;a++)
            diag_[i*navir_+a] = vire_[a] - occe_[i];

    diagonalize(Ecis_, Bcis_, num_roots_, false, false);
    
    int HOMO = naocc_;
    int LUMO = naocc_+1;
    for(int I = 0;I < num_roots_;I++) {
        fprintf(outfile, "->\tState %3d, Omega : %14.10f (a.u.), %11.7f (eV)\n", I+1, Ecis_[I], Ecis_[I]*_hartree2ev);
        fprintf(outfile, "\t--------------------------------------------------------------\n");
        for(int dim = 0;dim < naocc_*navir_;dim++)
            if(fabs(Bcis_[I][dim]) > pow(10, -cutoff_)) {
                int occ = dim / navir_+1;
                int vir = dim % navir_+1+naocc_;
                fprintf(outfile, "\t %5d -> %5d %14.10f", occ, vir, Bcis_[I][dim]);
                if(occ == HOMO) fprintf(outfile, " from HOMO");
                if(vir == LUMO) fprintf(outfile, " 2 LUMO");
                fprintf(outfile, "\n");
            }
        fprintf(outfile, "\n");
        // Compute oscillator strength
        double f = 0.0;
        double *mu = init_array(3);
        for(int X = 0;X < 3;X++) {
            mu[X] = C_DDOT(naocc_*(ULI)navir_,Bcis_[I], 1, dipole_ints_[X]->pointer()[0], 1);
            f += 4.0 / 3.0 * mu[X] * mu[X] * Ecis_[I];
        }
        fprintf(outfile, "\tOscillator strength : %10.7f\n", f);
        fprintf(outfile, "\t[X: %10.7f, Y: %10.7f, Z: %10.7f]\n\n", mu[0]*mu[0], mu[1]*mu[1], mu[2]*mu[2]);
        free(mu);
        fflush(outfile);
    }
    
    fprintf(outfile, "\t==> GS energy <==\n");
    
    double E_MP2J = 0;
    double E_MP2K = 0;
    double iajb, ibja;
    
    formInvSqrtJ(J_mhalf);            
    formDFtensor(occCa_, virCa_, J_mhalf, Icol, Qiap);
    double *V = (double *)malloc(navir_*navir_*sizeof(double));
    for(int i = 0;i < naocc_;i++){
        for(int j = 0;j < naocc_;j++){
            C_DGEMM('N', 'T', navir_, navir_, nQ, 1.0, Qiap[i*navir_], nQ, Qiap[j*navir_], nQ, 0.0, V, navir_);
            #pragma omp parallel for reduction(+:E_MP2J,E_MP2K)
            for(int a = 0;a < navir_;a++){
                for(int b = 0;b < navir_;b++){
                    iajb = V[a*navir_+b];
                    ibja = V[b*navir_+a];
                    E_MP2J -= 2.0 * iajb * iajb / (vire_[a]+vire_[b]-occe_[i]-occe_[j]);
                    E_MP2K +=       iajb * ibja / (vire_[a]+vire_[b]-occe_[i]-occe_[j]);
                }
            }
        }
    }
    free(V);
    free_block(Qiap);
    
    double E_DFMP2 = E_MP2J + E_MP2K;
    fprintf(outfile, "\n\tDF-MP2  J energy   = %14.10f\n", E_MP2J);
    fprintf(outfile, "\tDF-MP2  K energy   = %14.10f\n", E_MP2K);
    fprintf(outfile, "->\tDF-MP2 energy      = %14.10f\n", E_DFMP2);
    fflush(outfile);

    // Form X_{ij} and the symmetrized A_{ij} tensors
    formDFtensor(occCa_, virCa_, J_mhalf, Irow, Qpia);
    double **Viajb = block_matrix(naocc_*(ULI)navir_, naocc_*(ULI)navir_);
    
    C_DGEMM('T', 'N', naocc_*navir_, naocc_*navir_, nQ, 1.0, Qpia[0], naocc_*navir_, Qpia[0], naocc_*navir_, 0.0, Viajb[0], naocc_*navir_);
    
    free_block(Qpia);
        
    Kiajb_ = block_matrix(naocc_*(ULI)navir_, naocc_*(ULI)navir_);
    C_DCOPY(naocc_*navir_*naocc_*(ULI)navir_, Viajb[0], 1, Kiajb_[0], 1);
    C_DSCAL(naocc_*navir_*naocc_*(ULI)navir_,  2.0, Kiajb_[0], 1);
    #pragma omp parallel for num_threads(nthread)
    for(int i = 0;i < naocc_;i++)
        for(int a = 0;a < navir_;a++)
            for(int j = 0;j < naocc_;j++)
                for(int b = 0;b < navir_;b++){
                    Kiajb_[i*navir_+a][j*navir_+b] -= Viajb[i*navir_+b][j*navir_+a];
                    Kiajb_[i*navir_+a][j*navir_+b] /= - vire_[a] - vire_[b] + occe_[i] + occe_[j];
                }
                            
    // X_{ij} <-- \sum_{kab} D_{iajb} (2V_{iakb}-V_{ikab}) V_{jakb}
    double **Xij = block_matrix((ULI)naocc_, (ULI)naocc_);
    C_DGEMM('N', 'T', naocc_, naocc_, navir_*naocc_*navir_, 1.0, Kiajb_[0], navir_*naocc_*navir_, Viajb[0], navir_*naocc_*navir_, 0.0, Xij[0], naocc_);
    
    free_block(Viajb);
 
    // Symetrize Xij to form Aij
    Aij_ = block_matrix((ULI)naocc_, (ULI)naocc_);
    #pragma omp parallel for num_threads(nthread)
    for(int i = 0;i < naocc_;i++)
        for(int j = 0;j < naocc_;j++)
            Aij_[i][j] = 0.5 * (Xij[i][j] + Xij[j][i]);
    
    free_block(Xij);
    
    // Form X_{ab} and the symmetrized A_{ab} tensors
    double **Vaibj = block_matrix(naocc_*(ULI)navir_, naocc_*(ULI)navir_);
    formDFtensor(virCa_, occCa_, J_mhalf, Irow, Qpai);
    
    C_DGEMM('T', 'N', naocc_*navir_, naocc_*navir_, nQ, 1.0, Qpai[0], naocc_*navir_, Qpai[0], naocc_*navir_, 0.0, Vaibj[0], naocc_*navir_);
    
    free_block(Qpai);
        
    double **Kaibj = block_matrix(naocc_*(ULI)navir_, naocc_*(ULI)navir_);
    C_DCOPY(naocc_*navir_*naocc_*(ULI)navir_, Vaibj[0], 1, Kaibj[0], 1);
    C_DSCAL(naocc_*navir_*naocc_*(ULI)navir_,  2.0, Kaibj[0], 1);
    #pragma omp parallel for num_threads(nthread)
    for(int a = 0;a < navir_;a++)
        for(int i = 0;i < naocc_;i++)
            for(int b = 0;b < navir_;b++)
                for(int j = 0;j < naocc_;j++) {
                    Kaibj[a*naocc_+i][b*naocc_+j] -= Vaibj[b*naocc_+i][a*naocc_+j];
                    Kaibj[a*naocc_+i][b*naocc_+j] /= - vire_[a] - vire_[b] + occe_[i] + occe_[j];
                }
    
    // X_{ab} <-- \sum_{ijc} D_{aicj} (2V_{aicj}-V_{acij}) V_{bicj}
    double **Xab = block_matrix((ULI)navir_, (ULI)navir_);
    C_DGEMM('N', 'T', navir_, navir_, naocc_*navir_*naocc_, 1.0, Kaibj[0], naocc_*navir_*naocc_, Vaibj[0], naocc_*navir_*naocc_, 0.0, Xab[0], navir_);
    
    free_block(Kaibj);
    free_block(Vaibj);
    
    // Symmetrize X_{ab} to form A_{ab} tensor
    Aab_ = block_matrix((ULI)navir_, (ULI)navir_);
    #pragma omp parallel for num_threads(nthread)
    for(int a = 0;a < navir_;a++)
        for(int b = 0;b < navir_;b++)
            Aab_[a][b] = 0.5 * (Xab[a][b]+Xab[b][a]);

    free_block(Xab);
    
    int iter;
    bool is_first;
    double denom, *Ohms, **Xs;
    
    fprintf(outfile, "\n\t==> DF-ADC(2) Computation <==\n\n");
    omega_ = Ecis_[0];
    for(int nroot = 0;nroot < num_roots_;nroot++){
        omega_ = Ecis_[nroot];
        is_first = true;
        for(iter = 1;iter <= pole_max_;iter++){
            diagonalize(Ohms, Xs, nroot+1, is_first, true);
            denom = 1 - accelerate(Xs[nroot]); //printf("%14.10f\n", 1/denom);
            is_first = false;
            double diff = (omega_ - Ohms[nroot]) / denom;
            if(fabs(diff) > pow(10, -conv_)) omega_ -= diff;
            else break;
        } // End loop over iter
        if(iter == pole_max_) 
            fprintf(outfile, "\t Pole %3d did not converge!", nroot);
        else{
            fprintf(outfile, "->\tState %3d, Omega   : %14.10f (a.u.), %11.7f (eV)\n", nroot+1, omega_, omega_*_hartree2ev);
            fprintf(outfile, "\tNon-iterative value: %14.10f (a.u.), %11.7f (eV)\n", omega_ps_, omega_ps_*_hartree2ev);
            fprintf(outfile, "\t--------------------------------------------------------------\n");
            for(int dim = 0;dim < naocc_*navir_;dim++)
                if(fabs(Xs[nroot][dim]) > pow(10, -cutoff_)) {
                    int occ = dim / navir_+1;
                    int vir = dim % navir_+1+naocc_;
                    fprintf(outfile, "\t %5d -> %5d %14.10f", occ, vir, Xs[nroot][dim]);
                    if(occ == HOMO) fprintf(outfile, " from HOMO");
                    if(vir == LUMO) fprintf(outfile, " 2 LUMO");
                    fprintf(outfile, "\n");
                }
            fprintf(outfile, "\n");
            fprintf(outfile, "\tConverged in %3d iterations.\n", iter);
            fprintf(outfile, "\tSquared norm of the singles component   : %10.7f\n", 1/denom);
            // Calculate rotation angle from the CIS vector
            double theta = acos(C_DDOT(naocc_*(ULI)navir_, Xs[nroot], 1, Bcis_[nroot], 1)) * 180.0 / _pi;
            if((180.0-fabs(theta)) < theta) theta = 180.0 - fabs(theta);
            fprintf(outfile, "\tThe rotation angle of the singles vector: %11.7f\n", theta);
            // Compute *ZEROTH ORDER* oscillator strength
            double f = 0.0;
            double *mu = init_array(3);
            for(int X = 0;X < 3;X++) {
                mu[X] = C_DDOT(naocc_*(ULI)navir_,Xs[nroot], 1, dipole_ints_[X]->pointer()[0], 1);
                f += 4.0 / 3.0 * mu[X] * mu[X] * omega_ * denom;
            }
            fprintf(outfile, "\tZeroth order oscillator strength        : %10.7f\n", f);
            fprintf(outfile, "\t\t[X: %10.7f, Y: %10.7f, Z: %10.7f]\n\n", mu[0]*mu[0], mu[1]*mu[1], mu[2]*mu[2]);
            free(mu);
            fflush(outfile);
        } // End loop over dim
    } // End loop over nroot
    
    release_mem();
    
    return E_DFMP2;
}
        
}} // End Namespaces
