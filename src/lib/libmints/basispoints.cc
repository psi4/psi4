#include "molecule.h"
#include "vector3.h"
#include "integral.h"
#include "basisset.h"
#include "basispoints.h"
#include "gridblock.h"
#include "gshell.h"
#include <cstdio>
#include <cmath>
#include <float.h>

#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#define DANGER_TOL 1.0E-25

using namespace psi;
using namespace boost;

namespace psi {

BasisPoints::BasisPoints(boost::shared_ptr<BasisSet> bas, int block_size)
{
    basis_ = bas;
    block_size_ = block_size;
    have_points_ = false;
    have_gradients_ = false;
    have_hessians_ = false;
    have_laplacians_ = false;
    do_points_ = true;
    do_gradients_ = false;
    do_hessians_ = false;
    do_laplacians_ = false;
    
    // Setup the default cutoff radii squared (infinity)
    cutoff_radii_2_ = init_array(basis_->nshell());
    computeCutoffRadii2();

    sig_shells_ = (bool*) malloc(basis_->nshell()*sizeof(bool));
    false_shells_ = (bool*) malloc(basis_->nshell()*sizeof(bool));

    rel2abs_shells_ = init_int_array(basis_->nshell());
    abs2rel_shells_ = init_int_array(basis_->nshell());
    rel2abs_functions_ = init_int_array(basis_->nbf());
    abs2rel_functions_ = init_int_array(basis_->nbf());
    nsig_shells_ = basis_->nshell();
    nsig_functions_ = basis_->nbf();
    block_efficiency_ = 1.0;

    spherical_ = basis_->has_puream();
    max_am_ = basis_->max_am();
    max_carts_ = (max_am_ + 1)*(max_am_ + 2) / 2; 
    max_prim_ = basis_->max_nprimitive();
    nao_ = basis_->nao();
    nso_ = basis_->nbf();
    nshell_ = basis_->nshell();
    nprim_per_shell_ = init_int_array(basis_->nshell());
    ncart_per_shell_ = init_int_array(basis_->nshell());
    nso_per_shell_ = init_int_array(basis_->nshell());
    am_per_shell_ = init_int_array(basis_->nshell());
    a_ = init_int_array(nao_);
    b_ = init_int_array(nao_);
    c_ = init_int_array(nao_);
    xc_ = init_array(nshell_);
    yc_ = init_array(nshell_);
    zc_ = init_array(nshell_);
    coef_ = block_matrix(nao_, max_prim_);
    alpha_ = block_matrix(nshell_, max_prim_);
    prims_ = block_matrix(nao_, max_prim_);
    ex_ = init_array(max_prim_);
    dx_pow_ = init_array(max_am_ + 1);
    dy_pow_ = init_array(max_am_ + 1);
    dz_pow_ = init_array(max_am_ + 1);

    sotrans_count_ = init_int_array(max_am_ + 1);
    sotrans_ao_ = (int**) malloc((max_am_ + 1)*sizeof(int*));
    sotrans_so_ = (int**) malloc((max_am_ + 1)*sizeof(int*));
    sotrans_coef_ = (double**) malloc((max_am_ + 1)*sizeof(double*));

    allocate();
    init_basis_data();
}
BasisPoints::~BasisPoints()
{
    do_points_ = false;
    do_gradients_ = false;
    do_hessians_ = false;
    do_laplacians_ = false;
    release();

    free(cutoff_radii_2_);
    free(sig_shells_);
    free(false_shells_);

    free(rel2abs_shells_);
    free(abs2rel_shells_);
    free(rel2abs_functions_);
    free(abs2rel_functions_);
    
    free(nprim_per_shell_);
    free(ncart_per_shell_);
    free(nso_per_shell_);
    free(am_per_shell_);

    free(a_);
    free(b_);
    free(c_);
    free(xc_);
    free(yc_);
    free(zc_);
    free_block(coef_);
    free_block(alpha_);
    free_block(prims_);
    free(ex_);
    free(dx_pow_);
    free(dy_pow_);
    free(dz_pow_);

    for (int L = 0; L <= max_am_; L++) {
        free(sotrans_so_[L]);
        free(sotrans_ao_[L]);
        free(sotrans_coef_[L]);
    }
    free(sotrans_so_); 
    free(sotrans_ao_); 
    free(sotrans_coef_); 
    free(sotrans_count_); 
}
void BasisPoints::init_basis_data()
{
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(basis_,basis_,basis_,basis_));

    int offset = 0;
    for (int P = 0; P < nshell_; P++) {
        
        // Grab the shell
        boost::shared_ptr<GaussianShell> shell = basis_->shell(P);
        int l = shell->am(); 
        int nao = shell->ncartesian();
        int nprim = shell->nprimitive();
        
        // Shell center
        Vector3 v = shell->center();
        xc_[P] = v[0];
        yc_[P] = v[1];
        zc_[P] = v[2];

        // Cartesians and Primitives in this shell
        nprim_per_shell_[P] = shell->nprimitive();
        ncart_per_shell_[P] = shell->ncartesian();
        nso_per_shell_[P] = shell->nfunction();
        am_per_shell_[P] = shell->am();

        // Cartesian Exponents for each ao within the shell
        int ao = 0;
        for (int i=0; i<=l; ++i) {
            int x = l-i;
            for (int j=0; j<=i; ++j) {
                int y = i-j;
                int z = j;
                a_[ao + offset] = x;
                b_[ao + offset] = y;
                c_[ao + offset] = z;
                ao++;
            }
        }

        // Contraction coefficients, standard normalization, and am normalization
        for (int L = 0; L < nao; L++)
            for (int K = 0; K < nprim; K++) {
                coef_[offset + L][K] =  shell->coef(K) * \
                    shell->normalize(a_[offset + L], b_[offset + L], c_[offset + L]);
            }

        // Alphas for each primitive in the shell
        for (int K = 0; K < nprim; K++) {
            alpha_[P][K] = shell->exp(K);
        } 

        offset += nao;   
    }
    for (int L = 0; L <= max_am_; L++) {
        SphericalTransformIter* trans(fact->spherical_transform_iter(L));
        int G = 0;
        for (trans->first(); !trans->is_done();trans->next()) {
            //double trans_coef = trans.coef();
            //int ind_ao = trans.cartindex();
            //int ind_so = trans.pureindex();
            //fprintf(outfile,"  -L = %d, so = %d, ao = %d, val = %14.10f\n", L, ind_so, ind_ao, trans_coef);
            G++;
        }
        sotrans_count_[L] = G;
        sotrans_so_[L] = init_int_array(G); 
        sotrans_ao_[L] =  init_int_array(G); 
        sotrans_coef_[L] =  init_array(G); 

        G = 0;
        for (trans->first(); !trans->is_done();trans->next()) {
            sotrans_so_[L][G] = trans->pureindex();
            sotrans_ao_[L][G] = trans->cartindex();
            sotrans_coef_[L][G] = trans->coef();
            G++;
        }
        //fprintf(outfile, "  L = %d\n", L);
        //print_mat(so_transforms_[L],2*L+1, (L+1)*(L+2)/2,outfile);
    }
}
void BasisPoints::allocate()
{
    int lmax = basis_->max_am();
    if (do_points_ && !have_points_) {
        if (spherical_) {
            ao_points_ = init_array(nao_);
        }
        points_ = block_matrix(block_size_,basis_->nbf());
        have_points_ = true;
    }
    if (do_gradients_ && !have_gradients_) {
        if (spherical_) {
            ao_gradX_ = init_array(nao_);
            ao_gradY_ = init_array(nao_);
            ao_gradZ_ = init_array(nao_);
        }
        gradientsX_ = block_matrix(block_size_,basis_->nbf());
        gradientsY_ = block_matrix(block_size_,basis_->nbf());
        gradientsZ_ = block_matrix(block_size_,basis_->nbf());
        have_gradients_ = true;
    }
    if (do_hessians_ && !have_hessians_) {
        if (spherical_) {
            ao_hessXY_ = init_array(nao_);
            ao_hessXZ_ = init_array(nao_);
            ao_hessYZ_ = init_array(nao_);
            ao_hessXX_ = init_array(nao_);
            ao_hessYY_ = init_array(nao_);
            ao_hessZZ_ = init_array(nao_);
        }
        hessiansXY_ = block_matrix(block_size_,basis_->nbf());
        hessiansXZ_ = block_matrix(block_size_,basis_->nbf());
        hessiansYZ_ = block_matrix(block_size_,basis_->nbf());
        hessiansXX_ = block_matrix(block_size_,basis_->nbf());
        hessiansYY_ = block_matrix(block_size_,basis_->nbf());
        hessiansZZ_ = block_matrix(block_size_,basis_->nbf());
        have_hessians_ = true;
    }
    if (do_laplacians_ && !have_laplacians_) {
        if (spherical_) {
            ao_laplac_ = init_array((lmax+1)*(lmax+2)>>1);
        }
        laplacians_ = block_matrix(block_size_,basis_->nbf());
        have_laplacians_ = true;
    }

}
void BasisPoints::release()
{
    if (have_points_ && !do_points_) {
        if (spherical_) {
            free(ao_points_);
        }
        free_block(points_);
        have_points_ = false;
    }
    if (have_gradients_ && !do_gradients_) {
        if (spherical_) {
            free(ao_gradX_);
            free(ao_gradY_);
            free(ao_gradZ_);
        }
        free_block(gradientsX_);
        free_block(gradientsY_);
        free_block(gradientsZ_);
        have_gradients_ = false;
    }
    if (have_hessians_ && !do_hessians_) {
        if (spherical_) {
            free(ao_hessXY_);
            free(ao_hessXZ_);
            free(ao_hessYZ_);
            free(ao_hessXX_);
            free(ao_hessYY_);
            free(ao_hessZZ_);
        }
        free_block(hessiansXY_);
        free_block(hessiansXZ_);
        free_block(hessiansYZ_);
        free_block(hessiansXX_);
        free_block(hessiansYY_);
        free_block(hessiansZZ_);
        have_hessians_ = false;
    }
    if (have_laplacians_ && !do_laplacians_) {
        if (spherical_) {
            free(ao_laplac_);
        }
        free_block(laplacians_);
        have_laplacians_ = false;
    }
}
void BasisPoints::computePoints(SharedGridBlock grid)
{
    // NOTE: I have elected not to support mixed spherical/
    // cartesian basis sets at this time. I'm pretty sure libmints 
    // dies on those anyways. Also, C1 for now as usual.

    // Find out what's to be done first
    computeSignificantShells(grid);

    double *xg = grid->getX();
    double *yg = grid->getY();
    double *zg = grid->getZ();
    int ntrue = grid->getTruePoints();
    
    true_size_ = ntrue;    
  
    dx_pow_[0] = 1.0; 
    dy_pow_[0] = 1.0; 
    dz_pow_[0] = 1.0; 

    double x, y, z, dx, dy, dz, ang;
    double reg, regX, regY, regZ, preX, preY, preZ; 
    bool dangerX, dangerY, dangerZ;
    int index, offset, offset_so, P, K, L, G, nK, nL, nG, l; 

    int Prel;

    double trans_coef;
    int ind_ao, ind_so, ind_rel;

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
    //          PUREAM  ALGORITHM         //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//

    if (spherical_) {
    
        // Loop over grid points
        for (index = 0; index < ntrue; index++) {

           // Evaluate primitives
           x = xg[index];
           y = yg[index];   
           z = zg[index];  
           offset = 0; 
           for (Prel = 0; Prel < nsig_shells_; Prel++) {
                P = rel2abs_shells_[Prel];
                nK = nprim_per_shell_[P];                
                nL = ncart_per_shell_[P];                

                // Center displacement
                dx = x - xc_[P];
                dy = y - yc_[P];
                dz = z - zc_[P];
                
                // All possible powers of \vec dR
                for (L = 1; L < max_am_ + 1; L++) {
                    dx_pow_[L] = dx_pow_[L - 1] * dx;
                    dy_pow_[L] = dy_pow_[L - 1] * dy;
                    dz_pow_[L] = dz_pow_[L - 1] * dz;
                }

                // exponential
                for (K = 0; K < nK; K++) {
                    ex_[K] = exp(-alpha_[P][K]*(dx*dx+dy*dy+dz*dz));
                }
                // primitives
                for (L = offset; L < offset + nL; L++) {
                    ang = dx_pow_[a_[L]] * dy_pow_[b_[L]] * dz_pow_[c_[L]]; 
                    for (K = 0; K < nK; K++) {
                        prims_[L][K] = coef_[L][K] * ang * ex_[K];                         
                    } 
                }
                offset += ncart_per_shell_[P];   
            }
            
            // points
            if (do_points_) {

                offset = 0;
                for (Prel = 0; Prel < nsig_shells_; Prel++) {
                    P = rel2abs_shells_[Prel];
                    nK = nprim_per_shell_[P];                
                    nL = ncart_per_shell_[P];                
                    for (L = offset; L < offset + nL; L++) {
                        reg = prims_[L][0];
                        for (K = 1; K < nK; K++) {
                            reg += prims_[L][K];
                        }
                        ao_points_[L] = reg;
                    } 
                    offset += ncart_per_shell_[P];   
                }
            }           
 
            // gradients 
            if (do_gradients_) {

                offset = 0;
                for (Prel = 0; Prel < nsig_shells_; Prel++) {
                    P = rel2abs_shells_[Prel];
                    nK = nprim_per_shell_[P];                
                    nL = ncart_per_shell_[P];
                    
                    dx = x - xc_[P];
                    dy = y - yc_[P];
                    dz = z - zc_[P];

                    dangerX = fabs(dx) < DANGER_TOL;
                    dangerY = fabs(dy) < DANGER_TOL;
                    dangerZ = fabs(dz) < DANGER_TOL;
                    
                    for (L = offset; L < offset + nL; L++) {
                        if (dangerX || a_[L] == 0) 
                            preX = 0.0;
                        else 
                            preX = a_[L]/dx;
                        if (dangerY || b_[L] == 0) 
                            preY = 0.0;
                        else 
                            preY = b_[L]/dy;
                        if (dangerZ || c_[L] == 0) 
                            preZ = 0.0;
                        else 
                            preZ = c_[L]/dz;
                        regX = (preX - 2.0*alpha_[P][0]*dx) * prims_[L][0]; 
                        regY = (preY - 2.0*alpha_[P][0]*dy) * prims_[L][0]; 
                        regZ = (preZ - 2.0*alpha_[P][0]*dz) * prims_[L][0]; 
                        for (K = 1; K < nK; K++) {
                            regX += (preX - 2.0*alpha_[P][K]*dx) * prims_[L][K]; 
                            regY += (preY - 2.0*alpha_[P][K]*dy) * prims_[L][K]; 
                            regZ += (preZ - 2.0*alpha_[P][K]*dz) * prims_[L][K]; 
                        }
                        ao_gradX_[L] = regX;
                        ao_gradY_[L] = regY;
                        ao_gradZ_[L] = regZ;
                    } 
                    offset += ncart_per_shell_[P];   
                }
            }
            // TODO Laplacians 

            //AO -> SO code (sparse/small enough to not merit DGEMM) 
            if (do_points_) {
                memset((void*) &points_[index][0], '\0', nsig_functions_*sizeof(double)); 
            }
            if (do_gradients_) {
                memset((void*) &gradientsX_[index][0], '\0', nsig_functions_*sizeof(double)); 
                memset((void*) &gradientsY_[index][0], '\0', nsig_functions_*sizeof(double)); 
                memset((void*) &gradientsZ_[index][0], '\0', nsig_functions_*sizeof(double)); 
            }

            offset = 0;
            offset_so = 0;
            for (Prel = 0; Prel < nsig_shells_; Prel++) {
                P = rel2abs_shells_[Prel];

                l = am_per_shell_[P]; 
                nG = sotrans_count_[l];

                for (G = 0; G < nG; G++) {
                    ind_so = sotrans_so_[l][G];
                    ind_ao = sotrans_ao_[l][G];
                    trans_coef = sotrans_coef_[l][G];

                    if (do_points_) {
                        points_[index][ind_so + offset_so] += trans_coef*ao_points_[ind_ao + offset];
                    }
                    if (do_gradients_) {
                        gradientsX_[index][ind_so + offset_so] += trans_coef*ao_gradX_[ind_ao + offset];
                        gradientsY_[index][ind_so + offset_so] += trans_coef*ao_gradY_[ind_ao + offset];
                        gradientsZ_[index][ind_so + offset_so] += trans_coef*ao_gradZ_[ind_ao + offset];
                    }
                }
                
                offset += ncart_per_shell_[P];   
                offset_so += nso_per_shell_[P];   
            }
        }
    
    } else {

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
    //       CARTESIAN  ALGORITHM         //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//

        // Loop over grid points
        for (index = 0; index < ntrue; index++) {

           // Evaluate primitives
           x = xg[index];
           y = yg[index];   
           z = zg[index];  
           offset = 0; 
           for (Prel = 0; Prel < nsig_shells_; Prel++) {
                P = rel2abs_shells_[Prel];
                nK = nprim_per_shell_[P];                
                nL = ncart_per_shell_[P];                

                // Center displacement
                dx = x - xc_[P];
                dy = y - yc_[P];
                dz = z - zc_[P];
                
                // All possible powers of \vec dR
                for (L = 1; L < max_am_ + 1; L++) {
                    dx_pow_[L] = dx_pow_[L - 1] * dx;
                    dy_pow_[L] = dy_pow_[L - 1] * dy;
                    dz_pow_[L] = dz_pow_[L - 1] * dz;
                }

                // exponential
                for (K = 0; K < nK; K++) {
                    ex_[K] = exp(-alpha_[P][K]*(dx*dx+dy*dy+dz*dz));
                }
                // primitives
                for (L = offset; L < offset + nL; L++) {
                    ang = dx_pow_[a_[L]] * dy_pow_[b_[L]] * dz_pow_[c_[L]]; 
                    for (K = 0; K < nK; K++) {
                        prims_[L][K] = coef_[L][K] * ang * ex_[K];                         
                    } 
                }
                offset += ncart_per_shell_[P];   
            }
            
            // points
            if (do_points_) {

                offset = 0;
                for (Prel = 0; Prel < nsig_shells_; Prel++) {
                    P = rel2abs_shells_[Prel];
                    nK = nprim_per_shell_[P];                
                    nL = ncart_per_shell_[P];                
                    for (L = offset; L < offset + nL; L++) {
                        reg = prims_[L][0];
                        for (K = 1; K < nK; K++) {
                            reg += prims_[L][K];
                        }
                        points_[index][L] = reg;
                    } 
                    offset += ncart_per_shell_[P];   
                }
            }
            
            // gradients 
            if (do_gradients_) {

                offset = 0;
                for (Prel = 0; Prel < nsig_shells_; Prel++) {
                    P = rel2abs_shells_[Prel];
                    nK = nprim_per_shell_[P];                
                    nL = ncart_per_shell_[P];
                    
                    dx = x - xc_[P];
                    dy = y - yc_[P];
                    dz = z - zc_[P];
    
                    dangerX = fabs(dx) < DANGER_TOL;
                    dangerY = fabs(dy) < DANGER_TOL;
                    dangerZ = fabs(dz) < DANGER_TOL;
                    
                    for (L = offset; L < offset + nL; L++) {
                        if (dangerX || a_[L] == 0) 
                            preX = 0.0;
                        else 
                            preX = a_[L]/dx;
                        if (dangerY || b_[L] == 0) 
                            preY = 0.0;
                        else 
                            preY = b_[L]/dy;
                        if (dangerZ || c_[L] == 0) 
                            preZ = 0.0;
                        else 
                            preZ = c_[L]/dz;
                        regX = (preX - 2.0*alpha_[P][0]*dx) * prims_[L][0]; 
                        regY = (preY - 2.0*alpha_[P][0]*dy) * prims_[L][0]; 
                        regZ = (preZ - 2.0*alpha_[P][0]*dz) * prims_[L][0]; 
                        for (K = 1; K < nK; K++) {
                            regX += (preX - 2.0*alpha_[P][K]*dx) * prims_[L][K]; 
                            regY += (preY - 2.0*alpha_[P][K]*dy) * prims_[L][K]; 
                            regZ += (preZ - 2.0*alpha_[P][K]*dz) * prims_[L][K]; 
                        }
                        gradientsX_[index][L] = regX;
                        gradientsY_[index][L] = regY;
                        gradientsZ_[index][L] = regZ;
                    } 
                    offset += ncart_per_shell_[P];   
                }
    
            }
            // TODO Laplacians 
        }
    }
}
void BasisPoints::computeCutoffRadii2(double epsilon)
{
    cutoff_epsilon_ = epsilon;
    //printf("Epsilon is %14.10E\n", epsilon);

    if (epsilon == 0.0) {
        for (int P = 0; P < basis_->nshell(); P++)
            cutoff_radii_2_[P] = DBL_MAX;
        return;
    }

    // Approximate basis functions as S gaussians
    // This is a weak approximation
    // Valid for epsilon << 1
    for (int P = 0; P < basis_->nshell(); P++) {
        
        boost::shared_ptr<GaussianShell> shell = basis_->shell(P);

        double mean_alpha = 0;
        for (int K = 0; K < shell->nprimitive(); K++)
            mean_alpha += shell->exp(K);
        mean_alpha /= (double) shell->nprimitive();

        // Get some newton-raphson going
        double sigma = sqrt(1.0/(2.0*mean_alpha));

        double r_eps = 3.0*sigma; //Initial guess
        double r_eps_old = 3.0*sigma; //Initial old guess
        double phi = 0.0;
        double del_phi = 0.0;
        int iter = 0;
        int l = shell->am();
        double Nang = shell->normalize(l,0,0);

        do {
            iter = iter + 1;

            phi = 0.0;
            del_phi = 0.0;

            // Compute phi and del_phi
            for (int k = 0; k < shell->nprimitive(); k++) {
                phi += Nang*pow(r_eps,l)*shell->coef(k)*exp(-shell->exp(k)*r_eps*r_eps);
                del_phi += Nang*pow(r_eps,l)*shell->coef(k)*exp(-shell->exp(k)*r_eps*r_eps)*(-2.0*shell->exp(k)*r_eps);
                if (l > 0) {
                    del_phi +=  Nang*l*pow(r_eps,l - 1)*shell->coef(k)*exp(-shell->exp(k)*r_eps*r_eps);
                }
            }
            phi = phi - epsilon; // There's the tolerance

            // Do the newton raphson iteration
            r_eps_old = r_eps;
            r_eps = -phi/del_phi + r_eps;

            //printf("  Iteration %d: r_c = %14.10f, phi_c = %14.10f, r_c+1 = %14.10f\n", \
                iter, r_eps_old, phi, r_eps); 
            if (iter > 2000) {
                // Convergence failure (should never happen)
                printf("WARNING: Basis set cutoff radius computation did not converge.\n");
                r_eps = DBL_MAX;
                break;
            }
        } while (fabs(r_eps-r_eps_old) >= DBL_EPSILON*r_eps);
        //cutoff radii squares are stored on a per-shell basis
        cutoff_radii_2_[P] = r_eps*r_eps;
        // prevent infinity
        if (cutoff_radii_2_[P] > DBL_MAX)
            cutoff_radii_2_[P] = DBL_MAX; 
    }
}
void BasisPoints::computeSignificantShells(SharedGridBlock grid)
{
    double *xg = grid->getX();
    double *yg = grid->getY();
    double *zg = grid->getZ();

    Vector3 shell_center;        
    double R2, x, y, z;

    memset(sig_shells_,'\0',nshell_*sizeof(bool));
    for (int P = 0; P < nshell_; P++)
        false_shells_[P] = true; 

    //for (int Q = 0 ; Q < grid->getTruePoints(); Q ++)
    //    printf(" Q = %6d: (%14.10f, %14.10f, %14.10f)\n", Q, xg[Q], yg[Q],zg[Q]); 
    //for (int Q = 0 ; Q < nshell_; Q ++) 
    //    printf(" Q = %6d: (%14.10f, %14.10f, %14.10f)\n", Q, xc_[Q], yc_[Q],zc_[Q]); 

    // Activate significant shells in this grid block
    // The objective is to keep the trues the same for as many shells as possible
    // and to minimize the overall number of different trues within a block
    //
    //      01100000
    //      01110000
    //      11100000 is OK
    //      
    //      01100000
    //      00011000
    //      00000110 is not OK
    for (int Q = 0 ; Q < grid->getTruePoints(); Q ++) {
        for (int P = 0; P < nshell_; P++) {
            R2 = (xc_[P]-xg[Q])*(xc_[P]-xg[Q]) + (yc_[P]-yg[Q])*(yc_[P]-yg[Q]) + (zc_[P]-zg[Q])*(zc_[P]-zg[Q]);
            //printf(" Q = %d, P = %d, R2 = %14.10f\n", Q, P, R2);
            if (cutoff_radii_2_[P] == DBL_MAX || R2 < cutoff_radii_2_[P])
                sig_shells_[P] = true;
            else
                false_shells_[P] = false;
        }  
    }

    // Setup indexing and compute block efficiency
    int shell_counter = 0;
    int local_counter = 0;
    int function_counter = 0;
    int misses = 0;
    for (int P = 0; P < basis_->nshell(); P++ ) {
        if (sig_shells_[P]) {
            // Efficiency
            if (!false_shells_[P])
                misses++;

            // Indexing
            abs2rel_shells_[P] = shell_counter;
            rel2abs_shells_[shell_counter] = P;
            shell_counter++;
            local_counter = 0;
            for (int k = 0; k < basis_->shell(P)->nfunction(); k++) {
                abs2rel_functions_[basis_->shell(P)->function_index() + local_counter] = function_counter;
                rel2abs_functions_[function_counter] = basis_->shell(P)->function_index() + local_counter;
                local_counter++;
                function_counter++;
            }
        }
    }
    nsig_shells_ = shell_counter;
    nsig_functions_ = function_counter;
    block_efficiency_ = (nsig_shells_-misses)/(double)nsig_shells_; 
    if (nsig_shells_ == 0)
        block_efficiency_ = 1.0;
    
    //printf("n_sig_shells = %d, nsig_functions = %d, block_efficiency = %8.5f\n", nsig_shells_, \
        nsig_functions_, block_efficiency_);
    /**
    printf("Shells:\n");
    for (int P = 0; P < basis_->nshell(); P++ ) {
        printf("  P = %4d, sig = %3s, abs2rel = %4d, rel2abs = %4d\n", P, \
            (sig_shells_[P] ? "Yes" : "No"), abs2rel_shells_[P], rel2abs_shells_[P]);
    }
    printf("Functions:\n");
    for (int P = 0; P < basis_->nbf(); P++ ) {
        printf("  Function = %4d, abs2rel = %4d, rel2abs = %4d\n", P, \
            abs2rel_functions_[P], rel2abs_functions_[P]);
    }*/
}
int BasisPoints::nbf()
{
    return basis_->nbf();
}

}
