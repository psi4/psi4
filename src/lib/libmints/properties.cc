#include "mints.h"
#include <cstdio>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

using namespace psi;
using namespace boost;

namespace psi {

Properties::Properties(shared_ptr<BasisSet> _b, int _block_size): BasisPoints(_b, _block_size)
{
    // Default to do density only
    do_density_ = false;
    do_density_gradient_ = false;
    do_ke_density_ = false;
    setToComputeDensity(true);

    // Allocate D and C temporary arrays
    Da_ = block_matrix(_b->nbf(), _b->nbf());
    Ca_ = block_matrix(_b->nbf(), _b->nbf());
    Db_ = block_matrix(_b->nbf(), _b->nbf());
    Cb_ = block_matrix(_b->nbf(), _b->nbf());

    // Allocate temporary tensor
    temp_tens_ = block_matrix(_block_size, _b->nbf());
}
Properties::~Properties()
{
    // Free registers
    setToComputeDensity(false);
    setToComputeDensityGradient(false);
    setToComputeKEDensity(false);

    // Free D and C temporary arrays
    free_block(Da_);
    free_block(Ca_);
    free_block(Db_);
    free_block(Cb_);

    // Free temporary tensor
    free_block(temp_tens_);
}
void Properties::computeRKSProperties(shared_ptr<GridBlock> grid, shared_ptr<Matrix> D, shared_ptr<Matrix> C, int* docc)
{
    // All C1 for now (irrep information flow is a problem)

    // Compute significant basis points for this block
    timer_on("Points");
    computePoints(grid);
    timer_off("Points");

    // Get sizes
    int npoints = grid->getTruePoints();
    int nbf = basis_->nbf();
    int nsigf = nsig_functions_;

    // Build the reduced Da_ (always needed)
    for (int m = 0; m < nsigf; m++)
        for (int n = 0; n <= m; n++) {
            Da_[m][n] = D->get(0,rel2abs_functions_[m], rel2abs_functions_[n]);
            Da_[n][m] = Da_[m][n];
        }

    if (do_density_) {
        // rho_a_
        // rho_a^Q = phi_m^Q * Da_mn * phi_n^Q
        C_DGEMM('N', 'N', npoints, nsigf, nsigf, 1.0, &points_[0][0], nbf, &Da_[0][0], nbf, \
            0.0, &temp_tens_[0][0], nbf);

        for (int Q = 0; Q < npoints; Q++) { 
            rho_a_[Q] = C_DDOT(nsigf, &temp_tens_[Q][0], 1, &points_[Q][0], 1);
            //printf(" Q = %d, rho = %14.10E\n", Q, rho_a_[Q]);
        }
    }
    if (do_density_gradient_) {

        // rho_a_x_
        // rho_a_x^Q = Da_mn * (phi_m^Q * phi_x_n^Q + phi_m_x^Q * phi_n^Q) 
        C_DGEMM('N', 'N', npoints, nsigf, nsigf, 1.0, &points_[0][0], nbf, &Da_[0][0], nbf, \
            0.0, &temp_tens_[0][0], nbf);

        for (int Q = 0; Q < npoints; Q++) 
            rho_a_x_[Q] = 2.0*C_DDOT(nsigf, &temp_tens_[Q][0], 1, &gradientsX_[Q][0], 1);

        // rho_a_y_
        // rho_a_y^Q = Da_mn * (phi_m^Q * phi_y_n^Q + phi_m_y^Q * phi_n^Q) 
        C_DGEMM('N', 'N', npoints, nsigf, nsigf, 1.0, &points_[0][0], nbf, &Da_[0][0], nbf, \
            0.0, &temp_tens_[0][0], nbf);

        for (int Q = 0; Q < npoints; Q++) 
            rho_a_y_[Q] = 2.0*C_DDOT(nsigf, &temp_tens_[Q][0], 1, &gradientsY_[Q][0], 1);

        // rho_a_z_
        // rho_a_z^Q = Da_mn * (phi_m^Q * phi_z_n^Q + phi_m_z^Q * phi_n^Q) 
        C_DGEMM('N', 'N', npoints, nsigf, nsigf, 1.0, &points_[0][0], nbf, &Da_[0][0], nbf, \
            0.0, &temp_tens_[0][0], nbf);

        for (int Q = 0; Q < npoints; Q++) 
            rho_a_z_[Q] = 2.0*C_DDOT(nsigf, &temp_tens_[Q][0], 1, &gradientsZ_[Q][0], 1);

        // gamma_aa_^Q = | \nabla rho_a | ^ 2
        for (int Q = 0; Q < npoints; Q++)
            gamma_aa_[Q] = rho_a_x_[Q]*rho_a_x_[Q] + \
                           rho_a_y_[Q]*rho_a_y_[Q] + \
                           rho_a_z_[Q]*rho_a_z_[Q];
                            
    }
    if (do_ke_density_) {
        // Now we need to fill C 
        for (int m = 0; m < nsigf; m++)
            for (int i = 0; i < docc[0]; i++)
                Ca_[m][i] = C->get(0,rel2abs_functions_[m],i);

        // tau_a_
        // tau_a^Q = \sum_i(occ) | C_mi \del \phi_m | ^ 2
        C_DGEMM('N', 'N', npoints, docc[0], nsigf, 1.0, &gradientsX_[0][0], nbf, &Ca_[0][0], nbf, \
            0.0, &temp_tens_[0][0], nbf);

        for (int Q = 0; Q < npoints; Q++)
            tau_a_[Q] = C_DDOT(docc[0], &temp_tens_[Q][0], 1, &temp_tens_[Q][0], 1);
        
        C_DGEMM('N', 'N', npoints, docc[0], nsigf, 1.0, &gradientsY_[0][0], nbf, &Ca_[0][0], nbf, \
            0.0, &temp_tens_[0][0], nbf);

        for (int Q = 0; Q < npoints; Q++)
            tau_a_[Q] += C_DDOT(docc[0], &temp_tens_[Q][0], 1, &temp_tens_[Q][0], 1);
        
        C_DGEMM('N', 'N', npoints, docc[0], nsigf, 1.0, &gradientsZ_[0][0], nbf, &Ca_[0][0], nbf, \
            0.0, &temp_tens_[0][0], nbf);

        for (int Q = 0; Q < npoints; Q++)
            tau_a_[Q] += C_DDOT(docc[0], &temp_tens_[Q][0], 1, &temp_tens_[Q][0], 1);
    }
}
void Properties::setToComputeDensity(bool v)
{
    if (!do_density_ && v) {
            rho_a_ = init_array(block_size_);
            rho_b_ = init_array(block_size_);
        }
        if (do_density_ && !v) {
            free(rho_a_);
            free(rho_b_);
        }

        do_density_ = v;
    if (v)
        setToComputePoints(true);
}
void Properties::setToComputeDensityGradient(bool v)
{
    if (!do_density_gradient_ && v) {
            gamma_aa_ = init_array(block_size_);
            gamma_ab_ = init_array(block_size_);
            gamma_bb_ = init_array(block_size_);
            rho_a_x_  = init_array(block_size_); 
            rho_a_y_  = init_array(block_size_); 
            rho_a_z_  = init_array(block_size_); 
            rho_b_x_  = init_array(block_size_); 
            rho_b_y_  = init_array(block_size_); 
            rho_b_z_  = init_array(block_size_); 
        }
        if (do_density_gradient_ && !v) {
            free(gamma_aa_);
            free(gamma_ab_);
            free(gamma_bb_);
            free(rho_a_x_);
            free(rho_a_y_);
            free(rho_a_z_);
            free(rho_b_x_);
            free(rho_b_y_);
            free(rho_b_z_);
        }
    do_density_gradient_ = v;
    if (v) {
        setToComputeGradients(true);
        setToComputePoints(true);
    }
}
void Properties::setToComputeKEDensity(bool v)
{
    if (!do_ke_density_ && v) {
            tau_a_ = init_array(block_size_);
            tau_b_ = init_array(block_size_);
        }
        if (do_ke_density_ && !v) {
            free(tau_a_);
            free(tau_b_);
        }
    do_ke_density_ = v;
    if (v)
        setToComputeGradients(true);
}
shared_ptr<Properties> Properties::get_testbed()
{
    int npoints = 23;
    Properties * props = new Properties(BasisSet::zero_ao_basis_set(), npoints);
    props->setTrueSize(npoints);
    props->setToComputeDensity(true);
    props->setToComputeDensityGradient(true);
    props->setToComputeKEDensity(true);

    int i = 0;

    double* rho_a_ = props->rho_a_;    
    double* rho_b_ = props->rho_b_;    
    double* gamma_aa_ = props->gamma_aa_;    
    double* gamma_ab_ = props->gamma_ab_;    
    double* gamma_bb_ = props->gamma_bb_;    
    double* tau_a_ = props->tau_a_;    
    double* tau_b_ = props->tau_b_;    

    rho_a_[i] = 0.17E+01; rho_b_[i] = 0.17E+01; gamma_aa_[i] = 0.81E-11; gamma_ab_[i] = 0.81E-11; gamma_bb_[i] = 0.81E-11; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++;  
    rho_a_[i] = 0.17E+01; rho_b_[i] = 0.17E+01; gamma_aa_[i] = 0.17E+01; gamma_ab_[i] = 0.17E+01; gamma_bb_[i] = 0.17E+01; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.15E+01; rho_b_[i] = 0.15E+01; gamma_aa_[i] = 0.36E+02; gamma_ab_[i] = 0.36E+02; gamma_bb_[i] = 0.36E+02; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.88E-01; rho_b_[i] = 0.88E-01; gamma_aa_[i] = 0.87E-01; gamma_ab_[i] = 0.87E-01; gamma_bb_[i] = 0.87E-01; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.18E+04; rho_b_[i] = 0.18E+04; gamma_aa_[i] = 0.55E+00; gamma_ab_[i] = 0.55E+00; gamma_bb_[i] = 0.55E+00; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.18E+04; rho_b_[i] = 0.18E+04; gamma_aa_[i] = 0.86E+04; gamma_ab_[i] = 0.86E+04; gamma_bb_[i] = 0.86E+04; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.16E+04; rho_b_[i] = 0.16E+04; gamma_aa_[i] = 0.37E+10; gamma_ab_[i] = 0.37E+10; gamma_bb_[i] = 0.37E+10; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.26E+00; rho_b_[i] = 0.26E+00; gamma_aa_[i] = 0.28E+00; gamma_ab_[i] = 0.28E+00; gamma_bb_[i] = 0.28E+00; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.53E+05; rho_b_[i] = 0.53E+05; gamma_aa_[i] = 0.96E+05; gamma_ab_[i] = 0.96E+05; gamma_bb_[i] = 0.96E+05; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.47E+05; rho_b_[i] = 0.47E+05; gamma_aa_[i] = 0.29E+14; gamma_ab_[i] = 0.29E+14; gamma_bb_[i] = 0.29E+14; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.15E+00; rho_b_[i] = 0.15E+00; gamma_aa_[i] = 0.16E+00; gamma_ab_[i] = 0.16E+00; gamma_bb_[i] = 0.16E+00; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.35E+01; rho_b_[i] = 0.00E+00; gamma_aa_[i] = 0.46E-10; gamma_ab_[i] = 0.00E+00; gamma_bb_[i] = 0.00E+00; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.35E+01; rho_b_[i] = 0.00E+00; gamma_aa_[i] = 0.34E+01; gamma_ab_[i] = 0.00E+00; gamma_bb_[i] = 0.00E+00; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.30E+01; rho_b_[i] = 0.00E+00; gamma_aa_[i] = 0.20E+03; gamma_ab_[i] = 0.00E+00; gamma_bb_[i] = 0.00E+00; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.58E-01; rho_b_[i] = 0.00E+00; gamma_aa_[i] = 0.47E-01; gamma_ab_[i] = 0.00E+00; gamma_bb_[i] = 0.00E+00; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.82E+02; rho_b_[i] = 0.81E+02; gamma_aa_[i] = 0.49E+07; gamma_ab_[i] = 0.49E+07; gamma_bb_[i] = 0.49E+07; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.39E+02; rho_b_[i] = 0.38E+02; gamma_aa_[i] = 0.81E+06; gamma_ab_[i] = 0.82E+06; gamma_bb_[i] = 0.82E+06; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.13E+00; rho_b_[i] = 0.95E-01; gamma_aa_[i] = 0.15E+00; gamma_ab_[i] = 0.18E+00; gamma_bb_[i] = 0.22E+00; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.78E-01; rho_b_[i] = 0.31E-01; gamma_aa_[i] = 0.41E-02; gamma_ab_[i] = 0.38E-02; gamma_bb_[i] = 0.36E-02; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.50E+02; rho_b_[i] = 0.49E+02; gamma_aa_[i] = 0.11E+06; gamma_ab_[i] = 0.11E+06; gamma_bb_[i] = 0.11E+06; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.40E+02; rho_b_[i] = 0.40E+02; gamma_aa_[i] = 0.99E+05; gamma_ab_[i] = 0.98E+05; gamma_bb_[i] = 0.98E+05; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.12E+00; rho_b_[i] = 0.10E+00; gamma_aa_[i] = 0.12E+00; gamma_ab_[i] = 0.13E+00; gamma_bb_[i] = 0.14E+00; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 
    rho_a_[i] = 0.48E-01; rho_b_[i] = 0.25E-01; gamma_aa_[i] = 0.46E-02; gamma_ab_[i] = 0.44E-02; gamma_bb_[i] = 0.41E-02; tau_a_[i] = 0.00E+00; tau_b_[i] = 0.00E+00; i++; 

    return (shared_ptr<Properties>) props;
}

}
