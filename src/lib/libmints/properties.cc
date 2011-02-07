#include "mints.h"
#include <cstdio>
#include <cmath>
#include <libciomr/libciomr.h>

using namespace psi;
using namespace boost;

Properties::Properties(shared_ptr<BasisSet> _b, int _block_size): BasisPoints(_b, _block_size)
{
    do_mos_ = false;
    do_density_ = false;
    do_density_gradient_ = false;
    do_density_hessian_ = false;
    do_density_laplacian_ = false;
    do_ke_density_ = false;
    do_electrostatic_ = false;
    setToComputeDensity(true);
}
Properties::~Properties()
{
    setToComputeDensity(false);
    setToComputeDensityGradient(false);
    setToComputeDensityHessian(false);
    setToComputeDensityLaplacian(false);
    setToComputeKEDensity(false);
    setToComputeElectrostatic(false);
    int m[1];
    setToComputeMOs(false, m,0);
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
void Properties::computeProperties(SharedGridBlock grid, SharedMatrix D, SharedMatrix C)
{
    int* rows = D->rowspi();
    int nirreps = D->nirrep();
    computePoints(grid);

        double *xg = grid->getX();
        double *yg = grid->getY();
        double *zg = grid->getZ();

        shared_ptr<Molecule> mol = basis_->molecule();

        int ntrue = grid->getTruePoints();

        for (int grid_index = 0; grid_index<ntrue; grid_index++) {
        // << BEGIN OUTER LOOP OVER POINTS

            Vector3 v(xg[grid_index],yg[grid_index],zg[grid_index]);
        //Vector3 v(xg[grid_index],yg[grid_index],zg[grid_index]);

         //fprintf(oufile,"\nPoint:  <%14.10f,%14.10f,%14.10f>\n",v[0],v[1],v[2]);
    //for (int i=0; i<basis_->nbf(); i++)
    //{
        //fprintf(outfile,"Basis Function %d is %14.10f\n",i,points_[i]);
    //}
        if (do_electrostatic_) {

            double esp = 0.0;

            //Nuclear contribution
            int Z;
            Vector3 r;
            double R;
            for (int A = 0; A < mol->natom(); A++) {
                Z = mol->Z(A);
                r = v-mol->xyz(A);
                R = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
                esp += Z/R;
            }

            //Electronic contribution
            //TODO

            electrostatic_[grid_index] = esp;
        }

    if (do_mos_) {
        memset(mos_[grid_index],0,nmo_*sizeof(double));
        for (int index = 0 ; index<nmo_; index++)
        {
            for (int k = 0; k<rows[0]; k++)
                mos_[grid_index][index] += C->get(0,k,index)*points_[grid_index][k];
        }
        //TODO: Symmetrize
    }

    if (do_density_){
        for (int h = 0; h<nirreps; h++) {
            double temp = 0.0;
            for (int i = 0; i<rows[h]; i++)
                for (int j = 0; j<=i; j++)
                    temp+=((i==j)?1.0:2.0)*D->get(h,i,j)*points_[grid_index][i]*points_[grid_index][j];
            density_[grid_index] = temp;
            //fprintf(oufile, "Density at <%14.10f,%14.10f,%14.10f> is %14.10f\n",v[0],v[1],v[2],density_);
        }
    }
    if (do_density_gradient_){
        for (int h = 0; h<nirreps; h++) {
            double tempX = 0.0;
            double tempY = 0.0;
            double tempZ = 0.0;
            for (int i = 0; i<rows[h]; i++)
                for (int j = 0; j<=i; j++) {
                    tempX+=((i==j)?1.0:2.0)*D->get(h,i,j)*(gradientsX_[grid_index][i]*points_[grid_index][j]+points_[grid_index][i]*gradientsX_[grid_index][j]);
                    tempY+=((i==j)?1.0:2.0)*D->get(h,i,j)*(gradientsY_[grid_index][i]*points_[grid_index][j]+points_[grid_index][i]*gradientsY_[grid_index][j]);
                    tempZ+=((i==j)?1.0:2.0)*D->get(h,i,j)*(gradientsZ_[grid_index][i]*points_[grid_index][j]+points_[grid_index][i]*gradientsZ_[grid_index][j]);
                }
            densityX_[grid_index] = tempX;
            densityY_[grid_index] = tempY;
            densityZ_[grid_index] = tempZ;
            density_gradient_2_[grid_index] = tempX*tempX+tempY*tempY+tempZ*tempZ;
        }
    }
    if (do_density_hessian_){
        //TODO
    }
    if (do_density_laplacian_){
        //TODO
    }
    if (do_ke_density_){
        //TODO
    }
        // <<END OUTER LOOP OVER POINTS
        }
}
void Properties::setToComputeMOs(bool v, int* inds, int nmo)
{
    if (v == false && do_mos_ == true) {
        free_block(mos_);
        free(mo_inds_);
    }
    if (v == true && do_mos_ == false) {
        nmo_ = nmo;
        mos_ = block_matrix(block_size_,nmo);
        mo_inds_ = init_int_array(nmo);
        for (int k = 0; k< nmo; k++)
            mo_inds_[k] = inds[k];
    }
    //TODO Add case for reshape
    do_mos_ = v;
}
void Properties::setToComputeDensity(bool v)
{
    if (!do_density_ && v) {
            density_ = init_array(block_size_);
            rho_a_ = init_array(block_size_);
            rho_b_ = init_array(block_size_);
        }
        if (do_density_ && !v) {
            free(density_);
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
            densityX_ = init_array(block_size_);
            densityY_ = init_array(block_size_);
            densityZ_ = init_array(block_size_);
            density_gradient_2_ = init_array(block_size_);
            gamma_aa_ = init_array(block_size_);
            gamma_ab_ = init_array(block_size_);
            gamma_bb_ = init_array(block_size_);
        }
        if (do_density_gradient_ && !v) {
            free(densityX_);
            free(densityY_);
            free(densityZ_);
            free(density_gradient_2_);
            free(gamma_aa_);
            free(gamma_ab_);
            free(gamma_bb_);
        }
    do_density_gradient_ = v;
    if (v) {
        setToComputeGradients(true);
        setToComputePoints(true);
    }
}
void Properties::setToComputeDensityHessian(bool v)
{
    if (!do_density_hessian_ && v) {
            densityXX_ = init_array(block_size_);
            densityYY_ = init_array(block_size_);
            densityZZ_ = init_array(block_size_);
            densityXY_ = init_array(block_size_);
            densityXZ_ = init_array(block_size_);
            densityYZ_ = init_array(block_size_);
        }
        if (do_density_hessian_ && !v) {
            free(densityXX_);
            free(densityYY_);
            free(densityZZ_);
            free(densityXY_);
            free(densityXZ_);
            free(densityYZ_);
        }
    do_density_hessian_ = v;
    if (v) {
        setToComputeHessians(true);
        setToComputeGradients(true);
        setToComputePoints(true);
    }
}
void Properties::setToComputeDensityLaplacian(bool v)
{
    if (!do_density_laplacian_ && v) {
            density_laplacian_ = init_array(block_size_);
        }
        if (do_density_laplacian_ && !v) {
            free(density_laplacian_);
        }
    do_density_laplacian_ = v;
    if (v) {
        setToComputeGradients(true);
        setToComputeLaplacians(true);
        setToComputePoints(true);
    }
}
void Properties::setToComputeKEDensity(bool v)
{
    if (!do_ke_density_ && v) {
            ke_density_ = init_array(block_size_);
            tau_a_ = init_array(block_size_);
            tau_b_ = init_array(block_size_);
        }
        if (do_ke_density_ && !v) {
            free(ke_density_);
            free(tau_a_);
            free(tau_b_);
        }
    do_ke_density_ = v;
    if (v)
        setToComputeGradients(true);
}
void Properties::setToComputeElectrostatic(bool v)
{
    if (!do_electrostatic_ && v) {
            electrostatic_ = init_array(block_size_);
            IntegralFactory factory(basis_,basis_,basis_,basis_);
            std::vector<SphericalTransform> transformer = factory.spherical_transform();
            e_ints_ = shared_ptr<ElectrostaticInt>(new ElectrostaticInt(transformer,basis_,basis_));
        }
        if (do_electrostatic_ && !v) {
            free(ke_density_);
        }
    do_electrostatic_ = v;
}
