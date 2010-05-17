#include <libmints/basisset.h>
#include <libmints/properties.h>
#include <libmints/matrix.h>
#include <libmints/vector3.h>
#include <libmints/gridblock.h>
#include <cstdio>
#include <cmath>
//#include <psifiles.h>
#include <libciomr/libciomr.h>
//#include <libpsio/psio.h>
//#include <libchkpt/chkpt.hpp>
//#include <libipv1/ip_lib.h>
//#include <libiwl/iwl.hpp>
//#include <libqt/qt.h>
//#include <psifiles.h>

using namespace psi;

Properties::Properties(shared_ptr<BasisSet> _b, int _block_size): BasisPoints(_b, _block_size)
{
	do_mos_ = false;
	do_density_ = false;
	do_density_gradient_ = false;
	do_density_hessian_ = false;
	do_density_laplacian_ = false;
	do_ke_density_ = false;
        setToComputeDensity(true);
}
Properties::~Properties()
{
	setToComputeDensity(false);
        setToComputeDensityGradient(false);
        setToComputeDensityHessian(false);
        setToComputeDensityLaplacian(false);
        setToComputeKEDensity(false);
        int m[1];
        setToComputeMOs(false, m,0);
}
void Properties::get_RKS_GGA_Testbed()
{
    setToComputeDensity(true);
    setToComputeDensityGradient(true);

    density_[0] = 0.17E1; density_gradient_2_[0] = 0.81E-11; 
    density_[1] = 0.17E1; density_gradient_2_[1] = 0.17E1; 
    density_[2] = 0.15E1; density_gradient_2_[2] = 0.36E2; 
    density_[3] = 0.88E-1; density_gradient_2_[3]= 0.87E-1; 
    density_[4] = 0.18E4; density_gradient_2_[4] = 0.18E4; 
    density_[5] = 0.18E4; density_gradient_2_[5] = 0.86E4; 
    density_[6] = 0.16E4; density_gradient_2_[6] = 0.37E10; 
    density_[7] = 0.53E5; density_gradient_2_[7] = 0.96E5; 
    density_[8] = 0.47E5; density_gradient_2_[8] = 0.29E14;

    setTrueSize(RKS_GGA_TESTBED_SIZE_); 
    
}
void Properties::computeProperties(SharedGridBlock grid, SharedMatrix D, SharedMatrix C)
{
	int* rows = D->rowspi();
	int nirreps = D->nirreps();
	computePoints(grid);
    
        //double *xg = grid->getX();
        //double *yg = grid->getY();
        //double *zg = grid->getZ();
        int ntrue = grid->getTruePoints();

        for (int grid_index = 0; grid_index<ntrue; grid_index++) {
        // << BEGIN OUTER LOOP OVER POINTS
        
        //Vector3 v(xg[grid_index],yg[grid_index],zg[grid_index]); 
        
         //fprintf(oufile,"\nPoint:  <%14.10f,%14.10f,%14.10f>\n",v[0],v[1],v[2]);
	//for (int i=0; i<basis_->nbf(); i++)
	//{
		//fprintf(outfile,"Basis Function %d is %14.10f\n",i,points_[i]);
	//}
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
        }
        if (do_density_ && !v) {
            free(density_);
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
        }
        if (do_density_gradient_ && !v) {
            free(densityX_);
            free(densityY_);
            free(densityZ_);
            free(density_gradient_2_);
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
        }
        if (do_ke_density_ && !v) {
            free(ke_density_);
        }
	do_ke_density_ = v;
	if (v)
		setToComputeGradients(true);
}
