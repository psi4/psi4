#include <libutil/ref.h>
#include <libmints/basisset.h>
#include <libmints/properties.h>
#include <libmints/matrix.h>
#include <libmints/vector3.h>
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

Properties::Properties(shared_ptr<BasisSet> b): BasisPoints(b)
{
	do_mos_ = false;
	do_density_ = true;
	do_density_gradient_ = false;
	do_density_hessian_ = false;
	do_density_laplacian_ = false;
	do_ke_density_ = false;
}
Properties::~Properties()
{
	//Any frees
	if (do_mos_) {
		free(mos_);
		free(mo_inds_);
	}
}
void Properties::computeProperties(Vector3 v, SharedMatrix D, SharedMatrix C)
{
	int* rows = D->rowspi();
	int nirreps = D->nirreps();
	computePoints(v);
 //fprintf(oufile,"\nPoint:  <%14.10f,%14.10f,%14.10f>\n",v[0],v[1],v[2]);
	//for (int i=0; i<basis_->nbf(); i++)
	//{
		//fprintf(outfile,"Basis Function %d is %14.10f\n",i,points_[i]);
	//}
	if (do_mos_) {
		memset(mos_,0,nmo_*sizeof(double));
		for (int index = 0 ; index<nmo_; index++)
		{
			for (int k = 0; k<rows[0]; k++)
				mos_[index] += C->get(0,k,index)*points_[k];
		}
		//TODO: Symmetrize
	}

	if (do_density_){
		for (int h = 0; h<nirreps; h++) { 
			double temp = 0.0;
			for (int i = 0; i<rows[h]; i++)
				for (int j = 0; j<=i; j++)
					temp+=((i==j)?1.0:2.0)*D->get(h,i,j)*points_[i]*points_[j];
			density_ = temp;
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
					tempX+=((i==j)?1.0:2.0)*D->get(h,i,j)*(gradientsX_[i]*points_[j]+points_[i]*gradientsX_[j]);
					tempY+=((i==j)?1.0:2.0)*D->get(h,i,j)*(gradientsY_[i]*points_[j]+points_[i]*gradientsY_[j]);
					tempZ+=((i==j)?1.0:2.0)*D->get(h,i,j)*(gradientsZ_[i]*points_[j]+points_[i]*gradientsZ_[j]);
				}
			densityX_ = tempX;
			densityY_ = tempY;
			densityZ_ = tempZ;
			density_gradient_2_ = densityX_*densityX_+densityY_*densityY_+densityZ_*densityZ_;
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
}
void Properties::setToComputeMOs(bool v, int* inds, int nmo)
{
	if (v == false && do_mos_ == true) {
		free(mos_);
		free(mo_inds_);
	}
	if (v == true && do_mos_ == false) {
		nmo_ = nmo;
		mos_ = init_array(nmo);
		mo_inds_ = init_int_array(nmo);
		for (int k = 0; k< nmo; k++)
			mo_inds_[k] = inds[k];
	}
	//TODO Add case for reshape
	do_mos_ = v;
}
void Properties::setToComputeDensity(bool v)
{
	do_density_ = v;
	if (v)
		setToComputePoints(true);
}
void Properties::setToComputeDensityGradient(bool v)
{
	do_density_gradient_ = v;
	if (v) {
		setToComputeGradients(true);
		setToComputePoints(true);
	}
}
void Properties::setToComputeDensityHessian(bool v)
{
	do_density_hessian_ = v;
	if (v) {
		setToComputeHessians(true);
		setToComputeGradients(true);
		setToComputePoints(true);
	}
}
void Properties::setToComputeDensityLaplacian(bool v)
{
	do_density_laplacian_ = v;
	if (v) {
		setToComputeGradients(true);
		setToComputeLaplacians(true);
		setToComputePoints(true);
	}
}
void Properties::setToComputeKEDensity(bool v)
{
	do_ke_density_ = v;
	if (v)
		setToComputeGradients(true);
}
