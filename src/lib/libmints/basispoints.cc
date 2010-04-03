#include <libutil/ref.h>
#include <libmints/basisset.h>
#include <libmints/basispoints.h>
#include <libmints/gshell.h>
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

BasisPoints::BasisPoints(shared_ptr<BasisSet> bas)
{
	basis_ = bas;
	have_points_ = false;
	have_gradients_ = false;
	have_hessians_ = false;
	have_laplacians_ = false;
	do_points_ = true;
	do_gradients_ = false;
	do_hessians_ = false;
	do_laplacians_ = false;
	allocate();

	int lmax = basis_->max_am();
	Nam = init_array((lmax+1)*(lmax+2)>>1);
	ang = init_array((lmax+1)*(lmax+2)>>1);
	a = init_int_array((lmax+1)*(lmax+2)>>1);
	b = init_int_array((lmax+1)*(lmax+2)>>1);
	c = init_int_array((lmax+1)*(lmax+2)>>1);
		
	prims = init_array(basis_->max_nprimitive());

}
BasisPoints::~BasisPoints()
{
	do_points_ = false;
	do_gradients_ = false;
	do_hessians_ = false;
	do_laplacians_ = false;	
	release();

	free(Nam);
	free(ang);
	free(a);
	free(b);
	free(c);

	free(prims);
}
void BasisPoints::allocate()
{
	int lmax = basis_->max_am();
	if (do_points_ && !have_points_) {
		ao_points_ = init_array((lmax+1)*(lmax+2)>>1);
		points_ = init_array(basis_->nbf());
		have_points_ = true;
	}
	if (do_gradients_ && !have_gradients_)	{
		ao_gradX_ = init_array((lmax+1)*(lmax+2)>>1);
		ao_gradY_ = init_array((lmax+1)*(lmax+2)>>1);
		ao_gradZ_ = init_array((lmax+1)*(lmax+2)>>1);
		gradientsX_ = init_array(basis_->nbf());
		gradientsY_ = init_array(basis_->nbf());
		gradientsZ_ = init_array(basis_->nbf());
		have_gradients_ = true;
	}
	if (do_hessians_ && !have_hessians_) {
		ao_hessXY_ = init_array((lmax+1)*(lmax+2)>>1);
		ao_hessXZ_ = init_array((lmax+1)*(lmax+2)>>1);
		ao_hessYZ_ = init_array((lmax+1)*(lmax+2)>>1);
		ao_hessXX_ = init_array((lmax+1)*(lmax+2)>>1);
		ao_hessYY_ = init_array((lmax+1)*(lmax+2)>>1);
		ao_hessZZ_ = init_array((lmax+1)*(lmax+2)>>1);
		hessiansXY_ = init_array(basis_->nbf());
		hessiansXZ_ = init_array(basis_->nbf());
		hessiansYZ_ = init_array(basis_->nbf());
		hessiansXX_ = init_array(basis_->nbf());
		hessiansYY_ = init_array(basis_->nbf());
		hessiansZZ_ = init_array(basis_->nbf());
		have_hessians_ = true;
	}
	if (do_laplacians_ && !have_laplacians_) {
		ao_laplac_ = init_array((lmax+1)*(lmax+2)>>1);
		laplacians_ = init_array(basis_->nbf());
		have_laplacians_ = true;
	}
		
}
void BasisPoints::release()
{
	if (have_points_ && !do_points_) {
		free(ao_points_);
		free(points_);
		have_points_ = false;
	}
	if (have_gradients_ && !do_gradients_) {
		free(ao_gradX_);
		free(ao_gradY_);
		free(ao_gradZ_);
		free(gradientsX_);
		free(gradientsY_);
		free(gradientsZ_);
		have_gradients_ = false;
	}
	if (have_hessians_ && !do_hessians_) {
		free(ao_hessXY_);
		free(ao_hessXZ_);
		free(ao_hessYZ_);
		free(ao_hessXX_);
		free(ao_hessYY_);
		free(ao_hessZZ_);
		free(hessiansXY_);
		free(hessiansXZ_);
		free(hessiansYZ_);
		free(hessiansXX_);
		free(hessiansYY_);
		free(hessiansZZ_);
		have_hessians_ = false;
	}
	if (have_laplacians_ && !do_laplacians_) {
		free(ao_laplac_);
		free(laplacians_);
		have_laplacians_ = false;
	}
}
void BasisPoints::computePoints(Vector3 point)
{
	if (do_points_) {
		memset((void*)points_, 0, sizeof(double)*basis_->nbf());
	}
	if (do_gradients_) {
		memset((void*)gradientsX_, 0, sizeof(double)*basis_->nbf());
		memset((void*)gradientsY_, 0, sizeof(double)*basis_->nbf());
		memset((void*)gradientsZ_, 0, sizeof(double)*basis_->nbf());
	}
	if (do_hessians_) {
		memset((void*)hessiansXY_, 0, sizeof(double)*basis_->nbf());
		memset((void*)hessiansXZ_, 0, sizeof(double)*basis_->nbf());
		memset((void*)hessiansYZ_, 0, sizeof(double)*basis_->nbf());
		memset((void*)hessiansXX_, 0, sizeof(double)*basis_->nbf());
		memset((void*)hessiansYY_, 0, sizeof(double)*basis_->nbf());
		memset((void*)hessiansZZ_, 0, sizeof(double)*basis_->nbf());
	}
	if (do_laplacians_) {
		memset((void*)laplacians_, 0, sizeof(double)*basis_->nbf());
	}
	int lmax = basis_->max_am();
		

	//Run across all shells in the basis
	for (int P = 0; P<basis_->nshell(); P++)
	{
		if (do_points_) {
			memset((void*)ao_points_,0,sizeof(double)*(lmax+1)*(lmax+2)>>1);
		}
		if (do_gradients_) {
			memset((void*)ao_gradX_,0,sizeof(double)*(lmax+1)*(lmax+2)>>1);
			memset((void*)ao_gradY_,0,sizeof(double)*(lmax+1)*(lmax+2)>>1);
			memset((void*)ao_gradZ_,0,sizeof(double)*(lmax+1)*(lmax+2)>>1);
		}
		if (do_hessians_) {
			memset((void*)ao_hessXY_,0,sizeof(double)*(lmax+1)*(lmax+2)>>1);
			memset((void*)ao_hessXZ_,0,sizeof(double)*(lmax+1)*(lmax+2)>>1);
			memset((void*)ao_hessYZ_,0,sizeof(double)*(lmax+1)*(lmax+2)>>1);
			memset((void*)ao_hessXX_,0,sizeof(double)*(lmax+1)*(lmax+2)>>1);
			memset((void*)ao_hessYY_,0,sizeof(double)*(lmax+1)*(lmax+2)>>1);
			memset((void*)ao_hessZZ_,0,sizeof(double)*(lmax+1)*(lmax+2)>>1);
		}
		if (do_laplacians_) {
			memset((void*)ao_laplac_,0,sizeof(double)*(lmax+1)*(lmax+2)>>1);
		}
		//Get a reference to the current Gaussian Shell object
		shared_ptr<GaussianShell> shell = basis_->shell(P);
		//Get center of current shell
		Vector3 center = shell->center();
		//Compute displacement
		Vector3 r = point-center;
		//Compute distance squared
		double R2 = r.dot(r);
		
		//Evaluate the Gaussian part (gaussian normalization, coef and exponent) 
		//for each primitive
		for (int k = 0; k<shell->nprimitive(); k++) {
			prims[k] = shell->coef(k)*exp(-shell->exp(k)*R2);
			//printf("  Shell %d, Prim %d, Coef %14.10f, Alpha %14.10f, Val%14.10f\n",P,k,shell->coef(k),shell->exp(k),prims[k]);
		}
		int l = shell->am();


			//Get the AO exponents and angular norms for this shell
			//Multiply the angular part against the angular norm by m_c
			//TODO: Use the exp_ao array in BasisSet instead.
		
			int ao = 0;
			for (int i=0; i<=l; ++i) {
				int x = l-i;
				for (int j=0; j<=i; ++j) {
					int y = i-j;
					int z = j;
					a[ao] = x;
					b[ao] = y;
					c[ao] = z;
					Nam[ao] = shell->normalize(x,y,z);
					ang[ao] = pow(r[0],x)*pow(r[1],y)*pow(r[2],z);
					ao++;
				}
			}
			//Compute the AO basis functions in this shell 
			for (int mc = 0; mc<(l+1)*(l+2)>>1; mc++) {
				if (do_points_) {
					for (int k = 0; k<shell->nprimitive(); k++) {
						ao_points_[mc] += prims[k]; 	
					}
					ao_points_[mc] *= Nam[mc]*ang[mc];
					//fprintf(outfile, "AO basis function %d = %14.10f\n",mc,ao_points_[mc]); 				
				}
				if (do_gradients_) {
					for (int k = 0; k<shell->nprimitive(); k++) {
						if (a[mc] == 0)
							ao_gradX_[mc] += prims[k]*(-2.0*shell->exp(k)*r[0]*ang[mc]);
						else
							ao_gradX_[mc] += prims[k]*(a[mc]*pow(r[0],a[mc]-1)*pow(r[1],b[mc])*pow(r[2],c[mc])-2.0*shell->exp(k)*r[0]*ang[mc]);
						if (b[mc] == 0)
							ao_gradY_[mc] += prims[k]*(-2.0*shell->exp(k)*r[1]*ang[mc]);
						else
							ao_gradY_[mc] += prims[k]*(b[mc]*pow(r[0],a[mc])*pow(r[1],b[mc]-1)*pow(r[2],c[mc])-2.0*shell->exp(k)*r[1]*ang[mc]);
						if (c[mc] == 0)
							ao_gradZ_[mc] += prims[k]*(-2.0*shell->exp(k)*r[2]*ang[mc]);
						else
							ao_gradZ_[mc] += prims[k]*(c[mc]*pow(r[0],a[mc])*pow(r[1],b[mc])*pow(r[2],c[mc]-1)-2.0*shell->exp(k)*r[2]*ang[mc]);
					}
					ao_gradX_[mc] *= Nam[mc];
					ao_gradY_[mc] *= Nam[mc];
					ao_gradZ_[mc] *= Nam[mc];
					//fprintf(outfile, "AO basis gradient %d = <%14.10f,%14.10f,%14.10f>\n",mc,ao_gradX_[mc],ao_gradY_[mc],ao_gradZ_[mc]);
				}
				if (do_hessians_) {
					for (int k = 0; k<shell->nprimitive(); k++) {
						double q;
						q = 0.0;
						if (a[mc]>0 && b[mc]>0)
							q+= a[mc]*b[mc]*pow(r[0],a[mc]-1)*pow(r[1],b[mc]-1)*pow(r[2],c[mc]);
						if (a[mc]>0)
							q+= -2.0*shell->exp(k)*r[1]*a[mc]*pow(r[0],a[mc]-1)*pow(r[1],b[mc])*pow(r[2],c[mc]);
						if (b[mc]>0)
							q+= -2.0*shell->exp(k)*r[0]*b[mc]*pow(r[0],a[mc])*pow(r[1],b[mc]-1)*pow(r[2],c[mc]);
						q+= 4.0*(shell->exp(k))*(shell->exp(k))*r[0]*r[1]*ang[mc];
						ao_hessXY_[mc] += prims[k]*q;
					
						q = 0.0;
						if (a[mc]>0 && c[mc]>0)
							q+= a[mc]*c[mc]*pow(r[0],a[mc]-1)*pow(r[1],b[mc])*pow(r[2],c[mc]-1);
						if (a[mc]>0)
							q+= -2.0*shell->exp(k)*r[2]*a[mc]*pow(r[0],a[mc]-1)*pow(r[1],b[mc])*pow(r[2],c[mc]);
						if (c[mc]>0)
							q+= -2.0*shell->exp(k)*r[0]*c[mc]*pow(r[0],a[mc])*pow(r[1],b[mc])*pow(r[2],c[mc]-1);
						q+= 4.0*(shell->exp(k))*(shell->exp(k))*r[0]*r[2]*ang[mc];
						ao_hessXZ_[mc] += prims[k]*q;
					
						q = 0.0;
						if (b[mc]>0 && c[mc]>0)
							q+= b[mc]*c[mc]*pow(r[0],a[mc])*pow(r[1],b[mc]-1)*pow(r[2],c[mc]-1);
						if (b[mc]>0)
							q+= -2.0*shell->exp(k)*r[2]*b[mc]*pow(r[0],a[mc])*pow(r[1],b[mc]-1)*pow(r[2],c[mc]);
						if (c[mc]>0)
							q+= -2.0*shell->exp(k)*r[1]*c[mc]*pow(r[0],a[mc])*pow(r[1],b[mc])*pow(r[2],c[mc]-1);
						q+= 4.0*(shell->exp(k))*(shell->exp(k))*r[1]*r[2]*ang[mc];
						ao_hessYZ_[mc] += prims[k]*q;
					
						q = 0;
						if (a[mc]>1)
							q+= a[mc]*(a[mc]-1)*pow(r[0],a[mc]-2)*pow(r[1],b[mc])*pow(r[2],c[mc]);
						if (a[mc]>0)
							q+= -2.0*a[mc]*shell->exp(k)*ang[mc];
						q+= -2.0*shell->exp(k)*(a[mc]+1)*ang[mc]+4.0*(shell->exp(k))*(shell->exp(k))*r[0]*r[0]*ang[mc];
						ao_hessXX_[mc] += prims[k]*q;
					
						q = 0;
						if (b[mc]>1)
							q+= b[mc]*(b[mc]-1)*pow(r[0],a[mc])*pow(r[1],b[mc]-2)*pow(r[2],c[mc]);
						if (b[mc]>0)
							q+= -2.0*b[mc]*shell->exp(k)*ang[mc];
						q+= -2.0*shell->exp(k)*(b[mc]+1)*ang[mc]+4.0*(shell->exp(k))*(shell->exp(k))*r[1]*r[1]*ang[mc];
						ao_hessYY_[mc] += prims[k]*q;
				
						q = 0;
						if (c[mc]>1)
							q+= c[mc]*(c[mc]-1)*pow(r[0],a[mc])*pow(r[1],b[mc])*pow(r[2],c[mc]-2);
						if (a[mc]>0)
							q+= -2.0*c[mc]*shell->exp(k)*ang[mc];
						q+= -2.0*shell->exp(k)*(c[mc]+1)*ang[mc]+4.0*(shell->exp(k))*(shell->exp(k))*r[2]*r[2]*ang[mc];
						ao_hessZZ_[mc] += prims[k]*q;
						
					}
					ao_hessXY_[mc] *= Nam[mc];
					ao_hessXZ_[mc] *= Nam[mc];
					ao_hessYZ_[mc] *= Nam[mc];
					ao_hessXX_[mc] *= Nam[mc];
					ao_hessYY_[mc] *= Nam[mc];
					ao_hessZZ_[mc] *= Nam[mc];
					//fprintf(outfile, "AO basis Hessian %d = <%14.10f,%14.10f,%14.10f,%14.10f,%14.10f,%14.10f>\n",mc,ao_hessXY_[mc],ao_hessXZ_[mc],ao_jacobYZ[mc],ao_jacobXX[mc],ao_jacobYY[mc],ao_jacobZZ[mc]);
				}
			
				if (do_laplacians_) {
					if (do_hessians_)
						ao_laplac_[mc] = ao_hessXX_[mc]+ao_hessYY_[mc]+ao_hessZZ_[mc];
					else {
						for (int k = 0; k<shell->nprimitive(); k++) {
							ao_laplac_[mc] += prims[k]*((((a[mc]==0||r[0]==0.0)?0.0:-a[mc]/(r[0]*r[0]))-2.0*shell->exp(k))+(((a[mc]==0||r[0]==0.0)?0.0:a[mc]/r[0])-2.0*shell->exp(k)*r[0])*(((a[mc]==0||r[0]==0.0)?0.0:a[mc]/r[0])-2.0*shell->exp(k)*r[0]));
							ao_laplac_[mc] += prims[k]*((((b[mc]==0||r[1]==0.0)?0.0:-b[mc]/(r[1]*r[1]))-2.0*shell->exp(k))+(((b[mc]==0||r[1]==0.0)?0.0:b[mc]/r[1])-2.0*shell->exp(k)*r[1])*(((b[mc]==0||r[1]==0.0)?0.0:b[mc]/r[1])-2.0*shell->exp(k)*r[1]));
 	      						ao_laplac_[mc] += prims[k]*((((c[mc]==0||r[2]==0.0)?0.0:-c[mc]/(r[2]*r[2]))-2.0*shell->exp(k))+(((c[mc]==0||r[2]==0.0)?0.0:c[mc]/r[2])-2.0*shell->exp(k)*r[2])*(((c[mc]==0||r[2]==0.0)?0.0:c[mc]/r[2])-2.0*shell->exp(k)*r[2]));

						}
						ao_laplac_[mc] *= Nam[mc]*ang[mc];
					}
					//fprintf(outfile, "AO basis Laplacian %d = %14.10f\n",mc,ao_laplac_[mc]); 
				}
			}
			int start = shell->function_index();
			double trans_coef;
			int ind_ao;
			int ind_so;
			if (shell->is_pure()) {
				//AO -> SO (fairly easy)
				SphericalTransformIter trans(basis_->spherical_transform(shell->am()));
				for (trans.first(); trans.is_done();trans.next()) {
					trans_coef = trans.coef();
					ind_ao = trans.cartindex();
					ind_so = trans.pureindex();
					if (do_points_) {
						points_[ind_so+start] += trans_coef*ao_points_[ind_ao];
					}
					if (do_gradients_) {
						gradientsX_[ind_so+start] += trans_coef*ao_gradX_[ind_ao];
						gradientsY_[ind_so+start] += trans_coef*ao_gradY_[ind_ao];
						gradientsZ_[ind_so+start] += trans_coef*ao_gradZ_[ind_ao];
					}
					if (do_hessians_) {
						hessiansXY_[ind_so+start] += trans_coef*ao_hessXY_[ind_ao];
						hessiansXZ_[ind_so+start] += trans_coef*ao_hessXZ_[ind_ao];
						hessiansYZ_[ind_so+start] += trans_coef*ao_hessYZ_[ind_ao];
						hessiansXX_[ind_so+start] += trans_coef*ao_hessXX_[ind_ao];
						hessiansYY_[ind_so+start] += trans_coef*ao_hessYY_[ind_ao];
						hessiansZZ_[ind_so+start] += trans_coef*ao_hessZZ_[ind_ao];
					}
					if (do_laplacians_) {
						laplacians_[ind_so+start] += trans_coef*ao_laplac_[ind_ao];
					}
					//fprintf(out,"Transforming AO shell index %d to SO total index %d, c = %14.10f\n",ind_ao,ind_so+start,trans_coef); fflush(outfile);
				}			
			}
			else {
				//fprintf(outfile,"  AO -> AO transform\n");
				//AO -> AO (Easy)
				for (int i = 0; i<(l+1)*(l+2)>>1; i++) {
					if (do_points_) {
						points_[i+start] = ao_points_[i];
					}
					if (do_gradients_) {
						gradientsX_[i+start] += trans_coef*ao_gradX_[i];
						gradientsY_[i+start] += trans_coef*ao_gradY_[i];
						gradientsZ_[i+start] += trans_coef*ao_gradZ_[i];
					}
					if (do_hessians_) {
						hessiansXY_[i+start] += trans_coef*ao_hessXY_[i];
						hessiansXZ_[i+start] += trans_coef*ao_hessXZ_[i];
						hessiansYZ_[i+start] += trans_coef*ao_hessYZ_[i];
						hessiansXX_[i+start] += trans_coef*ao_hessXX_[i];
						hessiansYY_[i+start] += trans_coef*ao_hessYY_[i];
						hessiansZZ_[i+start] += trans_coef*ao_hessZZ_[i];
					}
					if (do_hessians_) {
						laplacians_[i+start] += trans_coef*ao_laplac_[i];
					}
				}
			}
		
	}
}
