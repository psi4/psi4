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

BasisPoints::BasisPoints(shared_ptr<BasisSet> b)
{
	basis_ = b;
	have_points_ = false;
	have_gradients_ = false;
	have_hessians_ = false;
	have_laplacians_ = false;
	do_points_ = true;
	do_gradients_ = false;
	do_hessians_ = false;
	do_laplacians_ = false;
	allocate();
}
BasisPoints::~BasisPoints()
{
	do_points_ = false;
	do_gradients_ = false;
	do_hessians_ = false;
	do_laplacians_ = false;	
	release();
}
void BasisPoints::allocate()
{
	if (do_points_ && !have_points_) {
		points_ = init_array(basis_->nbf());
		have_points_ = true;
	}
	if (do_gradients_ && !have_gradients_)	{
		gradientsX_ = init_array(basis_->nbf());
		gradientsY_ = init_array(basis_->nbf());
		gradientsZ_ = init_array(basis_->nbf());
		have_gradients_ = true;
	}
	if (do_hessians_ && !have_hessians_) {
		hessiansXY_ = init_array(basis_->nbf());
		hessiansXZ_ = init_array(basis_->nbf());
		hessiansYZ_ = init_array(basis_->nbf());
		hessiansXX_ = init_array(basis_->nbf());
		hessiansYY_ = init_array(basis_->nbf());
		hessiansZZ_ = init_array(basis_->nbf());
		have_hessians_ = true;
	}
	if (do_laplacians_ && !have_laplacians_) {
		laplacians_ = init_array(basis_->nbf());
		have_laplacians_ = true;
	}
}
void BasisPoints::release()
{
	if (have_points_ && !do_points_) {
		free(points_);
		have_points_ = false;
	}
	if (have_gradients_ && !do_gradients_) {
		free(gradientsX_);
		free(gradientsY_);
		free(gradientsZ_);
		have_gradients_ = false;
	}
	if (have_hessians_ && !do_hessians_) {
		free(hessiansXY_);
		free(hessiansXZ_);
		free(hessiansYZ_);
		free(hessiansXX_);
		free(hessiansYY_);
		free(hessiansZZ_);
		have_hessians_ = false;
	}
	if (have_laplacians_ && !do_laplacians_) {
		free(laplacians_);
		have_laplacians_ = false;
	}
}
void BasisPoints::computePoints(Vector3 point)
{
	if (do_points_) {
		memset(points_, 0.0, sizeof(double)*basis_->nbf());
	}
	if (do_gradients_) {
		memset(gradientsX_, 0.0, sizeof(double)*basis_->nbf());
		memset(gradientsY_, 0.0, sizeof(double)*basis_->nbf());
		memset(gradientsZ_, 0.0, sizeof(double)*basis_->nbf());
	}
	if (do_hessians_) {
		memset(hessiansXY_, 0.0, sizeof(double)*basis_->nbf());
		memset(hessiansXZ_, 0.0, sizeof(double)*basis_->nbf());
		memset(hessiansYZ_, 0.0, sizeof(double)*basis_->nbf());
		memset(hessiansXX_, 0.0, sizeof(double)*basis_->nbf());
		memset(hessiansYY_, 0.0, sizeof(double)*basis_->nbf());
		memset(hessiansZZ_, 0.0, sizeof(double)*basis_->nbf());
	}
	if (do_laplacians_) {
		memset(laplacians_, 0.0, sizeof(double)*basis_->nbf());
	}
	//Run across all shells in the basis
	for (int P = 0; P<basis_->nshell(); P++)
	{
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
		double* prims = init_array(shell->nprimitive());
		for (int k = 0; k<shell->nprimitive(); k++) {
			prims[k] = shell->coef(k)*exp(-shell->exp(k)*R2);
			//printf("  Shell %d, Prim %d, Coef %14.10f, Alpha %14.10f, Val%14.10f\n",P,k,shell->coef(k),shell->exp(k),prims[k]);
		}
		int l = shell->am();

		double* ao_points, *ao_gradX, *ao_gradY, *ao_gradZ;
		double* ao_hessXY, *ao_hessXZ, *ao_hessYZ, *ao_hessXX, *ao_hessYY, *ao_hessZZ, *ao_laplac; 
		if (do_points_)
			ao_points = init_array((l+1)*(l+2)>>1);
		if (do_gradients_) {
			ao_gradX = init_array((l+1)*(l+2)>>1);
			ao_gradY = init_array((l+1)*(l+2)>>1);
			ao_gradZ = init_array((l+1)*(l+2)>>1);
		}
		if (do_hessians_) {
			ao_hessXY = init_array((l+1)*(l+2)>>1);
			ao_hessXZ = init_array((l+1)*(l+2)>>1);
			ao_hessYZ = init_array((l+1)*(l+2)>>1);
			ao_hessXX = init_array((l+1)*(l+2)>>1);
			ao_hessYY = init_array((l+1)*(l+2)>>1);
			ao_hessZZ = init_array((l+1)*(l+2)>>1);
			}
		if (do_laplacians_) 
			ao_laplac = init_array((l+1)*(l+2)>>1);
		
			//Get the AO exponents and angular norms for this shell
			//Multiply the angular part against the angular norm by m_c
			//TODO: Use the exp_ao array in BasisSet instead.
			double* Nam = init_array((l+1)*(l+2)>>1);
			double* ang = init_array((l+1)*(l+2)>>1);
			int* a = init_int_array((l+1)*(l+2)>>1);
			int* b = init_int_array((l+1)*(l+2)>>1);
			int* c = init_int_array((l+1)*(l+2)>>1);
		
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
						ao_points[mc] += prims[k]; 	
					}
					ao_points[mc] *= Nam[mc]*ang[mc];
					//fprintf(outfile, "AO basis function %d = %14.10f\n",mc,ao_points[mc]); 				
				}
				if (do_gradients_) {
					for (int k = 0; k<shell->nprimitive(); k++) {
						if (a[mc] == 0)
							ao_gradX[mc] += prims[k]*(-2.0*shell->exp(k)*r[0]*ang[mc]);
						else
							ao_gradX[mc] += prims[k]*(a[mc]*pow(r[0],a[mc]-1)*pow(r[1],b[mc])*pow(r[2],c[mc])-2.0*shell->exp(k)*r[0]*ang[mc]);
						if (b[mc] == 0)
							ao_gradY[mc] += prims[k]*(-2.0*shell->exp(k)*r[1]*ang[mc]);
						else
							ao_gradY[mc] += prims[k]*(b[mc]*pow(r[0],a[mc])*pow(r[1],b[mc]-1)*pow(r[2],c[mc])-2.0*shell->exp(k)*r[1]*ang[mc]);
						if (c[mc] == 0)
							ao_gradZ[mc] += prims[k]*(-2.0*shell->exp(k)*r[2]*ang[mc]);
						else
							ao_gradZ[mc] += prims[k]*(c[mc]*pow(r[0],a[mc])*pow(r[1],b[mc])*pow(r[2],c[mc]-1)-2.0*shell->exp(k)*r[2]*ang[mc]);
					}
					ao_gradX[mc] *= Nam[mc];
					ao_gradY[mc] *= Nam[mc];
					ao_gradZ[mc] *= Nam[mc];
					//fprintf(outfile, "AO basis gradient %d = <%14.10f,%14.10f,%14.10f>\n",mc,ao_gradX[mc],ao_gradY[mc],ao_gradZ[mc]);
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
						ao_hessXY[mc] += prims[k]*q;
					
						q = 0.0;
						if (a[mc]>0 && c[mc]>0)
							q+= a[mc]*c[mc]*pow(r[0],a[mc]-1)*pow(r[1],b[mc])*pow(r[2],c[mc]-1);
						if (a[mc]>0)
							q+= -2.0*shell->exp(k)*r[2]*a[mc]*pow(r[0],a[mc]-1)*pow(r[1],b[mc])*pow(r[2],c[mc]);
						if (c[mc]>0)
							q+= -2.0*shell->exp(k)*r[0]*c[mc]*pow(r[0],a[mc])*pow(r[1],b[mc])*pow(r[2],c[mc]-1);
						q+= 4.0*(shell->exp(k))*(shell->exp(k))*r[0]*r[2]*ang[mc];
						ao_hessXZ[mc] += prims[k]*q;
					
						q = 0.0;
						if (b[mc]>0 && c[mc]>0)
							q+= b[mc]*c[mc]*pow(r[0],a[mc])*pow(r[1],b[mc]-1)*pow(r[2],c[mc]-1);
						if (b[mc]>0)
							q+= -2.0*shell->exp(k)*r[2]*b[mc]*pow(r[0],a[mc])*pow(r[1],b[mc]-1)*pow(r[2],c[mc]);
						if (c[mc]>0)
							q+= -2.0*shell->exp(k)*r[1]*c[mc]*pow(r[0],a[mc])*pow(r[1],b[mc])*pow(r[2],c[mc]-1);
						q+= 4.0*(shell->exp(k))*(shell->exp(k))*r[1]*r[2]*ang[mc];
						ao_hessYZ[mc] += prims[k]*q;
					
						q = 0;
						if (a[mc]>1)
							q+= a[mc]*(a[mc]-1)*pow(r[0],a[mc]-2)*pow(r[1],b[mc])*pow(r[2],c[mc]);
						if (a[mc]>0)
							q+= -2.0*a[mc]*shell->exp(k)*ang[mc];
						q+= -2.0*shell->exp(k)*(a[mc]+1)*ang[mc]+4.0*(shell->exp(k))*(shell->exp(k))*r[0]*r[0]*ang[mc];
						ao_hessXX[mc] += prims[k]*q;
					
						q = 0;
						if (b[mc]>1)
							q+= b[mc]*(b[mc]-1)*pow(r[0],a[mc])*pow(r[1],b[mc]-2)*pow(r[2],c[mc]);
						if (b[mc]>0)
							q+= -2.0*b[mc]*shell->exp(k)*ang[mc];
						q+= -2.0*shell->exp(k)*(b[mc]+1)*ang[mc]+4.0*(shell->exp(k))*(shell->exp(k))*r[1]*r[1]*ang[mc];
						ao_hessYY[mc] += prims[k]*q;
				
						q = 0;
						if (c[mc]>1)
							q+= c[mc]*(c[mc]-1)*pow(r[0],a[mc])*pow(r[1],b[mc])*pow(r[2],c[mc]-2);
						if (a[mc]>0)
							q+= -2.0*c[mc]*shell->exp(k)*ang[mc];
						q+= -2.0*shell->exp(k)*(c[mc]+1)*ang[mc]+4.0*(shell->exp(k))*(shell->exp(k))*r[2]*r[2]*ang[mc];
						ao_hessZZ[mc] += prims[k]*q;
						
					}
					ao_hessXY[mc] *= Nam[mc];
					ao_hessXZ[mc] *= Nam[mc];
					ao_hessYZ[mc] *= Nam[mc];
					ao_hessXX[mc] *= Nam[mc];
					ao_hessYY[mc] *= Nam[mc];
					ao_hessZZ[mc] *= Nam[mc];
					//fprintf(outfile, "AO basis Hessian %d = <%14.10f,%14.10f,%14.10f,%14.10f,%14.10f,%14.10f>\n",mc,ao_hessXY[mc],ao_hessXZ[mc],ao_jacobYZ[mc],ao_jacobXX[mc],ao_jacobYY[mc],ao_jacobZZ[mc]);
				}
			
				if (do_laplacians_) {
					if (do_hessians_)
						ao_laplac[mc] = ao_hessXX[mc]+ao_hessYY[mc]+ao_hessZZ[mc];
					else {
						for (int k = 0; k<shell->nprimitive(); k++) {
							ao_laplac[mc] += prims[k]*((((a[mc]==0||r[0]==0.0)?0.0:-a[mc]/(r[0]*r[0]))-2.0*shell->exp(k))+(((a[mc]==0||r[0]==0.0)?0.0:a[mc]/r[0])-2.0*shell->exp(k)*r[0])*(((a[mc]==0||r[0]==0.0)?0.0:a[mc]/r[0])-2.0*shell->exp(k)*r[0]));
							ao_laplac[mc] += prims[k]*((((b[mc]==0||r[1]==0.0)?0.0:-b[mc]/(r[1]*r[1]))-2.0*shell->exp(k))+(((b[mc]==0||r[1]==0.0)?0.0:b[mc]/r[1])-2.0*shell->exp(k)*r[1])*(((b[mc]==0||r[1]==0.0)?0.0:b[mc]/r[1])-2.0*shell->exp(k)*r[1]));
 	      						ao_laplac[mc] += prims[k]*((((c[mc]==0||r[2]==0.0)?0.0:-c[mc]/(r[2]*r[2]))-2.0*shell->exp(k))+(((c[mc]==0||r[2]==0.0)?0.0:c[mc]/r[2])-2.0*shell->exp(k)*r[2])*(((c[mc]==0||r[2]==0.0)?0.0:c[mc]/r[2])-2.0*shell->exp(k)*r[2]));

						}
						ao_laplac[mc] *= Nam[mc]*ang[mc];
					}
					//fprintf(outfile, "AO basis Laplacian %d = %14.10f\n",mc,ao_laplac[mc]); 
				}
			}
			//TODO: AO Transformation	
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
						points_[ind_so+start] += trans_coef*ao_points[ind_ao];
					}
					if (do_gradients_) {
						gradientsX_[ind_so+start] += trans_coef*ao_gradX[ind_ao];
						gradientsY_[ind_so+start] += trans_coef*ao_gradY[ind_ao];
						gradientsZ_[ind_so+start] += trans_coef*ao_gradZ[ind_ao];
					}
					if (do_hessians_) {
						hessiansXY_[ind_so+start] += trans_coef*ao_hessXY[ind_ao];
						hessiansXZ_[ind_so+start] += trans_coef*ao_hessXZ[ind_ao];
						hessiansYZ_[ind_so+start] += trans_coef*ao_hessYZ[ind_ao];
						hessiansXX_[ind_so+start] += trans_coef*ao_hessXX[ind_ao];
						hessiansYY_[ind_so+start] += trans_coef*ao_hessYY[ind_ao];
						hessiansZZ_[ind_so+start] += trans_coef*ao_hessZZ[ind_ao];
					}
					if (do_laplacians_) {
						laplacians_[ind_so+start] += trans_coef*ao_laplac[ind_ao];
					}
					//fprintf(out,"Transforming AO shell index %d to SO total index %d, c = %14.10f\n",ind_ao,ind_so+start,trans_coef); fflush(outfile);
				}			
			}
			else {
				//fprintf(outfile,"  AO -> AO transform\n");
				//AO -> AO (Easy)
				for (int i = 0; i<(l+1)*(l+2)>>1; i++) {
					if (do_points_) {
						points_[i+start] = ao_points[i];
					}
					if (do_gradients_) {
						gradientsX_[i+start] += trans_coef*ao_gradX[i];
						gradientsY_[i+start] += trans_coef*ao_gradY[i];
						gradientsZ_[i+start] += trans_coef*ao_gradZ[i];
					}
					if (do_hessians_) {
						hessiansXY_[i+start] += trans_coef*ao_hessXY[i];
						hessiansXZ_[i+start] += trans_coef*ao_hessXZ[i];
						hessiansYZ_[i+start] += trans_coef*ao_hessYZ[i];
						hessiansXX_[i+start] += trans_coef*ao_hessXX[i];
						hessiansYY_[i+start] += trans_coef*ao_hessYY[i];
						hessiansZZ_[i+start] += trans_coef*ao_hessZZ[i];
					}
					if (do_hessians_) {
						laplacians_[i+start] += trans_coef*ao_laplac[i];
					}
				}
			}

	
			//Free ao_buffers
			if (do_points_)
				free(ao_points);
			if (do_gradients_) {
				free(ao_gradX);
				free(ao_gradY);
				free(ao_gradZ);
			}
			if (do_hessians_) {
				free(ao_hessXY);
				free(ao_hessXZ);
				free(ao_hessYZ);
				free(ao_hessXX);
				free(ao_hessYY);
				free(ao_hessZZ);
			}
			if (do_laplacians_)
				free(ao_laplac);


		//Free AO angular exponents and norms
		free(a);
		free(b);
		free(c);
		free(Nam);
		free(ang);

		
		free(prims);
	}
}
