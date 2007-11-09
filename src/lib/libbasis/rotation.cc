/*! \file 
    \ingroup (BASIS)
    \brief Enter brief description of file here 
*/

#include <stdexcept>
#if HAVE_CMATH
# include <cmath>
#else
# include <math.h>
#endif
extern "C" {
#include <libciomr/libciomr.h>
}
#include "rotation.h"

RotationOp::RotationOp(BasisSet* bs) :
  gnorm_(bs->max_am()), scdata_(bs->max_am())
{
  bs_ = bs;
  maxam_ = bs_->max_am();
  
  Rl_ = new PSI_FLOAT**[maxam_ + 1];
  for(int l=0; l<=maxam_; l++) {
    int nao = (l+1)*(l+2)/2;
    Rl_[l] = block_matrix(nao,nao);
  }
}

RotationOp::~RotationOp()
{
  for(int l=0; l<=maxam_; l++)
    free_block(Rl_[l]);
}

void RotationOp::check_am(int l) const
{
  if (l < 0 || l > maxam_)
    throw std::runtime_error("ERROR: RotationOp::check_am -- angular momentum out of range");
  if (l > 1)
    throw std::runtime_error("ERROR: RotationOp::check_am -- can only handle up to p functions at the moment");
}

void RotationOp::init_Rl(PSI_FLOAT** R)
{
  // s-type shell
  Rl_[0][0][0] = 1.0;

  // Loop over each type of shell (p, d, etc.)
  for(int am=1; am<=maxam_; am++) {

    int nbf = (am+1)*(am+2)/2;
    int bf = 0;

    // Loop over non-rotated cartesian AOs
    for(int i=0; i<=am; i++) {
      int l = am - i;
      for(int j=0; j<=i; j++) {
	int m = i - j;
	int n = j;

	// Zero out coefficients
	for(int k=0; k<nbf; k++)
	  Rl_[am][k][bf] = 0.0;

	for(int lx=0; lx<=l; lx++) {
	  int lyz = l - lx;
	  for(int ly=0; ly<=lyz; ly++) {
	    int lz = lyz - ly;
	    double pfac_l = scdata_.binomc(l,lx) * scdata_.binomc(lyz,ly) *
	      pow(R[0][0],lx) * pow(R[0][1],ly) * pow(R[0][2],lz);

	    for(int mx=0; mx<=m; mx++) {
	      int myz = m - mx;
	      for(int my=0; my<=myz; my++) {
		int mz = myz - my;
		double pfac_m = scdata_.binomc(m,mx) * scdata_.binomc(myz,my) *
		  pow(R[1][0],mx) * pow(R[1][1],my) * pow(R[1][2],mz);

		for(int nx=0; nx<=n; nx++) {
		  int nyz = n - nx;
		  for(int ny=0; ny<=nyz; ny++) {
		    int nz = nyz - ny;
		    double pfac_n = scdata_.binomc(n,nx) * scdata_.binomc(nyz,ny) * 
		      pow(R[2][0],nx) * pow(R[2][1],ny) * pow(R[2][2],nz);
		    
		    // indices of the rotated basis function
		    int ll = lx + mx + nx;
		    int mm = ly + my + ny;
		    int nn = am - ll - mm;
		    int am_ll = am - ll;
		    int bf_rot = am_ll*(am_ll+1)/2 + nn;

		    double pfac = pfac_l * pfac_m * pfac_n;
		    Rl_[am][bf_rot][bf] += pfac * gnorm_.norm(am,bf)/gnorm_.norm(am,bf_rot);
		    //		    Rl_[am][bf_rot][bf] += pfac;
		  }
		}
	      }
	    }
	  }
	}

	bf++;

      }
    }

  }

}

PSI_FLOAT** RotationOp::rotation_mat(PSI_FLOAT** R, int l)
{
  check_am(l);
  init_Rl(R);

  return Rl_[l];
}

PSI_FLOAT** RotationOp::full_rotation_mat(PSI_FLOAT** R)
{
  int nao = bs_->num_ao();
  PSI_FLOAT** TM = block_matrix(nao,nao);

  init_Rl(R);

  int nshell = bs_->num_shells();
  int ao_offset = 0;
  for(int s=0; s<nshell; s++) {
    int am = bs_->shell(s).am();
    int nao = bs_->shell(s).num_ao();
    for(int ao1=0; ao1<nao; ao1++)
      for(int ao2=0; ao2<nao; ao2++)
	TM[ao_offset+ao1][ao_offset+ao2] = Rl_[am][ao1][ao2];
    ao_offset += nao;
  }

  return TM;
}
