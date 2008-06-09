/*! \file
    \ingroup BASIS
    \brief Enter brief description of file here 
*/

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <physconst.h>

#include "basisset.h"
#include "shell.h"
#include "overlap.h"

#define MAX(a,b) (a > b ? a : b)

namespace psi{

OverlapEngine::OverlapEngine(const BasisSet* bs1, const BasisSet* bs2) :
  bs1_(bs1), bs2_(bs2), overlap_recur_(bs1_->max_am(),bs2_->max_am()),
  gnorm_(MAX(bs1_->max_am(),bs2_->max_am()))
{
  int maxam1 = bs1_->max_am();
  int maxam2 = bs2_->max_am();
  
  int maxnao1 = (maxam1+1)*(maxam1+2)/2;
  int maxnao2 = (maxam2+1)*(maxam2+2)/2;
  buffer_ = new PSI_FLOAT[maxnao1*maxnao2];
}

OverlapEngine::~OverlapEngine()
{
  delete[] buffer_;
}

void OverlapEngine::compute_pair_(const GaussianShell& s1, const GaussianShell& s2)
{
  int am1 = s1.am();
  int am2 = s2.am();
  int nprim1 = s1.num_prims();
  int nprim2 = s2.num_prims();
  PSI_FLOAT A[3], B[3];
  A[0] = s1.origin(0);
  A[1] = s1.origin(1);
  A[2] = s1.origin(2);
  B[0] = s2.origin(0);
  B[1] = s2.origin(1);
  B[2] = s2.origin(2);

  // compute intermediates
  PSI_FLOAT AB2 = 0.0;
  for(int xyz=0; xyz<3; xyz++)
    AB2 += (A[xyz] - B[xyz]) * (A[xyz] - B[xyz]);

  // zero out buffer
  memset(buffer_, 0, s1.num_ao()*s2.num_ao()*sizeof(PSI_FLOAT));

  PSI_FLOAT** OIX = overlap_recur_.OIX();
  PSI_FLOAT** OIY = overlap_recur_.OIY();
  PSI_FLOAT** OIZ = overlap_recur_.OIZ();
  for(int p1=0; p1<nprim1; p1++) {
    PSI_FLOAT a1 = s1.exp(p1);
    PSI_FLOAT cc1 = s1.cc(0,p1);
    for(int p2=0; p2<nprim2; p2++) {
      PSI_FLOAT a2 = s2.exp(p2);
      PSI_FLOAT cc2 = s2.cc(0,p2);
      PSI_FLOAT gamma = a1 + a2;
      PSI_FLOAT oog = 1.0/gamma;

      PSI_FLOAT PA[3], PB[3];
      PSI_FLOAT P[3];

      P[0] = (a1*A[0] + a2*B[0])*oog;
      P[1] = (a1*A[1] + a2*B[1])*oog;
      P[2] = (a1*A[2] + a2*B[2])*oog;
      PA[0] = P[0] - A[0];
      PA[1] = P[1] - A[1];
      PA[2] = P[2] - A[2];
      PB[0] = P[0] - B[0];
      PB[1] = P[1] - B[1];
      PB[2] = P[2] - B[2];

      PSI_FLOAT over_pf = exp(-a1*a2*AB2*oog) * sqrt(M_PI*oog) * M_PI * oog * cc1 * cc2;
      
      // compute one-dimensional integrals
      overlap_recur_.compute(PA,PB,gamma,am1,am2);

      int ao12 = 0;
      for(int ii = 0; ii <= am1; ii++) {
	int l1 = am1 - ii;
	for(int jj = 0; jj <= ii; jj++) {
	  int m1 = ii - jj;
	  int n1 = jj;
	  /*--- create all am components of sj ---*/
	  for(int kk = 0; kk <= am2; kk++) {
	    int l2 = am2 - kk;
	    for(int ll = 0; ll <= kk; ll++) {
	      int m2 = kk - ll;
	      int n2 = ll;
	      
	      PSI_FLOAT x0 = OIX[l1][l2];
	      PSI_FLOAT y0 = OIY[m1][m2];
	      PSI_FLOAT z0 = OIZ[n1][n2];

	      buffer_[ao12++] += over_pf*x0*y0*z0;

	    }
	  }

	}
      }	       

    }
  }

  // Finish normalization
  int nao1 = s1.num_ao();
  int nao2 = s2.num_ao();
  int ao12 = 0;
  for(int ao1=0; ao1<nao1; ao1++)
    for(int ao2=0; ao2<nao2; ao2++)
      buffer_[ao12++] *= gnorm_.norm(am1,ao1)*gnorm_.norm(am2,ao2);

  // Done
}

void OverlapEngine::compute_shell_pair(int s1, int s2)
{
  compute_pair_(bs1_->shell(s1),bs2_->shell(s2));
}

PSI_FLOAT** OverlapEngine::compute_full_matrix()
{
  PSI_FLOAT** result = block_matrix(bs1_->num_ao(),bs2_->num_ao());

  int nshell1 = bs1_->num_shells();
  int nshell2 = bs2_->num_shells();

  for(int s1=0; s1<nshell1; s1++) {
    int ao_off1 = bs1_->first_ao(s1);
    int nao1 = bs1_->shell(s1).num_ao();

    for(int s2=0; s2<nshell2; s2++) {
      int ao_off2 = bs2_->first_ao(s2);
      int nao2 = bs2_->shell(s2).num_ao();

      compute_shell_pair(s1,s2);
      
      int ao12 = 0;
      for(int ao1=0; ao1<nao1; ao1++)
	for(int ao2=0; ao2<nao2; ao2++)
	  result[ao_off1+ao1][ao_off2+ao2] = buffer_[ao12++];

    }
  }

  return result;
}

}

