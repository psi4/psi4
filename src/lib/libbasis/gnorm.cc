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
#include "gnorm.h"

GaussianNormalization::GaussianNormalization(int am) :
  maxam_(am)
{
  int dfmax = 2*maxam_ + 1;
  df_ = new PSI_FLOAT[dfmax];
  df_[0] = 1.0;
  if (maxam_ > 0) {
    df_[1] = 1.0;
    df_[2] = 1.0;
  }
  for(int i=3; i<dfmax; i++)
    df_[i] = (i-1)*df_[i-2];

  norm_coeffs_ = new PSI_FLOAT*[maxam_+1];
  for(int l=0; l<=maxam_; l++) {
    int nbf = (l+1)*(l+2)/2;
    norm_coeffs_[l] = new PSI_FLOAT[nbf];

    int bf = 0;
    for(int i=0; i<=l; i++) {
      int l1 = l - i;
      for(int j=0; j<=i; j++) {
	int m1 = i-j;
	int n1 = j;
	if (use_cca_integrals_standard)
	  norm_coeffs_[l][bf++] = 1.0;
	else
	  norm_coeffs_[l][bf++] = sqrt(df_[2*l]/(df_[2*l1]*df_[2*m1]*df_[2*n1]));
      }
    }
  }
}

GaussianNormalization::~GaussianNormalization()
{
  for(int l=0; l<=maxam_; l++)
    delete[] norm_coeffs_[l];
  delete[] norm_coeffs_;

  delete[] df_;
}

PSI_FLOAT GaussianNormalization::norm(int l, int bf)
{
  if (l > maxam_)
    throw std::runtime_error("ERROR: GaussianNormalization::norm -- l out of range");
  int bfmax = (l+1)*(l+2)/2;
  if (bf > bfmax)
    throw std::runtime_error("ERROR: GaussianNormalization::norm -- bf out of range");

  return norm_coeffs_[l][bf];
}

PSI_FLOAT GaussianNormalization::norm(int x, int y, int z)
{
  int l = x + y + z;
  if (l > maxam_)
    throw std::runtime_error("ERROR: GaussianNormalization::norm -- l out of range");

  if (l <= 1)
    return 1.0;
  int l_x = l - x;
  int bf = (l_x+1)*l_x/2 + z;

  return norm_coeffs_[l][bf];
}
