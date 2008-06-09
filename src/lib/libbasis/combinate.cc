/*! \file
    \ingroup BASIS
    \brief Enter brief description of file here 
*/

#include <libciomr/libciomr.h>
#include <stdexcept>
#include "combinate.h"

namespace psi {

StatCombData::StatCombData(int imax) :
  imax_(imax)
{
  Fact_ = new PSI_FLOAT[imax+1];
  Fact_[0] = 1.0;
  for(int i=1; i<=imax; i++)
    Fact_[i] = i*Fact_[i-1];

  Fact2_ = new PSI_FLOAT[imax+1];
  Fact2_[0] = 1.0;
  Fact2_[1] = 1.0;
  for(int i=2; i<=imax; i++)
    Fact2_[i] = i*Fact_[i-2];

  BinomC_ = block_matrix(imax+1,imax+1);
  for(int n=0; n<=imax; n++)
    for(int m=0; m<=n; m++)
      BinomC_[n][m] = Fact_[n]/(Fact_[m]*Fact_[n-m]);
}

StatCombData::~StatCombData()
{
  delete[] Fact_;
  delete[] Fact2_;
  free_block(BinomC_);
}

}
