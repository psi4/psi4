/*! \file shell_block_matrix.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cstdlib>
#include<cstdio>
#include<libciomr/libciomr.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

namespace psi { namespace CINTS {

double ****init_shell_block_matrix()
{
  int si, sj;
  int ni, nj;
  double ****data;

  data = (double****) malloc(BasisSet.num_shells*sizeof(double***));
  for(si=0;si<BasisSet.num_shells;si++) {
    data[si] = (double***) malloc(BasisSet.num_shells*sizeof(double**));
    ni = ioff[BasisSet.shells[si].am];
    for(sj=0;sj<BasisSet.num_shells;sj++) {
      nj = ioff[BasisSet.shells[sj].am];
      data[si][sj] = block_matrix(ni,nj);
    }
  }

  return data;
}

void free_shell_block_matrix(double**** data)
{
  int si, sj;

  for(si=0;si<BasisSet.num_shells;si++) {
    for(sj=0;sj<BasisSet.num_shells;sj++) {
      free_block(data[si][sj]);
    }
    free(data[si]);
  }
  free(data);

  return;
}
  

void shell_block_to_block(double**** shell_block, double** block)
{
  int si, sj;
  int ni, nj;
  int i, j, ioffset, joffset;
  double **sb;

  for(si=0;si<BasisSet.num_shells;si++) {
    ni = ioff[BasisSet.shells[si].am];
    ioffset = BasisSet.shells[si].fao-1;
    for(sj=0;sj<BasisSet.num_shells;sj++) {
      nj = ioff[BasisSet.shells[sj].am];
      joffset = BasisSet.shells[sj].fao-1;
      sb = shell_block[si][sj];
      for(i=0;i<ni;i++)
	for(j=0;j<nj;j++)
	  block[i+ioffset][j+joffset] = sb[i][j];
    }
  }

  return;
}


void GplusGt(double**** G, double**** Gsym)
{
  int si, sj;
  int ni, nj;
  int i, j;
  double **sb, **sb_t;

  for(si=0;si<BasisSet.num_shells;si++) {
    ni = ioff[BasisSet.shells[si].am];
    for(sj=0;sj<BasisSet.num_shells;sj++) {
      nj = ioff[BasisSet.shells[sj].am];
      sb = G[si][sj];
      /*--- Symmetrize off-diagonal blocks ---*/
      if (1) {
	sb_t = G[sj][si];
	for(i=0;i<ni;i++)
	  for(j=0;j<nj;j++)
	    Gsym[si][sj][i][j] = sb[i][j] + sb_t[j][i];
      }
      /*--- Symmetrize the diagonal blocks ---*/
      else {
	for(i=0;i<ni;i++) {
	  for(j=0;j<i;j++)
	    Gsym[si][si][i][j] = sb[i][j] + sb[j][i];
	  Gsym[si][si][i][i] = sb[i][i];
	}
      }
    }
  }
  return;
}

};};
