/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "moinfo.h"
#include "mo_overlap.h"
#include "float.h"
#include "linalg.h"
#include "hfwfn.h"

namespace psi { namespace DBOC {

extern MOInfo_t MOInfo;
extern HFWavefunction* HFVectors[MAX_NUM_DISP];

extern void done(const char *);

double eval_uhf_derwfn_overlap(DisplacementIndex LDisp, DisplacementIndex RDisp)
{
  FLOAT **CSC_a = eval_S_alpha(LDisp,RDisp);
  FLOAT **CSC_b = eval_S_beta(LDisp,RDisp);

  int* clsdpi = HFVectors[LDisp]->clsdpi();
  int* openpi = HFVectors[LDisp]->openpi();
  int* orbspi = HFVectors[LDisp]->orbspi();
  int nirreps = HFVectors[LDisp]->nirreps();

#if USE_MOINFO
  int nalpha = MOInfo.nalpha;
  int nbeta = MOInfo.nbeta;
  int ndocc = nbeta;
#else
  int nalpha = HFVectors[LDisp]->nalpha();
  int nbeta = HFVectors[LDisp]->nbeta();
  int ndocc = nbeta;
#endif

  // Extract the alpha and beta blocks
  FLOAT **CSC_alpha = create_matrix(nalpha,nalpha);
  FLOAT **CSC_beta = create_matrix(nbeta,nbeta);
  int mo_offset1 = 0;
  int docc_offset1 = 0;
  int socc_offset1 = ndocc;
  for(int irrep1=0; irrep1<nirreps; irrep1++) {

    int ndocc1 = clsdpi[irrep1];
    int nsocc1 = openpi[irrep1];

    int mo_offset2 = 0;
    int docc_offset2 = 0;
    int socc_offset2 = ndocc;
    for(int irrep2=0; irrep2<nirreps; irrep2++) {

      int ndocc2 = clsdpi[irrep2];
      int nsocc2 = openpi[irrep2];

      for(int i=0;i<ndocc1;i++)
	for(int j=0;j<ndocc2;j++) {
	  CSC_alpha[i+docc_offset1][j+docc_offset2] = CSC_a[i+mo_offset1][j+mo_offset2];
	  CSC_beta[i+docc_offset1][j+docc_offset2] = CSC_b[i+mo_offset1][j+mo_offset2];
	}

      for(int i=0;i<ndocc1;i++)
	for(int j=0;j<nsocc2;j++) {
	  CSC_alpha[i+docc_offset1][j+socc_offset2] = CSC_a[i+mo_offset1][j+ndocc2+mo_offset2];
	}

      for(int i=0;i<nsocc1;i++)
	for(int j=0;j<ndocc2;j++) {
	  CSC_alpha[i+socc_offset1][j+docc_offset2] = CSC_a[i+mo_offset1+ndocc1][j+mo_offset2];
	}

      for(int i=0;i<nsocc1;i++)
	for(int j=0;j<nsocc2;j++) {
	  CSC_alpha[i+socc_offset1][j+socc_offset2] = CSC_a[i+mo_offset1+ndocc1][j+mo_offset2+ndocc2];
	}

      docc_offset2 += ndocc2;
      socc_offset2 += nsocc2;
      mo_offset2 += orbspi[irrep2];
    }

    docc_offset1 += ndocc1;
    socc_offset1 += nsocc1;
    mo_offset1 += orbspi[irrep1];
  }
  delete_matrix(CSC_a);
  delete_matrix(CSC_b);

  // Compute the overlap of alpha part
  int *tmpintvec = new int[nalpha];
  FLOAT sign;
  lu_decom(CSC_alpha, nalpha, tmpintvec, &sign);
  delete[] tmpintvec;
  FLOAT deter_a = 1.0;
  for(int i=0;i<nalpha;i++)
    deter_a *= CSC_alpha[i][i];
  deter_a = FABS(deter_a);
  delete_matrix(CSC_alpha);

  // Compute the overlap of beta part
  tmpintvec = new int[nbeta];
  lu_decom(CSC_beta, nbeta, tmpintvec, &sign);
  delete[] tmpintvec;
  FLOAT deter_b = 1.0;
  for(int i=0;i<nbeta;i++)
    deter_b *= CSC_beta[i][i];
  deter_b = FABS(deter_b);
  delete_matrix(CSC_beta);

  return (double)deter_a*deter_b;
}

}} // namespace psi::DBOC
