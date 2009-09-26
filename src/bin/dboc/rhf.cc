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
#include "defines.h"
#include "params.h"
#include "moinfo.h"
#include "mo_overlap.h"
#include "float.h"
#include "linalg.h"
#include "hfwfn.h"
#include <psi4-dec.h>

namespace psi { namespace DBOC {

extern MOInfo_t MOInfo;
extern Params_t Params;
extern HFWavefunction* HFVectors[MAX_NUM_DISP];

extern void done(const char *);

double eval_rhf_derwfn_overlap(DisplacementIndex LDisp, DisplacementIndex RDisp)
{
  FLOAT **CSC = eval_S_alpha(LDisp,RDisp);

  int* clsdpi = HFVectors[LDisp]->clsdpi();
  int* orbspi = HFVectors[LDisp]->orbspi();
  int nirreps = HFVectors[LDisp]->nirreps();
#if USE_MOINFO
  int ndocc = MOInfo.ndocc;
#else
  int ndocc = HFVectors[LDisp]->ndocc();
#endif

  // Extract the occupied block
  FLOAT **CSC_occ = create_matrix(ndocc,ndocc);
  int mo_offset1 = 0;
  int occ_offset1 = 0;
  for(int irrep1=0; irrep1<nirreps; irrep1++) {

    int nocc1 = clsdpi[irrep1];

    int mo_offset2 = 0;
    int occ_offset2 = 0;
    for(int irrep2=0; irrep2<nirreps; irrep2++) {

      int nocc2 = clsdpi[irrep2];

      for(int i=0;i<nocc1;i++)
	for(int j=0;j<nocc2;j++)
	  CSC_occ[i+occ_offset1][j+occ_offset2] = CSC[i+mo_offset1][j+mo_offset2];

      occ_offset2 += nocc2;
      mo_offset2 += orbspi[irrep2];
    }

    occ_offset1 += nocc1;
    mo_offset1 += orbspi[irrep1];
  }
  delete_matrix(CSC);

  if (Params.print_lvl > PrintLevels::print_contrib) {
    fprintf(outfile,"  +/- overlap in the basis of doubly-occupied MOs:\n");
    psi::DBOC::print_mat(CSC_occ, ndocc, ndocc, outfile);
  }

  // Compute the determinant
  int *tmpintvec = new int[ndocc];
  FLOAT sign;
  lu_decom(CSC_occ, ndocc, tmpintvec, &sign);
  delete[] tmpintvec;
  FLOAT deter1 = 1.0;
  for(int i=0;i<ndocc;i++)
    deter1 *= CSC_occ[i][i];
  deter1 = FABS(deter1);

  delete_matrix(CSC_occ);
  return (double)deter1*deter1;
}

}} // namespace psi::DBOC
