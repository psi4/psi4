/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "moinfo.h"
#include "mo_overlap.h"
#include "float.h"
#include "linalg.h"
#include <psi4-dec.h>

using namespace std;

// Wrap a,b indices into one composite index assuming S2 symmetry
#define INDEX2(a,b) ((a) > (b)) ? ( (((a)*(a+1)) >> 1) + (b) ) : ( (((b)*(b+1)) >> 1) + (a) )
// Wrap a>=b indices into one composite index assuming S2 symmetry
#define INDEX2_ORD(a,b) ( (((a)*(a+1))/2) + (b) )
// Wrap a>=b>=c indices into one composite index assuming S3 symmetry
#define INDEX3_ORD(a,b,c) ( ((a)*(((a)+4)*((a)-1)+6)/6) + (((b)*(b+1))/2) + (c) )

namespace psi { namespace DBOC {

extern MOInfo_t MOInfo;
extern void done(const char *);

double eval_rccsd_derwfn_overlap()
{
#if USE_MOINFO
  int ndocc = MOInfo.ndocc;
  int num_mo = MOInfo.num_mo;
#else
  abort();
#endif
  int nvirt = num_mo - ndocc;
  FLOAT **CSC_full = eval_S_alpha();
  FLOAT **CSC = create_matrix(ndocc,ndocc);
  int *tmpintvec = new int[ndocc];

  FLOAT **RmSp = create_matrix(ndocc,nvirt);
  FLOAT **RmDp = create_matrix(INDEX2_ORD(ndocc,0),INDEX2_ORD(nvirt,0));
  FLOAT **RmTp = create_matrix(INDEX3_ORD(ndocc,0,0),INDEX3_ORD(nvirt,0,0));

  //
  // Evaluate reference-reference overlap <Ref(-)|Ref(+)>
  //
  for(int i=0;i<ndocc;i++)
    for(int j=0;j<ndocc;j++)
      CSC[i][j] = CSC_full[i][j];
  //  C_DGETRF(ndocc,ndocc,&(CSC[0][0]),ndocc,tmpintvec);
  FLOAT sign;
  lu_decom(CSC, ndocc, tmpintvec, &sign);
  FLOAT deter_ref = 1.0;
  for(int i=0;i<ndocc;i++)
    deter_ref *= CSC[i][i];

  //
  // Evaluate all overlaps <Ref(-)|Singles(+)>
  //
  fprintf(outfile,"\n  -Overlap of Ref(-) with Singles(-):\n");
  for(int mo_i=0; mo_i<ndocc; mo_i++) {
    for(int mo_a=ndocc; mo_a<MOInfo.num_mo; mo_a++) {

      // Before generating the single restore the original (<Ref(-)|Ref(+)>) overlap matrix
      for(int i=0;i<ndocc;i++)
	for(int j=0;j<ndocc;j++)
	  CSC[i][j] = CSC_full[i][j];
      
      // Replace ith column with ath column
      for(int p=0; p<ndocc; p++)
	CSC[p][mo_i] = CSC_full[p][mo_a];

      // Compute the determinant
      //      C_DGETRF(ndocc,ndocc,&(CSC[0][0]),ndocc,tmpintvec);
      lu_decom(CSC, ndocc, tmpintvec, &sign);
      FLOAT deter1 = 1.0;
      for(int i=0;i<ndocc;i++)
	deter1 *= CSC[i][i];

      RmSp[mo_i][mo_a-ndocc] = deter_ref * deter1;
      fprintf(outfile,"  %d %d %20.15lf\n",mo_i,mo_a,RmSp[mo_i][mo_a-ndocc]);

    }
  }

  //
  // Evaluate all overlaps <Ref(-)|IJABDoubles(+)>
  //
  fprintf(outfile,"\n  -Overlap of Ref(-) with IJABDoubles(-):\n");
  for(int mo_i=0; mo_i<ndocc; mo_i++) {
    for(int mo_j=0; mo_j<mo_i; mo_j++) {
      for(int mo_a=ndocc, a=0; mo_a<MOInfo.num_mo; mo_a++,a++) {
	for(int mo_b=ndocc, b=0; mo_b<mo_a; mo_b++,b++) {

	  // Before generating the single restore the original (<Ref(-)|Ref(+)>) overlap matrix
	  for(int i=0;i<ndocc;i++)
	    for(int j=0;j<ndocc;j++)
	      CSC[i][j] = CSC_full[i][j];
	  
	  // Replace ith column with ath column
	  for(int p=0; p<ndocc; p++)
	    CSC[p][mo_i] = CSC_full[p][mo_a];
	  // Replace jth column with bth column
	  for(int p=0; p<ndocc; p++)
	    CSC[p][mo_j] = CSC_full[p][mo_b];
	  
	  // Compute the determinant
	  //	  C_DGETRF(ndocc,ndocc,&(CSC[0][0]),ndocc,tmpintvec);
	  lu_decom(CSC, ndocc, tmpintvec, &sign);
	  FLOAT deter1 = 1.0;
	  for(int i=0;i<ndocc;i++)
	    deter1 *= CSC[i][i];
	  
	  int ij = INDEX2_ORD(mo_i,mo_j);
	  int ab = INDEX2_ORD(a,b);
	  RmDp[ij][ab] = deter_ref * deter1;
	  fprintf(outfile,"  %d %d %d %d %20.15lf\n",mo_i,mo_j,mo_a,mo_b,RmDp[ij][ab]);
	  
	}
      }
    }
  }

  //
  // Evaluate all overlaps <Ref(-)|IJKABCTriples(+)>
  //
  fprintf(outfile,"\n  -Overlap of Ref(-) with IJKABCTriples(-):\n");
  for(int mo_i=0; mo_i<ndocc; mo_i++) {
    for(int mo_j=0; mo_j<mo_i; mo_j++) {
      for(int mo_k=0; mo_k<mo_j; mo_k++) {
	for(int mo_a=ndocc, a=0; mo_a<MOInfo.num_mo; mo_a++,a++) {
	  for(int mo_b=ndocc, b=0; mo_b<mo_a; mo_b++,b++) {
	    for(int mo_c=ndocc, c=0; mo_c<mo_b; mo_c++,c++) {

	      // Before generating the single restore the original (<Ref(-)|Ref(+)>) overlap matrix
	      for(int i=0;i<ndocc;i++)
		for(int j=0;j<ndocc;j++)
		  CSC[i][j] = CSC_full[i][j];
	  
	      // Replace ith column with ath column
	      for(int p=0; p<ndocc; p++)
		CSC[p][mo_i] = CSC_full[p][mo_a];
	      // Replace jth column with bth column
	      for(int p=0; p<ndocc; p++)
		CSC[p][mo_j] = CSC_full[p][mo_b];
	      // Replace jth column with bth column
	      for(int p=0; p<ndocc; p++)
		CSC[p][mo_k] = CSC_full[p][mo_c];
	  
	      // Compute the determinant
	      // C_DGETRF(ndocc,ndocc,&(CSC[0][0]),ndocc,tmpintvec);
	      lu_decom(CSC, ndocc, tmpintvec, &sign);
	      FLOAT deter1 = 1.0;
	      for(int i=0;i<ndocc;i++)
		deter1 *= CSC[i][i];
	      
	      int ijk = INDEX3_ORD(mo_i,mo_j,mo_k);
	      int abc = INDEX3_ORD(a,b,c);
	      RmTp[ijk][abc] = deter_ref * deter1;
	      fprintf(outfile,"  %d %d %d %d %d %d %20.15lf\n",mo_i,mo_j,mo_k,mo_a,mo_b,mo_c,RmTp[ijk][abc]);
	      
	    }
	  }
	}
      }
    }
  }

  delete[] tmpintvec;
  delete_matrix(CSC);
  delete_matrix(CSC_full);
  delete_matrix(RmTp);
  delete_matrix(RmDp);
  delete_matrix(RmSp);
  return (double)deter_ref*deter_ref;
}

}} // namespace psi::DBOC
