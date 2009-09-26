/*! \file
    \ingroup MVO
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <libiwl/iwl.h>
#include "MOInfo.h"
#include "params.h"
#include "globals.h"


/* First definitions of globals */
extern "C" {
  extern FILE *infile, *outfile;
}

namespace psi { namespace mvo {

extern int *ioff;
extern struct MOInfo moinfo;
extern struct Params params;

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void form_fock_full(double **F)
{
  double noei;   /* Number of one-electron integrals */
  double ntei;   /* Number of two-electron integrals */
  double *oei;   /* Array of one-electron integrals */
  double *tei;   /* Array of two-electron integrals */
  double *qts;   /* Orbital indexing in QTS ordering */
  double tval;   /* temporary variable */
  int irrep;     /* index for irrep */
  int p;         /* general index */
  int q;         /* general index */
  int pq;        /* compound general, general index */
  int d;         /* occupied orbital index (relative Pitzer) */
  int i;         /* occupied orbital index */
  int ii;
  int pi;
  int iq;
  int pqii;
  int piiq;
  
  noei = (moinfo.nmo*(moinfo.nmo+1))/2;
  oei = init_array(noei);

  ntei = (noei*(noei+1))/2;
  tei = init_array(ntei);
  
  //moinfo.F = init_matrix(moinfo.nmo, moinfo.nmo);

  /* The integrals are in Pitzer ordering */ 
  iwl_rdone(PSIF_OEI, PSIF_MO_FZC, oei, noei, 0, 0, outfile);
  iwl_rdtwo(PSIF_MO_TEI, tei, ioff, ntei, moinfo.nfzc, moinfo.nfzv, 0, outfile);

  for(p=0,pq=0; p<moinfo.nmo; p++) {
    for(q=0; q<=p; q++,pq++) {
      F[p][q] = oei[pq];
      tval = 0.0;
      for (irrep=0; irrep<moinfo.nirreps; irrep++) {
        for (d=0; d<moinfo.clsdpi[irrep]; d++) {
          i = moinfo.order[moinfo.first[irrep]+d]; 
	  ii = INDEX(i,i);
          pi = INDEX(p,i);
          iq = INDEX(i,q);
          pqii = INDEX(pq,ii);
          piiq = INDEX(pi,iq);
	
          tval +=  2.0 * tei[pqii] - tei[piiq];
	}
      }
      F[p][q] += tval;
      if (p!=q) F[q][p] = F[p][q];
    }
  }

  free(oei);
  free(tei);
}

}} // end namespace psi::mvo

