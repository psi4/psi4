/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libqt/slaterdset.h>
#include <psifiles.h>
#include "moinfo.h"
#include "params.h"
#include "float.h"
#include "linalg.h"
#include "mo_overlap.h"
#include "hfwfn.h"
#include "stringblocks.h"
#include "dets.h"
#include "ci_overlap.h"

using namespace std;
using namespace psi::dboc;

// Wrap a,b indices into one composite index assuming S2 symmetry
#define INDEX2(a,b) ((a) > (b)) ? ( (((a)*(a+1)) >> 1) + (b) ) : ( (((b)*(b+1)) >> 1) + (a) )
// Wrap a>=b indices into one composite index assuming S2 symmetry
#define INDEX2_ORD(a,b) ( (((a)*(a+1))/2) + (b) )
// Wrap a>=b>=c indices into one composite index assuming S3 symmetry
#define INDEX3_ORD(a,b,c) ( ((a)*(((a)+4)*((a)-1)+6)/6) + (((b)*(b+1))/2) + (c) )

// Set to 1 to use the new threaded code
#define USE_CI_OVERLAP 1
// Set to 1 to reduce I/O at the expense of more computation
#define LOOP_OVER_BLOCKS 1

namespace psi { namespace dboc {

extern MOInfo_t MOInfo;
extern Params_t Params;
extern char *CI_Vector_Labels[MAX_NUM_DISP];
extern HFWavefunction* HFVectors[MAX_NUM_DISP];
extern void done(const char *);
extern void mo_maps(short int**, short int**);
extern "C" FILE *outfile;

double eval_rci_derwfn_overlap(DisplacementIndex LDisp, DisplacementIndex RDisp)
{
#if USE_MOINFO
  int ndocc = MOInfo.ndocc;
#else
  int ndocc = HFVectors[LDisp]->ndocc();
#endif
  FLOAT **CSC_full = eval_S_alpha(LDisp,RDisp);
  FLOAT **CSC = create_matrix(ndocc,ndocc);
  int *tmpintvec = new int[ndocc];

  // Read in CI vectors
  SlaterDetVector *vecm, *vecp;
  slaterdetvector_read(PSIF_CIVECT,CI_Vector_Labels[RDisp],&vecm);
  slaterdetvector_read(PSIF_CIVECT,CI_Vector_Labels[LDisp],&vecp);

  // Compute overlap between strings for alpha spin case (beta is the same)
  StringSet *ssetm;
  ssetm = vecm->sdset->alphastrings;
  short int* fzc_occ = ssetm->fzc_occ;
  int nstr_a = ssetm->size;
  int nfzc = ssetm->nfzc;
  int nact = ndocc - nfzc;

  // Overlap matrix can be very large, especially for low excitation levels. 
  // Must handle a block of strings at a time.
  // How many strings can handle at once? Need to hold 2 blocks of overlap 
  // at the same time.
  const double nstr_per_block_float = 
    std::sqrt((double)Params.memory/(2*sizeof(FLOAT)));
  int nstr_per_block = (int)std::floor(nstr_per_block_float);
  {
    if (nstr_per_block > nstr_a)
      nstr_per_block = nstr_a;
    else {
      const int nblks = (nstr_a + nstr_per_block - 1)/ nstr_per_block;
      nstr_per_block = (nstr_a + nblks - 1) / nblks;
    }
  }
  StringBlocks strblks_a(nstr_a,nstr_per_block);
  const int nblks_a = strblks_a.nblocks();
  StringBlockedMatrix ovlp_a(&strblks_a,&strblks_a,"S_a");
  double** S_a = ovlp_a.buffer();

  // WARNING: Assume the order of strings is the same for - and + displacements!!!
  // Loop over blocks of + and - displacement strings
  for(int bp=0; bp<nblks_a; ++bp) {
    const int jp_begin = strblks_a.begin(bp);
    const int jp_end = strblks_a.end(bp);

    for(int bm=0; bm<nblks_a; ++bm) {
      const int im_begin = strblks_a.begin(bm);
      const int im_end = strblks_a.end(bm);

      // loop over strings in each block
      for(int jp=jp_begin; jp<=jp_end; jp++) {
	String *str_j = &ssetm->strings[jp];
	
	for(int im=im_begin; im<=im_end; im++) {
	  String *str_i = &ssetm->strings[im];

	  for(int j=0;j<nact;j++)
	    for(int i=0;i<nact;i++)
	      CSC[j+nfzc][i+nfzc] = CSC_full[str_j->occ[j]][str_i->occ[i]];

	  // frozen orbitals need to be mapped to pitzer order manually
	  for(int i=0; i<nfzc; i++) {
	    int ii = fzc_occ[i];
	    for(int j=0; j<nact; j++)
	      CSC[j+nfzc][i] = CSC_full[str_j->occ[j]][ii];
	  }

	  for(int j=0; j<nfzc; j++) {
	    int jj = fzc_occ[j];
	    for(int i=0; i<nact; i++)
	      CSC[j][i+nfzc] = CSC_full[jj][str_i->occ[i]];
	  }

	  for(int i=0;i<nfzc;i++) {
	    int ii = fzc_occ[i];
	    for(int j=0;j<nfzc;j++)
	      CSC[i][j] = CSC_full[ii][fzc_occ[j]];
	  }

	  // Compute the determinant
	  FLOAT sign;
	  lu_decom(CSC, ndocc, tmpintvec, &sign);
	  FLOAT deter1 = 1.0;
	  for(int i=0;i<ndocc;i++)
	    deter1 *= CSC[i][i];
	  
	  S_a[jp-jp_begin][im-im_begin] = sign*deter1;

	} // end of im loop
      } // end of jp loop

      // Write out the overlap block
      ovlp_a.write(bp,bm);

    } // end of bm loop
  } // end of bp loop

  int ndets = vecm->size;

  //
  // Sort list of determinants to increasing (Ia,Ib) order
  //
  // create vector of determinants to sort
  std::vector<DetI> dets(ndets);
  for(int I=0; I<ndets; ++I) {
    dets[I].first = I;
    SlaterDet *detI = vecm->sdset->dets + I;
    int Istra = detI->alphastring;
    int Istrb = detI->betastring;
    dets[I].second.Ia = strblks_a.block(Istra);
    dets[I].second.Ib = strblks_a.block(Istrb);
  }
  // sort
  std::sort(dets.begin(),dets.end(),detcomp);
  // print
  if (Params.print_lvl >= PrintLevels::print_everything)
    for(int I=0; I<ndets; ++I) {
      std::cout << "Det #" << I << ": "
                << dets[I].first << " ("
                << dets[I].second.Ia << ","
                << dets[I].second.Ib << ")" << std::endl;
    }

  //
  // Evaluate total overlap in the highest available precision
  //
  FLOAT S_tot = 0.0;

#if USE_CI_OVERLAP

  StringBlockedMatrix ovlp_b(ovlp_a);
  CIOverlap ciovlp(vecp,vecm,ovlp_a,ovlp_b,Params.num_threads);
  ciovlp.compute();
  S_tot = ciovlp.value();

#else

  // string overlaps will be held in these buffers
  // ovlp_a has been defined before
  StringBlocks& strblks_b = strblks_a;
  StringBlockedMatrix ovlp_b(&strblks_b,&strblks_b,"S_a");
  const int nblks_b = strblks_b.nblocks();
  FLOAT** S_b;
  // HACK: if nblocks == 1 then the overlap is held in memory
  if (strblks_a.nblocks() == 1)
    S_b = S_a;
  else
    S_b = ovlp_b.buffer();

#if LOOP_OVER_BLOCKS
  // loop over string overlap blocks and then loop over all determinant pairs whose overlap
  // can be computed using these string overlap blocks
  for(int blk_ma=0; blk_ma<nblks_a; ++blk_ma) {
    for(int blk_pa=0; blk_pa<nblks_a; ++blk_pa) {

      ovlp_a.read(blk_pa, blk_ma);

      for(int blk_mb=0; blk_mb<nblks_b; ++blk_mb) {
	for(int blk_pb=0; blk_pb<nblks_b; ++blk_pb) {
	  
	  ovlp_b.read(blk_pb, blk_mb);

	  for(int I=0; I<ndets; I++) {
	    // map determinant back to the order in which CI vector is expressed
	    const int II = dets[I].first;
	    SlaterDet *detI = vecm->sdset->dets + II;
	    int Istra = detI->alphastring;
	    int Istrb = detI->betastring;

	    // skip this determinant if these string overlap blocks don't contribute
	    if (blk_ma != strblks_a.block(Istra) ||
		blk_mb != strblks_a.block(Istrb))
	      continue;

	    FLOAT cI = vecm->coeffs[II];
	    const int istra = strblks_a.rel_to_block_begin(Istra);
	    const int istrb = strblks_a.rel_to_block_begin(Istrb);

	    for(int J=0; J<ndets; J++) {
	      const int JJ = dets[J].first;
	      SlaterDet *detJ = vecp->sdset->dets + JJ;
	      int Jstra = detJ->alphastring;
	      int Jstrb = detJ->betastring;

	      // skip this determinant if these string overlap blocks don't contribute
	      if (blk_pa != strblks_b.block(Jstra) ||
		  blk_pb != strblks_b.block(Jstrb))
		continue;

	      FLOAT cJ = vecp->coeffs[JJ];
	      const int jstra = strblks_a.rel_to_block_begin(Jstra);
	      const int jstrb = strblks_a.rel_to_block_begin(Jstrb);
	      
	      FLOAT S = S_a[jstra][istra] * S_b[jstrb][istrb];
	      
	      FLOAT contrib = cJ * S * cI;
	      S_tot += contrib;
	    }
	  }

	}
      }
    }
  }
	      
#else

  int blk_ma = 0, blk_mb = 0;     // overlap between these string blocks is held in memory
  for(int I=0; I<ndets; I++) {
    // map determinant back to the order in which CI vector is expressed
    const int II = dets[I].first;
    SlaterDet *detI = vecm->sdset->dets + II;
    int Istra = detI->alphastring;
    int Istrb = detI->betastring;
    FLOAT cI = vecm->coeffs[II];
    const int istra = strblks_a.rel_to_block_begin(Istra);
    const int istrb = strblks_a.rel_to_block_begin(Istrb);
    blk_ma = strblks_a.block(Istra);
    blk_mb = strblks_a.block(Istrb);

    int blk_pa = 0, blk_pb = 0;     // overlap between these string blocks is held in memory
    for(int J=0; J<ndets; J++) {
      const int JJ = dets[J].first;
      SlaterDet *detJ = vecp->sdset->dets + JJ;
      int Jstra = detJ->alphastring;
      int Jstrb = detJ->betastring;
      FLOAT cJ = vecp->coeffs[JJ];
      const int jstra = strblks_a.rel_to_block_begin(Jstra);
      const int jstrb = strblks_a.rel_to_block_begin(Jstrb);
      blk_pa = strblks_a.block(Jstra);
      blk_pb = strblks_a.block(Jstrb);

      std::cout << "dets = (" << I << "," << J << ") : ("
		<< blk_pa << "," << blk_ma << ")  ("
		<< blk_pb << "," << blk_mb << ")" << std::endl;
      ovlp_a.read(blk_pa, blk_ma);
      ovlp_b.read(blk_pb, blk_mb);

      FLOAT S = S_a[jstra][istra] * S_b[jstrb][istrb];

      FLOAT contrib = cJ * S * cI;
      S_tot += contrib;
      /* fprintf(outfile,"  %3d %3d %+15.10Le", I, J, cI);
      fprintf(outfile," %+15.10Le", cJ);
      fprintf(outfile," %+25.15Le", S);
      fprintf(outfile," %+25.15Le", contrib);
      fprintf(outfile," %+25.15Le\n", S_tot); */

    }
  }

#endif // if LOOP_OVER_BLOCKS
#endif // if USE_CI_OVERLAP

  // Cleanup
  slaterdetvector_delete_full(vecm);
  slaterdetvector_delete_full(vecp);
  delete[] tmpintvec;
  delete_matrix(CSC);
  delete_matrix(CSC_full);
  double S_tot_double = (double) S_tot;
  return fabs(S_tot_double);
}

}} // namespace psi::dboc