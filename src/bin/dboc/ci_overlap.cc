/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/

#include <pthread.h>
#include "ci_overlap.h"

using namespace psi::DBOC;

namespace {
  // Packages ptr to CIOverlap and thread id
  typedef std::pair<CIOverlap*,int> objptr_id_t;
  void*
  thread_compute(void* objptr_id_voidptr)
  {
    objptr_id_t* objptr_id = static_cast<objptr_id_t*>(objptr_id_voidptr);
    int tid = objptr_id->second;
    CIOverlap* obj = objptr_id->first;
    // call
    obj->thread_compute(tid);
    return 0;
  }
}

CIOverlap::CIOverlap(SlaterDetVector* vecbra, SlaterDetVector* vecket,
		     StringBlockedMatrix& ovlp_a, StringBlockedMatrix& ovlp_b,
		     unsigned int nthreads) :
  vecbra_(vecbra), vecket_(vecket), ovlp_a_(ovlp_a), ovlp_b_(ovlp_b),
  nthreads_(nthreads), evaluated_(false), S_((FLOAT)0.0), threadgrp_(nthreads)
{
}

CIOverlap::~CIOverlap()
{
}

void
CIOverlap::compute()
{
  //
  // Loop over string block quartets, inside the loop spawn threads, accumulate the result
  //

  const int nblksbra_a = ovlp_a_.strblk_bra()->nblocks();
  const int nblksbra_b = ovlp_b_.strblk_bra()->nblocks();
  const int nblksket_a = ovlp_a_.strblk_ket()->nblocks();
  const int nblksket_b = ovlp_b_.strblk_ket()->nblocks();
  
  // loop over string overlap blocks and then loop over all determinant pairs whose overlap
  // can be computed using these string overlap blocks
  for(int bk_a=0; bk_a<nblksket_a; ++bk_a) {
    for(int bb_a=0; bb_a<nblksbra_a; ++bb_a) {

      ovlp_a_.read(bb_a, bk_a);

      for(int bk_b=0; bk_b<nblksket_b; ++bk_b) {
	for(int bb_b=0; bb_b<nblksbra_b; ++bb_b) {
	  
	  ovlp_b_.read(bb_b, bk_b);

	  // prepare to spawn threads
	  threadgrp_.set_blocks(bb_a, bb_b, bk_a, bk_b);

	  // spawn threads
	  objptr_id_t* objids = new objptr_id_t[nthreads_];
	  pthread_t* threads = new pthread_t[nthreads_];
	  pthread_attr_t thread_attr; pthread_attr_init(&thread_attr);
	  pthread_attr_setscope(&thread_attr,
				PTHREAD_SCOPE_SYSTEM);
	  for(int t=1; t<nthreads_; ++t) {
	    objids[t].first = this;
	    objids[t].second = t;
	    pthread_create(&(threads[t]),&thread_attr,
			   &::thread_compute,(void *)&(objids[t]) );
	  }
	  this->thread_compute(0);
	  for(int t=1;t<nthreads_;++t)
	    pthread_join(threads[t], NULL);
	  delete[] threads;
	  delete[] objids;

	}
      }
    }
  }

  // Accumulate the result
  for(int t=0; t<nthreads_; ++t)
    S_ += threadgrp_.Sthr[t];
}

FLOAT
CIOverlap::value() const
{
  return S_;
}

void
CIOverlap::thread_compute(int tid)
{
  const int bb_a = threadgrp_.blkbra_a;
  const int bb_b = threadgrp_.blkbra_b;
  const int bk_a = threadgrp_.blkket_a;
  const int bk_b = threadgrp_.blkket_b;
  const StringBlocks& strblks_bra_a(*(ovlp_a_.strblk_bra()));
  const StringBlocks& strblks_ket_a(*(ovlp_a_.strblk_ket()));
  const StringBlocks& strblks_bra_b(*(ovlp_b_.strblk_bra()));
  const StringBlocks& strblks_ket_b(*(ovlp_b_.strblk_ket()));
  FLOAT** Sa = ovlp_a_.buffer();
  FLOAT** Sb = ovlp_b_.buffer();

  FLOAT S_tot = 0.0;

  const int ndets = vecbra_->sdset->size;
  int IJ = 0;
  for(int I=0; I<ndets; I++) {
    SlaterDet *detI = vecket_->sdset->dets + I;
    int Istra = detI->alphastring;
    int Istrb = detI->betastring;

    // skip this determinant if these string overlap blocks don't contribute
    if (bk_a != strblks_ket_a.block(Istra) ||
	bk_b != strblks_ket_b.block(Istrb))
      continue;

    FLOAT cI = vecket_->coeffs[I];
    const int istra = strblks_ket_a.rel_to_block_begin(Istra);
    const int istrb = strblks_ket_b.rel_to_block_begin(Istrb);

    for(int J=0; J<ndets; ++J, ++IJ) {
      if (IJ%nthreads_ != tid)
	continue;
      SlaterDet *detJ = vecbra_->sdset->dets + J;
      int Jstra = detJ->alphastring;
      int Jstrb = detJ->betastring;

      // skip this determinant if these string overlap blocks don't contribute
      if (bb_a != strblks_bra_a.block(Jstra) ||
	  bb_b != strblks_bra_b.block(Jstrb))
	continue;

      FLOAT cJ = vecbra_->coeffs[J];
      const int jstra = strblks_bra_a.rel_to_block_begin(Jstra);
      const int jstrb = strblks_bra_b.rel_to_block_begin(Jstrb);
	      
      FLOAT S = Sa[jstra][istra] * Sb[jstrb][istrb];
	      
      FLOAT contrib = cJ * S * cI;
      S_tot += contrib;
    }
  }

  threadgrp_.Sthr[tid] += S_tot;
}
