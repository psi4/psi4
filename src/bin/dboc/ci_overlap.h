/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/

#ifndef _psi3_bin_DBOC_cioverlap_h_
#define _psi3_bin_DBOC_cioverlap_h_

#include <libqt/slaterdset.h>
#include "stringblocks.h"

namespace psi{ namespace DBOC {

/// Computes overlap between 2 CI vectors. Uses precomputed string overlap matrix for alpha and beta spins.
class CIOverlap {
 public:
  CIOverlap(SlaterDetVector* vecbra, SlaterDetVector* vecket,
	    StringBlockedMatrix& ovlp_a, StringBlockedMatrix& ovlp_b,
	    unsigned int nthreads);
  ~CIOverlap();

  SlaterDetVector* vecbra() const { return vecbra_; }
  SlaterDetVector* vecket() const { return vecket_; }
  StringBlockedMatrix& ovlp_a() const { return ovlp_a_; }
  StringBlockedMatrix& ovlp_b() const { return ovlp_b_; }
  unsigned int nthreads() const { return nthreads_; }

  void compute();
  /// Returns the value of the overlap
  double value() const;
  // thread body -- computes overlaps for all determinants whose strings belong to blocks specified in threadgrp_
  void thread_compute(int tid);

 private:
  SlaterDetVector* vecbra_;
  SlaterDetVector* vecket_;
  StringBlockedMatrix& ovlp_a_;
  StringBlockedMatrix& ovlp_b_;
  unsigned int nthreads_;
  bool evaluated_;
  FLOAT S_;

  class ThreadGrp {
   public:
    ThreadGrp(int nthreads) :
      blkbra_a(0), blkbra_b(0), blkket_a(0), blkket_b(0),
      Sthr(new FLOAT[nthreads]) { for(int t=0; t<nthreads; ++t) Sthr[t] = 0.0; }
    ~ThreadGrp() { delete[] Sthr; }
    void set_blocks(int bb_a, int bb_b, int bk_a, int bk_b) {
      blkbra_a = bb_a;
      blkbra_b = bb_b;
      blkket_a = bk_a;
      blkket_b = bk_b;
    }
    /// quartet of string blocks
    int blkbra_a, blkbra_b;
    int blkket_a, blkket_b;
    // thread contributions
    FLOAT* Sthr;
  };
  ThreadGrp threadgrp_;
};

}} /* namespace psi::DBOC */

#endif
