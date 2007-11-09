/*! \file 
    \ingroup (BASIS)
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_lib_libbasis_overlap_h_
#define _psi_src_lib_libbasis_overlap_h_

#include <psitypes.h>
#include "basisset.h"
#include "shell.h"
#include "osrecur.h"
#include "gnorm.h"

class OverlapEngine {

  const BasisSet* bs1_;
  const BasisSet* bs2_;
  PSI_FLOAT* buffer_;
  OI_OSRecursor overlap_recur_;
  GaussianNormalization gnorm_;

  void compute_pair_(const GaussianShell&, const GaussianShell&);

 public:
  OverlapEngine(const BasisSet*, const BasisSet*);
  ~OverlapEngine();

  /// Compute overlap between 2 shells. Result is stored in buffer
  void compute_shell_pair(int, int);
  /// Return buffer which holds the result of compute_shell_pair()
  PSI_FLOAT* buffer() const { return buffer_; };
  /// Allocate and compute full overlap matrix
  PSI_FLOAT** compute_full_matrix();
};

#endif
