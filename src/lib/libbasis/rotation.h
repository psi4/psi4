/*! \file 
    \ingroup (BASIS)
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_lib_libbasis_rotation_h_
#define _psi_src_lib_libbasis_rotation_h_

#include "basisset.h"
#include "gnorm.h"
#include "combinate.h"

class RotationOp {

  BasisSet* bs_;
  GaussianNormalization gnorm_;
  StatCombData scdata_;

  int maxam_;     // maximum angular momentum this RotationOp can handle
  PSI_FLOAT*** Rl_;  // Transformation matrices for cartesian shells of each angular momentum

  void check_am(int l) const;
  void init_Rl(PSI_FLOAT** R);

  // no default constructor
  RotationOp();
  // no copy constructor
  RotationOp(const RotationOp&);
  // no assignment operator
  RotationOp& operator=(const RotationOp&);

 public:
  RotationOp(BasisSet* bs);
  ~RotationOp();

  /// Returns the full transformation matrix for the current cartesian Gaussian basis
  PSI_FLOAT** full_rotation_mat(PSI_FLOAT** R);
  /// Return the transformation matrix for a shell of a given angular momentum
  PSI_FLOAT** rotation_mat(PSI_FLOAT** R, int l);
};

#endif
