/*! \file
    \ingroup BASIS
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_lib_libbasis_basisset_h_
#define _psi_src_lib_libbasis_basisset_h_

#include <psitypes.h>
#include "shell.h"

namespace psi {

class BasisSet {

  int num_prims_;
  int num_shells_;
  int num_ao_;
  int num_bf_;
  int max_am_;
  bool puream_;

  GaussianShell** shells_;

  int *shell_fbf_;
  int *shell_fao_;
  int *shell_center_;
  
  // Cartesian coordinates of basis function centers  
  PSI_FLOAT** coords_;
  int ncenters_;

  void init_shells();
  /// Throw std::runtime_error if index out of bounds
  void check_shell_index(int si) const;

  // No default constructor
  BasisSet();
  // No assignment
  BasisSet& operator=(const BasisSet&);

 public:
  BasisSet(int chkptfile);
  BasisSet(const BasisSet&);
  ~BasisSet();

  /// Print out the basis set
  void print(char* id, FILE* outfile) const;

  int num_prims() const;
  int num_shells() const;
  int num_ao() const;
  int num_bf() const;
  int max_am() const;

  /// Return si'th gaussian shell
  GaussianShell& shell(int si) const;
  /// Return index of the first basis function from si'th shell
  int first_bf(int si) const;
  /// Return index of the first Cartesian Gaussian function from si'th shell
  int first_ao(int si) const;
  /// Return index of the center on which si'th shell is centered
  int center(int si) const;

  /// Set the coordinate of center ci to O
  void set_center(int ci, PSI_FLOAT[3]);
  /// Get i-th coordinate of center ci
  PSI_FLOAT get_center(int ci, int i);
};

} // end of namespace psi

#endif
