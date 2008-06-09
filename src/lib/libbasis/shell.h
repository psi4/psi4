/*! \file
    \ingroup BASIS
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_lib_libbasis_shell_h_
#define _psi_src_lib_libbasis_shell_h_

#include <psitypes.h>

namespace psi {

class BasisSet;

class GaussianShell {

  friend class BasisSet;
  
  PSI_FLOAT O_[3];
  PSI_FLOAT *exps_;
  PSI_FLOAT **ccoeffs_;   // contractions in columns
  int *am_;
  int max_am_;
  int min_am_;
  int num_bf_;
  int num_ao_;
  int num_prims_;
  int num_contr_;
  bool puream_;

  void check_contraction_index_(int ci) const;
  void check_primitive_index_(int pi) const;

  // no default constructor
  GaussianShell();
  // no assignment operator
  GaussianShell& operator=(const GaussianShell&);

  /// Set the origin
  void set_origin(PSI_FLOAT[3]);

 public:
  GaussianShell(int nprims, int ncontr, int *am, bool puream, PSI_FLOAT *exps, PSI_FLOAT **ccoeffs, PSI_FLOAT origin[3]);
  GaussianShell(const GaussianShell&);
  ~GaussianShell();

  /// Print out GaussianShell
  void print(int id, FILE* outfile) const;

  /// Return the number of Cartesian Gaussians in this shell
  int num_ao() const { return num_ao_; };
  /// Return the number of basis functions in this shell
  int num_bf() const { return num_bf_; };
  /// Return the number of primitive Gaussians in this shell
  int num_prims() const { return num_prims_; };
  /** Return the angular momentum of contractions (throw an exception if contractions of several
      angular momenta are in this shell) */
  int am() const;
  /// Return the angular momentum of ci'th contraction
  int am(int ci) const;
  /// Return the maximum angular momentum of any contraction in this shell
  int max_am() const;
  /// Return the minimum angular momentum of any contraction in this shell
  int min_am() const;
  /// Use spherical harmonics Gaussians?
  bool puream() const { return puream_; };
  /// Return i'th Cartesian coordinate of the origin
  PSI_FLOAT origin(int i) const { return O_[i]; };
  /// Return orbital exponent of i'th primitive
  PSI_FLOAT exp(int pi) const;
  /// Return coefficient of pi'th primitive in ci'th contraction
  PSI_FLOAT cc(int ci, int pi) const;

};

} // end of namespace psi

#endif
