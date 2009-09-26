/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here
*/

#ifndef _psi3_bin_DBOC_params_h_
#define _psi3_bin_DBOC_params_h_

namespace psi { namespace DBOC {

namespace PrintLevels {
  static const int print_intro = 1;
  static const int print_params = 1;
  static const int print_contrib = 2;
  static const int print_everything = 5;
}

typedef struct {

  enum RefType { rhf=1, rohf=2, uhf=3};

  /// Cartesian coordinate structure
  typedef struct {
    int index;    // index of the coordinate
    int atom;     // which atom
    int xyz;      // x (=0), y (=1), or z (=2)
    double coeff; // Degeneracy (number of equivalent coords)
    bool symm;    // Whether plus displacement is equivalent to minus displacement
  } Coord_t;

  char *label;
  char *wfn;
  RefType reftype;
  double delta;
  int disp_per_coord;
  int ncoord;
  Coord_t* coords;
  int nisotope;
  char** isotopes;
  bool ref_frame_wfn; // compute wave functions in the reference frame
                      // (i.e. avoid reorientation and COM shift associated with
                      //  wfn computation that takes advantage of symmetry)

  unsigned int num_threads;  /// number of threads
  size_t max_memory;         /// maximum available memory, in bytes
  size_t memory;             /// currently available memory, in bytes
  int print_lvl;

} Params_t;

}} /* namespace psi::DBOC */

#endif
