/*! \file
    \ingroup STABLE
    \brief Enter brief description of file here 
*/

namespace psi { namespace stable {

struct Params {
  int print_lvl;         /* Output level control */
  long int memory;       /* Memory available (in bytes) */
  int cachelev;
  int ref;
  int follow_instab;     /* follow a UHF->UHF instability of same symm? */
  int num_evecs_print;   /* print n lowest eigenvectors of MO hessian */
  int rotation_method;   /* 0 = by angles, 1 = by antisymmetric matrix */
  double scale;          /* scale factor for orbital rotation step */
};

}} // namespace psi::stable
