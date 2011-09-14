/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/

#include <string>

namespace psi { namespace ccsort {

struct Params {
  int print_lvl;         /* Output level control */
  int keep_TEIFile;      /* Should we keep the input two-elec. integrals? */
  int keep_OEIFile;      /* Should we keep the input one-elec. integrals? */
  double tolerance;      /* Cutoff value for integrals in IWL Buffers */
  long int memory;       /* Memory available (in bytes) */
  int cachelev;
  int ref;
  double *omega;         /* energy of applied field (a.u) for dynamic properties */
  int nomega;            /* number of field energies desired */
  std::string wfn;
  int make_abcd;
  int make_unpacked_abcd;
  int make_aibc;
  int dertype;
  std::string aobasis;
  int reset;            /* cmdline argument; if true, all CC-related files */
                        /* are deleted at the beginning of the run */
  int semicanonical;    /* semicanonical orbitals for perturbation theory */
  int local;
  std::string prop;
};

}} // namespace psi::ccsort
