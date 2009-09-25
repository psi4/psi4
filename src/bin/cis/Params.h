/*! \file
    \ingroup CIS
    \brief Enter brief description of file here 
*/

namespace psi { namespace cis {

/* Input parameters */
struct Params {
  long int memory;
  char *wfn;
  char *diag_method;
  double convergence;
  int maxiter;
  int ref;
  int cis_ref;
  int print;
  int *rpi;
  int local;
};


}} // namespace psi::cis
