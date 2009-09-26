/*! \file
    \ingroup CIS
    \brief Enter brief description of file here 
*/
#include <string> 

namespace psi { namespace cis {

/* Input parameters */
struct Params {
  long int memory;
  std::string wfn;
  std::string diag_method;
  double convergence;
  int maxiter;
  int ref;
  int cis_ref;
  int print;
  int *rpi;
  int local;
};


}} // namespace psi::cis
