/*! \file
    \ingroup CCTRIPLES
    \brief Enter brief description of file here 
*/
#include <string>

namespace psi { namespace cctriples {

struct Params {
  int ref;
  std::string wfn;
  int semicanonical;
  int nthreads;
  int dertype;
};

}} // namespace psi::CCTRIPLES
