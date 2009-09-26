/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here 
*/
#include <string>

namespace psi { namespace cchbar {

/* Input parameters for cchbar */
struct Params {
  long int memory;
  int cachelev;
  int ref;
  int print;
  std::string wfn;
  int dertype;
  int Tamplitude;
  int wabei_lowdisk;
};

}} // namespace psi::cchbar
