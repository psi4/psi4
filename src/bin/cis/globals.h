/*! \file
    \ingroup CIS
    \brief Enter brief description of file here
*/
#include <ccfiles.h>

namespace psi { namespace cis {

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
EXTERN struct Local local;

enum Spin {singlet, triplet, uhf};

}} // namespace psi::cis
