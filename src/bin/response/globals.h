/*! \file
    \ingroup RESPONSE
    \brief Enter brief description of file here 
*/
#include <ccfiles.h>

namespace psi { namespace response {

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN int *ioff;
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))

}} // namespace psi::response
