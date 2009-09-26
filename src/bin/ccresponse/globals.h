/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here 
*/
#include <ccfiles.h>
#include <psi4-dec.h>

namespace psi { namespace CCRESPONSE {

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
EXTERN struct Local local;

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))

}} // namespace psi::CCRESPONSE
