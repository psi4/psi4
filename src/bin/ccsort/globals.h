/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here
*/
#include <ccfiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

namespace psi {
extern FILE* outfile;
namespace ccsort {

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

}} // namespace psi::ccsort
