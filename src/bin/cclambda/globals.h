/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <ccfiles.h>
#include <libdpd/dpd.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

namespace psi {
extern FILE* outfile;
namespace cclambda {

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

/* #define EOM_DEBUG (1) */

EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
EXTERN struct L_Params *pL_params;
EXTERN struct Local local;
void check_sum(char *lbl, int L_irr);

}} // namespace psi::cclambda
