/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <ccfiles.h>
#include <libdpd/dpd.h>

namespace psi { namespace ccdensity {

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

/* #define DEBUG_XI (1)*/

extern "C" {
EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;
}
EXTERN struct MOInfo moinfo;
EXTERN struct Frozen frozen;
EXTERN struct Params params;
EXTERN struct RHO_Params *rho_params;
EXTERN struct TD_Params *td_params;

}} // namespace psi::ccdensity
