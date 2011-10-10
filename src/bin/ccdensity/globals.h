/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <ccfiles.h>
#include <libdpd/dpd.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

namespace psi {
extern FILE* outfile;
namespace ccdensity {

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

/* #define DEBUG_XI (1)*/

EXTERN struct MOInfo moinfo;
EXTERN struct Frozen frozen;
EXTERN struct Params params;
EXTERN struct RHO_Params *rho_params;
EXTERN struct TD_Params *td_params;

}} // namespace psi::ccdensity
