/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here
*/
#include <ccfiles.h>
#include <libdpd/dpd.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

namespace psi {
extern FILE* outfile;
namespace cchbar {

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN struct MOInfo moinfo;
EXTERN struct Params params;

}} // namespace psi::cchbar
