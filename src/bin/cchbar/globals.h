/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here 
*/
#include <psi4-dec.h>
#include <ccfiles.h>
#include <libdpd/dpd.h>
#include <libciomr/libciomr.>
#include <libpsio/psio.h>

namespace psi { namespace cchbar {

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
