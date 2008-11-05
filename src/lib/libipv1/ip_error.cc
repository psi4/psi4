/*! \file
    \ingroup IPV1
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include "tmpl.h"
#include "ip_types.h"
#include "ip_global.h"
#include <psifiles.h>

/* Cannot include ip_error.global due to xlc's handling of varargs. */
#include "ip_error.gbl"
#include "ip_error.lcl"
#include "ip_error.h"

#include "scan.gbl"

extern "C" {

/* Returns some text for an errcod. */
char *ip_error_message(int errcod)
{
  static char *ipe_ok = "No problem has been detected.";
  static char *ipe_key_not_found = "No match was found for the given keyword.";
  static char *ipe_out_of_bounds = "An array index is out of bounds.";
  static char *ipe_malloc = "Memory allocation failed.";
  static char *ipe_not_an_array = "An index was given for a scalar quantity.";
  static char *ipe_not_a_scalar = "Expected a scalar, but found an array.";
  static char *ipe_type = "The datum is not of the appropiate type.";
  static char *huh = "The nature of the problem is unknown.";

  if (errcod == IPE_OK) return ipe_ok;
  if (errcod == IPE_KEY_NOT_FOUND) return ipe_key_not_found;
  if (errcod == IPE_OUT_OF_BOUNDS) return ipe_out_of_bounds;
  if (errcod == IPE_MALLOC) return ipe_malloc;
  if (errcod == IPE_NOT_AN_ARRAY) return ipe_not_an_array;
  if (errcod == IPE_NOT_A_SCALAR) return ipe_not_a_scalar;
  if (errcod == IPE_TYPE) return ipe_type;
  return huh;
  }

void ip_error(char *msg, ...)
{
  va_list args;
  va_start(args,msg);
  fprintf(ip_out,"IP_ERROR: ");
  vfprintf(ip_out,msg,args);
  fprintf(ip_out,"\n");
  va_end(args);
  showpos();
  exit(PSI_RETURN_FAILURE);
  }

void ip_warn(char *msg, ...)
{
  va_list args;
  va_start(args,msg);
  fprintf(ip_out,"IP_WARN: ");
  vfprintf(ip_out,msg,args);
  fprintf(ip_out,"\n");
  va_end(args);
  }

} /* extern "C" */
