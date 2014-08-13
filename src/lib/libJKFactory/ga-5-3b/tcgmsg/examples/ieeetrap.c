#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*$Id: ieeetrap.c,v 1.2 1995-02-02 23:24:12 d3g681 Exp $*/
#include <floatingpoint.h>
#include <stdio.h>
#include <signal.h>

static void catchit()
{
  printf("!!  Floating point interrupt caught  !!\n");
  fflush(stdout);
  (void) signal(SIGIOT, SIG_DFL);
  abort();
}

void ieeetrap_()
{
 (void) ieee_handler("set","inexact", SIGFPE_IGNORE);
 (void) ieee_handler("set","underflow", SIGFPE_IGNORE);
 (void) ieee_handler("set","invalid", SIGFPE_IGNORE);

}
