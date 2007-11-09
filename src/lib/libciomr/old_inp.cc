#include "iomrparam.h"
#include "includes.h"
#include <libciomr/libciomr.h>

extern "C" {

int oldstyleinput(void)
{
  FILE *input;
  int ierr;

  input = fopen("input.dat","r");
  ierr = io_locate(input,"# FILES ##");
  fclose(input);
  if (ierr == 0) return(1);
  return 0;
  }

} /* extern "C" */
