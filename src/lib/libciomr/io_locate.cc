#include "includes.h"
#include "iomrparam.h"
#include <libciomr/libciomr.h>

extern "C" {

int
io_locate(FILE* input, char loc_token[MAX_STRING])
{
  char line[MAX_STRING];

  fseek(input,0L,0);

  for (;;) {
    if (io_getline(input,line) != 0) return(-1);
    if (!strncmp(loc_token,line,10)) {
      return(0);
      }
    }
  }

} /* extern "C" */
