#include <stdio.h>
#include "iomrparam.h"

extern "C" {

int
io_getline(FILE* input, char line[MAX_STRING])
{
  int i;

  for (i=0; i<MAX_STRING; i++) {
    if ((line[i]=fgetc(input)) == EOF) return(-1);
    if(feof(input)) return(-1);
    if (line[i] == '\n') {
      line[i] = '\0';
      return(0);
      }
    }

  fprintf(stdout,"io_getline: buffer size exceeded\n");
  fprintf(stderr,"io_getline: buffer size exceeded\n");
  return(-1);
  }

} /* extern "C" */
