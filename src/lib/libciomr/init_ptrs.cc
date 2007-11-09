#include <psifiles.h>
#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

extern "C" {

/* allocates memory for array of file pointers */

void init_ptrs(void)
{
  int num_ptrs = MAX_UNIT;

  ptr.wptr = (PSI_FPTR *) malloc(sizeof(PSI_FPTR)*num_ptrs);

  if (ptr.wptr == NULL) {
    fprintf(stderr,"trouble allocating memory for pointers!\n");
    exit(PSI_RETURN_FAILURE);
  }
}

void free_ptrs(void)
{
  free(ptr.wptr);
}

} /* extern "C" */
