#include "includes.h"
#include "iomrparam.h"
#include "types.h"

extern "C" {

r_async_t * r_async_ioopen(char *param, int unit)
{
  fprintf(stderr,"no r_async io yet\n");
  ioabort();
  }

void r_async_ioclos(r_async_t *ud, int status)
{
  fprintf(stderr,"no r_async io yet\n");
  ioabort();
  }

void r_async_iordr(r_async_t *ud, char *buffer, PSI_FPTR first, int length)
{
  fprintf(stderr,"no r_async io yet\n");
  ioabort();
  }

void r_async_iowrr(r_async_t *ud, char *buffer, PSI_FPTR first, int length)
{
  fprintf(stderr,"no r_async io yet\n");
  ioabort();
  }

PSI_FPTR r_async_iosize(r_async_t *ud)
{
  fprintf(stderr,"no r_async io yet\n");
  ioabort();
  }

} /* extern "C" */
