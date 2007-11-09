#include "includes.h"
#include "iomrparam.h"
#include "types.h"

extern "C" {

s_async_t *s_async_ioopen(char *param, int unit)
{
  fprintf(stderr,"no s_async io yet\n");
  ioabort();
  }

void s_async_ioclos(s_async_t *ud, int status)
{
  fprintf(stderr,"no s_async io yet\n");
  ioabort();
  }

void s_async_iordr(s_async_t *ud, char *buffer, PSI_FPTR first, int length)
{
  fprintf(stderr,"no s_async io yet\n");
  ioabort();
  }

void s_async_iowrr(s_async_t *ud, char *buffer, PSI_FPTR first, int length)
{
  fprintf(stderr,"no s_async io yet\n");
  ioabort();
  }

PSI_FPTR s_async_iosize(s_async_t *ud)
{
  fprintf(stderr,"no s_async io yet\n");
  ioabort();
  }

} /* extern "C" */
