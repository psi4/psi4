
#include "includes.h"
#include "iomrparam.h"
#include "types.h"

extern "C" {

ram_t * ram_ioopen(char *param, int unit)
{
  fprintf(stderr,"no ram io yet\n");
  ioabort();
  }

void ram_ioclos(ram_t *ud, int status)
{
  fprintf(stderr,"no ram io yet\n");
  ioabort();
  }

void ram_iordr(ram_t *ud, char *buffer, PSI_FPTR first, int length)
{
  fprintf(stderr,"no ram io yet\n");
  ioabort();
  }

void ram_iowrr(ram_t *ud, char *buffer, PSI_FPTR first, int length)
{
  fprintf(stderr,"no ram io yet\n");
  ioabort();
  }

PSI_FPTR ram_iosize(ram_t *ud)
{
  fprintf(stderr,"no ram io yet\n");
  ioabort();
  }

} /* extern "C" */
