#define MAIN

#include <psifiles.h>
#include "iomrparam.h"
#include "includes.h"
#include "types.h"
#include <libciomr/libciomr.h>

extern "C" {

ioFILE_t ud[MAX_UNIT];
unsigned int current_unit;

void
ioinit_()
{
  int i;

  for (i=0; i<MAX_UNIT; i++) {
    ud[i].method = PSIIO_UNOPENED;
    ud[i].ptr.sequential = NULL;
    }
  }

void
ioopen_(int* unit)
{
  char param[MAX_STRING];
  char methodparam[MAX_STRING];
  char method[MAX_STRING]; 
  char unitch[MAX_STRING];

  current_unit = *unit;

  sprintf(unitch,"%d",*unit);

  strcpy(param,"FILES:");
  strcat(param,unitch);
  strcat(param,":");

  strcpy(methodparam,param);

  if(oldstyleinput()) {
    strcat(methodparam,"method");
    if (get_param(methodparam,"%s",method) != 0) {
      strcpy(method,"sequential");
      }
    }
  else {
    strcat(methodparam,"METHOD");
    if (get_file_info(methodparam,"%s",method) != 0) {
      strcpy(method,"sequential");
      }
    }

  if (!strcmp(method,"sequential")) {
    ud[*unit].method = PSIIO_SEQUENTIAL;
    ud[*unit].ptr.sequential = sequential_ioopen(param,*unit);
    }
  else if (!strcmp(method,"r_async")) {
    ud[*unit].method = PSIIO_R_ASYNC;
    ud[*unit].ptr.r_async = r_async_ioopen(param,*unit);
    }
  else if (!strcmp(method,"s_async")) {
    ud[*unit].method = PSIIO_S_ASYNC;
    ud[*unit].ptr.s_async = s_async_ioopen(param,*unit);
    }
  else if (!strcmp(method,"ram")) {
    ud[*unit].method = PSIIO_RAM;
    ud[*unit].ptr.ram = ram_ioopen(param,*unit);
    }
  else {
    fprintf(stderr,"ioopen_: invalid method (%s) for unit %d\n",method,*unit);
    ioabort();
    }
  }

void
ioclos_(int* unit,int* status)
{

  current_unit = *unit;

  if (ud[*unit].method == PSIIO_SEQUENTIAL) {
    ud[*unit].method = PSIIO_UNOPENED;
    sequential_ioclos(ud[*unit].ptr.sequential, *status);
    }
  else if (ud[*unit].method == PSIIO_R_ASYNC) {
    ud[*unit].method = PSIIO_UNOPENED;
    r_async_ioclos(ud[*unit].ptr.r_async, *status);
    }
  else if (ud[*unit].method == PSIIO_S_ASYNC) {
    ud[*unit].method = PSIIO_UNOPENED;
    s_async_ioclos(ud[*unit].ptr.s_async, *status);
    }
  else if (ud[*unit].method == PSIIO_RAM) {
    ud[*unit].method = PSIIO_UNOPENED;
    ram_ioclos(ud[*unit].ptr.ram, *status);
    }
  else {
    fprintf(stderr,"ioclos_: invalid method (%d) for unit %d\n",
             ud[*unit].method,*unit);
    ioabort();
    }
  free(ud[*unit].ptr.sequential);
  ud[*unit].ptr.sequential = NULL;
  ud[*unit].method = PSIIO_UNOPENED;
  }

void
iowrr_(int* unit,char* buffer,PSI_FPTR* first,int* length)
{

  current_unit = *unit;

  if (*first + *length < *first) {
    fprintf(stderr,"iowrr_: *unit=%d, *first=%lu, *length=%lu\n",
            *unit, *first, *length);
    ioabort();
    }

  if (ud[*unit].method == PSIIO_SEQUENTIAL) {
    sequential_iowrr(ud[*unit].ptr.sequential,buffer,*first,*length);
    }
  else if (ud[*unit].method == PSIIO_R_ASYNC) {
    r_async_iowrr(ud[*unit].ptr.r_async,buffer,*first,*length);
    }
  else if (ud[*unit].method == PSIIO_S_ASYNC) {
    s_async_iowrr(ud[*unit].ptr.s_async,buffer,*first,*length);
    }
  else if (ud[*unit].method == PSIIO_RAM) {
    ram_iowrr(ud[*unit].ptr.ram,buffer,*first,*length);
    }
  else {
    fprintf(stderr,"iowrr_: invalid method (%d) for unit %d\n",
            ud[*unit].method,*unit);
    ioabort();
    }
  }

void
iordr_(int* unit, char* buffer,PSI_FPTR* first,int* length)
{

  current_unit = *unit;

  if (*first + *length < *first) {
    fprintf(stderr,"iordr_: *unit=%d, *first=%lu, *length=%lu\n",
            *unit, *first, *length);
    ioabort();
    }

  if (ud[*unit].method == PSIIO_SEQUENTIAL) {
    sequential_iordr(ud[*unit].ptr.sequential,buffer,*first,*length);
    }
  else if (ud[*unit].method == PSIIO_R_ASYNC) {
    r_async_iordr(ud[*unit].ptr.r_async,buffer,*first,*length);
    }
  else if (ud[*unit].method == PSIIO_S_ASYNC) {
    s_async_iordr(ud[*unit].ptr.s_async,buffer,*first,*length);
    }
  else if (ud[*unit].method == PSIIO_RAM) {
    ram_iordr(ud[*unit].ptr.ram,buffer,*first,*length);
    }
  else {
    fprintf(stderr,"iordr_: invalid method (%d) for unit %d\n",
            ud[*unit].method,*unit);
    ioabort();
    }

  }

void
ioabort()
{
  fprintf(stderr,"ioabort: current unit = %d\n",current_unit);
  exit(PSI_RETURN_FAILURE);
  }

PSI_FPTR
iosize_(int* unit)
{

   PSI_FPTR fsize=0;

   current_unit = *unit;

   if (ud[*unit].method == PSIIO_UNOPENED) {
      fprintf(stderr,"iosize_: Unit %d not opened\n", current_unit);
      return(0);
      }

   else if (ud[*unit].method == PSIIO_SEQUENTIAL) {
      fsize = sequential_iosize(ud[*unit].ptr.sequential); 
      }

   else if (ud[*unit].method == PSIIO_R_ASYNC) {
      fsize = r_async_iosize(ud[*unit].ptr.r_async); 
      }

   else if (ud[*unit].method == PSIIO_S_ASYNC) {
      fsize = s_async_iosize(ud[*unit].ptr.s_async); 
      }

   else if (ud[*unit].method == PSIIO_RAM) {
      fsize = ram_iosize(ud[*unit].ptr.ram); 
      }

   else {
      fprintf(stderr,"iosize_: Invalid method for unit %d\n", current_unit);
      return(0);
      }

   return(fsize);
}
      
} /* extern "C" */
