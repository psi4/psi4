#include "iomrparam.h"
#include "includes.h"
#include <libipv1/ip_lib.h>

extern "C" {

int get_file_info(char* token, char* format,void* val)
{
  int i,errcod;
  char *prog,unit[4],*junk;
  char ip_token[MAX_STRING];
  char *gprgid();

  prog = gprgid();
  junk = strchr(token,':');
  junk++;
  i=0;
  while(*junk != ':') {
    unit[i] = *junk;
    junk++;
    i++;
    }
  junk++;
  unit[i]='\0';

  sprintf(ip_token,":%s:FILES:FILE%s:%s",prog,unit,junk);
  errcod = ip_data(ip_token,format,val,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":%s:FILES:DEFAULT:%s",prog,junk);
  errcod = ip_data(ip_token,format,val,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":DEFAULT:FILES:FILE%s:%s",unit,junk);
  errcod = ip_data(ip_token,format,val,0);
  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:%s",junk);
  errcod = ip_data(ip_token,format,val,0);
  if(errcod == IPE_OK) return(0);

  return(-1);
  }

} /* extern "C" */
