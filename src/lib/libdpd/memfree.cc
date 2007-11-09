/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
/*
** Function to return number of double words available for allocation.
*/

#include <stdio.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

extern "C" {

long int dpd_memfree(void)
{
  return dpd_main.memory - (dpd_main.memused - 
			    dpd_main.memcache + 
			    dpd_main.memlocked);
}

void dpd_memset(long int memory)
{
  dpd_main.memory = memory;
}

} /* extern "C" */
