/*!
** \file
** \brief Get the amount of core memory available from input
** \ingroup CIOMR
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libipv1/ip_lib.h>
#include <psi4-dec.h>

namespace psi {

#define DEF_MAXCRR 32000000  /* default maxcor in doubles: used only if
                              * it can't be read in */

void set_memory(FILE *infile, FILE *outfile)
{
  char type[20];
  char *s;
  int count;
  long int maxcrr;            /* maxcor in real words */
  char *maxcrr_str;           /* string representation of maxcrr */
  double size;
  int errcod; 

  maxcrr = DEF_MAXCRR;        /* set maxcor to default first */

  if(ip_exist(const_cast<char*>("MEMORY"),0)) { /* check if the keyword exists */
    errcod = ip_count(const_cast<char*>("MEMORY"), &count, 0);

    if (errcod != IPE_OK) throw("Cannot read memory");

    else if (errcod == IPE_NOT_AN_ARRAY) { /* Scalar specification of MEMORY */
      errcod = ip_string(const_cast<char*>("MEMORY"), &maxcrr_str, 0);
      if (errcod != IPE_OK) throw("Cannot read memory");
      maxcrr = atol(maxcrr_str);
    }
    /* Array specification of MEMORY */
    else if (count == 1) {
      errcod = ip_string(const_cast<char*>("MEMORY"), &maxcrr_str, 0);
      if (errcod != IPE_OK) throw("Cannot read memory");
      maxcrr = atol(maxcrr_str);
    }
    else if (count == 2) {
      errcod = ip_data(const_cast<char*>("MEMORY"), const_cast<char*>("%lf"), &size, 1, 0);
      if (errcod != IPE_OK) throw("Cannot read memory");
      errcod = ip_data(const_cast<char*>("MEMORY"), const_cast<char*>("%s"), type, 1, 1);
      if (errcod != IPE_OK) throw("Cannot read memory");
      /* convert string to uppercase */
      for (s=type; *s!='\0'; s++) {
	    if (*s>='a' && *s <='z') *s = *s + 'A' - 'a';
      }
      if ((strcmp(type, "R")==0) || (strcmp(type, "REAL")==0))
        maxcrr = (long int) size;
      else if ((strcmp(type, "I")==0) || (strcmp(type, "INTEGER")==0)) 
        maxcrr = (long int) (size * sizeof(int) / sizeof(double));
      else if ((strcmp(type, "B")==0) || (strcmp(type, "BYTES")==0)) 
        maxcrr = (long int) (size / sizeof(double));
      else if ((strcmp(type, "KB")==0) || (strcmp(type, "KBYTES")==0)) 
        maxcrr = (long int) (1000.0 * size / sizeof(double));
      else if ((strcmp(type, "MB")==0) || (strcmp(type, "MBYTES")==0))
        maxcrr = (long int) (1000000.0 * size / sizeof(double));
      else if ((strcmp(type, "GB")==0) || (strcmp(type, "GBYTES")==0))
        maxcrr = (long int) (1000000000.0 * size / sizeof(double));
      else {
   	    fprintf(outfile, "bad data type, specify one of: \n") ;
	    fprintf(outfile, "REAL, INTEGER, BYTES, KBYTES, MBYTES, or GBYTES\n");
	    throw("Cannot read memory");
      }
    }
  }
  
  module.set_memory(maxcrr * sizeof(double));

  return;
}

}


