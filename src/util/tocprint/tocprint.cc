/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <psi4-def.h>

namespace psi{
  namespace tocprint{

  /* Function prototypes */
  void init_io(int argc, char *argv[]);
  void exit_io(void);
  int cc2unit(char *);
  
  } /* Namespace tocprint */
} /* Namespace psi */

using namespace psi::tocprint;
using namespace psi;

int main(int argc, char *argv[])
{
  int unit=0, i=0, get_chkpt_prefix=0;
  char *prefix;

  while(--argc > 0) {
       i++;
       if(!strcmp(argv[i], "-unit") || !strcmp(argv[i], "-u")) {
           unit = atoi(argv[i+1]);
           argc--;
         }
       else if(!strcmp(argv[i], "-cc")) {
           unit = cc2unit(argv[i+1]);
           argc--;
         }
       else if(!strcmp(argv[i], "-prefix")) {
           unit = PSIF_CHKPT;
           get_chkpt_prefix = 1;
           argc--;
         }
    }

  if(!unit) { printf("Bad unit number.\n"); exit(PSI_RETURN_FAILURE); }

  init_io(argc-1, argv+i);

  psio_open(unit,PSIO_OPEN_OLD);
  if (!get_chkpt_prefix)
    psio_tocprint(unit,stdout);
  else {
    chkpt_init(PSIO_OPEN_OLD);
    prefix = chkpt_rd_prefix();
    chkpt_close();
    printf("Checkpoint file prefix: %s\n", prefix);
    free(prefix);
  }
  psio_close(unit,1);

  exit_io();
  exit(PSI_RETURN_SUCCESS);
}

extern "C"{ char *gprgid() { char *prgid = "TOCPRINT"; return(prgid); } }

namespace psi {
  namespace tocprint {

  void init_io(int argc, char *argv[])
  {
   char *progid;
  
   progid = (char *) malloc(strlen(gprgid())+2);
   sprintf(progid, ":%s",gprgid());
  
   psi_start(&infile,&outfile,&psi_file_prefix,argc-1, argv+1, 0);
   ip_cwk_add(progid);
   free(progid);
  
   psio_init(); psio_ipv1_config();
  }

  void exit_io(void)
  {
    psio_done();
    psi_stop(infile,outfile,psi_file_prefix);
  }

  } /* Namespace tocprint */
} /* Namespace psi */
