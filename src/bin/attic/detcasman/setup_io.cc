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

/*! 
** \file
** \ingroup DETCASMAN
** \brief Set up and shut down I/O for detcasman
**
** C. David Sherrill
** University of California, Berkeley
** April 1998
*/

#include <cstdio>
#include <cstdlib>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "globals.h"

namespace psi { namespace detcasman {

/*!
** init_io(): Function opens input and output files
**
** \param argc   = number of arguments from main
** \param argv   = argument list from main
**
** Returns: none
** \ingroup detcasman
*/
void init_io(int argc, char *argv[])
{
  int i;
  int num_extra_args=0;
  char **extra_args;

  extra_args = (char **) malloc(argc*sizeof(char *));
  
  for (i=1; i<argc; i++) {
    extra_args[num_extra_args++] = argv[i];
  }
  
  psi_start(&infile,&outfile,&psi_file_prefix,num_extra_args, extra_args, 0);
  ip_cwk_add(":DETCASMAN");
  ip_cwk_add(":DETCAS"); 
  tstart(outfile);
  free(extra_args);
}



/*!
** close_io(): Function closes down I/O and exits
** 
** Returns: none
** \ingroup detcasman
*/
void close_io(void)
{
  tstop(outfile);
  psi_stop(infile,outfile,psi_file_prefix);
}



/*!
** check() acts as an abort function if the condition 'a' is not true;
**   it shuts down i/o and returns an error
**
** \param a  = integer which is 0 (false) or nonzero (true)
** \param errmsg = error message to write if a is 0 (false)
**
** Returns: none
** \ingroup detcasman
*/
void check(int a, const char *errmsg)
{
  if (!a) {
    psi::fprintf(outfile, "%s\n", errmsg);
    close_io();
    exit(1);
  }
}

}} // end namespace psi::detcasman

/*
** The gprgid() function required by the PSI I/O libraries
*/
extern "C" {
  const char *gprgid()
  {
    const char *prgid = "DETCASMAN";
    return(prgid);
  }
}


