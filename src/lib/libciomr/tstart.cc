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
** \brief Controls starting and stopping of timers
** \ingroup CIOMR
*/

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <string>
#include <ctime>
#define EXTERN

#include <sys/times.h>

namespace psi {

extern FILE *outfile;

time_t time_start, time_end;
time_t time_start_overall;
int running = 0;
double user_start, sys_start;
double user_start_overall, sys_start_overall;
double user_stop, sys_stop;

/*!
** tstart(): Starts a timer
**
** \param outfile = output file pointer
**
** \ingroup CIOMR
*/
void tstart()
{
  int error;
  char *name;
  struct tms total_tmstime;
  const long clk_tck = sysconf(_SC_CLK_TCK);
  times(&total_tmstime);



  /// host name info, needed?
  name = (char *) malloc(40 * sizeof(char));
  error = gethostname(name, 40);
  if(error != 0) strncpy(name,"nohostname", 11);

  /// start a global timer
  if(!running){
     time_start_overall = time(NULL);
     user_start_overall = ((double) total_tmstime.tms_utime)/clk_tck;
     sys_start_overall = ((double) total_tmstime.tms_stime)/clk_tck;
     running = 1;
  }

  /// start module timers
  time_start = time(NULL);
  user_start = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_start = ((double) total_tmstime.tms_stime)/clk_tck;


  fprintf(outfile,"\n*** tstart() called on %s\n", name);
  fprintf(outfile,"*** at %s\n",ctime(&time_start));

  free(name);
}

/*!
** tstop(): Stop timer
**
** \param outfile = output file pointer.
**
** \ingroup CIOMR
*/
void tstop()
{
  int error;
  time_t total_time;
  time_t total_time_overall;
  struct tms total_tmstime;
  char *name;
  double user_s, sys_s;

  name = (char *) malloc(40 * sizeof(char));
  error = gethostname(name, 40);
  if(error != 0) strncpy(name,"nohostname", 11);

  time_end = time(NULL);
  total_time = time_end - time_start;
  total_time_overall = time_end - time_start_overall;

  times(&total_tmstime);
  const long clk_tck = sysconf(_SC_CLK_TCK);
  user_stop = ((double) total_tmstime.tms_utime)/clk_tck;
  sys_stop = ((double) total_tmstime.tms_stime)/clk_tck;

  user_s = user_stop - user_start;
  sys_s = sys_stop - sys_start;


  fprintf(outfile,"\n*** tstop() called on %s at %s", name, ctime(&time_end));

  /// print all module timings
  fprintf(outfile,"Module time:\n");
  fprintf(outfile,"\tuser time   = %10.2f seconds = %10.2f minutes\n",
          user_s, user_s/60.0);
  fprintf(outfile,"\tsystem time = %10.2f seconds = %10.2f minutes\n",
          sys_s, sys_s/60.0);
  fprintf(outfile,"\ttotal time  = %10d seconds = %10.2f minutes\n",
          (int)total_time, ((double) total_time)/60.0);

  user_s = user_stop - user_start_overall;
  sys_s = sys_stop - sys_start_overall;

  /// print all overall timings
  fprintf(outfile,"Total time:\n");
  fprintf(outfile,"\tuser time   = %10.2f seconds = %10.2f minutes\n",
          user_s, user_s/60.0);
  fprintf(outfile,"\tsystem time = %10.2f seconds = %10.2f minutes\n",
          sys_s, sys_s/60.0);
  fprintf(outfile,"\ttotal time  = %10d seconds = %10.2f minutes\n",
          (int)total_time_overall, ((double) total_time_overall)/60.0);


  free(name);

}

}

