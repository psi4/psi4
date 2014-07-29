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

/*! \file
    \ingroup DETCI
    \brief DETCI-specific timing routines
*/

#include <unistd.h>
#include <sys/time.h>
#include <cstdlib>
#include <cstdio>
#include "structs.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

double
wall_time_new(void)
{
  struct timeval tod;
  gettimeofday(&tod,0);
  return (double) (tod.tv_sec + 0.000001 * tod.tv_usec);
}

void 
init_time_new(struct detci_timings time)
{
 time.s1_total_time = time.s1_before_time = time.s1_after_time = 0.0;
 time.s2_total_time = time.s2_before_time = time.s2_after_time = 0.0;
 time.s3_total_time = time.s3_before_time = time.s3_after_time = 0.0;
 time.write_total_time = time.write_after_time = time.write_before_time = 0.0;
 time.read_total_time = time.read_after_time = time.read_before_time = 0.0;
 time.Hd_total_time = time.Hd_before_time = time.Hd_after_time = 0.0;
 time.total_before_time = time.total_after_time = 0.0;
}

void
print_time_new(struct detci_timings time)
{
  psi::fprintf(outfile,"\n");
  psi::fprintf(outfile,"        Total Time (s)     %%Time 		%%Relative\n");
  psi::fprintf(outfile," -----------------------------------------------------\n");
  psi::fprintf(outfile," Read      %lf\n", time.read_total_time);
  psi::fprintf(outfile," Write     %lf\n", time.write_total_time);
  psi::fprintf(outfile," Sigma1    %lf\n", time.s1_total_time);
  psi::fprintf(outfile," Sigma2    %lf\n", time.s2_total_time);
  psi::fprintf(outfile," Sigma3    %lf\n", time.s3_total_time);
  psi::fprintf(outfile," S1 Thread %lf\n", time.s1_mt_total_time);
  psi::fprintf(outfile," S2 Thread %lf\n", time.s2_mt_total_time);
  psi::fprintf(outfile," S3 Thread %lf\n", time.s3_mt_total_time);
  psi::fprintf(outfile,"\n");
}

}} // namespace psi::detci

