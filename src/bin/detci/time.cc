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
  fprintf(outfile,"\n");
  fprintf(outfile,"        Total Time (s)     %%Time 		%%Relative\n");
  fprintf(outfile," -----------------------------------------------------\n");
  fprintf(outfile," Read      %lf\n", time.read_total_time);
  fprintf(outfile," Write     %lf\n", time.write_total_time);
  fprintf(outfile," Sigma1    %lf\n", time.s1_total_time);
  fprintf(outfile," Sigma2    %lf\n", time.s2_total_time);
  fprintf(outfile," Sigma3    %lf\n", time.s3_total_time);
  fprintf(outfile," S1 Thread %lf\n", time.s1_mt_total_time);
  fprintf(outfile," S2 Thread %lf\n", time.s2_mt_total_time);
  fprintf(outfile," S3 Thread %lf\n", time.s3_mt_total_time);
  fprintf(outfile,"\n");
}

}} // namespace psi::detci

