#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <mpi.h>

#include "ga.h"

#define SCF F77_FUNC(scf,SCF)
#define TIMER F77_FUNC_(timer,TIMER)
#define UTIL_TIME F77_FUNC_(util_time,UTIL_TIME)
#define UTIL_MTIME F77_FUNC_(util_mtime,UTIL_MTIME)

double UTIL_TIME()
{
    static int first_call = 1;
    static double first_time, last_time, cur_time;
    double diff;

    if (first_call) {
        first_time = MPI_Wtime();
        first_call = 0;
        last_time  = -1e-9;
    }

    cur_time = MPI_Wtime();
    diff = cur_time - first_time;

    /* address crappy MPI_Wtime: consectutive calls must be at least 1ns apart  */
    if(diff - last_time < 1e-9) {
        diff +=1e-9;
    }
    last_time = diff;

    return diff; /* Add logic here for clock wrap */
}


double TIMER()
{
    return UTIL_TIME()*0.01;
}


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    GA_Initialize();

    SCF();

    GA_Terminate();
    MPI_Finalize();
    return 0;
}
