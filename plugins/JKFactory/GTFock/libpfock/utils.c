#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mkl.h>
#include <mkl_lapack.h>
#include <mpi.h>

#include "config.h"
#include "utils.h"


void printmatRM (char *name, double *A, int nrows, int ncols)
{
    fprintf (stderr, "name[RM]: %s\n", name);
    int r;
    int c;
    for (r = 0; r < nrows; r++)
    {
        for (c = 0; c < ncols; c++)
        {
            fprintf (stderr, "%10.6e ", A[r * ncols + c]);
        }
        fprintf (stderr, "\n");
    }
}


void printmatCM (char *name, double *A, int nrows, int ncols)
{
    fprintf (stderr, "name[CM]: %s\n", name);
    int r;
    int c;
    for (c = 0; c < ncols; c++)
    {
        for (r = 0; r < nrows; r++)
        {
            fprintf (stderr, "%10.6e ", A[r * ncols + c]);
        }
        fprintf (stderr, "\n");
    }
}


void initomp (int nthreads, int verbose)
{
    char schedule[1024];

    if (verbose == 1)
    {
        sprintf (schedule, "KMP_AFFINITY=granularity=fine,compact,verbose");
    }
    else
    {
        sprintf (schedule, "KMP_AFFINITY=granularity=fine,compact");
    }
    kmp_set_defaults (schedule);
    mkl_set_num_threads (nthreads);
    omp_set_num_threads (nthreads);
}


void *my_malloc (size_t size, char *file, int line)
{
    void *p;
    int myrank;
    static double totalsize = 0.0;;

    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    totalsize += size/1024.0/1024.0;
    if (myrank == 0)
    {
        fprintf (stderr, "**** Allocated %.3lf MB at %s %d, totalsize %.3lf MB\n",
            size/1024.0/1024.0, file, line, totalsize);        
    }
    fflush (stdout);
    p = malloc (size);
    
    if (p == NULL)
    {
        fprintf (stderr, "**** Rank %d failed to allocate %.3lf MB at %s %d, totalsize %.3lf MB\n",
            myrank, size/1024.0/1024.0, file, line, totalsize);
    }

    return p;
}
