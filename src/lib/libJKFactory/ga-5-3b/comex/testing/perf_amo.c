/* Test Rmw Performance
 * The number of processes are increases from 2 to the number of 
 * processes present in the job */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include <mpi.h>

#include "comex.h"

static int me;
static int nproc;

#define FETCH_AND_ADD  0
#define FETCH_AND_ADD_LONG  1
#define SWAP 2
#define SWAP_LONG 3

#define MAX_MESSAGE_SIZE 1024
#define MEDIUM_MESSAGE_SIZE 8192
#define ITER_SMALL 10000
#define ITER_LARGE 10000

#define WARMUP 20
static void fill_array(double *arr, int count, int which);
static void rmw_test(size_t buffer_size, int op);

double dclock()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return(tv.tv_sec * 1.0e6 + (double)tv.tv_usec);
}

int main(int argc, char **argv)
{
    comex_init_args(&argc, &argv);
    comex_group_rank(COMEX_GROUP_WORLD, &me);
    comex_group_size(COMEX_GROUP_WORLD, &nproc);

    /* This test only works for two processes */

    if (0 == me) {
        printf("#Processes     avg time (us)\n");
        printf("\n\n");
    }

    if (0 == me) {
        printf("#PNNL armci Rmw-Fetch and Add Long Test\n");
        printf("\n\n");
    }
    rmw_test(MAX_MESSAGE_SIZE, FETCH_AND_ADD_LONG);

    if (0 == me)
        printf("\n\n");


    if (0 == me) {
        printf("#PNNL armci Rmw-Fetch and Add Test\n");
    }
    rmw_test(MAX_MESSAGE_SIZE, FETCH_AND_ADD);
    if (0 == me)
        printf("\n\n");
    
   
    if (0 == me) {
        printf("#PNNL armci Rmw-Swap Long Test\n");
    }
    rmw_test(MAX_MESSAGE_SIZE, SWAP_LONG);
   
    if (0 == me)
        printf("\n\n");
    
    
    if (0 == me) {
        printf("#PNNL armci Rmw-Swap Test\n");
    }
    rmw_test(MAX_MESSAGE_SIZE, SWAP);
    
    if (0 == me)
        printf("\n\n");
    comex_finalize();
    
    MPI_Finalize();

    return 0;
}


static void fill_array(double *arr, int count, int which)
{
    int i;

    for (i = 0; i < count; i++) {
        arr[i] = i * 8.23 + which * 2.89;
    }
}


static void rmw_test(size_t buffer_size, int op)
{
    void **dst_ptr;
    void **put_buf;
    void **get_buf;
    double *times;

    dst_ptr = (void*)malloc(nproc * sizeof(void*));
    put_buf = (void*)malloc(nproc * sizeof(void*));
    get_buf = (void*)malloc(nproc * sizeof(void*));
    times = (double*)malloc(nproc * sizeof(double));
    comex_malloc(dst_ptr, buffer_size, COMEX_GROUP_WORLD);
    comex_malloc(put_buf, buffer_size, COMEX_GROUP_WORLD);
    comex_malloc(get_buf, buffer_size, COMEX_GROUP_WORLD);

    /* initialize what we're putting */
    fill_array((double*)put_buf[me], buffer_size/sizeof(double), me);

    /* All processes perform Rmw on process 0*/
    int dst = 0;
    double t_start, t_end;

    int j;
    int iter = ITER_LARGE;

    int part_proc;

    for (part_proc = 2; part_proc <= nproc; part_proc *= 2) {
        if (me < part_proc) {
            for (j= 0; j < iter + WARMUP; ++j) {

                if (WARMUP == j) {
                    t_start = dclock();
                }

                switch (op) {
                    case FETCH_AND_ADD:
                        comex_rmw(COMEX_FETCH_AND_ADD,
                                put_buf[me], dst_ptr[dst], 1, dst, COMEX_GROUP_WORLD);
                        break;
                    case FETCH_AND_ADD_LONG:
                        comex_rmw(COMEX_FETCH_AND_ADD_LONG,
                                put_buf[me], dst_ptr[dst], 1, dst, COMEX_GROUP_WORLD);
                        break;
                    case SWAP:
                        comex_rmw(COMEX_SWAP,
                                put_buf[me], dst_ptr[dst], 1, dst, COMEX_GROUP_WORLD);
                        break;
                    case SWAP_LONG:
                        comex_rmw(COMEX_SWAP_LONG,
                                put_buf[me], dst_ptr[dst], 1, dst, COMEX_GROUP_WORLD);
                        break;
                    default:
                        comex_error("oops", 1);
                }
            }
        }
        comex_barrier(COMEX_GROUP_WORLD);
        /* calculate total time and average time */
        t_end = dclock();


        if (0 == me) {
            printf("%5d\t\t%6.2f\n",
                    part_proc,
                    ((t_end  - t_start))/iter);
        }
    }
    
    comex_free(dst_ptr[me], COMEX_GROUP_WORLD);
    comex_free(put_buf[me], COMEX_GROUP_WORLD);
    comex_free(get_buf[me], COMEX_GROUP_WORLD);
    free(dst_ptr);
    free(put_buf);
    free(get_buf);
    free(times);
}
