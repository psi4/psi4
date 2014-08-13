#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "comex.h"

static int me;
static int nproc;
static int size[] = {2,4,8,16,32,64,128,256,512,1024,0}; /* 0 is sentinal */

#define PUT_FORWARD  0
#define PUT_BACKWARD 1
#define GET_FORWARD  2
#define GET_BACKWARD 3

static void fill_array(double *arr, int count, int which);
static void shift(size_t buffer_size, int op);


int main(int argc, char **argv)
{
    int i;

    comex_init_args(&argc, &argv);
    comex_group_rank(COMEX_GROUP_WORLD, &me);
    comex_group_size(COMEX_GROUP_WORLD, &nproc);

    if (0 == me) {
        printf("msg size (bytes)     avg time (milliseconds)    avg b/w (bytes/sec)\n");
    }

    if (0 == me) {
        printf("shifting put forward\n");
    }
    for (i=0; size[i]!=0; ++i) {
        shift(size[i], PUT_FORWARD);
    }

    if (0 == me) {
        printf("shifting put backward\n");
    }
    for (i=0; size[i]!=0; ++i) {
        shift(size[i], PUT_BACKWARD);
    }

    if (0 == me) {
        printf("shifting get forward\n");
    }
    for (i=0; size[i]!=0; ++i) {
        shift(size[i], GET_FORWARD);
    }

    if (0 == me) {
        printf("shifting get backward\n");
    }
    for (i=0; size[i]!=0; ++i) {
        shift(size[i], GET_BACKWARD);
    }

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


static void shift(size_t buffer_size, int op)
{
    void **dst_ptr;
    void **put_buf;
    void **get_buf;
    int i=0;
    double *times;
    double *result;
    double total_time=0;
    MPI_Comm comm = MPI_COMM_NULL;

    comex_group_comm(COMEX_GROUP_WORLD, &comm);
    dst_ptr = (void*)malloc(nproc * sizeof(void*));
    put_buf = (void*)malloc(nproc * sizeof(void*));
    get_buf = (void*)malloc(nproc * sizeof(void*));
    times = (double*)malloc(nproc * sizeof(double));
    result = (double*)malloc(nproc * sizeof(double));
    comex_malloc(dst_ptr, buffer_size, COMEX_GROUP_WORLD);
    comex_malloc(put_buf, buffer_size, COMEX_GROUP_WORLD);
    comex_malloc(get_buf, buffer_size, COMEX_GROUP_WORLD);

    /* initialize what we're putting */
    fill_array((double*)put_buf[me], buffer_size/sizeof(double), me);

    /* initialize time keepers */
    (void)memset(times, 0, nproc*sizeof(double));
    (void)memset(result, 0, nproc*sizeof(double));
    times[me] = MPI_Wtime()*1.0e6;

    /* the shift */
    switch (op) {
        case PUT_FORWARD:
            for (i=1; i<nproc; ++i) {
                int dst = (me+i)%nproc;
                comex_put(put_buf[me], dst_ptr[dst], buffer_size, dst, COMEX_GROUP_WORLD);
                comex_barrier(COMEX_GROUP_WORLD);
            }
            break;
        case PUT_BACKWARD:
            for (i=1; i<nproc; ++i) {
                int dst = me<i ? me-i+nproc : me-i;
                comex_put(put_buf[me], dst_ptr[dst], buffer_size, dst, COMEX_GROUP_WORLD);
                comex_barrier(COMEX_GROUP_WORLD);
            }
            break;
        case GET_FORWARD:
            for (i=1; i<nproc; ++i) {
                int dst = (me+i)%nproc;
                comex_get(dst_ptr[dst], get_buf[me], buffer_size, dst, COMEX_GROUP_WORLD);
                comex_barrier(COMEX_GROUP_WORLD);
            }
            break;
        case GET_BACKWARD:
            for (i=1; i<nproc; ++i) {
                int dst = me<i ? me-i+nproc : me-i;
                comex_get(dst_ptr[dst], get_buf[me], buffer_size, dst, COMEX_GROUP_WORLD);
                comex_barrier(COMEX_GROUP_WORLD);
            }
            break;
        default:
            comex_error("oops", 1);
    }

    /* calculate total time and average time */
    times[me] = MPI_Wtime()*1.0e6 - times[me];
    MPI_Allreduce(times, result, nproc, MPI_DOUBLE, MPI_SUM, comm);
    for (i=0; i<nproc; ++i) {
        total_time += times[i];
    }
    if (0 == me) {
        printf("%5zu                %6.2f                   %10.2f\n",
                buffer_size,
                total_time/nproc*1000,
                buffer_size*(nproc-1)/total_time);
    }

    comex_free(dst_ptr[me], COMEX_GROUP_WORLD);
    comex_free(put_buf[me], COMEX_GROUP_WORLD);
    comex_free(get_buf[me], COMEX_GROUP_WORLD);
    free(dst_ptr);
    free(put_buf);
    free(get_buf);
    free(times);
}
