#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include <mpi.h>

#include "comex.h"

static int me;
static int nproc;

#define PUT  0
#define GET 1
#define ACC 2

#define MAX_MESSAGE_SIZE 1024*1024
#define MEDIUM_MESSAGE_SIZE 8192
#define ITER_SMALL 1000
#define ITER_LARGE 100

#define WARMUP 20
static void fill_array(double *arr, int count, int which);
static void contig_test(size_t buffer_size, int op);

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

    assert(nproc == 2);

    if (0 == me) {
        printf("msg size (bytes)     avg time (us)    avg b/w (MB/sec)\n");
    }
    printf("\n\n");

    if (0 == me) {
        printf("#PNNL armci Put Test\n");
    }
    contig_test(MAX_MESSAGE_SIZE, PUT);
    printf("\n\n");


    if (0 == me) {
        printf("#PNNL armci Get Test\n");
    }
    contig_test(MAX_MESSAGE_SIZE, GET);
   
   
    if (0 == me) {
        printf("#PNNL armci Accumulate Test\n");
    }
    contig_test(MAX_MESSAGE_SIZE, ACC);
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


static void contig_test(size_t buffer_size, int op)
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

    size_t msg_size;

    int dst = 1;
    double scale = 1.0;
    for (msg_size = 16; msg_size <= buffer_size; msg_size *= 2) {

        int j;
        int iter = msg_size > MEDIUM_MESSAGE_SIZE ? ITER_LARGE : ITER_SMALL;

        double t_start, t_end;
        if (0 == me) {
            for (j= 0; j < iter + WARMUP; ++j) {

                if (WARMUP == j) {
                    t_start = dclock();
                }

                switch (op) {
                    case PUT:
                        comex_put(put_buf[me], dst_ptr[dst], msg_size, dst, COMEX_GROUP_WORLD);
                        break;
                    case GET:
                        comex_get(dst_ptr[dst], get_buf[me], msg_size, dst, COMEX_GROUP_WORLD);
                        break;
                    case ACC:
                        comex_acc(COMEX_ACC_DBL, &scale, 
                                put_buf[me], dst_ptr[dst], msg_size, dst, COMEX_GROUP_WORLD);
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
            printf("%5zu\t\t%6.2f\t\t%10.2f\n",
                    msg_size,
                    ((t_end  - t_start))/iter,
                    msg_size*(nproc-1)*iter/((t_end - t_start)));
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
