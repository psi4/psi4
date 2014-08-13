#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "armci.h"

static int me;
static int nproc;

#define PUT  0
#define GET 1
#define ACC 2

#define MAX_MESSAGE_SIZE 1024*1024
#define MEDIUM_MESSAGE_SIZE 8192
#define ITER_SMALL 100
#define ITER_LARGE 100

#define WARMUP 10
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
    armci_msg_init(&argc,&argv);
    ARMCI_Init_args(&argc, &argv);
    me = armci_msg_me();
    nproc = armci_msg_nproc();

    /* This test only works for two processes */

    assert(nproc == 2);

    if (0 == me) {
        printf("msg size (bytes)     avg time (us)    avg b/w (MB/sec)\n");
    }

    if (0 == me) {
        printf("#PNNL comex Put Test\n");
    }
    contig_test(MAX_MESSAGE_SIZE, PUT);

    if (0 == me) {
        printf("#PNNL comex Get Test\n");
    }
    contig_test(MAX_MESSAGE_SIZE, GET);
   
    if (0 == me) {
        printf("#PNNL comex Accumulate Test\n");
    }
    contig_test(MAX_MESSAGE_SIZE, ACC);
    
    ARMCI_Finalize();
    armci_msg_finalize();

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
    ARMCI_Malloc(dst_ptr, buffer_size);
    ARMCI_Malloc(put_buf, buffer_size);
    ARMCI_Malloc(get_buf, buffer_size);

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
                        ARMCI_Put(put_buf[me], dst_ptr[dst], msg_size,
                                dst);
                        break;
                    case GET:
                        ARMCI_Get(dst_ptr[dst], get_buf[me], msg_size,
                                dst);
                        break;
                    case ACC:
                        ARMCI_Acc(ARMCI_ACC_DBL, &scale, 
                                put_buf[me], dst_ptr[dst], msg_size,
                                dst);
                        break;
                    default:
                        ARMCI_Error("oops", 1);
                }

            }
        }
        /* calculate total time and average time */
        t_end = dclock();
        ARMCI_Barrier();


        if (0 == me) {
            printf("%8zu\t\t%6.2f\t\t%10.2f\n",
                    msg_size,
                    ((t_end  - t_start))/iter,
                    msg_size*iter/((t_end - t_start)));
        }
    }
    ARMCI_Free(dst_ptr[me]);
    ARMCI_Free(put_buf[me]);
    ARMCI_Free(get_buf[me]);
    free(dst_ptr);
    free(put_buf);
    free(get_buf);
    free(times);
}
