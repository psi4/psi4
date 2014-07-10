#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 *                                Copyright (c) 2006
 *                      Pacific Northwest National Laboratory,
 *                           Battelle Memorial Institute.
 *                              All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met: - Redistributions of source code must retain the above
 * copyright notice, this list of conditions and the following disclaimer.
 * 
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 * - Neither the name of the Battelle nor the names of its contributors
 *   may be used to endorse or promote products derived from this software
 *   without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * $Id$
 */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STDARG_H
#   include <stdarg.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#elif HAVE_WINDOWS_H
#   include <windows.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_TIME_H
#   include <time.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif

#include <mpi.h>

#include "armci.h"
#include "message.h"

extern double exp2(double);
extern double round(double);
extern double log2(double);
#define NDEBUG
/*#define LOG2FILE*/

typedef int t_elem; /* type of an array element */
#define SIZE_ELEM   sizeof(t_elem)


#define STRIDE_OFF  (SIZE_ELEM * 4 - 1)

#define MIN_MSG_SIZE    8
#define MAX_MSG_SIZE    (1024 * 1024)
#define MSG_COUNT       20

int armci_error_code;
#define ARMCI_ASSERT(error_code) if ((armci_error_code = error_code)) {   \
        fprintf(stderr, "ARMCI error %d\n", armci_error_code);pause();         \
        ARMCI_Cleanup(); MPI_Abort(MPI_COMM_WORLD, armci_error_code); }

#define FIX_TIME(t) if (t < 0.0) t = 0.0;

int size, rank, second;

#define ITERS           18
#define ITER_STEPS      20
double iterations_times[ITERS];
int iterations[ITERS];

int *p_srcs, *p_dsts;

enum {CONT_PUT, CONT_GET,
    STRIDED_PUT, STRIDED_GET, STRIDED_ACC,
    VECTOR_PUT, VECTOR_GET, VECTOR_ACC};
#define OPS_COUNT   (STRIDED_ACC + 1)
#define NON_CONT(op)  (op > CONT_GET)

enum {NOWORK, TOTAL, OVERLAP};
#define STATS_COUNT (OVERLAP + 1)

/* prints formatted numbered message with processor's rank */
int log_debug(const char *fmt, ...)
{
    int r = 0;
#ifndef NDEBUG
    static int log_counter = 1;
    va_list ap;
    va_start(ap, fmt);

    printf("%03d@%1d: ", log_counter++, rank);
    r = vprintf(fmt, ap);

    va_end(ap);
#endif
    return r;
}

FILE *log_file = NULL;

void start_logging(const char *fname)
{
#ifdef  LOG2FILE
    char exe_name[255];
    char log_path[255];
    int i;
    char k;

    strcpy(exe_name, fname);
    if (exe_name[strlen(exe_name) - 2] == '.') /* remove .x */
        exe_name[strlen(exe_name) - 2] = 0;

    if (exe_name[0] == '/') { /* full path given */
        for (i = ((int)strlen(exe_name)) - 1, k = -1; i >= 0; i--)
            if (exe_name[i] == '/') {
                if (k == -1) k = i + 1;
                else {
                    exe_name[i] = 0;
                    break;
                }
            }
        log_debug("exe: path=%s, name=%s\n", exe_name, exe_name + k);
        sprintf(log_path, "%s/data/%s.dat", exe_name, exe_name + k);
    } else { /* only executable name */
        sprintf(log_path, "../data/%s.dat", exe_name);
    }
    log_debug("log: %s\n", log_path);

    log_file = fopen(log_path, "w");

    if (!log_file) {
        perror("cannot open log file");
        abort();
    }
#else
    log_file = stderr;
#endif
}

void finish_logging()
{
    fclose(log_file);
}

/* prints formatted message to ../data/<prog>.dat */
int log_printf(const char *fmt, ...)
{
    va_list ap;
    int r;
    
    va_start(ap, fmt);

    if (log_file)
        r = vfprintf(log_file, fmt, ap);
    else {
        fprintf(stderr, "warning: logging is not enabled for this process\n");
        r = vfprintf(stderr, fmt, ap);
    }

    va_end(ap);
    return r;
}


/* computes approximate time of n iterations for variable n */
void time_iterations()
{
        double time_start, time_after_start, time_stop;
        int i, j, k, l;

        for (i = 0, j = 1; i < ITERS; i++, j *= 2) {
                time_start = armci_timer();
                time_after_start = armci_timer();

                for (l = 0, k = rand(); l < j; l++) k *= rand();

                time_stop = armci_timer();
                iterations_times[i] = time_stop - time_after_start +
                    time_start - time_after_start;
                FIX_TIME(iterations_times[i]);
                iterations[i] = j;

                log_debug("it takes %.8f sec to iterate %d times\n",
                        iterations_times[i], iterations[i]);
        }
}


/* computes useful overlap time for contiguous/vector/strided arrays
 *  * op - operation
 *   * msg_size - size of a message/ 1st dimension (bytes)
 *    * size2 - not used for contiguous arrays
 *     *       - size of 2nd dimension for strided arrays (bytes)
 *      *       - # of vector segments for vectors
 *       * returns pointer to static array of stats (STATS_COUNT doubles)
 *        */
double * benchmark(int op, int msg_size, int size2)
{
    static double stats[STATS_COUNT]; /* return statistics in static array */

    void **array_ptrs;
    int stride_dist, block_sizes[2], scale = 2;
    int i=0, j=0, k=0, l=0, less=0, more=0;
    double time_start=0, time_after_start=0, time_after_call=0,
           time_after_work=0, time_after_wait=0;
    double time2call_nw=0, time2wait_nw = 1.0, time_total_nw=0;
    double time2call_fw, time2work_fw, time2wait_fw, time_total_fw;
    armci_hdl_t handle;

    array_ptrs = malloc(sizeof(void*)*size);

    log_debug("barrier O\n");
    armci_msg_barrier();
    /* initialize: obtain remote address and generate random array */
    switch (op) {
        case CONT_PUT:
        case CONT_GET:
            ARMCI_ASSERT(ARMCI_Malloc(array_ptrs, msg_size));
            for (i = 0; i < msg_size; i++)
                ((char *)array_ptrs[rank])[i] = (char)(rand() >> 24);
            break;

        /* 2D strided array of ints */
        case STRIDED_PUT:
        case STRIDED_GET:
        case STRIDED_ACC:
            block_sizes[0] = msg_size;
            block_sizes[1] = size2;
            stride_dist = STRIDE_OFF + msg_size;
            log_debug("strided: dim1 = %d (%d bytes), dim2 = %d, stride = %d\n",
                    msg_size / SIZE_ELEM, msg_size, size2, stride_dist);
            ARMCI_ASSERT(ARMCI_Malloc(array_ptrs,
                        (size2 - 1) * stride_dist + msg_size));

            for (i = 0; i < size2; i++)
                for (j = 0; j < (msg_size / ((int)SIZE_ELEM)); j++) {
                    l = stride_dist * i + SIZE_ELEM * i;
                    *(int *)((char *)array_ptrs[rank] + l) = rand();
                }
            break;
    }

    /* warm up call */
    log_debug("barrier A\n");
    armci_msg_barrier();
    if (second != -1) {
            log_debug("testing message size %d bytes\n", msg_size);
            switch (op) {
                case CONT_PUT:
                    ARMCI_INIT_HANDLE(&handle);
                    time_start = armci_timer();
                    time_after_start = armci_timer();

                    ARMCI_ASSERT(ARMCI_NbPut(array_ptrs[rank],
                                array_ptrs[second], msg_size,
                                second, &handle));
                    time_after_call = armci_timer();

                    ARMCI_ASSERT(ARMCI_Wait(&handle));
                    time_after_wait = armci_timer();
                    break;

                case CONT_GET:
                    ARMCI_INIT_HANDLE(&handle);
                    time_start = armci_timer();
                    time_after_start = armci_timer();

                    ARMCI_ASSERT(ARMCI_NbGet(array_ptrs[second],
                                array_ptrs[rank], msg_size,
                                second, &handle));
                    time_after_call = armci_timer();

                    ARMCI_ASSERT(ARMCI_Wait(&handle));
                    time_after_wait = armci_timer();
                    break;

                case STRIDED_PUT:
                    ARMCI_INIT_HANDLE(&handle);
                    time_start = armci_timer();
                    time_after_start = armci_timer();

                    ARMCI_ASSERT(ARMCI_NbPutS(array_ptrs[rank], &stride_dist,
                                array_ptrs[second], &stride_dist,
                                block_sizes, 1, second, &handle));
                    time_after_call = armci_timer();

                    ARMCI_ASSERT(ARMCI_Wait(&handle));
                    time_after_wait = armci_timer();
                    break;

                case STRIDED_GET:
                    ARMCI_INIT_HANDLE(&handle);
                    time_start = armci_timer();
                    time_after_start = armci_timer();

                    ARMCI_ASSERT(ARMCI_NbGetS(array_ptrs[second], &stride_dist,
                                array_ptrs[second], &stride_dist,
                                block_sizes, 1, second, &handle));
                    time_after_call = armci_timer();

                    ARMCI_ASSERT(ARMCI_Wait(&handle));
                    time_after_wait = armci_timer();
                    break;

                case STRIDED_ACC:
                    ARMCI_INIT_HANDLE(&handle);
                    time_start = armci_timer();
                    time_after_start = armci_timer();

                    ARMCI_ASSERT(ARMCI_NbAccS(ARMCI_ACC_INT, &scale,
                                array_ptrs[rank], &stride_dist,
                                array_ptrs[second], &stride_dist,
                                block_sizes, 1, second, &handle));
                    time_after_call = armci_timer();

                    ARMCI_ASSERT(ARMCI_Wait(&handle));
                    time_after_wait = armci_timer();

                    break;
            }

            time2call_nw = time_after_call - time_after_start + time_start -
                time_after_start;
            time2wait_nw = time_after_wait - time_after_call + time_start -
                time_after_start;
            time_total_nw = time_after_wait - time_after_start + time_start -
                time_after_start;

            log_debug("time (warm up): %.8f call, %.8f wait, %.8f total\n",
            time2call_nw, time2wait_nw, time_total_nw);
    }

    log_debug("barrier B\n");
    armci_msg_barrier();
    if (second != -1) {
            /* no work */
            ARMCI_INIT_HANDLE(&handle);
            switch (op) {
                case CONT_PUT:
                    time_start = armci_timer();
                    time_after_start = armci_timer();

                    ARMCI_ASSERT(ARMCI_NbPut(array_ptrs[rank],
                                array_ptrs[second], msg_size,
                                second, &handle));
                    time_after_call = armci_timer();

                    ARMCI_ASSERT(ARMCI_Wait(&handle));
                    time_after_wait = armci_timer();
                    break;

                case CONT_GET:
                    time_start = armci_timer();
                    time_after_start = armci_timer();

                    ARMCI_ASSERT(ARMCI_NbGet(array_ptrs[second],
                                array_ptrs[rank], msg_size,
                                second, &handle));
                    time_after_call = armci_timer();

                    ARMCI_ASSERT(ARMCI_Wait(&handle));
                    time_after_wait = armci_timer();
                    break;

                case STRIDED_PUT:
                    time_start = armci_timer();
                    time_after_start = armci_timer();

                    ARMCI_ASSERT(ARMCI_NbPutS(array_ptrs[rank], &stride_dist,
                                array_ptrs[second], &stride_dist,
                                block_sizes, 1, second, &handle));
                    time_after_call = armci_timer();

                    ARMCI_ASSERT(ARMCI_Wait(&handle));
                    time_after_wait = armci_timer();
                    break;

                case STRIDED_GET:
                    time_start = armci_timer();
                    time_after_start = armci_timer();

                    ARMCI_ASSERT(ARMCI_NbGetS(array_ptrs[second], &stride_dist,
                                array_ptrs[rank], &stride_dist,
                                block_sizes, 1, second, &handle));
                    time_after_call = armci_timer();

                    ARMCI_ASSERT(ARMCI_Wait(&handle));
                    time_after_wait = armci_timer();
                    break;

                case STRIDED_ACC:
                    time_start = armci_timer();
                    time_after_start = armci_timer();

                    ARMCI_ASSERT(ARMCI_NbAccS(ARMCI_ACC_INT, &scale,
                                array_ptrs[rank], &stride_dist,
                                array_ptrs[second], &stride_dist,
                                block_sizes, 1, second, &handle));
                    time_after_call = armci_timer();

                    ARMCI_ASSERT(ARMCI_Wait(&handle));
                    time_after_wait = armci_timer();
                    break;
            }

            time2call_nw = time_after_call - time_after_start + time_start -
                time_after_start;
            FIX_TIME(time2call_nw);
            time2wait_nw = time_after_wait - time_after_call + time_start -
                time_after_start;
            FIX_TIME(time2wait_nw);
            time_total_nw = time_after_wait - time_after_start + time_start -
                time_after_start;
            FIX_TIME(time_total_nw);

            log_debug("time (no work): %.8f call, %.8f wait, %.8f total\n",
                    time2call_nw, time2wait_nw, time_total_nw);
    }

    /* only perform tests if wait time is not 0 */
    if (time2wait_nw > 0.0) {
    /* time2wait_nw is always 1.0 on seconds (receiving nodes) */
        double overlaps[ITER_STEPS], totals[ITER_STEPS];
        if (second !=  -1) {
            /* compute approximate range of iterations */
            less = 0, more = iterations[ITERS - 1];
            assert(time2wait_nw < iterations_times[ITERS - 1]);

            for (i = 0; i < ITERS; i++)
                if (time2wait_nw > iterations_times[i])
                    less = iterations[i];
                else
                    break;
            for (i = 0; i < ITERS; i++)
                if (time2wait_nw < iterations_times[ITERS - i - 1])
                    more = iterations[ITERS - i - 1];
                else
                    break;

            log_debug("wait time (%.8f) is between %d and %d iterations\n",
                    time2wait_nw, less, more);
        }

        /* benchmark ITER_STEPS steps within computed range */
        for (i = 0, j = less; i < ITER_STEPS;
             i++, j += (more - less) / (ITER_STEPS - 1)) {
            /* time noneblocking call with j interations of fake work */
            log_debug("barrier C\n");
            armci_msg_barrier();
            if (second != -1) {
                ARMCI_INIT_HANDLE(&handle);
                switch (op) {
                    case CONT_PUT:
                        time_start = armci_timer();
                        time_after_start = armci_timer();

                        ARMCI_ASSERT(ARMCI_NbPut(array_ptrs[rank],
                                    array_ptrs[second], msg_size,
                                    second, &handle));
                        time_after_call = armci_timer();

                        for (l = 0, k = rand(); l < j; l++) k *= rand();
                        time_after_work = armci_timer();

                        ARMCI_ASSERT(ARMCI_Wait(&handle));
                        time_after_wait = armci_timer();
                        break;

                    case CONT_GET:
                        time_start = armci_timer();
                        time_after_start = armci_timer();

                        ARMCI_ASSERT(ARMCI_NbGet(array_ptrs[second],
                                    array_ptrs[rank], msg_size,
                                    second, &handle));
                        time_after_call = armci_timer();

                        for (l = 0, k = rand(); l < j; l++) k *= rand();
                        time_after_work = armci_timer();

                        ARMCI_ASSERT(ARMCI_Wait(&handle));
                        time_after_wait = armci_timer();
                        break;

                case STRIDED_PUT:
                        time_start = armci_timer();
                        time_after_start = armci_timer();

                       ARMCI_ASSERT(ARMCI_NbPutS(array_ptrs[rank], &stride_dist,
                                   array_ptrs[second], &stride_dist,
                                   block_sizes, 1, second, &handle));
                        time_after_call = armci_timer();

                        for (l = 0, k = rand(); l < j; l++) k *= rand();
                        time_after_work = armci_timer();

                        ARMCI_ASSERT(ARMCI_Wait(&handle));
                        time_after_wait = armci_timer();
                        break;

                case STRIDED_GET:
                        time_start = armci_timer();
                        time_after_start = armci_timer();

                        ARMCI_ASSERT(ARMCI_NbGetS(array_ptrs[second],
                                   &stride_dist, array_ptrs[rank], &stride_dist,
                                   block_sizes, 1, second, &handle));
                        time_after_call = armci_timer();

                        for (l = 0, k = rand(); l < j; l++) k *= rand();
                        time_after_work = armci_timer();

                        ARMCI_ASSERT(ARMCI_Wait(&handle));
                        time_after_wait = armci_timer();

                        break;

                case STRIDED_ACC:
                        time_start = armci_timer();
                        time_after_start = armci_timer();

                        ARMCI_ASSERT(ARMCI_NbAccS(ARMCI_ACC_INT, &scale,
                                    array_ptrs[rank], &stride_dist,
                                    array_ptrs[second], &stride_dist,
                                    block_sizes, 1, second, &handle));
                        time_after_call = armci_timer();

                        for (l = 0, k = rand(); l < j; l++) k *= rand();
                        time_after_work = armci_timer();

                        ARMCI_ASSERT(ARMCI_Wait(&handle));
                        time_after_wait = armci_timer();
                        break;
                }

                time2call_fw = time_after_call - time_after_start + time_start -
                    time_after_start;
                FIX_TIME(time2call_fw);
                time2work_fw = time_after_work - time_after_call + time_start -
                    time_after_start;
                FIX_TIME(time2work_fw);
                time2wait_fw = time_after_wait - time_after_work + time_start -
                    time_after_start;
                FIX_TIME(time2wait_fw);
                time_total_fw = time_after_wait - time_after_start +
                    time_start - time_after_start;
                FIX_TIME(time_total_fw);

                log_debug("time (%d iters): %.8f call, %.8f work, "
                        "%.8f wait %.8f total\n", j, time2call_fw, time2work_fw,
                        time2wait_fw, time_total_fw);

                overlaps[i] = time2work_fw;
                totals[i] = time_total_fw;
            }
        }

        /* pick overlap with closest total (less or equal) */
        if (second != -1) {
                double closest_total, closest_overlap;
                double smallest_total = totals[ITER_STEPS - 1],
                       smallest_overlap = overlaps[ITER_STEPS - 1];
                for (i = ITER_STEPS - 1; i >= 0; i--) {
                        closest_total = totals[i];
                        closest_overlap = overlaps[i];
                        if (closest_total < smallest_total) {
                            smallest_total = closest_total;
                            smallest_overlap = closest_overlap;
                        }
                        if (closest_total <= time_total_nw) break;
                }
                if (closest_total > time_total_nw) {
                    closest_total = smallest_total;
                    closest_overlap = smallest_overlap;
                }
                stats[NOWORK]   = time_total_nw;
                stats[TOTAL]    = closest_total;
                stats[OVERLAP]  = closest_overlap;
        }
    } else {
        if (second != -1) {
            for (i = 0; i < ITER_STEPS; i++) {
                log_debug("barrier C0\n");
                armci_msg_barrier();
            }
            stats[NOWORK]   = time_total_nw;
            stats[TOTAL]    = 0;
            stats[OVERLAP]  = 0;
        }
    }


    ARMCI_ASSERT(ARMCI_Free(array_ptrs[rank]));
    free(array_ptrs);

    log_debug("barrier D\n");
    armci_msg_barrier();

    return stats;
}



int main (int argc, char *argv[])
{
    int i, j, k, l;
    double u;
    char buf[255];

    int dist, pos, time_seed;

    int msg_sizes[MSG_COUNT], dim1_sizes[MSG_COUNT], dim2[MSG_COUNT], mul_elem;
    double *stats=NULL, *stats_all=NULL;
    double from_log = log2(MIN_MSG_SIZE);
    double to_log   = log2(MAX_MSG_SIZE);
    double step_log = (to_log - from_log) / (MSG_COUNT - 1);

    armci_msg_init(&argc, &argv);
    rank = armci_msg_me();
    size = armci_msg_nproc();
    assert((size & 1) ^ 1); /* works with even number of processors only */
    log_debug("Message passing initialized\n");

    ARMCI_ASSERT(ARMCI_Init());
    log_debug("ARMCI initialized\n");

    if (!rank) start_logging(argv[0]);

    /* generate MSG_COUNT message sizes MIN_MSG_SIZE thru MAX_MSG_SIZE */
    for (i = 0, u = from_log; i < MSG_COUNT; i++, u += step_log) {
        mul_elem = round(exp2(u));
        msg_sizes[i] = mul_elem % ((int)SIZE_ELEM)
            ? (mul_elem / ((int)SIZE_ELEM) + 1) * ((int)SIZE_ELEM)
            : mul_elem; /* multiple of SIZE_ELEM */
    }

    /* generate MSG_COUNT respective dim1 sizes and dim2 for strided */
    for (i = 0; i < MSG_COUNT; i++) {
        mul_elem = msg_sizes[i] / SIZE_ELEM;
        mul_elem = sqrt(2.0 * mul_elem);
        dim1_sizes[i] = mul_elem * SIZE_ELEM;
        dim2[i] = mul_elem / 2;
    }

    /* print msg_sizes and appropriate derivatives (debug mode only) */
    if (!rank) {
        log_debug("msg_sizes:\n");
        for (i = 0; i < MSG_COUNT; i++)
            log_debug("cont: %d bytes | strided: %d bytes X %d\n",
                    msg_sizes[i], dim1_sizes[i], dim2[i]);
    }

    /* inialize PRNG, use seed generated on processor 0 for uniform sequence */
    time_seed = time(NULL);
    MPI_Bcast(&time_seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    srand(time_seed); rand();
    log_debug("seed: %d\n", time_seed);

    /* generate random pairs of processors */
#define HALFSIZE    (size / 2)
    assert(p_srcs = malloc(sizeof(int) * size));
    for (i = 0; i < size; i++) p_srcs[i] = -1;
    p_dsts = p_srcs + HALFSIZE;

    for (i = 0, j = size - 1, pos = 0; i < size; i++, j--) {
        dist = round((double)rand() * j / RAND_MAX + 1); /* random 1..j */

        for (l = 0; l < dist; ) {
            pos = (pos + 1 == size) ? 0 : pos + 1;
            if ((p_srcs[pos] == -1) && (pos != i)) l++;
        }
        p_srcs[pos] = i;
    }

    for (i = 0, j = 0; i < HALFSIZE; i++)
        j += sprintf(buf + j, " %d->%d", p_srcs[i], p_dsts[i]);
    log_debug("random pairs:%s\n", buf);

    /* time interations: 1 thru ITERS */
    time_iterations();

    /* determine if processor initiates communication and where it sends to,
     *      * -1 for second(receiver) */
    second = -1;
    for (i = 0; i < HALFSIZE; i++) if (p_srcs[i] == rank) second = p_dsts[i];
    log_debug("second: %d\n", second);

    /* allocate memory for statisticis */
#define MSG_OFF (STATS_COUNT * size)
#define OPS_OFF (MSG_OFF * MSG_COUNT)
    assert(stats_all = malloc(OPS_COUNT * OPS_OFF * sizeof(double)));

    for (i = 0; i < OPS_COUNT; i++)
        for (j = 0; j < MSG_COUNT; j++) {
            switch (i) {
                case CONT_PUT:
                case CONT_GET:
                    stats = benchmark(i, msg_sizes[j], 0);
                    log_debug("stats: %8d | %.8f | %.8f | %.8f | %.2f\n",
                            msg_sizes[j], stats[NOWORK], stats[TOTAL],
                            stats[OVERLAP],
                            100.0 * stats[OVERLAP] / stats[TOTAL]);
                    break;

                case STRIDED_PUT:
                case STRIDED_GET:
                case STRIDED_ACC:
                     stats = benchmark(i, dim1_sizes[j], dim2[j]);
                     log_debug("stats: %8d | %.8f | %.8f | %.8f | %.2f\n",
                            dim1_sizes[j] * dim2[j], stats[NOWORK],
                            stats[TOTAL], stats[OVERLAP],
                            100.0 * stats[OVERLAP] / stats[TOTAL]);
                    break;
            }
            MPI_Gather(stats, STATS_COUNT, MPI_DOUBLE,
                        stats_all + i * OPS_OFF + j * MSG_OFF,
                        STATS_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

    if (!rank)
        for (l = 0; l < HALFSIZE; l++) { /* interate thru pairs */
            log_printf("for pair of processors %d -> %d:\n", p_srcs[l], p_dsts[l]);

            for (i = 0; i < OPS_COUNT; i++) { /* iterate thru operations */
                switch (i) {
                        case CONT_PUT:
                            log_printf("ARMCI_NbPut\n");
                            break;

                        case CONT_GET:
                            log_printf("ARMCI_NbGet\n");
                            break;

                        case STRIDED_PUT:
                            log_printf("ARMCI_NbPutS\n");
                            break;

                        case STRIDED_GET:
                            log_printf("ARMCI_NbGetS\n");
                            break;

                        case STRIDED_ACC:
                            log_printf("ARMCI_NbAccS\n");
                            break;
                }
                log_printf("msg size |   nowork   |    total   |   overlap  |"
                        " ratio\n");
                log_printf("---------+------------+------------+------------+"
                        "------\n");

                for (j = 0; j < MSG_COUNT; j++) { /* iterate thru msg sizes */
                    k = i * OPS_OFF + j * MSG_OFF + p_srcs[l] * STATS_COUNT;
                    log_printf("%8d | %.8f | %.8f | %.8f | %.2f\n",
                            NON_CONT(i) ? dim1_sizes[j] * dim2[j]: msg_sizes[j],
                            stats_all[k + NOWORK], stats_all[k + TOTAL],
                            stats_all[k + OVERLAP],
                            (stats_all[k + NOWORK] < stats_all[k + TOTAL]) ||
                            (stats_all[k + TOTAL] <= 0.0)
                            ? 0 : 100.0 * stats_all[k + OVERLAP] /
                            stats_all[k + TOTAL]);
                }
                log_printf("\n");
            }
        }

    if (!rank) finish_logging();

    ARMCI_Finalize();
    armci_msg_finalize();

    free(p_srcs);
    free(stats_all);

    return 0;
}
