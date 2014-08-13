#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 *                              Copyright (c) 2006
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
 * $Id: multidma.c,v 1.1.2.1 2007-06-20 17:41:46 vinod Exp $
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
#ifdef WINDOWS_H
#   include <windows.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
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
#if HAVE_TIME_H
#   include <time.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif

#include <mpi.h>

#include "armci.h"
#include "message.h"

#define NDEBUG
/*#define LOG2FILE*/

typedef int t_elem; /* type of an array element */
#define SIZE_ELEM   sizeof(t_elem)

/* # of simultaneous message for phase A*/
#define SML_MSGS    (size < 9 ? size - 1 : 8)

#define MIN_MSG_SIZE    8
#define MAX_MSG_SIZE    (1024 * 1024)
#define MSG_COUNT       20

static inline void ARMCI_ASSERT(int error_code) {
    if (error_code) {
        fprintf(stderr, "ARMCI error %d\n", error_code);
        pause();
        ARMCI_Cleanup();
        MPI_Abort(MPI_COMM_WORLD, error_code);
    }
}

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
    int r;
    va_list ap;
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

struct stats {
    int size;   /* size, including trailing space for arrays */
    double put;
    double get;
    double *nb_put_a;
    double *nb_get_a;
    double *nb_put_b;
    double *nb_get_b;
};

void benchmark(int msg_size, struct stats *st)
{
    void **out_ptrs, **in_ptrs;
    double time_start, time_after_start, time_after_call, time_after_wait;
    double time2put, time2get;
    int i, j;
    armci_hdl_t handle;

    out_ptrs = malloc(sizeof(void*)*size);
    in_ptrs = malloc(sizeof(void*)*size); 

    log_debug("testing message size %d bytes\n", msg_size);

    log_debug("barrier A\n");
    armci_msg_barrier();

    /* allocate buffers */
    ARMCI_ASSERT(ARMCI_Malloc(out_ptrs, msg_size));
    ARMCI_ASSERT(ARMCI_Malloc(in_ptrs, msg_size));
    for (i = 0; i < msg_size; i++) {
        ((char *)out_ptrs[rank])[i] = (char)(rand() >> 24);
        ((char *)in_ptrs[rank])[i] = (char)(rand() >> 24);
    }

    /* warmup call */
    log_debug("barrier B\n");
    ARMCI_Barrier();

    ARMCI_ASSERT(ARMCI_Put(in_ptrs[rank], out_ptrs[second], msg_size, second));

    log_debug("barrier C: warmed up\n");
    ARMCI_Barrier();

    /* phase A */
    if (!rank) {
        /* measure blocking Put 0 -> 1 */
        time_start = armci_timer();
        time_after_start = armci_timer();

        ARMCI_ASSERT(ARMCI_Put(out_ptrs[0], in_ptrs[1], msg_size, 1));

        time_after_call = armci_timer();

        time2put =
            time_after_call - time_after_start + time_start - time_after_start;
        log_debug("PutA(%d bytes) 0 -> 1: %.8f\n", msg_size, time2put);
        st->put = time2put;
    }

    log_debug("barrier D\n");
    ARMCI_Barrier();

    if (!rank) {
        /* measure blocking Get 0 <- 1 */
        time_start = armci_timer();
        time_after_start = armci_timer();

        ARMCI_ASSERT(ARMCI_Get(out_ptrs[1], in_ptrs[0], msg_size, 1));

        time_after_call = armci_timer();

        time2get =
            time_after_call - time_after_start + time_start - time_after_start;
        log_debug("GetA(%d bytes) 0 <- 1: %.8f\n", msg_size, time2get);
        st->get = time2get;
    }

    for (i = 1; i <= SML_MSGS; i++) {
        log_debug("barrier E\n");
        ARMCI_Barrier();

        if (!rank) {
            /* put i messages to procs 0..i-1 */
            time_start = armci_timer();
            time_after_start = armci_timer();

            for (j = 1; j <= i; j++)
                ARMCI_ASSERT(ARMCI_NbPut(out_ptrs[0], in_ptrs[j],
                            msg_size, j, NULL));
            ARMCI_WaitAll();

            time_after_wait = armci_timer();
            time2put = time_after_wait - time_after_start + time_start -
                time_after_start;
            log_debug("NbPutA(%d bytes) 0 -> 1..%d: %.8f\n",
                    msg_size, i, time2put);
            st->nb_put_a[i - 1] = time2put;
        }
    }

    for (i = 1; i <= SML_MSGS; i++) {
        log_debug("barrier F\n");
        ARMCI_Barrier();

        if (!rank) {
            /* get i messages from procs 0..i-1 */
            time_start = armci_timer();
            time_after_start = armci_timer();

            for (j = 1; j <= i; j++)
                ARMCI_ASSERT(ARMCI_NbGet(out_ptrs[j], in_ptrs[0],
                            msg_size, j, NULL));
            ARMCI_WaitAll();

            time_after_wait = armci_timer();
            time2get = time_after_wait - time_after_start + time_start -
                time_after_start;
            log_debug("NbGetA(%d bytes) 0 <- 1..%d: %.8f\n",
                    msg_size, i, time2get);
            st->nb_get_a[i - 1] = time2get;
        }
    }

    /* phase B */
    log_debug("barrier G\n");
    ARMCI_Barrier();

    ARMCI_INIT_HANDLE(&handle);
    time_start = armci_timer();
    time_after_start = armci_timer();

    ARMCI_ASSERT(ARMCI_NbPut(out_ptrs[rank], in_ptrs[second],
                msg_size, second, &handle));
    ARMCI_ASSERT(ARMCI_Wait(&handle));

    time_after_wait = armci_timer();
    time2put = time_after_wait - time_after_start + time_start -
        time_after_start;

    log_debug("NbPutB(%d bytes) %d -> %d: %.8f\n",
            msg_size, rank, second, time2put);
    st->nb_put_b[rank] = time2put;

    log_debug("barrier H\n");
    ARMCI_Barrier();

    ARMCI_INIT_HANDLE(&handle);
    time_start = armci_timer();
    time_after_start = armci_timer();

    ARMCI_ASSERT(ARMCI_NbGet(out_ptrs[second], in_ptrs[rank],
                msg_size, second, &handle));
    ARMCI_ASSERT(ARMCI_Wait(&handle));

    time_after_wait = armci_timer();
    time2get = time_after_wait - time_after_start + time_start -
        time_after_start;

    log_debug("NbGetB(%d bytes) %d <- %d: %.8f\n",
            msg_size, rank, second, time2get);
    st->nb_get_b[rank] = time2get;

    log_debug("barrier I\n");
    ARMCI_Barrier();

    log_debug("finished benchmarking\n");
    /* clean up */
    ARMCI_Free(out_ptrs[rank]);
    ARMCI_Free(in_ptrs[rank]);

    free(out_ptrs);
    free(in_ptrs);
}


int main (int argc, char *argv[])
{
    int i, j, l;
    char buf[255];

    int dist, pos, time_seed;
    struct stats *st;

#ifdef  USE_EXP2
    int msg_sizes[MSG_COUNT];
#else
    int msg_sizes[] = {8, 16, 28, 52, 96, 180, 332, 616, 1144, 2124, 3952,
        7344, 13652, 25384, 47196, 87748, 163144, 303332, 563972, 1048576};
#endif

    armci_msg_init(&argc, &argv);
    rank = armci_msg_me();
    size = armci_msg_nproc();
    assert((size & 1) ^ 1); /* works with even number of processors only */
    log_debug("Message passing initialized\n");

    if (!rank) start_logging(argv[0]);

    ARMCI_ASSERT(ARMCI_Init());
    log_debug("ARMCI initialized\n");

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

    /* determines a pair/second for current processor */
    for (i = 0; i < HALFSIZE; i++) {
        if (p_srcs[i] == rank) second = p_dsts[i];
        if (p_dsts[i] == rank) second = p_srcs[i];
    }
    log_debug("second: %d\n", second);

    /* allocate memory tatistics */
    l = sizeof(struct stats) + 2 * sizeof(double) * (size + SML_MSGS);
    assert(st = malloc(l));
    st->size = l;
    st->nb_put_a = (double *)(((char *)st) + sizeof(struct stats));
    st->nb_get_a = st->nb_put_a + SML_MSGS;
    st->nb_put_b = st->nb_get_a + SML_MSGS;
    st->nb_get_b = st->nb_put_b + size;

    for (i = 0; i < MSG_COUNT; i++) {
        benchmark(msg_sizes[i], st);
        if (rank != 0) {
            MPI_Gather(st->nb_put_b + rank, 1, MPI_DOUBLE,
                        st->nb_put_b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(st->nb_get_b + rank, 1, MPI_DOUBLE,
                        st->nb_get_b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    if (!rank) {
        double min_put, max_put, mean_put;
        double min_get, max_get, mean_get;
        log_printf("msg_size = %d:\n", msg_sizes[i]);
        log_printf("Put = %.8f Get = %.8f\n", st->put, st->get);

        for (j = 0; j < SML_MSGS; j++)
            log_printf("NbPutA(0 -> 1..%d) = %.8f\n",
                    j + 1, st->nb_put_a[j]);
        for (j = 0; j < SML_MSGS; j++)
            log_printf("NbGetA(0 -> 1..%d) = %.8f\n",
                    j + 1, st->nb_get_a[j]);

        /* determine min, max and mean for phase B */
        min_put = max_put = mean_put = st->nb_put_b[0];
        min_get = max_get = mean_get = st->nb_get_b[0];

        for (j = 1; j < size; j++) {
            if (min_put > st->nb_put_b[j]) min_put = st->nb_put_b[j];
            if (max_put < st->nb_put_b[j]) max_put = st->nb_put_b[j];
            mean_put += st->nb_put_b[j];

            if (min_get > st->nb_get_b[j]) min_get = st->nb_get_b[j];
            if (max_get < st->nb_get_b[j]) max_get = st->nb_get_b[j];
            mean_get += st->nb_get_b[j];
        }
        mean_put /= size;
        mean_get /= size;

        log_printf("NbPutB min = %.8f max = %.8f mean = %.8f\n",
                min_put, max_put, mean_put);
        log_printf("NbGetB min = %.8f max = %.8f mean = %.8f\n",
                min_get, max_get, mean_get);
        log_printf("\n");
    }
    }

    /* clean up */
    ARMCI_Finalize();
    armci_msg_finalize();

    free(st);
    free(p_srcs);

    if (!rank) finish_logging();

    return 0;
}
