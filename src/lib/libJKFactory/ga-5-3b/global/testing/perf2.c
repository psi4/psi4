/**
 * @author Jeff Daily, PNNL
 *
 * This is meant to directly compare against ARMCI's perf.x in order to
 * evaluate the overhead of GA for one-sided calls.
 */
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#elif HAVE_WINDOWS_H
#   include <windows.h>
#   define sleep(x) Sleep(1000*(x))
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif

#include "ga.h"
#include "armci.h"
#include "message.h"
#include "mp3.h"

#define SIZE 640 /**< must be >= biggest chunk + 128 */
#define CHUNK_NUM 28
#define TIMER armci_timer
#define MALLOC_LOC 1 /**< use ARMCI_Malloc_local instead of plain malloc */
#ifndef GA_ABS
#   define GA_ABS(a) ((a)>0? (a): -(a))
#endif
#define OP_GET 0
#define OP_PUT 1
#define OP_ACC 2
#define ENABLE_CLEANUP 0

static int CHECK_RESULT = 0;
static int chunk[CHUNK_NUM] = {
      1,   3,   4,   6,   9,  12,  16,  20,  24,  30,
     40,  48,  52,  64,  78,  91, 104, 128, 142, 171,
    210, 256, 300, 353, 400, 440, 476, 512 };
static int nproc = -1;
static int me = -1;
static int warn_accuracy = 0;

static void fill_array(double *arr, int count, int which);
static double time_op(int g_a, double *buf_, int chunk, int loop, int proc,
                      int ndim, int op);
static void test_1D();
static void test_2D();


int main(int argc, char **argv)
{
  /* initialize GA */
  MP_INIT(argc,argv);
  GA_Initialize_args(&argc, &argv);

  me = GA_Nodeid();
  nproc = GA_Nnodes();

  if (nproc < 2) {
    if (me == 0) {
      fprintf(stderr, "USAGE: 2 <= processes - got %d\n", nproc);
    }
    GA_Terminate();
    MP_FINALIZE();
    exit(0);
  }

  if (!me) {
    printf("\n             Performance of Basic Blocking Communication Operations\n");
  }
  GA_Sync();

  /* test 1 dimension array */
  if (!me) {
    printf("\n\t\t\tContiguous Data Transfer\n");
  }
  test_1D();

  /* test 2 dimension array */
  if (!me) {
    printf("\n\t\t\tStrided Data Transfer\n");
  }
  test_2D();

#if 0
  if (me == 0) {
    if (warn_accuracy) {
      printf("\nWARNING: Your timer does not have sufficient accuracy for this test (%d)\n", warn_accuracy);
    }
    printf("\n\n------------ Now we test the same data transfer for correctness ----------\n");
    fflush(stdout);
  }
#endif

  GA_Terminate();
  MP_FINALIZE();
  return(0);
}


double time_op(int g_a, double *buf_, int chunk, int loop, int proc,
               int ndim, int op)
{
  double start_time = 0;
  double stop_time = 0;
  double total_time = 0;
  int lo[2] = {-1,-1};
  int hi[2] = {-1,-1};
  int ld = -1;
  int i = 0;
  int bal = 0;
  double *buf = buf_;
  double alpha = 1;

  /* get the location within the g_a for the given proc */
  NGA_Distribution(g_a, proc, lo, hi);
  /* determine how much data to grab based on the chunk and dimensionality */
  if (ndim == 1) {
    hi[0] = lo[0] + chunk*chunk - 1;
  }
  else if (ndim == 2) {
    hi[0] = lo[0] + chunk - 1;
    hi[1] = lo[1] + chunk - 1;
    ld = chunk;
  }
  else {
    GA_Error("invalid ndim for time_op", ndim);
  }

  start_time = TIMER();
  for (i=0; i<loop; ++i) {
    switch (op) {
      case OP_GET:
        NGA_Get(g_a, lo, hi, buf, &ld);
        break;
      case OP_PUT:
        NGA_Put(g_a, lo, hi, buf, &ld);
        break;
      case OP_ACC:
        NGA_Acc(g_a, lo, hi, buf, &ld, &alpha);
        break;
      default:
        GA_Error("bad case value for op", op);
    }
    /* prepare next src location and dst ptr: avoid cache locality */
    if (bal == 0) {
      lo[0] += 128;
      lo[1] += 128;
      hi[0] += 128;
      hi[1] += 128;
      buf += 128;
      bal = 1;
    }
    else {
      lo[0] -= 128;
      lo[1] -= 128;
      hi[0] -= 128;
      hi[1] -= 128;
      buf -= 128;
      bal = 0;
    }
  }
  stop_time = TIMER();
  total_time = (stop_time - start_time);
  if (total_time == 0.0) {
    total_time = 0.000001; /* workaround for inaccurate timers */
    warn_accuracy++;
  }

  return(total_time / loop);
}


void test_1D()
{
  int i = 0;
  int dst = 0;
  int g_a = 0;
  int shape = SIZE * SIZE * nproc;
  int dist = SIZE * SIZE;
  int lo = 0;
  int hi = 0;
  double *buf = NULL;

  /* allocate the GA */
  g_a = NGA_Create(C_DBL, 1, &shape, "1d", &dist);
  NGA_Distribution(g_a, me, &lo, &hi);
  assert(hi-lo+1 == SIZE*SIZE);

  /* memory allocation */
  if (me == 0) {
#if MALLOC_LOC
    buf = (double *)ARMCI_Malloc_local(SIZE * SIZE * sizeof(double));
    assert(buf != NULL);
#else
    buf = (double *)malloc(SIZE * SIZE * sizeof(double));
    assert(buf != NULL);
#endif
  }

  /* only the proc 0 does the work */
  if (me == 0) {
    if (!CHECK_RESULT) {
      printf("  section               get                 put                 acc\n");
      printf("bytes   loop      usec      MB/s      usec      MB/s      usec      MB/s\n");
      printf("------- ------  --------  --------  --------  --------  --------  --------\n");
      fflush(stdout);
    }

    for (i=0; i<CHUNK_NUM; ++i) {
      int loop;
      int bytes = chunk[i] * chunk[i] * sizeof(double);
      double t_get = 0, t_put = 0, t_acc = 0;
      double latency_get, latency_put, latency_acc;
      double bandwidth_get, bandwidth_put, bandwidth_acc;

      loop = (SIZE * SIZE) / (chunk[i] * chunk[i]);
      loop = (int)sqrt((double)loop);
      if (loop < 2) {
        loop = 2;
      }

      for (dst=1; dst<nproc; ++dst) {
        /* contiguous get */
        fill_array(buf, SIZE * SIZE, me * 10);
        t_get += time_op(g_a, buf, chunk[i], loop, dst, 1, OP_GET);
        /* contiguous put */
        fill_array(buf, SIZE * SIZE, me * 10);
        t_put += time_op(g_a, buf, chunk[i], loop, dst, 1, OP_PUT);
        /* contiguous acc */
        fill_array(buf, SIZE * SIZE, me * 10);
        t_acc += time_op(g_a, buf, chunk[i], loop, dst, 1, OP_ACC);
      }

      latency_get = t_get / (nproc - 1);
      latency_put = t_put / (nproc - 1);
      latency_acc = t_acc / (nproc - 1);

      bandwidth_get = (bytes * (nproc - 1) * 1e-6) / t_get;
      bandwidth_put = (bytes * (nproc - 1) * 1e-6) / t_put;
      bandwidth_acc = (bytes * (nproc - 1) * 1e-6) / t_acc;

      /* print */
      if (!CHECK_RESULT) {
        printf("%d\t%d\t%.2e  %.2e  %.2e  %.2e  %.2e  %.2e\n",
               bytes, loop,
               latency_get / 1e-6, bandwidth_get,
               latency_put / 1e-6, bandwidth_put,
               latency_acc / 1e-6, bandwidth_acc);
      }
    }
  }
  else {
    sleep(3);
  }

  GA_Sync();

#if ENABLE_CLEANUP
  /* cleanup */
  GA_Destroy(g_a);
#endif

  if (me == 0) {
#if MALLOC_LOC
    ARMCI_Free_local(buf);
#else
    free(buf);
#endif
  }
}


void test_2D()
{
  int i = 0;
  int dst = 0;
  int g_a = 0;
  int shape[2] = {SIZE*nproc, SIZE};
  int dist[2] = {SIZE, SIZE};
  int lo[2] = {0,0};
  int hi[2] = {0,0};
  double *buf = NULL;

  /* allocate the GA */
  g_a = NGA_Create(C_DBL, 2, shape, "2d", dist);
  NGA_Distribution(g_a, me, lo, hi);
  assert(hi[0]-lo[0]+1 == SIZE);
  assert(hi[1]-lo[1]+1 == SIZE);

  /* memory allocation */
  if (me == 0) {
#if MALLOC_LOC
    buf = (double *)ARMCI_Malloc_local(SIZE * SIZE * sizeof(double));
    assert(buf != NULL);
#else
    buf = (double *)malloc(SIZE * SIZE * sizeof(double));
    assert(buf != NULL);
#endif
  }

  /* only the proc 0 does the work */
  if (me == 0) {
    if (!CHECK_RESULT) {
      printf("  section               get                 put                 acc\n");
      printf("bytes   loop      usec      MB/s      usec      MB/s      usec      MB/s\n");
      printf("------- ------  --------  --------  --------  --------  --------  --------\n");
      fflush(stdout);
    }

    for (i=0; i<CHUNK_NUM; ++i) {
      int loop;
      int bytes = chunk[i] * chunk[i] * sizeof(double);
      double t_get = 0, t_put = 0, t_acc = 0;
      double latency_get, latency_put, latency_acc;
      double bandwidth_get, bandwidth_put, bandwidth_acc;

      loop = (SIZE * SIZE) / (chunk[i] * chunk[i]);
      loop = (int)sqrt((double)loop);
      if (loop < 2) {
        loop = 2;
      }

      for (dst=1; dst<nproc; ++dst) {
        /* strided get */
        fill_array(buf, SIZE * SIZE, me * 10);
        t_get += time_op(g_a, buf, chunk[i], loop, dst, 2, OP_GET);
        /* strided put */
        fill_array(buf, SIZE * SIZE, me * 10);
        t_put += time_op(g_a, buf, chunk[i], loop, dst, 2, OP_PUT);
        /* strided acc */
        fill_array(buf, SIZE * SIZE, me * 10);
        t_acc += time_op(g_a, buf, chunk[i], loop, dst, 2, OP_ACC);
      }

      latency_get = t_get / (nproc - 1);
      latency_put = t_put / (nproc - 1);
      latency_acc = t_acc / (nproc - 1);

      bandwidth_get = (bytes * (nproc - 1) * 1e-6) / t_get;
      bandwidth_put = (bytes * (nproc - 1) * 1e-6) / t_put;
      bandwidth_acc = (bytes * (nproc - 1) * 1e-6) / t_acc;

      /* print */
      if (!CHECK_RESULT) {
        printf("%d\t%d\t%.2e  %.2e  %.2e  %.2e  %.2e  %.2e\n",
               bytes, loop,
               latency_get / 1e-6, bandwidth_get,
               latency_put / 1e-6, bandwidth_put,
               latency_acc / 1e-6, bandwidth_acc);
      }
    }
  }
  else {
    sleep(3);
  }

  GA_Sync();

#if ENABLE_CLEANUP
  /* cleanup */
  GA_Destroy(g_a);
#endif

  if (me == 0) {
#if MALLOC_LOC
    ARMCI_Free_local(buf);
#else
    free(buf);
#endif
  }
}


void fill_array(double *arr, int count, int which)
{
  int i;

  for (i = 0; i < count; i++) {
    arr[i] = i * 8.23 + which * 2.89;
  }
}
