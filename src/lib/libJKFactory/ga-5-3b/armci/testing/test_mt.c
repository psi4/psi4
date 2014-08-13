#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: test_mt.c,v 1.1.2.1 2007-07-02 05:34:13 d3p687 Exp $
 * test_mt.c
 *
 * Developed by Andriy Kot <andriy.kot@pnl.gov>
 * Copyright (c) 2007 Pacific Northwest National Laboratory
 *
 * Changelog:
 * 2007-02-17 - created
 *
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
#if HAVE_ASSERT_H
#   include <assert.h>
#endif

#include "armci.h"
#include "message.h"
#include "utils.h"

#define DEBUG /* note: a TBUFSIZE (per thread) buffer is used to print arrays */
#define TBUFSIZE 65535
/* prints debug information to files named test_mt.th<#th_rank>*/
/*#define LOG2FILE*/
#define NOTHREADS_ /* debug: does not spawn threads if set */

typedef double atype_t; /* type of the element of the array */
#define MAX_TPP 8 /* max treads per processor */
#define TPP     1//2 /* threads per processor */
#define ASIZE   4//5000 /* size of array */
#define ITERS   1//10 /* iterations per test */
enum {PUT = 101, GET, ACC};

int tpp = TPP;
int asize = ASIZE;
int iters = ITERS;
int delay = 0;
double scale = 2.0;

#define THREAD_OFF 1000.0
#define ITER_OFF 10.0
#define ELEM_INC 0.0003

/* each ARMCI mem block that allocated consists of
 * th_size blocks, one per each thread (systemwide), which consist of
 * iters blocks, one per each iteration, which consist of
 * asize doubles (or atype_t)
 */
#define ASIZE_BYTES (asize*sizeof(atype_t))
#define ASIZExITERS (asize*iters)
#define ASIZExITERS_BYTES (ASIZE_BYTES*iters)
#define ASIZExITERSxTH (asize*iters*th_size)
#define ASIZExITERSxTH_BYTES (ASIZE_BYTES*iters*th_size)

/* p_   - pointer (points to area owned by particula thread)
 * th_  - thread section
 * it_  - iteration
 * i_   - element offset (0 for beginning of array)
 */
#define AELEM(p_,th_,it_,i_) ((atype_t *)p_)[th_*ASIZExITERS+it_*asize+i_]
#define AELEM_VAL(th_,it_,i_) (THREAD_OFF*th_+ITER_OFF*it_+ELEM_INC*(i_+1))

int rank, size, th_size;
int th_rank[MAX_TPP];
unsigned time_seed;
int *pairs, *rnd_tgts, rnd_one;
void **ptrs1, **ptrs2;

FILE *dbg[MAX_TPP];
char fname[] = "test_mt.th000";
char cbuf[2048];
int cbufl;

#define PRINTF0 if(!rank)printf
#define PRINTF0T if(!TH_ME)printf
#define RND(l_,h_) (l_+(int)(((double)h_)*rand()/(RAND_MAX+((double)l_))))
#define TH_ME (th_rank[th_idx])

void prndbg(int th_idx, char *fmt, ...)
{
#ifdef DEBUG
  va_list ap;
  va_start(ap, fmt);
#ifdef LOG2FILE
#   define DSCR dbg[th_idx]
#else
#   define DSCR stdout
  printf("%3d: ", TH_ME);
#endif
  vfprintf(DSCR, fmt, ap);
  fflush(DSCR);
  va_end(ap);
#endif
}

void usage()
{
  if (!rank) {
    printf("Usage: test_mt, or \n");
    printf("       test_mt -tTHREADS_PER_PROC -sARRAY_SIZE -iITERATIONS_COUNT\n");
  }
  ARMCI_Barrier();
  armci_msg_finalize();
  exit(0);
}

void *thread_main(void *arg);
void zero_array(int th_idx, void *ptr);
void init_array(int th_idx, void *ptr);
void print_array(int th_idx, char *msg, atype_t *array);
void test_pairs(int th_idx); // deprecated?
void test_PutGetAcc(int th_idx, int tgt, int *rmt, int rmt_cnt);
void check_PutGetAcc(int th_idx, int rmt, int op, atype_t *array);

int main(int argc, char *argv[])
{
  int ch;
  extern char *optarg;
  int i, j, r;
  thread_t threads[MAX_TPP];

  /* init ARMCI */
  armci_msg_init(&argc, &argv);
  ARMCI_Init_args(&argc, &argv);
  size = armci_msg_nproc();
  rank = armci_msg_me();

  while ((ch = getopt(argc, argv, "t:s:i:d:h")) != -1) {
    switch (ch) {
      case 't': /* # of threads */
        tpp = atoi(optarg);
        if (tpp < 1 || tpp > MAX_TPP) {
          PRINTF0("\"%s\" is improper value for -t, should be a "
                  "number between 1 and %d(MAX_TPP)\n",
                  optarg, MAX_TPP);
          usage();
        }
        break;
      case 'i': /* # of iterations */
        iters = atoi(optarg);
        if (iters < 1) {
          PRINTF0("\"%s\" is improper value for -t, should be a "
                  "number equal or larger than 1\n", optarg);
          usage();
        }
        break;
      case 's': /* # of elements in the array */
        asize = atoi(optarg);
        if (iters < 1) {
          PRINTF0("\"%s\" is improper value for -s, should be a "
                  "number equal or larger than 1\n", optarg);
          usage();
        }
        break;
      case 'd':
        delay = atoi(optarg);
        break; /* delay before start */
      case 'h':
        usage();
        break; /* print usage info */
    }
  }
#ifdef NOTHREADS
  tpp = 1;
  PRINTF0("Warning: NOTHREADS debug symbol is set -- running w/o threads\n");
#endif
  th_size = size * tpp;
  PRINTF0("\nTest of multi-threaded capabilities:\n"
          "%d threads per process (%d threads total),\n"
          "%d array elements of size %d,\n"
          "%d iteration(s)\n\n", tpp, th_size, asize, sizeof(atype_t), iters);
  if (delay) {
    printf("%d: %d\n", rank, getpid());
    fflush(stdout);
    sleep(delay);
    ARMCI_Barrier();
  }
  TH_INIT(size, tpp);
  for (i = 0; i < tpp; i++) {
    th_rank[i] = rank * tpp + i;
  }

#if defined(DEBUG) && defined(LOG2FILE)
  for (i = 0; i < tpp; i++) {
    fname[10] = '0' + th_rank[i] / 100;
    fname[11] = '0' + th_rank[i] % 100 / 10;
    fname[12] = '0' + th_rank[i] % 10;
    dbg[i] = fopen(fname, "w");
  }
#endif
  for (i = 0; i < tpp; i++) {
    prndbg(i, "proc %d, thread %d(%d):\n", rank, i, th_rank[i]);
  }

  /* set global seed (to ensure same random sequence across procs) */
  time_seed = (unsigned)time(NULL);
  armci_msg_brdcst(&time_seed, sizeof(time_seed), 0);
  srand(time_seed);
  rand();
  prndbg(0, "seed = %u\n", time_seed);
  /* random pairs */
  pairs = calloc(th_size, sizeof(int));
  for (i = 0; i < th_size; i++) {
    pairs[i] = -1;
  }
  for (i = 0; i < th_size; i++) {
    if (pairs[i] != -1) {
      continue;
    }
    r = RND(0, th_size);
    while (i == r || pairs[r] != -1) {
      r = RND(0, th_size);
    }
    pairs[i] = r;
    pairs[r] = i;
  }
  for (i = 0, cbufl = 0; i < th_size; i++)
    cbufl += sprintf(cbuf + cbufl, " %d->%d|%d->%d",
                     i, pairs[i], pairs[i], pairs[pairs[i]]);
  prndbg(0, "random pairs:%s\n", cbuf);
  /* random targets */
  rnd_tgts = calloc(th_size, sizeof(int));
  for (i = 0, cbufl = 0; i < th_size; i++) {
    rnd_tgts[i] = RND(0, th_size);
    if (rnd_tgts[i] == i) {
      i--;
      continue;
    }
    cbufl += sprintf(cbuf + cbufl, " %d", rnd_tgts[i]);
  }
  prndbg(0, "random targets:%s\n", cbuf);
  /* random one */
  rnd_one = RND(0, th_size);
  prndbg(0, "random one = %d\n", rnd_one);

  assert(ptrs1 = calloc(th_size, sizeof(void *)));
  assert(ptrs2 = calloc(th_size, sizeof(void *)));
#ifdef NOTHREADS
  thread_main((void *)(long)0);
#else
  for (i = 0; i < tpp; i++) {
    THREAD_CREATE(threads + i, thread_main, (void *)(long)i);
  }
  for (i = 0; i < tpp; i++) {
    THREAD_JOIN(threads[i], NULL);
  }
#endif

  ARMCI_Barrier();
  PRINTF0("Tests Completed\n");

  /* clean up */
#if defined(DEBUG) && defined(LOG2FILE)
  for (i = 0; i < tpp; i++) {
    fclose(dbg[i]);
  }
#endif
  ARMCI_Finalize();
  TH_FINALIZE();
  armci_msg_finalize();

  return 0;
}


void *thread_main(void *arg)
{
  int th_idx, i;
  int tgt, *rmt, rmt_cnt;

  th_idx = (int)(long)arg;
  prndbg(th_idx, "thread %d(%d|%d) STARTED\n", TH_ME, rank, th_idx);

  assert(!ARMCI_MALLOC_MT(ptrs1, ASIZExITERSxTH_BYTES));
  assert(!ARMCI_MALLOC_MT(ptrs2, ASIZExITERSxTH_BYTES));
#if 0
  for (i = 0, cbufl = 0; i < th_size; i++) {
    cbufl += sprintf(cbuf + cbufl, " %p", ptrs1[i]);
  }
  prndbg(th_idx, "ptrs1: %s\n", cbuf);
  for (i = 0, cbufl = 0; i < th_size; i++) {
    cbufl += sprintf(cbuf + cbufl, " %p", ptrs2[i]);
  }
  prndbg(th_idx, "ptrs2: %s\n", cbuf);
#endif
#if 0
  init_array(th_idx, ptrs1[TH_ME]);
  init_array(th_idx, ptrs2[TH_ME]);
#endif

  assert(rmt = calloc(th_size, sizeof(int)));

  PRINTF0T("  TESTING GET/PUT/ACC\n\n");

  /* test pairs */
  PRINTF0T("> Testing pair-wise communication pattern ...\n");
  tgt = rmt[0] = pairs[TH_ME];
  rmt_cnt = 1;
  test_PutGetAcc(th_idx, tgt, rmt, rmt_cnt);
  PRINTF0T("  pair-wise is OK\n\n");
  //return NULL; // REMOVE WHEN DONE

  /* test random target */
  PRINTF0T("> Testing random target communication pattern ...\n");
  tgt = rnd_tgts[TH_ME];
  for (i = 0, rmt_cnt = 0; i < th_size; i++) if (rnd_tgts[i] == TH_ME) {
      rmt[rmt_cnt++] = i;
    }
  test_PutGetAcc(th_idx, tgt, rmt, rmt_cnt);
  PRINTF0T("  random target is OK\n\n");

  /* test all to one */
  PRINTF0T("> Testing hotspot(all to one) communication pattern ...\n");
  if (TH_ME == rnd_one) {
    tgt = -1;
    for (i = 0, rmt_cnt = 0; i < th_size; i++) if (i != TH_ME) {
        rmt[rmt_cnt++] = i;
      }
  }
  else {
    tgt = rnd_one;
    rmt_cnt = 0;
  }
  test_PutGetAcc(th_idx, tgt, rmt, rmt_cnt);
  PRINTF0T("  hotspot is OK\n\n");

  free(rmt);
}

void zero_array(int th_idx, void *ptr)
{
  int i, j, k;
  for (i = 0; i < th_size; i++)for (j = 0; j < iters; j++)for (k = 0; k < asize; k++) {
        AELEM(ptr, i, j, k) = 0.0;
      }
}

void init_array(int th_idx, void *ptr)
{
  int i, j, k;
  for (i = 0; i < th_size; i++)for (j = 0; j < iters; j++)for (k = 0; k < asize; k++) {
        AELEM(ptr, i, j, k) = AELEM_VAL(TH_ME, j, k);
      }
  /*AELEM(ptr, i, j) = THREAD_OFF*TH_ME + ITER_OFF*i + ELEM_INC*(j+1);*/

  print_array(th_idx, "initialized", ptr);
#if 0
#   if 1
  for (i = 0, cbufl = 0; i < th_size; i++) {
    for (j = 0; j < iters; j++) {
      cbufl += sprintf(cbuf + cbufl, "(%d,%d)%p:", i, j, &(((atype_t *)ptr)[i*ASIZExITERS+j*asize]));
      for (k = 0; k < asize; k++) {
        cbufl += sprintf(cbuf + cbufl, " %.4f", ((atype_t *)ptr)[i*ASIZExITERS+j*asize+k]);
      }
      cbufl += sprintf(cbuf + cbufl, "\n");
    }
    cbufl += sprintf(cbuf + cbufl, "\n");
  }
#   else
  for (i = 0, cbufl = 0; i < (th_size * iters * asize); i++) {
    cbufl += sprintf(cbuf + cbufl, " %.4f", ((atype_t *)ptr)[i]);
  }
#   endif
  prndbg(th_idx, "initialized:\n%s\n", cbuf);
#endif
}

void print_array(int th_idx, char *msg, atype_t *array)
{
#ifdef DEBUG
  int i, j, k, tbufl;
  char tbuf[TBUFSIZE];

  if (ASIZExITERSxTH_BYTES > TBUFSIZE / 2) {
    prndbg(th_idx, "%s:\n%s\n", msg, "array is too big to print");
  }

  for (i = 0, tbufl = 0; i < th_size; i++) {
    for (j = 0; j < iters; j++) {
      tbufl += sprintf(tbuf + tbufl, "(%d,%d)%p:", i, j, &AELEM(array, i, j, 0));
      for (k = 0; k < asize; k++) {
        tbufl += sprintf(tbuf + tbufl, " %.4f", AELEM(array, i, j, k));
      }
      tbufl += sprintf(tbuf + tbufl, "\n");
    }
    tbufl += sprintf(tbuf + tbufl, "\n");
  }
  prndbg(th_idx, "%s:\n%s\n", msg, tbuf);
#endif
}


/*
void print_array(int th_idx, char *msg, atype_t *array, int count) {
#ifdef DEBUG
    int i;
    for (i = 0, cbufl = 0; i < count; i++)
        cbufl+=sprintf(cbuf+cbufl, " %.4f", array[i]);
    prndbg(th_idx, "%s:%s\n", msg, cbuf);
#endif
}
*/

int check_result(atype_t *array, int th)
{
  int i, j, k, mismatch;
  for (i = 0, k = 0, mismatch; i < iters; i++) for (j = 0; j < asize; j++, k++) {
      if (array[k] != AELEM_VAL(th, i, j)) {
        printf("mismatch detected: th=%d, i=%d, j=%d, k=%d, elem=%d, array=%d\n",
               th, i, j, k, AELEM_VAL(th, i, j), array[k]);
        fflush(stdout);
        abort();
      }
    }
}

void test_pairs(int th_idx)
{
  int rem_th, rem_proc;
  int i, j;
  void *src, *dst;

  rem_th = pairs[TH_ME];
  rem_proc = TH2PROC(rem_th);

  prndbg(th_idx, "test_pair: %d<->%d(%d)\n", TH_ME, rem_th, rem_proc);

  MT_BARRIER();
#if 0
  print_array(th_idx, "before", &AELEM(ptrs2[TH_ME], rem_th, 0, 0), ASIZExITERS);
#endif
  for (i = 0; i < iters; i++) {
    /* src - addr of my thread block on remote proc/thread */
    src = &AELEM(ptrs1[rem_th], TH_ME, i, 0);
    /* src - addr of remote thread block on my proc/thread */
    dst = &AELEM(ptrs2[TH_ME], rem_th, i, 0);
    /* get from my pair */
    assert(!ARMCI_Get(src, dst, ASIZE_BYTES, rem_proc));
  }

  MT_BARRIER();
#if 0
  print_array(th_idx, "rcvd", &AELEM(ptrs2[TH_ME], rem_th, 0, 0), ASIZExITERS);
#endif
  /* check results */
  check_result(&AELEM(ptrs2[TH_ME], rem_th, 0, 0), rem_th);

}

/* test Put/Get/Acc sequence regardless of communication pattern
 *  tgt -- remote target for put/get/acc (none if -1)
 *  rmt -- list of remote thread that put/acc to here (correctness is cheked here)
 *  rmt_cnt -- # of threads in rmt
 */
void test_PutGetAcc(int th_idx, int tgt, int *rmt, int rmt_cnt)
{
  /* a - local thread, b - remote thread */
  int a, b, b_proc, stride[2], count[2];
  int i, j;
  void *src, *dst;
#ifdef DEBUG
  for (i = 0, cbufl = 0; i < rmt_cnt; i++) {
    cbufl += sprintf(cbuf + cbufl, " %d", rmt[i]);
  }
  prndbg(th_idx, "test_PutGetAcc: put/acc to %d, get from %d, check put/acc from %s\n",
         tgt, tgt, rmt_cnt ? cbuf : "none");
#endif
  a = TH_ME;
  stride[0] = ASIZE_BYTES;
  count[0] = ASIZE_BYTES;
  count[1] = 1;

  /* init arrays */
  init_array(th_idx, ptrs1[TH_ME]);
  init_array(th_idx, ptrs2[TH_ME]);
  MT_BARRIER();

  /* put - put a.ptrs1[b] into b.ptrs2[a] */
  if (tgt != -1) {
    b = tgt;
    b_proc = TH2PROC(b);
    for (i = 0; i < iters; i++) {
      src = &AELEM(ptrs1[a], b, i, 0); /* a.ptrs1[b] */
      dst = &AELEM(ptrs2[b], a, i, 0); /* b.ptrs2[a] */
      //            assert(!ARMCI_Put(src, dst, ASIZE_BYTES, b_proc));
      assert(!ARMCI_PutS(src, stride, dst, stride, count, 1, b_proc));
    }
    ARMCI_Fence(b_proc);
  }
  MT_BARRIER();
  print_array(th_idx, "PUT:ptrs1[TH_ME]", ptrs1[TH_ME]);
  print_array(th_idx, "PUT:ptrs2[TH_ME]", ptrs2[TH_ME]);
  MT_BARRIER();

  /* chk put(s) from b(s): a.ptrs2[b] */
  for (j = 0; j < rmt_cnt; j++) {
    b = rmt[j];
    b_proc = TH2PROC(b);
    check_PutGetAcc(th_idx, b, PUT, &AELEM(ptrs2[a], b, 0, 0));
  }
  //return; // REMOVE WHEN DONE

  /* init arrays */
  init_array(th_idx, ptrs1[TH_ME]);
  init_array(th_idx, ptrs2[TH_ME]);
  MT_BARRIER();

  /* get - get b.ptrs1[a] into a.ptrs2[b] */
  if (tgt != -1) {
    b = tgt;
    b_proc = TH2PROC(b);
    for (i = 0; i < iters; i++) {
      src = &AELEM(ptrs1[b], a, i, 0); /* b.ptrs1[a] */
      dst = &AELEM(ptrs2[a], b, i, 0); /* a.ptrs2[b] */
      assert(!ARMCI_GetS(src, stride, dst, stride, count, 1, b_proc));
    }
  }
  print_array(th_idx, "GET:ptrs1[TH_ME]", ptrs1[TH_ME]);
  print_array(th_idx, "GET:ptrs2[TH_ME]", ptrs2[TH_ME]);
  MT_BARRIER();

  /* chk get from b: a.ptrs2[b] */
  if (tgt != -1) {
    check_PutGetAcc(th_idx, b, GET, &AELEM(ptrs2[a], b, 0, 0));
  }

#if 1
  /* init arrays */
  init_array(th_idx, ptrs1[TH_ME]);
  init_array(th_idx, ptrs2[TH_ME]);
  MT_BARRIER();

  /* acc - acc a.ptrs1[b] * scale + b.ptrs2[a] into b.ptrs2[a] */
  if (tgt != -1) {
    b = tgt;
    b_proc = TH2PROC(b);
    for (i = 0; i < iters; i++) {
      src = &AELEM(ptrs1[a], b, i, 0); /* a.ptrs1[b] */
      dst = &AELEM(ptrs2[b], a, i, 0); /* b.ptrs2[a] */
      assert(!ARMCI_AccS(ARMCI_ACC_DBL, &scale, src, stride, dst, stride, count, 1, b_proc));
    }
    ARMCI_Fence(b_proc);
  }
  MT_BARRIER();
  print_array(th_idx, "ACC:ptrs1[TH_ME]", ptrs1[TH_ME]);
  print_array(th_idx, "ACC:ptrs2[TH_ME]", ptrs2[TH_ME]);
  MT_BARRIER();

  /* chk acc(s) from b(s): a.ptrs2[b] */
  for (j = 0; j < rmt_cnt; j++) {
    b = rmt[j];
    b_proc = TH2PROC(b);
    check_PutGetAcc(th_idx, b, ACC, &AELEM(ptrs2[a], b, 0, 0));
  }

#endif
  MT_BARRIER();
}

void check_PutGetAcc(int th_idx, int rmt, int op, atype_t *array)
{
  int i, j, k;
  double expected;

  for (i = 0, k = 0; i < iters; i++) for (j = 0; j < asize; j++, k++) {
      expected = op == ACC ? AELEM_VAL(TH_ME, i, j) + scale * AELEM_VAL(rmt, i, j) :
                 AELEM_VAL(rmt, i, j);
      if (array[k] != expected) {
        printf("mismatch detected: TM_ME=%d, rmt=%d, op=%d, i=%d, j=%d, "
               "k=%d, expected=%f, array=%f\n",
               TH_ME, rmt, op, i, j, k, expected, array[k]);
        fflush(stdout);
        sleep(5);
        abort();
      }
#if 0
      if (array[k] != AELEM_VAL(th, i, j)) {
        printf("mismatch detected: th=%d, i=%d, j=%d, k=%d, elem=%d, array=%d\n",
               th, i, j, k, AELEM_VAL(th, i, j), array[k]);
        fflush(stdout);
        sleep(10);
        abort();
      }
#endif
    }
}


