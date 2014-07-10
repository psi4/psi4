#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <vector>

#include <assert.h>
#include <math.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include <armci.h>
#include <ga.h>
#include <gfutex.h>

#include <tbb/atomic.h>
#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_queue.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/queuing_mutex.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/tbb_thread.h>
#include <tbb/tick_count.h>

#include "cscc.h"
#include "scf.h"
#include "twoel.h"

using namespace globalFutures;

namespace globalFutures {
  extern std::ofstream ofs;
}

static const int MaxTasks = 1000;
static const int MaxExecTasks = 10;

typedef double ChunkArray[ichunk][ichunk];

struct LocalBuffers {
  ChunkArray tf_ij, td_kl, tf_ik, td_jl, ts_ij, ts_kl;
};

struct MapProcs {
  int *map, *procs;

  MapProcs() : map(NULL), procs(NULL) { };
  ~MapProcs() { if (map) delete [] map; if (procs) delete [] procs; };
};

typedef tbb::enumerable_thread_specific<double> TIAccum;
typedef tbb::enumerable_thread_specific<LocalBuffers> TLocalBuffers;
typedef tbb::enumerable_thread_specific<MapProcs> TMapProcs;

struct RTaskData {
  int ntasks;
  long int tasks[MaxTasks];

  RTaskData() : ntasks(0) { };
};

typedef std::map<int, RTaskData> RTaskMap;
typedef tbb::enumerable_thread_specific<RTaskMap> TRTaskMap;

static void TwoElTaskExecute(int g_a, int ndims, int plo[], int phi[], void *arg);
static void TwoElTaskCheck(int g_a, int ndims, int plo[], int phi[], void *arg);
static void TwoElTaskSecondCheck(int g_a, int ndims, int plo[], int phi[], void *arg);
static void TwoElTaskSend(int g_a, int ndims, int plo[], int phi[], void *arg);
static void TwoElTaskRedist(int g_a, int ndims, int plo[], int phi[], void *arg);

struct TwoElExecuteArgs {
  double schwmax;
  int ntasks;
  long int tasks[MaxExecTasks];
};

struct TwoElCheckArgs {
  double schwmax;
  long int taskid;
  int lo[4], hi[4];
};

struct TwoElSecondCheckArgs {
  int ntasks;
  long int tasks[MaxTasks];

  TwoElSecondCheckArgs() : ntasks(0) { };
};

struct TwoElSendArgs {
  double schwmax;
  int ntasks;
  long int tasks[MaxTasks];
};

struct TwoElComputeArgs {
  double schwmax;
};

struct TwoElTaskRedistArgs {
  int ntasks;
  long int tasks[MaxTasks];
};

struct RedistSched {
  int dst;
  int ntasks;
};

typedef std::vector<RedistSched> RedistPlan;

typedef tbb::concurrent_queue<long int> TaskQueue;
typedef tbb::concurrent_hash_map<long int, bool> TaskFilter;

static TIAccum tiaccum(0.0);

static TLocalBuffers localBuffs;
static TMapProcs mapProcs;

static TRTaskMap ttmap;

static int nprocs, nid, me;

static GFHandle te_hndl, tch_hndl, t2ch_hndl, ts_hndl, trd_hndl;
static tbb::atomic<long int> execTasks;
static tbb::atomic<long long int> aicut1, aicut2, aicut3;

static bool filter_done = false;
static TaskQueue *taskQueue = NULL, *ftaskQueue = NULL;
static TaskFilter tfilter;

static bool is_task_local(int g_a, int *lo, int *hi, int me, int *map, int *procs);
static void add_rtask_and_send(RTaskMap &tmap, long int taskid, int dst, double schwmax);
static void flush_rtasks(RTaskMap &tmap, double schwmax);
static void flush_sctasks();
static void print_task_stats(int me);
static void redistribute_work(int me, RedistPlan &redistPlan);
static void send_redist_futures(int me, RedistPlan &redistPlan);
static void execute_tasks(int me, double schwmax);
static void execute_tasks_nf(int me, double schwmax);

int GFInitialize()
{
  int ret = globalFutures::GFInitialize();

  te_hndl = GFRegister(TwoElTaskExecute, sizeof(TwoElExecuteArgs));
  tch_hndl = GFRegister(TwoElTaskCheck, sizeof(TwoElCheckArgs));
  t2ch_hndl = GFRegister(TwoElTaskSecondCheck, sizeof(TwoElSecondCheckArgs));
  ts_hndl = GFRegister(TwoElTaskSend, sizeof(TwoElSendArgs));
  trd_hndl = GFRegister(TwoElTaskRedist, sizeof(TwoElTaskRedistArgs));

  execTasks = 0L;

  return ret;
} // GFInitialize

void GFFinalize()
{
  if (taskQueue)
    delete taskQueue;
  if (ftaskQueue)
    delete ftaskQueue;

  globalFutures::GFFinalize();
} // GFFinalize

void twoel_orig(double schwmax, double *etwo, int nproc)
{
  double f_ij[ichunk][ichunk], d_kl[ichunk][ichunk];
  double f_ik[ichunk][ichunk], d_jl[ichunk][ichunk];
  double s_ij[ichunk][ichunk], s_kl[ichunk][ichunk];
  double one;
      
  long long int ijcnt, klcnt, ijklcnt;
  int lo[4], hi[4], lo_ik[2], hi_ik[2], lo_jl[2], hi_jl[2];
  int i, j, k, l, iloc, jloc, kloc, lloc, ich, it, jt, kt, lt;
  int dotask, newtask, accum;
  long int itask, taskcnt = 0L;

  int ctask;
  double gg;

  ga_nbhdl_t f0, f1;

  // add in the two-electron contribution to the fock matrix;

  one = 1.00;
  ijcnt = icut1;
  klcnt = icut2;
  ijklcnt = icut3;

  ich = ichunk;

  me = GA_Nodeid();

  if (!ftaskQueue || ftaskQueue->empty()) {
    double t0, t1;

    assert(!taskQueue);

    t0 = MPI_Wtime();
    taskQueue = new TaskQueue;

    taskcnt = 0L;
    GA_Zero(g_counter);
    dotask = next_4chunk(lo, hi, &it, &jt, &kt, &lt, &itask);
    ctask = 0;
    newtask = 1;
    accum = 0;
      
    while (dotask) {
      lo_ik[0] = lo[0];
      lo_ik[1] = lo[2];
      hi_ik[0] = hi[0];
      hi_ik[1] = hi[2];
      lo_jl[0] = lo[1];
      lo_jl[1] = lo[3];
      hi_jl[0] = hi[1];
      hi_jl[1] = hi[3];

      taskcnt++;

      GF_CachedGet(g_schwarz, lo, hi, s_ij, &ich);

      clean_chunk(f_ij);
      clean_chunk(f_ik);

      for (i = lo[0]; i <= hi[0]; i++) {
	iloc = i - lo[0];
	for (j = lo[1]; j <= hi[1]; j++) {
	  jloc = j - lo[1];
	  if ((s_ij[iloc][jloc] * schwmax) < tol2e)
	    icut1 = icut1 + (hi[2] - lo[2] + 1) * (hi[3] - lo[3] + 1);
	  else {
	    ga_nbhdl_t s, d0, d1;

	    if (newtask) {
	      GF_CachedNbGet(g_schwarz, &lo[2], &hi[2], s_kl, &ich, &s);
	      GF_CachedNbGet(g_dens, &lo[2], &hi[2], d_kl, &ich, &d0);
	      GF_CachedNbGet(g_dens, lo_jl, hi_jl, d_jl, &ich, &d1);

	      GF_CachedNbWait(&s);
	    }

	    for (k = lo[2]; k <= hi[2]; k++) {
	      kloc = k - lo[2];
	      for (l = lo[3]; l <= hi[3]; l++) {
		lloc = l - lo[3];
		if (s_ij[iloc][jloc] * s_kl[kloc][lloc] < tol2e)
		  icut2 = icut2 + 1;
		else {
		  if (newtask) {
		    GF_CachedNbWait(&d0);
		    GF_CachedNbWait(&d1);

		    newtask = 0;
		    execTasks++;
		  }

		  g(&gg, i, j, k, l);
		  f_ij[iloc][jloc] = f_ij[iloc][jloc] + gg * d_kl[kloc][lloc];
		  f_ik[iloc][kloc] = f_ik[iloc][kloc] - 0.50 * gg * d_jl[jloc][lloc];
		  icut3 = icut3 + 1;
		  accum = 1;
		}
	      }
	    }

	    if (newtask) {
	      GF_CachedNbWait(&d0);
	      GF_CachedNbWait(&d1);

	      newtask = 0;
	    }
	  }
	}
      }
      if (accum) {
	// GF_NbAcc(g_fock, lo, hi, f_ij, &ich, &one, &f0);
	GF_CachedAcc(g_fock, lo, hi, f_ij, &ich, &one);
	// GF_NbAcc(g_fock, lo_ik, hi_ik, f_ik, &ich, &one, &f1);
	GF_CachedAcc(g_fock, lo_ik, hi_ik, f_ik, &ich, &one);

	taskQueue->push(itask);
      }

      dotask = next_4chunk(lo, hi, &it, &jt, &kt, &lt, &itask);

      if (accum) {
	// GF_NbWait(&f0);
	// GF_NbWait(&f1);
      }

      if (dotask) {
	newtask = 1;
	accum = 0;
      }

      ctask++;
    }

    ftaskQueue = new TaskQueue(*taskQueue);

    delete taskQueue;
    taskQueue = NULL;

    t1 = MPI_Wtime();

    TIAccum::reference ltiacc = tiaccum.local();

    ltiacc += (t1 - t0);
  }
  else {
    assert(!taskQueue);

    taskQueue = new TaskQueue(*ftaskQueue);

    aicut1 = icut1;
    aicut2 = icut2;
    aicut3 = icut3;

    execute_tasks_nf(me, schwmax); // Actual execution of enqueued tasks

    delete taskQueue;
    taskQueue = NULL;

    icut1 = aicut1;
    icut2 = aicut2;
    icut3 = aicut3;
  }

#if 1
  std::ostringstream ostr;

  ostr << "tasks.dat." << me;

  std::ofstream ofsl(ostr.str().c_str(), std::ios_base::out | std::ios_base::app);

  ofsl << "Number of tasks inspected: " << taskcnt << std::endl;
  ofsl << "Executed real tasks: " << execTasks << std::endl << std::endl;

  for (TIAccum::const_iterator iter = tiaccum.begin(); iter != tiaccum.end(); iter++)
    ofsl << *iter << " ";

  ofsl << std::endl << std::endl;

  tiaccum.clear();
  execTasks = 0L;
#endif

  // Flush Fock matrix accumulate cache
  GF_CacheAccFlush(g_fock);

  *etwo = 0.50 * contract_matrices(g_fock, g_dens);
  ijcnt = icut1 - ijcnt;
  klcnt = icut2 - klcnt;
  ijklcnt = icut3 - ijklcnt;
  icut4 = icut3;

  GF_CacheReadOnlyEmpty(g_dens);

  if (icut3 > 0)
    return;

  //    no integrals may be calculated if there is no work for;
  //    this node (ichunk too big), or, something is wrong;
  printf("no two-electron integrals computed by node %d\n", me);
  printf("\n"); 
} // twoel_orig

void twoel(double schwmax, double *etwo, int nproc)
{
  int dotask, itask, ndim;
  int lo[4], hi[4], it, jt, kt, lt;
  long long int ijcnt, klcnt, ijklcnt;
  long int taskid, tidlo, tidhi, currLim, incrT, tmpi;
  double tr0, tr1, tr2, tr3, teaccum, tqaccum, tgaccum, traccum;
  double tg0, tg1, te0, te1, teqaccum;
  long long int totTasks, procTasks, ppos;

  // add in the two-electron contribution to the fock matrix;
  ijcnt = icut1;
  klcnt = icut2;
  ijklcnt = icut3;

  me = GA_Nodeid();
  nid = GA_Cluster_nodeid();

#if 0
  if (!ofs.is_open()) {
    std::ostringstream ostr;

    ostr << "diags.dat." << me;
    ofs.open(ostr.str().c_str(), std::ios_base::out | std::ios_base::app);
  }
#endif

  nprocs = nproc;

  tmpi = ceil(static_cast<double>(nbfn) / static_cast<double>(ichunk));
  incrT = tmpi * tmpi;
  totTasks = incrT * incrT;
  
  procTasks = ceil(static_cast<double>(totTasks) / static_cast<double>(nproc));

  tidlo = me * procTasks;
  tidhi = (me + 1) * procTasks - 1;

  int *map = new int[2 * 2 * nproc];
  int *procs = new int[nproc];

  RTaskMap tmap;

  itask = 0;

  TwoElCheckArgs tcargs;

  tcargs.schwmax = schwmax;
  ndim = GA_Ndim(g_schwarz);

  aicut1 = icut1;
  aicut2 = icut2;
  aicut3 = icut3;

  teaccum = 0.0;
  tqaccum = 0.0;
  tgaccum = 0.0;
  traccum = 0.0;
  teqaccum = 0.0;

  execTasks = 0L;

  if (!ftaskQueue || ftaskQueue->empty()) {
    assert(!taskQueue);

    taskQueue = new TaskQueue;

    tr0 = MPI_Wtime();

    for (taskid = tidlo; taskid <= tidhi; taskid++) {
      dotask = translate_task(taskid, lo, hi, &it, &jt, &kt, &lt);
      tcargs.taskid = taskid;

      bool is_local = is_task_local(g_schwarz, lo, hi, me, map, procs);

      if (dotask && is_local) {
	memcpy(tcargs.lo, lo, sizeof(lo));
	memcpy(tcargs.hi, hi, sizeof(hi));
	GFExecute(tch_hndl, g_schwarz, ndim, lo, hi, &tcargs);
      }
      else if (dotask && !is_local)
	add_rtask_and_send(tmap, taskid, procs[0], schwmax);

      itask++;
    }
  
    flush_rtasks(tmap, schwmax);

    tr1 = MPI_Wtime();

    tr2 = MPI_Wtime();
    GFAllQuiesce();
    flush_sctasks();
    GFAllQuiesce();
    tr3 = MPI_Wtime();

    traccum += (tr1 - tr0);

    tqaccum += (tr3 - tr2);

    RedistPlan plan;

    tg0 = MPI_Wtime();
    redistribute_work(me, plan);
    send_redist_futures(me, plan);
    GFAllQuiesce();
    tg1 = MPI_Wtime();

    tgaccum += (tg1 - tg0);

#if 1
    size_t qsz = taskQueue->unsafe_size();
    size_t *rqsz = new size_t[nproc];

    memset(rqsz, 0, sizeof(size_t) * nproc);

    MPI_Allgather(&qsz, 1, MPI_UNSIGNED_LONG, rqsz, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

    if (me == 0) {
      for (int i = 0; i < nproc; i++)
	std::cout << "rqsz[" << i << "]: " << rqsz[i] << std::endl;

      std::cout << std::flush;
    }

    delete [] rqsz;
#endif

    ftaskQueue = new TaskQueue(*taskQueue);
  }
  else {
    assert(!taskQueue);
    
    taskQueue = new TaskQueue(*ftaskQueue);
  }

  te0 = MPI_Wtime();
  execute_tasks(me, schwmax); // Actual execution of enqueued tasks
  te1 = MPI_Wtime();
  tr2 = MPI_Wtime(); 
  GFAllQuiesce();
  filter_done = true;
  tr3 = MPI_Wtime();

  teaccum += (te1 - te0);
  teqaccum += (tr3 - tr2);

#if 1
  std::ostringstream ostr;

  ostr << "tasks.dat." << me;

  std::ofstream ofsl(ostr.str().c_str(), std::ios_base::out | std::ios_base::app);

  ofsl << "Number of tasks inspected: " << itask << ", exec. time: " <<
    teaccum << ", exec. quiesce time: " << teqaccum << ", check time: " <<
    traccum << ", check quiesce time: " << tqaccum << ", redistribution time: " <<
    tgaccum << std::endl;
  ofsl << "Executed real tasks: " << execTasks << std::endl << std::endl;

  for (TIAccum::const_iterator iter = tiaccum.begin(); iter != tiaccum.end(); iter++)
    ofsl << *iter << " ";

  ofsl << std::endl << std::endl;

  tiaccum.clear();
#endif

  icut1 = aicut1;
  icut2 = aicut2;
  icut3 = aicut3;

#if 0
  ofsl << "Cache hits g_schwarz: " << GF_CacheGetHits(g_schwarz) << ", misses: " <<
    GF_CacheGetMisses(g_schwarz) << std::endl;
  ofsl << "Cache hits g_dens: " << GF_CacheGetHits(g_dens) << ", misses: " <<
    GF_CacheGetMisses(g_dens) << std::endl;
  ofsl << "Cache hits g_fock: " << GF_CacheGetHits(g_fock) << ", misses: " <<
    GF_CacheGetMisses(g_fock) << std::endl;
#endif

  tr0 = MPI_Wtime();
  // Flush Fock matrix accumulate cache
  GF_CacheAccFlush(g_fock);
  tr1 = MPI_Wtime();

  ofsl << "Acc cache flush time: " << (tr1 - tr0) << std::endl << std::endl;

  *etwo = 0.50 * contract_matrices(g_fock, g_dens);
  ijcnt = icut1 - ijcnt;
  klcnt = icut2 - klcnt;
  ijklcnt = icut3 - ijklcnt;
  icut4 = icut3;

#if 0
  print_task_stats(me);
#endif

  delete [] map;
  delete [] procs;

  delete taskQueue;
  taskQueue = NULL;

  GF_CacheReadOnlyEmpty(g_dens);

  if (icut3 > 0)
    return;

  // no integrals may be calculated if there is no work for;
  // this node (ichunk too big), or, something is wrong;

  printf("no two-electron integrals computed by node %d\n", me);
  printf("\n"); 
} // twoel

void TwoElTaskExecute(int g_a, int ndims, int plo[], int phi[], void *arg)
{
  const double one = 1.0;

  TwoElExecuteArgs *teargs;
  double gg;

  int i, j, k, l, iloc, jloc, kloc, lloc, ich, it, jt, kt, lt;
  int lo[4], hi[4], lo_ik[2], hi_ik[2], lo_jl[2], hi_jl[2];
  int newtask, accum;

  ChunkArray *f_ij, *d_kl;
  ChunkArray *f_ik, *d_jl;
  ChunkArray *s_ij, *s_kl;

  ga_nbhdl_t f0, f1;

  unsigned int lTaskPos;

  long long int laicut1 = 0L, laicut2 = 0L, laicut3 = 0L;
  long int lexecTasks = 0L;

  teargs = reinterpret_cast<TwoElExecuteArgs *>(arg);

  TLocalBuffers::reference ltbuffs = localBuffs.local();

  f_ij = &ltbuffs.tf_ij;
  d_kl = &ltbuffs.td_kl;
  f_ik = &ltbuffs.tf_ik;
  d_jl = &ltbuffs.td_jl;
  s_ij = &ltbuffs.ts_ij;
  s_kl = &ltbuffs.ts_kl;

  ich = ichunk;

  tbb::tick_count t0, t1;

  t0 = tbb::tick_count::now();
  for (int pos = 0; pos < teargs->ntasks; pos++) {
    int dotask;
    long int taskid = teargs->tasks[pos];
    ga_nbhdl_t s;

    if (filter_done) {
      TaskFilter::const_accessor cacc;

      bool found = tfilter.find(cacc, taskid);

      if (!found)
        continue;
    }

    newtask = 1;
    accum = 0;

    dotask = translate_task(taskid, lo, hi, &it, &jt, &kt, &lt);

    assert(dotask);

    lo_ik[0] = lo[0];
    lo_ik[1] = lo[2];
    hi_ik[0] = hi[0];
    hi_ik[1] = hi[2];
    lo_jl[0] = lo[1];
    lo_jl[1] = lo[3];
    hi_jl[0] = hi[1];
    hi_jl[1] = hi[3];

    GF_CachedGet(g_schwarz, lo, hi, *s_ij, &ich);
    GF_CachedNbGet(g_schwarz, &lo[2], &hi[2], *s_kl, &ich, &s);
    
    clean_chunk(*f_ij);
    clean_chunk(*f_ik);

    for (i = lo[0]; i <= hi[0]; i++) {
      iloc = i - lo[0];
      for (j = lo[1]; j <= hi[1]; j++) {
	jloc = j - lo[1];
	if (((*s_ij)[iloc][jloc] * teargs->schwmax) < tol2e)
	  laicut1 += (hi[2] - lo[2] + 1) * (hi[3] - lo[3] + 1);
	else {
	  ga_nbhdl_t d0, d1;

	  if (newtask) {
	    GF_CachedNbGet(g_dens, &lo[2], &hi[2], *d_kl, &ich, &d0);
	    GF_CachedNbGet(g_dens, lo_jl, hi_jl, *d_jl, &ich, &d1);

	    GF_CachedNbWait(&s);
	  }

	  for (k = lo[2]; k <= hi[2]; k++) {
	    kloc = k - lo[2];
	    for (l = lo[3]; l <= hi[3]; l++) {
	      lloc = l - lo[3];
	      if ((*s_ij)[iloc][jloc] * (*s_kl)[kloc][lloc] < tol2e)
		laicut2++;
	      else {
		if (newtask) {
		  GF_CachedNbWait(&d0);
		  GF_CachedNbWait(&d1);

		  newtask = 0;
		  lexecTasks++;
		}

		g(&gg, i, j, k, l);
		(*f_ij)[iloc][jloc] = (*f_ij)[iloc][jloc] + gg *
		  (*d_kl)[kloc][lloc];
		(*f_ik)[iloc][kloc] = (*f_ik)[iloc][kloc] - 0.50 * gg *
		  (*d_jl)[jloc][lloc];
		laicut3++;
		accum = 1;
	      }
	    }
	  }

	  if (newtask) {
	    GF_CachedNbWait(&d0);
	    GF_CachedNbWait(&d1);

	    newtask = 0;
	  }
	}
      }
    }

    if (accum) {
      // GF_NbAcc(g_fock, lo, hi, *f_ij, &ich, const_cast<double *>(&one), &f0);
      // GF_NbAcc(g_fock, lo_ik, hi_ik, *f_ik, &ich, const_cast<double *>(&one), &f1);

      GF_CachedAcc(g_fock, lo, hi, *f_ij, &ich, const_cast<double *>(&one));
      GF_CachedAcc(g_fock, lo_ik, hi_ik, *f_ik, &ich, const_cast<double *>(&one));

      // GF_NbWait(&f0);
      // GF_NbWait(&f1);

      if (!filter_done) {
        TaskFilter::const_accessor cacc;

        tfilter.insert(cacc, taskid);
      }
    }
  }

  aicut1 += laicut1;
  aicut2 += laicut2;
  aicut3 += laicut3;

  execTasks += lexecTasks;

  t1 = tbb::tick_count::now();

  TIAccum::reference ltiacc = tiaccum.local();

  ltiacc += (t1 - t0).seconds();
} // TwoElTaskExecute

void TwoElTaskCheck(int g_a, int ndims, int plo[], int phi[], void *arg)
{
  TwoElCheckArgs *tcargs = reinterpret_cast<TwoElCheckArgs *>(arg);
  bool rte = false;
  ChunkArray *s_ij;
  double *rs_ij, *ors_ij;
  int ich = ichunk;
  int accs0;
  int *map, *procs;
  bool hget0 = false;

  TLocalBuffers::reference ltbuffs = localBuffs.local();
  TMapProcs::reference ltmapProcs = mapProcs.local();

  if (!ltmapProcs.map)
    ltmapProcs.map = new int[2 * 2 * nprocs];
  if (!ltmapProcs.procs)
    ltmapProcs.procs = new int[nprocs];

  map = ltmapProcs.map;
  procs = ltmapProcs.procs;

  int np = GF_Locate_region(g_schwarz, tcargs->lo, tcargs->hi, map, procs);

  if (np > 1 || procs[0] != me)
    hget0 = true;

  if (hget0) {
    s_ij = &ltbuffs.ts_ij;
    GF_CachedGet(g_schwarz, tcargs->lo, tcargs->hi, *s_ij, &ich);
    rs_ij = &(*s_ij)[0][0];
    ors_ij = rs_ij;
    accs0 = ich;
  }
  else {
    GF_Access(g_schwarz, tcargs->lo, tcargs->hi, &rs_ij, &accs0);
    ors_ij = rs_ij;
  }

  for (int il = tcargs->lo[0]; il <= tcargs->hi[0] && !rte; il++) {
    int iloc = il - tcargs->lo[0];

    rs_ij = (ors_ij + iloc * accs0);

    for (int jl = tcargs->lo[1]; jl <= tcargs->hi[1] && !rte; jl++) {
      int jloc = jl - tcargs->lo[1];

      if ((*rs_ij * tcargs->schwmax) >= tol2e)
	rte = true;

      rs_ij++;
    }
  }

  if (rte) {
#if 1
    TRTaskMap::reference lttmap = ttmap.local();

    np = GF_Locate_region(g_schwarz, &tcargs->lo[2], &tcargs->hi[2], map, procs);

    RTaskData &tdata = lttmap[procs[0]];

    tdata.tasks[tdata.ntasks] = tcargs->taskid;
    tdata.ntasks++;

    if (tdata.ntasks == MaxTasks) {
      TwoElSecondCheckArgs tscargs;

      tscargs.ntasks = tdata.ntasks;
      memcpy(tscargs.tasks, tdata.tasks, sizeof(tdata.tasks[0]) * tdata.ntasks);

      GFExecute(t2ch_hndl, g_proc, 1, &procs[0], &procs[0], &tscargs);

      tdata.ntasks = 0;
    }
#else
    taskQueue->push(tcargs->taskid);
#endif
  }

  if (!hget0)
    GF_Release(g_schwarz, tcargs->lo, tcargs->hi);
} // TwoElTaskCheck

void TwoElTaskSecondCheck(int g_a, int ndims, int plo[], int phi[], void *arg)
{
  TwoElSecondCheckArgs *tscargs = reinterpret_cast<TwoElSecondCheckArgs *>(arg);
  ChunkArray *s_kl;
  double *rs_kl, *ors_kl;
  int lo[4], hi[4], it, jt, kt, lt;
  long int taskid;
  int accs1;
  int ich = ichunk;
  int *map, *procs;

  TLocalBuffers::reference ltbuffs = localBuffs.local();
  TMapProcs::reference ltmapProcs = mapProcs.local();

  if (!ltmapProcs.map)
    ltmapProcs.map = new int[2 * 2 * nprocs];
  if (!ltmapProcs.procs)
    ltmapProcs.procs = new int[nprocs];

  map = ltmapProcs.map;
  procs = ltmapProcs.procs;

  for (int i = 0; i < tscargs->ntasks; i++) {
    bool rti = false;
    bool hget1 = false;

    taskid = tscargs->tasks[i];

    bool dotask = translate_task(taskid, lo, hi, &it, &jt, &kt, &lt);

    int np = GF_Locate_region(g_schwarz, &lo[2], &hi[2], map, procs);

    if (np > 1 || procs[0] != me)
      hget1 = true;

    if (hget1) {
      s_kl = &ltbuffs.ts_kl;
      GF_CachedGet(g_schwarz, &lo[2], &hi[2], *s_kl, &ich);
      rs_kl = &(*s_kl)[0][0];
      ors_kl = rs_kl;
      accs1 = ich;
    }
    else {
      GF_Access(g_schwarz, &lo[2], &hi[2], &rs_kl, &accs1);
      ors_kl = rs_kl;
    }

    for (int kl = lo[2]; kl <= hi[2] && !rti; kl++) {
      int kloc = kl - lo[2];

      rs_kl = (ors_kl + kloc * accs1);

      for (int ll = lo[3]; ll <= hi[3] && !rti; ll++) {
	int lloc = ll - lo[3];

	if (fabs(*rs_kl) >= tol2e)
	  rti = true;

	rs_kl++;
      }
    }

    if (rti)
      taskQueue->push(taskid);

    if (!hget1)
      GF_Release(g_schwarz, &lo[2], &hi[2]);
  }
} // TwoElTaskSecondCheck

void TwoElTaskSend(int g_a, int ndims, int plo[], int phi[], void *arg)
{
  TwoElSendArgs *tsargs = reinterpret_cast<TwoElSendArgs *>(arg);
  int lo[4], hi[4], it, jt, kt, lt, ndim;
  bool dotask;
  long int taskid;
  TwoElCheckArgs tcargs;

  tcargs.schwmax = tsargs->schwmax;
  ndim = GF_Ndim(g_schwarz);

  for (int i = 0; i < tsargs->ntasks; i++) {
    taskid = tsargs->tasks[i];
    dotask = translate_task(taskid, lo, hi, &it, &jt, &kt, &lt);
    tcargs.taskid = taskid;

    memcpy(tcargs.lo, lo, sizeof(lo));
    memcpy(tcargs.hi, hi, sizeof(hi));

    GFExecute(tch_hndl, g_schwarz, ndim, lo, hi, &tcargs);
  }
} // TwoElTaskSend

void TwoElTaskRedist(int g_a, int ndims, int plo[], int phi[], void *arg)
{
  TwoElTaskRedistArgs *trargs = reinterpret_cast<TwoElTaskRedistArgs *>(arg);
  const int ntasks = trargs->ntasks;

  for (int i = 0; i < ntasks; i++)
    taskQueue->push(trargs->tasks[i]);
} // TwoElTaskRedist

void print_task_stats(int me)
{
  std::ostringstream ostr;

  ostr << "tasks.dat." << me;

  std::ofstream ofs(ostr.str().c_str(), std::ios_base::out | std::ios_base::app);

  ofs << "Executed real tasks: " << execTasks << std::endl << std::endl;

  for (TIAccum::const_iterator iter = tiaccum.begin(); iter != tiaccum.end(); iter++)
    ofs << *iter << " ";

  ofs << std::endl << std::endl;
} // print_task_stats

bool is_task_local(int g_a, int *lo, int *hi, int me, int *map, int *procs)
{
  int np = GF_Locate_region(g_a, lo, hi, map, procs);

  return (procs[0] == me);
} // is_task_local

void add_rtask_and_send(RTaskMap &tmap, long int taskid, int dst, double schwmax)
{
  RTaskData &tdata = tmap[dst];

  tdata.tasks[tdata.ntasks] = taskid;
  tdata.ntasks++;

  if (tdata.ntasks == MaxTasks) {
    TwoElSendArgs tsargs;

    tsargs.schwmax = schwmax;
    memcpy(tsargs.tasks, tdata.tasks, sizeof(tdata.tasks[0]) * tdata.ntasks);
    tsargs.ntasks = tdata.ntasks;

    GFExecute(ts_hndl, g_proc, 1, &dst, &dst, &tsargs);

    tdata.ntasks = 0;
  }
} // add_rtask_and_send

void flush_rtasks(RTaskMap &tmap, double schwmax)
{
  TwoElSendArgs tsargs;

  tsargs.schwmax = schwmax;

  for (RTaskMap::iterator iter = tmap.begin(); iter != tmap.end(); iter++) {
    int dst = iter->first;

    memcpy(tsargs.tasks, iter->second.tasks,
	   sizeof(iter->second.tasks[0]) * iter->second.ntasks);
    tsargs.ntasks = iter->second.ntasks;

    GFExecute(ts_hndl, g_proc, 1, &dst, &dst, &tsargs);

    iter->second.ntasks = 0;
  }
} // flush_rtasks

void flush_sctasks()
{
  TwoElSecondCheckArgs tscargs;

  for (TRTaskMap::const_iterator iter = ttmap.begin(); iter != ttmap.end(); iter++)
    for (RTaskMap::const_iterator iiter = iter->begin(); iiter != iter->end(); iiter++) {
      int dst = iiter->first;

      memcpy(tscargs.tasks, iiter->second.tasks,
	     sizeof(iiter->second.tasks[0]) * iiter->second.ntasks);
      tscargs.ntasks = iiter->second.ntasks;

      GFExecute(t2ch_hndl, g_proc, 1, &dst, &dst, &tscargs);
    }

  ttmap.clear();
} // flush_sctasks

void redistribute_work(int me, RedistPlan &redistPlan)
{
  size_t qsz = taskQueue->unsafe_size();
  size_t *rqsz = new size_t[nprocs] /*, *crqsz = new size_t[nprocs] */;

  memset(rqsz, 0, sizeof(size_t) * nprocs);

  MPI_Allgather(&qsz, 1, MPI_UNSIGNED_LONG, rqsz, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

#if 1
  if (me == 0) {
    for (int i = 0; i < nprocs; i++)
      std::cout << "rqsz[" << i << "]: " << rqsz[i] << std::endl;

    std::cout << std::flush;
  }
#endif

  size_t tqsz = 0U, avgt, cavgt;

  for (int i = 0; i < nprocs; i++)
    tqsz += rqsz[i];

  // Integer floor
  avgt = tqsz / nprocs;
  cavgt = (tqsz + nprocs - 1) / nprocs;

#if 1
  if (me == 0) {
    std::cout << "avgt: " << avgt << ", cavgt: " << cavgt << std::flush << std::endl;
    std::cout << "tqsz: " << tqsz << std::endl;

    for (int i = 0; i < nprocs; i++)
      if (rqsz[i] > avgt)
	std::cout << "over avgt: " << i << ", " << rqsz[i] - avgt << std::endl;
  }
#endif

  // memcpy(crqsz, rqsz, sizeof(size_t) * nprocs);

  if (rqsz[me] < avgt)
    goto lend;

  for (int lme = 0; lme < nprocs; lme++) {
    for (int i = 0; i < nprocs; i++) {
      size_t tdd;

      if (i == me || i == lme)
	continue;

      if (rqsz[me] <= cavgt)
	break;

      if (rqsz[lme] <= cavgt)
	break;

      if (rqsz[i] < cavgt) {
	tdd = std::min(rqsz[lme] - cavgt, cavgt - rqsz[i]);

	if (tdd > 0) {
	  rqsz[lme] -= tdd;
	  rqsz[i] += tdd;
	}

	if (me == lme && tdd > 0) {
	  RedistSched sch;

	  sch.dst = i;
	  sch.ntasks = tdd;
	  redistPlan.push_back(sch);
	}
      }
    }
  }

 lend:
#if 0
  int *irqsz = new int[nprocs], *rrqsz = new int[nprocs];

  memset(irqsz, 0, sizeof(int) * nprocs);
  memset(rrqsz, 0, sizeof(int) * nprocs);

  for (RedistPlan::const_iterator iter = redistPlan.begin(); iter != redistPlan.end(); iter++) {
    irqsz[me] -= iter->ntasks;
    irqsz[iter->dst] += iter->ntasks;
  }

  MPI_Reduce(irqsz, rrqsz, nprocs, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (me == 0) {
    for (int i = 0; i < nprocs; i++)
      crqsz[i] += rrqsz[i];
  }

  delete [] rrqsz;
  delete [] irqsz;

  std::ostringstream ostr;

  ostr << "tasks.dat." << me;

  std::ofstream ofs(ostr.str().c_str(), std::ios_base::out | std::ios_base::app);

  if (me == 0) {
    size_t sqsz = 0U;

    for (int i = 0; i < nprocs; i++) {
      sqsz += crqsz[i];

      ofs << "crqsz[" << i << "]: " << crqsz[i] << std::endl;
    }

    ofs << "Total tasks: " << sqsz << std::endl;

    ofs << std::endl << std::endl;
  }
 
  for (RedistPlan::const_iterator iter = redistPlan.begin(); iter != redistPlan.end(); iter++)
    ofs << "Send: " << iter->ntasks << ", to proc: " << iter->dst << std::endl;
#endif

  delete [] rqsz;
  // delete [] crqsz;  
} // redistribute_work

void send_redist_futures(int me, RedistPlan &redistPlan)
{
  for (RedistPlan::iterator iter = redistPlan.begin(); iter != redistPlan.end(); iter++) {
    long int taskid;
    bool popped;
    int outer = (iter->ntasks + MaxTasks - 1) / MaxTasks;
    int lastBatch = iter->ntasks - (outer - 1) * MaxTasks;

    if (lastBatch < 0)
      lastBatch = iter->ntasks;

    for (int nit = 0; nit < outer; nit++) {
      TwoElTaskRedistArgs trargs;
      int pos = 0;
      int tasksToPop = (nit == (outer - 1)) ? lastBatch : MaxTasks;

      trargs.ntasks = tasksToPop;
      for (int i = 0; i < tasksToPop; i++) {
	popped = taskQueue->try_pop(taskid);

	assert(popped);
	trargs.tasks[pos] = taskid;
	pos++;
      }

      GFExecute(trd_hndl, g_proc, 1, &iter->dst, &iter->dst, &trargs);
    }
  }
} // send_redist_futures

void execute_tasks(int me, double schwmax)
{
  long int taskid;
  bool popped;
  int qsz = taskQueue->unsafe_size();

  int outer = (qsz + MaxExecTasks - 1) / MaxExecTasks;
  int lastBatch = qsz - (outer - 1) * MaxExecTasks;

  if (lastBatch < 0)
    lastBatch = qsz;

  for (int nit = 0; nit < outer; nit++) {
    TwoElExecuteArgs teargs;
    int pos = 0;
    int tasksToPop = (nit == (outer - 1)) ? lastBatch : MaxExecTasks;

    teargs.schwmax = schwmax;
    teargs.ntasks = tasksToPop;
    for (int i = 0; i < tasksToPop; i++) {
      popped = taskQueue->try_pop(taskid);

      assert(popped);
      teargs.tasks[pos] = taskid;
      pos++;
    }

    GFExecute(te_hndl, g_proc, 1, &me, &me, &teargs);
  }
} // execute_tasks

void execute_tasks_nf(int me, double schwmax)
{
  long int taskid;
  bool popped;
  int qsz = taskQueue->unsafe_size();

  int outer = (qsz + MaxExecTasks - 1) / MaxExecTasks;
  int lastBatch = qsz - (outer - 1) * MaxExecTasks;

  if (lastBatch < 0)
    lastBatch = qsz;

  for (int nit = 0; nit < outer; nit++) {
    TwoElExecuteArgs teargs;
    int pos = 0;
    int tasksToPop = (nit == (outer - 1)) ? lastBatch : MaxExecTasks;

    teargs.schwmax = schwmax;
    teargs.ntasks = tasksToPop;
    for (int i = 0; i < tasksToPop; i++) {
      popped = taskQueue->try_pop(taskid);

      assert(popped);
      teargs.tasks[pos] = taskid;
      pos++;
    }

    TwoElTaskExecute(g_proc, 1, &me, &me, &teargs);
  }
} // execute_tasks_nf
