#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <sched.h>
#include <stdint.h>
#include <unistd.h>

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <list>
#include <stdexcept>
#include <sstream>
#include <utility>

#include <mpi.h>

#include <ga.h>
#include <gpc.h>

#include <tbb/spin_mutex.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/tbb_thread.h>
#include <tbb/tick_count.h>

#include <boost/interprocess/sync/scoped_lock.hpp>

#include "gfutex.h"
#include "gfutex_internal.h"

namespace globalFutures_implementation {
  int nproc, me;
  int nodeme;
}

namespace globalFutures {
  using namespace globalFutures_implementation;

  typedef int (*Gpc_Func)();

  static const char *const msgq_name = "futures_queue";
  static const char *const tbb_nthreads = "TBB_NTHREADS";

  static const int maxCnt = INT_MAX / 100;
  static const int maxRets = 10000000;

  static const int msgqMaxNumMsgs = 51200;
  static const int msgqMaxMsgSz = 8192;

  static const unsigned int complete = 1U;
  static const unsigned int incomplete = 0U;

  static int msgLenOff;
  static int msgDatOff;

  static message_queue *msgq = NULL;
  static NodeQueues msgqs;
  static ProcIdMap pmap;

  static int nthreads;
  static int hgpced, hgpcxd;
  static atomic<unsigned int> listen;
  static GAMutex inGA; // Is any thread in the GA library?
  static FuncArray funcArray;
  static unsigned int *pretid_cnt;
  static char tmpdirname[] = "gfutex_XXXXXX";
  static MPI_Comm pcomm;
  static unsigned int volatile **retptrs;
  static unsigned int currRetPos;
  static RetMutex currRetLock;
  static DirectMap directActivity;

  static task_scheduler_init *tsched = NULL;

  static atomic<unsigned int> pbcnt, execnt, msgscnt, msgrcnt;

  static void send_gpc(GFHandle hndl, int g_a, int ndims, int lo[], int hi[], void *arg,
		       bool enqueue);
  static int gpc_enqueue_dispatcher(int to, int from, void *hdr, int hlen, void *data, int dlen,
				    void *rhdr, int rhlen, int *rhsize, void *rdata, int rdlen,
				    int *rdsize, int rtype);
  static int gpc_execute_dispatcher(int to, int from, void *hdr, int hlen, void *data, int dlen,
				    void *rhdr, int rhlen, int *rhsize, void *rdata, int rdlen,
				    int *rdsize, int rtype);
  static void msgq_listener();
  static void wait_remote(GFHandle hndl, int np, int procs[], unsigned int rdata[]);
  static void process_execmsg();
  static int put_value_long(unsigned int val, int offs, int target);
  static void print_debug_info();

  int GFInitialize()
  {
    char *retdir, *snthreads;
    int lendir;

    MPI_Comm_dup(MPI_COMM_WORLD, &pcomm); // Duplicate MPI_COMM_WORLD for compatibility

    me = GA_Nodeid();
    nproc = GA_Nnodes();

    snthreads = getenv(tbb_nthreads);

    if (snthreads)
      nthreads = atoi(snthreads) + 1;
    else
      nthreads = task_scheduler_init::default_num_threads();

    if (me == 0)
      std::cerr << "Number of TBB threads to use: " << nthreads << std::endl;

    int nnds = GA_Cluster_nnodes();

    if (nnds > 1)
      tsched = new task_scheduler_init(task_scheduler_init::deferred);
    else
      tsched = new task_scheduler_init(nthreads);

    hgpced = ARMCI_Gpc_register(reinterpret_cast<Gpc_Func>(gpc_enqueue_dispatcher));
    hgpcxd = ARMCI_Gpc_register(reinterpret_cast<Gpc_Func>(gpc_execute_dispatcher));

    retptrs = new unsigned int volatile *[nproc];

    ARMCI_Malloc(reinterpret_cast<void **>(const_cast<unsigned int **>(retptrs)),
		 sizeof(unsigned int) * maxRets);

    int nid = GA_Cluster_nodeid();
    int nodenproc = GA_Cluster_nprocs(nid);

    for (int i = 0; i < maxRets; i++)
      retptrs[me][i] = incomplete;

    for (int i = 0; i < nodenproc; i++)
      if (me == GA_Cluster_procid(nid, i)) {
	nodeme = i;
	break;
      }

    try {
      std::ostringstream ostr;
      std::string qname;

      if (nodeme == 0) {
	for (int i = 0; i < nodenproc; i++) {
	  int rpid = GA_Cluster_procid(nid, i);

	  ostr << msgq_name << rpid;

	  qname = ostr.str();

	  message_queue::remove(qname.c_str());

	  message_queue *nqueue = new message_queue(create_only, qname.c_str(),
						    msgqMaxNumMsgs, msgqMaxMsgSz);

	  QueueData qdata;

	  qdata.queue = nqueue;
	  qdata.qname = qname;

	  qdata.mname = qname + "_mutex";

	  named_mutex::remove(qdata.mname.c_str());

	  qdata.mutex = new named_mutex(create_only, qdata.mname.c_str());

	  msgqs.push_back(qdata);

	  ostr.str("");

	  pmap[rpid] = i;
	}

	msgq = msgqs.at(0).queue;

	MPI_Barrier(pcomm);
      }
      else {
	MPI_Barrier(pcomm);

	for (int i = 0; i < nodenproc; i++) {
	  int rpid = GA_Cluster_procid(nid, i);

	  ostr << msgq_name << rpid;

	  qname = ostr.str();

	  message_queue *nqueue = new message_queue(open_only, qname.c_str());

	  QueueData qdata;

	  qdata.queue = nqueue;
	  qdata.qname = qname;

	  qdata.mname = qname + "_mutex";

	  qdata.mutex = new named_mutex(open_only, qdata.mname.c_str());

	  msgqs.push_back(qdata);

	  if (rpid == me)
	    msgq = nqueue;

	  ostr.str("");

	  pmap[rpid] = i;
	}
      }
    }
    catch (const MsgqException &excp) {
      std::cerr << __LINE__ << " " << excp.what() << std::endl;
      exit(EXIT_FAILURE);
    }
    catch(interprocess_exception &ex){
      std::cerr << ex.what() << ": " << __LINE__ << std::endl;
      exit(EXIT_FAILURE);
    }

    msgLenOff = 1;
    msgDatOff = nproc + 1;

    listen = 1U;
    currRetPos = 0U;
    pbcnt = 0U;
    execnt = 0U;
    msgscnt = 0U;
    msgrcnt = 0U;

    MPI_Barrier(pcomm);

    return 1;
  } // GFInitialize

  void GFFinalize()
  {
    int retc;

    MPI_Barrier(pcomm);

    try {
      for (FuncArray::iterator iter = funcArray.begin(); iter != funcArray.end(); iter++) {
	RemoteFunction *func = *iter;

	if (func)
	  func->finalizeThreads();
      }

      listen = 0U;

      assert(tsched);

      delete tsched;

      for (NodeQueues::iterator iter = msgqs.begin(); iter != msgqs.end(); iter++) {
	message_queue *q = iter->queue;
	named_mutex *m = iter->mutex;

	assert(q);
	delete q;

	assert(m);
	delete m;

	if (nodeme == 0) {
	  message_queue::remove(iter->qname.c_str());
	  named_mutex::remove(iter->mname.c_str());
	}

	iter->queue = NULL;
	iter->mutex = NULL;
      }
    }
    catch (const MsgqException &excp) {
      std::cerr << __LINE__ << " " << excp.what() << std::endl;
      exit(EXIT_FAILURE);
    }

    MPI_Barrier(pcomm);

    for (FuncArray::iterator iter = funcArray.begin(); iter != funcArray.end(); iter++) {
      RemoteFunction *func = *iter;

      if (func)
	delete func;
    }

    ARMCI_Free(const_cast<unsigned int *>(retptrs[me]));

    delete [] retptrs;

    MPI_Barrier(pcomm);

    MPI_Comm_free(&pcomm);
  } // GFFinalize

  GFHandle GFRegister(RemoteFuncProto func, const size_t argSz, const size_t retSz)
  {
    if (argSz > msgqMaxMsgSz) {
      std::cerr << "GFRegister: Argument size exceeds maximum allowed!, argSz: " << argSz <<
	", max. allowed: " << msgqMaxMsgSz << std::endl;
      exit(EXIT_FAILURE);
    }

    RemoteFunction *ptr = new RemoteFunction(func, argSz, retSz);

    funcArray.push_back(ptr); // Elements are owned by the array now

    MPI_Barrier(pcomm);

    return funcArray.size() - 1;
  } // GFRegister

  void send_gpc(GFHandle hndl, int g_a, int ndims, int lo[], int hi[], void *arg,
		bool enqueue)
  {
    try {
      int *plo, *phi;
      void *parg;
      int np;
      int hgpc;

      RemoteFunction *ptr = funcArray.at(hndl);

      int *map = new int[2 * ndims * nproc];
      int *procs = new int[nproc];
      size_t tlen = sizeof(GPCArgsExec) + sizeof(lo[0]) * ndims + sizeof(hi[0]) * ndims +
	ptr->getArgSz();
      char *sndBuf = new char[tlen];
      GPCArgsExec *gpcArgs = reinterpret_cast<GPCArgsExec *>(sndBuf);
      bool direct;

      if (enqueue)
	hgpc = hgpced;
      else
	hgpc = hgpcxd;

      plo = reinterpret_cast<int *>(gpcArgs->dat);
      phi = plo + ndims;
      parg = phi + ndims;

      gpcArgs->g_a = g_a;
      gpcArgs->ndims = ndims;
      gpcArgs->src = me;

      direct = (task::self().parent() == NULL);

      memcpy(plo, lo, sizeof(lo[0]) * ndims);
      memcpy(phi, hi, sizeof(hi[0]) * ndims);
      memcpy(parg, arg, ptr->getArgSz());

      np = GF_Locate_region(g_a, lo, hi, map, procs);

      unsigned int *rdata = new unsigned int[np];

      for (int i = 0; i < np; i++) {
	GAMutex::scoped_lock glock;
	RetMutex::scoped_lock rlock;
	unsigned int retPos;
	unsigned int retDat;	
	
	rlock.acquire(currRetLock); // Lock position counter
	retPos = currRetPos++;

	if (currRetPos == maxRets)
	  currRetPos = 0U;

	rlock.release(); // Unlock position counter

	if (retptrs[me][retPos] == complete) {
	  std::cerr << "Error executing GPC for process: " << i <<
	    ", from process: " << me << ", too many concurrent GPCs ?" << std::endl;
	  std::cerr << "directActivity.size(): " << directActivity.size() << std::endl;
	  exit(EXIT_FAILURE);
	}

	gpcArgs->offs = retPos;
	// Address to indicate completion
	rdata[i] = retPos;

	{
	  DirectMap::accessor a;

	  directActivity.insert(a, retPos);
	  a->second = direct;
	}

	glock.acquire(inGA);
	ARMCI_Gpc_exec(hgpc, procs[i], &hndl, sizeof(hndl), gpcArgs, tlen, NULL, 0,
		       &retDat, sizeof(retDat), NULL);
	glock.release();
      }

      delete [] rdata;

      delete [] map;
      delete [] procs;
      delete [] sndBuf;
    }
    catch (const std::out_of_range &excp) {
      std::cerr << "Proc: " << me << ", function not registered: " << hndl <<
	", line:" << __LINE__ << ", " << excp.what() << std::endl;
      std::cerr << "funcArray.size(): " << funcArray.size() << std::endl;
      exit(EXIT_FAILURE);
    }
  } // send_gpc

  void GFEnqueue(GFHandle hndl, int g_a, int ndims, int lo[], int hi[], void *arg)
  {
    send_gpc(hndl, g_a, ndims, lo, hi, arg, true);
  } // GFEnqueue

  void GFExecute(GFHandle hndl, int g_a, int ndims, int lo[], int hi[], void *arg)
  {
    send_gpc(hndl, g_a, ndims, lo, hi, arg, false);
  } // GFExecute

  int GFMaxConcurrency()
  {
    return maxRets;
  } // GFMaxConcurrency

  void GFQuiesce(GFHandle hndl)
  {
    try {
      if (task::self().parent() != NULL)
	throw GFutException("GFQuiesce() must be called only from the main process!");

      bool allquiet;

      // Now wait for all activities spawned on other processes to complete
      do {
	allquiet = true;

	for (DirectMap::const_iterator iter = directActivity.begin();
	     iter != directActivity.end(); iter++) {	  
	  unsigned int pos = iter->first;

	  if (iter->second && retptrs[me][pos] != complete)
	    allquiet = false;
	}
      } while(!allquiet);

      // All directly spawned activities are quiesced for this process

      for (DirectMap::iterator iter = directActivity.begin(); iter != directActivity.end();
	   iter++) {
	unsigned int pos = iter->first;

	if (iter->second) {
	  assert(retptrs[me][pos] == complete);
	  retptrs[me][pos] = incomplete;

	  directActivity.erase(pos);
	}
      }
    }
    catch (const std::out_of_range &excp) {
      std::cerr << "Function not registered: " << __LINE__ << excp.what() << std::endl;
      exit(EXIT_FAILURE);
    }
    catch (const GFutException &excp) {
      std::cerr << excp.what() << std::endl;
      exit(EXIT_FAILURE);
    }
  } // GFQuiesce

  void GFAllQuiesce(GFHandle hndl)
  {
    try {
      if (task::self().parent() != NULL)
	throw GFutException("GFQuiesce() must be called only from the main process!");

      RemoteFunction *ptr = funcArray.at(hndl);

      int compcnt;
      int totcnt;
      int passcnt = 0;
      int maxRefCnt = 0;

#if 0
      tbb::tick_count t0, t1;
      double taccum = 0.0, caccum = 0.0, waccum = 0.0, raccum = 0.0;
#endif

      // Now wait for all activities spawned locally & on other processes to complete
#if 0
      t0 = tbb::tick_count::now();
#endif

      do {
	compcnt = 0;
	totcnt = 0;

#if 0
	tbb::tick_count t2, t3;

	t2 = tbb::tick_count::now();
#endif
	for (DirectMap::const_iterator iter = directActivity.begin();
	     iter != directActivity.end(); iter++) {
	  unsigned int pos = iter->first;

	  if (/* iter->second && */ retptrs[me][pos] != complete)
	    compcnt++;
	}

	msgq_listener();

#if 0
	t3 = tbb::tick_count::now();
	caccum += (t3 - t2).seconds();
#endif

	task *parentTask = ptr->getParentTask();
	int  refCnt = -1;

	if (!parentTask)
	  goto reduction;

#if 0
	t2 = tbb::tick_count::now();
#endif
	refCnt = parentTask->ref_count();

	if (refCnt > maxRefCnt)
	  maxRefCnt = refCnt;

	if (refCnt > 0) {
	  ptr->lockParent();

	  parentTask->wait_for_all();

	  parentTask->destroy(*parentTask);

	  parentTask = new (task::allocate_root()) empty_task;
	  parentTask->set_ref_count(1);

	  ptr->setParentTask(parentTask);
	  ptr->unlockParent();
	}

#if 0
	t3 = tbb::tick_count::now();
	waccum += (t3 - t2).seconds();
#endif
	
  reduction:
#if 0
        if (passcnt % 1000000 == 0 && me == 0) {
          std::cout << "Proc: " << me << ", passcnt: " << passcnt << ", refCnt: " << refCnt <<
            ", compcnt: " << compcnt << std::flush << std::endl;
          print_debug_info();
          std::cout << "Proc: " << me << ", directActivity.size(): " << directActivity.size() <<
            std::flush << std::endl;
        }
#endif

#if 0
	t2 = tbb::tick_count::now();
#endif

	MPI_Allreduce(&compcnt, &totcnt, 1, MPI_INT, MPI_SUM, pcomm);
#if 0
	t3 = tbb::tick_count::now();
	raccum += (t3 - t2).seconds();
#endif

	passcnt++;
      } while (totcnt > 0 /* && passcnt < maxCnt */);

#if 0
      t1 = tbb::tick_count::now();
      taccum += (t1 - t0).seconds();
#endif

      MPI_Barrier(pcomm);

#if 0
      if (me == 0) {
          for (int i = 0; i < directActivity.size(); i++)
            std::cout << "Proc: " << me << ", retptrs[" << me << "][" << i << "]: " <<
              (void *)&retptrs[me][i] << std::flush << std::endl;
        std::cout << "Proc: " << me << ", final passcnt: " << passcnt <<
          ", maxRefCnt: " << maxRefCnt << ", compcnt: " << compcnt <<
          ", totcnt: " << totcnt << std::flush << std::endl;
        print_debug_info();
        std::cout << "Proc: " << me << ", final directActivity.size(): " <<
          directActivity.size() << std::flush << std::endl;
      }
#endif

#if 0
      std::cout << "Proc: " << me << ", time spent checking: " << std::fixed <<
	taccum << ", passcnt: " << passcnt << std::endl;
      std::cout << "Proc: " << me << ", time spent in loop: " << caccum <<
	", time spent waiting for tasks: " << waccum << ", time spent in reduction: " <<
	raccum << std::flush << std::endl;
#endif

      // All spawned activities are quiesced

      for (DirectMap::iterator iter = directActivity.begin();
	   iter != directActivity.end(); iter++) {
	unsigned int pos = iter->first;

	if (retptrs[me][pos] != complete) {
	  std::cerr << "Proc: " << me << ", activity: " << pos << ", should be complete!" <<
	    std::flush << std::endl;
          std::cerr << "Proc: " << me << ", pos: " << pos << ", retptrs: " <<
            retptrs[me][pos] << ", direct: " << iter->second << std::flush << std::endl;
        }

	retptrs[me][pos] = incomplete;
      }
      directActivity.clear();

      // All indirectly spawned activities are quiesced at this point
      // Safe to call GA without lock
      MPI_Barrier(pcomm);
    }
    catch (const std::out_of_range &excp) {
      std::cerr << "Function not registered: " << __LINE__ << excp.what() << std::endl;
      exit(EXIT_FAILURE);
    }
  } // GFAllQuiesce

  void GFAllQuiesce()
  {
    for (GFHandle hndl = 0U; hndl < funcArray.size(); hndl++)
      GFAllQuiesce(hndl);
  } // GFAllQuiesce

  void GF_Access(int g_a, int lo[], int hi[], void *ptr, int ld[])
  {
    GAMutex::scoped_lock lock;

    lock.acquire(inGA);

    ::NGA_Access(g_a, lo, hi, ptr, ld);

    lock.release();
  } // GF_Access

  void GF_Release(int g_a, int lo[], int hi[])
  {
    GAMutex::scoped_lock lock;

    lock.acquire(inGA);

    ::NGA_Release(g_a, lo, hi);

    lock.release();
  } // GF_Release

  void GF_Get(int g_a, int lo[], int hi[], void *buf, int ld[])
  {
    GAMutex::scoped_lock lock;

    lock.acquire(inGA);

    ::NGA_Get(g_a, lo, hi, buf, ld);

    lock.release();
  } // GF_Get

  int GF_Locate_region(int g_a, int lo[], int hi[], int map[], int procs[])
  {
    GAMutex::scoped_lock lock;

    lock.acquire(inGA);

    int np = ::NGA_Locate_region(g_a, lo, hi, map, procs);

    lock.release();

    return np;
  } // GF_Locate_region

  void GF_Distribution(int g_a, int iproc, int lo[], int hi[])
  {
    GAMutex::scoped_lock lock;

    lock.acquire(inGA);

    ::NGA_Distribution(g_a, iproc, lo, hi);

    lock.release();
  } // GF_Distribution

  void GF_Acc(int g_a, int lo[], int hi[], void *buf, int ld[], void *alpha)
  {
    GAMutex::scoped_lock lock;

    lock.acquire(inGA);

    ::NGA_Acc(g_a, lo, hi, buf, ld, alpha);

    lock.release();
  } // GFAcc

  void GF_NbAcc(int g_a, int lo[], int hi[], void *buf, int ld[], void *alpha,
		ga_nbhdl_t *nbhandle)
  {
    GAMutex::scoped_lock lock;

    lock.acquire(inGA);

    ::NGA_NbAcc(g_a, lo, hi, buf, ld, alpha, nbhandle);

    lock.release();
  } // GF_NbAcc

  void GF_NbGet(int g_a, int lo[], int hi[], void *buf, int ld[], ga_nbhdl_t *nbhandle)
  {
    GAMutex::scoped_lock lock;

    lock.acquire(inGA);

    ::NGA_NbGet(g_a, lo, hi, buf, ld, nbhandle);

    lock.release();
  } // GF_NbGet

  void GF_NbPut(int g_a, int lo[], int hi[], void *buf, int ld[], ga_nbhdl_t *nbhandle)
  {
    GAMutex::scoped_lock lock;

    lock.acquire(inGA);

    ::NGA_NbPut(g_a, lo, hi, buf, ld, nbhandle);

    lock.release();
  } // GF_NbPut

  void GF_NbWait(ga_nbhdl_t *nbhandle)
  {
    GAMutex::scoped_lock lock;

    lock.acquire(inGA);

    ::NGA_NbWait(nbhandle);

    lock.release();
  } // GF_NbWait

  int GF_Ndim(int g_a)
  {
    GAMutex::scoped_lock lock;

    lock.acquire(inGA);

    int ret = ::GA_Ndim(g_a);

    lock.release();

    return ret;
  } // GF_Ndim

  int GF_Cluster_nodeid()
  {
    GAMutex::scoped_lock lock;

    lock.acquire(inGA);

    int ret = ::GA_Cluster_nodeid();

    lock.release();

    return ret;
  } // GF_Cluster_nodeid

  int GF_Cluster_proc_nodeid(int proc)
  {
    GAMutex::scoped_lock lock;

    lock.acquire(inGA);

    int ret = ::GA_Cluster_proc_nodeid(proc);

    lock.release();

    return ret;
  } // GF_Cluster_proc_nodeid

  void GF_Inquire(int g_a, int *type, int *ndim, int dims[])
  {
    GAMutex::scoped_lock lock;

    lock.acquire(inGA);

    ::NGA_Inquire(g_a, type, ndim, dims);

    lock.release();
  } // GF_Inquire

  // Internal member functions

  RemoteFunction::~RemoteFunction()
  {
  } // ~RemoteFunction
    
  void RemoteFunction::execute(void *buf, size_t bufSz)
  {
    // Very careful allocation
    char *rbuf = new char[bufSz + sizeof(ExecutorArgs) - sizeof(GPCArgsExec)];
    // Will be deleted as the user or library checks for completion
    ExecutorArgs *rargs = reinterpret_cast<ExecutorArgs *>(rbuf);

    rargs->func = func;
    memcpy(&rargs->gpcArgs, buf, bufSz);

    if (!tsched->is_active())
      tsched->initialize(nthreads);

    if (!parentTask) {
      lockParent();
      if (!parentTask) {
	parentTask = new (task::allocate_root()) empty_task;
	parentTask->set_ref_count(1);
      }
      unlockParent();
    }

    GFutTask *rec = new (parentTask->allocate_child()) GFutTask(rbuf);

    parentTask->increment_ref_count();
    pbcnt++;
#if 1
    task::enqueue(*rec);
#else
    parentTask->spawn(*rec);
#endif
  } // execute

  void RemoteFunction::finalizeThreads(const bool terminate)
  {
    if (!parentTask)
      return;

    if (parentTask->ref_count() > 0)
      parentTask->wait_for_all(); // Wait for all child tasks to finish

    parentTask->destroy(*parentTask);

    parentTask = NULL;
  } // finalizeThreads

  // Static function definitions

  int gpc_enqueue_dispatcher(int to, int from, void *hdr, int hlen, void *data, int dlen,
			     void *rhdr, int rhlen, int *rhsize,
			     void *rdata, int rdlen, int *rdsize, int rtype)
  {
    try {
      gfutmsglen msglen;
      char *alc = new char[sizeof(gfutmsgsnd) + hlen + dlen];
      gfutmsgsnd *msgbuf = reinterpret_cast<gfutmsgsnd *>(alc);

      msglen.len = hlen + dlen;

      memcpy(msgbuf->dat, hdr, hlen);
      memcpy(msgbuf->dat + hlen, data, dlen);

      ProcIdMap::const_iterator iter = pmap.find(to);

      if (iter == pmap.end())
	throw std::out_of_range("key not found in pmap!");

      QueueData &qdata = msgqs.at(iter->second);
      message_queue *mqueue = qdata.queue;

      assert(mqueue);

      scoped_lock<named_mutex> lock(*qdata.mutex);

      bool ret = mqueue->try_send(&msglen, sizeof(msglen), to);

      if (!ret)
	throw GFutException("Message queue is full for process!");
     
      ret = mqueue->try_send(msgbuf, hlen + dlen, to);

      if (!ret)
	throw GFutException("Message queue is full for process!");

      lock.unlock();

      msgscnt++;

      delete [] alc;
    }
    catch (const MsgqException &excp) {
      std::cerr << __LINE__ << " " << excp.what() << std::endl;
      exit(EXIT_FAILURE);
    }
    catch (const GFutException &excp) {
      std::cerr << excp.what() << std::endl;
      exit(EXIT_FAILURE);
    }
    catch(interprocess_exception &ex){
      std::cerr << ex.what() << ": " << __LINE__ << std::endl;
      exit(EXIT_FAILURE);
    }
    catch (const std::out_of_range &excp) {
      std::cerr << excp.what() << ": " << __LINE__ << std::endl;
      std::cerr << "Proc: " << me << ", to: " << to << ", pmap[to]: " << pmap[to] << std::endl;
      std::cerr << "Proc: " << me << ", msgqs.size(): " << msgqs.size() << std::endl;
      exit(EXIT_FAILURE);
    }

    return GPC_DONE;
  } // gpc_enqueue_dispatcher

  int gpc_execute_dispatcher(int to, int from, void *hdr, int hlen, void *data, int dlen,
			     void *rhdr, int rhlen, int *rhsize,
			     void *rdata, int rdlen, int *rdsize, int rtype)
  {
    try {
      GFHandle hndl;

      hndl = *(reinterpret_cast<GFHandle *>(hdr));

      RemoteFunction *func = funcArray.at(hndl);

      func->execute(data, dlen);
    }

    catch (const std::out_of_range &excp) {
      std::cerr << excp.what() << ": " << __LINE__ << std::endl;
      exit(EXIT_FAILURE);
    }

    return GPC_DONE;
  } // gpc_execute_dispatcher

  void msgq_listener()
  {
    unsigned int cont_listen;
    int cnt = 0;

    do {
      try {
	process_execmsg();

	cont_listen = listen;
	cnt++;
      }
      catch (const MsgqException &excp) {
	std::cerr << __LINE__ << " " << excp.what() << std::endl;
	exit(EXIT_FAILURE);
      }
      catch (const std::out_of_range &excp) {
	std::cerr << "Function not registered: " << __LINE__ << excp.what() << std::endl;
	exit(EXIT_FAILURE);
      }
      catch(interprocess_exception &ex){
	std::cerr << ex.what() << ": " << __LINE__ << std::endl;
	std::cerr << ex.get_native_error() << ", " << ex.get_error_code();
	std::cerr << ", " << size_error << std::endl;
	exit(EXIT_FAILURE);
      }
    } while ((msgq->get_num_msg() > 0) /* && (cnt < maxCnt) */);

  } // msgq_listener

  void GFutTask::msgq_executor(void *arg)
  {
    ExecutorArgs *earg;
    int *lo, *hi;
    void *rarg;

    earg = reinterpret_cast<ExecutorArgs *>(arg);

    lo = reinterpret_cast<int *>(earg->gpcArgs.dat);
    hi = lo + earg->gpcArgs.ndims;
    rarg = hi + earg->gpcArgs.ndims;

    earg->func(earg->gpcArgs.g_a, earg->gpcArgs.ndims, lo, hi, rarg);
    put_value_long(complete, earg->gpcArgs.offs, earg->gpcArgs.src);
    execnt++;
  } // msgq_executor

  void wait_remote(GFHandle hndl, int np, int procs[], unsigned int rdata[])
  {
    bool *status = new bool[np];
    bool flag;
    int cnt = 0;

    for (int i = 0; i < np; i++)
      status[i] = false;

    do {
      flag = true;

      for (int i = 0; i < np; i++) {
	int pos = rdata[i];

	if (!status[i])
	  status[i] = (retptrs[me][pos] == complete);

	this_tbb_thread::yield();
	flag = (flag && status[i]);
      }
      cnt++;
    } while (!flag && cnt < maxCnt);

    delete [] status;

    if (cnt >= maxCnt)
      throw GFutException("wait_remote() exceeded iteration count, response did not arrive");
  } // wait_remote

  void process_execmsg()
  {
    gfutmsglen msglen;
    size_t rsz = 0U;
    unsigned int prio = 0U;

    bool ret = msgq->try_receive(&msglen, msgqMaxMsgSz, rsz, prio);

    if (ret) {
      assert(rsz == sizeof(msglen));

      char *alc = new char[sizeof(gfutmsgsnd) + msglen.len]; // Length of function type is
      // included in the overall message length
      gfutmsgrcv *msgbuf = reinterpret_cast<gfutmsgrcv *>(alc);
      bool oret;

      do {
	oret = msgq->try_receive(msgbuf, msgqMaxMsgSz, rsz, prio);
      } while (!oret);

      assert(rsz == msglen.len);

      msgrcnt++;

      RemoteFunction *func = funcArray.at(msgbuf->func);

      func->execute(&msgbuf->gpcArgs, rsz - sizeof(msgbuf->func));

      delete [] alc;
    }
  } // process_execmsg

  int put_value_long(unsigned int val, int offs, int target)
  {
    GAMutex::scoped_lock lock;
    int ret = 0;

    lock.acquire(inGA);

    ret = ARMCI_PutValueInt(val, const_cast<unsigned int *>(&retptrs[target][offs]), target);

    ARMCI_Fence(target);

    lock.release();

    return ret;
  } // put_value_long

  void print_debug_info()
  {
    std::cout << "pbcnt: " << pbcnt << ", execnt: " << execnt << std::flush << std::endl;

    std::cout << "msgscnt: " << msgscnt << ", msgrcnt: " << msgrcnt << std::endl;
    std::cout << "number of messages in queue: " << msgq->get_num_msg() <<
      std::endl;

    std::cout << "listen: " << listen << std::flush << std::endl;
  } // print_debug_info
}
