#ifndef __GFUTEX_INTERNAL_H
#define __GFUTEX_INTERNAL_H

#include <exception>
#include <map>
#include <string>
#include <vector>

#include <tbb/concurrent_hash_map.h>
#include <tbb/queuing_mutex.h>
#include <tbb/task.h>

#include <boost/interprocess/ipc/message_queue.hpp>
#include <boost/interprocess/sync/named_mutex.hpp>

using namespace tbb;
using namespace boost::interprocess;

namespace globalFutures_implementation {
  extern int nproc, me;
  extern int nodeme;
}

namespace globalFutures {
  struct GPCArgsExec {
    int g_a;
    int ndims;
    int src;
    int offs; // Remote location to indicate completion of the execution
    char dat[];
  };

  struct gfutmsglen {
    int len;
  };

  struct gfutmsgsnd {
    char dat[];
  };

  struct gfutmsgrcv {
    GFHandle func; // Coming from GPC header
    GPCArgsExec gpcArgs; // Embedded structure with the rest of the data
  };

  struct ExecutorArgs {
    RemoteFuncProto func;
    int shep;
    GPCArgsExec gpcArgs;
  };

  class MsgqException : public std::exception
  {
  private:
    int errnum;

  public:
    MsgqException(int oerrnum) : errnum(oerrnum) { }

    virtual const char *what() const throw()
    {
      return strerror(errnum);
    }
  };

  class GFutException : public std::exception
  {
  private:
    const std::string str;
  public:
    GFutException(const std::string ostr) : str(ostr) { }
    virtual ~GFutException() throw() { }

    virtual const char *what() const throw()
    {
      return str.c_str();
    }
  };

  class RemoteFunction {
  public:
    RemoteFunction(RemoteFuncProto ofunc, const size_t oargSz, const size_t oretSz = 0U)
      : argSz(oargSz), retSz(oretSz), parentTask(NULL), func(ofunc) { }
    virtual ~RemoteFunction();
    
    virtual void execute(void *buf, size_t bufSz);

    void finalizeThreads(const bool terminate = true);

    size_t getArgSz() const { return argSz; }
    size_t getRetSz() const { return retSz; }
    RemoteFuncProto getRealFunc() { return func; }
    task *getParentTask() const { return parentTask; }
    void setParentTask(task *oParentTask) { parentTask = oParentTask; }
    void lockParent() { rlock.acquire(mutex); }
    void unlockParent() { rlock.release(); }

  protected:
    typedef queuing_mutex ParentMutex;

    const size_t argSz;
    const size_t retSz;

    task *parentTask;
    RemoteFuncProto func;
    ParentMutex mutex;
    ParentMutex::scoped_lock rlock;
  };

  struct QueueData {
    message_queue *queue;
    std::string qname;
    named_mutex *mutex;
    std::string mname;

    QueueData() : queue(NULL), mutex(NULL) { };
    QueueData(const QueueData &othis) : queue(othis.queue), qname(othis.qname),
					mutex(othis.mutex), mname(othis.mname) { };

    ~QueueData() { };
  };

  typedef std::vector<RemoteFunction *> FuncArray;
  typedef queuing_mutex GAMutex;
  typedef queuing_mutex RetMutex;
  typedef concurrent_hash_map<unsigned int, bool> DirectMap;
  typedef std::map<int, int> ProcIdMap;
  typedef std::vector<QueueData> NodeQueues;

  struct GFutTask : public task {
  private:
    void *rarg;

    static void msgq_executor(void *arg);

  public:
    GFutTask() : rarg(NULL) { }
    GFutTask(void *orarg) : rarg(orarg) { }
    explicit GFutTask(const GFutTask &othis)
      : rarg(othis.rarg) { }
    ~GFutTask() { if (rarg) delete rarg; }
    task *execute() { msgq_executor(rarg); return NULL; };
  };
}

#endif // __GFUTEX_INTERNAL_H
