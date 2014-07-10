#ifndef _GA_CLASSIC_PORT_H
#define _GA_CLASSIC_PORT_H

class GAClassicPort: public virtual ::classic::gov::cca::Port {

 public:
  virtual GlobalArray* createGA(int type, int ndim, int dims[],
				char *arrayname,    int chunk[]) = 0;
  virtual GlobalArray* createGA(int type, int ndim, int dims[],
				char *arrayname, int maps[], int block[]) = 0;
  virtual GlobalArray* createGA(const GlobalArray *g_b, char *arrayname)  = 0;
  virtual GlobalArray* createGA(const GlobalArray &g_b) = 0;
  virtual GlobalArray* createGA() = 0;
  virtual GlobalArray* createGA_Ghosts(int type, int ndim, int dims[], 
				       int width[], char *array_name, 
				       int chunk[]) = 0;
  virtual GlobalArray* createGA_Ghosts(int type, int ndim, int dims[], 
				       int width[], char *array_name, 
				       int map[], int nblock[]) = 0;

  virtual void brdcst(void *buf, int lenbuf, int root) = 0;
  virtual int clusterNnodes() = 0;
  virtual int clusterNodeid() = 0;
  virtual int clusterNprocs(int inode) = 0;
  virtual int clusterProcid(int inode, int iproc) = 0;
  virtual int createMutexes(int number) = 0;
  virtual int destroyMutexes() = 0;
  virtual void dgop(double x[], int n, char *op) = 0;
  virtual int duplicate(int g_a, char* array_name) = 0;
  virtual void error(const char *message, int code) = 0;
  virtual void fence() = 0;
  virtual void igop(int x[], int n, char *op) = 0;
  virtual void initFence() = 0;
  virtual size_t inquireMemory() = 0;
  virtual void lgop(long x[], int n, char *op) = 0;
  virtual void lock(int mutex) = 0;
  virtual void maskSync(int first, int last) = 0;
  virtual int memoryAvailable() = 0;
  virtual int memoryLimited() = 0;
  virtual int nodeid() = 0;
  virtual int nodes() = 0;
  virtual void printStats() = 0;
  virtual void setMemoryLimit(size_t limit) = 0;
  virtual void summarize(int verbose) = 0;
  virtual void sync() = 0;
  virtual void unlock(int mutex) = 0;
  virtual int usesMA() = 0;
  virtual int usesFAPI() = 0;
  virtual void setServices(::classic::gov::cca::Services *cc) = 0;
};

#endif // _GA_CLASSIC_PORT_H
