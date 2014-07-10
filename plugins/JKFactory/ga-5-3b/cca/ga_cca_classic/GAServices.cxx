/**              
 * module: GAServices.cc
 * Author: Manoj Kumar Krishnan, PNNL.
 */

#include "gacca.h"

#define GA_STACKSIZE 50000
#define GA_HEAPSIZE  50000


/**
 *  Constructor and Destructor of GAServices               
 */

GA::GAServices::GAServices() {

  svc = 0; /* services to NULL */
  
  // GA Initialization
  GA_Initialize();
  // later do it with parameter ports
  if(!MA_init(MT_F_DBL, GA_STACKSIZE, GA_HEAPSIZE))
    GA_Error((char *)"MA_init failed", GA_STACKSIZE+GA_HEAPSIZE);
}

GA::GAServices::~GAServices()  {
  svc = 0;
  
  // GA Termination
  GA_Terminate();
}

GA::GlobalArray *  
GA::GAServices::createGA(int type, int ndim, int dims[], 
			 char *arrayname, int chunk[]) {
  
  GA::GlobalArray * GA = new GA::GlobalArray(type, ndim, dims, arrayname, 
					     chunk);
  return GA;
}

GA::GlobalArray *  
GA::GAServices::createGA(int type, int ndim, int dims[], 
			 char *arrayname, int maps[], int block[]) {
  
  GA::GlobalArray * GA = new GA::GlobalArray(type, ndim, dims, arrayname, 
					     maps, block);
  return GA;
}

GA::GlobalArray *  
GA::GAServices::createGA(const GA::GlobalArray *g_b, char *arrayname) {
  GA::GlobalArray * GA = new GA::GlobalArray(*g_b, arrayname);
  return GA;
}

GA::GlobalArray *  
GA::GAServices::createGA(const GA::GlobalArray &g_b) {
  GA::GlobalArray * GA = new GA::GlobalArray(g_b);
  return GA;
}

GA::GlobalArray * 
GA::GAServices::createGA() {
  GA::GlobalArray * GA = new GA::GlobalArray();
  return GA;
}

GA::GlobalArray * 
GA::GAServices::createGA_Ghosts(int type, int ndim, int dims[], int width[],
				char *array_name, int chunk[]) {
  /* last argument is a dummy argument, just to increase the count of the 
     number of arguments, inorder to avoid conflict in # of args */
  GA::GlobalArray * GA = new GA::GlobalArray(type, ndim, dims, width, 
					     array_name, chunk, 'g');
  return GA;
}

GA::GlobalArray * 
GA::GAServices::createGA_Ghosts(int type, int ndim, int dims[], int width[], 
				char *array_name, int map[], int nblock[]) {
  GA::GlobalArray * GA = new GA::GlobalArray(type, ndim, dims, width,
					     array_name, map, nblock, 'g');
  return GA;
}

void 
GA::GAServices::brdcst(void *buf, int lenbuf, int root) {
  GA_Brdcst(buf, lenbuf, root);
}

int 
GA::GAServices::clusterNnodes() {
  return GA_Cluster_nnodes();
}

int 
GA::GAServices::clusterNodeid() {
  return GA_Cluster_nodeid();
}

int 
GA::GAServices::clusterNprocs(int inode) {
  return GA_Cluster_nprocs(inode) ;
}
  
int 
GA::GAServices::clusterProcid(int inode, int iproc) {
  return GA_Cluster_procid(inode, iproc);
}

int 
GA::GAServices::createMutexes(int number) {
  return GA_Create_mutexes(number);
}

int 
GA::GAServices::destroyMutexes() {
  return GA_Destroy_mutexes(); 
}

void 
GA::GAServices::dgop(double x[], int n, char *op) {
  GA_Dgop(x, n, op);
}

int 
GA::GAServices::duplicate(int g_a, char* array_name) {
  return GA_Duplicate(g_a, array_name);
}

void 
GA::GAServices::error(const char *message, int code) { 
  GA_Error((char *)message, code); 
}

void 
GA::GAServices::fence() {
  GA_Fence();
}

void 
GA::GAServices::igop(Integer x[], int n, char *op) {
  GA_Igop(x, n, op);
}

void 
GA::GAServices::initFence() {
  GA_Init_fence(); 
}

size_t 
GA::GAServices::inquireMemory() {
  return GA_Inquire_memory();
}

void 
GA::GAServices::lgop(long x[], int n, char *op) {
  GA_Lgop(x, n, op);
}

void 
GA::GAServices::lock(int mutex) {
  GA_Lock(mutex);
}

void 
GA::GAServices::maskSync(int first, int last) {
  GA_Mask_sync(first, last);
}

int 
GA::GAServices::memoryAvailable() {
   return GA_Memory_avail();
}

int 
GA::GAServices::memoryLimited() {
  return GA_Memory_limited();
}

int
GA::GAServices::nodeid() {
  return GA_Nodeid();
}

int
GA::GAServices::nodes() {
  return GA_Nnodes();
}

void 
GA::GAServices::printStats() {
  GA_Print_stats();
}

void 
GA::GAServices::setMemoryLimit(size_t limit) {
  GA_Set_memory_limit(limit);
}

void 
GA::GAServices::summarize(int verbose) {
  GA_Summarize(verbose);
}

void 
GA::GAServices::sync() { 
  GA_Sync(); 
}

void 
GA::GAServices::unlock(int mutex) {
  GA_Unlock(mutex);
}

int 
GA::GAServices::usesMA() {
  return GA_Uses_ma();
}

int 
GA::GAServices::usesFAPI() {
  return GA_Uses_fapi();
}

void
GA::GAServices::setServices(::classic::gov::cca::Services *cc){  

  int err = 0;
  classic::gov::cca::PortInfo* pInfo;
  
  IO_dn1("In DistArrayDescriptorFactory::setServices entry\n");
 
  // We're being shut down
  if (cc == 0) {
    // Are we shut down already?
    if (svc == 0) { return; }
 
    // Shutdown
    svc->removeProvidesPort("ga_classic_port");
    svc->removeProvidesPort("TemplateFactory");
    svc->removeProvidesPort("DescriptorFactory");
    svc = cc;
    return;
  }

  svc = cc;
 
  IO_dn1("In GA::GAServices::setServices entry\n");

  /****** Provide GA Classic Port ******/
  err = svc->addProvidesPort(dynamic_cast<GA::GAClassicPort *>(this),
			     svc->createPortInfo("ga_classic_port", "GAClassicPort", 0));  
  if ( err != 0) { 
    IO_dn1("In GA::GAServices::setServices addProvidesPort(ga_classic) failed\n");  
    ::abort(); 
  }
  pInfo = 0; // Current version of spec says we should destroy our copy
 
  
  /******  Provide DADF Template Factory Port *******/
  pInfo = svc->createPortInfo("TemplateFactory",
                              "DistArrayTemplFactoryPort", 0);
  err = svc->addProvidesPort(this, pInfo);
  if ( err != 0) { 
    IO_dn1("In GA::GAServices::setServices addProvidesPort(ga_dadf) failed\n");
    ::abort(); 
  }
  pInfo = 0; // Current version of spec says we should destroy our copy
  
  
  /******  Provide DADF Array Descriptor Factory Port *******/
  pInfo = svc->createPortInfo("DescriptorFactory",
                              "DistArrayDescrFactoryPort", 0);
  err = svc->addProvidesPort(this, pInfo);
  if ( err != 0) { 
    IO_dn1("In GA::GAServices::setServices addProvidesPort(ga_dadf) failed\n");
    ::abort(); 
  }
  pInfo = 0; // Current version of spec says we should destroy our copy
}
