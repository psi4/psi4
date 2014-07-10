#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "ga++.h"

GA::GlobalArray *  
GA::createGA(int type, int ndim, int dims[], char *arrayname, int chunk[]) {
  return new GA::GlobalArray(type, ndim, dims, arrayname, chunk);
}

GA::GlobalArray *  
GA::createGA(int type, int ndim, int dims[], char *arrayname, int block[],
        int maps[]) {
  return new GA::GlobalArray(type, ndim, dims, arrayname, block, maps);
}

GA::GlobalArray *  
GA::createGA(const GA::GlobalArray *g_b, char *arrayname) {
  return new GA::GlobalArray(*g_b, arrayname);
}

GA::GlobalArray *  
GA::createGA(const GA::GlobalArray &g_b) {
  return new GA::GlobalArray(g_b);
}

GA::GlobalArray * 
GA::createGA() {
  return new GA::GlobalArray();
}

GA::GlobalArray * 
GA::createGA_Ghosts(int type, int ndim, int dims[], int width[],
				char *array_name, int chunk[]) {
  /* last argument is a dummy argument, just to increase the count of the 
     number of arguments, inorder to avoid conflict in # of args */
  return new GA::GlobalArray(type, ndim, dims, width, array_name, chunk, 'g');
}

GA::GlobalArray * 
GA::createGA_Ghosts(int type, int ndim, int dims[], int width[], 
				char *array_name, int map[], int nblock[]) {
  return new GA::GlobalArray(type, ndim, dims, width, array_name, map, nblock, 'g');
}

int
GA::getDebug() {
    return GA_Get_debug();
}

void 
GA::brdcst(void *buf, int lenbuf, int root) {
  GA_Brdcst(buf, lenbuf, root);
}

int 
GA::clusterNnodes() {
  return GA_Cluster_nnodes();
}

int 
GA::clusterNodeid() {
  return GA_Cluster_nodeid();
}

int 
GA::clusterProcNodeid(int iproc) {
  return GA_Cluster_proc_nodeid(iproc);
}

int 
GA::clusterNprocs(int inode) {
  return GA_Cluster_nprocs(inode);
}
  
int 
GA::clusterProcid(int inode, int iproc) {
  return GA_Cluster_procid(inode, iproc);
}

int 
GA::createMutexes(int number) {
  return GA_Create_mutexes(number);
}

int
GA::deregisterType(int size) {
  return NGA_Deregister_type(size);
}

int 
GA::destroyMutexes() {
  return GA_Destroy_mutexes(); 
}

void 
GA::dgop(double x[], int n, char *op) {
  GA_Dgop(x, n, op);
}

int 
GA::duplicate(int g_a, char* array_name) {
  return GA_Duplicate(g_a, array_name);
}

void 
GA::error(const char *message, int code) { 
  GA_Error((char *)message, code); 
}

void 
GA::fence() {
  GA_Fence();
}

void
GA::gop(int x[], int n, char *op) {
    GA_Igop(x, n, op);
}

void
GA::gop(long x[], int n, char *op) {
    GA_Lgop(x, n, op);
}

void
GA::gop(float x[], int n, char *op) {
    GA_Fgop(x, n, op);
}

void
GA::gop(double x[], int n, char *op) {
    GA_Dgop(x, n, op);
}

void 
GA::igop(int x[], int n, char *op) {
  GA_Igop(x, n, op);
}

void 
GA::initFence() {
  GA_Init_fence(); 
}

size_t 
GA::inquireMemory() {
  return GA_Inquire_memory();
}

void 
GA::lgop(long x[], int n, char *op) {
  GA_Lgop(x, n, op);
}

void 
GA::lock(int mutex) {
  GA_Lock(mutex);
}

void 
GA::maskSync(int first, int last) {
  GA_Mask_sync(first, last);
}

int 
GA::memoryAvailable() {
   return GA_Memory_avail();
}

int 
GA::memoryLimited() {
  return GA_Memory_limited();
}

void GA::nbWait(GANbhdl *nbhandle) {
  NGA_NbWait(nbhandle);
}

int
GA::nodeid() {
  return GA_Nodeid();
}

int
GA::nodes() {
  return GA_Nnodes();
}

void 
GA::printStats() {
  GA_Print_stats();
}

int
GA::registerType(size_t size) {
  return NGA_Register_type(size);
}

void
GA::setDebug(int dbg) {
    return GA_Set_debug(dbg);
}

void 
GA::setMemoryLimit(size_t limit) {
  GA_Set_memory_limit(limit);
}

void 
GA::summarize(int verbose) {
  GA_Summarize(verbose);
}

void 
GA::sync() { 
  GA_Sync(); 
}

void 
GA::unlock(int mutex) {
  GA_Unlock(mutex);
}

int 
GA::usesMA() {
  return GA_Uses_ma();
}

int 
GA::usesFAPI() {
  return GA_Uses_fapi();
}

double GA::wtime() {
    return GA_Wtime();
}
