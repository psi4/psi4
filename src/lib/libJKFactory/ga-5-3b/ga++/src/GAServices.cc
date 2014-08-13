/**              
 * @file GAServices.cc
 * @author Manoj Kumar Krishnan, PNNL.
 */
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "ga++.h"

GA::GAServices::GAServices()
{
}

GA::GAServices::~GAServices() 
{
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
			 char *arrayname, int block[], int maps[]) {
  
  GA::GlobalArray * GA = new GA::GlobalArray(type, ndim, dims, arrayname, 
					     block, maps);
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

int
GA::GAServices::getDebug() {
    return GA_Get_debug();
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
GA::GAServices::clusterProcNodeid(int iproc) {
  return GA_Cluster_proc_nodeid(iproc);
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
GA::GAServices::deregisterType(int type) {
  return NGA_Deregister_type(type);
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
GA::GAServices::gop(int x[], int n, char *op) {
    GA_Igop(x, n, op);
}

void
GA::GAServices::gop(long x[], int n, char *op) {
    GA_Lgop(x, n, op);
}

void
GA::GAServices::gop(float x[], int n, char *op) {
    GA_Fgop(x, n, op);
}

void
GA::GAServices::gop(double x[], int n, char *op) {
    GA_Dgop(x, n, op);
}

void 
GA::GAServices::igop(int x[], int n, char *op) {
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

void GA::GAServices::nbWait(GANbhdl *nbhandle) {
  NGA_NbWait(nbhandle);
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

int
GA::GAServices::registerType(size_t size) {
  return NGA_Register_type(size);
}

void
GA::GAServices::setDebug(int dbg) {
    return GA_Set_debug(dbg);
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

double GA::GAServices::wtime() {
    return GA_Wtime();
}

GA::GAServices GA::SERVICES;
