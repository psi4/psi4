/* $Id: DistHashmap.cc,v 1.1.2.1 2007-06-20 17:41:53 vinod Exp $  */
/*
 * AUTHOR: Manojkumar Krishnan, PNNL
 */
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;

#define DEBUG 0

#ifdef WIN32
#  include <windows.h>
#  define sleep(x) Sleep(1000*(x))
#else
#  include <unistd.h>
#endif

#include "armci.h"
#include "message.h"
#define ARMCI_ENABLE_GPC_CALLS
#include "gpc.h"

#include "Hash_common.h"
#include "DistHashmap.h"
#include "Util.h"

#define BUF_LIMIT 8192 // in bytes. Make sure this value can fit in an
                       // integer variable

extern int gpc_disthashmap_handler(int to, int from, void *hdr,   int hlen,
                                   void *data,  int dlen,
                                   void *rhdr,  int rhlen, int *rhsize,
                                   void *rdata, int rdlen, int *rdsize,
                                   int rtype);
extern void gpc_disthashmap_exec(int hash_op, char *buf, size_t bufsize,
                                 int proc, int gpc_handle);
extern void gpc_disthashmap_exec_nb(int hash_op, char *buf, size_t bufsize,
                                    int proc, int gpc_handle, gpc_hdl_t *nbh);
extern void gpc_disthashmap_exec_wait(gpc_hdl_t *nbh);

extern int armci_master;

short int DistHashmap::sm_initialized=0;

DistHashmap::DistHashmap()
  : mGPCHandle(0), mGlobalIdMap(NULL), m_globalHashMapSize(0)
{
    MP_MYID(&m_me);

    int nproc;
    MP_PROCS(&nproc);
    
    mIndex = new int[nproc];
    for(int i=0; i<nproc; i++) mIndex[i] = 0;

    mBuf = new char*[nproc];    
    for(int i=0; i<nproc; i++) {
      mBuf[i] = new char[BUF_LIMIT];
    }
    mTmpBuf = new char[BUF_LIMIT];
}

DistHashmap::~DistHashmap()
{
    if(sm_initialized != 0) this->destroy();

    delete[] mIndex;

    int nproc;
    MP_PROCS(&nproc);
    for(int i=0; i<nproc; i++) delete[] mBuf[i];
    delete[] mBuf;
    delete[] mTmpBuf;
}

/**
 * This is a collective call
 */
void DistHashmap::create() 
{
    if(sm_initialized != 0) {
      ARMCI_Error("DistHashmap::create(): A Distributed Hashmap already exists. At a given time, only one distributed hashmap should exist. Multiple distributed hashmaps not yet supported", 0);
    }

    // to cache the globalIds locally
    mGlobalIdMap = new VocabIntMap();
    
    MP_BARRIER();
    mGPCHandle = ARMCI_Gpc_register((int (*)())gpc_disthashmap_handler);

    /* SMP master process creates the hashmap */
    if(m_me == armci_master) {
      gpc_disthashmap_exec(HASHMAP_CREATE, NULL, 0, m_me, mGPCHandle);
    }
    MP_BARRIER();

    sm_initialized=1;
}

/**
 * This is a collective call
 */
void DistHashmap::destroy() 
{

    if (mGlobalIdMap != NULL) delete mGlobalIdMap;
  
    MP_BARRIER();  
    /* SMP master process destroys the hashmap */
    if(m_me == armci_master) {
      gpc_disthashmap_exec(HASHMAP_DESTROY, NULL, 0, m_me, mGPCHandle);
    }
    
    MP_BARRIER();
    
    sm_initialized=0;
    if(mGPCHandle != 0) ARMCI_Gpc_release(mGPCHandle);
}

/**
 * str - string to be inserted into the distributed hashmap
 * size - size of strlen array
 * Insert() is just to ensure that there is a new entry to be inserted into
 * the distributed hash map, and this entry is marked in the local buffer as
 * "to be inserted". The reason is, it will be expensive to do the insert for
 * each and every new element. In order to avoid those latency costs, insert
 * aggregates the entries into a single buffer of size BUF_LIMIT and only
 * completes the insertion into the distributed hashh map, if and only if
 * this aggregate local buffer is full. In order to complete this aggregation
 * call (even though the buffer is not full), commit() can be used to flush the
 * local buffer into the distributed hash map. So it is a good practice to
 * call commit after all insert or whenever consistency is required. 
 */
void DistHashmap::insert(string str) {

  // identify the destination process to insert (using hash function)
  int nproc;
  MP_PROCS(&nproc);
  int dst_proc = armci_djb2_hash((unsigned char*)str.c_str(), nproc);
  int padding = 2*sizeof(int); // string length + int-byte alignment
  
  /*
   * The structure of this internal buffer is:
   *  ----------------------------------------------------------------------
   * | # of Strings | 1st string length | 1st string | 2nd string length | ...
   *  ----------------------------------------------------------------------
   * For example:
   *       ---------------------------------------
   *      | 3 | 2 | Hello | 3 | Sir | 7 | Welcome |
   *       ---------------------------------------
   */

  if(str.length() > BUF_LIMIT)
    ARMCI_Error("String length is greater than internal buffer. Increase BUF_LIMIT", 0);
  
  // if the buffer is full, trigger an insert into the GPC hashmap
  if(mIndex[dst_proc] + str.length() + padding  > BUF_LIMIT) this->commit(dst_proc);

  // initalize the "number of strings" field in the buffer
  int *numstrings = (int*)mBuf[dst_proc];
  if(mIndex[dst_proc] == 0) {
    *numstrings = 0;
    mIndex[dst_proc] += sizeof(int);
  }
  (*numstrings)++;

  // insert into the internal buffer
  int index = mIndex[dst_proc];
  mIndex[dst_proc] += armci_hashmap_pack(&mBuf[dst_proc][index], str);
}

/**
 * Call commit() to complete the data insertion into the distributed hashmap
 * for consistency purposes. This is a collective call.
 */
void DistHashmap::commit() {

  int nproc; MP_PROCS(&nproc);

  for(int i=0; i<nproc; i++) this->commit(i);
  this->rehash();
}

/**
 * Call commit() to complete the data insertion into the distributed hashmap
 * for consistency purposes (of a specified process).
 */
void DistHashmap::commit(int proc) {

  // complete the hashmap insertion if the buffer still has data
  if(mIndex[proc] > 0) {
    char *buf = &mBuf[proc][0];
    int bufsize  = mIndex[proc];
#   if DEBUG
    printf("%d: Commit() to %d: numstrings=%d bufsize=%d\n",
           m_me, proc, *((int*)buf), bufsize);fflush(stdout);
#   endif

    memcpy(mTmpBuf, (const char*)buf, (size_t)bufsize);
    int tmpBufsize = bufsize;

    // insert into the distributed hashmap
    gpc_disthashmap_exec(HASHMAP_INSERT, buf, bufsize, proc, mGPCHandle);
    mIndex[proc] = 0;

    // verify the returned data
    if( *((int*)buf) != *((int*)mTmpBuf) )
      ARMCI_Error("DistHashmap::commit() failed", 0);

    // i.e. save the global ids in the locally cached hashmap (mGlobalIdMap)
    int *globalTermIds = (int*) (buf + sizeof(int));
    
    armci_hashmap_insert2(mGlobalIdMap, (const char*)mTmpBuf, tmpBufsize,
                          globalTermIds, HASHMAP_INSERT);
  }
}

/**
 * This is a collective call to rehash the distributed hashmap
 */
void DistHashmap::rehash() 
{

    int nclus = armci_domain_count(ARMCI_DOMAIN_SMP);
    int clus_me = armci_domain_my_id(ARMCI_DOMAIN_SMP);

    // no need to rehash if #of clusters is <= 1
    if(nclus <=1) {
      m_globalHashMapSize = mGlobalIdMap->size();
      return;
    }

    // no GPC support for single node multiple processes
    if(nclus == 1) {
      int nproc; MP_PROCS(&nproc);
      if(nproc > 1) ARMCI_Error("DistHashmap::rehash(): no GPC support for single node multiple processes", 0);
    }
    
    int size=0;
    int *globalHashOffset = new int[nclus];
    int *offsetMap = new int[nclus];

    if(globalHashOffset == NULL || offsetMap == NULL)
      ARMCI_Error("DistHashmap::rehash() new alloc failed", 0);
    for(int i=0; i<nclus; i++) globalHashOffset[i]=0;
    
    MP_BARRIER();
    ARMCI_AllFence();
    
    /* SMP master process get the size of its GPC Hashmap */
    if(m_me == armci_master) {
      gpc_disthashmap_exec(HASHMAP_REHASH, (char*)&size, sizeof(int), m_me,
                           mGPCHandle);

      if(size <0) ARMCI_Error("DistHashmap::rehash(): invalid size", 0);

      // get the global hasp map size
      m_globalHashMapSize = size;
      armci_msg_gop_scope(SCOPE_MASTERS, &m_globalHashMapSize, 1, "+",
                          ARMCI_INT);
      
      
      // masters only gop
      globalHashOffset[clus_me] = size;
      armci_msg_gop_scope(SCOPE_MASTERS, globalHashOffset, nclus, "+",
                          ARMCI_INT);
    }
    
    MP_BARRIER();  

    // cluster bcast to processes in the same node 
    armci_msg_clus_brdcst(&m_globalHashMapSize, 1*sizeof(int));
    armci_msg_clus_brdcst(globalHashOffset, nclus*sizeof(int));
      
    if(m_me==0) {
      for(int i=0; i<nclus; i++)
        printf("globalHashOffset[%d] = %d\n", i, globalHashOffset[i]);
      fflush(stdout);
    }
    
    // global offset in a "distribution map" fashion
    offsetMap[0] = 0;
    for(int i=1; i<nclus; i++)
      offsetMap[i] = offsetMap[i-1] + globalHashOffset[i-1];

    if(m_me==0) {
      for(int i=0; i<nclus; i++)
        printf("offsetMap[%d] = %d\n", i, offsetMap[i]); fflush(stdout);
    }

    // Now, every processes have the hashmap size of master processes who
    // maintain the dist hashmap. So offset is done accordingly, and the
    // term ids in the local hashmap is updated with the global ids.
    VocabIntMap::iterator iter;
    string termStr;
    for(iter=mGlobalIdMap->begin(); iter!=mGlobalIdMap->end(); iter++) {
      termStr = (*iter).first;
      int nproc; MP_PROCS(&nproc);
      int dst_proc = armci_djb2_hash((unsigned char*)termStr.c_str(), nproc);
      int dst_clus = armci_domain_id(ARMCI_DOMAIN_SMP, dst_proc);
      (*iter).second += offsetMap[dst_clus];
    }

    delete[] offsetMap;
    delete[] globalHashOffset;
}

/**
 * to check if a distributed hashmap exists
 * This is a local call
 */
bool DistHashmap::isCreated() 
{
  if (sm_initialized >0) return true;
  return false;
}

/**
 * This is a collective call
 */
void DistHashmap::print() 
{
    MP_BARRIER();
    
    /* SMP master process destroys the hashmap */
    if(m_me == armci_master) {
      gpc_disthashmap_exec(HASHMAP_PRINT, NULL, 0, m_me, mGPCHandle);
    }
    
    MP_BARRIER();
}

/**
 * This is a local call
 */
void DistHashmap::print2() 
{
  printf("%d: Locally cached Hashmap[%d:%ld]\n", m_me, 1, mGlobalIdMap->size());

  VocabIntMap::const_iterator iter;
  for(iter=mGlobalIdMap->begin();
      iter != mGlobalIdMap->end();
      iter++) {
    cout << iter->second  << "\t: " << iter->first << endl;
  }
  cout << endl;
  
}

const VocabIntMap *
DistHashmap:: getLocalMapPtr() 
{
  return mGlobalIdMap;
}

int
DistHashmap:: getGlobalHashMapSize() 
{
  return m_globalHashMapSize;
}

/*
  TODO:
  1. put the gpc exec functions in GPCHashmapHandler.cc as private
  methods in this class
  2. for better performance, use ARMCI_Malloc instead of new for data transfer GPC buffers
*/

