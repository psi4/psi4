/* $Id:  */

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <string>
using std::string;
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
#define ARMCI_ENABLE_GPC_CALLS
#include "gpc.h"

/***************************** macros ************************/
extern "C" {
#  if defined(__ia64)
#    if defined(__GNUC__) && !defined (__INTEL_COMPILER)
#       define MEM_FENCE __asm__ __volatile__ ("mf" ::: "memory");
#    else /* Intel Compiler */
        extern void _armci_ia64_mb();
#       define MEM_FENCE _armci_ia64_mb();
#    endif
#  endif
}

#include "Hash_common.h"
#include "GPCHashmap.h"
#include "Util.h"

short int GPCHashmap::sm_initialized=0;

GPCHashmap::GPCHashmap()
  : mVocabMap(NULL)
{
}

GPCHashmap::~GPCHashmap()
{
    if(sm_initialized != 0) this->destroy();
}


void GPCHashmap::create() 
{
    if(sm_initialized != 0) {
      ARMCI_Error("GPCHashmap::create(): Hashmap already exists. At a given time, only one distributed hashmap should exist. Multiple distributed hashmaps not yet supported", 0);
    }
  
    mVocabMap = new VocabIntMap();
    sm_initialized=1;
}

void GPCHashmap::destroy() 
{
    if (mVocabMap != NULL) delete mVocabMap;
    sm_initialized=0;
}

/**
 * buf - character array
 * size - size of strlen array
 */
void GPCHashmap::insert(const char *buf, size_t bufsize) {

    armci_hashmap_insert(mVocabMap, buf, bufsize);
    
#ifdef MEM_FENCE
    MEM_FENCE;
#else
    ARMCI_Error("gpc_insert_handler: MEM_FENCE not defined", 0);
#endif
    
}

/**
 * get the global term IDs
 */
void GPCHashmap::getGlobalIds(const char *buf, size_t bufsize,
                              int *globalTermIds) {

  globalTermIds[0] = *((int*)buf);
  armci_hashmap_insert2(mVocabMap, buf, bufsize, &globalTermIds[1], HASHMAP_GET);
  
#ifdef MEM_FENCE
    MEM_FENCE;
#else
    ARMCI_Error("gpc_insert_handler: MEM_FENCE not defined", 0);
#endif
    
}

/**
 * prints the hashmap in this server process
 */ 
void GPCHashmap::print() 
{
    int me; MP_MYID(&me);
    printf("%d: Hashmap[%d:%ld]\n", me, 1, mVocabMap->size());

    VocabIntMap::const_iterator iter;
    for(iter=mVocabMap->begin();
        iter != mVocabMap->end();
        iter++) {
      cout << iter->second  << "\t: " << iter->first << endl;
    }
    cout << endl;
}

// to check if a hashmap exists
void GPCHashmap::rehash(int *size) 
{
  *size = (int) mVocabMap->size();

#ifdef MEM_FENCE
    MEM_FENCE;
#else
    ARMCI_Error("gpc_insert_handler: MEM_FENCE not defined", 0);
#endif
    
}

// to check if a hashmap exists
bool GPCHashmap::isCreated() 
{
  if (sm_initialized >0) return true;
  return false;
}
