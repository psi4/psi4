#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "mp3.h"

#include "ga++.h"

void
GA::Initialize(int argc, char *argv[], size_t limit) {
  MP_INIT(argc, argv);

  // GA Initialization
  if(limit == 0) 
    GA_Initialize();
  else 
    GA_Initialize_ltd(limit);
}

void 
GA::Initialize(int argc, char *argv[], unsigned long heapSize, 
	   unsigned long stackSize, int type, size_t limit) {
  MP_INIT(argc, argv);
  
  // GA Initialization
  if(limit == 0) 
    GA_Initialize();
  else 
    GA_Initialize_ltd(limit);

  
  //if(GA_Uses_ma()) {
  
  int nProcs = GA_Nnodes();
  
  // Initialize memory allocator
  heapSize /= ((unsigned long) nProcs);
  stackSize /= ((unsigned long) nProcs);
  
  if(!MA_init(type, stackSize, heapSize)) 
    GA_Error((char *)"MA_init failed",stackSize+heapSize);
  // }
}

void 
GA::Terminate()
{
  
  /* Terminate GA */
  GA_Terminate();    
  
  MP_FINALIZE();
}
