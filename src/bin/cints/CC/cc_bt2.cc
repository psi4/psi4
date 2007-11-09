/*! \file 
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<memory.h>
#include<stdlib.h>
#include<pthread.h>
#include<cstdio>
#include<cstring>
#include<cmath>
#include<libciomr/libciomr.h>
#include<libint/libint.h>
#include<libqt/qt.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"Tools/int_fjt.h"

#define USE_SYMM_CODE 0

namespace psi { namespace CINTS {
  /*-------------------------------
    Explicit function declarations
    -------------------------------*/
  extern void *cc_bt2_thread(void *);
  extern void *cc_bt2_thread_symm(void *);
  
  /*!-------------------------------------------------------
    Algorithm
    
    ***Split into threads. Each thread does the following:
    
    Loop over all symmetry-unique shells UP, UQ<=UP
    ***if this (UP UQ| is to be handled by this thread - proceed, else skip to next one
    Find all symmetry-distinct shell doublets resulting from (UP UQ|UR US)
    Loop over all resulting shell quartets (PQ|RS)
    
        Loop over p in P, q in Q, r in R, s in S, ij
	bt2(pq,ij) +=  (pq|rs) * t2(rs,ij)
        End p, q, r, s, ij loop
	
	End P,Q,R,S loop
	End UP, UQ loop
	
	***Merge all threads
	
  -------------------------------------------------------*/
  void cc_bt2()
  {
    pthread_attr_t thread_attr;
    pthread_t *thread_id;
    
    long int libint_memory;
    int max_bf_per_shell;
    int max_num_prim_comb;
    int i;
    
    /*---------------
      Initialization
      ---------------*/
    init_fjt(BasisSet.max_am*4+1);
    init_libint_base();
    timer_init();
    
    /*-------------------------
      Allocate data structures
      -------------------------*/
    max_bf_per_shell = ioff[BasisSet.max_am];
    /*--- Use this dirty trick to get how much memory integrals library needs ---*/
    max_num_prim_comb = (BasisSet.max_num_prims*
			 BasisSet.max_num_prims)*
      (BasisSet.max_num_prims*
       BasisSet.max_num_prims);
    libint_memory = libint_storage_required(BasisSet.max_am-1,max_num_prim_comb);
    UserOptions.memory -= libint_memory*UserOptions.num_threads;
    
    /*---
      Assuming 1 pass for now
      ---*/
    
    /*--------------------------
      Start compute threads now
      --------------------------*/
    thread_id = (pthread_t *) malloc(UserOptions.num_threads*sizeof(pthread_t));
    pthread_attr_init(&thread_attr);
    pthread_attr_setscope(&thread_attr,
			  PTHREAD_SCOPE_SYSTEM);
#if USE_SYMM_CODE
    for(long int i=0;i<UserOptions.num_threads-1;i++)
      pthread_create(&(thread_id[i]),&thread_attr,
		     cc_bt2_thread_symm,(void *)i);
    cc_bt2_thread_symm( (void *) (UserOptions.num_threads - 1) );
#else
    for(long int i=0;i<UserOptions.num_threads-1;i++)
      pthread_create(&(thread_id[i]),&thread_attr,
		     cc_bt2_thread,(void *)i);
    cc_bt2_thread( (void *) (UserOptions.num_threads - 1) );
#endif
    for(i=0;i<UserOptions.num_threads-1;i++)
      pthread_join(thread_id[i], NULL);
    free(thread_id);
    
    /*---------
      Clean-up
      ---------*/
    free_fjt();
    timer_done();
    
    return;
  }
};
};
