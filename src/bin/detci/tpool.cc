/*! \file
 *  \ingroup DETCI
 *  \brief Thread pools
 *
 * An example source module to accompany...
 *
 * "Using POSIX Threads: Programming with Pthreads"
 *     by Brad nichols, Dick Buttlar, Jackie Farrell
 *     O'Reilly & Associates, Inc.
 *
 ********************************************************
 * tpool.c -- 
 * 
 * Example thread pooling library
 */

#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <sys/types.h>
#include <cstring>

#include <pthread.h>
#include "tpool.h"
#include "structs.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

void *tpool_thread(void *);

void tpool_init(tpool_t   *tpoolp,
		int       num_worker_threads, 
		int       max_queue_size,
		int       do_not_block_when_full)
{
  int i, rtn;
  tpool_t tpool;
  std::string str;
   
  /* allocate a pool data structure */ 
  if ((tpool = (tpool_t )malloc(sizeof(struct tpool))) == NULL)
    throw PsiException("malloc error",__FILE__,__LINE__);

  /* initialize th fields */
  tpool->num_threads = num_worker_threads;
  tpool->max_queue_size = max_queue_size;
  tpool->do_not_block_when_full = do_not_block_when_full;
  if ((tpool->threads = 
       (pthread_t *)malloc(sizeof(pthread_t)*num_worker_threads)) 
      == NULL)
    throw PsiException("malloc error",__FILE__,__LINE__);
  tpool->cur_queue_size = 0;
  tpool->queue_head = NULL; 
  tpool->queue_tail = NULL;
  tpool->queue_closed = 0;  
  tpool->shutdown = 0; 
  tpool->threads_awake = 0;
  if ((rtn = pthread_mutex_init(&(tpool->queue_lock), NULL)) != 0){
    str = "pthread_mutex_init ";
    str += strerror(rtn);
    throw PsiException(str,__FILE__,__LINE__);
  }
  if ((rtn = pthread_cond_init(&(tpool->queue_not_empty), NULL)) != 0){
    str = "pthread_cond_init ";
    str += strerror(rtn);
    throw PsiException(str,__FILE__,__LINE__);
  }
  if ((rtn = pthread_cond_init(&(tpool->queue_not_full), NULL)) != 0){
    str = "pthread_cond_init ";
    str += strerror(rtn);
    throw PsiException(str,__FILE__,__LINE__);
  }
  if ((rtn = pthread_cond_init(&(tpool->queue_empty), NULL)) != 0){
    str = "pthread_cond_init ";
    str += strerror(rtn);
    throw PsiException(str,__FILE__,__LINE__);
  }
  if ((rtn = pthread_cond_init(&(tpool->all_work_done), NULL)) != 0){
    str = "pthread_cond_init ";
    str += strerror(rtn);
    throw PsiException(str,__FILE__,__LINE__);
  }

  /* create threads */
  for (i=0; i<num_worker_threads; i++) {
    if ((rtn = pthread_create( &(tpool->threads[i]),
			      NULL,
			      tpool_thread,
			      (void *)tpool)) != 0){
      str = "pthread_create ";
      str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
      throw PsiException(str,__FILE__,__LINE__);
  }
    tpool->threads_awake++;
  }

  *tpoolp = tpool;
}

int tpool_add_work(
		   tpool_t          tpool,
		   void             (*routine)(void *),
		   void             *arg)
{
  std::string str;
  int rtn;
  tpool_work_t *workp;

  if ((rtn = pthread_mutex_lock(&(tpool->queue_lock))) != 0){
    str = "pthread_mutex_lock ";
    str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
    throw PsiException(str,__FILE__,__LINE__);
  }

  /* no space and this caller doesn't want to wait */
  if ((tpool->cur_queue_size == tpool->max_queue_size) &&
      tpool->do_not_block_when_full) {
    if ((rtn = pthread_mutex_unlock(&(tpool->queue_lock))) != 0){
      str = "pthread_mutex_unlock ";
      str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
      throw PsiException(str,__FILE__,__LINE__);
    }

    return -1;
  }

  while( (tpool->cur_queue_size == tpool->max_queue_size) &&
	(!(tpool->shutdown || tpool->queue_closed))  ) {

    if ((rtn = pthread_cond_wait(&(tpool->queue_not_full),
				 &(tpool->queue_lock))) != 0){
      str = "pthread_cond_wait ";
      str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
      throw PsiException(str,__FILE__,__LINE__);
    }

  }

  /* the pool is in the process of being destroyed */
  if (tpool->shutdown || tpool->queue_closed) {
    if ((rtn = pthread_mutex_unlock(&(tpool->queue_lock))) != 0){
      str = "pthread_mutex_unlock ";
      str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
      throw PsiException(str,__FILE__,__LINE__);
    }
 
    return -1;
  }


  /* allocate work structure */
  if ((workp = (tpool_work_t *)malloc(sizeof(tpool_work_t))) == NULL){
    throw PsiException("pthread malloc error",__FILE__,__LINE__);
  }
  workp->routine = routine;
  workp->arg = arg;
  workp->next = NULL;

   if (tpool->cur_queue_size == 0) {
    tpool->queue_tail = tpool->queue_head = workp;


    if ((rtn = pthread_cond_broadcast(&(tpool->queue_not_empty))) != 0){
      str = "pthread_cond_signal ";
      str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
      throw PsiException(str,__FILE__,__LINE__);
    }
  } else {
    tpool->queue_tail->next = workp;
    tpool->queue_tail = workp;
  }

  tpool->cur_queue_size++; 
  if ((rtn = pthread_mutex_unlock(&(tpool->queue_lock))) != 0){
    str = "pthread_mutex_unlock ";
    str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
    throw PsiException(str,__FILE__,__LINE__);
  }
  return 1;
}

int tpool_destroy(tpool_t          tpool,
		  int              finish)
{
  std::string str;
  int          i,rtn;
  tpool_work_t *cur_nodep;
  

  if ((rtn = pthread_mutex_lock(&(tpool->queue_lock))) != 0){
    str = "pthread_mutex_lock ";
    str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
    throw PsiException(str,__FILE__,__LINE__);
  }

  /* Is a shutdown already in progress? */
  if (tpool->queue_closed && tpool->shutdown) {
    if ((rtn = pthread_mutex_unlock(&(tpool->queue_lock))) != 0){
      str = "pthread_mutex_unlock ";
      str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
      throw PsiException(str,__FILE__,__LINE__);
    }
    return 0;
  }

  tpool->queue_closed = 1;

  /* If the finish flag is set, wait for workers to 
     drain queue */ 
  if (finish == 1) {
    while (tpool->cur_queue_size != 0) {
      if ((rtn = pthread_cond_wait(&(tpool->queue_empty),
				   &(tpool->queue_lock))) != 0){
        str = "pthread_cond_wait ";
        str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
        throw PsiException(str,__FILE__,__LINE__);
      }
    }
  }

  tpool->shutdown = 1;

  if ((rtn = pthread_mutex_unlock(&(tpool->queue_lock))) != 0){
    str = "pthread_mutex_unlock ";
    str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
    throw PsiException(str,__FILE__,__LINE__);
  }


  /* Wake up any workers so they recheck shutdown flag */
  if ((rtn = pthread_cond_broadcast(&(tpool->queue_not_empty))) != 0){
    str = "pthread_cond_broadcast ";
    str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
    throw PsiException(str,__FILE__,__LINE__);
  }
  if ((rtn = pthread_cond_broadcast(&(tpool->queue_not_full))) != 0){
    str = "pthread_cond_broadcast ";
    str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
    throw PsiException(str,__FILE__,__LINE__);
  }


  /* Wait for workers to exit */
  for(i=0; i < tpool->num_threads; i++) {
      if ((rtn = pthread_join(tpool->threads[i],NULL)) != 0){
          str = "pthread_join ";
          str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
          throw PsiException(str,__FILE__,__LINE__);
      }
    }

  /* Now free pool structures */
  free(tpool->threads);
  while(tpool->queue_head != NULL) {
    cur_nodep = tpool->queue_head->next; 
    tpool->queue_head = tpool->queue_head->next;
    free(cur_nodep);
  }
  free(tpool); 
  return(1);
}
void tpool_queue_open(tpool_t tpool)
{
  pthread_mutex_lock(&tpool->queue_lock);
  tpool->queue_closed = 0;
  tpool->threads_awake = 0;
  pthread_mutex_unlock(&tpool->queue_lock);
}
  
void tpool_queue_close(tpool_t tpool, int finish)
{
  std::string str;
  int rtn;
  
  pthread_mutex_lock(&tpool->queue_lock);
  tpool->queue_closed = 1;
  
  if (finish) {
      if (tpool->cur_queue_size !=0 || tpool->threads_awake != 0) {
          if ((rtn = pthread_cond_wait(&(tpool->all_work_done), &(tpool->queue_lock))) !=0){
              str = "pthread_cond_wait ";
              str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
              throw PsiException(str,__FILE__,__LINE__);
          }
        }
    
    }
  
  pthread_mutex_unlock(&tpool->queue_lock);
  
}
  
void *tpool_thread(void *arg)
{
  std::string str;
  tpool_t tpool = (tpool_t)arg; 
  int rtn;
  tpool_work_t	*my_workp;
	
  for(;;) {

    /* Check queue for work */ 
    if ((rtn = pthread_mutex_lock(&(tpool->queue_lock))) != 0){
      str = "pthread_mutex_lock ";
      str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
      throw PsiException(str,__FILE__,__LINE__);
    }

    while ((tpool->cur_queue_size == 0) && (!tpool->shutdown)) {

        tpool->threads_awake--;
                 
        if (tpool->threads_awake == 0 && tpool->queue_closed) {
            if ((rtn = pthread_cond_signal(&(tpool->all_work_done))) != 0){
                str = "pthread_cond_signal ";
                str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
                throw PsiException(str,__FILE__,__LINE__);
            }
          }
      
        if ((rtn = pthread_cond_wait(&(tpool->queue_not_empty),
                                     &(tpool->queue_lock))) != 0){
            str = "pthread_cond_wait ";
            str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
            throw PsiException(str,__FILE__,__LINE__);
        }
        
        tpool->threads_awake++;
      }
 
 
    /* Has a shutdown started while i was sleeping? */
    if (tpool->shutdown == 1) {
        tpool->threads_awake--;
        if ((rtn = pthread_mutex_unlock(&(tpool->queue_lock))) != 0){
            str = "pthread_mutex_unlock ";
            str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
            throw PsiException(str,__FILE__,__LINE__);
        }
        pthread_exit(NULL);
      }


    /* Get to work, dequeue the next item */ 
    my_workp = tpool->queue_head;
    tpool->cur_queue_size--;
    if (tpool->cur_queue_size == 0)
      tpool->queue_head = tpool->queue_tail = NULL;
    else
      tpool->queue_head = my_workp->next;
 
    /* Handle waiting add_work threads */
    if ((!tpool->do_not_block_when_full) &&
	(tpool->cur_queue_size ==  (tpool->max_queue_size - 1))) 

      if ((rtn = pthread_cond_broadcast(&(tpool->queue_not_full))) != 0){
        str = "pthread_cond_broadcast ";
        str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
        throw PsiException(str,__FILE__,__LINE__);
      }

    /* Handle waiting destroyer threads */
    if (tpool->cur_queue_size == 0)
      if ((rtn = pthread_cond_signal(&(tpool->queue_empty))) != 0){
        str = "pthread_cond_signal ";
        str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
        throw PsiException(str,__FILE__,__LINE__);
      }

    if ((rtn = pthread_mutex_unlock(&(tpool->queue_lock))) != 0){
      str = "pthread_mutex_unlock ";
      str += static_cast<std::ostringstream*>( &(std::ostringstream() << rtn) )->str();
      throw PsiException(str,__FILE__,__LINE__);
    }
      
    /* Do this work item */
    (*(my_workp->routine))(my_workp->arg);
    free(my_workp);

  } 
  return(NULL);            
}

}} // namespace psi::detci

