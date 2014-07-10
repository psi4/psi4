#ifndef __tascel_SplitQueue_h__
#define __tascel_SplitQueue_h__

#include "TerminationDetector.h"

namespace tascel {

  /**
   * A set of task queues shared among processes.
   *
   * Each process has a local queue with a maximum number of possible tasks.
   * A process can steal tasks from any other process's queue and have
   * tasks stolen from its own queue.
   *
   * This implementation splits the queue into a local and shared portion.
   * This allows the local process to modify its queue without locking.
   */
  class SplitQueue {
    private:
      /**
       * The local queue state which other procs can modify.
       */
      struct sq_state_t {
        int dirty;      /**< whether a steal has occurred */
        int split;      /**< the split between head and tail */
        int tail;       /**< the tail of the queue */
        int size_shared;/**< the number of shared tasks in the queue */
      };

      const int max_ntsks;    /**< max number of tasks allowed in queue */
      const int tsk_size;     /**< size of a single task */
      char **q;               /**< addresses of queue data on all procs */
      sq_state_t **sq_state;  /**< addresses of queue state on all procs */
      int *head;              /**< pointer to this procs head_val */
      int *tail;              /**< pointer to this procs tail state */
      int *size_local;        /**< pointer to this procs size_local_val */
      int *size_shared;       /**< pointer to this procs size_shared state */
      int *split;             /**< pointer to this procs split state */
      int *dirty;             /**< pointer to this procs dirt state */
      int size_local_val;     /**< size of local portion of local queue */
      int head_val;           /**< index of the local queue's head */
      TerminationDetector td; /**< the termination detector */

      /**
       * Moves tasks from the shared portion to the local portion.
       *
       * @param[in] _nacquire number of tasks to move
       */
      int acquireFromShared(int _nacquire);

      /**
       * Moves tasks from the local portion to the shared portion.
       *
       * @param[in] ndonate number of tasks to move
       */
      void releaseToShared(int ndonate);

    public:
      /**
       * Constructs the SplitQueue instance.
       *
       * @param[in] tsk_size size of a single task
       * @param[in] max_ntsks max number of tasks allowed in a queue
       *
       * @pre tsk_size must be the same value on all procs
       * @pre max_ntsks must be the same value on all procs
       */
      SplitQueue(int tsk_size, int max_ntsks);

      /**
       * Destroys the SplitQueue instance.
       */
      ~SplitQueue();

      /**
       * Returns true if the local queue is empty.
       */
      bool empty() const;

      /**
       * Returns true if the local queue is full.
       */
      bool full() const;

      /**
       * Retrieves a task from the private portion of the local queue.
       *
       * @param[out] dscr the retrieved task
       * @param[in] dlen size of the retrieved task
       *
       * @pre dscr is not NULL
       * @pre dlen == tsk_size
       */
      bool getTask(void *dscr, int dlen);

      /**
       * Adds a task to the private portion of the local queue.
       *
       * @param[in] dscr the task to add
       * @param[in] dlen size of the added task
       *
       * @pre dscr is not NULL
       * @pre dlen == tsk_size
       */
      void addTask(void *dscr, int dlen);
      
      /**
       * Steals one or more tasks from the given proc's shared portion.
       *
       * @param[in] proc to steal from
       */
      bool steal(int proc);
      
      /**
       * Returns true if all of the tasks have been processed.
       */
      bool hasTerminated();
      
      /**
       * Returns the number of tasks in the entire local queue.
       */
      int numTasks() const {
        return *size_local + *size_shared;
      }
      
      /**
       * Runs the termination detector.
       */
      void td_progress();

  }; /* SplitQueue */
}; /* tascel */

#endif /*__tascel_SplitQueue_h__*/

