#ifndef __tascel_SharedQueue_h__
#define __tascel_SharedQueue_h__

#include "TerminationDetector.h"

namespace tascel {

  /**
   * A set of task queues shared among processes.
   *
   * Each process has a local queue with a maximum number of possible tasks.
   * A process can steal tasks from any other process's queue and have
   * tasks stolen from its own queue.
   *
   * This implementation forces the entire queue on a given proc to be locked
   * during any modification by itself or by any other proc.
   */
  class SharedQueue {
    private:
      /**
       * The local queue state which other procs can modify.
       */
      struct sq_state_t {
        int dirty; /**< whether a steal has occurred */
        int  head; /**< the head of the queue */
        int  tail; /**< the tail of the queue */
        int  size; /**< the number of tasks in the queue */
      };

      const int max_ntsks;    /**< max number of tasks allowed in queue */
      const int tsk_size;     /**< size of a single task */
      char **q;               /**< addresses of queue data on all procs */
      sq_state_t **sq_state;  /**< addresses of queue state on all procs */
      int *head;              /**< pointer to this procs head state */
      int *tail;              /**< pointer to this procs tail state */
      int *size;              /**< pointer to this procs size state */
      int *dirty;             /**< pointer to this procs dirty state */
      TerminationDetector td; /**< the TerminationDetector */

    public:
      /**
       * Constructs the SharedQueue instance.
       *
       * @param[in] tsk_size size of a single task
       * @param[in] max_ntsks max number of tasks in a queue
       *
       * @pre tsk_size must be the same value on all procs
       * @pre max_ntsks must be the same value on all procs
       */
      SharedQueue(int tsk_size, int max_ntsks);
      
      /**
       * Destroys the SharedQueue instance.
       */
      ~SharedQueue();
      
      /**
       * Returns true if the local queue is empty.
       *
       * @return true if the local queue is empty.
       */
      bool empty() const;
      
      /**
       * Retrieves a task from the local queue.
       *
       * @param[out] dscr the retrieved task
       * @param[in] dlen size of the retrieved task
       *
       * @pre dscr is not NULL
       * @pre dlen == tsk_size
       *
       * @return true if task was successfully retrieved
       */
      bool getTask(void *dscr, int dlen);
      
      /**
       * Adds a task to the local queue.
       *
       * @param[in] dscr the task to add
       * @param[in] dlen size of the added task
       *
       * @pre dscr is not NULL
       * @pre dlen == tsk_size
       *
       * @return true if task was successfully added
       */
      void addTask(void *dscr, int dlen);
      
      /**
       * Steals one or more tasks from the given proc.
       *
       * @param[in] proc to steal from
       */
      bool steal(int proc);
      
      /**
       * Returns true if all of the tasks have been processed.
       *
       * @return true if all of the tasks have been processed.
       */
      bool hasTerminated();
      
      /**
       * Runs the termination detector.
       */
      void td_progress();

  }; /* SharedQueue */
}; /* tascel */

#endif /*__tascel_SharedQueue_h__*/

