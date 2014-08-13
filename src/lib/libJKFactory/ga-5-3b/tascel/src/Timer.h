#ifndef __tascel_Timer_h__
#define __tascel_Timer_h__

#include <mpi.h>

namespace tascel {
  class Timer {
  private:
    double ttime;   /**< Total time so far  */
    double stime;   /**< Time at last start */

    /**
     * Current time. A simple (portable) implementation for now.
     *
     * @return Current time stamps (in seconds)
     */
    inline double timestamp() {
      return MPI_Wtime();
    }

  public:
    
    /**
     * Constructor
     */
  Timer() : ttime(0) {}

    /**
     * Start the timer. Overwrites previous start
     */
    inline void startTimer() { 
      stime = timestamp(); 
    }

    /**
     * Stop timer. Adds time since previous start to internal state.
     */
    inline void stopTimer() { 
      double etime_ = timestamp(); 
      ttime += etime_ - stime;      
    }
    
    /**
     * Returns sum of all times in intervals between start and stop of
     * times.
     */
    inline double time() const { 
      return ttime; 
    }
  }; /* Timer */
}; /* tascel */

#endif /*__tascel_Timer_h__*/

