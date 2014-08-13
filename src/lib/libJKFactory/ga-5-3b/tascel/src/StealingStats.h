#ifndef __tascel_StealingStats_h__
#define __tascel_StealingStats_h__

#include "Counter.h"
#include "Timer.h"

namespace tascel {
  class StealingStats {
  public:
    Counter numStealAttempts;   /**< Num. steal attempts */
    Counter numSteals;          /**< Num steals that got work */
    Counter numTasks;           /**< Num. tasks executed by this proc*/

    Timer stealTime;            /**< Time spent stealing */
    Timer taskTime;             /**< Time spend executing tasks */
    Timer tdTime;               /**< Time spent in termination detection */

    /**
     * Print all stats
     */
    void print() const;

  }; /* StealingStats */
}; /* tascel */

#endif /* __tascel_StealingStats_h__ */

