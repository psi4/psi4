#ifndef __tascel_TerminationDetector_h__
#define __tascel_TerminationDetector_h__

#include <mpi.h>

namespace tascel {

  /**
   * Termination Detector.
   */
  class TerminationDetector {
    private:
      int parent;           /**< id of this proc's parent */
      int child[2];         /**< id of this proc's children */
      int cvote[2];         /**< children's votes */
      MPI_Request creq[2];  /**< children's request handles */
      MPI_Request preq;     /**< parent's request handle */
      int myvote;           /**< value of this proc's vote */
      int pdecision;        /**< cumulative decision (myvote and vote of parent/children) */
      bool dirUp;           /**< voting direction */

      /**
       * Initiates receive of messages from children and resets vote to 'yes'.
       */
      void resetUpPhase();

      /**
       * Initiates recieve of message from parent.
       */
      void resetDownPhase();

      /**
       * Integrates votes from children and sends to parent.
       */
      void propogateUp();

      /**
       * Sends decision to children.
       */
      void propogateDown();

    public:
      /**
       * Constructs the TerminationDetector.
       */
      TerminationDetector();

      /**
       * Destroys the TerminationDetector.
       */
      ~TerminationDetector();

      /**
       * Assigns the value of myvote.
       */
      void vote(bool val);

      /**
       * Propogates votes up or down depending on the current mode.
       */
      void progress(bool vote);

      /**
       * Returns true if termination has been decided.
       *
       * @return true if termination has been decided.
       */
      bool hasTerminated();

  }; /*TerminationDetector*/

}; /*tascel*/

#endif /*__tascel_TerminationDetector_h__*/

