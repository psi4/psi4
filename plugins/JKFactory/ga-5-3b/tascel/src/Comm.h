#ifndef __tascel_Comm_h__
#define __tascel_Comm_h__

namespace tascel {

  /**
   * Encapsulation of common communication functions.
   *
   * Much of the communication operations will be abstracted here. At
   * this point, the operations are not group-aware.  
   */
  namespace comm {

    /** 
     * barrier (collective).
     *
     * Equivalent of MPI's MPI_Barrier(MPI_COMM_WORLD)
     */
    void barrier();

    /**
     * Number of processes.
     *
     * Equivalent of nproc in MPI_Comm_size(MPI_COMM_WORLD, &nproc)
     */
    int nproc();

    /**
     * Rank of invoking process.
     *
     * Equivalent of me in MPI_Comm_rank(MPI_COMM_WORLD, &me)
     */
    int me();
  }; /*comm*/

}; /*tascel*/

#endif /*__tascel_Comm_h__*/

