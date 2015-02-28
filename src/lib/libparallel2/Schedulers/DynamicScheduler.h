/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */
#ifndef DYNAMICSCHEDULER_H_
#define DYNAMICSCHEDULER_H_

#include <vector>
#include "MPIScheduler.h"

namespace psi{
namespace LibParallel{
class ToolSet;
class Communicator;
/** \brief The base class for algorithms that are a master/slave model
 *
 *   For dynamic schedulers our MPI comm tree starts with the communicator
 *   that created the MPIJob, which we call Comm1.  Comm1 is then split into
 *   two communicators, iff the user won't let us use all the processes. We
 *   call that Comm2, realizing it may be equal to Comm1.  Management of
 *   Comm2 is done by MPIJob (Comm1 proceeded our object creation and is
 *   not our concern).  The first step of a dynamic scheduler is to determine
 *   the number of master processes (currently hard-coded to 1, but set-up
 *   to take this into consideration in the future).  If we have say three
 *   masters, then Comm2 is split into 3 comms, each with a master, and an
 *   approximately equal amount of slaves.  This comm is Comm3.  Now the
 *   scheduler is called and it splits each process in Comm3 into its own
 *   comm, called Comm4.  Scheduler is responsible for Comm4, this class
 *   is responsible for Comm3.
 */
class DynamicScheduler:public MPIScheduler{
   private:
      ///Function that makes sure each process knows who the Masters are
      void SynchMasters();

   protected:
      ///The number of master processes
      int NMasters_;

      ///The identities of the masters on Comm2
      std::vector<int> MasterIDs_;

      ///True if the current MPI process is a Master Task
      bool AmMaster_;

      /** \brief This is the MPIState we make by splitting IState_ (Comm2)
       *
       *   Comm2_ is all the processes we were allowed by the user. Comm3_
       *   is Comm2_ divided into NMaster_ groups of processes, as evenly
       *   as possible.  In the event that NMaster_=1, Comm3_ is effectively
       *   Comm2_, but derived from it, i.e. we need to free Comm3_ when
       *   finished.
       */
      boost::shared_ptr<Communicator> Comm3_;

      ///These are little operations needed to facilitate master/slave tasks
      boost::shared_ptr<ToolSet> ToolKit_;

   public:
      ///Passes IState (Comm2) to base class, then sets up Comm3
      DynamicScheduler(boost::shared_ptr<const Communicator>& IState,
            const int NMasters=1);

      /** \brief Returns a task queue
       *
       *  If the MPI process is a Master then this function will in turn
       *  return a queue with size 0.  Otherwise it returns a pointer
       *  to the queue where your future tasks will be dumped, don't
       *  loose it...(if you do just call this function again)
       */
      boost::shared_ptr<MPITaskQueue> GetTaskQueue();

      /** \brief This is where the magic happens
       *
       */
      bool Done();

      ///Adds to debugging printing
      void PrintOut()const;
};

}}//End namespaces


#endif /* DYNAMICSCHEDULER_H_ */
