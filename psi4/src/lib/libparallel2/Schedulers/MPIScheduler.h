/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */
#ifndef MPISCHEDULER_H_
#define MPISCHEDULER_H_
#include <boost/shared_ptr.hpp>
#include "../Algorithms.h"
#include "../LibParallelBase.h"
namespace psi {
namespace LibParallel {
class MPITaskQueue;
class Communicator;
/** \brief This class is responsible for assigning tasks
 *
 *   All schedulers are predicated on the calls occurring in the order:
 *   Begin, Done, Next, Done, Next, ..., Next , Done.  This is the
 *   order both a for loop and a while loop (with begin right before the
 *   loop, Done in the termination condition, and next at the end) will
 *   call them in.  Note that these Begin, Done, and Next calls emminate
 *   from the MPIJobGuts class.
 *
 *   For future schedulers we assume several things.  First, Begin
 *   must return an integer that will not segfault when mapped to
 *   a label.  Currently there are two reasons why Begin may fail to do this.
 *   1)The current MPI process is not active, and thus doesn't have a
 *   queue or 2) That queue is empty.  In both cases we solve this by
 *   returning 0, which is a valid label, except if the user has given us
 *   no tasks to parallelize, a different problem entirely...
 *
 *   The Done call will first check if the process is active, if it isn't
 *   that process leaves the loop at this point.   If it's active we check
 *   if the queue is Done (true if the last task it gave away is equal to
 *   the last task it holds).  Namely a queue is done immediately if it is
 *   empty.  We exploit this for letting the Master process get lost in
 *   the Done call (a master process should always have an empty queue).
 *   Next we check if the scheduler agrees that the process is done.  For
 *   static schedulers, this condition is always true.  For dynamic
 *   schedulers, this condition is true only if the ToolKit says it's done.
 *   The toolkit says a Master process is done if it has given out all
 *   of it's tasks, and told each process it is free.  Slaves are free
 *   only after the Master has told them they are.  The Master process
 *   spools in the done function until it is done, at which point it returns
 *   true.
 *
 *   The Next call is pretty simple
 *
 */
class MPIScheduler:public LibParallelBase {
   protected:
      ///This is the state before this scheduler started doing anything
      boost::shared_ptr<const Communicator> IState_;

      ///This is the algorithm we actually are
      SchedAlgorithm Alg_;
   public:
      ///Returns the algorithm we are using
      SchedAlgorithm Algorithm() const;

      ///Sets IState_ to State
      MPIScheduler(boost::shared_ptr<const Communicator> State);

      ///Nothing to free up
      virtual ~MPIScheduler() {}

      ///Returns a task queue
      virtual boost::shared_ptr<MPITaskQueue> GetTaskQueue()=0;

      ///Returns true if the current MPIProcess is free to go
      virtual bool Done()=0;

      ///Debugging print
      virtual void PrintOut() const;

};
}} //End namespaces

#endif /* MPISCHEDULER_H_ */