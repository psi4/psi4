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
#ifndef MPISTATICSCHEDULER_H_
#define MPISTATICSCHEDULER_H_

#include "MPIScheduler.h"
namespace psi{
namespace LibParallel{
class MPITaskQueue;
class Communicator;

/** \brief The abstract base class for a static MPI scheduler
 *
 *  Static MPI schedulers assign all tasks at object creation.
 *  Thus which tasks each process is going to perform is already known.
 */
class StaticScheduler: public MPIScheduler{
   protected:
      ///The task queue the current MPI process is responsible for
      boost::shared_ptr<MPITaskQueue> Tasks_;

      ///This state is made here, and puts all comms on their own communicator
      boost::shared_ptr<Communicator> FState_;

   public:

      ///Splits IState_ into FState_
      StaticScheduler(boost::shared_ptr<const Communicator> IState);

      ///Frees FState
      ~StaticScheduler();

      ///Static schedulers may not add tasks as they go, so it's done
      bool Done();

      ///Simply returns our TaskQueue, again it's in its finished state
      boost::shared_ptr<MPITaskQueue> GetTaskQueue();
};

}}




#endif /* MPISTATICSCHEDULER_H_ */
