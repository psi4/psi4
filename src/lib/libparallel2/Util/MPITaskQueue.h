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
#ifndef MPITASKQUEUE_H_
#define MPITASKQUEUE_H_

#include <vector>
#include <boost/shared_ptr.hpp>

namespace psi{
namespace LibParallel{
class MPIScheduler;
/** \brief Class designed to hold a list of tasks a MPI processes must
 *         perform.
 *
 *   The MPIJob class is given one of these by the scheduler, on each
 *   MPI process.  When the user requests the next task of the MPIJob
 *   it checks with this object, which returns the next in its queue.
 *   When MPIJob Done checks for done status, this object returns false
 *   if it's queue is non-empty.  If it's queue is non-empty, it checks
 *   with the MPIScheduler to see if there are more tasks to perform.
 *
 */
class MPITaskQueue{
   private:
      ///The tasks I need to perform
      std::vector<int> Tasks_;

      ///The value of the last task I ran
      int LastTask_;

      ///The scheduler I report to
      boost::shared_ptr<MPIScheduler> Scheduler_;

      ///Makes Deep Copy of Tasks_ & LastTask_, shallow of Scheduler
      void Copy(const MPITaskQueue& other);
   public:
      MPITaskQueue();

      int size()const;

      ///Calls copy
      MPITaskQueue(const MPITaskQueue& other);

      ///Calls copy iff this!=other
      const MPITaskQueue& operator=(const MPITaskQueue& other);

      ///Combines two task queues
      void operator+=(const MPITaskQueue& other);

      ///Combines two task queues and returns the result
      MPITaskQueue operator+(const MPITaskQueue& other);

      ///Mechanism for inserting tasks into the queue
      MPITaskQueue& operator<<(const int NewTask);

      ///Returns a particular value of the queue, no updating of LastTask
      int operator[](const int Task)const;

      ///Returns the next task to perform, updates LastTask
      int Next();

      ///Returns the first task, updates LastTask, resets if necessary
      int Begin();

      ///Returns true if all tasks have been run
      bool Done()const;

};

}
}



#endif /* MPITASKQUEUE_H_ */
