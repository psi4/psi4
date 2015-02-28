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
#ifndef MPISETTERUPPER_H_
#define MPISETTERUPPER_H_

#include <vector>
#include <boost/shared_ptr.hpp>

#include "../Algorithms.h"
#include "../MPITask.h"
#include "../TaskJobGuts/MPITaskGuts.h"

namespace psi{
namespace LibParallel{
class MPIScheduler;
class Communicator;
/** \brief The class charged with turning an array of MPITask<T>s into
 *       a working MPIJob.
 *
 *  When it makes a scheduler you may assume your tasks are already sorted
 */
class MPISetterUpper{
   private:

      ///The algorithm this job will be using
      SchedAlgorithm ChoosenAlg_;

      ///The state this SetterUpper is working on
      boost::shared_ptr<const Communicator> State_;

      /** \brief The tasks the user wants done, stripped down to
       * std::pairs<int,int> of the form (priority,original task number)
       *
       */
      boost::shared_ptr<std::vector<MPITaskGuts> > Tasks_;

      ///The scheduler
      boost::shared_ptr<MPIScheduler> Sched_;

      ///Sorts the Tasks in descending order or priority
      void SortTasks();

      ///Chooses the algorithm
      void ChooseAlgorithm(bool ForceDynamic);

   public:
      ///Changes your tasks to MPITaskGuts, and determines the best Algorithm
      template<typename T>
      MPISetterUpper(const std::vector<MPITask<T> >& Tasks,
            boost::shared_ptr<const Communicator> IState,
            bool ForceDynamic);

      /** \brief Based on the value of ChoosenAlg_ returns your scheduler
       *
       *   If the current MPI process is not active you will get a NULL
       *   pointer back, so make sure you check it.
       */
      boost::shared_ptr<MPIScheduler> GetScheduler()const;



};

template<typename T>
MPISetterUpper::MPISetterUpper(const std::vector<MPITask<T> >& Tasks,
      boost::shared_ptr<const Communicator> IState,
      bool ForceDynamic):Tasks_(new std::vector<MPITaskGuts>),
      State_(IState){
      for(int i=0;i<Tasks.size();i++)
         Tasks_->push_back(MPITaskGuts(Tasks[i].Priority(),i));
      this->ChooseAlgorithm(ForceDynamic);
}

}}//End namespaces



#endif /* MPISETTERUPPER_H_ */
