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

#include "TaskQueue.h"
#include "psi4-dec.h"
#include "MPIProperty.h"
namespace psi{

void SimpleScheduler::SetValues(){
   Remainder_=NTasks_%NProcs_;
   BatchSize_=(NTasks_-Remainder_)/NProcs_;
   Me_=WorldComm->me();
   MyStart_=Me_*BatchSize_;
   MyEnd_=(Me_+1)*BatchSize_;
}

boost::shared_ptr<MPIQueue> SimpleScheduler::Schedule(){
   SetValues();
   std::vector<boost::shared_ptr<MPITask> > TempQueue;
   for(int i=MyStart_;i<MyEnd_;i++){
      boost::shared_ptr<MPITask> Task(new MPITask(i));
      TempQueue.push_back(Task);
   }
   if(Me_<Remainder_){
      int AdditionalTask=NProcs_*BatchSize_+Me_;
      boost::shared_ptr<MPITask> Task(new MPITask(AdditionalTask));
      TempQueue.push_back(Task);
   }
   boost::shared_ptr<MPIQueue> Queue(new MPIQueue(TempQueue));
   return Queue;
}

template<typename T>
boost::shared_ptr<MPIQueue> SimpleScheduler::SynchImpl(std::vector<MPIProperty<T> >)


}

