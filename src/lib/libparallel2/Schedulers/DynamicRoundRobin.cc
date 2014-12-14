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

#include "DynamicRoundRobin.h"
#include "../TaskJobGuts/MPITaskGuts.h"
namespace psi{
namespace LibParallel{
bool DynamicRoundRobin::Done(){
   bool done=DynamicScheduler::Done();
   if(done)Final_->FreeComm();
   return done;
}
void DynamicRoundRobin::FillMasterQueues(
      boost::shared_ptr<std::vector<MPITaskGuts> >& Tasks){
   NTasks_=Tasks->size();
   Remainder_=NTasks_%NMasters_;
   EvenTasks_=NTasks_-Remainder_;
   ///Load the Tasks
   if(AmMaster_){
      boost::shared_ptr<MasterTools> temp=
            boost::dynamic_pointer_cast<MasterTools>(ToolKit_);
      for(int task=0;task<EvenTasks_;){
         for(int i=0;i<NMasters_;i++){
            TaskMap_.AddTask(task,(*Tasks)[task].Number(),i);
            if(IState_->Me()==MasterIDs_[i]){
               (*temp)<<(*Tasks)[task];
            }
            task++;
         }
      }
      for(int task=EvenTasks_;task<NTasks_;task++)
         if(IState_->Me()==MasterIDs_[task-EvenTasks_])
            (*temp)<<(*Tasks)[task];
   }
}
void DynamicRoundRobin::PrintOut()const{
   for(int i=0;i<NMasters_;i++){
      if(IState_->Me()==MasterIDs_[i]){
         (*this)<<"I am Master "<<i<<" and am responsible for:\n";
         ToolKit_->PrintOut();
      }
   }
}

DynamicRoundRobin::DynamicRoundRobin(
      boost::shared_ptr<std::vector<MPITaskGuts> >& Tasks,
      boost::shared_ptr<const Communicator>& IState,const int NMasters):
            DynamicScheduler(IState,NMasters),TaskMap_(NMasters){
   Alg_=DYNAMICRR;
   FillMasterQueues(Tasks);
   if(!AmMaster_)ToolKit_->GiveTasks();
   Final_=IState->MakeComm(IState_->Me());
}


}}///End namespaces
