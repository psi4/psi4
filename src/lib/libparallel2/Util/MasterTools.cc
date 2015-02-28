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

#include "MasterTools.h"
#include "MPITaskQueue.h"
namespace psi{
namespace LibParallel{

void MasterTools::PrintOut()const{
   for(int i=0;i<Tasks_->size();i++){
      this->LibParallelBase::operator<<((*Tasks_)[i].Number())<<" ";
   }
   this->LibParallelBase::operator<<("\n");
}

bool MasterTools::AllTasksComplete()const{
   return(CurrentTask_==Tasks_->size());
}

bool MasterTools::Done(){
   //Can't release myself
   bool AllReleased=(Released_.size()==(State_->NProc()-1));
   return AllReleased;
}

MasterTools::MasterTools(boost::shared_ptr<const Communicator> State):
      CurrentTask_(0),TaskChunk_(1),Tasks_(new std::vector<MPITaskGuts>()),
      ToolSet(State),TaskMap_(State->NProc()-1){
}

void MasterTools::GiveTasks(){
   if(!this->Done()){
      int sender=State_->Probe(ANY,NEXT);
      //Actually get the message...
      State_->Receive(sender,NEXT,(int*)NULL,0);
      if(!AllTasksComplete()){
         int TasksLeft=Tasks_->size()-CurrentTask_;
         int size=(TasksLeft>=TaskChunk_?TaskChunk_:TasksLeft);
         State_->Send(sender,TASKSIZE,&size,1);
         std::vector<int> buffer(size);
         for(int i=0;i<size;i++){
            int AbsTask=(*Tasks_)[CurrentTask_].Number();
            ///Minus one because the first slave is MPI process 1
            TaskMap_.AddTask(CurrentTask_++,AbsTask,sender-1);
            buffer[i]=AbsTask;
         }
         State_->Send(sender,TASKS,&buffer[0],size);
      }
      else{
         int Done=0;
         State_->Send(sender,TASKSIZE,&Done,1);
         Released_.push_back(sender);
      }
   }
}

MasterTools& MasterTools::operator<<(const MPITaskGuts& Task){
   Tasks_->push_back(Task);
   return *this;
}

}}

