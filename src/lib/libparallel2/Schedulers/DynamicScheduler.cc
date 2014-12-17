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

#include "DynamicScheduler.h"
#include "../Util/MasterTools.h"
#include "../Communicator.h"
#include "../Util/MPITaskQueue.h"
#include "../Util/SlaveTools.h"
namespace psi{
namespace LibParallel{

DynamicScheduler::DynamicScheduler(
      boost::shared_ptr<const Communicator>& IState,
      const int NMasters):
      MPIScheduler(IState),NMasters_(NMasters),AmMaster_(false){
   int NProcs=IState->NProc();
   int Me=IState->Me();
   //Each Master needs a slave
   if(NProcs/NMasters_<2.0)
      Error("Not enough processes for a dynamic scheduler.");
   //Number of extra processes
   int Remainder=NProcs%NMasters_;
   //Number of processes per master/slave team
   int BatchSize=(NProcs-Remainder)/NMasters_;
   int MyBatch=-1;
   for(int i=0;i<NMasters_&&MyBatch==-1;i++)
      if(Me<i*BatchSize)MyBatch=i;
   Comm3_=IState->MakeComm(MyBatch);
   AmMaster_=(Comm3_->Me()==0);
   //Now let each process know who the Masters are
   SynchMasters();
   if(AmMaster_)ToolKit_=boost::shared_ptr<MasterTools>(new MasterTools(Comm3_));
   else ToolKit_=boost::shared_ptr<SlaveTools>(new SlaveTools(Comm3_));
}

void DynamicScheduler::PrintOut()const{
   MPIScheduler::PrintOut();
   std::cout<<"Dynamic Scheduler using: "<<NMasters_<<" masters, which are:\n";
   for(int i=0;i<NMasters_;i++)std::cout<<MasterIDs_[i]<<" ";
   std::cout<<"Comm3 is:"<<std::endl;
   Comm3_->PrintOut();
}

void DynamicScheduler::SynchMasters(){
   int EffMasters=NMasters_;
   if(AmMaster_){
      if(IState_->Me()==0){
         //We know one of the masters
         EffMasters--;
         MasterIDs_.push_back(0);
      }
      else{
         int value2send=IState_->Me();
         IState_->Send(0,1,&value2send,1);
      }
   }
   if(IState_->Me()==0){
      for(int i=0;i<EffMasters;i++){
         int value2recv=0;
         IState_->Receive(-1,-1,&value2recv,1);
         MasterIDs_.push_back(value2recv);
      }
   }
   else MasterIDs_.reserve(NMasters_);
   IState_->Bcast(&MasterIDs_[0],NMasters_,0);
}

boost::shared_ptr<MPITaskQueue> DynamicScheduler::GetTaskQueue(){
   boost::shared_ptr<MPITaskQueue> Return;
   if(AmMaster_)
      Return=boost::shared_ptr<MPITaskQueue>(new MPITaskQueue());
   else
      Return=boost::dynamic_pointer_cast<SlaveTools>(ToolKit_)->Tasks();
   return Return;
}

bool DynamicScheduler::Done(){
   bool AmDone=false;
   if(AmMaster_){
      while(!ToolKit_->Done()){
         ToolKit_->GiveTasks();
      }
      AmDone=true;
   }
   else AmDone=ToolKit_->Done();
   if(AmDone)Comm3_->FreeComm();
   return AmDone;
}
}}//End namespaces

