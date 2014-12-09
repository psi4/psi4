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

#include "MPIJobGuts.h"
#include "../Schedulers/MPIScheduler.h"
#include "psi4-dec.h"
namespace psi{
namespace LibParallel{

MPIJobGuts::~MPIJobGuts(){
   //Let the scheduler clean-up int's communicators first
   Scheduler_.reset();
   if(FState_){
      if(FState_->Active())FState_->FreeComm();
   }
}

bool MPIJobGuts::EnoughProcs(const int Procs)const{
   bool enough=!(Procs<0);
   if(enough){
      int NProcs=IState_->NProc();
      enough=(Procs<NProcs);
   }
   return enough;
}

bool MPIJobGuts::EnoughProcs()const{return Enough_;}

int MPIJobGuts::Begin(){
   return (Active_&&Queue_->size()>0?Queue_->Begin():0);
}

int MPIJobGuts::Next(){
   int next=0;
   if(Active_){
      //If the queue returns 0 we are done with our queue
      next=Queue_->Next();
      if(next==-1){
         //See if scheduler wants to update our queue
         if(!Scheduler_->Done())next=Queue_->Next();
         else next=0;
      }
   }
   return next;
}

bool MPIJobGuts::Done(){
   bool QD=false,SD=false;
   if(Active_){
      //True if queue is done
      QD=Queue_->Done();
      //If the queue is done see if the scheduler is done
      if(QD)SD=(Scheduler_?Scheduler_->Done():true);
      if(QD && SD && FState_)FState_->FreeComm();
   }
   return (Active_?(QD && SD):true);
}

boost::shared_ptr<Communicator> MPIJobGuts::InitialState(const int MaxProcs) {
   int Me=IState_->Me();
   boost::shared_ptr<Communicator> state;
   int MyValue=Me;
   if (MaxProcs==0)Error(
         "I don't know what it means to restrict the "
               "number of MPI processes to 0");
   else if (MaxProcs>0) {
      MyValue=(Me<MaxProcs ? 0 : MaxProcs);
      if (MyValue==MaxProcs) {
         //These MPI processes aren't actually participating
         Active_=false;
      }
      state=IState_->MakeComm(MyValue);
   }
   return state;
}

void MPIJobGuts::Wait()const{
   (FState_?FState_->Barrier():IState_->Barrier());
}
}}

