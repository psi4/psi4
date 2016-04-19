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

#include "SlaveTools.h"
#include "../Communicator.h"
#include "MPITaskQueue.h"
namespace psi{
namespace LibParallel{

void SlaveTools::PrintOut()const{
   for(int i=0;i<MyTasks_->size();i++){
      (*this)<<(*MyTasks_)[i]<<" ";
   }
   (*this)<<"\n";
}

SlaveTools::SlaveTools(boost::shared_ptr<const Communicator> State):
      ToolSet(State),MyTasks_(new MPITaskQueue()),Released_(false){}

bool SlaveTools::Done(){
   if(MyTasks_->Done()&&!Released_)GiveTasks();
   return MyTasks_->Done();
}

void SlaveTools::GiveTasks(){
   State_->Send(0,NEXT,(int*)NULL,0);
   int ntasks=0;
   State_->Receive(0,TASKSIZE,&ntasks,1);
   if(ntasks!=0){
      std::vector<int> Tasks(ntasks);
      State_->Receive(0,TASKS,&Tasks[0],ntasks);
      for(int i=0;i<ntasks;i++)(*MyTasks_)<<Tasks[i];
   }
   else Released_=true;
}

boost::shared_ptr<MPITaskQueue> SlaveTools::Tasks(){return MyTasks_;}

}}
