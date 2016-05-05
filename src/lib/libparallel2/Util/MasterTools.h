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
#ifndef MASTERTOOLS_H_
#define MASTERTOOLS_H_

#include "../Communicator.h"
#include "../TaskJobGuts/MPITaskGuts.h"
#include "../Algorithms.h"
#include "TaskMap.h"
#include "ToolSet.h"
namespace psi{
namespace LibParallel{

class MasterTools:public ToolSet{
   private:
      ///The task we just gave out
      unsigned int CurrentTask_;

      ///The number of tasks we are giving away at a time
      unsigned int TaskChunk_;

      ///Our Task queue
      boost::shared_ptr<std::vector<MPITaskGuts> > Tasks_;

      ///A map of who I gave, what task to (0 is Tasks_[0].Number())
      TaskMap TaskMap_;

      ///A vector of the slaves I have released
      std::vector<int> Released_;

      bool AllTasksComplete()const;
   public:
      ///Debugging printing
      void PrintOut()const;
      ///Implements the function responsible for assigning tasks
      void GiveTasks();
      MasterTools& operator<<(const MPITaskGuts& Task);
      MasterTools(boost::shared_ptr<const Communicator> State);
      ///True if we have given out all tasks
      bool Done();


      template<typename T>
      std::vector<T> Reorder(const std::vector<T>& UnOrder,const int N);

      template<typename T>
      std::vector<T> SynchImpl(const std::vector<T>& Local, const int N,
            boost::shared_ptr<const Communicator> Comm3_);

      template<typename T>
      std::vector<T> ReduceImpl(const std::vector<T>& Local, const int N,
            const MPIOperation& op,boost::shared_ptr<const Communicator> Comm3_);
};

template<typename T>
std::vector<T> MasterTools::ReduceImpl(
      const std::vector<T>& Local,
      const int N,
      const MPIOperation& op,
      boost::shared_ptr<const Communicator> Comm3_){
   std::vector<T> Results(N);
   boost::shared_ptr<Communicator> Comm=Comm3_->MakeComm(true);
   Comm3_->Receive(1,1,&Results[0],N);
   Comm->FreeComm();
   return Results;
}

template<typename T>
std::vector<T> MasterTools::SynchImpl(
      const std::vector<T>&,
      const int N,
      boost::shared_ptr<const Communicator> Comm3_){
   int NTasks=N*Tasks_->size();
   std::vector<T> Temp;
   for(int i=1;i<Comm3_->NProc();i++){
      int size=0;
      Comm3_->Receive(i,1,&size,1);
      std::vector<T> temp2(size);
      Comm3_->Receive(i,1,&temp2[0],size);
      Temp.reserve(Temp.size()+temp2.size());
      Temp.insert(Temp.end(),temp2.begin(),temp2.end());
   }
   //Now we have the entire array
   //Next task: get it back in the right order
   std::vector<T> Results=Reorder(Temp,N);
   Comm3_->Bcast(&NTasks,1,0);
   Comm3_->Bcast(&Results[0],NTasks,0);
   return Results;

}

template<typename T>
std::vector<T> MasterTools::Reorder(const std::vector<T>& UnOrder,
      const int N){
   int NTasks=Tasks_->size();
   std::vector<T> Results(NTasks*N);
   for(int i=0;i<NTasks;i++){
      int Task=(*Tasks_)[i].Number();
      int slave=TaskMap_.Slave(Task);
      int index=TaskMap_.Index(Task);
      int offset=TaskMap_.Offset(slave);
      for(int j=0;j<N;j++){
         int OverallI=i*N+j;
         int TempI=N*offset+N*index+j;
         Results[OverallI]=UnOrder[TempI];
      }
   }
   return Results;
}
}}



#endif /* MASTERTOOLS_H_ */