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
#ifndef DYNAMICROUNDROBIN_H_
#define DYNAMICROUNDROBIN_H_

#include "DynamicScheduler.h"
#include "../Util/MasterTools.h"
#include "../Util/SlaveTools.h"
#include "../Util/TaskMap.h"
#include "../Algorithms.h"
#include "../Communicator.h"
namespace psi{
namespace LibParallel{
class MPITaskGuts;
class DynamicRoundRobin:public DynamicScheduler{
   private:
      int NTasks_;
      int Remainder_;
      int EvenTasks_;
      TaskMap TaskMap_;
      void FillMasterQueues(boost::shared_ptr<std::vector<MPITaskGuts> >&);

      ///A comm with each processor on its own communicator
      boost::shared_ptr<Communicator> Final_;
   public:
      ///Debug printing
      void PrintOut()const;

      ///Free Final
      bool Done();


      DynamicRoundRobin(boost::shared_ptr<std::vector<MPITaskGuts> >& Tasks,
                        boost::shared_ptr<const Communicator>& IState,
                        const int NMasters=1);
      template<typename T>
      std::vector<T> SynchImpl(const std::vector<T>& Local,const int N);

      template<typename T>
      std::vector<T> ReduceImpl(const std::vector<T>& Local,const int N,
            const MPIOperation& op);
};

template<typename T>
std::vector<T> DynamicRoundRobin::ReduceImpl(
      const std::vector<T>& Local,
      const int N,
      const MPIOperation& op){
   std::vector<T> Results(N);
   if(!AmMaster_)Results=boost::dynamic_pointer_cast<SlaveTools>
                           (ToolKit_)->ReduceImpl(Local,N,op,Comm3_);
   else Results=boost::dynamic_pointer_cast<MasterTools>
                  (ToolKit_)->ReduceImpl(Local,N,op,Comm3_);
   //We are synched across each master at this point
   if(NMasters_>1){
      boost::shared_ptr<Communicator> Comm=IState_->MakeComm(AmMaster_);
      if(AmMaster_){
         std::vector<T> temp(N);
         Comm->AllReduce(&Results[0],N,&temp[0],op);
         Results=temp;
      }
      Comm->FreeComm();
      IState_->Bcast(&Results[0],N,MasterIDs_[0]);
   }
   return Results;
}

template<typename T>
std::vector<T> DynamicRoundRobin::SynchImpl(const std::vector<T>& Local,
      const int N){
   std::vector<T> Results;
   if(!AmMaster_)Results=boost::dynamic_pointer_cast<SlaveTools>
                           (ToolKit_)->SynchImpl(Local,N,Comm3_);
   else Results=boost::dynamic_pointer_cast<MasterTools>
                           (ToolKit_)->SynchImpl(Local,N,Comm3_);
   //We are synched across each master at this point
   if(NMasters_>1){
     int NData=N*NTasks_;
     if(AmMaster_){
        if(IState_->Me()==MasterIDs_[0]){
           std::vector<std::vector<T> >temp;
           temp.push_back(Results);
           for(int i=1;i<MasterIDs_.size();i++){
              int size=0;
              IState_->Receive(MasterIDs_[i],1,&size,1);
              std::vector<T> temp2(size);
              IState_->Receive(MasterIDs_[i],1,&temp2[0],size);
              temp.push_back(temp2);
           }
           std::vector<T> Temp(NData);
           for(int i=0;i<NTasks_;i++){
                 int Task=i;
                 int slave=TaskMap_.Slave(Task);
                 int index=TaskMap_.Index(Task);
                 for(int j=0;j<N;j++){
                    int OverallI=i*N+j;
                    int TempI=N*index+j;
                    Temp[OverallI]=temp[slave][TempI];
                 }
              }
           Results=Temp;
        }
        else{
           int size=Results.size();
           IState_->Send(MasterIDs_[0],1,&size,1);
           IState_->Send(MasterIDs_[0],1,&Results[0],size);
        }
     }
     if(IState_->Me()!=MasterIDs_[0]){
        std::vector<T> Temp(NData);
        Results=Temp;
     }
     IState_->Bcast(&Results[0],NData,MasterIDs_[0]);
   }
   return Results;
}
}}//End namespaces



#endif /* DYNAMICROUNDROBIN_H_ */
