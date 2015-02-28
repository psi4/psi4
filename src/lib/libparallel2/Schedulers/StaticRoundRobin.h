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
#ifndef SRC_LIB_LIBPARALLEL2_STATICROUNDROBIN_H_
#define SRC_LIB_LIBPARALLEL2_STATICROUNDROBIN_H_

#include <vector>
#include "../Communicator.h"
#include "MPIStaticScheduler.h"
#include "../Util/TaskMap.h"

namespace psi{
namespace LibParallel{
class MPITaskGuts;

class StaticRoundRobin: public StaticScheduler{
   private:
      boost::shared_ptr<TaskMap> Map_;
      int NTasks_;
   public:
      StaticRoundRobin(boost::shared_ptr<std::vector<MPITaskGuts> > Tasks,
            boost::shared_ptr<const Communicator>& State);
      template<typename T>
      std::vector<T> ReduceImpl(const std::vector<T>& LocalValues,
            const int N,
            const MPIOperation& op);
      template<typename T>
      std::vector<T> SynchImpl(const std::vector<T>& LocalValues,
            const int N);
};

template<typename T>
std::vector<T> StaticRoundRobin::ReduceImpl(
      const std::vector<T>& LocalValues,
      const int N,const MPIOperation& op){
   std::vector<T> Result(N);
   IState_->AllReduce(&LocalValues[0],N,&Result[0],op);
   return Result;
}

template<typename T>
std::vector<T> StaticRoundRobin::SynchImpl(
      const std::vector<T>& LocalValues,
      const int N){
   int TotSize=N*NTasks_;
   std::vector<T> Result(TotSize);
   int NProc=IState_->NProc();
   int Remainder=NTasks_%NProc;
   int Even=(NTasks_-Remainder)/NProc;
   int Me=IState_->Me();
   IState_->AllGather(&LocalValues[0],Even*N,&Result[0]);
   if(Remainder>0){
      for(int i=0;i<N;i++)Result[NProc*N*Even+i]=LocalValues[N*Even+i];
   }
   for(int i=1;i<Remainder;i++){
      if(Me==0){
         IState_->Receive(i,1,&Result[NProc*N*Even+i*N],N);
      }
      if(Me==i){
         IState_->Send(0,1,&LocalValues[N*Even],N);
      }
   }
   std::vector<T> Temp(TotSize);
   if(Me==0){
      for(int AbsTask=0;AbsTask<NTasks_;AbsTask++){
         int Owner=Map_->Slave(AbsTask);
         int index=Map_->Index(AbsTask);
         int offset=(index<Even?Owner*N*Even+index*N:NProc*N*Even+Owner*N);
         for(int i=0;i<N;i++){
            Temp[AbsTask*N+i]=Result[offset+i];
         }
      }
   }
   IState_->Bcast(&Temp[0],TotSize,0);
   Result=Temp;
   return Result;
}

}}//End namespaces



#endif /* SRC_LIB_LIBPARALLEL2_STATICROUNDROBIN_H_ */
