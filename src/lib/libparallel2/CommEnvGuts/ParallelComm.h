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
#ifndef SRC_LIB_LIBPARALLEL_PARALLELCOMM_H_
#define SRC_LIB_LIBPARALLEL_PARALLELCOMM_H_
#include "CommBase.h"
#include "../Algorithms.h"
#include <boost/mpi.hpp>
#include <functional>
namespace psi {
namespace LibParallel {
class ParallelComm:public CommBase<ParallelComm> {
   private:
      typedef CommBase<ParallelComm> BaseType_;
      friend class CommBase<ParallelComm>;
      boost::shared_ptr<boost::mpi::communicator> Comm_;

      template <typename T>
      void AllGatherImpl(const T* LocalData, const int NElem, T* Target) const;

      template <typename T>
      void AllReduceImpl(const T* LocalData, const int NElem, T* Target,
            const MPIOperation& Op) const;

      template<typename T>
      void GatherImpl(const T* LocalData, const int NElem,T* Target,
            const int Root)const;

      template<typename T>
      void BcastImpl(T* Data, const int NElem, const int Broadcaster)const;

      template<typename T>
      void SendImpl(const int Receiver,const int MessageTag,T* Message=NULL,
                  const int Length=0,bool Block=true)const;

      template<typename T>
      void ReceiveImpl(const int Sender, const int MessageTag,T* Message=NULL,
            const int Length=0,bool Block=true)const;

   public:
      ParallelComm();
      boost::shared_ptr<ParallelComm> MakeComm(const int Color) const;
      void Barrier() const;
      int Probe(const int Sender, const int MessageTag, const bool Block) const;
      int Me()const;
      int NProc()const;
};

template <typename T>
void ParallelComm::AllGatherImpl(const T* LocalData, const int NElem,
      T* Target) const {
   boost::mpi::all_gather((*Comm_), LocalData, NElem, Target);
}

template<typename T>
void ParallelComm::GatherImpl(const T* LocalData, const int NElem,T* Target,
      const int Root)const{
   boost::mpi::gather((*Comm_),LocalData,NElem,Target,Root);
}

template<typename T>
void ParallelComm::BcastImpl(T* Data, const int NElem,
      const int Broadcaster)const{
   boost::mpi::broadcast((*Comm_),Data,NElem,Broadcaster);
}

template<typename T>
void ParallelComm::SendImpl(const int Receiver,const int MessageTag,
      T* Message,const int Length,bool Block)const{
   int SCopy=(Receiver<0?boost::mpi::any_source:Receiver);
   int TCopy=(MessageTag<0?boost::mpi::any_tag:MessageTag);
   int LCopy=Length;
   if(Message!=NULL && Length!=0){
      if(Length!=1){
         if(!Block)Comm_->isend(SCopy,TCopy,Message,LCopy);
         else Comm_->send(SCopy,TCopy,Message,LCopy);
      }
      else{
         if(!Block)Comm_->isend(SCopy,TCopy,Message[0]);
         else Comm_->send(SCopy,TCopy,Message[0]);
      }
   }
   else{
      if(!Block)Comm_->isend(SCopy,TCopy);
      else Comm_->send(SCopy,TCopy);
   }
}

template<typename T>
void ParallelComm::ReceiveImpl(const int Sender, const int MessageTag,
      T* Message,const int Length,bool Block)const{
   int SCopy=(Sender<0?boost::mpi::any_source:Sender);
   int TCopy=(MessageTag<0?boost::mpi::any_tag:MessageTag);
   int LCopy=Length;
   if(Message!=NULL && Length!=0){
      if(Length!=1){
         if(!Block)Comm_->irecv(SCopy,TCopy,Message,LCopy);
         else Comm_->recv(SCopy,TCopy,Message,LCopy);
      }
      else{
         if(!Block)Comm_->irecv(SCopy,TCopy,Message[0]);
         else Comm_->recv(SCopy,TCopy,Message[0]);
      }
   }
   else{
      if(!Block)Comm_->irecv(SCopy,TCopy);
      else Comm_->recv(SCopy,TCopy);
   }
}

template <typename T>
void ParallelComm::AllReduceImpl(const T* LocalData, const int NElem, T* Target,
      const MPIOperation& Op) const {
   switch (Op) {
      case (MULTIPLY): {
         boost::mpi::all_reduce((*Comm_), LocalData, NElem, Target,
               std::multiplies<T>());
         break;
      }
      case (DIVIDE): {
         boost::mpi::all_reduce((*Comm_), LocalData, NElem, Target,
               std::divides<T>());
         break;
      }
      case (ADD): {
         boost::mpi::all_reduce((*Comm_), LocalData, NElem, Target,
               std::plus<T>());
         break;
      }
      case (SUBTRACT): {
         boost::mpi::all_reduce((*Comm_), LocalData, NElem, Target,
               std::minus<T>());
         break;
      }
      /*Can't do modulus for doubles, don't know why no other compiler
       * complained till now...
       case (MODULUS): {
         boost::mpi::all_reduce((*Comm_), LocalData, NElem, Target,
               std::modulus<T>());
         break;
      }*/
      default:{
         throw PSIEXCEPTION("Unrecognized operation in mpiwrapper.h");
         break;
      }
   }
}

}} //End namespaces

#endif /* SRC_LIB_LIBPARALLEL_PARALLELCOMM_H_ */
