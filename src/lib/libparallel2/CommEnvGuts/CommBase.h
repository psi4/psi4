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
#ifndef SRC_LIB_LIBPARALLEL2_COMMBASE_H_
#define SRC_LIB_LIBPARALLEL2_COMMBASE_H_
#include <boost/shared_ptr.hpp>
#include "../LibParallelBase.h"
#include "../Algorithms.h"
namespace psi{
namespace LibParallel{

template<typename DerivedType>
class CommBase: public LibParallelBase{
   private:
      DerivedType* Derived_;
      typedef CommBase<DerivedType> MyType_;
   protected:
      CommBase(DerivedType* Derived):Derived_(Derived){}
   public:
      virtual ~CommBase(){}
      virtual boost::shared_ptr<DerivedType>
      MakeComm(const int Color)const{
         return Derived_->MakeComm(Color);
      }
      virtual void Barrier()const{Derived_->Barrier();}
      virtual int Probe(const int Sender,const int MessageTag,
                        const bool Block)const{
         return Derived_->Probe(Sender,MessageTag,Block);
      }
      virtual int Me()const{
         return Derived_->Me();
      }
      virtual int NProc()const{
         return Derived_->NProc();
      }
      template<typename T>
      void AllGather(const T* LocalData, const int NElem,
            T* Target)const{
         Derived_->AllGatherImpl(LocalData,NElem,Target);
      }
      template<typename T>
      void AllReduce(const T* LocalData, const int NElem,T* Target,
                     const MPIOperation& Op)const{
         Derived_->AllReduceImpl(LocalData,NElem,Target,Op);
      }
      template<typename T>
      void Gather(const T* LocalData, const int NElem,T* Target,
            const int Root)const{
         Derived_->GatherImpl(LocalData,NElem,Target,Root);
      }

      template<typename T>
      void Bcast(T* Data, const int NElem, const int Broadcaster)const{
         Derived_->BcastImpl(Data,NElem,Broadcaster);
      }

      template<typename T>
      void Send(const int Receiver,const int MessageTag,T* Message=NULL,
                  const int Length=0,bool Block=true)const{
         Derived_->SendImpl(Receiver,MessageTag,Message,Length,Block);
      }

      template<typename T>
      void Receive(const int Sender, const int MessageTag,T* Message=NULL,
            const int Length=0,bool Block=true)const{
         Derived_->ReceiveImpl(Sender,MessageTag,Message,Length,Block);
      }
};


}}



#endif /* SRC_LIB_LIBPARALLEL2_COMMBASE_H_ */
