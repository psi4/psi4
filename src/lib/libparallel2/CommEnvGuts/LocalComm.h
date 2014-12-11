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
#ifndef SRC_LIB_LIBPARALLEL2_LOCALCOMM_H_
#define SRC_LIB_LIBPARALLEL2_LOCALCOMM_H_
#include "CommBase.h"
#include <boost/shared_ptr.hpp>
#include "../Algorithms.h"
namespace psi{
namespace LibParallel{

class LocalComm:public CommBase<LocalComm>{
   private:
      typedef CommBase<LocalComm> BaseType_;
      friend class CommBase<LocalComm>;
      template<typename T>
      void AllGatherImpl(const T* LocalData, const int NElem,
            T* Target)const{
         if(Target!=LocalData)memcpy(Target,LocalData,NElem*sizeof(T));
      }
      template<typename T>
      void AllReduceImpl(const T* LocalData, const int NElem,T* Target,
                     const MPIOperation& Op)const{
         if(Target!=LocalData)memcpy(Target,LocalData,NElem*sizeof(T));
      }
      template<typename T>
      void GatherImpl(const T* LocalData, const int NElem,T* Target,
            const int Root)const{
         if(Target!=LocalData)memcpy(Target,LocalData,NElem*sizeof(T));
      }
      template<typename T>
      void BcastImpl(T* /*Data*/, const int /*NElem*/, const int /*Broadcaster*/)const{
      }

      template<typename T>
      void SendImpl(const int /*Receiver*/,const int /*MessageTag*/,T* /*Message*/,
                  const int /*Length*/,bool /*Block*/)const{
      }

      template<typename T>
      void ReceiveImpl(const int /*Sender*/, const int /*MessageTag*/,T* /*Message*/,
            const int /*Length*/,bool /*Block*/)const{
      }
   public:
      LocalComm();
      boost::shared_ptr<LocalComm> MakeComm(const int Color)const;
      void Barrier()const{}
      int Probe(const int /*Sender*/,const int /*MessageTag*/,
            const bool /*Block*/)const{return 0;}
      int Me()const{return 0;}
      int NProc()const{return 1;}


};

}}//End namespaces


#endif /* SRC_LIB_LIBPARALLEL2_LOCALCOMM_H_ */
