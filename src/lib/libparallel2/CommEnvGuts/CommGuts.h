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
#ifndef SRC_LIB_LIBPARALLEL2_COMMENVGUTS_COMMGUTS_H_
#define SRC_LIB_LIBPARALLEL2_COMMENVGUTS_COMMGUTS_H_
#include "../LibParallelBase.h"
#include "../Algorithms.h"
#include <boost/shared_ptr.hpp>
#include <string>
#if HAVE_MPI
  #include "ParallelComm.h"
#else
  #include "LocalComm.h"
#endif
namespace psi{
namespace LibParallel{
namespace Comm{
#if HAVE_MPI
    typedef ParallelComm CommType_;
#else
    typedef LocalComm CommType_;
#endif
}
class ParallelEnvironmentGuts;
class Communicator;
class CommGuts: public LibParallelBase{
   private:
      boost::shared_ptr<Comm::CommType_> DaComm_;
      ParallelEnvironmentGuts* Env_;
      std::string Name_;
      void Copy(const CommGuts& other);
      bool Active_;
   protected:
      CommGuts(const std::string& Name,ParallelEnvironmentGuts* Env);
      CommGuts(const CommGuts& other);
      const CommGuts& operator=(const CommGuts& other);
      void RegisterComm(boost::shared_ptr<Communicator> Comm)const;
   public:
      virtual ~CommGuts();
      virtual bool Active()const{return Active_;}
      virtual void FreeComm();
      std::string Name()const{return Name_;}

      /** \brief Makes a communicator
       *
       *   The user ultimately needs a Communicator object back, not
       *   a CommGuts object.  The layer above this class allocates the
       *   object, and then passes it to this class to have it set-up.
       *   That's the first argument, the second is how we actually split
       *   it.
       */
      boost::shared_ptr<CommGuts> MakeComm(const std::string& Name,
            const int Color)const;

      virtual void Barrier()const;

      virtual int Probe(const int Sender, const int MessageTag,
                        const bool Block)const;

      virtual int Me()const;

      virtual int NProc()const;

      template<typename T>
      void AllGatherImpl(const T* LocalData, const int NElem,
            T* Target)const{
         DaComm_->AllGather(LocalData,NElem,Target);
      }

      template<typename T>
      void AllReduceImpl(const T* LocalData, const int NElem,T* Target,
                     const MPIOperation& Op)const{
         DaComm_->AllReduce(LocalData,NElem,Target,Op);
      }

      template<typename T>
      void GatherImpl(const T* LocalData, const int NElem,T* Target,
            const int Root)const{
         DaComm_->Gather(LocalData,NElem,Target,Root);
      }

      template<typename T>
      void BcastImpl(T* Data, const int NElem, const int Broadcaster)const{
         DaComm_->Bcast(Data,NElem,Broadcaster);
      }

      template<typename T>
      void SendImpl(const int Receiver,const int MessageTag,T* Message=NULL,
                  const int Length=0,bool Block=true)const{
         DaComm_->Send(Receiver,MessageTag,Message,Length,Block);
      }

      template<typename T>
      void ReceiveImpl(const int Sender, const int MessageTag,T* Message=NULL,
            const int Length=0,bool Block=true)const{
         DaComm_->Receive(Sender,MessageTag,Message,Length,Block);
      }
};


}}//End namespaces



#endif /* SRC_LIB_LIBPARALLEL2_COMMENVGUTS_COMMGUTS_H_ */