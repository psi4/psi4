/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

/*
 * File:   parallel.h
 * Author: jturney, jjwilke
 *
 * Created on December 11, 2009, 3:34 PM
 */

#ifndef _psi_src_lib_libparallel_parallel_h_
#define	_psi_src_lib_libparallel_parallel_h_

#include <vector>
#include <string>
#include "psi4/libparallel/process.h"
#include "libparallel.h"

// Use this to silence "warning: unused parameter" messages from the compiler.
// Sometimes you have to define the parameter for Doxygen regardless if you
// use it or not.
#ifndef UNUSED
#define UNUSED(expr) (void)(expr)
#endif

namespace psi {

template <typename DerivedType>
class Parallel {
   public:
      ///This typedef is the type of this class
      typedef Parallel<DerivedType> ThisType;
   protected:
      ///Array of our communicator names, in the order they are derived
      std::vector<std::string> CurrentComm;
   public:
      ///Sets current comm to COMM_WORLD
      Parallel() {
         CurrentComm.push_back("COMM_WORLD");
      }

      ///No memory to free-up, does nothing
      virtual ~Parallel() {

      }

      ///Provides MPI barrier functionality
      virtual void sync(const std::string& CommName="NONE") const {
          UNUSED(CommName);
      }

      /** \brief Performs and all reduce
       *
       *   For each MPI process:
       *   \f[
       *     Target_i=\bigoplus_j Local_{ji},
       *   \f]
       *   where the direct sum is meant as a general binary operation, and
       *   \f$Local_{ji}\f$ is the i-th element held on the j-th MPI
       *   process
       *
       *  \param[in] localdata The local data that is being reduced
       *  \param[in] nelem     The length of the vector begin combined
       *  \param[in] target    Preallocated place for putting the result
       *  \param[in] op        The combining operation
       *  \param[in] Name      The name of the communicator we are combining
       *                       over
       */
      template<typename T>
      void all_reduce(const T* localdata, const int nelem,T* target,
           const MPIOperation& op,
            const std::string& Name)const{
         static_cast<const DerivedType*>(this)->AllReduceImpl(localdata,
               nelem,target,op, Name);
      }

      /** \brief Function for gathering a vector whose elements lie on
       *         different MPI processes.
       *
       *
       *   \param[in] localdata Chunk of data held locally
       *   \param[in] nelem     How long the local chunk is
       *   \param[in] target    A preallocated place for the gathered data
       *   \param[in] CommName  MPI communicator this gather will be performed
       *                        over (default=Current Comm)
       */
      template <class T>
      void all_gather(const T* localdata, const int nelem, T* target,
            const std::string& CommName="NONE") const {
         static_cast<const DerivedType*>(this)->all_gatherImpl(localdata, nelem,
               target, CommName);
      }

      ///Same as all_gather, except only root has copy
      template <class T>
      void gather(const T* localdata, const int nelem, T* target,const int Root,
            const std::string& CommName="NONE") const {
         static_cast<const DerivedType*>(this)->gatherImpl(localdata, nelem,
               target, Root,CommName);
      }

      /** \brief Function for broadcasting data held by one process to
       *         all other processes.
       *
       *
       *   \param[in] data A pre-allocated space for the data (if receiving)
       *                   else, the data to send
       *   \param[in] nelem The length of the data that is being sent
       *   \param[in] broadcaster Which MPI process is sending the data
       *   \param[in] CommName The name of the MPI communicator that the
       *                     broadcast is over (default=current)
       */
      template <class T>
      void bcast(T* data, const int nelem, const int broadcaster,
            const std::string& CommName="NONE") const {
         static_cast<const DerivedType*>(this)->bcastImpl(data, nelem,
               broadcaster, CommName);
      }

      ///See other bcast, difference is this is for serializable things
      ///like std::vector (untested!!!)
      template <class T>
      void bcast(T& data, const int broadcaster,
            const std::string& CommName="NONE") const {
         static_cast<const DerivedType*>(this)->bcastImpl(data,
               broadcaster, CommName);
      }

      /** \brief Sees if a message is available
       *
       *  \param[in] Sender The value of the MPI process who we want a message
       *                    from.  Negative int for any source.
       *  \param[in] MessageTag The label of the message.  Negative for any
       *                     messsage
       *
       *  \return The identity of the sender, or -1 if there is no message
       */
      virtual int Iprobe(const int Sender,const int MessageTag,
            const std::string& Comm="NONE") const {
          UNUSED(Sender);
          UNUSED(MessageTag);
          UNUSED(Comm);
         return true;
      }

      ///Same as Iprobe, except blocking
      virtual int probe(const int Sender,const int MessageTag,
            const std::string& Comm="NONE")const{
          UNUSED(Sender);
          UNUSED(MessageTag);
          UNUSED(Comm);
         return true;
      }


      /** \brief Non-blocking send
       *
       *   \param[in] destination The identity of the the message receiver,
       *                     negative values for any_source
       *   \param[in] tag    The message label to send, negative value
       *                     for any_tag
       *   \param[out] message A preallocated position for the message,NULL
       *                      for just a tag
       *   \param[in]  length The length of the message, 0 for just a tag
       *
       */
      template<typename T>
      void Isend(const int destination, const int tag, T* message=NULL,
            const int length=0,const std::string& Comm="NONE")const{
         static_cast<const DerivedType*>(this)->
               IsendImpl(destination,tag,message,length,Comm);
      }

      /** \brief blocking send
       *
       *   \param[in] destination The identity of the the message receiver,
       *                     negative values for any_source
       *   \param[in] tag    The message label to send, negative value
       *                     for any_tag
       *   \param[out] message A preallocated position for the message, NULL
       *                      for just a tag
       *   \param[in]  length The length of the message, 0 for just a tag
       *
       */
      template<typename T>
      void send(const int destination, const int tag, T* message=NULL,
            const int length=0,const std::string& Comm="NONE")const{
         static_cast<const DerivedType*>(this)->
               sendImpl(destination,tag,message,length,Comm);
      }

      /** \brief Non-blocking receive
       *
       *   \param[in] source The identity of the the message sender,
       *                     negative values for any_source
       *   \param[in] tag    The message label to receive, negative value
       *                     for any_tag
       *   \param[out] message A preallocated position for the message,NULL
       *                     if we are just sending a tag
       *   \param[in]  length The length of the message, 0 if we are just
       *                      sending a tag
       *
       */
      template<typename T>
      void Irecv(const int source, const int tag, T* message=NULL,
            const int length=0,const std::string& Comm="NONE")const{
         static_cast<const DerivedType*>(this)->
               IrecvImpl(source,tag,message,length,Comm);
      }

      /** \brief Blocking receive
       *
       *   \param[in] source The identity of the the message sender,
       *                     negative values for any_source
       *   \param[in] tag    The message label to receive, negative value
       *                     for any_tag
       *   \param[out] message A preallocated position for the message,
       *                     NULL for
       *   \param[in]  length The length of the message
       *
       */
      template<typename T>
      void recv(const int source, const int tag, T* message=NULL,
            const int length=0,const std::string& Comm="NONE")const{
         static_cast<const DerivedType*>(this)->
               recvImpl(source,tag,message,length,Comm);
      }

      ///Returns the current MPI process number
      virtual int me(const std::string& CommName="NONE") const {
          UNUSED(CommName);
         return 0;
      }

      ///Returns the current number of MPI processes
      virtual int nproc(const std::string& CommName="NONE") const {
          UNUSED(CommName);
         return 1;
      }

      ///Needed for legacy compatibility.  Don't use.
      virtual int nthread() const {
         return Process::environment.get_n_threads();
      }

      ///Returns the current communicator
      std::string communicator() const {
         return CurrentComm.back();
      }

      int thread_id(const pthread_t&) {
         return 0;
      }

      virtual void MakeComm(const std::string& /*Name*/, const int /*Color*/,
            const std::string& /*Comm2Split*/="NONE"){}

      virtual void FreeComm(const std::string& /*Name*/="NONE"){}
};// End Parallel base class
}//End namespace psi

#endif  /* _psi_src_lib_libparallel_parallel_h_ */
