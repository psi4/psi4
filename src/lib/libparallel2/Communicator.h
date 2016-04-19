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
#ifndef SRC_LIB_LIBPARALLEL2_COMMUNICATOR_H_
#define SRC_LIB_LIBPARALLEL2_COMMUNICATOR_H_
#include "CommEnvGuts/CommGuts.h"
#include "Algorithms.h"
#include <boost/shared_ptr.hpp>

namespace psi{
namespace LibParallel{
class ParallelEnvironment;

/** \brief The public interface to my communicator objects
 *
 *  All of the local vs. MPI stuff is buried within this wrapper.
 *  This object works regardless of whether MPI is running or not.
 */
class Communicator: private CommGuts{
   public:
      ///Initializes this communicator as a copy of other
      Communicator(const Communicator& other);

      ///Makes this comm a copy of other by assignment
      const Communicator& operator=(const Communicator& other);

      ///Creates a comm w/ name "Name", part of environment "Env"
      Communicator(const std::string& Name,ParallelEnvironmentGuts* Env);

      ///True if MPI processes of current communicator are still in use
      bool Active()const;

      ///Frees the MPI processes of the current communicator
      void FreeComm();

      /** \brief Splits the current communicator based on Color
       *
       *   \param[in] Color Which subcommunicator I'm going to.
       *   \param[in] Name The name of the resulting communicator. If set
       *                   to "" a random name will be generated.
       */
      boost::shared_ptr<Communicator> MakeComm(const int Color,
            const std::string& Name="")const;

      ///Halts program progress till all processes on current comm arrive
      void Barrier()const;

      /** \brief Function to see if a message is available
       *
       *  In Master/Slave models the master needs to see if messages
       *  are available from the slaves.  This is what the probe routine
       *  is for.
       *
       *  \param[in] Sender What MPI process are we checking for a message
       *                    (default=-1, which is any process)
       *  \param[in] MessageTag What is the tag of the message we are
       *                        looking for? (default=-1, any tag)
       *  \param[in] Block Is this a blocking operation (default true)
       *  \return The number of the MPI process that has a message.  A
       *          negative number means no message is available
       *
       */
      int Probe(const int Sender=-1,const int MessageTag=-1,bool Block=true)const;

      ///Returns the identity of the current process
      int Me()const;

      ///Returns the number of processes on the current comm
      int NProc()const;

      /** \brief Performs an all gather operation on the current comm
       *
       *  If the current communicator has N processes, and each process
       *  is holding onto M data points, this will gather the M data
       *  points from each process into a vector of length N*M, where
       *  elements (i-1)*M, through i*M-1 are the M elements that came
       *  from process i (i runs from 0 to N-1).  Note each process
       *  must send exactly M elements.  Additional note, Target is
       *  expected to be preallocated.
       *
       *  \param[in] LocalData The M data points this process has collected
       *  \param[in] NElem   M, from the description above
       *  \param[out] Target  The N*M long array that is returned (must be
       *                     allocated prior to this call).
       */
      template <typename T>
      void AllGather(const T* LocalData, const int NElem,
            T* Target)const{
         AllGatherImpl(LocalData,NElem,Target);
      }

      /** \brief Performs an all reduce operation on the current comm
       *
       *  If we have N communicators, that generated M data points each,
       *  each, process i has an array, \f$A_i\f$, that looks like:
       *  \f$A_i=[a_0,a_1,\ldots,a_M]\f$, where the various \f$a_j\f$
       *  elements are the M data points generated on that process.
       *  What this function does is creates an array B, on each
       *  process, whose j-th element is given by:
       *  \f[
       *    B_j=\bigotimes_i^N A_{ij},
       *  \f]
       *  where \f$\bigotimes\f$ is a generic operation, such as
       *  addition, subtraction, multiplication, etc, and \f$A_{ij}\f$
       *  is the j-th element of \f$A_i\f$.  Note that  A complete
       *  list of available operations is available in Algorithms.h,
       *  and this function assumes B is preallocated.
       *
       *  \param[in] LocalData A_i from the description
       *  \param[in] NElem  M from the description
       *  \param[out] Target The array B from the description,
       *                     preallocated
       *  \param[in] Op The operation used for the reduction, see
       *                Algorithms.h for a complete list.
       */
      template<typename T>
      void AllReduce(const T* LocalData, const int NElem,T* Target,
                     const MPIOperation& Op)const{
         AllReduceImpl(LocalData,NElem,Target,Op);
      }

      /** \brief Analogous to AllGather, except the result is not
       *         given to each process, instead it only exists
       *         on a root process.
       *
       *   See AllGather for a complete description.
       *
       *   \param[in] LocalData The data the current process has
       *                  collected
       *   \param[in] NElem    The number of data points the
       *                        in LocalData
       *   \param[out] Target  Where the data will go on the root
       *                       process (ignored on all others)
       *   \param[in]  Root    Which process is the root process
       */
      template<typename T>
      void Gather(const T* LocalData, const int NElem,T* Target,
            const int Root)const{
         GatherImpl(LocalData,NElem,Target,Root);
      }

      template<typename T>
      void Bcast(T* Data, const int NElem, const int Broadcaster)const{
         BcastImpl(Data,NElem,Broadcaster);
      }

      template<typename T>
      void Send(const int Receiver,const int MessageTag,T* Message=NULL,
                  const int Length=0,bool Block=true)const{
         SendImpl(Receiver,MessageTag,Message,Length,Block);
      }

      template<typename T>
      void Receive(const int Sender, const int MessageTag,T* Message=NULL,
            const int Length=0,bool Block=true)const{
         ReceiveImpl(Sender,MessageTag,Message,Length,Block);
      }

      ///Useful debugging information
      void PrintOut()const;

};

}}



#endif /* SRC_LIB_LIBPARALLEL2_COMMUNICATOR_H_ */