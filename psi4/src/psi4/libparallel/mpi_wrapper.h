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
#ifndef _psi_src_lib_libparallel_mpi_wrapper_h_
#define _psi_src_lib_libparallel_mpi_wrapper_h_

#include <string>
#include <functional>
 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include "parallel.h"
#include "libparallel.h"
#include <map>

#if HAVE_MPI
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <mpi.h>
namespace psi {
   class MPICommunicator:public Parallel<MPICommunicator> {
      public:
      //Convenient typedef of base type
      typedef Parallel<MPICommunicator> Base;

      ///Starts an MPI session, really should only be used in the beginning of Psi's execution
      MPICommunicator(int &argc, char **argv);

      ///Memory management is handled via smart pointer, so does nothing
      ~MPICommunicator() {}

      /** \brief Basic assignment operator, checks for self-assignment.
       *
       * As is standard practice, this MPICommunicator calls the base assingment
       * operator before copying MPICommunicator member variables.
       *
       * \param[in] other The MPICommunicator we are copying
       */
      MPICommunicator& operator=(const MPICommunicator& other);

      /** \brief Provides MPI_Barrier functionality
       *
       * The point of MPI_Barrier is to halt further execution of a program until
       * all processes in the communicator reach that point.  In that sense it
       * synchronizes the processes.  Note in particular that this is specific
       * to the communicator that it is called on.
       *
       * \param[in] Comm The communicator we are erecting a barrier on.  Defaults to current Comm.
       */
      void sync(const std::string& CommName="NONE") const;

      int me(const std::string& CommName="NONE") const;

      int nproc(const std::string& CommName="NONE") const;

      void MakeComm(const std::string& Name, const int Color,
            const std::string& Comm2Split="NONE");

      int Iprobe(const int Sender,const int MessageTag,
            const std::string& Comm)const;

      int probe(const int Sender,const int MessageTag,
            const std::string& Comm)const;



      ///Frees the processes in the communicator
      void FreeComm(const std::string& Name="NONE");

      ///Returns the MPI_Comm associated with the current communicator
      MPI_Comm GetMPIComm()const {return (MPI_Comm)GetComm("NONE");}

      private:
      /**\brief This is Boost's MPI environment object.
       *
       * In its creation Boost calls MPI_Initialize and in its destruction
       *  MPI_Finalize is called.  Consequentially, it is imperative that
       *  only one of these objects exist, and it is the first thing
       *  created by Psi and the last thing destroyed by Psi.
       *
       */
      std::shared_ptr<boost::mpi::environment> Env;

      ///This is a mapping between communicator names and communicators.
      std::map<std::string, boost::mpi::communicator> Communicators;

      /** \brief Helper function that returns the appropriate communicator by name.
       *
       * If Comm is "", then the default communicator is used, otherwise whatever
       * is passed in as Comm is used as the key value.
       *
       * \param[in] Comm The name of the communicator we want
       */
      boost::mpi::communicator GetComm(const std::string& CommName) const;

      ///So that the base class can access the implementation
      friend class Parallel<MPICommunicator>;

      /* \brief Wrapper to MPI's Allgatherv
       *
       * Allgatherv parallelizes the process of calculating the elements of some vector V by P
       * processes.  Unlike Allgather, Allgatherv allows each process to do different sized chunks
       * of V.  At the end of the computation Allgatherv makes sure each process has its own copy
       * of V.
       *
       * Let's assume we have P processes, that are splitting up the task of computing the
       * elements of a vector V.  The i'th process, is incharge of n_i data points in V.
       *
       * \param[in]  data   This is the data calculated on i
       * \param[in]  nbyte  This is the number of bytes that is to be transferred
       * \param[out] target This is V
       * \param[in]  counts This is a P element long array, where P_i=n_i
       * \param[in]  displs This is the offset in V, where i's elements go
       *
       *
       */
      template <typename T>
      void all_gatherImpl(const T* data, int nelem, T* target,
            const std::string& Comm="NONE") const {
         boost::mpi::all_gather(GetComm(Comm), data, nelem, target);
      }

      template <typename T>
      void gatherImpl(const T* data, int nelem, T* target, const int Root,
            const std::string& Comm="NONE") const {
         boost::mpi::gather(GetComm(Comm), data, nelem, target,Root);
      }

      template <typename T>
      void bcastImpl(T* data, int nelem,const int broadcaster,
            const std::string& Comm="NONE") const {
         boost::mpi::broadcast(GetComm(Comm), data, nelem, broadcaster);
      }

      template <typename T>
      void bcastImpl(T& data,const int broadcaster,
            const std::string&Comm="NONE")const {
         boost::mpi::broadcast(GetComm(Comm), data, broadcaster);
      }

      template<typename T>
      void IrecvImpl(const int source, const int tag, T* message,
            const int length,const std::string& Comm)const{
         int SCopy=(source<0?boost::mpi::any_source:source);
         int TCopy=(tag<0?boost::mpi::any_tag:tag);
         int LCopy=length;
         if(message!=NULL&&length!=0){
            if(length!=1)GetComm(Comm).irecv(SCopy,TCopy,message,LCopy);
            else GetComm(Comm).irecv(SCopy,TCopy,message[0]);
         }
         else
            GetComm(Comm).irecv(SCopy,TCopy);
      }

      template<typename T>
      void recvImpl(const int source, const int tag, T* message,
            const int length,const std::string& Comm)const{
         int SCopy=(source<0?boost::mpi::any_source:source);
         int TCopy=(tag<0?boost::mpi::any_tag:tag);
         int LCopy=length;
         if(message!=NULL&&length!=0){
            if(length!=1)GetComm(Comm).recv(SCopy,TCopy,message,LCopy);
            else GetComm(Comm).recv(SCopy,TCopy,message[0]);
         }
         else
            GetComm(Comm).recv(SCopy,TCopy);
      }

      template<typename T>
      void IsendImpl(const int source, const int tag,const T* message,
            const int length,const std::string& Comm)const{
         int SCopy=(source<0?boost::mpi::any_source:source);
         int TCopy=(tag<0?boost::mpi::any_tag:tag);
         int LCopy=length;
         if(message!=NULL&&length!=0){
            if(length!=1)GetComm(Comm).isend(SCopy,TCopy,message,LCopy);
            else GetComm(Comm).isend(SCopy,TCopy,message[0]);
         }
         else
            GetComm(Comm).isend(SCopy,TCopy);
      }

      template<typename T>
      void sendImpl(const int source, const int tag,const T* message,
            const int length,const std::string& Comm)const{
         int SCopy=(source<0?boost::mpi::any_source:source);
         int TCopy=(tag<0?boost::mpi::any_tag:tag);
         int LCopy=length;
         if(message!=NULL&&length!=0){
            if(length!=1)GetComm(Comm).send(SCopy,TCopy,message,LCopy);
            else GetComm(Comm).send(SCopy,TCopy,message[0]);
         }
         else GetComm(Comm).send(SCopy,TCopy);
      }

      template<typename T>
      void AllReduceImpl(const T* localdata,const int nelem,T* target,
            const MPIOperation& op,
            const std::string& Name)const {
         switch(op) {
            case(MULTIPLY): {
               boost::mpi::all_reduce(GetComm(Name),localdata,
                     nelem, target,std::multiplies<T>());
               break;
            }
            case(DIVIDE): {
               boost::mpi::all_reduce(GetComm(Name),localdata,
                     nelem, target,std::divides<T>());
               break;
            }
            case(ADD): {
               boost::mpi::all_reduce(GetComm(Name),localdata,
                     nelem, target,std::plus<T>());
               break;
            }
            case(SUBTRACT): {
               boost::mpi::all_reduce(GetComm(Name),localdata,
                     nelem, target,std::minus<T>());
               break;
            }
            case(MODULUS): {
               boost::mpi::all_reduce(GetComm(Name),localdata,
                     nelem, target,std::modulus<T>());
               break;
            }
            default:
            throw PSIEXCEPTION("Unrecognized operation in mpiwrapper.h");
            break;
         }
      }
   };
} //End namespace psi
//End mpiwrapper class
#endif //End on HAVE_MPI

#endif //End on header guard
