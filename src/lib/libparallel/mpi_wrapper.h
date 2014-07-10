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
#ifndef _psi_src_lib_libparallel_mpi_wrapper_h_
#define _psi_src_lib_libparallel_mpi_wrapper_h_

#include <string>

#include <boost/shared_ptr.hpp>
#include "parallel.h"
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
      void sync(const std::string& CommName="NONE") const{
         boost::mpi::communicator comm=GetComm(CommName);
         comm.barrier();
      }

      int me(const std::string& CommName="NONE") const{
         boost::mpi::communicator comm=GetComm(CommName);
         return comm.rank();
      }

      int nproc(const std::string& CommName="NONE") const{
         boost::mpi::communicator comm=GetComm(CommName);
         return comm.size();
      }

      void MakeComm(const std::string& Name, const int Color,
            const std::string& Comm2Split="NONE");

      ///Frees the processes in the communicator
      void FreeComm(const std::string& Name="NONE");

      ///Returns the MPI_Comm associated with the current communicator
      MPI_Comm GetMPIComm()const{return (MPI_Comm)GetComm("NONE");}

   private:
      /**\brief This is Boost's MPI environment object.
       *
       * In its creation Boost calls MPI_Initialize and in its destruction
       *  MPI_Finalize is called.  Consequentially, it is imperative that
       *  only one of these objects exist, and it is the first thing
       *  created by Psi and the last thing destroyed by Psi.
       *
       */
      boost::shared_ptr<boost::mpi::environment> Env;

      ///This is a mapping between communicator names and communicators.
      std::map<std::string, boost::mpi::communicator> Communicators;

      /** \brief Helper function that returns the appropriate communicator by name.
       *
       * If Comm is "", then the default communicator is used, otherwise whatever
       * is passed in as Comm is used as the key value.
       *
       * \param[in] Comm The name of the communicator we want
       */
      boost::mpi::communicator GetComm(const std::string& CommName) const {
         std::map<std::string, boost::mpi::communicator>DaCopy(Communicators);
         return DaCopy[(CommName!="NONE" ? CommName : CurrentComm.back())];
      }

      ///So that the base class can access the implementation
      friend Parallel<MPICommunicator> ;

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
      void bcastImpl(T* data, int nelem,const int broadcaster,
            const std::string& Comm="NONE") const {
         boost::mpi::broadcast(GetComm(Comm), data, nelem, broadcaster);
      }

      template <typename T>
      void bcastImpl(T& data,const int broadcaster,
            const std::string&Comm="NONE")const{
         boost::mpi::broadcast(GetComm(Comm), data, broadcaster);
      }

};
}//End namespace psi
//End mpiwrapper class
#endif //End on HAVE_MPI


#endif //End on header guard
