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

#include "mpi_wrapper.h"

#include <omp.h>
#if HAVE_MPI
namespace psi {
   void MPICommunicator::sync(const std::string& CommName) const {
      boost::mpi::communicator comm=GetComm(CommName);
      comm.barrier();
   }

   int MPICommunicator::me(const std::string& CommName) const {
      boost::mpi::communicator comm=GetComm(CommName);
      return comm.rank();
   }

   int MPICommunicator::nproc(const std::string& CommName) const {
      boost::mpi::communicator comm=GetComm(CommName);
      return comm.size();
   }

   boost::mpi::communicator MPICommunicator::GetComm(const std::string& CommName)
   const {
      std::string Comm2Check=(CommName=="NONE"?CurrentComm.back():CommName);
      if(Communicators.count(Comm2Check)!=1)
         throw PSIEXCEPTION("Comm: "+Comm2Check+" doesn't exist!!!");
      return (const_cast<MPICommunicator*>(this))->
      Communicators[Comm2Check];
   }

   MPICommunicator::MPICommunicator(int &argc,char **argv):
   Env(new boost::mpi::environment(argc,argv)) {
      boost::mpi::communicator world;
      Communicators["COMM_WORLD"]=world;
      //The next three lines are what the old local comm did
      //the code breaks if I do not include them
      //The way I understand this is that the number of openmp threads
      //is getting hard-coded to 1, and somewhere in the code people are
      //counting on this behavior...
      omp_set_nested(0);
      if (Process::environment("OMP_NUM_THREADS") == "")
      Process::environment.set_n_threads(1);
   }

   void MPICommunicator::MakeComm(const std::string& Name,const int Color,const std::string& Comm2Split) {
      boost::mpi::communicator comm=GetComm(Comm2Split);
      if(Communicators.count(Name)==1)
         throw PSIEXCEPTION("Communicator already exists.  Be more original with"
               " your name.  If you are using an MPIJob object, then boy are you"
               " unlucky (lucky?) you choose the same random number twice.");
      Communicators[Name]=comm.split(Color);
      CurrentComm.push_back(Name);
   }

   void MPICommunicator::FreeComm(const std::string& Name) {
      //Refuse to free mpi_comm_world
      if(CurrentComm.size()>1&&Name!="COMM_WORLD") {
         Communicators.erase(Name);
      }
   }

   int MPICommunicator::Iprobe(const int Sender,const int MessageTag,
               const std::string& Comm)const{
            int Sendcopy=(Sender<0?boost::mpi::any_source:Sender);
            int Messagecopy=(MessageTag<0?boost::mpi::any_tag:MessageTag);
            boost::optional<boost::mpi::status> status=
                  GetComm(Comm).iprobe(Sendcopy,Messagecopy);
            return (status?status.get().source():-1);
         }

   int MPICommunicator::probe(const int Sender,const int MessageTag,
               const std::string& Comm)const{
            int Sendcopy=(Sender<0?boost::mpi::any_source:Sender);
            int Messagecopy=(MessageTag<0?boost::mpi::any_tag:MessageTag);
            boost::optional<boost::mpi::status> status=
                  GetComm(Comm).probe(Sendcopy,Messagecopy);
            return (status?status.get().source():-1);
         }

}        //End namespace psi

#endif
