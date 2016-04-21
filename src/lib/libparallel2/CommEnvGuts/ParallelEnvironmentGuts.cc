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

#include "ParallelEnvironmentGuts.h"
#include "Communicator.h"
#if HAVE_MPI
   #include <boost/mpi/environment.hpp>
   #include <boost/mpi/communicator.hpp>
#endif
namespace psi{
namespace LibParallel{

ParallelEnvironmentGuts::ParallelEnvironmentGuts(int argc, char* argv[])
#if HAVE_MPI
   :Env_(new boost::mpi::environment(argc,argv))
#endif
{
   Comms_.push_back(
         boost::shared_ptr<Communicator>(
               new Communicator("COMM_WORLD",this)
   ));
}

boost::shared_ptr<const Communicator> ParallelEnvironmentGuts::GetComm()const{
   return Comms_.back();
}

void ParallelEnvironmentGuts::UpdateComms(){
   int MyTry=Comms_.size()-1;
   bool good=false;
   while(!good &&MyTry>0){
      if(!Comms_[MyTry]->Active()){
         Comms_.pop_back();
      }
      else good=true;
      MyTry--;
   }
}

void ParallelEnvironmentGuts::AddComm(boost::shared_ptr<Communicator> Comm){
  Comms_.push_back(Comm);
}

int ParallelEnvironmentGuts::Original()const{
   return Comms_[0]->Me();
}

void ParallelEnvironmentGuts::PrintOut()const{
   int size=Comms_.size();
   for(int i=0;i<size-1;i++){
      Comms_[i]->PrintOut();
   }
   (*this)<<"Current Communicator--> ";
   Comms_[size-1]->PrintOut();
}
}}