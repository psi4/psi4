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
#ifndef SLAVETOOLS_H_
#define SLAVETOOLS_H_

#include "../Communicator.h"
#include "../Algorithms.h"
#include "ToolSet.h"
namespace psi{
namespace LibParallel{
class MPITaskQueue;
class SlaveTools:public ToolSet{
   protected:
      boost::shared_ptr<MPITaskQueue> MyTasks_;
      bool Released_;
   public:
      ///Used for debug printing
      void PrintOut()const;
      SlaveTools(boost::shared_ptr<const Communicator> State);
      boost::shared_ptr<MPITaskQueue> Tasks();
      void GiveTasks();
      bool Done();

      template<typename T>
      std::vector<T> SynchImpl(
            const std::vector<T>& Local,
            const int N,
            boost::shared_ptr<const Communicator> Comm3_);

      template<typename T>
      std::vector<T> ReduceImpl(
            const std::vector<T>& Local,
            const int N,
            const MPIOperation& op,
            boost::shared_ptr<const Communicator> Comm3_);
};

template<typename T>
std::vector<T> SlaveTools::ReduceImpl(
      const std::vector<T>& Local,
      const int N,
      const MPIOperation& op,boost::shared_ptr<const Communicator> Comm3_){
   std::vector<T> Results(N);
   boost::shared_ptr<Communicator> Comm=Comm3_->MakeComm(false);
   Comm->AllReduce(&Local[0],N,&Results[0],op);
   if(Comm3_->Me()==1)Comm3_->Send(0,1,&Results[0],N);
   Comm->FreeComm();
   return Results;
}

template<typename T>
std::vector<T> SlaveTools::SynchImpl(const std::vector<T>& Local, const int,
      boost::shared_ptr<const Communicator> Comm3_){
   int size=Local.size();
   Comm3_->Send(0,1,&size,1);
   Comm3_->Send(0,1,&Local[0],size);
   Comm3_->Bcast(&size,1,0);
   std::vector<T> Results(size);
   Comm3_->Bcast(&Results[0],size,0);
   return Results;

}

}}




#endif /* SLAVETOOLS_H_ */