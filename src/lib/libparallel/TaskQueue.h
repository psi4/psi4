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
#ifndef TASKQUEUE_H_
#define TASKQUEUE_H_

#include "psi4-dec.h"
#include <time.h>
#include <stdlib.h>

namespace psi{

class MPITask{
   private:
      int TaskNumber_;
   public:
      MPITask(const int TaskNumber):TaskNumber_(TaskNumber){}
      virtual ~MPITask(){}
};

class MPIQueue{
   protected:
      typedef std::vector<boost::shared_ptr<MPITask> > QueueType;
      QueueType Queue_;
   public:
      MPIQueue(const QueueType& Queue):Queue_(Queue){}
      typedef std::vector<boost::shared_ptr<MPITask> >::iterator ItrType;
      ItrType Begin(){return Queue_.begin();}
      ItrType End(){return Queue_.end();}

};

/** \brief Base for an MPI Task scheduler
 *
 *  Each MPI task, on a communicator, calls the Schedule function and is
 *  then given an MPIQueue of tasks to run.
 */
template<typename T>
class MPIScheduler{
   private:
      T* Base_;
   protected:
      ///Comm name, uniqueness based on random number
      std::string Comm_;
      ///The number of processors Comm_ has
      int NProcs_;
      ///My index
      int Me_;
      MPIScheduler(T* Base):
         Base_(Base),NProcs_(WorldComm->nproc()),Me_(WorldComm->me()){

         srand(time(NULL));
         std::stringstream Name;
         Name<<"Comm"<<rand()<<std::endl;
      }
   public:
      virtual boost::shared_ptr<MPIQueue> Schedule(int Tasks2Schedule)=0;
      template<typename T1>
      void Synch(std::vector<MPIProperty<T1> >& prop){
         Base_->SynchImpl(prop);
      }
      virtual ~MPIScheduler(){}
};
}//End namespace



#endif /* TASKQUEUE_H_ */
