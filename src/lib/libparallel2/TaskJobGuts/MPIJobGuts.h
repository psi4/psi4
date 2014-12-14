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

#include "../Schedulers/MPIScheduler.h"
#include "../Util/MPISetterUpper.h"
#include "../MPITask.h"
#include "../Util/MPITaskQueue.h"
#include "../Communicator.h"
#include "../LibParallelBase.h"
#include "../Util/UpCastWrapper.h"
#include "../ParallelEnvironment.h"
namespace psi {
namespace LibParallel {
/** \brief Class meant to serve as the implementation of a MPI Job
 *
 */
class MPIJobGuts:public LibParallelBase {
   protected:
      ///The MPI state when we started
      boost::shared_ptr<const Communicator> IState_;

      /** \brief The MPI state after removing any processors we can't use
       *
       *  Within the advanced interface, the user is allowed to tell
       *  us that there are processes we can't use.  These are taken
       *  out of IState_, and the result is FState_.  For all intents
       *  and puropses all schedulers work on FState_.  If we are allowed
       *  to use all processors FState_ will be NULL
       */
      boost::shared_ptr<Communicator> FState_;

      ///The object that is assigning tasks
      boost::shared_ptr<MPIScheduler> Scheduler_;

      ///The current MPI process's job queue
      boost::shared_ptr<MPITaskQueue> Queue_;

      ///Did we split the comm
      bool Enough_;

      /** \brief True if this MPI process is actually participating in
       the Job
       *
       *  An MPI process may not be participating in a job if the
       *  advanced interface was used.  In this case Enough_ will
       *  be true, and Active_ will be true for all MPI processes
       *  that will participate, and false for the remainder.
       */
      bool Active_;

      /** \brief Determines the initial MPI state
       *
       *   If MaxProcs<0 then the initial MPI state is:
       *      - The current communicator
       *      - The current number of MPI processes
       *
       *   If it's 0, i'm just going to throw an error for now, but
       *   if someone can think of why this should be an allowed state
       *   they are welcome to code that up.
       *
       *   Otherwise, the initial state will be a new communicator, that we
       *   call NewComm, that is derived from the current communicator,
       *   that we call OldComm.  NewComm is comprised of processes
       *   0 through MaxProcs-1.
       */
      boost::shared_ptr<Communicator> InitialState(const int MaxProcs);

      /** \brief Tells you if there are enough processes to support Procs
       *
       *  If Procs<0 returns false.  If NProcs<=Procs returns false.
       *  If NProcs>Procs returns true
       */
      bool EnoughProcs(const int Procs) const;

   public:

      ///Returns whether enough processes were available upon creation
      bool EnoughProcs() const;

      template <typename T>
      MPIJobGuts(const std::vector<MPITask<T> >& Tasks, int MaxProcs=-1,
            bool ForceDynamic=false);

      ///Calls Queue_'s Next() if MPI process is active, else returns 0
      int Next();

      ///Calls Queue_'s Begin() if active, else returns 0
      int Begin();

      ///Calls Queue_'s Done(). If true, then calls Scheduler_'s Done()
      bool Done();

      ///If we made a comm on account of MaxProcs, we free it
      ~MPIJobGuts();

      template <typename T>
      std::vector<T> Synch(const std::vector<T>& LocalValues, const int N);

      template <typename T2>
      std::vector<T2> Reduce(const std::vector<T2>& LocalValues, const int N,
            const MPIOperation& op) const;

      ///Calls MPI barrier on FState_
      void Wait() const;

};

template <typename T>
MPIJobGuts::MPIJobGuts(const std::vector<MPITask<T> >& Tasks, int MaxProcs,
      bool ForceDynamic) :
      Active_(true), Enough_(false),
            IState_(Env_->GetComm()) {
   Enough_=MPIJobGuts::EnoughProcs(MaxProcs);
   if(Enough_)FState_=InitialState(MaxProcs);
   if (Active_) { ///Only set-up if the MPI process is participating
      MPISetterUpper SetUp(Tasks, (FState_?FState_:IState_),
                           ForceDynamic);
      Scheduler_=SetUp.GetScheduler();
      if (Scheduler_) Queue_=Scheduler_->GetTaskQueue();
   }
}

template <typename T>
std::vector<T> MPIJobGuts::Synch(const std::vector<T>& LocalValues,
      const int N) {
   std::vector<T> Result;
   if (Active_) {
      UpCastWrapper temp(Scheduler_,Scheduler_->Algorithm());
      Result=temp.Synch(LocalValues, N);
      if (IState_->Me()==0&&Enough_) {
         int size=Result.size();
         IState_->Bcast(&size, 1, 0);
         IState_->Bcast(&Result[0], size, 0);
      }
   }
   if (IState_->Me()!=0&&Enough_) {
      int size=0;
      IState_->Bcast(&size, 1, 0);
      std::vector<T> temp(size);
      IState_->Bcast(&temp[0], size, 0);
      Result=temp;
   }
   return Result;
}

template <typename T2>
std::vector<T2> MPIJobGuts::Reduce(const std::vector<T2>& LocalValues,
      const int N, const MPIOperation& op) const {
   std::vector<T2> Result;
   if (Active_) {
      UpCastWrapper temp(Scheduler_,Scheduler_->Algorithm());
      Result=temp.Reduce(LocalValues, N, op);
      if (IState_->Me()==0&&Enough_) {
         IState_->Bcast(&Result[0], N, 0);
      }
   }
   if (IState_->Me()!=0&&Enough_) {
      std::vector<T2> temp(N);
      IState_->Bcast(&temp[0], N, 0);
      Result=temp;
   }
   return Result;
}

}} //End namespaces

