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
#ifndef SIMPLESCHED_H_
#define SIMPLESCHED_H_

#include "../Communicator.h"
#include "MPIStaticScheduler.h"
namespace psi {
namespace LibParallel{
/** \brief An MPI scheduler designed to schedule a series of tasks that
 *         are all of approximately the same complexity.
 *
 *  This scheduler takes the number of tasks it has to schedule, N, and
 *  the number of MPI processes available to it, n, and assigns the tasks
 *  in consecutive batches of size, B=(N-R)/n, where R is the remainder,
 *  given by R=N%n.  Therefore, MPI process i is in charge of all tasks
 *  in the range: [i*B,i*B+B) where "[" is inclusive, and ")" is
 *  non-inclusive.  Additionally, MPI processes, i, such that i is between
 *  [0,R) is also responsible for one additional task, that with
 *  number (N-R)+i.
 *
 */
class SimpleScheduler:public StaticScheduler{
   private:
      ///Fills in the values of the member variables above and splits comm
      void SetValues();
   protected:
      ///The number of tasks we are scheduling
      int NTasks_;

      ///NTasks_ modulo the number of MPI processes, i.e. the remainder
      int Remainder_;

      ///NTasks_ less the remainder, i.e. how many tasks each process has
      int BatchSize_;

      ///For convenience where this process's tasks start
      int MyStart_;

      ///Same as MyStart_, but for the non-inclusive ending point
      int MyEnd_;

   public:

      ///Calls SetValues() to initialize the object
      SimpleScheduler(const int NTasks,
            boost::shared_ptr<const Communicator> IState);

      ///No Memory to free up
      virtual ~SimpleScheduler() {}

      template<typename T>
      std::vector<T> SynchImpl(const std::vector<T>& LocalValues,
            const int N)const;

      ///Actually reduces our data
      template <typename T>
      std::vector<T> ReduceImpl(const std::vector<T>&LocalValues,
            const int N,
            const MPIOperation& op) const;
};

template <typename T>
std::vector<T> SimpleScheduler::ReduceImpl(
      const std::vector<T>&LocalValues,
      const int N,
      const MPIOperation& op) const {
   std::vector<T> Result(N);
   IState_->AllReduce(&LocalValues[0], N, &Result[0], op);
   return Result;
}

template <typename T>
std::vector<T> SimpleScheduler::SynchImpl(
      const std::vector<T>& LocalValues,
      const int N)const {
   int Me=IState_->Me();
   std::vector<T> ReturnVec(NTasks_*N);
   IState_->AllGather(&LocalValues[0],N*BatchSize_,&ReturnVec[0]);
   //At this point we are good up to the remainder
   for (int i=0; i<Remainder_; i++) {
      int offset=(NTasks_-Remainder_)*N+i*N;
      if (Me==i) {
         //I know I'm broadcasting, and not receiving, so this
         //const_cast is safe
         T* Temp=const_cast<T*>(&LocalValues[BatchSize_*N]);
         IState_->Bcast(Temp, N, Me);
         size_t length=N*sizeof(T);
         memcpy(&ReturnVec[offset], &LocalValues[BatchSize_*N], length);
      }
      else IState_->Bcast(&ReturnVec[offset], N, i);
   }
   IState_->Barrier();
   return ReturnVec;
}
}}//End namespaces

#endif /* SIMPLESCHED_H_ */