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
#include "MPISetterUpper.h"

#include "../Schedulers/DynamicRoundRobin.h"
#include "../Schedulers/MPIScheduler.h"
#include "../TaskJobGuts/MPITaskGuts.h"
#include "../Schedulers/SimpleScheduler.h"
#include "../Schedulers/StaticRoundRobin.h"
#include "TaskStatistics.h"
namespace psi {
namespace LibParallel {
typedef boost::shared_ptr<MPIScheduler> SharedSched;

bool Sorter(const MPITaskGuts& left,const MPITaskGuts& right){
   return left>right;
}

void MPISetterUpper::SortTasks() {
   std::sort(Tasks_->begin(), Tasks_->end(),Sorter);
}

void MPISetterUpper::ChooseAlgorithm(bool ForceDynamic) {
   TaskStats Stats(Tasks_,ForceDynamic,State_);
   this->SortTasks();
   ChoosenAlg_=Stats.SuggestedAlgorithm();
   switch (ChoosenAlg_) {
      case (SIMPLE): {
         Sched_=boost::shared_ptr<SimpleScheduler>(
               new SimpleScheduler(Tasks_->size(),State_));
         break;
      }
      case(STATICRR):{
         Sched_=boost::shared_ptr<StaticRoundRobin>(
               new StaticRoundRobin(Tasks_,State_));
         break;
      }
      case (DYNAMICRR): {
         Sched_=boost::shared_ptr<DynamicRoundRobin>(
          new DynamicRoundRobin(Tasks_,State_,2));
         break;
      }
      default: {
         throw PSIEXCEPTION("Scheduling algorithm not recognized");
         break;
      }
   }

}

SharedSched MPISetterUpper::GetScheduler()const { return Sched_;}

}} //End namespaces

