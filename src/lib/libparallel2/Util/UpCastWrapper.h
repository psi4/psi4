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
#ifndef SRC_LIB_LIBPARALLEL2_UTIL_UPCASTWRAPPER_H_
#define SRC_LIB_LIBPARALLEL2_UTIL_UPCASTWRAPPER_H_
#include <vector>
#include <boost/shared_ptr.hpp>
#include "../Algorithms.h"
#include "../Schedulers/MPIScheduler.h"
#include "../Schedulers/DynamicRoundRobin.h"
#include "../Schedulers/StaticRoundRobin.h"
#include "../Schedulers/SimpleScheduler.h"
namespace psi{
namespace LibParallel{
class UpCastWrapper{
   private:
      boost::shared_ptr<MPIScheduler> Obj_;
      SchedAlgorithm Alg_;
   public:
      UpCastWrapper(boost::shared_ptr<MPIScheduler> object,
            const SchedAlgorithm& Alg):
         Obj_(object),Alg_(Alg){}

      template <typename T2>
      std::vector<T2> Reduce(
            const std::vector<T2>& Local,
            const int size,
            const MPIOperation& op) const;

      template <typename T2>
      std::vector<T2> Synch(
            const std::vector<T2>& Local,
            const int size) const;
};

template <typename T>
std::vector<T> UpCastWrapper::Synch(
      const std::vector<T>& Local,
      const int size) const{
   std::vector<T> Return;
   switch (Alg_) {
      case (SIMPLE): {
         Return=
               boost::dynamic_pointer_cast<SimpleScheduler>(Obj_)->
               SynchImpl(Local, size);
         break;
      }
      case (STATICRR):{
         Return=
               boost::dynamic_pointer_cast<StaticRoundRobin>(Obj_)->
               SynchImpl(Local, size);
         break;
      }
      case (DYNAMICRR):{
         Return=
               boost::dynamic_pointer_cast<DynamicRoundRobin>(Obj_)->
               SynchImpl(Local, size);
         break;
      }
   }
   return Return;
}



template <typename T>
std::vector<T> UpCastWrapper::Reduce(
      const std::vector<T>& Local,
      const int size,
      const MPIOperation& op) const{
   std::vector<T> Return;
   switch (Alg_) {
      case (SIMPLE): {
         Return=
               boost::dynamic_pointer_cast<SimpleScheduler>(Obj_)->
               ReduceImpl(Local, size, Op);
         break;
      }
      case (STATICRR):{
         Return=
               boost::dynamic_pointer_cast<StaticRoundRobin>(Obj_)->
               ReduceImpl(Local, size, Op);
         break;
      }
      case (DYNAMICRR):{
         Return=
               boost::dynamic_pointer_cast<DynamicRoundRobin>(Obj_)->
               ReduceImpl(Local, size, Op);
         break;
      }
   }
   return Return;
}


}}



#endif /* SRC_LIB_LIBPARALLEL2_UTIL_UPCASTWRAPPER_H_ */
