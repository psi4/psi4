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

/*
 * File:   parallel.h
 * Author: jturney, jjwilke
 *
 * Created on December 11, 2009, 3:34 PM
 */

#ifndef _psi_src_lib_libparallel_parallel_h_
#define	_psi_src_lib_libparallel_parallel_h_

#include <boost/shared_ptr.hpp>
#include <vector>
#include <string>
#include "process.h"

namespace psi {

template <typename DerivedType>
class Parallel {
   public:
      ///This typedef is the type of this class
      typedef Parallel<DerivedType> ThisType;
   protected:
      ///Array of our communicator names, in the order they are derived
      std::vector<std::string> CurrentComm;
   public:
      Parallel() {
         CurrentComm.push_back("COMM_WORLD");
      }
      virtual ~Parallel() {
      }
      virtual void sync(const std::string& /*CommName*/="NONE") const {
      }

      template <class T>
      void all_gather(const T* localdata, const int nelem, T* target,
            const std::string& CommName="NONE") const {
         static_cast<const DerivedType*>(this)->all_gatherImpl(localdata, nelem,
               target, CommName);
      }

      template <class T>
      void bcast(T* data, const int nelem, const int broadcaster,
            const std::string& CommName="NONE") const {
         static_cast<const DerivedType*>(this)->bcastImpl(data, nelem,
               broadcaster, CommName);
      }

      template <class T>
      void bcast(T& data, const int broadcaster,
            const std::string& CommName="NONE") const {
         static_cast<const DerivedType*>(this)->bcastImpl(data,
               broadcaster, CommName);
      }

      virtual int me(const std::string& /*CommName*/="NONE") const {
         return 0;
      }

      virtual int nproc(const std::string& /*CommName*/="NONE") const {
         return 1;
      }

      virtual int nthread() const {
         return Process::environment.get_n_threads();
      }

      std::string communicator() const {
         return CurrentComm.back();
      }

      int thread_id(const pthread_t&) {
         return 0;
      }

      virtual void MakeComm(const std::string& /*Name*/, const int /*Color*/,
            const std::string& /*Comm2Split*/="NONE"){}

      virtual void FreeComm(const std::string& /*Name*/="NONE"){}
};// End Parallel base class
}//End namespace psi
#endif  /* _psi_src_lib_libparallel_parallel_h_ */
