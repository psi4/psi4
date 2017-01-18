/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

#ifndef _psi_src_lib_libparallel_local_h_
#define	_psi_src_lib_libparallel_local_h_

#include "parallel.h"
#include "libparallel.h"
#include <cstring>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

class LocalCommWrapper:public Parallel<LocalCommWrapper> {
   private:
      template <typename type>
      void bcastImpl(const type* /*data*/,
                     const int /*nelem*/,
                     const int /*broadcaster*/,
                     const std::string& /*CommName*/="None") const
      {}

      template <typename T>
      void bcastImpl(T& data,const int broadcaster,
            const std::string&Comm="NONE") const {
          UNUSED(data);
          UNUSED(broadcaster);
          UNUSED(Comm);
      }
      template<typename type>
      void AllReduceImpl(const type* localdata,const int nelem,
            type* target,
            const MPIOperation& operation,
            const std::string& Name){
          UNUSED(operation);
          UNUSED(Name);
         memcpy(target,localdata,nelem);
      }
      template<typename T>
      void IrecvImpl(const int source, const int tag, T* message,
            const int length)const{
          UNUSED(source);
          UNUSED(tag);
          UNUSED(message);
          UNUSED(length);
      }
      template<typename T>
      void recvImpl(const int source, const int tag, T* message,
            const int length)const{
          UNUSED(source);
          UNUSED(tag);
          UNUSED(message);
          UNUSED(length);
      }
      template<typename T>
      void IsendImpl(const int source, const int tag, T* message,
            const int length)const{
          UNUSED(source);
          UNUSED(tag);
          UNUSED(message);
          UNUSED(length);
      }
      template<typename T>
      void sendImpl(const int source, const int tag, T* message,
            const int length)const{
          UNUSED(source);
          UNUSED(tag);
          UNUSED(message);
          UNUSED(length);
      }
      template <typename type>
      void all_gatherImpl(const type* localdata,
                          const int nelem,
                          type* target,
                          const std::string& /*CommName*/="NONE") const
      {
         if (localdata!=target)
            std::memcpy(const_cast<type*>(localdata), target,
                  sizeof(type)*nelem);
      }
      template <typename type>
      void gatherImpl(const type* localdata, const int nelem, type* target,
            const int Root,const std::string& CommName="NONE") const {
          UNUSED(Root);
          UNUSED(CommName);
         if (localdata!=target)
            std::memcpy(const_cast<type*>(localdata), target,
                  sizeof(type)*nelem);
      }
      friend class Parallel<LocalCommWrapper> ;
   public:
      LocalCommWrapper(const int &/*argc*/, char **/*argv*/) {
         //The next three lines are what the old local comm did
         //the code breaks if I do not include them
         //The way I understand this is that the number of openmp threads
         //is getting hard-coded to 1, and somewhere in the code people are
         //counting on this behavior...
#ifdef _OPENMP
         omp_set_nested(0);
#endif
         if (Process::environment("OMP_NUM_THREADS")=="")
            Process::environment.set_n_threads(1);
      }
};
// End of LocalCommWrapper class

}// End of namespace psi

#endif // End of _psi_src_lib_libparallel_local_h_