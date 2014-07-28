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

#ifndef _psi_src_lib_libparallel_local_h_
#define	_psi_src_lib_libparallel_local_h_

#include "parallel.h"
#include <cstring>
#include <string>
#include <omp.h>
namespace psi {

class LocalCommWrapper:public Parallel<LocalCommWrapper> {
   private:
      template <typename type>
      void bcastImpl(const type* data, const int nelem, const int broadcaster,
            const std::string& CommName="None") const {
      }
      template <typename T>
      void bcastImpl(T& data,const int broadcaster,
            const std::string&Comm="NONE") const {
      }

      template <typename type>
      void all_gatherImpl(const type* localdata, const int nelem, type* target,
            const std::string& CommName="NONE") const {
         if (localdata!=target)
            std::memcpy(const_cast<type*>(localdata), target,
                  sizeof(type)*nelem);
      }
      friend Parallel<LocalCommWrapper> ;
   public:
      LocalCommWrapper(const int &argc, char **argv) {
         //The next three lines are what the old local comm did
         //the code breaks if I do not include them
         //The way I understand this is that the number of openmp threads
         //is getting hard-coded to 1, and somewhere in the code people are
         //counting on this behavior...
         omp_set_nested(0);
         if (Process::environment("OMP_NUM_THREADS")=="")
            Process::environment.set_n_threads(1);
      }
};
// End of LocalCommWrapper class

}// End of namespace psi

#endif // End of _psi_src_lib_libparallel_local_h_
