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
#ifndef SRC_LIB_LIBPARALLEL2_LIBPARALLELHELPER_H_
#define SRC_LIB_LIBPARALLEL2_LIBPARALLELHELPER_H_
#include <boost/python.hpp>
#include "MPITask.h"

#include <vector>
namespace psi{
template <typename T> class MPIJob;
namespace LibParallel{


///Class designed to wrap the libparallel interface for python
class LibParallelHelper{
   private:
      std::vector<MPITask<boost::python::str> > Tasks_;
      MPIJob<boost::python::str>* Job_;
   public:
      ///Just sets Job_ to NULL
      LibParallelHelper();
      ///Frees Job_, iff it was allocated
      ~LibParallelHelper();

      void AddTask(boost::python::str& name, const int Prior);

      void MakeJob();

      boost::python::str Begin();

      bool Done();

      boost::python::str Next();

      boost::python::list Synch(boost::python::list& LocalValues,
            const int N=1);
};


}}//End namespaces

#endif /* SRC_LIB_LIBPARALLEL2_LIBPARALLELHELPER_H_ */