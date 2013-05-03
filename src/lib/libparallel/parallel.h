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
#include <psiconfig.h>
#include <cstdio>

// Define the default parallel environment. In general, developers should just use this type.
#if defined(HAVE_MPI)
#   include <libparallel/mpi_wrapper.h>
typedef psi::MPICommunicator            worldcomm;
#else
#   include <libparallel/local.h>
typedef psi::LocalCommWrapper           worldcomm;
#endif

#if 0
#if !defined(HAVE_MADNESS) && !defined(HAVE_ELEMENTAL)
#include <libparallel/local.h>
typedef     psi::LocalCommWrapper       worldcomm;
#elif !defined(HAVE_MADNESS) && defined(HAVE_ELEMENTAL)
#include <libparallel/elem.h>
typedef     psi::ElemCommWrapper        worldcomm;
#elif defined(HAVE_MADNESS)
#include <libparallel/mad.h>
typedef     psi::MADNESSCommWrapper     worldcomm;
#endif
#endif

namespace psi {

    extern FILE *outfile;
    extern boost::shared_ptr<worldcomm> WorldComm;

    // A templated version of init comunicator.
    template <typename comm_type>
    static comm_type* init_specific_communicator(int &argc, char **argv) {
        return new comm_type(argc, argv);
    }

    // Create a communicator from Comm typedef'ed above.
    worldcomm* Init_Communicator(int &argc, char **argv);

}

//#include "threaded_storage.h"

#endif  /* _psi_src_lib_libparallel_parallel_h_ */
