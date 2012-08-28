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
#if !defined(HAVE_MADNESS) && !defined(HAVE_ELEMENTAL)
#include <libparallel/local.h>
typedef     psi::LocalCommWrapper       worldcomm;
#elif !defined(HAVE_MADNESS) && defined(HAVE_ELEMENTAL)
#include <libparallel/elem.h>
typedef     psi::ElemCommWrapper             worldcomm;
#elif defined(HAVE_MADNESS)
#include <libparallel/mad.h>
typedef     psi::MADNESSCommWrapper     worldcomm;
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

#endif  /* _psi_src_lib_libparallel_parallel_h_ */
