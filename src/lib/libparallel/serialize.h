/*
 * File:   parallel.h
 * Author: jturney
 *
 * Created on May 28, 2010
 */

#ifndef _psi_src_lib_libparallel_serialize_h_
#define	_psi_src_lib_libparallel_serialize_h_

#include <boost/shared_ptr.hpp>
#include <psiconfig.h>
#include <cstdio>

namespace psi {

    class Communicator;

    class Serializable
    {
    public:
        virtual void send(Communicator* comm)  = 0;
        virtual void recv(Communicator* comm)  = 0;
        virtual void bcast(Communicator* comm, int broadcaster) = 0;
    };

}

#endif // _psi_src_lib_libparallel_serialize_h_

