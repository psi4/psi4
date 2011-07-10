/*
 * File:   parallel.h
 * Author: jturney
 *
 * Created on May 28, 2010
 */

#ifndef _psi_src_lib_libparallel_serialize_h_
#define	_psi_src_lib_libparallel_serialize_h_

namespace psi {

    class Communicator;

    class Serializable
    {
    public:
        virtual void send()  = 0;
        virtual void recv()  = 0;
        virtual void bcast(int broadcaster) = 0;
        virtual void sum() = 0;
    };

}

#endif // _psi_src_lib_libparallel_serialize_h_
