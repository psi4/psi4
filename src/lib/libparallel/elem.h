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

#ifndef _psi_src_lib_libparallel_elem_h_
#define	_psi_src_lib_libparallel_elem_h_

#include <psi4-dec.h>
#include <cstring>
#include <libparallel/openmp.h>

#if defined(HAVE_ELEMENTAL)

#ifdef _OPENMP
#include <omp.h>
#endif

#include "elemental.hpp"

namespace psi {

class ElemCommWrapper {
private:
    int me_;
    int nproc_;
    int nthread_;
    std::string communicator_;
    static const int DEFAULT_SEND_RECV_TAG = 1000;
    elem::mpi::Comm comm_;
    MPI_Request req_;
    MPI_Status stat_;
    typedef unsigned char byte;

    void bcast_vector_string(std::vector<std::string> &data, int broadcaster=0)
    {
        int vsize = data.size();
        bcast(&vsize, 1, broadcaster);

        if (me_ != broadcaster) data.resize(vsize);

        for (int i=0; i < data.size(); i++) {
            int ssize = data[i].size();
            bcast(&ssize, 1, broadcaster);

            if (me_ != broadcaster) data[i].resize(ssize);
            bcast((char *) &(data[i][0]), data[i].size(), broadcaster);
        }
    }

public:

    ElemCommWrapper (const int &argc, char **argv)
    {
        init_openmp();
        elem::Initialize(argc, argv);
        comm_ = elem::mpi::COMM_WORLD;
        me_ = elem::mpi::CommRank(comm_);
        nproc_ = elem::mpi::CommSize(comm_);
        nthread_ = Process::environment.get_n_threads();
        communicator_ = "ELEMENTAL";
    }

    ElemCommWrapper(const ElemCommWrapper  &copy) { *this = copy; }

    ~ElemCommWrapper ()
    { }

    ElemCommWrapper & operator=(const ElemCommWrapper & other)
    {
        me_   = other.me_;
        nproc_ = other.nproc_;
        nthread_ = other.nthread_;
        return *this;
    }

    inline int thread_id(const pthread_t &thread) { return 0; }

    inline void sync() { elem::mpi::Barrier(comm_); }

    inline void print(std::string OutFileRMR) const
    {
        if (me_ == 0) {
            printer->Printf( "\n    Using ElemCommWrapper  (Number of procs = %d)\n", nproc_);
            printer->Printf( "                          (Number of threads in pool = %d)\n\n", nthread_);
        }
    }

    inline void finalize()
    {
        elem::Finalize( );
    }

    template<typename type>
    inline void send(type *data, int nelem, int target, int tag=DEFAULT_SEND_RECV_TAG) const
    {
        elem::mpi::ISend(data, nelem, target, tag, comm_, req_);
    }

    template<typename type>
    inline void recv(type *data, int nelem, int sender, int tag=DEFAULT_SEND_RECV_TAG) const
    {
        elem::mpi::IRecv(data, nelem, sender, tag, comm_, req_);
    }

    template<typename type>
    inline void bcast(type* data, size_t nelem, int broadcaster=0)
    {
        elem::mpi::Broadcast((byte *)data, nelem*sizeof(type), broadcaster, comm_);
    }

    template<typename type>
    inline void bcast_serializable(type &data, int broadcaster=0)
    {
        // this is a temporary hack.  Need to serialize vector, strings, etc.
        bcast_vector_string(data, broadcaster);
    }

    template<typename type>
    inline void sum(type *data, int nelem, type *receive_buffer=0, int target=-1)
    {
        if (target == -1) {
            if (receive_buffer == 0) {
                receive_buffer = new type[nelem];
                elem::mpi::AllReduce(&data[0], receive_buffer, nelem, elem::mpi::SUM, comm_);
                std::memcpy(data, receive_buffer, nelem*sizeof(type));
                delete receive_buffer;
            }
            else
                elem::mpi::AllReduce(&data[0], receive_buffer, nelem, elem::mpi::SUM, comm_);
        }
        else
            elem::mpi::Reduce(&data[0], receive_buffer, nelem, elem::mpi::SUM, target, comm_);
    }

    inline int me() const { return me_; }
    inline int nproc() const { return nproc_; }
    inline int nthread() const { return nthread_; }
    inline std::string communicator() const {return communicator_;}

}; // End of ElemCommWrapper class

} // End of psi namespace

#endif // End of HAVE_ELEMENTAL

#endif // End of _psi_src_lib_libparallel_elem_h_
