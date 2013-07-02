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

#ifndef _psi_src_lib_libparallel_mad_h_
#define	_psi_src_lib_libparallel_mad_h_

#include <psi4-dec.h>
#include <cstring>
#include <libparallel/openmp.h>

#if defined(HAVE_MADNESS)

#if defined(HAVE_ELEMENTAL)
#include "elemental.hpp"
#endif

#ifdef WORLD_INSTANTIATE_STATIC_TEMPLATES
#undef WORLD_INSTANTIATE_STATIC_TEMPLATES
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#else
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#endif

#include <world/world.h>
#include <world/worldobj.h>
#include <world/worlddc.h>
#include <world/archive.h>
#include <world/safempi.h>
#include <tensor/tensor.h>

namespace psi {

typedef boost::shared_ptr<madness::Spinlock> SharedLock;
typedef boost::shared_ptr<madness::MutexFair> SharedMutex;
typedef boost::shared_ptr<madness::World> SharedMadWorld;

class MADNESSCommWrapper {
private:
    int me_;
    int nproc_;
    int nthread_;
    std::string communicator_;
    SharedMadWorld madworld_;
    std::map<pthread_t, int> thread_id_;
    static const int DEFAULT_SEND_RECV_TAG = 1000;

public:

    MADNESSCommWrapper(const int &argc, char **argv)
    {
        init_openmp();

        // Madness is initialized in the communicator constructor
        // Initialize madness and create a madness world communicator
        madness::initialize(argc, argv);

        madworld_ = SharedMadWorld(new madness::World(MPI::COMM_WORLD));

        me_ = madworld_->rank();
        nproc_ = madworld_->nproc();
        nthread_ = madworld_->taskq.nthread();
        Process::environment.set_n_threads(nthread_);
        communicator_ = "MADNESS";


#if defined(HAVE_ELEMENTAL)
        elem::Initialize( argc, argv );
#endif

        thread_id_.clear();

        // If we have threads in the thread pool, this gets all of the pthread ids from madness, and
        // creates a map of pthread ids to a unique integer value (i.e. 0, 1, 2, ...)
        int i = 0;
        if (Process::environment("MAD_NUM_THREADS") != "1" && Process::environment("POOL_NTHREAD") != "0") {
            pthread_t *id = madworld_->taskq.thread_id();
            for (i=0; i < nthread_; i++) {
                thread_id_.insert( std::pair<pthread_t, int>(id[i], i) );
            }
            free(id);
        }
        // We need to add the master thread because it can perform tasks that are in the task queue
        thread_id_.insert( std::pair<pthread_t, int>(pthread_self(), i));
        // Add one to the total number of threads for the master thread
        nthread_++;

#ifdef _OPENMP
        // We need to disable OpenMP threads since MADNESS has its own.
        omp_set_num_threads(1);
#endif

    }

    MADNESSCommWrapper(const MADNESSCommWrapper &copy) { *this = copy; }

    ~MADNESSCommWrapper()
    { }

    MADNESSCommWrapper& operator=(const MADNESSCommWrapper& other)
    {
        me_   = other.me_;
        nproc_ = other.nproc_;
        nthread_ = other.nthread_;
        madworld_ = other.madworld_;
        thread_id_ = other.thread_id_;
        return *this;
    }

    inline int thread_id(const pthread_t &thread) { return thread_id_[thread]; }

    inline void sync() { madworld_->gop.fence(); }

    inline void print(FILE *out=outfile) const
    {
        if (me_ == 0) {
            fprintf(out, "\n    Using MadCommunicator (Number of procs = %d)\n", nproc_);
            fprintf(out, "                          (Number of threads in pool = %d)\n\n", nthread_);
        }
    }

    inline void finalize()
    {
        madness::finalize();

#if defined(HAVE_ELEMENTAL)
        elem::Finalize( );
#endif
    }

    template<typename type>
    inline void send(type *data, int nelem, int target, int tag=DEFAULT_SEND_RECV_TAG) const
    {
        madworld_->mpi.Isend(data, nelem, target, tag);
    }

    template<typename type>
    inline void recv(type *data, int nelem, int sender, int tag=DEFAULT_SEND_RECV_TAG) const
    {
        madworld_->mpi.Irecv(data, nelem, sender, tag);
    }

    template<typename type>
    inline void bcast(type* data, size_t nelem, int broadcaster=0)
    {
        madworld_->gop.broadcast(data, nelem, broadcaster);
    }

    template<typename type>
    inline void bcast_serializable(type &data, int broadcaster=0)
    {
        madworld_->gop.broadcast_serializable(data, broadcaster);
    }

    template<typename type>
    inline void sum(type *data, int nelem, type *receive_buffer=0, int target=-1)
    {
        madworld_->gop.sum(data, nelem);
    }

    inline int me() const { return me_; }
    inline int nproc() const { return nproc_; }
    inline int nthread() const { return nthread_; }
    inline std::string communicator() const {return communicator_;}

    inline SharedLock spinlock() { return SharedLock(new madness::Spinlock()); }

    inline SharedMutex mutex() { return SharedMutex(new madness::MutexFair()); }

    inline SharedMadWorld get_madworld() { return madworld_; }

}; // End of MADNESSCommWrapper class

} // End of psi namespace

#endif // End of HAVE_MADNESS

#endif // End of _psi_src_lib_libparallel_mad_h_
