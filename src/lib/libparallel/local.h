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

#include <psi4-dec.h>
#include <cstring>
#include <libparallel/openmp.h>

namespace psi {

class LocalCommWrapper {

private:
    int me_;
    int nproc_;
    int nthread_;
    std::string communicator_;

public:

    LocalCommWrapper(int &argc, char **argv)
    {
        init_openmp();
        me_ = 0;
        nproc_ = 1;
        nthread_ = Process::environment.get_n_threads();
        communicator_ = "LOCAL";
    }

    ~LocalCommWrapper() { }

    /**
     * Returns the local processor number in the communicator.
     * @return Integer between 0 and n, the number of processors
     */
    inline int me() const { return me_; }

    /**
     * Returns the number of processors in the communicator
     * @return The number of processors
     */
    inline int nproc() const { return nproc_; }

    /**
     * Return the total number of threads
     * @return The number of threads
     */
    inline int nthread() const { return nthread_; }

    /**
     * Return the communicator type
     * @return The specific communicator type
     */
    inline std::string communicator() const {return communicator_;}

    /**
     * Return a unique integer value for a given pthread
     * @param thread The thread id of a given pthread
     */
    inline int thread_id(const pthread_t &thread) { return 0; }

    /**
     * This function is a barrier that makes sure that all communication is complete
     * before continuing.
     */
    inline void sync() { }

    /**
     * Send array of data to remote node.
     * @param target The node number to send to
     * @param data The array of data to send
     * @param nelem The size of the array.  This is the number of elements, not the number of bytes.
     */
    template<typename type>
    inline void send(const type data, int nelem, int target) { }

    /**
     * Receive array of data from remote node. You must allocate memory before calling this function.
     * @param sender The node number that sent it
     * @param data The destination array to receive data
     * @param nelem The size of the array.  This is the number of elements, not the number of bytes.
     */
    template<typename type>
    inline void recv(type data, int nelem, int sender) { }


    /**
     * Broadcast data to all nodes. At the end of the call all nodes have a complete copy of the data sent.
     * @param data The array of data to be broadcasted.
     * @param nelem The size of the array. This is the number of elements, not the number of bytes.
     * @broadcaster The node sending the data.
     */
    template<typename type>
    inline void bcast(type data, int nelem, int broadcaster=0) { }

    /**
     * Broadcast serializable data to all nodes (i.e. vector, string).
     * At the end of the call all nodes have a complete copy of the data sent.
     * @param data The array of data to be broadcasted.
     * @param nelem The size of the array. This is the number of elements, not the number of bytes.
     * @broadcaster The node sending the data.
     */
    template<typename type>
    inline void bcast_serializable(type &data, int broadcaster=0) { }

    /**
     * Performs element-by-element sum of all data from all nodes.  The sum will either appear
     * in a new buffer or will overwrite the original data.
     * @param data The array of data to be summed.  If no receive buffer is given, this is overwritten with summed data.
     * @param n The size of the array. This is the number of elements, not the number of bytes.
     * @param receive_buffer Optional receive buffer. If given, summed data is placed in this array.
     *                       If used, this must be allocated before entering method.
     * @param target Data is summed on target. By default all nodes receive data.
     */
    template<typename type>
    inline void sum(type data, int nelem, type *receive_buffer=0, int target=-1)
    {
        if (receive_buffer != 0)
            ::memcpy(receive_buffer, data, sizeof(type) * nelem);
    }

    /**
     * Print some communicator information
     * @param out The location where the communicator information is written.
     */
    inline void print(FILE *out) const
    {
        fprintf(out, "\n    Using LocalCommunicator (Number of processes = 1)\n\n");
    }

    /**
     * Call the correct finalize for a given communicator type
     */
    inline void finalize()
    { }

//#ifdef HAVE_MADNESS
//    SharedLock get_spinlock() { }
//    SharedMutex get_mutex() { }
//    SharedMadWorld get_madworld() { }
//#endif

}; // End of LocalCommWrapper class


} // End of namespace psi4

#endif // End of _psi_src_lib_libparallel_local_h_
