/* 
 * File:   parallel.h
 * Author: jturney
 *
 * Created on December 11, 2009, 3:34 PM
 */

#ifndef _psi_src_lib_libutil_parallel_h_
#define	_psi_src_lib_libutil_parallel_h_

#include <boost/shared_ptr.hpp>
#include <mpi.h>

namespace psi {


    class Communicator {
    protected:
        int me_;
        int nproc_;
    public:

        Communicator();
        virtual ~Communicator();

        static boost::shared_ptr<Communicator> world;

        /**
         * Returns the local processor number in the communicator.
         * @return Integer between 0 and n, the number of processors
         */
        int me() const { return me_; }

        /**
         * Returns the number of processors in the communicator
         * @return The number of processors
         */
        int nproc() const { return nproc_; }

        virtual void sync() = 0;

        /**
         * Send array of data to remote node.
         * @param target The node number to send to
         * @param data The array of data to send
         * @param ndata The size of the array.  This is the number of elements, not the number of bytes.
         */
        virtual void send(int target, const double *data, int ndata);
        virtual void send(int target, const unsigned int *data, int ndata);
        virtual void send(int target, const int *data, int ndata);
        virtual void send(int target, const char *data, int ndata);
        virtual void send(int target, const long *data, int ndata);
        void send(int target, double data);
        void send(int target, int data);
        virtual void raw_send(int target, const void *data, int nbyte) = 0;

        /**
         * Receive array of data from remote node. You must allocate memory before calling this function.
         * @param sender The node number that sent it
         * @param data The destination array to receive data
         * @param ndata The size of the array.  This is the number of elements, not the number of bytes.
         */
        virtual void recv(int sender, double *data, int ndata);
        virtual void recv(int sender, unsigned int *data, int ndata);
        virtual void recv(int sender, int *data, int ndata);
        virtual void recv(int sender, char *data, int ndata);
        virtual void recv(int sender, long *data, int ndata);
        void recv(int sender, double& data);
        void recv(int sender, int& data);
        virtual void raw_recv(int sender, void *data, int nbyte) = 0;

        /**
         * Broadcast data to all nodes. At the end of the call all nodes have a complete copy of the data sent.
         * @param data The array of data to be broadcasted.
         * @param ndata The size of the array. This is the number of elements, not the number of bytes.
         * @broadcaster The node sending the data.
         */
        virtual void bcast(double *data, int ndata, int broadcaster=0);
        virtual void bcast(unsigned int *data, int ndata, int broadcaster=0);
        virtual void bcast(int *data, int ndata, int broadcaster=0);
        virtual void bcast(char *data, int ndata, int broadcaster=0);
        virtual void bcast(long *data, int ndata, int broadcaster=0);
        void bcast(double &data, int broadcaster=0);
        void bcast(int &data, int broadcaster=0);
        virtual void raw_bcast(void *data, int nbyte, int broadcaster=0);

        /**
         * Performs element-by-element sum of all data from all nodes.  The sum will either appear
         * in a new buffer or will overwrite the original data.
         * @param data The array of data to be summed.  If no receive buffer is given, this is overwritten with summed data.
         * @param n The size of the array. This is the number of elements, not the number of bytes.
         * @param receive_buffer Optional receive buffer. If given, summed data is placed in this array.
         *                       If used, this must be allocated before entering method.
         * @param target Data is summed on target. By default all nodes receive data.
         */
        virtual void sum(double *data, int n, double *receive_buffer=0, int target=-1);
        virtual void sum(unsigned int *data, int n, unsigned int *receive_buffer=0, int target=-1);
        virtual void sum(int *data, int n, int *receive_buffer=0, int target=-1);
        virtual void sum(char *data, int n, char *receive_buffer=0, int target=-1);
        virtual void sum(long *data, int n, long *receive_buffer=0, int target=-1);
        void sum(double &data);
        void sum(int &data);

    };

    class MPICommunicator : public Communicator {
        MPI_Comm comm_;

    public:
        MPICommunicator(MPI_Comm comm);
        virtual ~MPICommunicator();

        virtual void sync();
        
        virtual void raw_send(int target, const void *data, int nbyte);
        virtual void raw_recv(int sender, void *data, int nbyte);
        virtual void raw_bcast(void *data, int nbyte, int broadcaster=0);

        virtual void sum(double *data, int n, double *receive_buffer=0, int target=-1);
        virtual void sum(unsigned int *data, int n, unsigned int *receive_buffer=0, int target=-1);
        virtual void sum(int *data, int n, int *receive_buffer=0, int target=-1);
        virtual void sum(char *data, int n, char *receive_buffer=0, int target=-1);
        virtual void sum(long *data, int n, long *receive_buffer=0, int target=-1);
    };

    class LocalCommunicator : public Communicator {
    };

}

#endif	/* _PARALLEL_H */

