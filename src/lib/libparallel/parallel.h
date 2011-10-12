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
//#include <psi4-dec.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_MADNESS

#ifdef WORLD_INSTANTIATE_STATIC_TEMPLATES
#undef WORLD_INSTANTIATE_STATIC_TEMPLATES
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#else
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#endif

#include <world/world.h>
#include <world/worldobj.h>
#include <world/archive.h>
#include <world/safempi.h>

#endif


namespace psi {

    extern FILE *outfile;

    class Serializable;

#ifdef HAVE_MADNESS
    typedef boost::shared_ptr<madness::Spinlock> SharedLock;
    typedef boost::shared_ptr<madness::MutexFair> SharedMutex;
    typedef boost::shared_ptr<madness::World> SharedMadWorld;
#endif

    class Communicator {
    protected:
        int me_;
        int nproc_;
        int nthread_;
        std::string communicator_;

        /**
         * The raw send, recv, bcast, and sums will be overwritten by
         * a given communicator.
         */
        virtual void raw_send(const void* data, int nbyte, int target) = 0;
        virtual void raw_recv(void* data, int nbyte, int sender) = 0;
        virtual void raw_bcast(void* data, int nbyte, int broadcaster) = 0;
        virtual void raw_sum(double *data, int n, double *receive_buffer=0, int target=-1) = 0;
        virtual void raw_sum(unsigned int *data, int n, unsigned int *receive_buffer=0, int target=-1) = 0;
        virtual void raw_sum(int *data, int n, int *receive_buffer=0, int target=-1) = 0;
        virtual void raw_sum(char *data, int n, char *receive_buffer=0, int target=-1) = 0;
        virtual void raw_sum(long *data, int n, long *receive_buffer=0, int target=-1) = 0;
//        virtual void raw_put(const void* data, int nbyte, int target) = 0;
//        virtual void raw_get(const void* data, int nbyte, int target) = 0;


    public:

        Communicator();
        virtual ~Communicator();

        static boost::shared_ptr<Communicator> world;


        /**
         * Returns the local processor number in the communicator.
         * @return Integer between 0 and n, the number of processors
         */
        virtual int me() const { return me_; }

        /**
         * Returns the number of processors in the communicator
         * @return The number of processors
         */
        virtual int nproc() const { return nproc_; }

        /**
         * Return the total number of threads
         * @return The number of threads
         */
        virtual int nthread() const { return nthread_; }

        /**
         * Return the communicator type
         * @return The specific communicator type
         */
        virtual std::string communicator() const {return communicator_;}

        /**
         * Return a unique integer value for a given pthread
         * @param thread The thread id of a given pthread
         */
        virtual int thread_id(const pthread_t &thread) = 0;


        /**
         * This function is a barrier that makes sure that all communication is complete
         * before continuing.
         */
        virtual void sync() = 0;

        /**
         * Send array of data to remote node.
         * @param target The node number to send to
         * @param data The array of data to send
         * @param nelem The size of the array.  This is the number of elements, not the number of bytes.
         */
        template<typename type>
        void send(const type *data, int nelem, int target) {
            raw_send(static_cast<const void *>(data), nelem*sizeof(type), target);
        }

        /**
         * Receive array of data from remote node. You must allocate memory before calling this function.
         * @param sender The node number that sent it
         * @param data The destination array to receive data
         * @param nelem The size of the array.  This is the number of elements, not the number of bytes.
         */
        template<typename type>
        void recv(type *data, int nelem, int sender) {
            raw_recv(static_cast<void *>(data), nelem*sizeof(type), sender);
        }


        /**
         * Broadcast data to all nodes. At the end of the call all nodes have a complete copy of the data sent.
         * @param data The array of data to be broadcasted.
         * @param nelem The size of the array. This is the number of elements, not the number of bytes.
         * @broadcaster The node sending the data.
         */
        template<typename type>
        void bcast(type *data, int nelem, int broadcaster=0) {
            raw_bcast(static_cast<void *>(data), nelem*sizeof(type), broadcaster);
        }

        /**
         * Broadcast data to all nodes. At the end of the call all nodes have a complete copy of the data sent.
         * @param data The array of data to be broadcasted.
         * @param nelem The size of the array. This is the number of elements, not the number of bytes.
         * @broadcaster The node sending the data.
         */
        void bcast(std::string& data, int broadcaster=0);

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
        void sum(type *data, int nelem, type *receive_buffer=0, int target=-1) {
            raw_sum(data, nelem, receive_buffer, target);
        }

        /**
         * Print some communicator information
         * @param out The location where the communicator information is written.
         */
        virtual void print(FILE *out=outfile) const = 0;

        /**
         * Call the correct finalize for a given communicator type
         */
        virtual void finalize() = 0;

//        /**
//         * Put an array of data to remote node.
//         * @param target The node number to send to
//         * @param data The array of data to send
//         * @param nelem The size of the array.  This is the number of elements, not the number of bytes.
//         */
//        template<typename type>
//        void put(const type *data, int nelem, int target) {
//            raw_put(static_cast<const void *>(data), nelem*sizeof(type), target);
//         }


//        /**
//         * Get array of data to remote node.
//         * @param target The node number to send to
//         * @param data The array of data to send
//         * @param nelem The size of the array.  This is the number of elements, not the number of bytes.
//         */
//        template<typename type>
//        void get(const type *data, int nelem, int target) {
//            raw_get(static_cast<const void *>(data), nelem*sizeof(type), target);
//         }


#ifdef HAVE_MADNESS
        virtual SharedLock get_spinlock() { };
        virtual SharedMutex get_mutex() { };
        virtual SharedMadWorld get_madworld() = 0;
#endif
    };

    typedef boost::shared_ptr<Communicator> SharedComm;

//    void p_fprintf(FILE * __restrict __stream, const char * __restrict __format, ...);

    class LocalCommunicator : public Communicator {
    protected:
#ifdef HAVE_MADNESS
        SharedMadWorld madworld_;
#endif

        virtual void raw_send(const void* data, int nbyte, int target);
        virtual void raw_recv(void* data, int nbyte, int sender);
        virtual void raw_bcast(void* data, int nbyte, int broadcaster);
        virtual void raw_sum(double *data, int n, double *receive_buffer=0, int target=-1);
        virtual void raw_sum(unsigned int *data, int n, unsigned int *receive_buffer=0, int target=-1);
        virtual void raw_sum(int *data, int n, int *receive_buffer=0, int target=-1);
        virtual void raw_sum(char *data, int n, char *receive_buffer=0, int target=-1);
        virtual void raw_sum(long *data, int n, long *receive_buffer=0, int target=-1);
//        virtual void raw_put(const void* data, int nbyte, int target) { throw PSIEXCEPTION("don't do that"); }
//        virtual void raw_get(const void* data, int nbyte, int target) { throw PSIEXCEPTION("don't do that"); }


    public:
        LocalCommunicator(const int &argc, char **argv);
        LocalCommunicator(const LocalCommunicator &copy);
        virtual ~LocalCommunicator();

        LocalCommunicator& operator=(const LocalCommunicator& other);

        virtual int me() const { return me_; }
        virtual int nproc() const { return nproc_; }
        virtual int nthread() const { return nthread_; }
        virtual std::string communicator() const {return communicator_;}

        virtual int thread_id(const pthread_t &thread) { return 0; }

        virtual void sync();

        virtual void print(FILE *out=outfile) const;

        virtual void finalize();

#ifdef HAVE_MADNESS
        virtual SharedMadWorld get_madworld() { return madworld_; }
#endif
    };

#ifdef HAVE_MADNESS

    class MadCommunicator : public Communicator
    {
    protected:
        SharedMadWorld madworld_;
        std::map<pthread_t, int> thread_id_;

        virtual void raw_send(const void* data, int nbyte, int target);
        virtual void raw_recv(void* data, int nbyte, int sender);
        virtual void raw_bcast(void* data, int nbyte, int sender);
        virtual void raw_sum(double *data, int n, double *receive_buffer=0, int target=-1);
        virtual void raw_sum(unsigned int *data, int n, unsigned int *receive_buffer=0, int target=-1);
        virtual void raw_sum(int *data, int n, int *receive_buffer=0, int target=-1);
        virtual void raw_sum(char *data, int n, char *receive_buffer=0, int target=-1);
        virtual void raw_sum(long *data, int n, long *receive_buffer=0, int target=-1);
//        virtual void raw_put(const void* data, int nbyte, int target);
//        virtual void raw_get(const void* data, int nbyte, int target);


    public:
        MadCommunicator(const int &argc, char **argv);
        MadCommunicator(const MadCommunicator &copy);
        virtual ~MadCommunicator();

        MadCommunicator& operator=(const MadCommunicator& other);

        virtual int thread_id(const pthread_t &thread) { return thread_id_[thread]; }

        virtual void sync();

        virtual void print(FILE *out=outfile) const;

        virtual void finalize() { madness::finalize(); }

        virtual SharedLock get_spinlock() { return SharedLock(new madness::Spinlock()); }

        virtual SharedMutex get_mutex() { return SharedMutex(new madness::MutexFair()); }

        virtual SharedMadWorld get_madworld() { return madworld_; }

    };

#endif

}

#endif  /* _PARALLEL_H */
