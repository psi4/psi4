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

#if HAVE_MPI == 1
#include <mpi.h>
#endif

#if HAVE_MADNESS == 1

#ifdef WORLD_INSTANTIATE_STATIC_TEMPLATES
#undef WORLD_INSTANTIATE_STATIC_TEMPLATES
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#else
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#endif

#include <world/world.h>

#endif


namespace psi {

    extern FILE *outfile;

    class Serializable;

    class Communicator {
    protected:
        int me_;
        int nproc_;
        const std::string communicator_;

    public:

        Communicator(const std::string &communicator);
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
        };

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

        virtual void raw_send(const void* data, int nbyte, int target) = 0;
        virtual void raw_recv(void* data, int nbyte, int sender) = 0;
        virtual void raw_bcast(void* data, int nbyte, int broadcaster) = 0;
        virtual void raw_sum(double *data, int n, double *receive_buffer=0, int target=-1) = 0;
        virtual void raw_sum(unsigned int *data, int n, unsigned int *receive_buffer=0, int target=-1) = 0;
        virtual void raw_sum(int *data, int n, int *receive_buffer=0, int target=-1) = 0;
        virtual void raw_sum(char *data, int n, char *receive_buffer=0, int target=-1) = 0;
        virtual void raw_sum(long *data, int n, long *receive_buffer=0, int target=-1) = 0;

        const std::string get_comm() {return communicator_;}
        virtual int nthread() {return 1;};
        virtual void print(FILE *out=outfile) const {};

#if HAVE_MADNESS == 1
        virtual boost::shared_ptr<madness::World> get_madworld() {};
        virtual madness::Spinlock* get_mutex() {};
#endif
        virtual int get_threadid(pthread_t thr_) { return 0; };


    };
//    void p_fprintf(FILE * __restrict __stream, const char * __restrict __format, ...);

    class LocalCommunicator : public Communicator {
    private:
        const std::string communicator_;

    public:
        LocalCommunicator(const std::string &communicator);
        LocalCommunicator(const LocalCommunicator &copy);
        virtual ~LocalCommunicator();

        LocalCommunicator& operator=(const LocalCommunicator& other);

        virtual void raw_send(const void* data, int nbyte, int target);
        virtual void raw_recv(void* data, int nbyte, int sender);
        virtual void raw_bcast(void* data, int nbyte, int broadcaster);
        virtual void raw_sum(double *data, int n, double *receive_buffer=0, int target=-1);
        virtual void raw_sum(unsigned int *data, int n, unsigned int *receive_buffer=0, int target=-1);
        virtual void raw_sum(int *data, int n, int *receive_buffer=0, int target=-1);
        virtual void raw_sum(char *data, int n, char *receive_buffer=0, int target=-1);
        virtual void raw_sum(long *data, int n, long *receive_buffer=0, int target=-1);

        virtual void sync();

        //virtual void runTask () { if (where == me_) task_->runTask(); }
        virtual int nthread() {return 1;};
        virtual void print(FILE *out=outfile) const;

    };

#if HAVE_MPI == 1
    class MPICommunicator : public Communicator {
        MPI_Comm comm_;
        const std::string communicator_;

    public:
        MPICommunicator(MPI_Comm comm, const std::string &communicator);
        MPICommunicator(const MPICommunicator &copy);
        virtual ~MPICommunicator();

        MPICommunicator& operator=(const MPICommunicator& other);

        virtual void raw_send(const void* data, int nbyte, int target);
        virtual void raw_recv(void* data, int nbyte, int sender);
        virtual void raw_bcast(void* data, int nbyte, int broadcaster);
        virtual void raw_sum(double *data, int n, double *receive_buffer=0, int target=-1);
        virtual void raw_sum(unsigned int *data, int n, unsigned int *receive_buffer=0, int target=-1);
        virtual void raw_sum(int *data, int n, int *receive_buffer=0, int target=-1);
        virtual void raw_sum(char *data, int n, char *receive_buffer=0, int target=-1);
        virtual void raw_sum(long *data, int n, long *receive_buffer=0, int target=-1);

        virtual void sync();

        //virtual void runTask () { if (where == me_) task_->runTask(); }
        virtual int nthread() { return 1; }
        virtual void print(FILE *out=outfile) const;


    };
#endif

#if HAVE_MADNESS == 1

    class MadCommunicator : public Communicator {
    private:
        boost::shared_ptr<madness::World> madworld_;
        int nthread_;
        std::map<pthread_t, int> thread_id_;

        madness::Void set_thread_id();

        boost::shared_ptr<madness::Spinlock> mutex_;
        int thread_count;

    public:
        MadCommunicator(boost::shared_ptr<madness::World> madness_world, const std::string &communicator);
        virtual ~MadCommunicator();

        virtual void raw_send(const void* data, int nbyte, int target);
        virtual void raw_recv(void* data, int nbyte, int sender);
        virtual void raw_bcast(void* data, int nbyte, int sender);
        virtual void raw_sum(double *data, int n, double *receive_buffer=0, int target=-1);
        virtual void raw_sum(unsigned int *data, int n, unsigned int *receive_buffer=0, int target=-1);
        virtual void raw_sum(int *data, int n, int *receive_buffer=0, int target=-1);
        virtual void raw_sum(char *data, int n, char *receive_buffer=0, int target=-1);
        virtual void raw_sum(long *data, int n, long *receive_buffer=0, int target=-1);

        virtual void sync();

        virtual void print(FILE *out=outfile) const;
        //virtual void runTask () { task_->madTask(madworld_); }

        virtual int nthread() {return nthread_; };

        virtual madness::Spinlock* get_mutex() { return new madness::Spinlock(); };
        virtual boost::shared_ptr<madness::World> get_madworld() { return madworld_; };
        virtual int get_threadid(pthread_t thr_) { return thread_id_[thr_]; };

    };

#endif

}

#endif  /* _PARALLEL_H */
