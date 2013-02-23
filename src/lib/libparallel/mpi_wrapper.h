#if defined(HAVE_MPI)

#include <mpi.h>

namespace psi {

class MPICommunicator {
    int me_;
    int nproc_;
    int nthread_;
    std::string communicator_;
    MPI_Comm comm_;
public:
    MPICommunicator(const int &argc, char **argv)
    {
        MPI_Init(const_cast<int*>(&argc), &argv);

        comm_ = MPI_COMM_WORLD;
        MPI_Comm_rank(comm_, &me_);
        MPI_Comm_size(comm_, &nproc_);
        nthread_ = 1;
        communicator_ = "MPI";
    }

    MPICommunicator(const MPICommunicator &copy) :
        me_(copy.me_), nproc_(copy.nproc_),
        nthread_(copy.nthread_), communicator_(copy.communicator_),
        comm_(copy.comm_)
    { }

    ~MPICommunicator()
    { }

    MPICommunicator& operator =(const MPICommunicator& other)
    {
        if (this != &other) {
            me_   = other.me_;
            nproc_ = other.nproc_;
            nthread_ = other.nthread_;
            comm_ = other.comm_;
        }
        return *this;
    }

    void sync()
    {
        MPI_Barrier(comm_);
    }


    void raw_send(const void* data, int nbyte, int target)
    {
        MPI_Send(const_cast<void*>(data), nbyte, MPI_BYTE, target, 0, comm_);
    }

    void raw_recv(void* data, int nbyte, int sender)
    {
        MPI_Status status;
        MPI_Recv(data, nbyte, MPI_BYTE, sender, 0, comm_, &status);
    }

    void raw_bcast(void* data, int nbyte, int broadcaster)
    {
        MPI_Bcast(data, nbyte, MPI_BYTE, broadcaster, comm_);
    }

    template<typename type>
    inline void bcast(type *data, int nelem, int broadcaster) {
        raw_bcast(data, nelem * sizeof(type), broadcaster);
    }

    template<typename type>
    inline void bcast_serializable(type& data, int broadcaster) {

    }

    void sum(double* data, size_t nelem) {
        double *receive_buffer = new double[nelem];
        MPI_Allreduce(static_cast<void*>(data), static_cast<void*>(receive_buffer), nelem, MPI_DOUBLE, MPI_SUM, comm_);
        ::memcpy(static_cast<void*>(data), static_cast<void*>(receive_buffer), sizeof(double)*nelem);
        delete[] receive_buffer;
    }

#define SUMMEMBER(T, M) \
    void raw_sum(T *data, int n, T *receive_buffer, int target) \
    { \
        bool alloc = false; \
        if (receive_buffer == NULL) { \
            alloc = true; \
            receive_buffer = new T[n]; \
        } \
     \
        if (target >= 0) \
            MPI_Reduce(static_cast<void*>(data), static_cast<void*>(receive_buffer), n, M, MPI_SUM, target, comm_); \
        else \
            MPI_Allreduce(static_cast<void*>(data), static_cast<void*>(receive_buffer), n, M, MPI_SUM, comm_); \
     \
        if (alloc) { \
            ::memcpy(static_cast<void*>(data), static_cast<void*>(receive_buffer), sizeof(T)*n); \
            delete[] receive_buffer; \
        } \
    }

    SUMMEMBER(double, MPI_DOUBLE)
    SUMMEMBER(unsigned int, MPI_INT)
    SUMMEMBER(int, MPI_INT)
    SUMMEMBER(char, MPI_CHAR)
    SUMMEMBER(long, MPI_LONG)

    void print(FILE *out) const
    {
        if (me() == 0) {
            fprintf(out, "\n    Using MPICommunicator (Number of processes = %d)\n\n", nproc());
        }
    }

    void finalize() {
        MPI_Finalize();
    }

    inline int me() const { return me_; }
    inline int nproc() const { return nproc_; }
    inline int nthread() const { return nthread_; }
    inline int thread_id(pthread_t) const { return 0; }
    inline const std::string& communicator() const { return communicator_; }
};

}

#endif
