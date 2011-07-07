
namespace MPI {

    struct Status {};

    struct Request {
        bool Test() {
            return false;
        }

        static bool Testany(int n, Request* request, int& ind) {
            return false;
        }

        static int Testsome(int n, Request* request, int* ind, MPI::Status* status) {
            return false;
        }

        static int Testsome(int n, Request* request, int* ind) {
            return false;
        }
    };

    struct Intracomm {
        int Get_rank() const {
            return 0;
        }

        int Get_size() const {
            return 1;
        }

        Request Isend(const void* buf, size_t count, const MPI::Datatype& datatype, int dest, int tag) const {
            throw "not implemented";
        }

        Request Irecv(void* buf, size_t count, const MPI::Datatype& datatype, int src, int tag) const {
            throw "not implemented";
        }

        void Send(const void* buf, size_t count, const MPI::Datatype& datatype, int dest, int tag) const {
            throw "not implemented";
        }

        void Recv(void* buf, int count, const MPI::Datatype& datatype, int source, int tag, MPI::Status& status) const {
            throw "not implemented";
        }

        void Recv(void* buf, int count, const MPI::Datatype& datatype, int source, int tag) const {
            throw "not implemented";
        }

        void Bcast(void* buf, size_t count, const MPI::Datatype& datatype, int root) const {
            return;
        }

        void Reduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype, const MPI::Op& op, int root) const {
            memcpy(recvbuf, sendbuf, count*sizeof the damn data type);
        }

        void Allreduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype, const MPI::Op& op) const {
            memcpy(recvbuf, sendbuf, count*sizeof the damn data type);
        }

        void Get_attr(int key, void* value) const {
            throw "not implemented";
        }

        void Abort(int code=1) const {
            exit(code);
        }

        void Barrier() const {
            return;
        }
    };
}
