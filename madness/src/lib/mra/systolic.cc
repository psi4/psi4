/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/
#include <world/world.h>

#include <utility>
#include <tensor/tensor.h>

namespace madness {

    /// Manages data associated with a row/column/block distributed array

    /// The class itself provides limited functionality for accessing the
    /// data and is primarily intended to provide base functionality for
    /// use by matrix algorithms and other matrix classes.
    ///
    /// The constructor is deliberately simple.  Factory functions
    /// are expected to be the main construction tool.
    template <typename T>
    class DistributedMatrix {
        World& _world;
        const int64_t P;                //< No. of processors
        const ProcessID rank;           //< My processor rank
        const int64_t n;                //< Column dimension of A(n,m)
        const int64_t m;                //< Row dimension of A(n,m)
        const int64_t tilen;            //< Tile size for column
        const int64_t tilem;            //< Tile size for row
        const int64_t Pcoldim;          //< Column dimension of processor grid
        const int64_t Prowdim;          //< Row dimension of processor grid
        const int64_t Pcol;             //< Column of processor grid for this processor
        const int64_t Prow;             //< Row of processor grid for this processor
        const int64_t ilo,ihi;          //< Range of column indices on this processor
        const int64_t jlo,jhi;          //< Range of row indices on this processor
        const int64_t idim,jdim;        //< Dimension of data on this processor

        Tensor<T> t;                    //< The data

    public:
        DistributedMatrix(World& world, int64_t n, int64_t m, int64_t coltile, int64_t rowtile)
            : _world(world)
            , P(world.size())
            , rank(world.rank())
            , n(n)
            , m(m)
            , tilen(coltile)
            , tilem(rowtile)
            , Pcoldim((n-1)/tilen+1)
            , Prowdim((m-1)/tilem+1)
            , Pcol(rank/Prowdim)
            , Prow(rank - Pcol*Prowdim)
            , ilo(Pcol*tilen)
            , ihi(std::min(ilo+tilen-1,n-1))
            , jlo(Prow*tilem)
            , jhi(std::min(jlo+tilem-1,m-1))
            , idim(std::max(ihi-ilo+1,int64_t(0)))
            , jdim(std::max(jhi-jlo+1,int64_t(0)))
        {
            //print("DM: dims", n, m, "tiles", tilen, tilem, "ilo ihi jlo jhi", ilo, ihi, jlo, jhi, "idim jdim", idim, jdim);

            if (idim>0 && jdim>0) t = Tensor<T>(idim,jdim);
        }

        virtual ~DistributedMatrix() {}

        /// Returns the column dimension of the matrix ... i.e., n for A(n,m)
        int64_t coldim() const {
            return n;
        }

        /// Returns the row dimension of the matrix ... i.e., m for A(n,m)
        int64_t rowdim() const {
            return m;
        }

        /// Returns the column tile size
        int64_t coltile() const {
            return tilen;
        }

        /// Returns the row tile size
        int64_t rowtile() const {
            return tilem;
        }

        /// Returns the no. of processors in the column dimension
        int64_t process_coldim() const {return Pcoldim;}

        /// Returns the no. of processors in the row dimension
        int64_t process_rowdim() const {return Prowdim;}

        /// Returns the total no. of elements stored on this processor
        int64_t local_size() const {return idim*jdim;}

        /// Returns the no. of column elements stored on this processor
        int64_t local_coldim() const {return idim;}

        /// Returns the no. of row elements stored on this processor
        int64_t local_rowdim() const {return jdim;}

        /// Returns the inclusive range of column indices on this processor

        /// If there is no data on this processor it returns ilow=0 and ihigh=-1
        void local_colrange(int64_t& ilow, int64_t& ihigh) const {
            if (ilo <= ihi) {
                ilow = ilo;
                ihigh = ihi;
            }
            else {
                ilow = 0;
                ihigh = -1;
            }
        }

        /// Returns the inclusive range of row indices on this processor

        /// If there is no data on this processor it returns jlow=0 and jhigh=-1
        void local_rowrange(int64_t& jlow, int64_t& jhigh) const {
            if (jlo <= jhi) {
                jlow = jlo;
                jhigh = jhi;
            }
            else {
                jlow = 0;
                jhigh = -1;
            }
        }

        /// Returns the inclusive range of column indices on processor p

        /// If there is no data on this processor it returns ilow=0 and ihigh=-1
        void get_colrange(int p, int64_t& ilow, int64_t& ihigh) const {
            int pi = p/Prowdim;
            int pj = p - pi*Prowdim;
            if (pi >= process_coldim() || pj >= process_rowdim()) {
                ilow = 0;
                ihigh = -1;
            }
            else {
                ilow = pi*tilen;
                ihigh= std::min(ilow+tilen-1,n-1);
            }
            return;
        }

        /// Returns associated world
        World& get_world() {return _world;}

        /// Returns reference to data
        Tensor<T>& data() {return t;}

        /// Returns const reference to data
        const Tensor<T>& data() const {return t;}

        /// Returns true if the matrix is column distributed (i.e., row dimension not distributed)
        bool is_column_distributed() const {return process_rowdim()==1;}

        /// Returns true if the matrix is row distributed (i.e., column dimension not distributed)
        bool is_row_distributed() const {return process_coldim()==1;}

        /// Given the full matrix(n,m), copy in the data range local to this processor
        void copyin(const Tensor<T>& s) {
            if (local_size() > 0) t(___) = s(Slice(ilo,ihi),Slice(jlo,jhi));
        }

        /// Given the full matrix s(n,m), copy out the data range local to this processor
        void copyout(Tensor<T>& s) const {
            if (local_size() > 0) s(Slice(ilo,ihi),Slice(jlo,jhi)) = t(___);
        }
    };


    /// Generates an (n,m) matrix distributed by columns (row dimension is not distributed)

    /// Quietly forces an even column tile size for ease of use in the systolic matrix algorithms
    template <typename T>
    DistributedMatrix<T> column_distributed_matrix(World& world, int64_t n, int64_t m, int64_t coltile=0) {
        if (world.size()*coltile < n) coltile = (n-1)/world.size() + 1;
        if ((coltile&0x1)) ++coltile;
        coltile = std::min(coltile,n);

        return DistributedMatrix<T>(world, n, m, coltile, m);
    }


    /// Generates an (n,m) matrix distributed by rows (column dimension is not distributed)
    template <typename T>
    DistributedMatrix<T> row_distributed_matrix(World& world, int64_t n, int64_t m, int64_t rowtile=0) {
        if (world.size()*rowtile < n) rowtile = (n-1)/world.size() + 1;
        rowtile = std::min(rowtile,m);

        return DistributedMatrix<T>(world, n, m, n, rowtile);
    }


    /// Generates a distributed matrix with rows of \c a and \c b interleaved

    /// I.e., the even rows of the result will be rows of \c a , and the
    /// odd rows those of \c b .
    ///
    /// The matrices a and b must have the same dimensions and be
    /// identically distributed.  The result will have a doubled column
    /// dimension and column tile size.  The row dimension is unchanged.
    template <typename T>
    DistributedMatrix<T> interleave_rows(const DistributedMatrix<T>& a, const DistributedMatrix<T>& b) {
        MADNESS_ASSERT(a.rowdim()==b.rowdim() && a.coldim()==b.coldim() && a.coltile()==b.coltile() && a.rowtile()==b.rowtile());

        DistributedMatrix<T> c(a.get_world(), a.coldim()*2, a.rowdim(), a.coltile()*2, a.rowtile());
        c.data()(Slice(0,-1,2),_) = a.data()(___);
        c.data()(Slice(1,-1,2),_) = b.data()(___);
    }


    /// Generates a distributed matrix with rows of \c a and \c b contatenated

    /// I.e., c[i,j] = a[i,j] if n<na or b[i,j] if j>=na
    ///
    /// The matrices a and b must have the same column size (i.e., the
    /// same number of rows) and be column distributed with the same
    /// column tilesze.  The result is also column distributed with
    /// the same column tilesize as the input matrices.
    template <typename T>
    DistributedMatrix<T> concatenate_rows(const DistributedMatrix<T>& a, const DistributedMatrix<T>& b) {
        MADNESS_ASSERT(a.coldim()==b.coldim() && a.coltile()==b.coltile() && a.is_column_distributed() && b.is_column_distributed());

        int64_t ma = a.rowdim();
        int64_t mb = b.rowdim();

        DistributedMatrix<T> c(a.get_world(), a.coldim(), ma+mb, a.coltile(), ma+mb);
        c.data()(_,Slice(0,ma-1)) = a.data()(___);
        c.data()(_,Slice(ma,-1))  = b.data()(___);
    }


    /// Base class for parallel algorithms that employ a systolic loop to generate all row pairs in parallel
    template <typename T>
    class SystolicMatrixAlgorithm : public TaskInterface {
    private:
        DistributedMatrix<T>& A;
        const int64_t nproc;            //< No. of processes with rows of the matrix (not size of world)
        const int64_t coldim;           //< A(coldim,rowdim)
        const int64_t rowdim;           //< A(coldim,rowdim)
        const int64_t nlocal;           //< No. of local pairs
        const ProcessID rank;           //< Rank of current process
        const int tag;                  //< MPI tag to be used for messages
        std::vector<T*> iptr, jptr;     //< Indirection for implementing cyclic buffer !! SHOULD BE VOLATILE ?????
        std::vector<int64_t> map;       //< Used to keep track of actual row indices

        using PoolTaskInterface::run;

        void iteration(const TaskThreadEnv& env) {
            start_iteration_hook(env);
            env.barrier();

            int64_t ilo, ihi;
            A.local_colrange(ilo, ihi);

            int neven = coldim + (coldim&0x1);

            int pairlo = rank*A.coltile()/2;

            int threadid = env.id();
            int nthread = env.nthread();

            for (int loop=0; loop<(neven-1); ++loop) {

                // This loop is parallelized over threads
                for (int pair=env.id(); pair<nlocal; pair+=nthread) {

                    int rp = neven/2-1-(pair+pairlo);
                    int iii = (rp+loop)%(neven-1);
                    int jjj = (2*neven-2-rp+loop)%(neven-1);
                    if (rp == 0) jjj = neven-1;

                    iii = map[iii];
                    jjj = map[jjj];

                    if (jptr[pair]) {
                        kernel(iii, jjj, iptr[pair], jptr[pair]);
                    }
                }
                env.barrier();

                if (threadid == 0) cycle();

                env.barrier();
            }
        }

        /// Call this after iterating to restore correct order of rows in original matrix

        /// At the end of each iteration the matrix rows are logically back in
        /// their correct order.  However, due to indirection to reduce data motion,
        /// if the local column dimension is not a factor of the number of cycles
        /// the underlying data may be in a different order.  This restores sanity.
        ///
        /// Only one thread should invoke this routine
        void unshuffle() {
            if (nlocal <= 0) return;
            Tensor<T>& t = A.data();
            Tensor<T> tmp(2L, t.dims(), false);
            T* tp = tmp.ptr();
            for (int64_t i=0; i<nlocal; ++i) {
                memcpy(tp+i*rowdim, iptr[i], rowdim*sizeof(T));
                if (jptr[i]) {
                    memcpy(tp+(i+nlocal)*rowdim, jptr[i], rowdim*sizeof(T));
                }
                iptr[i] = &t(i,0);
                jptr[i] = &t(i+nlocal,0);
            }
            memcpy(t.ptr(), tmp.ptr(), t.size()*sizeof(T));

            if (rank==(nproc-1) && (coldim&0x1)) jptr[nlocal-1] = 0;
        }

        /// Cycles data around the loop ... only one thread should invoke this
        void cycle() {
            if (coldim <= 2) return; // No cycling necessary
            if (nlocal <= 0) {       // Nothing local
                MADNESS_ASSERT(rank >= nproc);
                return;
            }

            // Check assumption that tiling put incomplete tile at the end
            MADNESS_ASSERT(A.local_coldim() == A.coltile()  ||  rank == (nproc-1));

            const ProcessID left = rank-1; //Invalid values are not used
            const ProcessID right = rank+1;

            /*
              Consider matrix (10,*) distributed with coltile=4 over
              three processors.

              .   0 1 2 3      4 5 6 7      8 9

              This is divided up as follows into this initial
              configuration for the loop

              .            P=0          P=1         P=2
              .                  msg          msg
              .   i    -->0-->1  -->   4-->5  -->    8  -->
              .       ^                                   |  msg
              .       |                         <---------
              .   j    <--2<--3  <--   6<--7  <--|   9
              msg          msg

              The first and last processes in the loop have to wrap ... others
              just pass left and right.  Note that 9 stays put.

              Note that the algorithm is assuming distribution puts equal
              amount of data on all nodes except the last.

              The i data is considered as flowing to the right.
              The j data is considered as flowing to the left.


              Hence, we should explore the pairs in this order
              (n-1 sets of n/2 pairs)

              .          P=0         P=1        P=2
              .          0  1        4  5       8
              .          2  3        6  7       9

              .          2  0        1  4       5
              .          3  6        7  8       9

              .          3  2        0  1       4
              .          6  7        8  5       9

              .          6  3        2  0       1
              .          7  8        5  4       9

              .          7  6        3  2       0
              .          8  5        4  1       9

              .          8  7        6  3       2
              .          5  4        1  0       9

              .          5  8        7  6       3
              .          4  1        0  2       9

              .          4  5        8  7       6
              .          1  0        2  3       9

              .          1  4        5  8       7
              .          0  2        3  6       9
            */

            // Copy end elements before they are overwritten
            T* ilast  = iptr[nlocal-1];
            T* jfirst = jptr[0];

            // Cycle local pointers
            for (int64_t i=0; i<nlocal-1; ++i) {
                iptr[nlocal-i-1] = iptr[nlocal-i-2];
                jptr[i] = jptr[i+1];
            }

            World& world = A.get_world();

            // Lazily program unsafely assuming enough MPI buffering available
            if (nproc == 1) {
                iptr[0] = jfirst;
                jptr[nlocal-2] = ilast;
            }
            else if (rank == 0) {
                iptr[0] = jfirst;
                world.mpi.Send(ilast, rowdim, right, tag);
                jptr[nlocal-1] = ilast;
                world.mpi.Recv(ilast, rowdim, right, tag);
            }
            else if (rank == (nproc-1)) {
                if (nlocal > 1) {
                    iptr[0] = jfirst;
                    jptr[nlocal-2] = ilast;
                }
                world.mpi.Send(iptr[0], rowdim, left, tag);
                world.mpi.Recv(iptr[0], rowdim, left, tag);
            }
            else {
                world.mpi.Send( ilast, rowdim, right, tag);
                world.mpi.Send(jfirst, rowdim,  left, tag);
                world.mpi.Recv( ilast, rowdim, right, tag);
                world.mpi.Recv(jfirst, rowdim,  left, tag);
                iptr[0] = jfirst;
                jptr[nlocal-1] = ilast;
            }
        }


    public:
        /// A must be a column distributed matrix with an even column tile >= 2
        SystolicMatrixAlgorithm(DistributedMatrix<T>& A, int tag, int nthread=ThreadPool::size()+1)
            : A(A)
            , nproc(A.process_coldim()*A.process_rowdim())
            , coldim(A.coldim())
            , rowdim(A.rowdim())
            , nlocal((A.local_coldim()+1)/2)
            , rank(A.get_world().rank())
            , tag(tag)
            , iptr(nlocal)
            , jptr(nlocal)
            , map(coldim+(coldim&0x1))
        {
            TaskInterface::set_nthread(nthread);

            MADNESS_ASSERT(A.is_column_distributed() && (nproc==1 || (A.coltile()&0x1)==0));

            // Initialize vectors of pointers to matrix rows
            Tensor<T>& t = A.data();

            //madness::print(nproc, coldim, rowdim, nlocal, rank, tag);

            for (int64_t i=0; i<nlocal; ++i) {
                iptr[i] = &t(i,0);
                jptr[i] = &t(i+nlocal,0);
            }

            // If no. of rows is odd, last process should have an empty last row
            if (rank==(nproc-1) && (coldim&0x1)) jptr[nlocal-1] = 0;

            // Initialize map from logical index order to actual index order

            int neven = (coldim+1)/2;
            int ii=0;
            for (ProcessID p=0; p<nproc; ++p) {
                int64_t lo, hi;
                A.get_colrange(p, lo, hi);
                int p_nlocal = (hi - lo + 2)/2;
                //print("I think process",p,"has",lo,hi,p_nlocal);
                for (int i=0; i<p_nlocal; ++i) {
                    map[ii+i] = lo+i;
                    //map[coldim-ii-nlocal+i] = lo+i+nlocal;
                    map[ii+i+neven] = lo+i+p_nlocal;
                }
                ii += p_nlocal;
            }

            std::reverse(map.begin(),map.begin()+neven);

            //print("MAP", map);
        }

        virtual ~SystolicMatrixAlgorithm() {}

        /// Threadsafe routine to apply the operation to rows i and j of the matrix
        virtual void kernel(int i, int j, T* rowi, T* rowj) = 0;


        /// Invoked simultaneously by all threads after each sweep to test for convergence

        /// There is a thread barrier before and after the invocation of this routine
        virtual bool converged(const TaskThreadEnv& env) const = 0;


        /// Invoked by all threads at the start of each iteration

        /// There is a thread barrier before and after the invocation of this routine
        virtual void start_iteration_hook(const TaskThreadEnv& env) {}


        /// Invoked by the task queue to run the algorithm with multiple threads
        void run(World& world, const TaskThreadEnv& env) {
            if (nlocal <= 0) return; // Nothing to do

            do {
                env.barrier();
                iteration(env);
                env.barrier();
            } while (!converged(env));

            if (env.id() == 0) unshuffle();

            env.barrier();
        }


        /// Invoked by the user to run the algorithm with one thread

        /// This is a collective call ... all processes in world should call
        /// this routine, though processes without data will immediately return
        /// without any synchronization.
        void solve() {
            run(A.get_world(), TaskThreadEnv(1,0,0));
        }


        /// Returns length of row
        int get_rowdim() const {return rowdim;}


        /// Returns length of column
        int get_coldim() const {return coldim;}

        /// Returns a reference to the world
        World& get_world() const {
            return A.get_world();
        }

        /// Returns rank of this process in the world
        ProcessID get_rank() const {
            return rank;
        }
    };


    template <typename T>
    class TestSystolicMatrixAlgorithm : public SystolicMatrixAlgorithm<T> {
        volatile int niter;
    public:
        TestSystolicMatrixAlgorithm(DistributedMatrix<T>& A, int tag)
            : SystolicMatrixAlgorithm<T>(A, tag)
            , niter(0)
        {
            madness::print("Testing SystolicMatrixAlgorithm ",
                           SystolicMatrixAlgorithm<T>::get_coldim(),
                           SystolicMatrixAlgorithm<T>::get_rowdim());
        }

        void kernel(int i, int j, T* rowi, T* rowj) {
            for (int k=0; k < SystolicMatrixAlgorithm<T>::get_rowdim(); ++k) {
                MADNESS_ASSERT(rowi[k] == i);
                MADNESS_ASSERT(rowj[k] == j);
            }
        }

        void start_iteration_hook(const TaskThreadEnv& env) {
            int id = env.id();
            if (id == 0) {
                ++niter;
            }
        }

        bool converged(const TaskThreadEnv& env) const {
            if (niter >= 3) {
                if (env.id() == 0) {
                    madness::print("    done!");
                }
                return true;
            }
            else {
                return false;
            }
        }
    };

}

using namespace madness;

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    madness::World world(MPI::COMM_WORLD);

    redirectio(world);

    try {
        for (int64_t n=1; n<100; ++n) {
            int64_t m = 2*n;
            DistributedMatrix<double> A = column_distributed_matrix<double>(world, n, m);
            int64_t ilo, ihi;
            A.local_colrange(ilo, ihi);
            for (int i=ilo; i<=ihi; ++i) A.data()(i-ilo,_) = i;

            world.taskq.add(new TestSystolicMatrixAlgorithm<double>(A, 3333));
            world.taskq.fence();

            for (int i=ilo; i<=ihi; ++i) {
                for (int k=0; k<m; ++k) {
                    MADNESS_ASSERT(A.data()(i-ilo,k) == i);
                }
            }
        }
    }
    catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    }
    catch (char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    }
    catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    MPI::Finalize();
}
