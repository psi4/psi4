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
#ifndef MADNESS_SYSTOLIC_H
#define MADNESS_SYSTOLIC_H
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
        int64_t coldim() const { return n; }
        
        /// Returns the row dimension of the matrix ... i.e., m for A(n,m)
        int64_t rowdim() const { return m; }
        
        /// Returns the column tile size
        int64_t coltile() const { return tilen; }
        
        /// Returns the row tile size
        int64_t rowtile() const { return tilem; }
        
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
        World& get_world() const {return _world;}
        
        /// Returns reference to data
        Tensor<T>& data()  {return t;}
        
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

    private:
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
        
    };
    
    
    /// Generates an (n,m) matrix distributed by columns (row dimension is not distributed)

    /// Quietly forces an even column tile size for ease of use in the systolic matrix algorithms
    template <typename T>
    DistributedMatrix<T> column_distributed_matrix(World& world, int64_t n, int64_t m, int64_t coltile=0) {
        /* 
           if given coltile is not enough to hold entire data of matrix or used default coltile = 0,
           recalculate it.
        */
        if (world.size()*coltile < n) coltile = (n-1)/world.size() + 1;
        if ((coltile&0x1)) coltile++;
        coltile = std::min(coltile,n);
        
        return DistributedMatrix<T>(world, n, m, coltile, m);
    }
    
    
    /// Generates an (n,m) matrix distributed by rows (column dimension is not distributed)
    template <typename T>
    DistributedMatrix<T> row_distributed_matrix(World& world, int64_t n, int64_t m, int64_t rowtile=0) {
        if (world.size()*rowtile < m) rowtile = (n-1)/world.size() + 1;
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
    
    
    /// Generates a distributed matrix with rows of \c a and \c b concatenated
    
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

        if(a.local_size() > 0) c.data()( _ , Slice(0,ma-1) ) = a.data()(___);
        if(b.local_size() > 0) c.data()( _ , Slice(ma, -1) ) = b.data()(___);

        return c;
    }

    template <typename T>
    DistributedMatrix<T> concatenate_rows( const DistributedMatrix<T>& a, const DistributedMatrix<T>& b, const DistributedMatrix<T>& c, const DistributedMatrix<T>& d) {
        MADNESS_ASSERT(a.coldim()==b.coldim() && b.coldim()==c.coldim() && c.coldim()==d.coldim());
        MADNESS_ASSERT(a.coltile()==b.coltile() && b.coltile()==c.coltile() && c.coltile()==d.coltile());
        MADNESS_ASSERT(a.is_column_distributed() && b.is_column_distributed() && c.is_column_distributed() && d.is_column_distributed());
        
        int64_t ma = a.rowdim();
        int64_t mb = b.rowdim();
        int64_t mc = c.rowdim();
        int64_t md = d.rowdim();
        
        DistributedMatrix<T> result(a.get_world(), a.coldim(), ma+mb+mc+md, a.coltile(), ma+mb+mc+md);

        if(a.local_size() > 0) result.data()( _ , Slice(0,ma-1) ) = a.data()(___);
        if(b.local_size() > 0) result.data()( _ , Slice(ma, ma+mb-1) ) = b.data()(___);
        if(c.local_size() > 0) result.data()( _ , Slice(ma+mb, ma+mb+mc-1) ) = c.data()(___);
        if(d.local_size() > 0) result.data()( _ , Slice(ma+mb+mc, -1) ) = d.data()(___);

        return result;
    }


    template <typename T>
    DistributedMatrix<T> concatenate_columns(const DistributedMatrix<T>& a, const DistributedMatrix<T>& b) {
        MADNESS_ASSERT(a.rowdim()==b.rowdim() && a.rowtile()==b.rowtile() && a.is_row_distributed() && b.is_row_distributed());

        int64_t ma = a.coldim();
        int64_t mt = ma + b.coldim();

        DistributedMatrix<T> c(a.get_world(), mt, a.rowdim(), b.rowtile(), mt);

        if(a.local_size() > 0) c.data()( Slice(0,ma-1), _ ) = a.data()(___);
        if(a.local_size() > 0) c.data()( Slice(ma,-1), _ ) = b.data()(___);

        return c;
    }
    
    /// make identity matrix with same propaties of A
    template <typename T>
    DistributedMatrix<T> idMatrix(const DistributedMatrix<T>& A){
	    int64_t n = A.coldim();
	    int64_t m = A.rowdim();
	    MADNESS_ASSERT( n == m );

	    DistributedMatrix<T> result( A.get_world(), n, m, A.coltile(), A.rowtile() );

	    int64_t ilo, ihi;
	    result.local_colrange(ilo,ihi);
	    for(int64_t i=ilo; i<=ihi; ++i) {
		    result.data()(i-ilo, i) = 1;
	    }

	    return result;
    }

    /// Base class for parallel algorithms that employ a systolic loop to generate all row pairs in parallel
    template <typename T>
    class SystolicMatrixAlgorithm : public TaskInterface {
    private:
        DistributedMatrix<T>& A;        //< internal data
        int64_t nproc;                  //< No. of processes with rows of the matrix (not size of world)
        int64_t coldim;                 //< A(coldim,rowdim)
        int64_t rowdim;                 //< A(coldim,rowdim)
        int64_t nlocal;                 //< No. of local pairs
        const ProcessID rank;           //< Rank of current process
        const int tag;                  //< MPI tag to be used for messages
        std::vector<T*> iptr, jptr;     //< Indirection for implementing cyclic buffer !! SHOULD BE VOLATILE ?????
        std::vector<int64_t> map;       //< Used to keep track of actual row indices

        void iteration(const TaskThreadEnv& env) {
            env.barrier();
            start_iteration_hook(env);
            env.barrier();

            int64_t ilo, ihi;
            A.local_colrange(ilo, ihi);
            
            int neven = coldim + (coldim&0x1);

            int pairlo = rank*A.coltile()/2;

            int threadid = env.id();
            int nthread = env.nthread();

            for (int loop=0; loop<(neven-1); loop++) {

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
            
            end_iteration_hook(env);
            env.barrier();
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
            //Tensor<T> tmp(2L, t.ndim(), false);
            Tensor<T> tmp(2L, t.dims(), false);
            T* tp = tmp.ptr();
            for (int64_t i=0; i<nlocal; i++) {
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
            for (int64_t i=0; i<nlocal-1; i++) {
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
            
            for (int64_t i=0; i<nlocal; i++) {
                iptr[i] = &t(i,0);
                jptr[i] = &t(i+nlocal,0);
            }

            // If no. of rows is odd, last process should have an empty last row
            if (rank==(nproc-1) && (coldim&0x1)) jptr[nlocal-1] = 0;

            // Initialize map from logical index order to actual index order

            int neven = (coldim+1)/2;
            int ii=0;
            for (ProcessID p=0; p<nproc; p++) {
                int64_t lo, hi;
                A.get_colrange(p, lo, hi);
                int p_nlocal = (hi - lo + 2)/2;
                //print("I think process",p,"has",lo,hi,p_nlocal);
                for (int i=0; i<p_nlocal; i++) {
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


        /// Invoked by all threads at the end of each iteration

        /// There is a thread barrier before and after the invocation of this routine
        virtual void end_iteration_hook(const TaskThreadEnv& env) {}


        /// Invoked by the task queue to run the algorithm with multiple threads
        void run(World& world, const TaskThreadEnv& env) {
            //if (nlocal <= 0) return; // Nothing to do
            
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
        void solve() { run(A.get_world(), TaskThreadEnv(1,0,0)); }


        /// Returns length of row
        int get_rowdim() const {return rowdim;}

        /// Returns length of column
        int get_coldim() const {return coldim;}

        /// Returns a reference to the world
        World& get_world() const { return A.get_world(); }

        /// Returns rank of this process in the world
        ProcessID get_rank() const{ return rank; }

    };

    template <typename T>
    class LocalizeBoys : public SystolicMatrixAlgorithm<T>
    {
    public:
        LocalizeBoys<T>( DistributedMatrix<T> &M, const std::vector<int>& set, long nmo, int tag,
                const double thresh = 1e-6, const double thetamax = 0.5, const bool randomized = true, const bool doprint = false):
            SystolicMatrixAlgorithm<T>(M, tag),
            M(M), // distributed marix of ( U, X, Y, Z)
            world(M.get_world()),
            set(set), // set of orbital
            nmo(nmo), // number of molecule
            niter(0), // number of iteration
            ndone(0), // number of rotation
            nrot(0), 
            thetamax(thetamax),
            thresh(thresh),
            tol(thetamax),
            randomized(false),
            doprint(doprint)
        {
            if (doprint) madness::print("Start boys localization\n");
        }

        virtual ~LocalizeBoys() {}   

        void start_iteration_hook(const TaskThreadEnv& env);
        void kernel(int i, int j, T* rowi, T* rowj);
        void end_iteration_hook(const TaskThreadEnv& env);
        bool converged(const TaskThreadEnv& env) const;

    private:
        DistributedMatrix<T> M;
        World& world;
        std::vector<int> set;
        long nmo, niter, ndone, ndone_iter; // number of molecule, number of iteration, number of rotation, number of rotation in 1 iteration
        volatile int64_t nrot; // number of rotation in 1 iteration
        const double thetamax, thresh;
        double tol; // current error
        volatile double maxtheta; // max rotation angle
        const bool randomized, doprint;

        void drot(T restrict a[], T restrict b[], double sin, double cos); 
        inline T DIP(const T ij[], const T kl[]) const;
        inline T inner(const T a[], const T b[] ) const;
    };
    template <typename T>
    void LocalizeBoys<T>::start_iteration_hook(const TaskThreadEnv& env)
    {
        int id = env.id();
        if(doprint) {
            T sum = 0.0;
            int64_t ilo, ihi;
            M.local_colrange(ilo, ihi);
            for(int64_t i=0; i <=(ihi-ilo); i++){
                T ii[] = { M.data()(i,i+ilo+nmo), M.data()(i,i+ilo+2*nmo), M.data()(i,i+ilo+3*nmo) };
                sum += DIP(ii, ii);
            }
            env.barrier();
            if (env.id() == 0) world.gop.sum(sum);
            env.barrier();

            // print a result of previous iteration
            printf("\titeration %ld sum=%.4f ndone=%ld tol=%.2e, maxtheta=%.2e\n", niter, sum, ndone, tol, maxtheta);
        }

        if( id == 0 ){
            this->nrot = 0; /// number of rotation in this iteration
            this->maxtheta = 0.0; /// maximum rotation angle in this iteration
        }
        env.barrier();
    }
    template <typename T>
    void LocalizeBoys<T>::kernel(int i, int j, T rowi[], T rowj[])
    {
        if(set[i] != set[j]) return;

        // make rowi and rowj since we're using one-sided jacobi
        T *ui = rowi;
        T *uj = rowj;
        T *ri[] = { rowi + nmo, rowi + 2*nmo, rowi + 3*nmo };
        T *rj[] = { rowj + nmo, rowj + 2*nmo, rowj + 3*nmo };
        T ii[] = { inner(ui, ri[0]), inner(ui, ri[1]), inner(ui, ri[2]) };
        T ij[] = { inner(ui, rj[0]), inner(ui, rj[1]), inner(ui, rj[2]) };
        T jj[] = { inner(uj, rj[0]), inner(uj, rj[1]), inner(uj, rj[2]) };
        bool doit = false;

        double g = DIP(ij, jj) - DIP(ij, ii);
        double h = 4.0*DIP(ij, ij) + 2.0*DIP(ii, jj) - DIP(ii, ii) - DIP(jj, jj);

        if (h >= 0.0) {
            if (doprint) print("\t\tforcing negative h", i, j, h);
            doit = true;
            h = -1.0;
        }
        double theta = -g / h;

        this->maxtheta = std::max<double>(fabs(theta), (double)maxtheta); // this doesn't seem to thread safe

        /// restriction
        if (fabs(theta) > thetamax){
            if (doprint) print("\t\trestricting", i,j);

            if (g < 0) theta = -thetamax;
            else theta = thetamax * 0.8;
            doit = true;
        }

        // randomize will be implemented here
        // double sij = DIP(ij, ij); // this line will be used by randomize

        // rotation
        if (fabs(theta) >= tol || randomized || doit){
            if (doprint) print("\t\trotating", i,j, theta);
            this->nrot++;

            double c = cos(theta);
            double s = sin(theta);
            drot(ri[0], rj[0], s, c); 
            drot(ri[1], rj[1], s, c);
            drot(ri[2], rj[2], s, c);
            drot(ui, uj, s, c);
        }
    }
    template <typename T>
    void LocalizeBoys<T>::end_iteration_hook(const TaskThreadEnv& env)
    {
        if (env.id() == 0) {
            this->niter++;
            this->tol = std::max(0.1 * maxtheta, thresh);
            this->ndone_iter = nrot;
            world.gop.max(tol);

            world.gop.sum(ndone_iter); // get total number of rotation whole processes
            this->ndone += ndone_iter;
        }
        env.barrier();
    }
    template <typename T>
    bool LocalizeBoys<T>::converged(const TaskThreadEnv& env) const
    {
        if( ndone_iter == 0 && tol <= thresh){
            if (doprint) madness::print("\tBoys localization converged in", ndone, "steps.");
            return true;
        }
        else if(niter >= 300){
            if(doprint) madness::print("\tWARNING!! Boys localization did not fully converged in", niter, "iteration!\n");
            return true;
        }
        else return false;
        env.barrier();
    }
    /// rotate matrix using sin and cos
    template <typename T>
    void LocalizeBoys<T>::drot(T restrict a[], T restrict b[], double sin, double cos)
    {
        for ( long i=0; i<this->nmo; i++ ) {
            T aa = a[i]*cos - b[i]*sin;
            T bb = a[i]*sin + b[i]*cos;
            a[i] = aa;
            b[i] = bb;
        }
    }
    template <typename T>
    inline T LocalizeBoys<T>::DIP(const T ij[], const T kl[]) const
    {
        return ij[0] * kl[0] + ij[1] * kl[1] + ij[2] * kl[2];
    }
    template <typename T>
    inline T LocalizeBoys<T>::inner(const T a[], const T b[] ) const
    {
        T s=0;
        for(int64_t i=0; i<nmo; i++){
            s += a[i] * b[i];
        }
        return s;
    }
}    
#endif
