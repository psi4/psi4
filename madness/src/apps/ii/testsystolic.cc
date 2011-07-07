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
/* \file testsystolic2.cc
 * systolic example of eigen solver using one-sided Jacobi method.
 */
#include <world/world.h>
#include <utility>
#include <tensor/tensor.h>
#include <ii/systolic.h>
#include <ctime>

using namespace madness;

template <typename T>
class SystolicEigensolver : public SystolicMatrixAlgorithm<T> {
public:
    SystolicEigensolver<T>(DistributedMatrix<T>& AV, int tag);

    void kernel(int i, int j, T* rowi, T* rowj) {
        /// get elements of A and V from concatenated row
        T *ai = rowi;
        T *aj = rowj;
        T *vi = rowi + size;
        T *vj = rowj + size;

        T aii = inner(vi, ai);
        T ajj = inner(vj, aj);
        T aij = inner(vi, aj);
        T daij = fabs(aij);

        maxd = std::max<T>( std::max<T>( fabs(aii), fabs(ajj) ), maxd );
        maxdaij = std::max<T>( maxdaij, daij/maxd );

        if( daij < maxd*tol ) return;

        T s = ajj-aii;
        T ds = fabs(s);

        if( daij > ds*tolmin ){
            nrot++;
            T c,t,u;
            /// make prameters of rotation matrix
            if( tolmin*daij > ds ) c = s = 1/sqrt(2.0);
            else{
                t = aij/s;
                u = 0.25/sqrt(0.25+t*t);
                c = sqrt(0.5+u);
                s = 2.0*t*u/c;
            }

            /// update all elements
            for (int k=0; k < size; k++) {
                t = ai[k];
                u = aj[k];
                ai[k] = c*t - s*u;
                aj[k] = c*u + s*t;

                t = vi[k];
                u = vj[k];
                vi[k] = c*t - s*u;
                vj[k] = c*u + s*t;
            }
        }
    }

    void start_iteration_hook(const TaskThreadEnv& env) {
        int id = env.id();

        if (id == 0) world.gop.max(maxdaij);
        env.barrier();

        if (id == 0){
            tol = std::min( tol, std::min(maxdaij*0.1, maxdaij*maxdaij) ); 
            tol = std::max( tol, 5.0e-16 );
            niter++;

            maxdaij = 0;
            nrot = 0; // don't move this line to above
        }
    }

    void end_iteration_hook(const TaskThreadEnv& env) {

        int id = env.id();

        if (id == 0) world.gop.sum(nrot);
        env.barrier();

        nrotsum += nrot;
    }

    bool converged(const TaskThreadEnv& env) const; 

    Tensor<T> get_eval() const;

    DistributedMatrix<T> get_evec() const;

private:
    /** constant members */
    static const T tolmin = (T)5.0e-16; ///threshld

    DistributedMatrix<T>& AV; /// concatnated two matrix A and V. V will holds eigen vector after calcuration
    World& world;
    volatile int niter;
    int nrot, nrotsum, size; 
    T tol, maxd, maxdaij;

    inline T inner(const T* a, const T* b ) const{
        T s=0;
        for(int64_t i=0; i<size; i++){
            s += a[i] * b[i];
        }
        return s;
    }

};
template <typename T>
SystolicEigensolver<T>::SystolicEigensolver(DistributedMatrix<T>& AV, int tag):
    SystolicMatrixAlgorithm<T>( AV, tag ), AV(AV),
    world(AV.get_world()),
    niter(0),
    nrot(0), nrotsum(0), size(AV.rowdim()/2), 
    tol((T)1e-2), maxd(0), maxdaij(1e3) // just I want a very big value
{
    MADNESS_ASSERT(AV.is_column_distributed());
    MADNESS_ASSERT(AV.coldim()*2 == AV.rowdim());
    print("One-sided Jacobi start");
}
template <typename T>
DistributedMatrix<T> SystolicEigensolver<T>::get_evec() const{
    int64_t m = AV.local_coldim();
    int64_t n = AV.local_rowdim()/2;
    DistributedMatrix<T> result = column_distributed_matrix<T>(world, m, n);
    result.data() = AV.data()(_,Slice(size,-1));

    return result;
}
template <typename T>
Tensor<T> SystolicEigensolver<T>::get_eval() const{
    long int lsize = AV.local_coldim();
    Tensor<T> result(lsize);

    for(int64_t i=0; i<lsize; i++){
        Tensor<T> ai= AV.data()(i, Slice(0, size-1)); 
        Tensor<T> vi= AV.data()(i, Slice(size, -1));
        result[i] = madness::inner(vi, ai)(0,0);
    }

    return result;
}
template <typename T>
bool SystolicEigensolver<T>::converged(const TaskThreadEnv& env) const {
    int id = env.id();

    if(nrot == 0 && tol <= tolmin){
        if (id == 0) {
            madness::print("    Converged! ", size, "\n");
        }
        return true;
    } else if (niter >= 50) {
        if (id == 0) {
            madness::print("    Did not converged in 50 iteration!", "\n");
        }
        return true;
    } else
        return false; 
}
/* trial function.
template <typename T>
SystolicEigensolver<T> void systolic_eigensolver(DistributedMatrix<T>& A, int tag )
{
    MADNESS_ASSERT(A.is_column_distributed() == true);
    /// initialize V as identity matrix of size(A)
    DistributedMatrix<T> V = column_distributed_matrix<T>( A.get_world(), A.coldim(), A.rowdim() );

    int64_t ilo, ihi;
    V.local_colrange(ilo, ihi);
        for(int i=ilo; i<=ihi; i++){
        V.data()(i-ilo,i) = 1.0;
    }

    DistributedMatrix<T> A_V = concatenate_rows(A,V);

    A.get_world().taskq.add(new SystolicEigensolver<T>(A_V, tag));
    A.get_world().taskq.fence();

} */

int main(int argc, char** argv) {
    initialize(argc, argv);

    madness::World world(MPI::COMM_WORLD);

    redirectio(world); /* redirect a results to file "log.<number of processor>" */

    try {
        print("Test of testsystolic2.cc\n");
        print("result: size time eig_val");
        for (int64_t n=10; n>1; n-=2) {
            DistributedMatrix<double> A = column_distributed_matrix<double>(world, n, n);

            int64_t ilo, ihi, jlo, jhi;
            A.local_rowrange(ilo, ihi); // local row range is equal to global row range
            A.local_colrange(jlo, jhi); /* get column range of A */
            for (int j=jlo; j<=jhi; j++) {
                for (int i=ilo; i<=ihi; i++) {
                    A.data()(j-jlo, i-ilo) = (i + j) * 0.1 ; //in this way, matrix is symmetric
                }
                A.data()(j-jlo,j) = ( A.data()(j-jlo,j)+0.5 ) * n;  
            }

            DistributedMatrix<double>  V = idMatrix(A);
            DistributedMatrix<double> AV = concatenate_rows(A, V);
            SystolicEigensolver<double> sA(AV, 3334);

            double t = cpu_time();
            sA.solve();
            print("result:", n, cpu_time()-t, sA.get_eval());

            if(world.size() == 1){
                /* check the answer*/
                Tensor<double> eigvec = sA.get_evec().data();
                print("eigen vector\n", eigvec);

                /* U^T * U = I */
                print("U^T * U\n", inner(transpose(eigvec), eigvec));
                //print("U^T * U\n", mxm2(transpose(eigvec), eigvec)); // must be identity matrix

                /* A * U = lambda * U */
                print("eval\n", sA.get_eval());
                print("A * U\n", inner(A.data(), transpose(eigvec)));
                //print("A * U\n", mxm2(A.data(), transpose(eigvec)));

                ///test for mTxm ... O.K. 
                /*
                   Tensor<double> t(2,2);
                   t(0,0) = 1; t(0,1) = 2; t(1,0) = 3; t(1,1) = 4;
                   print(t);
                   print( t(1,Slice()), t(Slice(),0)); /// result ... (3,4) , (1,3)
                /// result must be ( 10 , 14, 14, 20)
                print(mxm2(transpose(t), t)); //... O.K.
                print(transpose(t).emul(t)); //... N.G. emul makes multiply of each element
                 */
            }

            /*
               print("matrix A");
               print(A.data());

               DistributedMatrix<double> B = column_distributed_matrix<double>(world, n, l);
               B.local_rowrange(ilo, ihi);
               B.local_colrange(jlo, jhi); // get column range of B
               for (int i=ilo; i<=ihi; i++) {
               for (int j=jlo; j<=jhi; j++) {
               B.data()(j-jlo,i-ilo) = j + i*100 ;
               }
               }
               print("matrix B");
               print(B.data());

               DistributedMatrix<double> D = concatenate_rows(A, B);

               print("matrix D");
               print(D.data()); //... O.K

               world.taskq.add(new TestSystolicMatrixAlgorithm<double>(C, 3333));
               world.taskq.fence();

               world.taskq.add(new SystolicEigensolver<double>(A, 3334));
               world.taskq.fence();
             */
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
    catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (char* s) {
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

    finalize();
}

