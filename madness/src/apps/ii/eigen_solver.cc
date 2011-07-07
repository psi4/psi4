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
/// \file testsystolic.cc
/// systolic example of eigen solver using one-sided Jacobi method.

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/world.h>
#include <utility>
#include <tensor/tensor.h>
#include <ii/systolic.h>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace madness;

template <typename T>
class SystolicEigensolver : public SystolicMatrixAlgorithm<T> {
public:
    SystolicEigensolver<T>(DistributedMatrix<T> &AV, int tag);
    virtual ~SystolicEigensolver(){}

    void start_iteration_hook(const TaskThreadEnv& env);
    void kernel(int i, int j, T* rowi, T* rowj);
    void end_iteration_hook(const TaskThreadEnv& env);
    bool converged(const TaskThreadEnv& env) const;

private:
    DistributedMatrix<T>& AV; /// concatnated two matrix A and V. V will holds eigen vector after calcuration
    World& world;
    const T tolmin; ///threshld
    int niter; /// number of iteration
    int nrot; /// number of rotation in one iteration
    int nrotsum; /// sum of rotarion for all iteration
    int size; /// size of A
    T tol; /// current threshold
    T maxd; /// maximum value of diagonal element
    T maxdaij; /// maximum value of off diagonal element

    inline T inner(const T* a, const T* b ) const{
        T s=0;
        for(int64_t i=0; i<size; i++) s += a[i] * b[i];
        return s;
    }

};
template <typename T>
SystolicEigensolver<T>::SystolicEigensolver(DistributedMatrix<T>& AV, int tag):
    SystolicMatrixAlgorithm<T>( AV, tag ), AV(AV),
    world(AV.get_world()),
    tolmin(1.0e-6),
    niter(0),
    nrot(0), nrotsum(0), size(AV.rowdim()/2),
    tol(1.0e-2), maxd(0), maxdaij(1.0e-1)
{
    MADNESS_ASSERT(AV.is_column_distributed());
    MADNESS_ASSERT(AV.coldim()*2 == AV.rowdim());
    print("One-sided Jacobi start");
}
template <typename T>
void SystolicEigensolver<T>::start_iteration_hook(const TaskThreadEnv& env) {
    if ( env.id() == 0){

        // calculate threshold using parameters from this iteration
        tol = std::min<T>( tol, std::min<T>(maxdaij*1.0e-1, maxdaij*maxdaij) );
        tol = std::max<T>( tol, tolmin );

        // clear some paremeters for a new iteration
        niter++;
        nrot = 0;
        maxdaij = 0;
    }
}
template <typename T>
void SystolicEigensolver<T>::kernel(int i, int j, T* rowi, T* rowj) {
    /// get elements of A and V from concatenated row
    T *ai = rowi;
    T *aj = rowj;
    T *vi = rowi + size;
    T *vj = rowj + size;

    T aii = inner(vi, ai);
    T ajj = inner(vj, aj);
    T aij = inner(vi, aj);
    T daij = std::abs<T>(aij);
    T s = ajj-aii;
    T ds = std::abs<T>(s);

    maxd = std::max<T>(maxd, std::max<T>( std::abs<T>(aii), std::abs<T>(ajj) ) );
    maxdaij = std::max<T>( maxdaij, daij/maxd ); // maximum value of ratio off diagonal element with diagonal element

    if( daij < ds*tol ) return; // if off diagonal elements much smaller than diagonal elements skip this step

    nrot++;

    T c,t,u;
    /// make prameters of rotation matrix
    if( ds < daij*tolmin ) c = s = sqrt(0.5); // if two diagonal elements are almost same, then rotation angle is pi/4.
    else{
        //print("trial 2"); // not good
        u = s / (2.0*aij);
        if( u > 0 ) t = 1 / (u + sqrt( u*u + 1 ));
        else t = 1 / (u - sqrt( u*u + 1 ));
        c = 1 / sqrt( t*t + 1 );
        s = c*t;
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
template <typename T>
void SystolicEigensolver<T>::end_iteration_hook(const TaskThreadEnv &env) {
    if( env.id() == 0 ) {
        world.gop.max(maxdaij);
        world.gop.max(maxd);
        world.gop.sum(nrot);
        nrotsum += nrot;
    }
}
template <typename T>
bool SystolicEigensolver<T>::converged(const TaskThreadEnv& env) const {
    int id = env.id();

    if(nrot == 0 && tol <= tolmin) {
        if (id == 0) madness::print("\tConverged! ", niter, "iteration.", size);
        return true;
    }
    else if (niter >= 5000) {
        if (id == 0) madness::print("\tDid not converged in 5000 iteration!");
        return true;
    }
    else {
        return false;
    }
}

void START_TIMER(World& world) {
}

void END_TIMER(World& world, const char* msg) {
}

double ttt, sss;
int main(int argc, char** argv) {
    initialize(argc, argv);

    madness::World world(MPI::COMM_WORLD);

    //redirectio(world); /* redirect a results to file "log.<number of processor>" */

    try {
        if( world.rank() == 0 ) {
            print("Test of Eigen solver");
            print("result: size parallel cpu_time wall_time");
            print("check: size max(off-diagonal) max(diagonal) max(V^TAV-lambdaE)");
        }
        std::srand(time(NULL));
        int64_t pow[2];
        for (int64_t n=2; n<=256; n*=2) {
            // make symmetolic matrix, then distribute it all processes
            Tensor<double> sym_tensor( n, n*2 ); // right n*n is for identity matrix, left for symmetolic matrix A
            pow[0] = (int64_t)(std::rand() % 5 - 2);
            pow[1] = pow[0] + (int64_t)(std::rand() % 3 );
            if (world.rank() == 0) {
                sym_tensor.fillrandom();
                sym_tensor -= 0.5;
                for(int i=0; i<n; i++){
                    // all diagonal elements tensor(i,i) must be begger than sum_j tensor(i,j)
                    double abs_sum=0;
                    for(int j=0; j<i; j++){
                        if (i != j) {
                            sym_tensor(i, j) = sym_tensor(j, i) *= std::pow<double>(10, pow[0]);
                            sym_tensor(i, j+n) = sym_tensor(j, i+n) = 0.0;
                        }
                        abs_sum += std::abs<double>( sym_tensor(i,j) );
                    }
                    sym_tensor(i, i+n) = 1.0;
                    if( sym_tensor(i,i) > 0 )
                        sym_tensor(i,i) = (sym_tensor(i,i) + 2.0*abs_sum) * std::pow( 10, pow[1]);
                    else
                        sym_tensor(i,i) = (sym_tensor(i,i) - 2.0*abs_sum) * std::pow( 10, pow[1]);
                }
            }
            world.gop.broadcast(sym_tensor.ptr(), sym_tensor.size(), 0);
            world.gop.fence(); // end of making symmetolic matrix
            ttt=wall_time(); sss=cpu_time();
            // broad cast matrix element then each processor copys elements
            DistributedMatrix<double> AV = column_distributed_matrix<double>(world, n, n*2);
            AV.copyin(sym_tensor);
            //if( n < 4 ) print(AV.data()); // for check

            world.taskq.add(new SystolicEigensolver<double>(AV, 3333));
            world.taskq.fence();

            // gather all data on whole porcessors
            // NOTE: since I wanted to use each colmuns as a sequence, both of the results are transposed
            Tensor<double> eigens( n, n*2 );
            AV.copyout(eigens);
            world.gop.sum(eigens.ptr(), eigens.size());
            world.gop.fence();
            ttt=wall_time()-ttt; sss=cpu_time()-sss;
            if (world.rank()==0) printf("timer: %05d %04d %8.2fs %8.2fs\n", (int)n, (int)world.size(), sss, ttt);

            Tensor<double> eigvec = transpose( eigens(_, Slice(n, -1)) );
            Tensor<double> eig_val = inner( eigens(_, Slice(0, n-1)), eigvec );

            // check the answer
            if(world.rank() == 0){
                // V^T A V - \ E : this should be zero
                Tensor<double> error = inner( transpose(eigvec), inner( sym_tensor(_, Slice(0,n-1)), eigvec )) - eig_val;

                double maxdaij=0; // maximum value of off diagonal element in error matrix. should be zero
                double maxdaii=0; // maximum value of diagonal element in error matrix. should be zero
                double max_off_diag=0; // muximum off diagonal element value of \ I
                for( int64_t i=0; i<error.dim(0); i++ ){
                    for( int64_t j=0; j<error.dim(1); j++){
                        // max ( [ abs(aij) | aij <-{A}ij, i!=j ] )
                        if( i!=j ){
                            maxdaij = std::max<double>( maxdaij, std::abs<double>( error(i,j) ));
                            max_off_diag = std::max<double>( max_off_diag, std::abs<double>( eig_val(i,j) ));
                        }
                        // max ( [ abs(aij) | aij <-{A}ij, i!=j ] )
                        else maxdaii = std::max<double>( maxdaii, std::abs<double>( error(i,j) ));
                    }
                }
                print("check:", n, maxdaij, maxdaii, max_off_diag);
            }
            print("\n");
            world.gop.fence();
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

