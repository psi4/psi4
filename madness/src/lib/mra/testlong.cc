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


  $Id: tests-hqi.cc 257 2007-06-25 19:09:38Z HartmanBaker $
*/

/// \file testlong.cc
/// \brief The loadbal test suite


#include <mra/mra.h>

const double PI = 3.1415926535897932384;

using namespace madness;

template <typename T, std::size_t NDIM>
class GaussianFunctor : public FunctionFunctorInterface<T,NDIM> {
private:
    typedef Vector<double,NDIM> coordT;
    const std::vector<coordT> center;
    const std::vector<double> exponent;
    const std::vector<T> coefficient;

public:
    GaussianFunctor(const std::vector<coordT>& center, std::vector<double> exponent, std::vector<T> coefficient)
            : center(center), exponent(exponent), coefficient(coefficient) {};

    GaussianFunctor(const coordT& center, double exponent, T coefficient)
            : center(std::vector<coordT>(1,center)), exponent(std::vector<double>(1,exponent)), coefficient(std::vector<T>(1,coefficient)) {};

    T operator()(const coordT& x) const {
        T retval = 0;
        for (unsigned int j=0; j<center.size(); ++j) {
            double sum = 0.0;
            for (int i=0; i<NDIM; ++i) {
                double xx = center[j][i]-x[i];
                sum += xx*xx;
            }
            retval += coefficient[j]*exp(-exponent[j]*sum);
        }
        return retval;
    };

    GaussianFunctor operator+ (const GaussianFunctor& other) const {
        std::vector<coordT> newcenter = this->center;
        std::vector<double> newexponent = this->exponent;
        std::vector<T> newcoefficient = this->coefficient;
        newcenter.insert(newcenter.end(), other.center.begin(), other.center.end());
        newexponent.insert(newexponent.end(), other.exponent.begin(), other.exponent.end());
        newcoefficient.insert(newcoefficient.end(), other.coefficient.begin(), other.coefficient.end());
        return (GaussianFunctor(newcenter, newexponent, newcoefficient));
    };

    GaussianFunctor operator- (const GaussianFunctor& other) const {
        std::vector<T> newcoefficient;
        //	for (std::vector<T>::iterator it = other.coefficient.begin(); it != other.coefficient.end(); ++it) {
        //	  newcoefficient.push_back(-1*(*it));
        //	}
        int size = other.coefficient.size();
        for (int i = 0; i < size; ++i) {
            newcoefficient.push_back(-1*other.coefficient[i]);
            return (this + GaussianFunctor(other.center, other.exponent, newcoefficient));
        }
    };

    void print() const {
        madness::print("Sum of", center.size(), "gaussians:");
        for (unsigned int i = 0; i < center.size(); ++i) {
            madness::print("   g[", i, "] : =", coefficient[i], "* exp(", -exponent[i], "(", center[i], "- x )^2 )");
        }
    };
};


/// Returns a new functor combining two functors via operation op(left,right)

template <typename resultT, typename L, typename R, typename opT, std::size_t NDIM>
class BinaryOp : public FunctionFunctorInterface<resultT,NDIM> {
    typedef Vector<double,NDIM> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<L,NDIM> > functorL;
    typedef std::shared_ptr< FunctionFunctorInterface<R,NDIM> > functorR;

    functorL left;
    functorR right;
    opT op;

public:
    BinaryOp(functorL& left, functorR& right, opT& op)
            : left(left), right(right), op(op) {};

    resultT operator()(const coordT& x) const {
        return op((*left)(x),(*right)(x));
    };
};

double ttt, sss;
#define START_TIMER world.gop.fence(); ttt=wall_time(); sss=cpu_time()
#define END_TIMER(msg) ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s   %.6e   %.6e\n", msg, sss, ttt)


template <typename T, std::size_t NDIM>
Cost a_cost_function(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) {
    return 1;
}

template <typename T, std::size_t NDIM>
void test_loadbal(World& world) {
    typedef Vector<double,NDIM> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0) print("at beginning of test_loadbal");

    for (int i=0; i<NDIM; ++i) {
        FunctionDefaults<NDIM>::cell(i,0) = -10.0;
        FunctionDefaults<NDIM>::cell(i,1) =  10.0;
    }

    int nspikes = 2;
    //int nspikes = 10;
    std::vector<coordT> vorigin(nspikes);
    const double expnt = 64.0;
    const double expnt1 = 4096;
    std::vector<double> vexpnt(nspikes);
    Vector<double, NDIM> dcell, avgcell;
    for (int i = 0; i < NDIM; ++i) {
        dcell[i] = FunctionDefaults<NDIM>::cell(i,1) - FunctionDefaults<NDIM>::cell(i,0);
        avgcell[i] = (FunctionDefaults<NDIM>::cell(i,0) + FunctionDefaults<NDIM>::cell(i,1))/2;
    }
    for (int i = 0; i < nspikes; ++i) {
        Vector<double, NDIM> v(0);
        for (int j = 0; j < NDIM; ++j) {
            v[j] = 0.2;
            //v[j] = ((double) rand()/(double) RAND_MAX)*dcell[j] - FunctionDefaults<NDIM>::cell(j,1);
        }
        if (i%2) {
            vexpnt[i] = expnt1;
        }
        else {
            vexpnt[i] = expnt;
        }
        vorigin[i] = v;
    }
    const double coeff = pow(2.0/PI,0.25*NDIM);
    std::vector<double> vcoeff(nspikes,coeff);
    if (world.rank() == 0) print("about to make gf");
    GaussianFunctor<T,NDIM> *gf = new GaussianFunctor<T,NDIM>(vorigin, vexpnt, vcoeff);
    if (world.rank() == 0) print("about to make functor");
    functorT functor(gf);
    if (world.rank() == 0) {
        gf->print();
    }

    for (int k=5; k >=5; k-=2) {
        FunctionDefaults<NDIM>::pmap.reset(new MyPmap<NDIM>(world));

        int n = 2;
        double thresh = 1e-4;
//	double thresh = 1e-12;
//	double thresh = 5e-4;
        if (world.rank() == 0) {
            madness::print("for f with k =", k, "n =", n, "and thresh =",
                           thresh, "running on", world.nproc(), "processors:");
        }
        START_TIMER;
        Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor).thresh(thresh).initial_level(n).k(k);
        END_TIMER("H: project");
        double err2 = f.err(*functor);
        std::size_t size = f.size();
        std::size_t tree_size = f.tree_size();
        std::size_t maxsize = f.max_nodes();
        std::size_t minsize = f.min_nodes();
        std::size_t maxdepth = f.max_depth();
        START_TIMER;
        if (world.rank() == 0) {
            printf("   n=%d err=%.2e #coeff=%.2e tree_size=%.2e max_depth=%.2e max_nodes=%.2e min_nodes=%.2e log(err)/(n*k)=%.2e\n",
                   n, err2, double(size), double(tree_size), double(maxdepth), double(maxsize), double(minsize), abs(log(err2)/n/k));
        }
        world.gop.fence();
        END_TIMER("H: diagnostics");
        START_TIMER;
        test_ops(f);
        END_TIMER("H: test_ops");
        //	if (world.rank() == 0) print("ABOUT TO CREATE LOADBAL");
        START_TIMER;
        LoadBalImpl<NDIM> lb(f,a_cost_function<T,NDIM>);
        //	if (world.rank() == 0) print("ABOUT TO FENCE");
        world.gop.fence();
        END_TIMER("create LB");
        //	if (world.rank() == 0) print("ABOUT TO DO LOADBAL");
        START_TIMER;
        FunctionDefaults<NDIM>::pmap = lb.load_balance();
        //	if (world.rank() == 0) print("ABOUT TO FENCE");
        world.gop.fence();
        END_TIMER("do LB");
        START_TIMER;
        Function<T,NDIM> g = copy(f, FunctionDefaults<NDIM>::pmap);
        world.gop.fence();
        END_TIMER("copy H->LB");
        START_TIMER;
        Function<T,NDIM> h = FunctionFactory<T,NDIM>(world).functor(functor).thresh(thresh).initial_level(n).k(k);
        END_TIMER("new LB fn");
        START_TIMER;
        double err21 = h.err(*functor);
        std::size_t size1 = h.size();
        std::size_t tree_size1 = h.tree_size();
        std::size_t maxsize1 = h.max_nodes();
        std::size_t minsize1 = h.min_nodes();
        std::size_t maxdepth1 = h.max_depth();
        if (world.rank() == 0) {
            printf("   n=%d err=%.2e #coeff=%.2e tree_size=%.2e max_depth=%.2e max_nodes=%.2e min_nodes=%.2e log(err)/(n*k)=%.2e\n",
                   n, err21, double(size1), double(tree_size1), double(maxdepth1), double(maxsize1), double(minsize1), abs(log(err21)/n/k));
        }
        world.gop.fence();
        END_TIMER("LB: diagnostics");
        START_TIMER;
        test_ops(h);
        world.gop.fence();
        END_TIMER("LB: test_ops");
    }

    if (world.rank() == 0) print("test loadbal OK\n\n");
}


/*
void print_results(Vector<double,20> t) {
  double compress, reconstruct, r_copy, r_mult, r_add, r_mineq;
  double compress2, c_copy, c_mult, c_add, c_mineq, c_gaxpy;
  double reconstruct2, compress3, comp_avg, recon_avg, total_time;

  int counter = 0, counter1 = 1;

  compress = t[counter1++] - t[counter++];
  reconstruct = t[counter1++] - t[counter++];
  r_copy = t[counter1++] - t[counter++];
  r_mult = t[counter1++] - t[counter++];
  r_add = t[counter1++] - t[counter++];
  r_mineq = t[counter1++] - t[counter++];
  compress2 = t[counter1++] - t[counter++];
  c_copy = t[counter1++] - t[counter++];
  c_mult = t[counter1++] - t[counter++];
  c_add = t[counter1++] - t[counter++];
  c_mineq = t[counter1++] - t[counter++];
  c_gaxpy = t[counter1++] - t[counter++];
  reconstruct2 = t[counter1++] - t[counter++];
  compress3 = t[counter1] - t[counter];

  comp_avg = (compress + compress2 + compress3)/3.0;
  recon_avg = (reconstruct + reconstruct2)/2.0;
  total_time = t[counter1] - t[0];

  madness::print("OPERATION      | RECONSTRUCTED  | COMPRESSED ");
  madness::print("---------------+----------------+----------------");
  madness::print("Compress (1)   |", compress,   "| ---");
  madness::print("Compress (2)   |", compress2,  "| ---");
  madness::print("Compress (3)   |", compress3,  "| ---");
  madness::print("Reconstruct(1) | ---            |", reconstruct);
  madness::print("Reconstruct(2) | ---            |", reconstruct2);
  madness::print("Copy           |", r_copy, "|", c_copy);
  madness::print("Multiply       |", r_mult, "|", c_mult);
  madness::print("Min-eq (-=)    |", r_mineq, "|", c_mineq);
  madness::print("Gaxpy          | ---            |", c_gaxpy);
  madness::print("Average compress time    =", comp_avg);
  madness::print("Average reconstruct time =", recon_avg);
  madness::print("Total time for test      =", total_time);
}



template <typename T, std::size_t NDIM>
void test_ops(Function<T,NDIM>& f) {
  Vector<double, 20> t(0);
  int counter = 0;
  t[counter++] = MPI::Wtime();
  f.compress();
  t[counter++] = MPI::Wtime();
  Function<T,NDIM> c0 = f; // c = compressed
  t[counter++] = MPI::Wtime();
  f.reconstruct();
  t[counter++] = MPI::Wtime();
  Function<T,NDIM> r0 = f; // r = reconstructed
  t[counter++] = MPI::Wtime();
  Function<T,NDIM> r1 = f*r0;
  t[counter++] = MPI::Wtime();
  Function<T,NDIM> r2 = 2.0*r1;
  t[counter++] = MPI::Wtime();
  //Function<T,NDIM> r3 = r1 + PI;
  t[counter++] = MPI::Wtime();
  //r2 -= r3;
  //r2 -= r1;
  t[counter++] = MPI::Wtime();
  // Now work on compressed
  f.compress();
  t[counter++] = MPI::Wtime();
  Function<T,NDIM> c1 = f*c0;
  t[counter++] = MPI::Wtime();
  Function<T,NDIM> c2 = 2.0*c1;
  t[counter++] = MPI::Wtime();
  //Function<T,NDIM> c3 = c1 + PI;
  t[counter++] = MPI::Wtime();
  //c2 -= c3;
  //c2 -= c1;
  c2.compress();
  f.compress();
  t[counter++] = MPI::Wtime();
  // the big one:
  //c3.gaxpy(-3.46, f, 73.216);
  c2.gaxpy(-3.46, f, 73.216);
  t[counter++] = MPI::Wtime();
  //c3.reconstruct();
  c2.reconstruct();
  t[counter++] = MPI::Wtime();
  //c3.compress();
  c2.compress();
  t[counter++] = MPI::Wtime();
  if (f.world().rank() == 0) {
    print_results(t);
  }
}
*/

template <typename T, std::size_t NDIM>
void test_ops(Function<T,NDIM>& f) {
    World& world=f.world();
    START_TIMER;
    f.compress();
    END_TIMER("compress");
    START_TIMER;
    Function<T,NDIM> c0 = f; // copy
    END_TIMER("copy (C)");
    START_TIMER;
    f.truncate();
    END_TIMER("truncate (C)");
    START_TIMER;
    f.reconstruct();
    END_TIMER("reconstruct");
    for (int axis=0; axis < NDIM; ++axis) {
        START_TIMER;
        Function<T,NDIM> dfdx = diff(f,axis);
        END_TIMER("differentiate (R)");
    }
    START_TIMER;
    Function<T,NDIM> c1 = c0;
    END_TIMER("copy (C)");
    START_TIMER;
    c1.reconstruct();
    END_TIMER("reconstruct");
    START_TIMER;
    Function<T,NDIM> r0 = c1;
    END_TIMER("copy (R)");
    START_TIMER;
    Function<T,NDIM> r1 = r0*c1;
    END_TIMER("multiply (R)");
    START_TIMER;
    c1.compress();
    END_TIMER("compress");
    START_TIMER;
    Function<T,NDIM> c2 = c1*c0;
    END_TIMER("multiply");
    START_TIMER;
    c2.compress();
    END_TIMER("compress");
    START_TIMER;
    c1.compress();
    END_TIMER("compress");
    START_TIMER;
    c2.gaxpy(-3.46, c1, 73.216);
    END_TIMER("gaxpy (C)");
    START_TIMER;
    c2.reconstruct();
    END_TIMER("reconstruct");
}



int main(int argc, char**argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);

    try {
        if (world.rank() == 0) print("before startup");
        startup(world,argc,argv);
        if (world.rank() == 0) print("before test_loadbal");
        test_loadbal<double,3>(world);
    }
    catch (const MPI::Exception& e) {
        print(e);
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
        error("caught a string exception");
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

    //    print("entering final fence");
    world.gop.fence();
    //    print("done with final fence");
    MPI::Finalize();

    return 0;
}
