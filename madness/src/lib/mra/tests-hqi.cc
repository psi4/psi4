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

/// \file tests-hqi.cc
/// \brief The loadbal test suite


#include <mra/mra.h>

const double PI = 3.1415926535897932384;

using namespace madness;

template <typename T, std::size_t NDIM>
struct lbcost {
    double leaf_value;
    double parent_value;
    lbcost(double leaf_value=1.0, double parent_value=1.0) : leaf_value(leaf_value), parent_value(parent_value) {}

    double operator()(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) const {
        if (key.level() <=1) return 1000.0;
        else if (node.is_leaf()) {
            return leaf_value;
        }
        else {
            return parent_value;
        }
    }
};

template <typename T, int NDIM>
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

    void print() const {
        madness::print("Sum of", center.size(), "gaussians:");
        for (unsigned int i = 0; i < center.size(); ++i) {
            madness::print("   g[", i, "] : =", coefficient[i], "* exp(", -exponent[i], "(", center[i], "- x )^2 )");
        }
    };
};


template <typename T, int NDIM>
void test_loadbal(World& world) {
    typedef Vector<double,NDIM> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0) print("at beginning of test_loadbal");

    //    for (int i=0; i<NDIM; ++i) {
    //        FunctionDefaults<NDIM>::cell(i,0) = -10.0;
    //        FunctionDefaults<NDIM>::cell(i,1) =  10.0;
    //    }

    FunctionDefaults<NDIM>::set_cubic_cell(-10.0, 10.0);

    int nspikes = 2;
    //int nspikes = 10;
    std::vector<coordT> vorigin(nspikes);
    const double expnt = 64.0;
    const double expnt1 = 4096;
    std::vector<double> vexpnt(nspikes);
    Vector<double, NDIM> dcell, avgcell;
    for (int i = 0; i < NDIM; ++i) {
        double cell0 = FunctionDefaults<NDIM>::get_cell()(i,0);
        double cell1 = FunctionDefaults<NDIM>::get_cell()(i,1);
        dcell[i] = cell0 - cell1;
        avgcell[i] = (cell0 + cell1)/2;
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

    for (int k=9; k >=9; k-=2) {
        //MyPmap<NDIM> temppmap(world);
        FunctionDefaults<NDIM>::set_pmap(std::shared_ptr<MyPmap<NDIM> >(new MyPmap<NDIM>(world)));
        /*
        if (world.rank() == 0) {
            FunctionDefaults<NDIM>::pmap->print();
        }
        */
        int n = 2;
        //int n = 7;
        //double thresh = 1e-12;
        //double thresh = 5e-4;
        double thresh = 1e-3;
        if (world.rank() == 0) {
            printf("k=%d:\n", k);
        }
        double t0 = MPI::Wtime();
        //Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor).norefine().initial_level(n).k(k);
        Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor).thresh(thresh).initial_level(n).k(k);
        if (world.rank() == 0) {
            print("just made function f");
        }
        double t1 = MPI::Wtime();
        double err2 = f.err(*functor);
        std::size_t size = f.size();
        std::size_t tree_size = f.tree_size();
        std::size_t maxsize = f.max_nodes();
        std::size_t minsize = f.min_nodes();
        std::size_t maxdepth = f.max_depth();
        double t2 = MPI::Wtime();
        if (world.rank() == 0) {
            printf("   n=%d err=%.2e #coeff=%.2e tree_size=%.2e max_depth=%.2e max_nodes=%.2e min_nodes=%.2e log(err)/(n*k)=%.2e\n",
                   n, err2, double(size), double(tree_size), double(maxdepth), double(maxsize), double(minsize), abs(log(err2)/n/k));
        }
        world.gop.fence();
        double t3 = MPI::Wtime();
        Function<T,NDIM> g = copy(f);
        world.gop.fence();
        double t4 = MPI::Wtime();
        if (world.rank() == 0) print("ABOUT TO COMPRESS");
        f.compress(true);
        if (world.rank() == 0) print("ABOUT TO FENCE");
        world.gop.fence();
        double t5 = MPI::Wtime();
        if (world.rank() == 0) print("ABOUT TO RECON");
        f.reconstruct(true);
        if (world.rank() == 0) print("ABOUT TO FENCE");
        world.gop.fence();
        double t6 = MPI::Wtime();
        if (world.rank() == 0) print("ABOUT TO CREATE LOADBAL");
        //	typedef Cost (*my_default_fun<3,double>) (const Key<3>&, const FunctionNode<double,3>&);
        LoadBalImpl<NDIM> lb(g, lbcost<T,NDIM>(1.0, 1.0));
        //       	LoadBalImpl<NDIM> lb(g, &my_default_fun<T,NDIM>(Key<NDIM>, FunctionNode<T,NDIM>));
        //	LoadBalImpl<NDIM> lb(g, &f_ptr);
        if (world.rank() == 0) print("ABOUT TO FENCE");
        world.gop.fence();
        double t7 = MPI::Wtime();
        if (world.rank() == 0) print("ABOUT TO DO LOADBAL");
        FunctionDefaults<NDIM>::set_pmap(lb.load_balance());
        /*
        if (world.rank() == 0) {
            FunctionDefaults<NDIM>::pmap->print();
        }
        */
        if (world.rank() == 0) print("ABOUT TO FENCE");
        world.gop.fence();
        Function<T,NDIM> h = FunctionFactory<T,NDIM>(world).functor(functor).thresh(thresh).initial_level(n).k(k);
        double t8 = MPI::Wtime();
        f = copy(g,FunctionDefaults<NDIM>::get_pmap());
        double t9 = MPI::Wtime();
        double err21 = h.err(*functor);
        std::size_t size1 = h.size();
        std::size_t tree_size1 = h.tree_size();
        std::size_t maxsize1 = h.max_nodes();
        std::size_t minsize1 = h.min_nodes();
        std::size_t maxdepth1 = h.max_depth();
        double t10 = MPI::Wtime();
        if (world.rank() == 0) {
            printf("   n=%d err=%.2e #coeff=%.2e tree_size=%.2e max_depth=%.2e max_nodes=%.2e min_nodes=%.2e log(err)/(n*k)=%.2e\n",
                   n, err21, double(size1), double(tree_size1), double(maxdepth1), double(maxsize1), double(minsize1), abs(log(err21)/n/k));
        }
        world.gop.fence();
        double t11 = MPI::Wtime();
        h.compress(true);
        world.gop.fence();
        double t12 = MPI::Wtime();
        h.reconstruct(true);
        world.gop.fence();
        double t13 = MPI::Wtime();


        if (world.rank() == 0) {
            madness::print("for f with k =", k, "n =", n, "and thresh =",
                           thresh, "running on", world.nproc(), "processors:");
            madness::print("Routine            |  Time");
            madness::print("-------------------+--------------");
            madness::print("Init f             | ", t1-t0);
            madness::print("Diagnostics for f  | ", t2-t1);
            madness::print("Copy f to g        | ", t4-t3);
            madness::print("Compress f         | ", t5-t4);
            madness::print("Reconstruct f      | ", t6-t5);
            madness::print("Create LoadBalImpl | ", t7-t6);
            madness::print("Load Balance       | ", t8-t7);
            madness::print("Copy g to f        | ", t9-t8);
            madness::print("Diagnostics for g  | ", t10-t9);
            madness::print("Compress g         | ", t12-t11);
            madness::print("Reconstruct g      | ", t13-t12);
            madness::print("-------------------+--------------");
            madness::print("Total Time         | ", t13-t0);
            madness::print("");
        }
    }

    if (world.rank() == 0) print("test loadbal OK\n\n");
}


int main(int argc, char**argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);

    try {
        if (world.rank() == 0) print("before startup");
        startup(world,argc,argv);
        if (world.rank() == 0) print("before test_loadbal");
        test_loadbal<double, 3>(world);
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
