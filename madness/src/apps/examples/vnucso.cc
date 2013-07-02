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

//// \file examples/vnucso.cc
/// \brief Solves the spin-orbit nuclear potential problem

/*!

  \file vnucso.cc
  \brief Solves the Hartree-Fock equation for the 2-cosh potential with spin-orbit in Nuclear
  Density Functional Theory witough assumption on spatial symmetry.
  \ingroup examples

  Points of interest

  - Forming a Hamiltonian and a Fock matrix
  - Forming and solving a generalized eigensystem problem using LAPACK
  - Refining the representation of real and complex functions
  - Vectors of functions and operators, inner-product, gaxpy
  - Application of the Helmholtz bound-state Green function as vector of operators and functions
  - Projection and change of representation of functions from multiwavelets of degree k to k+1

  This is a more involved example than the Hydrogen and Helium.
  This example couples the traditional diagonalization approach with that of the integral equation approach.
  The details are described in:
  G. I. Fann , J. Pei, R. J. Harrison1, J. Jia, J. Hill1 , M. Ou, W. Nazarewicz, W. A. Shelton
  and N. Schunck, "Fast multiresolution methods for density functional theory in nuclear physics,"
  Journal of Physics, 180 (2009) 012080.


*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/vmra.h>
#include <mra/operator.h>
#include <constants.h>

using namespace madness;
using namespace std;

typedef Vector<double,3> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > real_functorT;
typedef std::shared_ptr< FunctionFunctorInterface<double_complex,3> > complex_functorT;
typedef Function<double,3> real_functionT;
typedef Function<double_complex,3> complex_functionT;
typedef FunctionFactory<double,3> real_factoryT;
typedef FunctionFactory<double_complex,3> complex_factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;

static const double L = 120.0;     // box is [-L,L]
static const double zeta = 7.5;   // potential wells at +/-zeta
static const double R1 = 2.0;     // potential parameter
static const double R2 = 2.0;     // potential parameter
static const double a1 = 1.0;     // potential parameter
static const double a2 = 1.0;     // potential parameter
static const double reduced = 0.04825964488415279478;
static const double V1 = -50.0*reduced;   // potential parameter
static const double V2 = -50.0*reduced;   // potential parameter
static const double lambda_correct = 0.0026608048208104861/reduced; // SO potential parameter

//static const double lambda_fudge = lambda_correct*100.0;

static       double lambda = lambda_correct;
static const double fac1=exp(-R1/a1);
static const double fac2=exp(-R2/a2);

struct Guess : FunctionFunctorInterface<double_complex,3> {
    const double z;        // z-coordinate center
    const double exponent; // exponent
    const int nx, ny, nz;  // powers of x, y, z ... only 0 or 1 supported!

    Guess(double z, double exponent, int nx, int ny, int nz)
        : z(z), exponent(exponent), nx(nx), ny(ny), nz(nz) {}

    double_complex operator()(const coordT& r) const {
        double rsq = r[0]*r[0] + r[1]*r[1] + (r[2]-z)*(r[2]-z);
        double psi = exp(-exponent*rsq);
        if (nx) psi *= r[0];
        if (ny) psi *= r[1];
        if (nz) psi *= (r[2]-z);
        return double_complex(psi,0.0);
    }
};

static double V(const coordT& r)
{
    const double x=r[0], y=r[1], z=r[2];
    const double zp=(z+zeta), zm = (z-zeta);
    const double rp = sqrt(x*x + y*y + zp*zp);
    const double rm = sqrt(x*x + y*y + zm*zm);

    return V1/(1.0 + fac1*cosh(rp/a1)) + V2/(1.0 + fac2*cosh(rm/a2));
}

void moments(World& world, const vector<complex_functionT>& u, const vector<complex_functionT>& v) {
  FunctionDefaults<3>::set_autorefine(true);
  complex_functionT rho = complex_factoryT(world);
  rho.compress();
  reconstruct(world, u);
  reconstruct(world, v);
  for (unsigned int i=0; i<u.size(); i++) {
    complex_functionT psisq = conj(u[i])*u[i] + conj(v[i])*v[i];
    psisq.compress();
    rho.gaxpy(1.0, psisq, 1.0);
  }
  double_complex moment1 = rho.trace();
  complex_functionT rhosq = (rho*rho).scale(2.);
  double_complex moment2a = rhosq.trace();
  double_complex moment2b = inner(rho,rho);
  double_complex moment3 = inner(rho,rhosq);
  if (world.rank() == 0) {
    print("  moment1", moment1);
    print("  moment2", moment2a, moment2b);
    print("  moment3", moment3);
  }
  FunctionDefaults<3>::set_autorefine(false);
}

void gaxpy1(World& world,
	    const double_complex alpha,
	    vector<complex_functionT>& a,
	    const double_complex beta,
	    const vector<complex_functionT>& b,
	    bool fence=true)
{
  MADNESS_ASSERT(a.size() == b.size());

  for (unsigned int i=0; i<a.size(); i++) {
    a[i] = alpha*a[i] + beta*b[i];
  }
  if (fence) world.gop.fence();
}

vector<poperatorT> make_bsh_operators(World& world, const Tensor<double>& evals, double tol)
{
    int n = evals.dim(0);
    vector<poperatorT> ops(n);
    for (int i=0; i<n; i++) {
        double eps = evals(i);
        if (eps > 0) eps = -0.05;
        double lo = 0.1*tol; // heuristic
        ops[i] = poperatorT(BSHOperatorPtr3D(world, sqrt(-eps), lo, tol));
    }
    return ops;
}

Tensor<double_complex> hamiltonian_matrix(World& world,
                                          const vector<complex_functionT>& u,
                                          const vector<complex_functionT>& v,
                                          const vector<complex_functionT>& Vu,
                                          const vector<complex_functionT>& Vv,
                                          const vector<complex_functionT> du[3],
                                          const vector<complex_functionT> dv[3])
{
    reconstruct(world, u);
    reconstruct(world, v);
    int n = u.size();
    Tensor<double_complex> r(n,n);

    for (int axis=0; axis<3; axis++) {
      r += matrix_inner(world, du[axis], du[axis], true);
      r += matrix_inner(world, dv[axis], dv[axis], true);
    }
    return r + matrix_inner(world, u, Vu, true) + matrix_inner(world, v, Vv, true);
}

template <typename T, int NDIM>
Cost lbcost(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) {
  return 1;
}

void apply_potential(World& world,
                     const real_functionT& V0,
                     const real_functionT& V0x,
                     const real_functionT& V0y,
                     const real_functionT& V0z,
                     const vector<complex_functionT>& u,
                     const vector<complex_functionT>& v,
                     vector<complex_functionT>& Vu,
                     vector<complex_functionT>& Vv,
                     vector<complex_functionT> du[3],
                     vector<complex_functionT> dv[3],
                     bool doso)
{
    const double_complex lam(lambda,0.0);
    const double_complex one(1.0,0.0);
    const double_complex I(0.0,1.0);
    const int x=0, y=1, z=2; // Attempt to make SO term more readable

    reconstruct(world, u);
    reconstruct(world, v);
    V0.reconstruct();
    V0x.reconstruct();
    V0y.reconstruct();
    V0z.reconstruct();

    for (int axis=0; axis<3; axis++) {
        complex_derivative_3d D = free_space_derivative<double_complex,3>(world, axis);
        du[axis] = apply(world, D, u);
        dv[axis] = apply(world, D, v);
    }

    Vu = mul(world, V0, u);
    Vv = mul(world, V0, v );

    if (!doso) return;

    gaxpy(world, one, Vu, -I*lam, mul(world, V0y, du[x]));
    gaxpy(world, one, Vu,  I*lam, mul(world, V0x, du[y]));
    gaxpy(world, one, Vu,    lam, mul(world, V0z, dv[x]));
    gaxpy(world, one, Vu, -I*lam, mul(world, V0z, dv[y]));
    gaxpy(world, one, Vu,   -lam, mul(world, V0x, dv[z]));
    gaxpy(world, one, Vu,  I*lam, mul(world, V0y, dv[z]));


    gaxpy(world, one, Vv,  I*lam, mul(world, V0y, dv[x]));
    gaxpy(world, one, Vv, -I*lam, mul(world, V0x, dv[y]));
    gaxpy(world, one, Vv,   -lam, mul(world, V0z, du[x]));
    gaxpy(world, one, Vv, -I*lam, mul(world, V0z, du[y]));
    gaxpy(world, one, Vv,    lam, mul(world, V0x, du[z]));
    gaxpy(world, one, Vv,  I*lam, mul(world, V0y, du[z]));

}

void normalize2(World& world, vector<complex_functionT>& u, vector<complex_functionT>& v) {
  vector<double> unorm = norm2s(world,u);
  vector<double> vnorm = norm2s(world,v);
  vector<double> normu(u.size());

  for (unsigned int i=0; i<u.size(); i++) {
    normu[i] = sqrt(unorm[i]*unorm[i] + vnorm[i]*vnorm[i]);
  }
  for (unsigned int i=0; i<u.size(); i++) {
    u[i].scale(double_complex(1.0/normu[i],0.));
    v[i].scale(double_complex(1.0/normu[i],0.));
  }
}



void doit(World& world) {
    // We are not using time-reversal symmetry and hence are doing
    // about 2x too much work.  With a little bit more complexity we
    // could readily exploit it.  Time reversal symmetry implies that
    // if (u,v) (two-component wave function) is an eigen function
    // then (-v*, u*) is a degenerate solution.  We could use this
    // symmetry by ensuring that after updating (and before
    // diagonalization) the higher energy states are orthogonal to
    // lower states and their time-reversed form.

    if (world.rank() == 0) print("entered solver at time", wall_time());
    long k = 5;              // wavelet order
    double thresh = 1e-3;    // precision for wave function
    coordT zero;

    zero[0]=0.; zero[1]=0.; zero[2]=0.;

    // FunctionDefaults<3>::k = k;
    //    FunctionDefaults<3>::thresh = thresh;
    // FunctionDefaults<3>::refine = true;
    // FunctionDefaults<3>::autorefine = false;
    // FunctionDefaults<3>::initial_level = 2;
    // FunctionDefaults<3>::truncate_mode = 1;

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_initial_level(2);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_cubic_cell(-L,L);

    // for (int i=0; i<3; i++) {
    // FunctionDefaults<3>::cell(i,0) = -L;
    // FunctionDefaults<3>::cell(i,1) =  L;
    // }
    if (world.rank() == 0) print("Making guesses");
    vector<complex_functionT> u;
    vector<complex_functionT> v;
    vector<complex_functionT> w;

    // Initial guess is as for SO-free but first with (u,0) and then (0,v)

    // Could this get any more verbose?  This needs to be fixed.
    u.push_back(complex_factoryT(world).functor(complex_functorT(new Guess( zeta, 0.15, 0, 0, 0))).nofence());  // s
    v.push_back(complex_functionT(complex_factoryT(world)));

    u.push_back(complex_factoryT(world).functor(complex_functorT(new Guess(-zeta, 0.15, 0, 0, 0))).nofence());  // s
    v.push_back(complex_functionT(complex_factoryT(world)));

    u.push_back(complex_factoryT(world).functor(complex_functorT(new Guess( zeta, 0.15, 0, 0, 1))).nofence());  // pz
    v.push_back(complex_functionT(complex_factoryT(world)));

    u.push_back(complex_factoryT(world).functor(complex_functorT(new Guess(-zeta, 0.15, 0, 0, 1))).nofence());  // pz
    v.push_back(complex_functionT(complex_factoryT(world)));
    u.push_back(complex_factoryT(world).functor(complex_functorT(new Guess( zeta, 0.15, 0, 1, 0))).nofence());  // py
    v.push_back(complex_functionT(complex_factoryT(world)));
    u.push_back(complex_factoryT(world).functor(complex_functorT(new Guess(-zeta, 0.15, 0, 1, 0))).nofence());  // py
    v.push_back(complex_functionT(complex_factoryT(world)));
    u.push_back(complex_factoryT(world).functor(complex_functorT(new Guess( zeta, 0.15, 1, 0, 0))).nofence());  // px
    v.push_back(complex_functionT(complex_factoryT(world)));

    u.push_back(complex_factoryT(world).functor(complex_functorT(new Guess(-zeta, 0.15, 1, 0, 0))).nofence());  // px
    v.push_back(complex_functionT(complex_factoryT(world)));
    u.push_back(complex_factoryT(world).functor(complex_functorT(new Guess( zeta, 0.03, 0, 0, 0))).nofence());  // s
    v.push_back(complex_functionT(complex_factoryT(world)));
    u.push_back(complex_factoryT(world).functor(complex_functorT(new Guess(-zeta, 0.03, 0, 0, 0))).nofence());  // s
    v.push_back(complex_functionT(complex_factoryT(world)));
    //
    v.push_back(complex_factoryT(world).functor(complex_functorT(new Guess( zeta, 0.15, 0, 0, 0))).nofence());  // s
    u.push_back(complex_functionT(complex_factoryT(world)));
    v.push_back(complex_factoryT(world).functor(complex_functorT(new Guess(-zeta, 0.15, 0, 0, 0))).nofence());  // s
    u.push_back(complex_functionT(complex_factoryT(world)));
    v.push_back(complex_factoryT(world).functor(complex_functorT(new Guess( zeta, 0.15, 0, 0, 1))).nofence());  // pz
    u.push_back(complex_functionT(complex_factoryT(world)));
    v.push_back(complex_factoryT(world).functor(complex_functorT(new Guess(-zeta, 0.15, 0, 0, 1))).nofence());  // pz
    u.push_back(complex_functionT(complex_factoryT(world)));
    v.push_back(complex_factoryT(world).functor(complex_functorT(new Guess( zeta, 0.15, 0, 1, 0))).nofence());  // py
    u.push_back(complex_functionT(complex_factoryT(world)));
    v.push_back(complex_factoryT(world).functor(complex_functorT(new Guess(-zeta, 0.15, 0, 1, 0))).nofence());  // py
    u.push_back(complex_functionT(complex_factoryT(world)));
    v.push_back(complex_factoryT(world).functor(complex_functorT(new Guess( zeta, 0.15, 1, 0, 0))).nofence());  // px
    u.push_back(complex_functionT(complex_factoryT(world)));
    v.push_back(complex_factoryT(world).functor(complex_functorT(new Guess(-zeta, 0.15, 1, 0, 0))).nofence());  // px
    u.push_back(complex_functionT(complex_factoryT(world)));
    v.push_back(complex_factoryT(world).functor(complex_functorT(new Guess( zeta, 0.03, 0, 0, 0))).nofence());  // s
    u.push_back(complex_functionT(complex_factoryT(world)));
    v.push_back(complex_factoryT(world).functor(complex_functorT(new Guess(-zeta, 0.03, 0, 0, 0))).nofence());  // s
    u.push_back(complex_functionT(complex_factoryT(world)));
    world.gop.fence();
    normalize2(world, u, v);


    int nvec=u.size();
    Tensor<double_complex> H1(nvec,nvec);
    Tensor<double_complex> S1(nvec,nvec);

    Tensor<double_complex> c;
    Tensor<double> e, e1(nvec);
    double maxerr;
    real_functionT V0, V0x, V0y, V0z;
    //double shift=0.;
    bool doso = false; // turned on once converged to 1e-4
    print(" u size ",  u.size(), "v size ", v.size(), "\n");

    doitagain:

    while (thresh > 0.9e-10) {
      if (world.rank() == 0) {
	print("\n Solving with thresh",thresh,"k",k,"at time", wall_time(), "\n");
            printf("  iter   root        energy         err      time\n");
            printf("  ----   ----  ------------------ -------   ------\n");
      }

      V0 = real_factoryT(world).f(V).thresh(thresh);
      V0x = real_derivative_3d(world,0)(V0);
      V0y = real_derivative_3d(world,1)(V0);
      V0z = real_derivative_3d(world,2)(V0);

      // If it does not converge in 10 iters it probably needs more precision.
      for (int iter=0; iter<10; iter++) {

	// loadbalance(world, u, v, V0, V0x, V0y, V0z);

	vector<complex_functionT> Vu, Vv, du[3], dv[3];

	apply_potential(world, V0, V0x, V0y, V0z, u, v, Vu, Vv, du, dv, doso);
	world.gop.fence(); // free memory

	Tensor<double_complex> H = hamiltonian_matrix(world, u, v, Vu, Vv, du, dv);
	Tensor<double_complex> S = matrix_inner(world, u, u) + matrix_inner(world, v, v);

	if ( iter==0) {
	  for (int iii = 0; iii < nvec; iii++ )
	    for (int jjj = 0; jjj < nvec; jjj++ ){
	      H1(jjj,iii) = double_complex(0.5,0.0)*(H(iii,jjj)+H(jjj,iii));
	    }
	}
	else
	  for (int iii = 0; iii < nvec; iii++ )
	    for (int jjj = 0; jjj < nvec; jjj++ ){
	      H1(jjj,iii) = H(jjj,iii);
	    }

	if ( iter==0) {
	  for (int iii = 0; iii < nvec; iii++ )
	    for (int jjj = 0; jjj < nvec; jjj++ ){
	      S1(jjj,iii) = double_complex(0.5,0.0)*(S(iii,jjj)+S(jjj,iii));
	    }
	}
	else
	  for (int iii = 0; iii < nvec; iii++ )
	    for (int jjj = 0; jjj < nvec; jjj++ ){
	      S1(jjj,iii) = S(jjj,iii);
	    }

	world.gop.fence(); // free memory

	sygv(H1, S1, 1, c, e);

	for (int axis=0; axis<3; axis++) {
	  du[axis].clear();
	  dv[axis].clear();
	}

	world.gop.fence(); // free memory
	e = e(Slice(0,nvec-1));
	u  = transform(world, u,  c(_,Slice(0,nvec-1)));
	v  = transform(world, v,  c(_,Slice(0,nvec-1)));

	Vu = transform(world, Vu, c(_,Slice(0,nvec-1)));
	Vv = transform(world, Vv, c(_,Slice(0,nvec-1)));
	world.gop.fence(); // free memory

	truncate(world, u);
	truncate(world, v);
	truncate(world, Vu);
	truncate(world, Vv);
	world.gop.fence();

	normalize2(world, u, v);

	world.gop.fence(); // free memory

	for (int iii = 0; iii < nvec; iii++ )
	  e1[iii] = e[iii]; // +shift;


	for (int i=0; i<nvec; i++) {
	  Vu[i].refine();
	  Vv[i].refine();
	}

	vector<poperatorT> ops = make_bsh_operators(world, e1, thresh);

	vector<complex_functionT> u_new = apply(world, ops, Vu);
	vector<complex_functionT> v_new = apply(world, ops, Vv);

	normalize2(world, u_new, v_new);

	Vu.clear();
	Vv.clear();
	world.gop.fence();

	vector<double> rnormu = norm2s(world,add(world, u, u_new));
	vector<double> rnormv = norm2s(world,add(world, v, v_new));
	vector<double> rnorm(nvec);

	for (int i=0; i<nvec; i++)
	  rnorm[i] = sqrt(rnormu[i]*rnormu[i] + rnormv[i]*rnormv[i]);

	maxerr= 0.;
	for (int i=0; i<nvec; i++) {
	  if (world.rank() == 0) printf("  %3d    %3d  %18.12f  %.1e  %7.1f\n",
					iter, i, e[i]/reduced, rnorm[i], wall_time());
	  maxerr = max(maxerr,rnorm[i]);
	}

	u = u_new;
	v = v_new;

	u_new.clear();
	v_new.clear();

	if (maxerr < 1.e-3 && !doso) {
	  // We initially converged to a threshold of 1e-4 without SO
	  // and now we repeat the 1e-4 iteration with SO and continue
	  // with SO on.
	  if (world.rank() == 0) print("\n Turning on the SO interaction \n");
	  doso = true;
	  goto doitagain;  // Gotta love it
	  break;
	}
      }

      moments(world, u, v);

      thresh *= 1e-1;
      k += 1;
      FunctionDefaults<3>::set_k(k);
      FunctionDefaults<3>::set_thresh(thresh);
      reconstruct(world,u);
      reconstruct(world,v);
      for (unsigned int i=0; i<u.size(); i++) {
	u[i] = madness::project(u[i], k, thresh);
	v[i] = madness::project(v[i], k, thresh);
      }
    }
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);

//     cpu_set_t mask;
//     CPU_ZERO(&mask);
//     CPU_SET(MPI::COMM_WORLD.Get_rank(), &mask);
//     if( sched_setaffinity( 0, sizeof(mask), &mask ) == -1 ) {
//         printf("WARNING: Could not set CPU Affinity, continuing...\n");
//     }

    try {
        doit(world);
    } catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    } catch (char* s) {
        print(s);
        error("caught a string exception");
    } catch (const char* s) {
        print(s);
        error("caught a string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }

    finalize();
    return 0;
}
