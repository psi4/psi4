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

/*!
	\file tdse1d.cc
	\brief Example propagation of TDSE (translating atom) using various propagators
	\defgroup exampletdse1d Solves a 1D time-dependent Schr&ouml;dinger equation using splitting and semi-group approaches with the free-particle propagator.
  \ingroup examples

  The source is <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/tdse1d.cc >here</a>.

  \par Points of interest
  - convolution with the Green's function/free-particle propagator
  - use free-particle propagators construct different schemes to find the solution
  - collocation methods on quadrature points

  \par Background
  This illustrates solution of a time-dependent Schr&ouml;dinger equation.

  We solve the following PDE
  \f[
  -\nabla^2 \psi(x,t) + V(x,t) \psi(x,t) = i\psi_t(x,t)
  \f]
  where the potential is
  \f[
  V(x,t) = -8 \exp((x - v \cdot t)^2)
  \f]
  and the velocity \f$ v \f$ is given in the code.

  \par Implementation

  Splitting based schemes, such as Trotter and Chin-Chen, can be found in existing literatures.


  The quadrature collocation method is based on the semi-group form of the equation
  \f[
     \psi(x,t) = \psi(x,0) * G(x,t) - \int_0^t V(x,s) \cdot \psi(x,s) * G(x,t-s) ds
  \f]
  where \f$ G \f$ is the Green's function/free-particle propagator
  \f[
     \left( - \frac{d^2}{dx^2} - i \frac{d}{dt} \right) G(x) = \delta(x)).
  \f]

  To find \f$ \psi(x,t) \f$, all temporal integrals are computed by a \f$n\f$ point Gauss-Legendre quadrature rule and we need to
  calculate employ a simple fixed-point iteration to the self-consistent
  solutions at the \f$n\f$ quadrature points on the intervel \f$ [0,t] \f$.
  All \f$ \psi \f$ involved in computing the integrals over the subintervels
  are interpolated by the \f$n\f$ values on the largest intervel \f$ [0,t] \f$.

  The fixed-point iteration is applied to the correction term of the semi-group formulation,
 \f[
     \psi^{m+1}(x,t) - \psi^{m}(x,t) = - \int_0^t V(x,s) \cdot ( \psi^{m}(x,s) - \psi^{m-1}(x,s)) * G(x,t-s) ds
 \f]
 or
   \f[
     \delta^{m+1}(x,t) = - \int_0^t V(x,s) \cdot  \delta^{m}(x,s) * G(x,t-s) ds
 \f]
 where \f$ \delta^m \f$ is the \f$ m_{th} \f$ correction term.

  A much more efficient scheme would involve use of a non-linear
  equation solver instead of simple iteration.

 Once we have the solutions at the \f$n\f$ quadrature points on \f$ [0,t] \f$, quadrature rule is used to construct the solution at \f$t\f$.

 \par Reference
  Preprint.

*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/qmprop.h>
#include <mra/operator.h>
#include <constants.h>
#include <string>

using namespace madness;

// This to test various implementaitons of the 1d propagator
#ifdef USE_NSFORM
  typedef SeparatedConvolution<double_complex,1> complex_operatorT;
  typedef std::shared_ptr<complex_operatorT> pcomplex_operatorT;
  #define APPLY(G,psi) apply(*G,psi)
  #define MAKE_PROPAGATOR(world, t) qm_free_particle_propagatorPtr<1>(world, k, c, t)

#else
  typedef Convolution1D<double_complex> complex_operatorT;
  typedef std::shared_ptr<complex_operatorT> pcomplex_operatorT;

  complex_function_1d APPLY(complex_operatorT* q1d, const complex_function_1d& psi) {
      psi.reconstruct();

      psi.broaden();
      psi.broaden();
      psi.broaden();
      psi.broaden();

      psi.broaden();
      // psi.broaden();
      // psi.broaden();
      // psi.broaden();

      complex_function_1d r = apply_1d_realspace_push(*q1d, psi, 0);
      r.sum_down();
      return r;
  }
  #define MAKE_PROPAGATOR(world, t) qm_1d_free_particle_propagator(k, c, t, 2.0*L)

#endif


// Simulation parameters
static const double L = 100.0; // Simulation in [-L,L]
static const double x0 = -L + 10.0; // Initial position of the atom
static const double energy_exact = -6.188788775728796797594788; // From Maple
static const long k = 20;        // wavelet order
static const double thresh = 1e-8; // precision
static const double velocity = 3.0;
//static const double eshift = energy_exact - 0.5*velocity*velocity; // Use this value to remove rotating phase
static const double eshift = 0.0;

// Estimate the bandwidth and largest practical time step using G0
static double ctarget = 20.0; // Estimated from FT of exact solution ... was 20
static double c = 1.86*ctarget;
static double tcrit = 2*constants::pi/(c*c);

static double current_time = 0.0; // Lazy but easier than making functors for everything

/////////////////////////////////// For quadrature collocations ///////////////////////////////
// global vars for the laziness
static const double_complex I = double_complex(0,1);
int np; // number of quadrature pts
Tensor<double> B, tc;
pcomplex_operatorT G;
std::vector<pcomplex_operatorT> Gs, Gss, Gtrs;
const int maxiter = 20;
double fix_iter_tol = thresh*10.0;

struct refop {
    bool operator()(FunctionImpl<double_complex,1>* impl, const Key<1>& key, const Tensor<double_complex>& t) const {
        double tol = impl->truncate_tol(impl->get_thresh(), key);
        double lo, hi;
        impl->tnorm(t, &lo, &hi);
        return hi > tol;;
    }
    template <typename Archive> void serialize(Archive& ar) {}
};


// Position of atom at current time
double atom_position() {
    return x0 + velocity*current_time;
}

// Exact solution ... (H-E)psi is accurate to 2e-7 or better inside Maple
double_complex psi_exact(const coord_1d& r) {
    const double x = r[0] - atom_position();

    if (fabs(x) > 9.0) return 0.0;

    const double xsq = x*x;

    // Horner's form for stability ... yes, it is 70-order polyn ... don't panic ... all is OK
    const double psi = exp(-1.30*xsq)*(-1.02151632756275513018719494826+(.522210612113707231536059971069+(-.378478352719362210706386739834+(.183732263756009855463432582593+(-0.866826311335724362186706464428e-1+(0.364601910940641762284283418688e-1+(-0.144289291226899801775738667691e-1+(0.536464813679927807734520149659e-2+(-0.188945345474975346004237994967e-2+(0.628725522158030920219217207400e-3+(-0.195986657875763402765072721284e-3+(0.563993909330309881048156461300e-4+(-0.147273758530730072646826354339e-4+(0.343202525037692780236348165792e-5+(-7.03765391498970506694108123682e-7+(1.25577395245191089173671652503e-7+(-1.93270666918809750376161513191e-8+(2.54624395753990033043923369519e-9+(-2.84968491109847694778452937732e-10+(2.68398879018928855922002181223e-11+(-2.09811331054703124038096336404e-12+(1.32869596552058891546970129944e-13+(-6.47453843054578193348503832223e-15+(2.08146181239264250219157355910e-16+(-8.27692933283146365185698371000e-19+(-4.21076100567125673420604632290e-19+(3.34873683682880953223636900553e-20+(-1.62840449033428572820374528252e-21+(5.80781234060753332169593869292e-23+(-1.58754665363960354625758075480e-24+(3.34792415609741263184593450566e-26+(-5.37727687523701580992550859153e-28+(6.37706272777548158355242766766e-30+(-5.27154378837014394526580093435e-32+(2.71430533655911926094030668033e-34-6.55694230766452931026127498540e-37*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq)*xsq);

    // Galilean translation factor
    const double arg = x*velocity - (energy_exact - velocity*velocity*0.5 - eshift)*current_time;
    const double_complex tranfac = exp(double_complex(0,arg));
    return psi*tranfac;
}

// Time-dependent potential for translating atom
double V(const coord_1d& r) {
    const double x = r[0] - atom_position();
    if (fabs(x) > 6.2) return 0.0;

    return -8.0*exp(-x*x) - eshift;
}

// (dV/dr)^2
double dVsq(const coord_1d& r) {
    const double x = r[0] - atom_position();
    if (fabs(x) > 4.5) return 0.0;

    double dv = 16.0*x*exp(-x*x);
    return dv*dv;
}


template<typename T, int NDIM>
struct unaryexp {
    void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
        UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = exp(*_p0););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};

// exp(-i (vcoeff*V + dcoeff*dV^2))
complex_function_1d expV(World& world, double vcoeff, double dcoeff) {
    real_function_1d potn = real_factory_1d(world).f(V);
    potn.scale(vcoeff);
    if (dcoeff) {
        real_function_1d d = real_factory_1d(world).f(dVsq).initial_level(10);
        potn.compress(); d.compress();
        potn.gaxpy(1.0, d, dcoeff);   //delta = 1.5*t;/
    }
    complex_function_1d expV = double_complex(0.0,-1.0)*potn;
    expV.unaryop(unaryexp<double_complex,1>());
    expV.truncate();
    return expV;
}

// Evolve forward one time step using symplectic gradient 4-th-order local accurate
// G should G0(tstep/2)
// CC is xi=0.0, chi=1.0/72.0
// optimal is xi=-17/18000 chi=71/4500
complex_function_1d sympgrad4(World& world,
                            const complex_function_1d& psi0,
                            const double tstep,
                            const double xi,
                            const double chi) {

    static complex_operatorT* G = 0;
    if (!G) G = MAKE_PROPAGATOR(world, tstep*0.5);

    const double lambda = 1.0/6.0;
    complex_function_1d psi;

    psi = expV(world, tstep*lambda, -xi*tstep*tstep*tstep)*psi0;           psi.truncate();
    psi = APPLY(G, psi);                                                   psi.truncate();

    current_time += 0.5*tstep;

    psi = expV(world, tstep*(1.0-2.0*lambda), -chi*tstep*tstep*tstep)*psi; psi.truncate();

    current_time += 0.5*tstep;

    psi = APPLY(G,psi);                                                    psi.truncate();
    psi = expV(world, tstep*lambda, -xi*tstep*tstep*tstep)*psi;            psi.truncate();

    return psi;
}

complex_function_1d sympgrad6(World& world,
                            const complex_function_1d& psi0,
                            const double tstep) {

    const double rho = 0.1097059723948682e+00;
    const double theta = 0.4140632267310831e+00;
    const double nu = 0.2693315848935301e+00;
    const double lambda = 0.1131980348651556e+01;
    const double chi = -0.1324638643416052e-01;
    const double mu =  0.8642161339706166e-03;

    static complex_operatorT *Grho=0, *Gtheta=0, *Gmid=0;

    if (Grho == 0) {
        Grho   = MAKE_PROPAGATOR(world, tstep*rho);
        Gtheta = MAKE_PROPAGATOR(world, tstep*theta);
        Gmid   = MAKE_PROPAGATOR(world, tstep*(1.0-2.0*(theta+rho))/2.0);
    }

    complex_function_1d psi;

    psi = APPLY(Grho, psi0);                                                psi.truncate();
    current_time += rho*tstep;
    psi = expV(world, tstep*nu, -mu*tstep*tstep*tstep)*psi;                  psi.truncate();
    psi = APPLY(Gtheta, psi);                                               psi.truncate();
    current_time += theta*tstep;
    psi = expV(world, tstep*lambda, 0.0)*psi;                                psi.truncate();
    psi = APPLY(Gmid, psi);                                                 psi.truncate();
    current_time += tstep*(1.0-2.0*(theta+rho))/2.0;
    psi = expV(world, tstep*(1.0-2.0*(lambda+nu)), -chi*tstep*tstep*tstep)*psi; psi.truncate();
    psi = APPLY(Gmid, psi);                                                 psi.truncate();
    current_time += tstep*(1.0 - 2.0*(theta+rho))/2.0;
    psi = expV(world, tstep*lambda, 0.0)*psi;                                psi.truncate();
    psi = APPLY(Gtheta, psi);                                               psi.truncate();
    current_time += theta*tstep;
    psi = expV(world, tstep*nu, -mu*tstep*tstep*tstep)*psi;                  psi.truncate();
    psi = APPLY(Grho, psi);                                                 psi.truncate();
    current_time += rho*tstep;

    return psi;
}


// Evolve forward one time step using Trotter ... G = G0(tstep/2)
complex_function_1d trotter(World& world, const complex_function_1d& psi0, const double tstep, complex_operatorT* G0 = 0) {
    static complex_operatorT* G = 0;
	if (G0) G = G0;
    if (!G) G = MAKE_PROPAGATOR(world, tstep*0.5);

    complex_function_1d psi = APPLY(G, psi0);    psi.truncate();

    current_time += 0.5*tstep;

    psi = expV(world, tstep, 0.0)*psi;         psi.truncate();

    current_time += 0.5*tstep;

    psi = APPLY(G,psi);                        psi.truncate();

    return psi;
}

void print_info(World& world, const complex_function_1d& psi, int step) {
    real_function_1d potn = real_factory_1d(world).f(V).truncate_on_project();
    complex_derivative_1d D = free_space_derivative<double_complex,1>(world,0);
    complex_function_1d dpsi = D(psi);
    double ke = inner(dpsi,dpsi).real() * 0.5;
    double pe = psi.inner(psi*potn).real();
    double norm = psi.norm2();

    complex_function_1d psiX = complex_factory_1d(world).f(psi_exact);
    double err = (psiX - psi).norm2();

    if (world.rank() > 0) return;
    if ((step%40) == 0) {
        printf("\n");
        printf(" step    time      atom x              norm               kinetic              potential             energy              err norm     depth   size  \n");
        printf("------ -------- ------------     ----------------     ----------------     ----------------     ----------------     ---------------- ----  --------\n");
    }
    printf("%6d %8.4f %12.8f %20.13f %20.13f %20.13f %20.13f %20.13f %4d %9ld\n",
           step, current_time, atom_position(), norm, ke, pe, ke+pe, err, int(psi.max_depth()), long(psi.size()));
}


static void readin(int np) {
    Tensor<double> BB(np), tctc(np);
    gauss_legendre(np, 0, 1, tctc.ptr(), BB.ptr());
    B=BB; tc=tctc;
}

//generate j_th lagrange interpolating coefficients
double icoeff(const int np, const int j, const double t) {
    double dum = 1;
    for (int i=  0; i< j; ++i) dum *= (t-tc[i])/(tc[j]-tc[i]);
    for (int i=j+1; i<np; ++i) dum *= (t-tc[i])/(tc[j]-tc[i]);
    return dum;
}

//evaluate the interpolating polynomail of ps at t
template<typename T> T myp(const std::vector<T>& ps, const double t) {
    int np = ps.size();
    T p = ps[0]*icoeff(np, 0, t);
    for (int j=1; j<np; ++j) p += ps[j]*icoeff(np, j, t);
    return p;
}

// Evolve forward one time step using quadrature rules
complex_function_1d q_c(World& world, const int np, const complex_function_1d psi0, const double tstep) {
    //can be more vectorized.

    std::vector<complex_function_1d> ps(np), ps1(np);
    //    for (int i=0; i<np; ++i) ps[i] = copy(psi0);

    double tdum = current_time;
    complex_function_1d pdum;
    std::vector<complex_function_1d> qs(np); //qs stores current solution values at quadrature points
    std::vector<real_function_1d> Vs(np); //Vs stores the function values of V at quadrature points
    std::vector< std::vector<real_function_1d> > Vss(np);  //Vas stores the function values of V at SUB quadrature points
    for (int i=0; i<np; ++i) {
        current_time = tdum + tstep*tc[i];
        Vs[i] = real_factory_1d(world).f(V).truncate_on_project();
        Vs[i].truncate();
        Vss[i].resize(np);
        for (int k=0; k<np; ++k) {
            current_time = tdum + tstep*tc[i]*tc[k];
            Vss[i][k]=real_factory_1d(world).f(V).truncate_on_project();
            Vss[i][k].truncate();
        }
        //~ ps[i] = APPLY(Gs[np - i - 1].get(), psi0).truncate();
        //~ qs[i]=copy(ps[i]);
        ps1[i] = APPLY(Gs[np - i - 1].get(), psi0).truncate();  //ps1 temporarily holds the solutions after first correction. Part 1/2
        current_time = tdum + (i?tstep*tc[i-1]:0);
        //Use trotter to generate initial guess at quadrature points; note that current_time is implicitly changed.
        qs[i] =  trotter(world, (i?qs[i-1]:psi0), tstep*(tc[i]-(i?tc[i-1]:0)),Gtrs[i].get());
    }

    compress(world, ps1);
    for (int i=0; i<np; ++i) for (int k=0; k<np; ++k) {
		//ps1[i].gaxpy(1.0,APPLY(Gss[i*np+k].get(), ).compress(), -tstep*tc[i]*I*B[k]);
		complex_function_1d tmp = (Vss[i][k]*myp(qs,tc[i]*tc[k])).scale(-tstep*tc[i]*I*B[k]).truncate();
		ps1[i].gaxpy(1.0,APPLY(Gss[i*np+k].get(), tmp).compress(), 1.0);  //ps1 temporarily holds the solutions after first correction. Part 2/2
    }

	// Preparation for one time step is done now, and preceed to iterations.

    ps = sub(world, ps1, qs);  // initialize ps, which is the correction term
    qs = ps1;  // update the solutions with the first correction.

    // now the fix-pt iterations for \psi's or qs' on the quadrature pts
    fix_iter_tol = thresh / tstep;
    if (world.rank() == 0) printf("fix iters ");
    double err = norm2(world, ps1)/np ;

    //~ ps=copy(world, qs);
    for (int j=0; j<maxiter; ++j) {
        if (world.rank() == 0) printf(" %6.0e",err);
        if (err <= fix_iter_tol) break;
        err = 0;
        ps1 = zero_functions<double_complex,1>(world, np); // ps1 holds the newly computed correction term
        compress(world,ps1);
        for (int i=0; i<np; ++i) {
            for (int k=0; k<np; ++k) { // compute the correction term
                //ps1[i].gaxpy(1.0, APPLY(Gss[i*np+k].get(), Vss[i][k]*myp(ps,tc[i]*tc[k])).compress(), -tstep*tc[i]*I*B[k]);
                complex_function_1d tmp = (Vss[i][k]*myp(ps,tc[i]*tc[k])).scale(-tstep*tc[i]*I*B[k]).truncate();
                ps1[i].gaxpy(1.0, APPLY(Gss[i*np+k].get(), tmp).compress(), 1.0);
            }
            err += ps1[i].truncate().norm2();
        }
        err /= np;
        gaxpy(world, 1.0, qs, 1.0, ps1);  // accumulate the corrections into the solutions.
        ps = ps1; //copy(world, ps1);
    }
    if (world.rank() == 0) printf("\n");

    // hopefully the iteration converged correctly and now use the semi-group formula with quadrature rule applied to the integral to compute the result.
    pdum = APPLY(G.get(), psi0);
    pdum.compress();
    for (int k=0; k<np; ++k) {
        //pdum.gaxpy(1.0, APPLY(Gs[k].get(), Vs[k]*qs[k]).compress(), -I*B[k]*tstep);
        complex_function_1d tmp = (Vs[k]*qs[k]).scale(-I*B[k]*tstep);
        pdum.gaxpy(1.0, APPLY(Gs[k].get(), tmp).compress(), 1.0);
    }

    current_time = tdum + tstep;

    return pdum.truncate();
}


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(8);

    FunctionDefaults<1>::set_k(k);                 // Wavelet order
    FunctionDefaults<1>::set_thresh(thresh);       // Accuracy
    FunctionDefaults<1>::set_autorefine(false);
    FunctionDefaults<1>::set_cubic_cell(-L,L);
    FunctionDefaults<1>::set_initial_level(8);     // Initial projection level
    FunctionDefaults<1>::set_truncate_mode(0);

    complex_function_1d psi = complex_factory_1d(world).f(psi_exact);
    psi.truncate();

    double tstep = 0.0;
    int selection = -1;

    if (world.rank() == 0) {
	print("   0: Trotter");
	print("   1: Symplectic-grad-4 Chin-Chen");
	print("   2: Symplectic-grad-4 optimal");
	print("   3: Symplectic-grad-6");
	print("   4: Quadrature 1pt");
	print("   5: Quadrature 2pt");
	print("   6: Quadrature 3pt");
	print("   7: Quadrature 4pt");
	print("   8: Quadrature 5pt");
	print("   9: Quadrature 6pt");
	print(" 10+: Quadrature 7pt+ extremely experimental, try at your own risk:)");

	while (selection < 0 || selection > 20) {
	  std::cout << " Select propagation method (0-20):";
	  std::cout.flush();
	  std::cin >> selection;
	}

	print("Critical time step is", tcrit, "\n");

	std::cout << " Enter time step: ";
	std::cout.flush();
	std::cin >> tstep;
    }

    world.gop.broadcast(selection);
    world.gop.broadcast(tstep);

    int nstep = (velocity==0) ? 100 : int((L - 10 - x0)/velocity/tstep);

      if (world.rank() == 0) {
    print("");
    print(" Wavelet order", k);
    print("     Threshold", thresh);
    print("      Velocity", velocity);
    print("     Time step", tstep);
    print("      No.steps", nstep);
    print("        Method", selection);
      }

    if (selection == 0) {
        for (int step=0; step<nstep; step++) {
            print_info(world, psi,step);
            psi = trotter(world, psi, tstep);
        }
    }
    else if (selection == 1) {
        for (int step=0; step<nstep; step++) {
            print_info(world, psi,step);
            psi = sympgrad4(world, psi, tstep, 0.0, 1.0/72.0); // CC
        }
    }
    else if (selection == 2) {
        for (int step=0; step<nstep; step++) {
            print_info(world, psi,step);
            psi = sympgrad4(world, psi, tstep, -17.0/18000.0, 71.0/4500.0); // Optimal
        }
    }
    else if (selection == 3) {
        for (int step=0; step<nstep; step++) {
            print_info(world, psi,step);
            psi = sympgrad6(world, psi, tstep);
        }
    }
    else {
        np = selection - 3;
        if (world.rank() == 0) print(" No. quad. pt.", np);

	readin(np);
	G = pcomplex_operatorT(MAKE_PROPAGATOR(world, tstep));
        for (int i=0; i<np; ++i)
            Gs.push_back(pcomplex_operatorT(MAKE_PROPAGATOR(world, (1-tc[i])*tstep))); //Generate the free-space FS kernel at quadrature points

        for (int j=0; j<np; ++j)  {    //Generate the free-space FS kernel at SUB quadrature points
            //~ for (int i=0; i<np; ++i) Gss.push_back(pcomplex_operatorT(MAKE_PROPAGATOR(world, (1-tc[i])*tstep*tc[j]))); //[j*np+i] , 1-tc[i] = tc[np-i-1]
            for (int i=0; i<np-j; ++i) Gss.push_back(pcomplex_operatorT(MAKE_PROPAGATOR(world, (1-tc[i])*tstep*tc[j]))); //[j*np+i] , 1-tc[i] = tc[np-i-1]
			for (int i=np-j; i<np; ++i)  Gss.push_back(Gss[(np-1-i)*np+(np-1-j)]); //[j*np+i] , 1-tc[i] = tc[np-i-1] // make use of symmetry
		}

        for (int i=0; i<1+(np>>1); ++i) Gtrs.push_back(pcomplex_operatorT(MAKE_PROPAGATOR(world, (tc[i]-(i?tc[i-1]:0))*tstep/2)));  //Generate the free-space FS kernel for Trotter to construct initial guess
		for (int i=1+(np>>1); i<np; ++i) Gtrs.push_back(Gtrs[np-i]);   // make use of symmetry

        for (int step=0; step<nstep; step++) {
            print_info(world, psi, step);
            psi = q_c(world, np, psi, tstep);
        }
    }

    print("CURRENT TIME",current_time);
    complex_function_1d exact = complex_factory_1d(world).f(psi_exact);
    coord_1d lo, hi;
    lo[0] = -L; hi[0] = L;
    plot_line("psiplot.dat",30001,lo,hi,psi,exact);

    finalize();
    return 0;
}

