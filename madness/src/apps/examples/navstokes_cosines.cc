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
  \file navstokes_cosines.cc
  \brief Example Solving the Navier-Stokes equations
  \defgroup examplense Solves a Navier-Stokes equation
  \ingroup examples

  The source is <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/navstokes_cosines.cc >here</a>.

  \par Points of interest
  - convolution with periodic Green's function (Possion/Coulomb kernel and Modified Helmholtz/Bound State Helmholtz/Yukawa kernel)
  - output data for ParaView

  \par Background
  This illustrates the solution of a Navier-Stokes equation for incompressible flows,
  \f{eqnarray*}{
  u_t - u \cdot \nabla u + \nabla p & = &\mu \Delta u + f \\
  \nabla \cdot u & = & 0
  \f}
  where the force and the viscocity  \f$ f \f$ and \f$ \mu \f$ are given in the code.

  \par Implementation

    Step 1.  Calculate the pressure at time \f$ n+1 \f$ explicitly.
    \f[
    \Delta p = \nabla \cdot (f - u_{n} \cdot \nabla u_{n} )
    \f]
    Everything in the RHS is either given or known; thus \f$p\f$ can be obtained by applying a Coulomb operator.

    Step 2.  Calculate the velocity at time n+1.
    \f[
    (\frac{1}{ \delta t \mu } - \Delta) u_{n+1} = \frac {f - \nabla p - u_n \cdot \nabla u_n }{ \mu } + \frac {u_n}{ \delta t \mu }
    \f]
    Again, \f$u_{n+1}\f$ is calculated by applying the BSH operator to the RHS.

    The resulting method is a first order in time scheme and can be extended by Spectral/Krylov deferred corrections to construct higher order methods.
    Particularly, the construction of a second order scheme under this frame is easy and similar to the Crank-Nicolson technique, which is also demonstrated by the example.

    \par Reference
    Jia, J.; Hill, J.; Fann, G. & Harrison, R. J.
    MULTIRESOLUTION FAST METHODS FOR A PERIODIC 3-D NAVIER-STOKES SOLVER
    Proceedings of the Eighth International Symposium on Distributed Computing and Applications to Business, Engineering and Science, Publishing House of Electronics Industry, 2009


*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/vmra.h>
#include <constants.h>

using namespace madness;

typedef Vector<double, 3> coordT3d,coordT;
typedef Vector<double, 1> coordT1d;
typedef Function<double, 3> functionT;
typedef std::vector<functionT> functT;
const double pi=madness::constants::pi;

const double L = 2*pi; //Cell length
const double N = 8.0;

const double mu = 1; // Effective Viscosity
const double deltaT = pi*0.0001; // Size of time step
const int Nts = L/deltaT+10; // Number of time steps
const int k = 10; // Wavelet order (usually precision + 2)
const double pthresh = 1.e-9; // Precision
const double pthresh1 = 1e-10;// * pthresh;
const double uthresh = pthresh; // Precision
const double uthresh1 = pthresh1;


double mytime = 0.0; // Global variable for the current time
// This should be passed in thru the class or app context
const double cc = 0;// L/(deltaT*Nts)/2;

//wrapper, but no longer needed.
struct FunctorInterfaceWrapper : public FunctionFunctorInterface<double,3> {
    double (*f)(const coordT&);

    FunctorInterfaceWrapper(double (*f)(const coordT&)) : f(f) {}

    double operator()(const coordT& x) const {return f(x);}
};

template<typename T, int NDIM> //used to simplifying the mirgrating to new bc. use with caution if you change your bc at runtime.
static Function<T,NDIM> diff(const Function<T,NDIM>& f, int axis) {
	static Vector<std::shared_ptr<Derivative<T,NDIM> >,NDIM> df;
	if (df[axis] == NULL) df[axis] = std::shared_ptr<Derivative<T,NDIM> >(new Derivative<T,NDIM>(f.world(), axis));
	return (*df[axis])(f);
}

//*****************************************************************************
static double init_zero(const coordT3d& r) {
	return 0.0;
}
//*****************************************************************************
//*****************************************************************************

static double uxexact(const coordT3d& r) {
	const double x=r[0]+cc*mytime, y = r[1], z = r[2];
	double t = mytime;

	return cos(.5*t) * sin(x) * sin(x) * (sin(2. * y) * sin(z) * sin(z) - sin(y)
			* sin(y) * sin(2. * z));

}
//*****************************************************************************
//*****************************************************************************
static double uyexact(const coordT3d& r) {
	const double x=r[0]+cc*mytime, y = r[1], z = r[2];
	double t = mytime;

	return cos(.5*t) * sin(y) * sin(y) * (sin(2. * z) * sin(x) * sin(x) - sin(z)
			* sin(z) * sin(2. * x));
}
//*****************************************************************************
//*****************************************************************************
static double uzexact(const coordT3d& r) {
	const double x=r[0]+cc*mytime, y = r[1], z = r[2];
	double t = mytime;

	return cos(.5*t) * sin(z) * sin(z) * (sin(2. * x) * sin(y) * sin(y) - sin(x)
			* sin(x) * sin(2. * y));
}
//*****************************************************************************
//*****************************************************************************
static double fxexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

 return -sin(t / 0.2e1) * pow(sin(x + cc * t), 0.2e1) * (sin(0.2e1 * y) * pow(sin(z), 0.2e1) - pow(sin(y), 0.2e1) * sin(0.2e1 * z)) / 0.2e1 + 0.2e1 * cos(t / 0.2e1) * sin(x + cc * t) * (sin(0.2e1 * y) * pow(sin(z), 0.2e1) - pow(sin(y), 0.2e1) * sin(0.2e1 * z)) * cos(x + cc * t) * cc - cos(t / 0.2e1) * sin(x + cc * t) * sin(y) * cos(z) + 0.2e1 * pow(cos(t / 0.2e1), 0.2e1) * pow(sin(x + cc * t), 0.3e1) * pow(sin(0.2e1 * y) * pow(sin(z), 0.2e1) - pow(sin(y), 0.2e1) * sin(0.2e1 * z), 0.2e1) * cos(x + cc * t) + pow(cos(t / 0.2e1), 0.2e1) * pow(sin(y), 0.2e1) * (sin(0.2e1 * z) * pow(sin(x + cc * t), 0.2e1) - pow(sin(z), 0.2e1) * sin(0.2e1 * x + 0.2e1 * cc * t)) * pow(sin(x + cc * t), 0.2e1) * (0.2e1 * cos(0.2e1 * y) * pow(sin(z), 0.2e1) - 0.2e1 * sin(y) * sin(0.2e1 * z) * cos(y)) + pow(cos(t / 0.2e1), 0.2e1) * pow(sin(z), 0.2e1) * (sin(0.2e1 * x + 0.2e1 * cc * t) * pow(sin(y), 0.2e1) - pow(sin(x + cc * t), 0.2e1) * sin(0.2e1 * y)) * pow(sin(x + cc * t), 0.2e1) * (0.2e1 * sin(0.2e1 * y) * sin(z) * cos(z) - 0.2e1 * pow(sin(y), 0.2e1) * cos(0.2e1 * z)) - mu * (0.2e1 * cos(t / 0.2e1) * pow(cos(x + cc * t), 0.2e1) * (sin(0.2e1 * y) * pow(sin(z), 0.2e1) - pow(sin(y), 0.2e1) * sin(0.2e1 * z)) - 0.2e1 * cos(t / 0.2e1) * pow(sin(x + cc * t), 0.2e1) * (sin(0.2e1 * y) * pow(sin(z), 0.2e1) - pow(sin(y), 0.2e1) * sin(0.2e1 * z)) + cos(t / 0.2e1) * pow(sin(x + cc * t), 0.2e1) * (-0.4e1 * sin(0.2e1 * y) * pow(sin(z), 0.2e1) - 0.2e1 * pow(cos(y), 0.2e1) * sin(0.2e1 * z) + 0.2e1 * pow(sin(y), 0.2e1) * sin(0.2e1 * z)) + cos(t / 0.2e1) * pow(sin(x + cc * t), 0.2e1) * (0.2e1 * sin(0.2e1 * y) * pow(cos(z), 0.2e1) - 0.2e1 * sin(0.2e1 * y) * pow(sin(z), 0.2e1) + 0.4e1 * pow(sin(y), 0.2e1) * sin(0.2e1 * z)));
}
//*****************************************************************************
//*****************************************************************************
static double fyexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

  return -sin(t / 0.2e1) * pow(sin(y), 0.2e1) * (sin(0.2e1 * z) * pow(sin(x + cc * t), 0.2e1) - pow(sin(z), 0.2e1) * sin(0.2e1 * x + 0.2e1 * cc * t)) / 0.2e1 + cos(t / 0.2e1) * pow(sin(y), 0.2e1) * (0.2e1 * sin(0.2e1 * z) * sin(x + cc * t) * cos(x + cc * t) * cc - 0.2e1 * pow(sin(z), 0.2e1) * cos(0.2e1 * x + 0.2e1 * cc * t) * cc) + cos(t / 0.2e1) * cos(x + cc * t) * cos(y) * cos(z) + pow(cos(t / 0.2e1), 0.2e1) * pow(sin(x + cc * t), 0.2e1) * (sin(0.2e1 * y) * pow(sin(z), 0.2e1) - pow(sin(y), 0.2e1) * sin(0.2e1 * z)) * pow(sin(y), 0.2e1) * (0.2e1 * sin(0.2e1 * z) * sin(x + cc * t) * cos(x + cc * t) - 0.2e1 * pow(sin(z), 0.2e1) * cos(0.2e1 * x + 0.2e1 * cc * t)) + 0.2e1 * pow(cos(t / 0.2e1), 0.2e1) * pow(sin(y), 0.3e1) * pow(sin(0.2e1 * z) * pow(sin(x + cc * t), 0.2e1) - pow(sin(z), 0.2e1) * sin(0.2e1 * x + 0.2e1 * cc * t), 0.2e1) * cos(y) + pow(cos(t / 0.2e1), 0.2e1) * pow(sin(z), 0.2e1) * (sin(0.2e1 * x + 0.2e1 * cc * t) * pow(sin(y), 0.2e1) - pow(sin(x + cc * t), 0.2e1) * sin(0.2e1 * y)) * pow(sin(y), 0.2e1) * (0.2e1 * cos(0.2e1 * z) * pow(sin(x + cc * t), 0.2e1) - 0.2e1 * sin(z) * sin(0.2e1 * x + 0.2e1 * cc * t) * cos(z)) - mu * (cos(t / 0.2e1) * pow(sin(y), 0.2e1) * (0.2e1 * sin(0.2e1 * z) * pow(cos(x + cc * t), 0.2e1) - 0.2e1 * sin(0.2e1 * z) * pow(sin(x + cc * t), 0.2e1) + 0.4e1 * pow(sin(z), 0.2e1) * sin(0.2e1 * x + 0.2e1 * cc * t)) + 0.2e1 * cos(t / 0.2e1) * pow(cos(y), 0.2e1) * (sin(0.2e1 * z) * pow(sin(x + cc * t), 0.2e1) - pow(sin(z), 0.2e1) * sin(0.2e1 * x + 0.2e1 * cc * t)) - 0.2e1 * cos(t / 0.2e1) * pow(sin(y), 0.2e1) * (sin(0.2e1 * z) * pow(sin(x + cc * t), 0.2e1) - pow(sin(z), 0.2e1) * sin(0.2e1 * x + 0.2e1 * cc * t)) + cos(t / 0.2e1) * pow(sin(y), 0.2e1) * (-0.4e1 * sin(0.2e1 * z) * pow(sin(x + cc * t), 0.2e1) - 0.2e1 * pow(cos(z), 0.2e1) * sin(0.2e1 * x + 0.2e1 * cc * t) + 0.2e1 * pow(sin(z), 0.2e1) * sin(0.2e1 * x + 0.2e1 * cc * t)));

}
//*****************************************************************************
//*****************************************************************************
static double fzexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

  return -sin(t / 0.2e1) * pow(sin(z), 0.2e1) * (sin(0.2e1 * x + 0.2e1 * cc * t) * pow(sin(y), 0.2e1) - pow(sin(x + cc * t), 0.2e1) * sin(0.2e1 * y)) / 0.2e1 + cos(t / 0.2e1) * pow(sin(z), 0.2e1) * (0.2e1 * cos(0.2e1 * x + 0.2e1 * cc * t) * cc * pow(sin(y), 0.2e1) - 0.2e1 * sin(x + cc * t) * sin(0.2e1 * y) * cos(x + cc * t) * cc) - cos(t / 0.2e1) * cos(x + cc * t) * sin(y) * sin(z) + pow(cos(t / 0.2e1), 0.2e1) * pow(sin(x + cc * t), 0.2e1) * (sin(0.2e1 * y) * pow(sin(z), 0.2e1) - pow(sin(y), 0.2e1) * sin(0.2e1 * z)) * pow(sin(z), 0.2e1) * (0.2e1 * cos(0.2e1 * x + 0.2e1 * cc * t) * pow(sin(y), 0.2e1) - 0.2e1 * sin(x + cc * t) * sin(0.2e1 * y) * cos(x + cc * t)) + pow(cos(t / 0.2e1), 0.2e1) * pow(sin(y), 0.2e1) * (sin(0.2e1 * z) * pow(sin(x + cc * t), 0.2e1) - pow(sin(z), 0.2e1) * sin(0.2e1 * x + 0.2e1 * cc * t)) * pow(sin(z), 0.2e1) * (0.2e1 * sin(0.2e1 * x + 0.2e1 * cc * t) * sin(y) * cos(y) - 0.2e1 * pow(sin(x + cc * t), 0.2e1) * cos(0.2e1 * y)) + 0.2e1 * pow(cos(t / 0.2e1), 0.2e1) * pow(sin(z), 0.3e1) * pow(sin(0.2e1 * x + 0.2e1 * cc * t) * pow(sin(y), 0.2e1) - pow(sin(x + cc * t), 0.2e1) * sin(0.2e1 * y), 0.2e1) * cos(z) - mu * (cos(t / 0.2e1) * pow(sin(z), 0.2e1) * (-0.4e1 * sin(0.2e1 * x + 0.2e1 * cc * t) * pow(sin(y), 0.2e1) - 0.2e1 * pow(cos(x + cc * t), 0.2e1) * sin(0.2e1 * y) + 0.2e1 * pow(sin(x + cc * t), 0.2e1) * sin(0.2e1 * y)) + cos(t / 0.2e1) * pow(sin(z), 0.2e1) * (0.2e1 * sin(0.2e1 * x + 0.2e1 * cc * t) * pow(cos(y), 0.2e1) - 0.2e1 * sin(0.2e1 * x + 0.2e1 * cc * t) * pow(sin(y), 0.2e1) + 0.4e1 * pow(sin(x + cc * t), 0.2e1) * sin(0.2e1 * y)) + 0.2e1 * cos(t / 0.2e1) * pow(cos(z), 0.2e1) * (sin(0.2e1 * x + 0.2e1 * cc * t) * pow(sin(y), 0.2e1) - pow(sin(x + cc * t), 0.2e1) * sin(0.2e1 * y)) - 0.2e1 * cos(t / 0.2e1) * pow(sin(z), 0.2e1) * (sin(0.2e1 * x + 0.2e1 * cc * t) * pow(sin(y), 0.2e1) - pow(sin(x + cc * t), 0.2e1) * sin(0.2e1 * y)));
}

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
static double pexact(const coordT3d& r) {
	const double x=r[0]+cc*mytime, y = r[1], z = r[2];
	double t = mytime;

	return cos(.5*t) * cos(x) * sin(y) * cos(z);
}

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

// Compute the divergence
inline functionT div(const functT& uint) {
	return diff(uint[0], 0) + diff(uint[1], 1) + diff(uint[2], 2);
}

// Compute the Laplacian
inline functionT lap(const functionT& uint) {
	return diff(diff(uint, 0),0) + diff(diff(uint,1), 1) + diff(diff(uint,2), 2);
}

World *pworld;
#define myfun std::vector< Function<T,NDIM> >

// Compute the advection of \c uu and store it in \c advu for return
template<typename T, int NDIM> void adv(const myfun& uu, myfun& advu) {
	for (int i=0; i < 3; ++i)  advu[i] = diff(uu[0]*uu[i],0) + diff(uu[1]*uu[i],1) + diff(uu[2]*uu[i],2);
}

template<typename T, int NDIM> inline myfun operator-(const myfun& l, const myfun& r) { return sub(*pworld, l, r); }

void testNavierStokes(int argc, char**argv) {
	initialize(argc, argv);
	try {
	World world(MPI::COMM_WORLD);

	pworld = &world;
	startup(world, argc, argv);

	// Function defaults
	FunctionDefaults<3>::set_k(k);
	FunctionDefaults<3>::set_cubic_cell(0.0, L);
	FunctionDefaults<3>::set_thresh(pthresh);
	FunctionDefaults<3>::set_bc(BC_PERIODIC);

	// construct the periodic Coulomb operator for later use
	SeparatedConvolution<double, 3> op = CoulombOperator (world, pthresh1, pthresh1);

	// construct the periodic BSH operator for later use
	double const dum = 1 / deltaT / mu;
	SeparatedConvolution<double, 3> op1 = BSHOperator<3>(world,	sqrt(dum), uthresh1, uthresh1);

	
	// construct the periodic BSH operator for later use
	Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
//	SeparatedConvolution<double, 3> op11 = PeriodicBSHOp<double, 3> (world, sqrt(dum), k, uthresh1, uthresh1, cellsize);

	// Initialize the old solution and print out to vts files
	mytime = 0.0;

	functT u(3);
	functT rhs(3);
	functT f(3);
	u[0] = FunctionFactory<double, 3> (world).f(uxexact  ) .truncate_on_project();
	u[1] = FunctionFactory<double, 3> (world).f(uyexact  ) .truncate_on_project();
	u[2] = FunctionFactory<double, 3> (world).f(uzexact  ) .truncate_on_project();


	print("col error",op(lap(u[0])).scale(-1. / (4. * pi)).err(uxexact));
	print("bsh error", op1(dum*u[0] - lap(u[0])).err(uxexact));


	Function<double, 3> divu = div(u); // computer the divergence
	double divun=divu.norm2();
	int dd=divu.max_depth();
	if (world.rank()==0) print("initial div, depth:", divun, dd); //print out some information

	for (int t =  0; t < Nts; t++) {
		mytime = deltaT*(t+1);
		//if (world.rank()==0) print("current time: ", mytime);
		// Step 1.  Calculate the pressure at time t+1.
		//            Laplace p = div (f-u grad u)
		f[0] = FunctionFactory<double, 3> (world).f(fxexact).truncate_on_project();
		f[1] = FunctionFactory<double, 3> (world).f(fyexact).truncate_on_project();
		f[2] = FunctionFactory<double, 3> (world).f(fzexact).truncate_on_project();


		adv(u, rhs);
		//adv works like: for (int i=0; i < 3; ++i)  rhs[i] = u[0]*diff(u[i],0) + u[1]*diff(u[i],1) + u[2]*diff(u[i],2);

		functionT divf = div(f-rhs);

		Function<double,3> p = op(divf); // apply the Coulomb operator to compute the pressure \c p
		p.scale(-1. / (4. * pi));

		// Step 2.  Calculate the velocity at time t+1.
		//            (1/(deltaT mu) - Laplace) u_t+1 = (f - grad p - u grad u)/mu + u_t/(deltaT mu)

		// do the following calculation
		//~ rhs[0] = (f[0] - diff(p, 0) -rhs[0])*(1. / mu) + u[0]*dum;
		//~ rhs[1] = (f[1] - diff(p, 1) -rhs[1])*(1. / mu) + u[1]*dum;
		//~ rhs[2] = (f[2] - diff(p, 2) -rhs[2])*(1. / mu) + u[2]*dum;
		f[0] -= diff(p,0);
		f[1] -= diff(p,1);
		f[2] -= diff(p,2);
		gaxpy(world, 1, rhs, -1, f);
		gaxpy(world, -1./mu, rhs, dum, u);

		functT ue = apply(world, op1, rhs); // apply the BSH operator to update the velocity \c ue

		//u = ue;  // use this line for first order/mixed Euler's method

		gaxpy(world,-1,u,2,ue); ++t; mytime += deltaT; //for (int i=0; i < 3; ++i) u[i] = 2.0*ue[i] - u[i];// += (mu*lap(ue[i])-(ue[0]*diff(ue[i],0) + ue[1]*diff(ue[i],1) + ue[2]*diff(ue[i],2)) - diff(p,i) + f[i])*(2*deltaT);
		//use the above line for second order/Crank-Nicolson like scheme
		//note that in this case, the time-step size is instead 2*deltaT


		if ( (t%10)==0 && world.rank()==0) { // output the current status to vts files every 10 steps
                   char filename[100];
                   sprintf(filename, "data-%02d.vts", t);
                   Vector<double, 3> plotlo, plothi;
                   Vector<long, 3> npts;
                   for(int i = 0; i < 3; ++i) {
                           plotlo[i] = 0.0;
                           plothi[i] = 1.0;
                           npts[i] = 21;
                   }
                   plotvtk_begin(world, filename, plotlo, plothi, npts);
                   plotvtk_data(u[0], "u", world, filename, plotlo, plothi, npts);
                   plotvtk_data(u[1], "v", world, filename, plotlo, plothi, npts);
                   plotvtk_data(u[2], "w", world, filename, plotlo, plothi, npts);
                   plotvtk_data(p   , "p", world, filename, plotlo, plothi, npts);
                   plotvtk_end<3>(world, filename);
		}

		{
		//~ Function<double, 3> du = FunctionFactory<double, 3> (world).f(uxexact).truncate_on_project();
		//~ du -= u[0];
		//~ Function<double, 3> dv = FunctionFactory<double, 3> (world).f(uyexact).truncate_on_project();
		//~ dv -= u[1];
		//~ Function<double, 3> dw = FunctionFactory<double, 3> (world).f(uzexact).truncate_on_project();
		//~ dw -= u[2];

		double  a=div(u).norm2(), b=u[0].err(uxexact), c=u[1].err(uyexact),d=u[2].err(uzexact);
		if (world.rank()==0)  print(t+1, mytime, a,b,c,d); // print out some information
		}
	}

//	RMI::end();
	} catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    } catch (char* s) {
        print(s);
        error("caught a c-string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    };
    finalize();
}
//*****************************************************************************

//*****************************************************************************
int main(int argc, char**argv) {
	testNavierStokes(argc, argv);
	return 0;
}
//*****************************************************************************

