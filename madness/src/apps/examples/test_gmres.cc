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
/// A set of routines for testing the GMRES solver.
///
/// Combinations of real/complex MADNESS Vectors/Functions are tested using
/// the child AbstractVectorSpace classes in gmres.h.  Other operators are
/// defined in this file for the tests.
///
/// For each type of number system and class of object, there are three tests:
///   1) x is already equal to the solution (converge on the start)
///   2) x is 0, but A is the identity operator (converge after 1 iteration)
///   3) an arbitrary test that requires more than one iteration.
///
/// NOTE: The final test (complex Function, #3) may sometimes report FAILURE
/// if the final residual is close (but still lower than) the threshold.
/// Tightening up the default tolerance of the MADNESS functions makes this
/// go away.

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <linalg/gmres.h>

using namespace madness;

/// the identity operator for real vectors
class RealVecIdentOp : public Operator<Vector<double, 3> > {
	protected:
		void action(const Vector<double, 3> &invec, Vector<double, 3> &outvec)
			const {

			outvec[0] = 1.0*invec[0] + 0.0*invec[1] + 0.0*invec[2];
			outvec[1] = 0.0*invec[0] + 1.0*invec[1] + 0.0*invec[2];
			outvec[2] = 0.0*invec[0] + 0.0*invec[1] + 1.0*invec[2];
		}
};

/// an arbitrary operator for real vectors
class RealVecLinearOp : public Operator<Vector<double, 3> > {
	protected:
		void action(const Vector<double, 3> &invec, Vector<double, 3> &outvec)
			const {

			outvec[0] = 1.0*invec[0] + 0.0*invec[1] + 1.0*invec[2];
			outvec[1] = 3.0*invec[0] + 1.0*invec[1] + 0.0*invec[2];
			outvec[2] = 2.0*invec[0] + 1.0*invec[1] + 1.0*invec[2];
		}
};

/// the identity operator for complex vectors
class ComplexVecIdentOp : public Operator<Vector<std::complex<double>, 3> > {
	protected:
		void action(const Vector<std::complex<double>, 3> &invec,
			Vector<std::complex<double>, 3> &outvec) const {

			outvec = invec;
		}
};

/// an arbitrary operator for complex vectors
class ComplexVecLinearOp : public Operator<Vector<std::complex<double>, 3> > {
	protected:
		void action(const Vector<std::complex<double>, 3> &invec,
			Vector<std::complex<double>, 3> &outvec) const {

			outvec[0] = std::complex<double>(37.0, 36.0)*invec[0] +
             std::complex<double>(-47.0, 12.0)*invec[1] +
             std::complex<double>(-50.0, 0.0)*invec[2];
			outvec[1] = std::complex<double>(23.0, 12.0)*invec[0] +
             std::complex<double>(0.0, 0.0)*invec[1] +
             std::complex<double>(19.0, 0.0)*invec[2];
			outvec[2] = std::complex<double>(0.0, -3.0)*invec[0] +
             std::complex<double>(17.0, -10.0)*invec[1] +
             std::complex<double>(0.0, 4.0)*invec[2];
		}
};

/// the identity operator for a real function
class RealFuncIdentOp : public Operator<Function<double, 3> > {
	protected:
		void action(const Function<double, 3> &invec,
			Function<double, 3> &outvec) const {

			outvec = copy(invec);
		}
};

/// an arbitrary operator for a real function
/// assumes the function b is never zero in the domain
class RealFuncLinearOp : public Operator<Function<double, 3> > {
	protected:
		const Function<double, 3> &b;

		void action(const Function<double, 3> &invec,
			Function<double, 3> &outvec) const {

			outvec = b * invec;
			invec.compress();
			outvec.truncate();
		}

	public:
		RealFuncLinearOp(const Function<double, 3> &_b) : b(_b) {}
};

/// the identity operator for a complex function
class ComplexFuncIdentOp : public Operator<Function<std::complex<double>, 3> > {
	protected:
		void action(const Function<std::complex<double>, 3> &invec,
			Function<std::complex<double>, 3> &outvec) const {

			outvec = copy(invec);
		}
};

/// an arbitrary operator for a complex function
/// assumes the function b is never zero in the domain
class ComplexFuncLinearOp : public Operator<Function<std::complex<double>, 3> > {
	protected:
		const Function<std::complex<double>, 3> &b;

		void action(const Function<std::complex<double>, 3> &invec,
			Function<std::complex<double>, 3> &outvec) const {

			outvec = b * invec;
			invec.compress();
			outvec.truncate();
		}

	public:
		ComplexFuncLinearOp(const Function<std::complex<double>, 3> &_b) : b(_b) {}
};


/// test functions: true for success, false for failure
const int NTESTS = 12;
bool realvec0();
bool realvec1();
bool realvec2();
bool cplxvec0();
bool cplxvec1();
bool cplxvec2();
bool realfunc0();
bool realfunc1();
bool realfunc2();
bool cplxfunc0();
bool cplxfunc1();
bool cplxfunc2();

/// pointer to the world
World *worldptr;

/// function pointer typedef (for making an array)
typedef bool (*testptr)(void);

/// main routine: execute the tests
int main(int argc, char **argv) {
	int i, passed = 0;

	testptr tests[NTESTS];
	char names[NTESTS][80];
	tests[0] = realvec0;
	sprintf(names[0], "Testing real vectors, 0-step convergence");
	tests[1] = realvec1;
	sprintf(names[1], "Testing real vectors, 1-step convergence");
	tests[2] = realvec2;
	sprintf(names[2], "Testing real vectors, >1-step convergence");
	tests[3] = cplxvec0;
	sprintf(names[3], "Testing complex vectors, 0-step convergence");
	tests[4] = cplxvec1;
	sprintf(names[4], "Testing complex vectors, 1-step convergence");
	tests[5] = cplxvec2;
	sprintf(names[5], "Testing complex vectors, >1-step convergence");
	tests[6] = realfunc0;
	sprintf(names[6], "Testing real functions, 0-step convergence");
	tests[7] = realfunc1;
	sprintf(names[7], "Testing real functions, 1-step convergence");
	tests[8] = realfunc2;
	sprintf(names[8], "Testing real functions, >1-step convergence");
	tests[9] = cplxfunc0;
	sprintf(names[9], "Testing complex functions, 0-step convergence");
	tests[10] = cplxfunc1;
	sprintf(names[10], "Testing complex functions, 1-step convergence");
	tests[11] = cplxfunc2;
	sprintf(names[11], "Testing complex functions, >1-step convergence");

   initialize(argc, argv);
	World world(MPI::COMM_WORLD);
	worldptr = &world;
	startup(world,argc,argv);

	// Function defaults
	FunctionDefaults<3>::set_k(6);
	FunctionDefaults<3>::set_cubic_cell(-2.0, 2.0);
	FunctionDefaults<3>::set_thresh(1.0e-4);
	FunctionDefaults<3>::set_max_refine_level(4);

	// run the tests
	for(i = 0; i < NTESTS; ++i) {
		printf("*** %s:\n", names[i]);
		if(tests[i]()) {
			printf("   --------- PASSED\n\n");
			++passed;
		}
		else
			printf("   --------- FAILED\n\n");
	}

	printf("%d of %d tests passed\n", passed, NTESTS);

   finalize();

	return 0;
}

/// test real vectors, converge on zeroth step
bool realvec0() {
	RealVecIdentOp lo;
	VectorSpace<double, 3> space(*worldptr);
	Vector<double, 3> x, b;
	double resid_thresh = 5.0e-4;
	double update_thresh = 1.0e-10;
	int maxiters = 10;

	b[0] = 1.0;
	b[1] = 2.0;
	b[2] = 3.0;
	x = b;

	GMRES(space, lo, b, x, maxiters, resid_thresh, update_thresh, true);
	// solution is [1,2,3]
	return space.norm(x-b) < 5.0e-4;
}

/// test real vectors, converge after 1 step
bool realvec1() {
	RealVecIdentOp lo;
	VectorSpace<double, 3> space(*worldptr);
	Vector<double, 3> x, b;
	double resid_thresh = 5.0e-4;
	double update_thresh = 1.0e-10;
	int maxiters = 10;

	b[0] = 1.0;
	b[1] = 2.0;
	b[2] = 3.0;
	x = 0.0;

	GMRES(space, lo, b, x, maxiters, resid_thresh, update_thresh, true);
	// solution is [1,2,3]
	return space.norm(x-b) < 5.0e-4;
}

/// test real vectors, converge after >1 steps
bool realvec2() {
	RealVecLinearOp lo;
	VectorSpace<double, 3> space(*worldptr);
	Vector<double, 3> x, b;
	double resid_thresh = 5.0e-4;
	double update_thresh = 1.0e-10;
	int maxiters = 10;

	b[0] = 1.0;
	b[1] = 2.0;
	b[2] = 3.0;
	x = 0.0;

	GMRES(space, lo, b, x, maxiters, resid_thresh, update_thresh, true);
	// solution is [0,2,1]
	b[0] = 0.0; b[1] = 2.0; b[2] = 1.0;
	return space.norm(x-b) < 5.0e-4;
}

/// test complex vectors, converge after 0 steps
bool cplxvec0() {
	ComplexVecIdentOp lo;
	VectorSpace<std::complex<double>, 3> space(*worldptr);
	Vector<std::complex<double>, 3> x, b;
	double resid_thresh = 5.0e-4;
	double update_thresh = 1.0e-10;
	int maxiters = 10;

	b[0] = std::complex<double>(1.0, 2.0);
	b[1] = std::complex<double>(2.0, 1.0);
	b[2] = std::complex<double>(3.0, 0.0);
	x = b;

	GMRES(space, lo, b, x, maxiters, resid_thresh, update_thresh, true);
	// solution is [(1,2), (2,1), (3,0)]
	return space.norm(x-b) < 5.0e-4;
}

/// test complex vectors, converge after 1 step
bool cplxvec1() {
	ComplexVecIdentOp lo;
	VectorSpace<std::complex<double>, 3> space(*worldptr);
	Vector<std::complex<double>, 3> x, b;
	double resid_thresh = 5.0e-4;
	double update_thresh = 1.0e-10;
	int maxiters = 10;

	b[0] = std::complex<double>(1.0, 2.0);
	b[1] = std::complex<double>(2.0, 2.0);
	b[2] = std::complex<double>(3.0, 0.0);
	x = 0.0;

	GMRES(space, lo, b, x, maxiters, resid_thresh, update_thresh, true);
	// solution is [(1,2), (2,2), (3,0)]
	return space.norm(x-b) < 5.0e-4;
}

/// test complex vectors, converge after >1 steps
bool cplxvec2() {
	ComplexVecLinearOp lo;
	VectorSpace<std::complex<double>, 3> space(*worldptr);
	Vector<std::complex<double>, 3> x, b;
	double resid_thresh = 5.0e-4;
	double update_thresh = 1.0e-10;
	int maxiters = 10;

	b[0] = std::complex<double>(1.0, 2.0);
	b[1] = std::complex<double>(2.0, 2.0);
	b[2] = std::complex<double>(3.0, 0.0);
	x = 0.0;

	GMRES(space, lo, b, x, maxiters, resid_thresh, update_thresh, true);
	// solution is in the following vector
	b[0] = std::complex<double>(142537429.0, 13617861.0) / 1013086721.0;
	b[1] = std::complex<double>(114173858.0, 105798015.0) / 1013086721.0;
	b[2] = std::complex<double>(-57303847.0, 132289.0) / 1013086721.0;
	return space.norm(x-b) < 5.0e-4;
}

/// Some functions for testing MADNESS Functions...
static double magfunc(const Vector<double, 3> &pt) {
	return sqrt(pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2]);
}

static std::complex<double> zmagfunc(const Vector<double, 3> &pt) {
	return std::complex<double>(pt[0], sqrt(pt[1]*pt[1] + pt[2]*pt[2]));
}

// inverts a madness function... assumes the function is never 0
template <typename T>
inline static void invert(const Key<3> &key, Tensor<T> &t) {
	UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = 1.0 / *_p0);
}


/// test real functions, converge after 0 steps
bool realfunc0() {
	RealFuncIdentOp lo;
	FunctionSpace<double, 3> space(*worldptr);
	Function<double, 3> x, b;
	double resid_thresh = 1.0e-3;
	double update_thresh = 1.0e-10;
	int maxiters = 10;

	b = FunctionFactory<double, 3>(*worldptr).f(magfunc);
	b.truncate();
	x = copy(b);

	GMRES(space, lo, b, x, maxiters, resid_thresh, update_thresh, true);
	return space.norm(x-b) < 1.0e-3;
}

/// test real functions, converge after 1 step
bool realfunc1() {
	RealFuncIdentOp lo;
	FunctionSpace<double, 3> space(*worldptr);
	Function<double, 3> x, b;
	double resid_thresh = 1.0e-3;
	double update_thresh = 1.0e-10;
	int maxiters = 10;

	b = FunctionFactory<double, 3>(*worldptr).f(magfunc);
	b.truncate();
	x = FunctionFactory<double, 3>(*worldptr); // zero function
	x.truncate();

	GMRES(space, lo, b, x, maxiters, resid_thresh, update_thresh, true);
	return space.norm(x-b) < 1.0e-3;
}

/// test real functions, converge after >1 steps
bool realfunc2() {
	FunctionSpace<double, 3> space(*worldptr);
	Function<double, 3> x, b, mult;
	double resid_thresh = 1.0e-3;
	double update_thresh = 1.0e-10;
	int maxiters = 10;

	b = FunctionFactory<double, 3>(*worldptr).f(magfunc);
	b.truncate();
	x = FunctionFactory<double, 3>(*worldptr); // zero function
	x.truncate();
	mult = copy(b);
	mult.add_scalar(1.0);
	mult.compress();

	RealFuncLinearOp lo(mult);

	GMRES(space, lo, b, x, maxiters, resid_thresh, update_thresh, true);

	// compute the actual solution
	mult.unaryop(&invert<double>);
	b = mult * b;
	b.truncate();
	return space.norm(x-b) < 1.0e-3;
}

/// test complex functions, converge after 0 steps
bool cplxfunc0() {
	ComplexFuncIdentOp lo;
	FunctionSpace<std::complex<double>, 3> space(*worldptr);
	Function<std::complex<double>, 3> x, b;
	double resid_thresh = 1.0e-3;
	double update_thresh = 1.0e-10;
	int maxiters = 10;

	b = FunctionFactory<std::complex<double>, 3>(*worldptr).f(zmagfunc);
	b.truncate();
	x = copy(b);

	GMRES(space, lo, b, x, maxiters, resid_thresh, update_thresh, true);
	return space.norm(x-b) < 1.0e-3;
}

/// test complex functions, converge after 1 step
bool cplxfunc1() {
	ComplexFuncIdentOp lo;
	FunctionSpace<std::complex<double>, 3> space(*worldptr);
	Function<std::complex<double>, 3> x, b;
	double resid_thresh = 1.0e-3;
	double update_thresh = 1.0e-10;
	int maxiters = 10;

	b = FunctionFactory<std::complex<double>, 3>(*worldptr).f(zmagfunc);
	b.truncate();
	x = FunctionFactory<std::complex<double>, 3>(*worldptr); // zero function
	x.truncate();

	GMRES(space, lo, b, x, maxiters, resid_thresh, update_thresh, true);
	return space.norm(x-b) < 1.0e-3;
}

/// test complex functions, converge after >1 steps
bool cplxfunc2() {
	FunctionSpace<std::complex<double>, 3> space(*worldptr);
	Function<std::complex<double>, 3> x, b, mult;
	double resid_thresh = 2.0e-3;
	double update_thresh = 1.0e-10;
	int maxiters = 25;

	b = FunctionFactory<std::complex<double>, 3>(*worldptr).f(zmagfunc);
	b.truncate();
	x = FunctionFactory<std::complex<double>, 3>(*worldptr); // zero function
	x.truncate();
	mult = copy(b);
	mult.add_scalar(std::complex<double>(1.0, 1.0));
	mult.truncate();

	ComplexFuncLinearOp lo(mult);

	GMRES(space, lo, b, x, maxiters, resid_thresh, update_thresh, true);

	// compute the actual solution
	mult.unaryop(&invert<double_complex>);
	b = mult * b;
	b.truncate();
	return space.norm(x-b) < 2.0e-3;
}
