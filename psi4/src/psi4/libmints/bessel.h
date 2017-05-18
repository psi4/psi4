/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2016-2017 Robert A. Shaw.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/* 	class BesselFunction defines the weights and abscissae required for 
   	calculating the modified spherical Bessel function of the first kind,
   	needed for ECP integrals.

   	By Robert A. Shaw 2016 

	REFERENCES:
	R. Flores-Moreno et al., J. Comput. Chem. 27 (2005), 1009
	L.E. McMurchie, E. Davidson, J. Comp. Phys. 44 (1981), 289 
*/

#ifndef BESSEL_FUNCTION_HEAD
#define BESSEL_FUNCTION_HEAD

#include <vector>

namespace psi {

const double SMALL = 1.0E-7;
const int TAYLOR_CUT = 5;

/*! \ingroup MINTS
 *  \class BesselFunction
 *  \brief Computes a modified spherical Bessel function of the first kind.
 *  
 *  Uses pretabulation to calculate the Bessel function up to a given maximum
 *  angular momentum. Values are interpolated using local Taylor series.
 */
class BesselFunction 
{
private:
	/// Maximum angular momentum
	int lMax;
	/// Number of abscissae
	int N; 
	/// Order to which the Bessel series is expanded
	int order;
	
	/// Bessel function values
	double **K;
	/// Coefficients for derivatives of Bessel function
	double *C; 
	
	/**
	  * Pretabulates the Bessel function to a given accuracy.
	  * @param accuracy - the tolerance at which a value is considered converged
	  * @return zero if successful, -1 if not converged
	  */
	int tabulate(const double accuracy);
	
public:
	/// Default constructor. Creates a blank object.
	BesselFunction();
	/// Specified constructor. Will call init with given arguments.
	BesselFunction(int lMax, int N, int order, const double accuracy);
	~BesselFunction();
	
	/**
	  * Initialises and pretabulates the BesselFunction up to the given angular momentum. 
	  * @param lMax - the maximum angular momentum needed
	  * @param N - the maximum number of points to be used in pretabulation, suggested 1600
	  * @param order - the order at which the expansion is cut off, suggested 200
	  * @param accuracy - the tolerance below which a value is considered converged
	  */
	void init(int lMax, int N, int order, const double accuracy);
	
	/**
	  * Calculates the Bessel function values at a given point up to a given angular momentum
	  * @param z - point at which to evaluate
	  * @param maxL - maximum angular momentum needed; must be <= lMax for object
	  * @param values - reference to vector in which to put the values for l = 0 to maxL
	  */
	void calculate(const double z, int maxL, std::vector<double> &values);
};

}
#endif
