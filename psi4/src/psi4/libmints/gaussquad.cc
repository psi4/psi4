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

/* 
	Implements gaussquad.hpp
	Robert A. Shaw 2016
 */

#include "gaussquad.h"
#include <cmath>
#include <iostream>

namespace psi {

// Constructor
GCQuadrature::GCQuadrature() {
	// Currently does nothing
}


GCQuadrature::GCQuadrature(const GCQuadrature &other) {
	maxN = other.maxN;
	M = other.M;
	I = other.I;
	start = other.start;
	end = other.end; 
	t = other.t;
	x = other.x;
	w = other.w;
}

// Initialise the quadrature grid
// As described in both Perez92 and Perez93
void GCQuadrature::initGrid(int points, GCTYPE _t) {
	t = _t;
	
	// Initialise parameters for grid
	int p;
	if (t == ONEPOINT) { // Perez92 one point method
		// We need the number of points to be of the form
		// 2^p - 1 for some power p. 
		p = (int) floor(log(points + 1)/log(2));
		maxN = pow(2, p) - 1;
	} else if (t == TWOPOINT) { // Perez93 two point method
		// Here we need it instead to be of the form
		// 3 * 2^p - 1 for some p.
		p = (int) floor(log((points + 2)/3.0)/log(2));
		maxN = 3*pow(2, p) - 1;
	}
	M = (maxN-1)/2; // Midpoint
	start = 0;
	end = maxN - 1;
	
	// initialise arrays
	x.assign(maxN, 0.0);
	w.assign(maxN, 0.0);
	
	// At the midpoint, M, x[M] = 0 and w[M] = 1
	x[M] = 0.0; w[M] = 1.0;
	// The rest of the abscissae and weights are then given by:
	// z_i = i*Pi / (maxN + 1), s_i = sin(z_i), c_i = cos(z_i)
	// x_i = 1 + 2/(3*pi) * [ (3 + 2*s_i^2) * c_i * s_i - 3z_i]
	// (3(maxN + 1)/16) * w_i = s_i^4
	// We then note that s_(i+1) = c_1 s_i + s_1 c_i
	// and c_(i+1) = c_1 c_i - s_1 s_i
	// with z_(i+1) = z_i + z_1
	// Clearly s_(maxN + 1 - i) = s_i
	// and c_(maxN + 1 - i) = -c_i
	// Therefore x_(maxN + 1 - i) = -x_i
	// and w_(maxN + 1 - i) = w_i
	// Meaning that we only have to calculate half the weights and abscissae
	double z1 = M_PI / ((double)(maxN + 1));
	double c1 = cos(z1); double s1 = sin(z1);
	double zi, si, ci, zi1, si1, ci1; //z_i, s_i, c_i, z_(i+1), s_(i+1), c_(i+1)
	zi1 = z1; si1 = s1; ci1 = c1;
	double o23pi = 2.0 / (3.0 * M_PI); // Convenient
	double s2; //si * si
	for (int n = 0; n < M; n++) {
		// First update zi, si, ci
		zi = zi1;
		si = si1;
		ci = ci1;
		s2 = si * si;
		
		// Now determine the w and x values
		w[maxN - 1 - n] = w[n] = s2 * s2;
		x[n] = 1 + o23pi * ( (3.0 + 2.0 * s2) * ci * si - 3.0*zi );
		x[maxN - 1 - n] = x[n];
		x[n] = -x[n];
		
		// Then update zi1, si1, ci1
		zi1 = zi + z1;
		si1 = c1 * si + s1 * ci;
		ci1 = c1 * ci - s1 * si;
	}
	
	/*std::cout << maxN << " " << M << " " << start << " " << end << "\n";
	for (int q = 0; q < maxN; q++) std::cout << x[q] << " " << w[q] << "\n";*/
}

// Perform the GC integration on the function f
bool GCQuadrature::integrate(std::function<double(double, double*, int)> &f, double *params, const double tolerance) {
	bool converged = false; // 0 for converged, -1 for not converged
	
	if (t == ONEPOINT) {
		// Perez92 Case
		// Integration proceeds in the sequence T_1, T_3, T_7, ..., T_{maxN}
		// where T_m = (3(m+1)/16)I_m
		// by using the fact that T_{2m + 1} = T_{m} + sum_{k = 0}^m w_{2k+1}f(x_{2k+1})
		// The indices in terms of the maxN indices are given by 
		// 2k + 1 = (2k + 1) * M / 2^n = (2k + 1) * p
		// and checking convergence via whether
		// (T_{2m + 1} - 2T_m)^2 <= |T_{2m+1} - 4T_{(m-1)/2}| x tolerance
		double Tn, T2n1, Tn12; // T_n, T_{2n+1} and 4T_{(n-1)/2}
		
		// Initialise values, 
		// Single point integration would use midpoint, M
		Tn = w[M]*f(x[M], params, M);
		Tn12 = 2.0 * Tn;
		
		// Main loop
		int n = 1;
		double dT; // T_{2n+1} - 2T_n
		int ix; // Index needs to be calculated to know which points to use
		int p = (M+1) / 2; // M / 2^n 
		while (n < maxN && !converged) {
			// Compute T_{2n+1}
		 	T2n1 = Tn + sumTerms(f, params, n, p, 2);
			
			// Check convergence
			dT = T2n1 - 2.0*Tn;
			n = 2*n + 1;
			if (dT*dT <= std::fabs(T2n1 - Tn12)*tolerance) {
				converged = true;  
			} else {
				Tn12 = 4.0 * Tn; 
				Tn = T2n1;
				p /= 2; 
			}
		}
		// Finalise the integral
		I = 16.0 * T2n1 / (3.0*(n + 1.0));
		
	} else if (t == TWOPOINT) {
		// Perez93 case
		// We instead proceed along T_2, T_5, T_11, ..., T_{2m + 1} where maxN = 2m+1
		// but also compute T_1, T_3, ..., T_n, T_{2n+1} etc. as before,
		// where m + 1 = 3/2(n+1), so as to get better error control
		// To do this, we use that 
		// T_{2m+1} = T_m + T_n - T_{(n-1)/2} + sum_{i=0}^{(m-2)/3} [w_{6i+1}f(x_{6i+1}) + w_{6i+5}f(x_{6i+5})]
		// along with the same results as before. 
		// The algorithm proceeds by calculating one in the two-point sequence,
		// using an error of |I_{2m+1} - I_m|, then calculates one in the one-point sequence
		// and uses an error of |I_m - I_n|, to check convergence.
		double Tn, Tm, T2n1, T2m1, Tn12;
		
		// Initialise values
		Tn12 = 0.0; 
		Tn = w[M]*f(x[M], params, M);
		int M2 = (maxN - 2)/3; //Index of first point in twopoint sequence
		Tm = w[M2]*f(x[M2], params, M2) + w[maxN - M2 - 1]*f(x[maxN - M2 - 1], params, maxN - M2 - 1);
		int p = (M+1) / 2; // as before
		M2 = (M2 + 1)/2; 
		int ix; 
		int n = 1; int m = 2;
		double error;
		 
		while(m < maxN && !converged) {
			// Propagate the two-point sequence first 
		  T2m1 = Tm + Tn - Tn12 + sumTerms(f, params, (2*m - 1)/3, M2, 3);
			
			// Check convergence
			error = 16.0 * std::fabs(0.5*T2m1 - Tm) / (3.0 * (m + 1)); 
			if (error > tolerance) {
				// Propagate the one-point sequence
			  T2n1 = Tn + sumTerms(f, params, n, p, 2); 
				
				// Check convergence again
				error = 16.0 * std::fabs(2.0*T2m1 - 3.0*T2n1) / (18.0 * (n+1) );
				m = 2 * m + 1;
				n = 2 * n + 1;
				if ( error < tolerance) {
					converged = true; 
				} else {
					Tn12 = Tn;
					Tn = T2n1;
					Tm = T2m1; 
					p /= 2;
					M2 /= 2; 
				}
			} else {
				m = 2 * m + 1;
				converged = true; 
			}
		}
		// Finalise the integral
		I = 16.0 * T2m1 / (3.0 * (m + 1.0));
	}
	
	return converged;
}

// Worker function to do the additional sum terms when going from I_n to I_{2n+1}
double GCQuadrature::sumTerms(std::function<double(double, double*, int)> &f, double *p, int limit, int shift, int skip) {
	double value = 0.0;
	int ix; 
	for (int i = 0; i <= limit; i+=2) {	
		ix = (skip*i + 1)*shift - 1;
		if (ix >= start)
		  value += w[ix] * f(x[ix], p, ix);
		
		ix = maxN - ix - 1; 
		if (ix <= end)
		  value += w[ix] * f(x[ix], p, ix);
	}
	return value;
}

// The GC integrations above are over the interval [-1, 1] and thus need to be transformed
// to the interval[0, infty), or [rmin, rmax]. We do this by the logarithmic transformation from Krack98
// or the linear mapping of Flores06, respectively.  
void GCQuadrature::transformZeroInf() {
	double ln2 = log(2.0);
	double xt;
	
	for (int i = 0; i < maxN; i++) {
		xt = 1.0 - log(1.0-x[i])/ln2;
		w[i] = w[i]/(ln2 * (1.0 - x[i]));
		x[i] = xt;
	}
}

void GCQuadrature::transformRMinMax(double z, double p) {
	double osz = 1.0 / sqrt(z);
	
	// Determine interval
	double rmin = p - 7.0 * osz;
	rmin = rmin > 0 ? rmin : 0.0;
	double rmax = p + 9.0 * osz;
	
	// Find the relative and absolute midpoints 
	double rmid = 0.5*(rmax - rmin); // Midpoint of interval relative to rmin
	double amid = rmid + rmin; // Midpoint of interval
	
	// Transform weights and abscissae by linearly transforming
	// both are scaled by rmid, and the abscissae are translated by amid
	for (int i = 0; i < maxN; i++) {
		x[i] = rmid * x[i] + amid;
		w[i] *= rmid;
	}
}
}
