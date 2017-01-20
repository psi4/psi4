/* Implements ecp.hpp 

   Robert A. Shaw 2016
*/

#include "ecp.h"

#include <cmath>
#include <iostream>
#include <algorithm>

namespace psi {

// GaussianECP constructor and copy constructor
GaussianECP::GaussianECP() : n(0), l(0), a(0), d(0) {}
GaussianECP::GaussianECP(int _n, int _l, double _a, double _d) : n(_n), l(_l), a(_a), d(_d) {}
GaussianECP::GaussianECP(const GaussianECP& other) : n(other.n), l(other.l), a(other.a), d(other.d) {}


// class ECP

ECP::ECP() : N(0), L(-1) {}
ECP::ECP(const double *_center) : N(0), L(-1), center_(_center) {}
ECP::ECP(const ECP &other) {
	gaussians = other.gaussians;
	N = other.N;
	L = other.L;
	center_ = other.center_;
}

void ECP::addPrimitive(int n, int l, double a, double d, bool needSort) {
	GaussianECP newEcp(n, l, a, d);
	gaussians.push_back(newEcp);
	N++;
	L = l > L ? l : L;
	if (needSort) sort();
}

void ECP::sort() {
	std::sort(gaussians.begin(), gaussians.end(),
	 [&] (const GaussianECP& g1, const GaussianECP& g2) {return (g1.l < g2.l);});
}

// Evaluate U_l(r), assuming that gaussians sorted by angular momentum
double ECP::evaluate(double r, int l) {
	double value = 0.0;
	int am = 0;
	double r2 = r*r;
	for (int i = 0; i < N; i++) {
		if (gaussians[i].l == l) 
		  value += pow(r, gaussians[i].n) * gaussians[i].d * exp(-gaussians[i].a * r2);
	} 
	return value; 
}

ECPBasis::ECPBasis() : N(0), maxL(-1) {}

void ECPBasis::addECP(ECP &U) {
	basis.push_back(U);
	N++;
	maxL = U.getL() > maxL ? U.getL() : maxL;
}

ECP& ECPBasis::getECP(int i) { return basis[i]; }

}
