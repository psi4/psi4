/* class ECP contains the contracted expansion of primitive GaussianECPs that define a particular ECP
   class GaussianECP is simply a data structure for the primitive gaussian parameters.  
   class ECPBasis is just a glorified list of all ECPs being used
   These are just skeletons to be replaced at a later date. 

	Robert A. Shaw 2016

	REFERENCES:
	L.E. McMurchie, E.R. Davidson, J. Comput. Phys. 44 (1981), 289 - 301
	R. Flores-Moreno et al., J. Comp. Chem. 27 (2006), 1009 - 1019
 */

#ifndef ECP_HEAD
#define ECP_HEAD

#include <vector>

namespace psi {
// Object describing a Gaussian of angular momentum l of the form
// d r^n e^{-ax^2}
struct GaussianECP {
	int n, l;
	double a, d;
	
	GaussianECP();
	GaussianECP(int n, int l, double a, double d);
	GaussianECP(const GaussianECP& other);
};

class ECP {
private:
	std::vector<GaussianECP> gaussians; // All the primitives in the ECP expansion
	int N, L; // # of Gaussians and maximum angular momentum
	const double *center_;
	
public:
	ECP();
	ECP(const double *_center);
	ECP(const ECP &other);
	
	void addPrimitive(int n, int l, double a, double d, bool needSort = true);
	const double* center() const { return center_; }
	void sort(); // Sort primitives according to angular momentum
	
	// Evaluate U_l(r)
	double evaluate(double r, int l);
  
 	int getL() const { return L; }
	
};

class ECPBasis {
private:
	std::vector<ECP> basis;
	int N, maxL;
	
public:
	ECPBasis();
	
	void addECP(ECP &U);
	ECP& getECP(int i);
	int getMaxL() const { return maxL; }
	int getN() const { return N; }
};

}

#endif
