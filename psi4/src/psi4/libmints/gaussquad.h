/* 	class GCQuadrature contains the information and methods needed to carry out
   	Gauss-Chebyshev quadrature of the second kind. 

   	Robert A. Shaw 2016	

   	REFERENCES:
   	(Perez92) J.M. Perez-Jorda et al., Comput. Phys. Comm. 70 (1992), 271-284
   	(Perez93) J.M. Perez-Jorda et al., Comput. Phys. Comm. 77 (1993), 46-56
	(Krack98) M. Krack, A.M. Koster, J. Chem. Phys. 108 (1998), 3226 - 3234
	(Flores06) R. Flores-Moreno et al., J. Comput. Chem. 27 (2006), 1009-1019
*/

#ifndef GC_QUAD_HEAD
#define GC_QUAD_HEAD

namespace psi {

#include <vector>

enum GCTYPE {
	ONEPOINT, // Described in Perez92
	TWOPOINT // Described in Perez93
};

/** 
  * \ingroup MINTS
  * \class GCQuadrature
  * \brief Performs adaptive Gauss-Chebyshev quadrature for any given function.
  * 
  * Stores the weights and abscissae for the quadrature, and provides two different methods to integrate on [-1, 1] 
  * Also contains means to transform the region of integration to [0, infinity) and [rmin, rmax]
  */
class GCQuadrature {
private:
	/// Maximum number of points to use in quadrature
	int maxN;
	/// Index of midpoint
	int M;
	
	/// Weights
	std::vector<double> x;
	/// Abscissae
	std::vector<double> w;
	/// Integration value
	double I;
	
	/// Algorithm to be used
	GCTYPE t;
	
   /// Worker function for integration routines, should not be called directly.	
   double sumTerms(double (*f)(double, double*, int), double *p, int limit, int shift, int skip);

public:
	/// Start and endpoints of integration, used for prescreening
	int start, end;
	
	/// Default constructor, creates empty object
	GCQuadrature();
	/// Copy constructor, carbon copies all members
	GCQuadrature(const GCQuadrature &other);
	
	/**
 	  * Intialises the integration grid to the given number of points, and integration type. 
	  * ONEPOINT will choose N = 2^n - 1 closest to the given number of points, whilst
	  * TWOPOINT will choose N= 3*2^n - 1 in the same way.
	  *
  	  * @param points - maximum number of quadrature points to be used
  	  * @param t - the algorithm to be used (ONEPOINT / TWOPOINT)
  	  */
	void initGrid(int points, GCTYPE t);
	
	/**
  	  * Integrates the given function (over [-1, 1] by default) to within the given tolerance. 
  	  * @param f - the function to be integrated
	  * @param params - array of parameters for the function to be integrated
	  * @param tolerance - change below which convergenced is considered to be achieved
  	  * @return true if integration converged, false otherwise
  	  */
	bool integrate(double (*f)(double, double*, int), double *params, const double tolerance);
	
	/**
	  * Transforms the region of integration to [0, inf) using the logarithmic transformation of Krack98
	  */
	void transformZeroInf();
	/**
	  * Transforms region of integration to [rmin, rmax] using the linear transformation from Flores06, assuming 
	  * a Gaussian envelope. rmin/rmax are the distances from the centre of the envelope such that the integrand is effectively zero.
	  * @param z - the exponent of the Gaussian envelope
	  * @param p - the centre of the Gaussian envelope
	  */
	void transformRMinMax(double z, double p);  // Transfromation from [-1, 1] to [rmin, rmax] from Flores06
	
	/// Returns the calculated integral value - must have called integrate first
	double getI() const { return I; }
	
	/// Returns the maximum number of quadrature points
	int getN() const { return maxN; }
	
	/// Returns a reference to the abscissae
	std::vector<double>& getX() { return x; }
};
}

#endif
