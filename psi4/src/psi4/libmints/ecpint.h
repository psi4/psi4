/* 	class ECPInt calculates one-body ECP integrals
	class AngularIntegral calculates and stores the angular integrals needed for the ECP integration
 	class RadialIntegral abstracts the calculation of the radial integrals needed for the ECP integration, 
	such that if a different approach was desired later, this could be done with minimal alterations to ECPInt.

   	Robert A. Shaw 2016	

   	REFERENCES:
	(Flores06) R. Flores-Moreno et al., J. Comput. Chem. 27 (2006), 1009-1019
	(MM81) L. E. McMurchie and E. R. Davidson, J. Comp. Phys. 44 (1981), 289 - 301
*/

#ifndef ECPINT_HEAD
#define ECPINT_HEAD

#include <vector>
#include "psi4/libmints/multiarr.h"
#include "psi4/libmints/gaussquad.h"
#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/ecp.h"
#include "psi4/libmints/bessel.h"

namespace psi {

    class BasisSet;
    class GaussianShell;
    class IntegralFactory;
    class SphericalTransform;

/** 
  * Calculates real spherical harmonics S_lm(theta, phi) for all l, m up to lmax
  * @param lmax - the maximum angular momentum needed
  * @param x - cos(theta), where theta is the polar angle in spherical coordinates
  * @param phi - the azimuth angle in spherical coordinates
  * @return a matrix S(l, l+m) of the spherical harmonic values
  */
static TwoIndex<double> realSphericalHarmonics(int lmax, double x, double phi);  

/**
  * \ingroup MINTS
  * \class ShellPairData
  * \brief Stores the (shifted) angular momenta, number of cartesians in a shell pair, and shifted centers
  */
struct ShellPairData {
	int LA, LB, maxLBasis;
	int ncartA, ncartB;
	double A[3], B[3];
	double A2, Am, B2, Bm, RAB2, RABm;
};

/** 
  * \ingroup MINTS
  * \class AngularIntegral
  * \brief Calculates and stores the angular integrals needed for ECP integration.
  * 
  * This should not usually be created directly, it is instead owned by an ECPIntegral object,
  * so that integrals can be performed over multiple ECP centers without duplicating work.
  */
class AngularIntegral 
{
private: 
	/// Maximum angular momentum of orbital basis and ECP basis, respectively
	int LB, LE; 
	/// Limits for the w-integral indices, and angular momentum indices
	int wDim, maxL; 
	
	/// Stores the type 1 angular integrals
	FiveIndex<double> W; 
	/// Stores the type 2 angular integrals
	SevenIndex<double> omega; 
	
	/// Worker functions for calculating terms in the USP to spherical transformation coefficients
	double calcG(int l, int m) const;
	double calcH1(int i, int j, int l, int m) const;
	double calcH2(int i, int j, int k, int m) const;
	
	/**
	  * Calculates all possible USP to spherical transformation coefficients for a given angular momentum
	  * @param lam - the angular momentum
	  * @param mu - the subshell
	  * @return ThreeIndex of the values U_lam,mu(k, l, m)
	  */
	ThreeIndex<double> uklm(int lam, int mu) const;
	/**
	  * Calculates the polynomial integrals, int (x^i y^j z^k dOm) where dOm is the solid angle differential
	  * @param maxI - the maximum power, i, to determine
	  * @return ThreeIndex of polynomials P(i, j, k), where strictly i >= j >= k
	  */
	ThreeIndex<double> Pijk(int maxI) const; 
	
	/**
	  * Builds the USP to spherical transformation coefficients for use in calculating the type 1 and 2 integrals
	  * @return FiveIndex of the coefficients U(lam, lam+mu, k, l, m)
	  */
	FiveIndex<double> makeU();
	/**
	  * Builds the type 1 angular integrals
	  * @param U - the USP to spherical transformation coefficients
	  */
	void makeW(FiveIndex<double> &U);
	/** 
	  * Builds the type 2 angular integrals
	  * @param U - the USP to spherical transformation coefficients
	  */
	void makeOmega(FiveIndex<double> &U);
	
public:
	
	/// Default constructor creates empty object
	AngularIntegral(); 
	/// Specified constructor calls init with given arguments
	AngularIntegral(int LB, int LE); 
	/**
	  * Initialises the object, must be called before anything else if default constructor was used.
	  * @param LB - the maximum angular momentum of the orbital basis
	  * @param LE - the maximum angular momentum of the ECP basis
	  */
	void init(int LB, int LE);
	/**
	  * Computes the type 1 and 2 angular integrals
	  */
	void compute();
	
	/// TODO: Clears the W and omega arrays
	void clear();
	
	/**
	  * Returns the type 1 angular integral W(k, l, m, lam, mu)
	  * @param k - x index
	  * @param l - y index
	  * @param m - z index
	  * @param lam - angular momentum
	  * @param mu - subshell
	  * @return value of type 1 angular integral
	  */
	double getIntegral(int k, int l, int m, int lam, int mu) const; 
	/**
  	  * Returns the type 2 angular integral Omega(k, l, m, lam, mu, rho, sigma)
  	  * @param k - x index
  	  * @param l - y index
  	  * @param m - z index
  	  * @param lam - angular momentum of current ECP shell
  	  * @param mu - subshell of lam
	  * @param rho - angular momentum of current basis shell
	  * @param sigma - subshell of rho
  	  * @return value of type 2 angular integral
  	  */
	double getIntegral(int k, int l, int m, int lam, int mu, int rho, int sigma) const;

 	/// is W(k, l, m, lam, mu) zero to within a given tolerance?
	bool isZero(int k, int l, int m, int lam, int mu, double tolerance) const;
	/// is Omega(k, l, m, lam, mu, rho, sigma) zero to within a given tolerance?
	bool isZero(int k, int l, int m, int lam, int mu, int rho, int sigma, double tolerance) const;	
};


/** 
  * \ingroup MINTS
  * \class RadialIntegral
  * \brief Abstracts the calculation of radial integrals for ECP integration.
  * 
  * This should not be used directly, and is owned by ECPIntegral.
  * It provides the interface to the adaptive quadrature algorithms used to calculate the type 1 and 2 radial integrals. 
  */
class RadialIntegral
{
private:
	/// The larger integration grid for type 1 integrals, and for when the smaller grid fails for type 2 integrals 
	GCQuadrature bigGrid;
	/// The smaller integration grid, default for the type 2 integrals
    GCQuadrature smallGrid;
	/// Modified spherical Bessel function of the first kind
	BesselFunction bessie;
	
	/// Matrices of parameters needed in both type 1 and 2 integrations
	TwoIndex<double> p, P, P2, K;
	
	/// Tolerance for change below which an integral is considered converged
	double tolerance;
	
	/// This integrand simply returns the pretabulated integrand values stored in p given an index ix
	static double integrand(double r, double *p, int ix);

	/**
	  * Builds a matrix of Bessel at the given points up to the given maximum angular momentum. 
	  * @param r - vector of points to evaluate at
	  * @param nr - number of points in r (for convenience)
	  * @param maxL - the maximum angular momentum needed
	  * @param values - TwoIndex<double> to store the values in
	  * @param weight - factor to weight r by (defaults to 1)
	  */
	void buildBessel(std::vector<double> &r, int nr, int maxL, TwoIndex<double> &values, double weight = 1.0);
	
	double calcKij(double Na, double Nb, double zeta_a, double zeta_b, double R2) const;
	
	/**
	  * Tabulate r^{N+2} times the ECP for all quadrature points. 
	  * @param U - the ECP to be pretabulated
	  * @param l - the angular momentum shell of the ECP to be used
	  * @param N - the power of r to weight the ECP by
	  * @param grid - the quadrature grid to be used
	  * @param Utab - the array to put the values into.
	  */
    void buildU(const GaussianShell &U, int l, int N, GCQuadrature &grid, double *Utab);
	
	/**
	  * Tabulate the F function values for the default mode of calculating type 2 integrals.
	  * @param shell - the shell of orbital basis functions to tabulate over
	  * @param lstart - the lowest angular momentum needed
	  * @param lend - the maximum angular momentum needed
	  * @param r - quadrature grid points
	  * @param nr - the number of grid points in r
	  * @param start - the grid point to start at
	  * @param end - the grid point to stop at
	  * @param F - the matrix to put the values in
	  */
	void buildF(const GaussianShell &shell, double A, int lstart, int lend, std::vector<double> &r, int nr, int start, int end, TwoIndex<double> &F);
	
	/**
	  * Performs the integration given the pretabulated integrand values. 
	  * @param maxL - the maximum angular momentum needed
	  * @param gridSize - the number of quadrature points
	  * @param intValues - the TwoIndex<double> of pretabulated integrand values for each angular momentum needed
	  * @param grid - the quadrature grid
	  * @param values - the vector to put the resulting integrals into
	  * @param offset - the angular momentum to start at (defaults to 0)
	  * @param skip - the steps of angular momentum to go up in (defaults to 1)
	  */
	int integrate(int maxL, int gridSize, TwoIndex<double> &intValues, GCQuadrature &grid, std::vector<double> &values, int offset = 0, int skip = 1);

public:
	/// Default constructor creates an empty object
	RadialIntegral();
	
	/**
	  * Initialises the object, in turn intialising the quadrature grids and BesselFunction
	  * @param maxL - the maximum angular momentum of integral needed
	  * @param tol - the tolerance for convergence of integrals (defaults to 1e-15)
	  * @param small - the maximum number of quadrature points for the small integration grid (default 256, minimum recommended)
	  * @param large - the maximum number of quadrature points for the large integration grid (default 1024, minimum recommended)
	  */
	void init(int maxL, double tol = 1e-15, int small = 256, int large = 1024);
	
	/**
	  * Given two GaussianShells, builds the parameters needed by both kind of integral. 
	  * @param shellA - the first GaussianShell
	  * @param shellB - the second GaussianShell
	  * @param A - position vector (relative to the ECP center) of shell A
	  * @param B - position vector (relative to the ECP center) of shell B
	  */
	void buildParameters(const GaussianShell &shellA, const GaussianShell &shellB, ShellPairData &data);
	
	/**
	  * Calculates all type 1 radial integrals over two Gaussian shells up to the given maximum angular momentum.
	  * @param maxL - the maximum angular momentum
	  * @param N - the power of r that the integrand is weighted by
	  * @param offset - the starting angular momentum
	  * @param U - the ECP to be integrated over
	  * @param shellA - the first GaussianShell
	  * @param shellB - the second GaussianShell
	  * @param A - position vector (relative to the ECP center) of shell A
	  * @param B - position vector (relative to the ECP center) of shell B
	  * @param values - the matrix to return the integrals in
	  */
    void type1(int maxL, int N, int offset, const GaussianShell &U, const GaussianShell &shellA, const GaussianShell &shellB, ShellPairData &data, TwoIndex<double> &values);
	
    /**
      * Calculates all type 2 radial integrals over two Gaussian shells for the given ECP angular momentum l
      * @param lam - the ECP shell angular momentum to be calculated over
	  * @param l1start - the angular momentum to start on for the first shell
	  * @param l1end - the angular momentum to stop at for the first shell
	  * @param l2start - the angular momentum to start on for the second shell
	  * @param l2end - the angular momentum to stop at for the second shell
      * @param N - the power of r that the integrand is weighted by
      * @param U - the ECP to be integrated over
      * @param shellA - the first GaussianShell
      * @param shellB - the second GaussianShell
      * @param A - position vector (relative to the ECP center) of shell A
      * @param B - position vector (relative to the ECP center) of shell B
      * @param values - the matrix to return the integrals in
      */
    void type2(int lam, int l1start, int l1end, int l2start, int l2end, int N, const GaussianShell &U, const GaussianShell &shellA, const GaussianShell &shellB, ShellPairData &data, TwoIndex<double> &values);
};


/** 
  * \ingroup MINTS
  * \class ECPIntegral
  * \brief Calculates ECP integrals. 
  * 
  * Given an ECP basis, and orbital bases, this will calculate the ECP integrals over all ECP centers. 
  * TODO: Implement derivatives (identical to normal integrals, but with shifted angular momenta)
  */
class ECPInt : public OneBodyAOInt
{
private:
	/// The interface to the radial integral calculation
	RadialIntegral radInts;
	/// The angular integrals, which can be reused over all ECP centers
	AngularIntegral angInts;
	/// The ECP basis
    std::shared_ptr<BasisSet> basis;
	
	/// Worker functions for calculating binomial expansion coefficients
	double calcC(int a, int m, double A) const;
	void makeC(FiveIndex<double> &C, int L, double *A);
	
	/// Calculates the type 1 integrals for the given ECP center over the given shell pair
    void type1(const GaussianShell& U, const GaussianShell &shellA, const GaussianShell &shellB, ShellPairData &data, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &values);
	/// Calculates the type 2 integrals for the given ECP center over the given shell pair
    void type2(int l, const GaussianShell& U, const GaussianShell &shellA, const GaussianShell &shellB, ShellPairData &data, FiveIndex<double> &CA, FiveIndex<double> &CB, ThreeIndex<double> &values);
	
	/// Overridden shell-pair integral calculation over all ECP centers
	void compute_pair(const GaussianShell &shellA, const GaussianShell &shellB);
	
	/// Computes the overall ECP integrals over the given ECP center and shell pair
    void compute_shell_pair(const GaussianShell &U, const GaussianShell &shellA, const GaussianShell &shellB, TwoIndex<double> &values, int shiftA = 0, int shiftB = 0);
	
public:
	/**
	  * Sets the reference to the ECP basis and initialises the radial and angular integrals
	  * @param basis - reference to the ECP basis set
	  * @paramm maxLB - the maximum angular momentum in the orbital basis
	  */
    ECPInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv = 0);
    virtual ~ECPInt();
	
};
    
class ECPSOInt : public OneBodySOInt
{
    int natom_;
public:
    ECPSOInt(const std::shared_ptr<OneBodyAOInt>& , const std::shared_ptr<IntegralFactory> &);
    ECPSOInt(const std::shared_ptr<OneBodyAOInt>& , const IntegralFactory*);
};


}
#endif
