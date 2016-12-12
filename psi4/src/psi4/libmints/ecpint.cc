/* Implements ecpint.hpp */

#include "ecpint.h"
#include "gshell.h"
#include <iostream>
#include <functional>
#include <cmath>

// Compute all the real spherical harmonics Slm(theta, phi) for l,m up to lmax
// x = cos (theta)
static TwoIndex<double> realSphericalHarmonics(int lmax, double x, double phi, std::vector<double> &fac, std::vector<double> &dfac){
	TwoIndex<double> rshValues(lmax+1, 2*lmax+1, 0.0);

	if (lmax > 0) {
		// First calculate the associated Legendre polynomials, Plm(cos theta), using the recursion relation
		// (l-m)Plm = x(2l - 1)P{l-1}m - (l+m-1)P{l-2}m
		// along with the zeroth order term
		// Pmm = (-1)^m (2m-1)!!(1-x^2)^{m/2}
		double x2 = x * x;
		double Plm[lmax+1][lmax+1]; 
		// First get all Pmm terms
		Plm[0][0] = 1.0;
		double sox2 = sqrt(1.0 - x2);
		double ox2m = 1.0;
		for (int m = 1; m <= lmax; m++) {
			ox2m *= -sox2;
			Plm[m][m] = ox2m * dfac[2*m-1];
		}
		
		// Then increment l for each m
		Plm[1][0] = x;
		Plm[0][1] = 0.0;
		for (int l = 2; l <= lmax; l++) {
			ox2m = x * (2*l - 1);
			for (int m = 0; m < l; m++) {
				Plm[l][m] = ox2m * Plm[l-1][m] - (l + m - 1)*Plm[l-2][m];
				Plm[l][m] /= ((double) (l -m));
			}
			Plm[l-1][l] = 0.0;
		}
		
		// Now we compute the spherical harmonics via
		// Slm(theta, phi) = Clm * Plm(cos(theta)) * cos(m * phi), m > 0
		// Sl{-m}(theta, phi) = Clm * Plm(cos(theta)) * sin(m * phi)
		// Sl0(theta, phi) = sqrt(2) * Cl0 * Pl0(cos(theta))
		// where Clm^2 = (2l + 1)*(l - m)! / (8*pi * (l+m)!)
		double osq4pi = 1.0 / sqrt(4.0 * M_PI); 
		int sign;
		for (int l = 0; l <= lmax; l++) {
			rshValues(l, l) = osq4pi * sqrt(2.0 * l + 1.0) * Plm[l][0];
			sign = -1;
			for (int m = 1; m <= l; m++) {
				ox2m = (2.0 * l + 1.0) * fac[l-m] / fac[l+m];
				ox2m = sign * osq4pi * sqrt(2.0 * ox2m) * Plm[l][m];
				rshValues(l, l+m) = ox2m * cos(m * phi);
				rshValues(l, l-m) = ox2m * sin(m * phi);
				sign *= -1;
			}
		}
		
	} else {
		rshValues(0, 0) = 1.0 / sqrt(4.0 * M_PI);
	}
		
	return rshValues;
}

double AngularIntegral::calcG(int l, int m, std::vector<double> &fac) const {
	double value = 0.0;
	double value1 = pow(2.0, l) * fac[l];
	value1 = 1.0 / value1; 
	double value2 = (2.0 * l + 1) * fac[l - m] / (2.0 * M_PI * fac[l + m]);
	value2 = sqrt(value2); 
	value = value1 * value2;
	return value;
} 

double AngularIntegral::calcH1(int i, int j, int l, int m, std::vector<double> &fac) const {
	double value = 0.0; 

	value = fac[l]/(fac[j]*fac[l - i]*fac[i-j]);
	value *= (1 - 2*(i%2)) * fac[2*(l - i)] / (fac[l - m - 2*i]);

	return value;
}

double AngularIntegral::calcH2(int i, int j, int k, int m, std::vector<double> &fac) const {
	double value = 0.0; 
	int ki2 = k - 2*i;
	if ( m >= ki2 && ki2 >= 0 ) {
		value = fac[j]*fac[m]/(fac[i] * fac[j-i] * fac[ki2] * fac[m-ki2]);
		int p = (m - k + 2*i)/2;
		value *= (1.0 - 2.0*(p%2));
	}
	return value;
}


ThreeIndex<double> AngularIntegral::uklm(int lam, int mu, std::vector<double> &fac) const {
	ThreeIndex<double> values(lam+1, lam+1, 2);
	 
  	double or2 = 1.0/sqrt(2.0);
  	double u = 0.0;
	double um = 0.0;
	double g = calcG(lam, mu, fac);

  	double u1, h1, h2;
  	int j;
  	for (int k = 0; k <= lam; k++) {
  	  for (int l = 0; l <= lam - k; l++) {
		u = um = 0.0;
	  	j = k + l - mu;
		if (j % 2 == 0 && j > -1) { 
			u1 = 0.0;
			j/=2;
			for (int i = j; i <= (lam - mu)/2; i++) u1 += calcH1(i, j, lam, mu, fac);
			
			u = g * u1;
			u1 = 0;
			for (int i = 0; i <= j; i++) u1 += calcH2(i, j, k, mu, fac);
			u *= u1;
			um = u;
			
			j = l % 2;
			u *= (1 - j);
			um *= j;
			if (mu == 0) {
				u *= or2;
				um = u;
			} 
		}
		values(k, l, 0) = u;
		values(k, l, 1) = um;
	  }
	}
	return values;						
}


ThreeIndex<double> AngularIntegral::Pijk(int maxI) const {
	int dim = maxI+1;
	ThreeIndex<double> values(dim, dim, dim);
	double pi4 = 4.0*M_PI;
	
	values(0, 0, 0) = pi4;
	for (int i = 1; i <= maxI; i++) {
		values(i, 0, 0) = pi4 / ((double) (2*i+1));
		
		for (int j = 1; j <= i; j++) {
			values(i, j, 0) = values(i, j-1, 0) * (2.0*j - 1.0) / (2.0 * ((double)(i + j)) + 1.0);
			
			for (int k = 1; k <= j; k++)
				values(i, j, k) = values(i, j, k-1) * (2.0*k - 1.0) / (2.0 * ((double)(i + j + k)) + 1.0);
			
		}
	}
	return values;
}

FiveIndex<double> AngularIntegral::makeU(std::vector<double> &fac) {
	int dim = maxL + 1;

	FiveIndex<double> values(dim, dim, dim, dim, 2);
	for (int lam = 0; lam <= maxL; lam++) {
		for (int mu = 0; mu <= lam; mu++) {
			ThreeIndex<double> Uij = uklm(lam, mu, fac);
			for (int i = 0; i <= lam; i++) {
				for (int j = 0; j <= lam - i; j++){
					values(lam, mu, i, j, 0) = Uij(i, j, 0);
					values(lam, mu, i, j, 1) = Uij(i, j, 1);
				}
			}
		}
	}
	
	return values;
}

void AngularIntegral::makeW(std::vector<double> &fac, FiveIndex<double> &U) {
	int LB2 = 2*LB;
	int dim = wDim;
	int maxI = (maxL + dim)/2;
	int maxLam = maxL;
	
	FiveIndex<double> values{dim+1, dim+1, dim+1, maxLam+1, 2*(maxLam + 1)};
	ThreeIndex<double> pijk = Pijk(maxI);
	
	int plam, pmu;
	double smu, w;
	std::vector<int> ix(3);
	for (int k = 0; k <= dim; k++) {	
		for (int l = 0; l <= dim; l++) {	
			for(int m = 0; m <= dim; m++) {
				plam = (k + l + m)%2;
				
				int limit = maxLam > k+l+m ? k+l+m : maxLam;
				for(int lam = plam; lam <= limit; lam += 2){
					smu = 1 - 2*(l%2);
					pmu = (k+l) % 2;
					
					for (int mu = pmu; mu <= lam; mu+=2) {
						w = 0.0;
						for (int i = 0; i <= lam; i++) {
							for (int j = 0; j <= lam - i; j++) {
								ix[0] = k+i;
								ix[1] = l+j;
								ix[2] = m + lam - i - j; 
								
								if (ix[0]%2 + ix[1]%2 + ix[2]%2 == 0){
									std::sort(ix.begin(), ix.end()); 
									w += U(lam, mu, i, j, (1 - (int)(smu))/2)*pijk(ix[2]/2, ix[1]/2, ix[0]/2);
								}
								
							}
						}
						
						values(k, l, m, lam, lam+(int)(smu*mu)) = w;
					}
				}	
			}	
		}	
	}
	W = values;
}

void AngularIntegral::makeOmega(FiveIndex<double> &U) {
	
	int lamDim = LE + LB; 
	int muDim = 2*lamDim + 1;
	SevenIndex<double> values{LB+1, LB+1, LB+1, lamDim+1, muDim+1, lamDim+1, muDim+1};
	
	double om_plus=0.0, om_minus=0.0;
	double wval; 
	bool test1, test2, test3;
	for (int k = 0; k <= LB; k++) {
		for (int l = 0; l <= LB; l++) {
			for (int m = 0; m <= LB; m++) {
					
				for (int rho = 0; rho <= lamDim; rho++ ) {
					for (int sigma = -rho; sigma <= rho; sigma++) {
						
						for (int lam = 0; lam <= rho; lam++) {
							test1 = (k+l+m+lam) % 2 == rho % 2;
	
							for (int mu = 0; mu <= lam; mu++) {
								
								om_plus = om_minus = 0.0;
								for (int i = 0; i<= lam; i++ ) {
									for (int j = 0; j <= lam - i; j++) {
											wval = W(k+i, l+j, m+lam-i-j, rho, rho+sigma);
											om_plus += U(lam, mu, i, j, 0) * wval;
											om_minus += U(lam, mu, i, j, 1) * wval;
									}
								}
								if (mu == 0) om_minus = om_plus;
								values(k, l, m, rho, sigma+rho, lam, lam+mu) = om_plus;
								values(k, l, m, lam, lam+mu, rho, sigma+rho) = om_plus;
								values(k, l, m, rho, sigma+rho, lam, lam-mu) = om_minus;
								values(k, l, m, lam, lam-mu, rho, sigma+rho) = om_minus;
								
							}
						}
						
					}
				}
					
			}
		}
	}
	
	omega = values;
}

AngularIntegral::AngularIntegral() { init(0, 0); }
AngularIntegral::AngularIntegral(int _LB, int _LE) { init(_LB, _LE); }
void AngularIntegral::init(int _LB, int _LE ) {
	LB = _LB;
	LE = _LE;
	wDim = 4*LB > 3*LB + LE ? 4*LB : 3*LB + LE;
	maxL = 2*LB > LB + LE ? 2*LB : LB+LE;
	
}

void AngularIntegral::compute() {
	std::vector<double> fac = facArray(2*wDim);
	
	FiveIndex<double> U = makeU(fac);
	makeW(fac, U);
	makeOmega(U);
}

void AngularIntegral::clear() {}

double AngularIntegral::getIntegral(int k, int l, int m, int lam, int mu) const { return W(k, l, m, lam, lam+mu); }
double AngularIntegral::getIntegral(int k, int l, int m, int lam, int mu, int rho, int sigma) const { return omega(k, l, m, lam, lam+mu, rho, rho+sigma); }

bool AngularIntegral::isZero(int k, int l, int m, int lam, int mu, double tolerance) const {
	if (wDim > 0) return fabs(W(k, l, m, lam, lam+mu)) < tolerance;
	else return true;
}
bool AngularIntegral::isZero(int k, int l, int m, int lam, int mu, int rho, int sigma, double tolerance) const {
	if (wDim > 0) return fabs(omega(k, l, m, lam, lam+mu, rho, rho+sigma)) < tolerance;
	else return true;
}

//****************************************** RADIAL INTEGRAL *********************************************

RadialIntegral::RadialIntegral() {}

void RadialIntegral::init(int maxL, double tol, int small, int large) {
	bigGrid.initGrid(large, ONEPOINT);
	smallGrid.initGrid(small, TWOPOINT);
	smallGrid.transformZeroInf();
	
	bessie.init(maxL, 1600, 200, tol);
	
	tolerance = tol;
}

void RadialIntegral::buildBessel(std::vector<double> &r, int nr, int maxL, TwoIndex<double> &values, double weight) {
	std::vector<double> besselValues;
	for (int i = 0; i < nr; i++) {
		bessie.calculate(weight * r[i], maxL, besselValues);
		for (int l = 0; l <= maxL; l++) values(l, i) = besselValues[l];
	}
}

double RadialIntegral::calcKij(double Na, double Nb, double zeta_a, double zeta_b, double *A, double *B) const {
	double muij = zeta_a * zeta_b / (zeta_a + zeta_b);
	double R[3] = {A[0] - B[0], A[1] - B[1], A[2] - B[2]};
	double R2 = R[0] * R[0] + R[1] * R[1] + R[2] * R[2];
	return Na * Nb * exp(-muij * R2);
}

// Assumes that p is the pretabulated integrand at the abscissae
double RadialIntegral::integrand(double r, double *p, int ix) {
	return p[ix];
}

void RadialIntegral::buildParameters(const GaussianShell &shellA, const GaussianShell &shellB, double *A, double *B) {
	int npA = shellA.nprimitive();
	int npB = shellB.nprimitive();

	p.assign(npA, npB, 0.0);
	P.assign(npA, npB, 0.0);
	P2.assign(npA, npB, 0.0);
	K.assign(npA, npB, 0.0);

	double Pvec[3];
	double zetaA, zetaB;
	for (int a = 0; a < npA; a++) {
		zetaA = shellA.exp(a);
		
		for (int b = 0; b < npB; b++) {
			zetaB = shellB.exp(b);
			
			p(a, b) = zetaA + zetaB;
			for (int n = 0; n < 3; n++) 
				Pvec[n] = (zetaA * A[n] + zetaB * B[n])/p(a, b);
			
			P2(a, b) = Pvec[0]*Pvec[0] + Pvec[1]*Pvec[1] + Pvec[2]*Pvec[2];
			P(a, b) = sqrt(P2(a, b));
			K(a, b) = calcKij(1.0, 1.0, zetaA, zetaB, A, B);
			
		}
	}
}

void RadialIntegral::buildU(ECP &U, int l, int N, GCQuadrature &grid, double *Utab) {
	int gridSize = grid.getN();
	std::vector<double> &gridPoints = grid.getX();
	
	// Tabulate weighted ECP values
	double r;
	bool foundStart = false;
	double startTolerance = tolerance / 1000.0; // Distributions are left-skewed
	for (int i = 0; i < gridSize; i++) {
		r = gridPoints[i];
		Utab[i] = pow(r, N+2) * U.evaluate(r, l);
		
		// Set pre-screening limits
		if(Utab[i] > startTolerance && !foundStart) {
			 grid.start = i;
			 foundStart = true;
		}
		if(Utab[i] < tolerance && foundStart) {
			grid.end = i-1;
			foundStart = false;
		}
	}
}

int RadialIntegral::integrate(int maxL, int gridSize, TwoIndex<double> &intValues, GCQuadrature &grid, std::vector<double> &values, int offset, int skip) {
	std::function<double(double, double*, int)> intgd = integrand; 
	values.assign(maxL+1, 0.0);
	int test;
	double params[gridSize];
	for (int i = 0; i < grid.start; i++) params[i] = 0.0;
	for (int i = grid.end+1; i < gridSize; i++) params[i] = 0.0;
	for (int l = offset; l <= maxL; l+=skip) {
		for (int i = grid.start; i <= grid.end; i++) params[i] = intValues(l, i); 
		test = grid.integrate(intgd, params, tolerance);
		values[l] = grid.getI();
		if (test == 0) break;
	}
	return test;
}

void RadialIntegral::type1(int maxL, int N, int offset, ECP &U, const GaussianShell &shellA, const GaussianShell &shellB, double *Avec, double *Bvec, TwoIndex<double> &values) {
	int npA = shellA.nprimitive();
	int npB = shellB.nprimitive();
	
	buildParameters(shellA, shellB, Avec, Bvec);
	
	int gridSize = bigGrid.getN();

	// Now pretabulate integrand
	TwoIndex<double> intValues(maxL+1, gridSize, 0.0);
	// and bessel function
	TwoIndex<double> besselValues(maxL+1, gridSize, 0.0);
	// Calculate type1 integrals
	double da, db, za, zb, val;
	double A = Avec[0]*Avec[0] + Avec[1]*Avec[1] + Avec[2]*Avec[2];
	double B = Bvec[0]*Bvec[0] + Bvec[1]*Bvec[1] + Bvec[2]*Bvec[2];
	A = sqrt(A); B = sqrt(B);
	std::vector<double> tempValues;
	values.assign(maxL+1, 2*maxL + 1, 0.0);
	
	std::vector<double> fac = facArray(2*maxL);
	std::vector<double> dfac = dfacArray(2*maxL);
	
	// Tabulate integrand
	double x, phi, Px, Py;
	for (int a = 0; a < npA; a++) {
		da = shellA.coef(a);
		za = shellA.exp(a);
		
		for (int b = 0; b < npB; b++) {
			db = shellB.coef(b);
			zb = shellB.exp(b);
			
			// Reset grid starting points
			GCQuadrature newGrid = bigGrid;
			newGrid.transformRMinMax(p(a, b), (za * A + zb * B)/p(a, b));
			std::vector<double> &gridPoints = newGrid.getX();
			newGrid.start = 0;
			newGrid.end = gridSize-1;
			
			// Build U and bessel tabs
			double Utab[gridSize];
			buildU(U, U.getL(), N, newGrid, Utab);
			buildBessel(gridPoints, gridSize, maxL, besselValues, 2.0*p(a,b)*P(a,b));
			
			for (int i = newGrid.start; i <= newGrid.end; i++) {
				val = -p(a, b) * (gridPoints[i]*(gridPoints[i] - 2*P(a, b)) + P2(a, b));
				val = exp(val);
				for (int l = offset; l <= maxL; l+=2)
					intValues(l, i) = Utab[i] * val * besselValues(l, i);
			}

			int test = integrate(maxL, gridSize, intValues, newGrid, tempValues, offset, 2);
			if (test == 0) std::cout << "Failed to converge: \n";
				
			// Calculate real spherical harmonic
			x = fabs(P(a, b)) < 1e-12 ? 0.0 : (za * Avec[2] + zb * Bvec[2]) / (p(a, b) * P(a, b));
			Py = (za * Avec[1] + zb * Bvec[1]) / p(a, b);
			Px = (za * Avec[0] + zb * Bvec[0]) / p(a, b);
			phi = atan2(Py, Px);

			TwoIndex<double> harmonics = realSphericalHarmonics(maxL, x, phi, fac, dfac);
			for (int l = offset; l <= maxL; l+=2) {
				for (int mu = -l; mu <= l; mu++)
					values(l, l+mu) += da * db * harmonics(l, l+mu) * K(a, b) * tempValues[l];
			}
		}
	}
}

// F_a(lam, r) = sum_{i in a} d_i K_{lam}(2 zeta_a A r)*exp(-zeta_a(r - A)^2)
void RadialIntegral::buildF(const GaussianShell &shell, double *Avec, int lstart, int lend, std::vector<double> &r, int nr, int start, int end, TwoIndex<double> &F) {
	int np = shell.nprimitive();
	double A = Avec[0]*Avec[0] + Avec[1]*Avec[1] + Avec[2]*Avec[2];
	A = sqrt(A); 
		
	double weight, zeta, c;
	TwoIndex<double> besselValues(lend+1, nr, 0.0);
	
	F.assign(lend + 1, nr, 0.0);
	for (int a = 0; a < np; a++) {
		zeta = shell.exp(a);
		c = shell.coef(a);
		weight = 2.0 * zeta * A;
		
		buildBessel(r, nr, lend, besselValues, weight);
		
		for (int i = start; i <= end; i++) {
			weight = r[i] - A;
			weight = c * exp(-zeta * weight * weight);
			
			for (int l = lstart; l <= lend; l+=2) 
				F(l, i) += weight * besselValues(l, i); 
		}
	}
}

void RadialIntegral::type2(int l, int l1start, int l1end, int l2start, int l2end, int N, ECP &U, const GaussianShell &shellA, const GaussianShell &shellB, double *Avec, double *Bvec, TwoIndex<double> &values) {
	int npA = shellA.nprimitive();
	int npB = shellB.nprimitive();
	
	//buildParameters(shellA, shellB, Avec, Bvec);
	
	// Start with the small grid
	// Pretabulate U
	int gridSize = smallGrid.getN();
	std::vector<double> &gridPoints = smallGrid.getX();
	
	// Reset grid starting points
	smallGrid.start = 0;
	smallGrid.end = gridSize-1;
	
	double Utab[gridSize];
	buildU(U, l, N, smallGrid, Utab);
	
	// Build the F matrices
	TwoIndex<double> Fa;
	TwoIndex<double> Fb;
	buildF(shellA, Avec, l1start, l1end, gridPoints, gridSize, smallGrid.start, smallGrid.end, Fa);
	buildF(shellB, Bvec, l2start, l2end, gridPoints, gridSize, smallGrid.start, smallGrid.end, Fb);
	
	// Build the integrals
	TwoIndex<double> intValues(l2end +1, gridSize, 0.0);
	std::vector<int> tests(l1end +1);
	std::vector<double> tempValues;
	bool failed = false;
	values.assign(l1end+1, l2end+1, 0.0);
	for (int l1 = l1start; l1 <= l1end; l1+=2) {
		for (int i = smallGrid.start; i <= smallGrid.end; i++) {
			for (int l2 = l2start; l2 <= l2end; l2+=2) 
				intValues(l2, i) = Utab[i] * Fa(l1, i) * Fb(l2, i);
		}
		tests[l1] = integrate(l2end, gridSize, intValues, smallGrid, tempValues, l2start, 2);
		failed = failed || (tests[l1] == 0);
		for (int l2 = l2start; l2 <= l2end; l2+=2) values(l1, l2) = tempValues[l2];
	}
	
	if (failed) {
		std::cerr << "Failed at first attempt\n";
		// Not converged, switch to big grid
		double zeta_a, zeta_b, c_a, c_b, weight, XA, XB;
		double A = Avec[0]*Avec[0] + Avec[1]*Avec[1] + Avec[2]*Avec[2];
		double B = Bvec[0]*Bvec[0] + Bvec[1]*Bvec[1] + Bvec[2]*Bvec[2];
		A = sqrt(A); B = sqrt(B);
		
		gridSize = bigGrid.getN();
		intValues.assign(l2end+1, gridSize, 0.0);
		Fa.assign(l2end+1, gridSize, 0.0);
		Fb.assign(l2end+1, gridSize, 0.0);
		
		for (int l1 = l1start; l1 <= l1end; l1+=2) {
			if (tests[l1] == 0) { 
				for (int l2 = l2start; l2 <= l2end; l2+=2) values(l1, l2) = 0.0;
			
				for (int a = 0; a < npA; a++) {
					zeta_a = shellA.exp(a);
					c_a = shellA.coef(a);
					
					for (int b = 0; b < npB; b++) {
						zeta_b = shellB.exp(b);
						c_b = shellB.coef(b); 
					
						// Set up grid
						GCQuadrature newGrid = bigGrid;
						std::vector<double> &gridPoints2 = newGrid.getX();
						newGrid.start = 0;
						newGrid.end = gridSize-1;
						newGrid.transformRMinMax(p(a,b), (zeta_a * A + zeta_b * B)/p(a, b));
				
						// Build the U tab
						double Utab2[gridSize];
						buildU(U, l, N, newGrid, Utab2);	
					
						// Build bessel function values
						weight = 2.0 * zeta_a * A;
						buildBessel(gridPoints2, gridSize, l2end, Fa, weight);
						weight = 2.0 * zeta_b * B;
						buildBessel(gridPoints2, gridSize, l2end, Fb, weight);
						for (int i = newGrid.start; i <= newGrid.end; i++) {
							XA = gridPoints2[i] - A;
							XA = exp(-zeta_a * XA * XA);
							XB = gridPoints2[i] - B;
							XB = exp(-zeta_b * XB * XB);
							for (int l2 = l2start; l2 <= l2end; l2+=2)
								intValues(l2, i) = Utab2[i] * XA * XB * Fa(l2, i) * Fb(l2, i);
						}
					
						if(integrate(l2end, gridSize, intValues, newGrid, tempValues, l2start, 2) == 0) {
							std::cerr << " Failed at second attempt!\n";
						}
						for (int l2 = l2start; l2 <= l2end; l2+=2) values(l1, l2) += c_a*c_b*tempValues[l2];
					}
				}
			}
		}
	}
	
}

//***************************************** ECP INTEGRAL ***********************************************

ECPIntegral::ECPIntegral(ECPBasis &_basis, int maxLB) : basis(_basis) {
	// Initialise angular and radial integrators
	angInts.init(maxLB, U.getL());
	angInts.compute();
	radInts.(2*maxLB + U.getL());
}

double ECPIntegral::calcC(int a, int m, double A, std::vector<double> &fac) const {
	double value = 1.0 - 2*((a-m) % 2);
	value *= pow(A, a-m);
	value *= fac[a]/(fac[m] * fac[a-m]);
	return value;
}

void ECPIntegral::makeC(ThreeIndex<double> &C, int L, double *A, std::vector<double> &fac) {
	int z; double Ck, Cl;
	for (int x = 0; x <= L; x++) {
		for (int y = 0; y <= L - x; y++) {
			z = L - x - y;
			
			for (int k = 0; k<= x; k++) {
				Ck = calcC(x, k, A[0], fac);
				for (int l = 0; l <= y; l++) {
					Cl = calcC(y, l, A[1], fac);
					for (int m = 0; m <= z; m++) C(k, l, m) = Ck * Cl * calcC(z, m, A[2], fac);
				}
			}
		}
	}
}

void ECPIntegral::type1(ECP &U, const GaussianShell &shellA, const GaussianShell &shellB, double *A, double *B, ThreeIndex<double> &CA, ThreeIndex<double> &CB, TwoIndex<double> &values) { 
	
	int LA = shellA.am(); int LB = shellB.am();
	int maxLBasis = LA > LB ? LA : LB;
	
	// Build radial integrals
	int L = LA + LB + U.getL();
	TwoIndex<double> temp;
	ThreeIndex<double> radials(L+1, L+1, 2*L+1);
	for (int ix = 0; ix <= L; ix++) {
		radInts.type1(ix, ix, ix % 2, U, shellA, shellB, A, B, temp);
		for(int l = 0; l <= ix; l++) {
			for (int m = -l; m <= l; m++) radials(ix, l, l+m) = temp(l, l+m);
		}
	}
	
	values.assign(shellA.ncartesian(), shellB.ncartesian(), 0.0);
	
	// Unpack positions
	double Ax = A[0]; double Ay = A[1]; double Az = A[2];
	double Bx = B[0]; double By = B[1]; double Bz = B[2];
	
	// Calculate chi_ab for all ab in shells
	int z1, z2, lparity, mparity, msign, ix, k, l, m;
	double C;
	int na = 0, nb = 0;
	for (int x1 = 0; x1 <= LA; x1++) {
		for (int y1 = 0; y1 <= LA - x1; y1++) {
			z1 = LA - x1 - y1;
			nb = 0;
			
			for (int x2 = 0; x2 <= LB; x2++) {
				for (int y2 = 0; y2 <= LB - x2; y2++) {
					z2 = LB - x2 - y2;
					
					for (int k1 = 0; k1 <= x1; k1++) {
						for (int k2 = 0; k2 <= x2; k2++) {
							k = k1 + k2;
							
							for (int l1 = 0; l1 <= y1; l1++) {
								for (int l2 = 0; l2 <= y2; l2++) {
									l = l1 + l2;
									
									for (int m1 = 0; m1 <= z1; m1++) {
										for (int m2 = 0; m2 <= z2; m2++){
											m = m1 + m2;
											C = CA(k1, l1, m1) * CB(k2, l2, m2);

											if ( fabs(C) > 1e-14 ) {
												// Build radial integrals
												ix = k + l + m;
												lparity = ix % 2;
												msign = 1 - 2*(l%2);
												mparity = (lparity + m) % 2;
												
												for (int lam = lparity; lam <= ix; lam+=2) {
													for (int mu = mparity; mu <= lam; mu+=2) 
														values(na, nb) += C * angInts.getIntegral(k, l, m, lam, msign*mu) * radials(ix, lam, lam+msign*mu);
												}
								
											}
										}
									}
								}
							}
						}
					}
					
					values(na, nb) *= 4.0 * M_PI;
					nb++;
				}
			}
			
			na++;
		}
	}
	
}

void ECPIntegral::type2(int lam, ECP& U, const GaussianShell &shellA, const GaussianShell &shellB, double *A, double *B, ThreeIndex<double> &CA, ThreeIndex<double> &CB, ThreeIndex<double> &values) {
	double prefac = 16.0 * M_PI * M_PI;
	int LA = shellA.am();
	int LB = shellB.am();
	int L = LA + LB;	
	int maxLBasis = LA > LB ? LA : LB;
	
	// Build radial integrals
	int lparity = lam % 2;
	int l1start,l2start;
	ThreeIndex<double> radials(L+1, lam + LA + 1, lam + LB + 1);
	TwoIndex<double> temp;
	for (int N1 = 0; N1 <= LA; N1++) {
		l1start = abs(lam - N1);
		
		for (int N2 = 0; N2 <= LB; N2++) {
			l2start = abs(lam - N2);
			
			radInts.type2(lam, l1start, lam + N1, l2start, lam + N2, N1 + N2, U, shellA, shellB, A, B, temp);
			for (int l1 = l1start; l1 <= lam+N1; l1+=2) {
				for (int l2 = l2start; l2 <= lam+N2; l2+=2) radials(N1 + N2, l1, l2) = temp(l1, l2);
			}
		}
	}
	
	// Loop over all basis functions in shell
	std::vector<double> fac = facArray(2*(lam + maxLBasis) + 1);
	std::vector<double> dfac = dfacArray(2*(lam + maxLBasis) + 1);
	
	// Unpack positions
	double Ax = A[0]; double Ay = A[1]; double Az = A[2];
	double Bx = B[0]; double By = B[1]; double Bz = B[2];
	double A2 = Ax * Ax + Ay * Ay + Az * Az;
	double B2 = Bx * Bx + By * By + Bz * Bz;
	double Am = sqrt(A2); double Bm = sqrt(B2);
	
	// Build spherical harmonics
	double xA = Am > 0 ? Az / Am : 0.0;
	double xB = Bm > 0 ? Bz / Bm : 0.0;
	double phiA = atan2(Ay, Ax);
	double phiB = atan2(By, Bx);
	TwoIndex<double> SA = realSphericalHarmonics(lam+LA, xA, phiA, fac, dfac);
	TwoIndex<double> SB = realSphericalHarmonics(lam+LB, xB, phiB, fac, dfac);
	
	// Calculate chi_ab for all ab in shells
	int z1, z2, ix, N1, N2;
	double val, C;
	int na = 0, nb = 0;
	for (int x1 = 0; x1 <= LA; x1++) {
		for (int y1 = 0; y1 <= LA - x1; y1++) {
			z1 = LA - x1 - y1;
			nb = 0;
			
			for (int x2 = 0; x2 <= LB; x2++) {
				for (int y2 = 0; y2 <= LB - x2; y2++) {
					z2 = LB - x2 - y2;
					
					for (int k1 = 0; k1 <= x1; k1++) {
						for (int k2 = 0; k2 <= x2; k2++) {
							for (int l1 = 0; l1 <= y1; l1++) {
								
								for (int l2 = 0; l2 <= y2; l2++) {
									for (int m1 = 0; m1 <= z1; m1++) {
										for (int m2 = 0; m2 <= z2; m2++){
											C = CA(k1, l1, m1) * CB(k2, l2, m2);
											
											N1 = k1 + l1 + m1;
											N2 = k2 + l2 + m2;
											l1start = abs(lam - N1);
											l2start = abs(lam - N2);
											ix = N1 + N2;
											// Sum over rho/kappa/sigma/tau
											for (int rho = l1start; rho <= lam + N1; rho+=2) {
												for (int sigma = -rho; sigma <= rho; sigma++) {
													for (int kappa = l2start; kappa <= lam + N2; kappa += 2) {
														for (int tau = -kappa; tau <= kappa; tau++) {
															val = C * SA(rho, rho+sigma) * SB(kappa, kappa+tau) * radials(ix, rho, kappa);
								
															for (int mu = -lam; mu <= lam; mu++) 
															 	values(na, nb, lam+mu) += val * angInts.getIntegral(k1, l1, m1, lam, mu, rho, sigma) * angInts.getIntegral(k2, l2, m2, lam, mu, kappa, tau);
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				
					for (int mu = -lam; mu <= lam; mu++) values(na, nb, lam+mu) *= prefac;
					nb++;
				}
			}
			
			na++;
		}
	}
}

void ECPIntegral::compute_shell_pair(ECP &U, const GaussianShell &shellA, const GaussianShell &shellB, TwoIndex<double> &values) {
	
	// Shift A and B to be relative to U
	const double* C = U.center();
	double A[3] = {shellA.center()[0] - C[0], shellA.center()[1] - C[1], shellA.center()[2] - C[2]};
	double B[3] = {shellB.center()[0] - C[0], shellB.center()[1] - C[1], shellB.center()[2] - C[2]};
	
	int LA = shellA.am(); int LB = shellB.am();
	int maxLBasis = LA > LB ? LA : LB;
	
	std::vector<double> fac = facArray(maxLBasis);
	
	// Construct coefficients 
	ThreeIndex<double> CA(LA+1, LA+1, LA+1);
	ThreeIndex<double> CB(LB+1, LB+1, LB+1);
	makeC(CA, LA, A, fac);
	makeC(CB, LB, B, fac);
	
	// Calculate type1 integrals
	type1(U, shellA, shellB, A, B, CA, CB, values);
	
	// Now all the type2 integrals
	ThreeIndex<double> t2vals(shellA.ncartesian(), shellB.ncartesian(), 2*U.getL() + 1);
	for (int l = 0; l < U.getL(); l++) {
		type2(l, U, shellA, shellB, A, B, CA, CB, t2vals);
		for (int m = -l; m <= l; m++) {
			for(int na = 0; na < shellA.ncartesian(); na++) {
				for (int nb = 0; nb < shellB.ncartesian(); nb++) values(na, nb) += t2vals(na, nb, l+m); 
			}
		}
		t2vals.fill(0.0);
	}
	
}

void ECPIntegral::compute_pair(const GaussianShell &shellA, const GaussianShell &shellB) {
	memset(buffer_, 0, shellA.ncartesian() * shellB.ncartesian() * sizeof(double));
	
	TwoIndex<double> tempValues;
	int ao12;
	for (int i = 0; i < basis.getN(); i++) {
		 compute_shell_pair(basis.getECP(i), shellA, shellB, tempValues);
		 ao12 = 0;
		 for (int a = 0; a < shellA.ncartesian(); a++) {
			 for (int b = 0; b < shellB.ncartesian(); b++) {
				 buffer_[ao12++] += tempValues(a, b);
			 }
		 }
	}
}


