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

/* Implements ecpint.hpp */

#include "psi4/libmints/ecpint.h"
#include "psi4/libmints/gshell.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/potential.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/libmints/molecule.h"

#include "psi4/libciomr/libciomr.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>

namespace psi {

// Compute all the real spherical harmonics Slm(theta, phi) for l,m up to lmax
// x = cos (theta)
static TwoIndex<double> realSphericalHarmonics(int lmax, double x, double phi) {
    TwoIndex<double> rshValues(lmax + 1, 2 * lmax + 1, 0.0);

    if (lmax > 0) {
        // First calculate the associated Legendre polynomials, Plm(cos theta), using the recursion relation
        // (l-m)Plm = x(2l - 1)P{l-1}m - (l+m-1)P{l-2}m
        // along with the zeroth order term
        // Pmm = (-1)^m (2m-1)!!(1-x^2)^{m/2}
        double x2 = x * x;
        std::vector<std::vector<double>> Plm(lmax + 1, std::vector<double>(lmax + 1));
        // First get all Pmm terms
        Plm[0][0] = 1.0;
        // Make sure that 1-x^2 doesn't go below 0, due to roundoff
        double sox2 = sqrt(std::max(0.0, 1.0 - x2));
        double ox2m = 1.0;
        for (int m = 1; m <= lmax; m++) {
            ox2m *= -sox2;
            Plm[m][m] = ox2m * df[2 * m];
        }

        // Then increment l for each m
        Plm[1][0] = x;
        Plm[0][1] = 0.0;
        for (int l = 2; l <= lmax; l++) {
            ox2m = x * (2 * l - 1);
            for (int m = 0; m < l; m++) {
                Plm[l][m] = ox2m * Plm[l - 1][m] - (l + m - 1) * Plm[l - 2][m];
                Plm[l][m] /= ((double)(l - m));
            }
            Plm[l - 1][l] = 0.0;
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
                ox2m = (2.0 * l + 1.0) * fac[l - m] / fac[l + m];
                ox2m = sign * osq4pi * sqrt(2.0 * ox2m) * Plm[l][m];
                rshValues(l, l + m) = ox2m * cos(m * phi);
                rshValues(l, l - m) = ox2m * sin(m * phi);
                sign *= -1;
            }
        }

    } else {
        rshValues(0, 0) = 1.0 / sqrt(4.0 * M_PI);
    }

    return rshValues;
}

double AngularIntegral::calcG(int l, int m) const {
    double value = 0.0;
    double value1 = pow(2.0, l) * fac[l];
    value1 = 1.0 / value1;
    double value2 = (2.0 * l + 1) * fac[l - m] / (2.0 * M_PI * fac[l + m]);
    value2 = sqrt(value2);
    value = value1 * value2;
    return value;
}

double AngularIntegral::calcH1(int i, int j, int l, int m) const {
    double value = 0.0;

    value = fac[l] / (fac[j] * fac[l - i] * fac[i - j]);
    value *= (1 - 2 * (i % 2)) * fac[2 * (l - i)] / (fac[l - m - 2 * i]);

    return value;
}

double AngularIntegral::calcH2(int i, int j, int k, int m) const {
    double value = 0.0;
    int ki2 = k - 2 * i;
    if (m >= ki2 && ki2 >= 0) {
        value = fac[j] * fac[m] / (fac[i] * fac[j - i] * fac[ki2] * fac[m - ki2]);
        int p = (m - k + 2 * i) / 2;
        value *= (1.0 - 2.0 * (p % 2));
    }
    return value;
}

ThreeIndex<double> AngularIntegral::uklm(int lam, int mu) const {
    ThreeIndex<double> values(lam + 1, lam + 1, 2);

    double or2 = 1.0 / sqrt(2.0);
    double u = 0.0;
    double um = 0.0;
    double g = calcG(lam, mu);

    double u1, h1, h2;
    int j;
    for (int k = 0; k <= lam; k++) {
        for (int l = 0; l <= lam - k; l++) {
            u = um = 0.0;
            j = k + l - mu;
            if (j % 2 == 0 && j > -1) {
                u1 = 0.0;
                j /= 2;
                for (int i = j; i <= (lam - mu) / 2; i++) u1 += calcH1(i, j, lam, mu);

                u = g * u1;
                u1 = 0;
                for (int i = 0; i <= j; i++) u1 += calcH2(i, j, k, mu);
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
    int dim = maxI + 1;
    ThreeIndex<double> values(dim, dim, dim);
    double pi4 = 4.0 * M_PI;

    values(0, 0, 0) = pi4;
    for (int i = 1; i <= maxI; i++) {
        values(i, 0, 0) = pi4 / ((double)(2 * i + 1));

        for (int j = 1; j <= i; j++) {
            values(i, j, 0) = values(i, j - 1, 0) * (2.0 * j - 1.0) / (2.0 * ((double)(i + j)) + 1.0);

            for (int k = 1; k <= j; k++)
                values(i, j, k) = values(i, j, k - 1) * (2.0 * k - 1.0) / (2.0 * ((double)(i + j + k)) + 1.0);
        }
    }
    return values;
}

FiveIndex<double> AngularIntegral::makeU() {
    int dim = maxL + 1;

    FiveIndex<double> values(dim, dim, dim, dim, 2);
    for (int lam = 0; lam <= maxL; lam++) {
        for (int mu = 0; mu <= lam; mu++) {
            ThreeIndex<double> Uij = uklm(lam, mu);
            for (int i = 0; i <= lam; i++) {
                for (int j = 0; j <= lam - i; j++) {
                    values(lam, mu, i, j, 0) = Uij(i, j, 0);
                    values(lam, mu, i, j, 1) = Uij(i, j, 1);
                }
            }
        }
    }

    return values;
}

void AngularIntegral::makeW(FiveIndex<double> &U) {
    int LB2 = 2 * LB;
    int dim = wDim;
    int maxI = (maxL + dim) / 2;
    int maxLam = maxL;

    FiveIndex<double> values{dim + 1, dim + 1, dim + 1, maxLam + 1, 2 * (maxLam + 1)};
    ThreeIndex<double> pijk = Pijk(maxI);

    int plam, pmu;
    double smu, w;
    std::vector<int> ix(3);
    for (int k = 0; k <= dim; k++) {
        for (int l = 0; l <= dim; l++) {
            for (int m = 0; m <= dim; m++) {
                plam = (k + l + m) % 2;

                int limit = maxLam > k + l + m ? k + l + m : maxLam;
                for (int lam = plam; lam <= limit; lam += 2) {
                    smu = 1 - 2 * (l % 2);
                    pmu = (k + l) % 2;

                    for (int mu = pmu; mu <= lam; mu += 2) {
                        w = 0.0;
                        for (int i = 0; i <= lam; i++) {
                            for (int j = 0; j <= lam - i; j++) {
                                ix[0] = k + i;
                                ix[1] = l + j;
                                ix[2] = m + lam - i - j;

                                if (ix[0] % 2 + ix[1] % 2 + ix[2] % 2 == 0) {
                                    std::sort(ix.begin(), ix.end());
                                    w += U(lam, mu, i, j, (1 - (int)(smu)) / 2) * pijk(ix[2] / 2, ix[1] / 2, ix[0] / 2);
                                }
                            }
                        }

                        values(k, l, m, lam, lam + (int)(smu * mu)) = w;
                    }
                }
            }
        }
    }
    W = values;
}

void AngularIntegral::makeOmega(FiveIndex<double> &U) {
    int lamDim = LE + LB;
    int muDim = 2 * lamDim + 1;
    SevenIndex<double> values{LB + 1, LB + 1, LB + 1, lamDim + 1, muDim + 1, lamDim + 1, muDim + 1};

    double om_plus = 0.0, om_minus = 0.0;
    double wval;
    bool test1, test2, test3;
    for (int k = 0; k <= LB; k++) {
        for (int l = 0; l <= LB; l++) {
            for (int m = 0; m <= LB; m++) {
                for (int rho = 0; rho <= lamDim; rho++) {
                    for (int sigma = -rho; sigma <= rho; sigma++) {
                        for (int lam = 0; lam <= rho; lam++) {
                            test1 = (k + l + m + lam) % 2 == rho % 2;

                            for (int mu = 0; mu <= lam; mu++) {
                                om_plus = om_minus = 0.0;
                                for (int i = 0; i <= lam; i++) {
                                    for (int j = 0; j <= lam - i; j++) {
                                        wval = W(k + i, l + j, m + lam - i - j, rho, rho + sigma);
                                        om_plus += U(lam, mu, i, j, 0) * wval;
                                        om_minus += U(lam, mu, i, j, 1) * wval;
                                    }
                                }
                                if (mu == 0) om_minus = om_plus;
                                values(k, l, m, rho, sigma + rho, lam, lam + mu) = om_plus;
                                values(k, l, m, lam, lam + mu, rho, sigma + rho) = om_plus;
                                values(k, l, m, rho, sigma + rho, lam, lam - mu) = om_minus;
                                values(k, l, m, lam, lam - mu, rho, sigma + rho) = om_minus;
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
void AngularIntegral::init(int _LB, int _LE) {
    LB = _LB;
    LE = _LE;
    wDim = 4 * LB > 3 * LB + LE ? 4 * LB : 3 * LB + LE;
    maxL = 2 * LB > LB + LE ? 2 * LB : LB + LE;
}

void AngularIntegral::compute() {
    FiveIndex<double> U = makeU();
    makeW(U);
    makeOmega(U);
}

void AngularIntegral::clear() {}

double AngularIntegral::getIntegral(int k, int l, int m, int lam, int mu) const { return W(k, l, m, lam, lam + mu); }
double AngularIntegral::getIntegral(int k, int l, int m, int lam, int mu, int rho, int sigma) const {
    return omega(k, l, m, lam, lam + mu, rho, rho + sigma);
}

bool AngularIntegral::isZero(int k, int l, int m, int lam, int mu, double tolerance) const {
    if (wDim > 0)
        return std::fabs(W(k, l, m, lam, lam + mu)) < tolerance;
    else
        return true;
}
bool AngularIntegral::isZero(int k, int l, int m, int lam, int mu, int rho, int sigma, double tolerance) const {
    if (wDim > 0)
        return std::fabs(omega(k, l, m, lam, lam + mu, rho, rho + sigma)) < tolerance;
    else
        return true;
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

double RadialIntegral::calcKij(double Na, double Nb, double zeta_a, double zeta_b, double R2) const {
    double muij = zeta_a * zeta_b / (zeta_a + zeta_b);
    return Na * Nb * exp(-muij * R2);
}

// Assumes that p is the pretabulated integrand at the abscissae
double RadialIntegral::integrand(double r, double *p, int ix) { return p[ix]; }

void RadialIntegral::buildParameters(const GaussianShell &shellA, const GaussianShell &shellB, ShellPairData &data) {
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
            for (int n = 0; n < 3; n++) Pvec[n] = (zetaA * data.A[n] + zetaB * data.B[n]) / p(a, b);

            P2(a, b) = Pvec[0] * Pvec[0] + Pvec[1] * Pvec[1] + Pvec[2] * Pvec[2];
            P(a, b) = sqrt(P2(a, b));
            K(a, b) = calcKij(1.0, 1.0, zetaA, zetaB, data.RAB2);
        }
    }
}

void RadialIntegral::buildU(const GaussianShell &U, int l, int N, GCQuadrature &grid, double *Utab) {
    int gridSize = grid.getN();
    std::vector<double> &gridPoints = grid.getX();
    // Tabulate weighted ECP values
    double r;
    for (int i = 0; i < gridSize; i++) {
        r = gridPoints[i];
        Utab[i] = pow(r, N) * U.evaluate(r, l);
    }
}

int RadialIntegral::integrate(int maxL, int gridSize, TwoIndex<double> &intValues, GCQuadrature &grid,
                              std::vector<double> &values, int offset, int skip) {
    std::function<double(double, double *, int)> intgd = integrand;
    values.assign(maxL + 1, 0.0);
    int test;
    std::vector<double> params(gridSize);
    for (int i = 0; i < grid.start; i++) params[i] = 0.0;
    for (int i = grid.end + 1; i < gridSize; i++) params[i] = 0.0;
    for (int l = offset; l <= maxL; l += skip) {
        for (int i = grid.start; i <= grid.end; i++) params[i] = intValues(l, i);
        test = grid.integrate(intgd, params.data(), tolerance);
        values[l] = grid.getI();
        if (test == 0) break;
    }
    return test;
}

void RadialIntegral::type1(int maxL, int N, int offset, const GaussianShell &U, const GaussianShell &shellA,
                           const GaussianShell &shellB, ShellPairData &data, TwoIndex<double> &values) {
    int npA = shellA.nprimitive();
    int npB = shellB.nprimitive();

    buildParameters(shellA, shellB, data);

    int gridSize = bigGrid.getN();

    // Now pretabulate integrand
    TwoIndex<double> intValues(maxL + 1, gridSize, 0.0);
    // and bessel function
    TwoIndex<double> besselValues(maxL + 1, gridSize, 0.0);
    // Calculate type1 integrals
    double da, db, za, zb, val;
    double A = data.Am;
    double B = data.Bm;
    std::vector<double> tempValues;
    values.assign(maxL + 1, 2 * maxL + 1, 0.0);

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
            newGrid.transformRMinMax(p(a, b), (za * A + zb * B) / p(a, b));
            std::vector<double> &gridPoints = newGrid.getX();
            newGrid.start = 0;
            newGrid.end = gridSize - 1;

            // Build U and bessel tabs
            std::vector<double> Utab(gridSize);
            buildU(U, U.am(), N, newGrid, Utab.data());
            buildBessel(gridPoints, gridSize, maxL, besselValues, 2.0 * p(a, b) * P(a, b));

            // Start building intvalues, and prescreen
            bool foundStart = false, tooSmall = false;
            for (int i = 0; i < gridSize; i++) {
                for (int l = offset; l <= maxL; l += 2) {
                    intValues(l, i) = Utab[i] * besselValues(l, i);
                    tooSmall = intValues(l, i) < tolerance;
                }
                if (!tooSmall && !foundStart) {
                    foundStart = true;
                    newGrid.start = i;
                }
                if (tooSmall && foundStart) {
                    newGrid.end = i - 1;
                    break;
                }
            }

            for (int i = newGrid.start; i <= newGrid.end; i++) {
                val = -p(a, b) * (gridPoints[i] * (gridPoints[i] - 2 * P(a, b)) + P2(a, b));
                val = exp(val);
                for (int l = offset; l <= maxL; l += 2) intValues(l, i) *= val;
            }

            int test = integrate(maxL, gridSize, intValues, newGrid, tempValues, offset, 2);
            if (test == 0) std::cout << "Failed to converge: \n";

            // Calculate real spherical harmonic
            x = std::fabs(P(a, b)) < 1e-12 ? 0.0 : (za * data.A[2] + zb * data.B[2]) / (p(a, b) * P(a, b));
            Py = (za * data.A[1] + zb * data.B[1]) / p(a, b);
            Px = (za * data.A[0] + zb * data.B[0]) / p(a, b);
            phi = atan2(Py, Px);

            TwoIndex<double> harmonics = realSphericalHarmonics(maxL, x, phi);
            for (int l = offset; l <= maxL; l += 2) {
                for (int mu = -l; mu <= l; mu++)
                    values(l, l + mu) += da * db * harmonics(l, l + mu) * K(a, b) * tempValues[l];
            }
        }
    }
}

// F_a(lam, r) = sum_{i in a} d_i K_{lam}(2 zeta_a A r)*exp(-zeta_a(r - A)^2)
void RadialIntegral::buildF(const GaussianShell &shell, double A, int lstart, int lend, std::vector<double> &r, int nr,
                            int start, int end, TwoIndex<double> &F) {
    int np = shell.nprimitive();

    double weight, zeta, c;
    TwoIndex<double> besselValues(lend + 1, nr, 0.0);

    F.assign(lend + 1, nr, 0.0);
    for (int a = 0; a < np; a++) {
        zeta = shell.exp(a);
        c = shell.coef(a);
        weight = 2.0 * zeta * A;

        buildBessel(r, nr, lend, besselValues, weight);

        for (int i = start; i <= end; i++) {
            weight = r[i] - A;
            weight = c * exp(-zeta * weight * weight);

            for (int l = lstart; l <= lend; ++l) F(l, i) += weight * besselValues(l, i);
        }
    }
}

void RadialIntegral::type2(int l, int l1start, int l1end, int l2start, int l2end, int N, const GaussianShell &U,
                           const GaussianShell &shellA, const GaussianShell &shellB, ShellPairData &data,
                           TwoIndex<double> &values) {
    std::function<double(double, double *, int)> intgd = integrand;
    int npA = shellA.nprimitive();
    int npB = shellB.nprimitive();

    double A = data.Am;
    double B = data.Bm;

    // Start with the small grid
    // Pretabulate U
    int gridSize = smallGrid.getN();
    std::vector<double> &gridPoints = smallGrid.getX();

    // Reset grid starting points
    smallGrid.start = 0;
    smallGrid.end = gridSize - 1;

    std::vector<double> Utab(gridSize);
    buildU(U, l, N, smallGrid, Utab.data());
    values.assign(l1end + 1, l2end + 1, 0.0);

    // Build the F matrices
    // If shell is on same center as ECP, only l = 0 will be nonzero
    if (A < 1e-15) l1end = 0;
    if (B < 1e-15) l2end = 0;
    TwoIndex<double> Fa;
    TwoIndex<double> Fb;
    buildF(shellA, data.Am, l1start, l1end, gridPoints, gridSize, smallGrid.start, smallGrid.end, Fa);
    buildF(shellB, data.Bm, l2start, l2end, gridPoints, gridSize, smallGrid.start, smallGrid.end, Fb);

    // Build the integrals
    bool foundStart, tooSmall;
    std::vector<int> tests((l1end + 1) * (l2end + 1));
    std::vector<double> params(gridSize);
    bool failed = false;
    int ix = 0;
    for (int l1 = 0; l1 <= l1end; l1++) {
        for (int l2 = 0; l2 <= l2end; l2++) {
            for (int i = 0; i < gridSize; i++) params[i] = Utab[i] * Fa(l1, i) * Fb(l2, i);
            tests[ix] = smallGrid.integrate(intgd, params.data(), tolerance);
            failed = failed || (tests[ix] == 0);
            values(l1, l2) = tests[ix] == 0 ? 0.0 : smallGrid.getI();
            ix++;
        }
    }

    if (failed) {
        // Not converged, switch to big grid
        double zeta_a, zeta_b, c_a, c_b;

        gridSize = bigGrid.getN();
        Fa.assign(l1end + 1, gridSize, 0.0);
        Fb.assign(l2end + 1, gridSize, 0.0);

        for (int a = 0; a < npA; a++) {
            c_a = shellA.coef(a);
            zeta_a = shellA.exp(a);

            for (int b = 0; b < npB; b++) {
                c_b = shellB.coef(b);
                zeta_b = shellB.exp(b);

                GCQuadrature newGrid = bigGrid;
                newGrid.transformRMinMax(p(a, b), (zeta_a * A + zeta_b * B) / p(a, b));
                std::vector<double> &gridPoints2 = newGrid.getX();
                newGrid.start = 0;
                newGrid.end = gridSize - 1;

                // Build U and bessel tabs
                std::vector<double> Utab2(gridSize);
                buildU(U, l, N, newGrid, Utab2.data());
                buildBessel(gridPoints2, gridSize, l1end, Fa, 2.0 * zeta_a * A);
                buildBessel(gridPoints2, gridSize, l2end, Fb, 2.0 * zeta_b * B);

                std::vector<double> Xvals(gridSize);
                double ria, rib;
                for (int i = 0; i < gridSize; i++) {
                    ria = gridPoints2[i] - A;
                    rib = gridPoints2[i] - B;
                    Xvals[i] = exp(-zeta_a * ria * ria - zeta_b * rib * rib) * Utab2[i];
                }

                std::vector<double> params2(gridSize);
                int test;
                ix = 0;
                for (int l1 = 0; l1 <= l1end; l1++) {
                    for (int l2 = 0; l2 <= l2end; l2++) {
                        if (tests[ix] == 0) {
                            for (int i = 0; i < gridSize; i++) params2[i] = Xvals[i] * Fa(l1, i) * Fb(l2, i);
                            test = newGrid.integrate(intgd, params2.data(), tolerance);
                            if (test == 0) std::cerr << "Failed at second attempt" << std::endl;
                            values(l1, l2) += c_a * c_b * newGrid.getI();
                        }
                        ix++;
                    }
                }
            }
        }
    }
}

//***************************************** ECP INTEGRAL ***********************************************

ECPInt::ECPInt(std::vector<SphericalTransform> &st, std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
               int deriv)
    : OneBodyAOInt(st, bs1, bs2, deriv) {
    // Initialise angular and radial integrators
    int maxam1 = bs1->max_am();
    int maxam2 = bs2->max_am();
    int maxLB = maxam1 > maxam2 ? maxam1 : maxam2;
    int maxLU = bs1_->max_ecp_am();
    angInts.init(maxLB + deriv, maxLU);
    angInts.compute();
    radInts.init(2 * (maxLB + deriv) + maxLU);

    int maxnao1 = INT_NCART(maxam1);
    int maxnao2 = INT_NCART(maxam2);
    buffer_ = new double[maxnao1 * maxnao2];
}

ECPInt::~ECPInt() { delete[] buffer_; }

double ECPInt::calcC(int a, int m, double A) const {
    double value = 1.0 - 2 * ((a - m) % 2);
    value *= pow(A, a - m);
    value *= fac[a] / (fac[m] * fac[a - m]);
    return value;
}

void ECPInt::makeC(FiveIndex<double> &C, int L, double *A) {
    int z;
    double Ck, Cl;
    int na = 0;
    for (int x = L; x >= 0; x--) {
        for (int y = L - x; y >= 0; y--) {
            z = L - x - y;

            for (int k = 0; k <= x; k++) {
                Ck = calcC(x, k, A[0]);
                for (int l = 0; l <= y; l++) {
                    Cl = calcC(y, l, A[1]);
                    for (int m = 0; m <= z; m++) C(0, na, k, l, m) = Ck * Cl * calcC(z, m, A[2]);
                }
            }

            na++;
        }
    }
}

void ECPInt::type1(const GaussianShell &U, const GaussianShell &shellA, const GaussianShell &shellB,
                   ShellPairData &data, FiveIndex<double> &CA, FiveIndex<double> &CB, TwoIndex<double> &values) {
    int LA = data.LA;
    int LB = data.LB;
    int maxLBasis = data.maxLBasis;

    // Build radial integrals
    int L = LA + LB;
    TwoIndex<double> temp;
    ThreeIndex<double> radials(L + 1, L + 1, 2 * L + 1);
    for (int ix = 0; ix <= L; ix++) {
        radInts.type1(ix, ix, ix % 2, U, shellA, shellB, data, temp);
        for (int l = 0; l <= ix; l++) {
            for (int m = -l; m <= l; m++) {
                radials(ix, l, l + m) = temp(l, l + m);
            }
        }
    }

    // Unpack positions
    double Ax = data.A[0];
    double Ay = data.A[1];
    double Az = data.A[2];
    double Bx = data.B[0];
    double By = data.B[1];
    double Bz = data.B[2];

    // Calculate chi_ab for all ab in shells
    int z1, z2, lparity, mparity, msign, ix, k, l, m;
    double C;
    int na = 0, nb = 0;
    for (int x1 = LA; x1 >= 0; x1--) {
        for (int y1 = LA - x1; y1 >= 0; y1--) {
            z1 = LA - x1 - y1;
            nb = 0;

            for (int x2 = LB; x2 >= 0; x2--) {
                for (int y2 = LB - x2; y2 >= 0; y2--) {
                    z2 = LB - x2 - y2;

                    for (int k1 = 0; k1 <= x1; k1++) {
                        for (int k2 = 0; k2 <= x2; k2++) {
                            k = k1 + k2;

                            for (int l1 = 0; l1 <= y1; l1++) {
                                for (int l2 = 0; l2 <= y2; l2++) {
                                    l = l1 + l2;

                                    for (int m1 = 0; m1 <= z1; m1++) {
                                        for (int m2 = 0; m2 <= z2; m2++) {
                                            m = m1 + m2;
                                            C = CA(0, na, k1, l1, m1) * CB(0, nb, k2, l2, m2);

                                            if (std::fabs(C) > 1e-14) {
                                                ix = k + l + m;
                                                lparity = ix % 2;
                                                msign = 1 - 2 * (l % 2);
                                                mparity = (lparity + m) % 2;

                                                for (int lam = lparity; lam <= ix; lam += 2) {
                                                    for (int mu = mparity; mu <= lam; mu += 2)
                                                        values(na, nb) +=
                                                            C * angInts.getIntegral(k, l, m, lam, msign * mu) *
                                                            radials(ix, lam, lam + msign * mu);
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

void ECPInt::type2(int lam, const GaussianShell &U, const GaussianShell &shellA, const GaussianShell &shellB,
                   ShellPairData &data, FiveIndex<double> &CA, FiveIndex<double> &CB, ThreeIndex<double> &values) {
    double prefac = 16.0 * M_PI * M_PI;
    int LA = data.LA;
    int LB = data.LB;
    int L = LA + LB;
    int maxLBasis = data.maxLBasis;

    ThreeIndex<double> radials(L + 1, lam + LA + 1, lam + LB + 1);
    TwoIndex<double> temp;
    for (int N = 0; N < L + 1; N++) {
        radInts.type2(lam, 0, lam + LA, 0, lam + LB, N, U, shellA, shellB, data, temp);
        for (int l1 = 0; l1 < lam + LA + 1; l1++)
            for (int l2 = 0; l2 < lam + LB + 1; l2++) radials(N, l1, l2) = temp(l1, l2);
    }

    double Ax = data.A[0];
    double Ay = data.A[1];
    double Az = data.A[2];
    double Bx = data.B[0];
    double By = data.B[1];
    double Bz = data.B[2];
    double Am = data.Am;
    double Bm = data.Bm;
    double xA = Am > 0 ? Az / Am : 0.0;
    double xB = Bm > 0 ? Bz / Bm : 0.0;
    double phiA = atan2(Ay, Ax);
    double phiB = atan2(By, Bx);
    TwoIndex<double> SA = realSphericalHarmonics(lam + LA, xA, phiA);
    TwoIndex<double> SB = realSphericalHarmonics(lam + LB, xB, phiB);

    int z1, z2;
    double C, val1, val2;
    int na = 0;
    for (int x1 = LA; x1 >= 0; x1--) {
        for (int r1 = LA - x1; r1 >= 0; r1--) {
            z1 = LA - x1 - r1;

            int nb = 0;
            for (int x2 = LB; x2 >= 0; x2--) {
                for (int y2 = LB - x2; y2 >= 0; y2--) {
                    z2 = LB - x2 - y2;

                    for (int alpha_x = 0; alpha_x <= x1; alpha_x++) {
                        for (int alpha_y = 0; alpha_y <= r1; alpha_y++) {
                            for (int alpha_z = 0; alpha_z <= z1; alpha_z++) {
                                int alpha = alpha_x + alpha_y + alpha_z;

                                for (int beta_x = 0; beta_x <= x2; beta_x++) {
                                    for (int beta_y = 0; beta_y <= y2; beta_y++) {
                                        for (int beta_z = 0; beta_z <= z2; beta_z++) {
                                            int beta = beta_x + beta_y + beta_z;
                                            int N = alpha + beta;
                                            C = CA(0, na, alpha_x, alpha_y, alpha_z) *
                                                CB(0, nb, beta_x, beta_y, beta_z);

                                            for (int lam1 = 0; lam1 <= lam + alpha; lam1++) {
                                                for (int lam2 = 0; lam2 <= lam + beta; lam2++) {
                                                    val1 = prefac * C * radials(N, lam1, lam2);

                                                    for (int mu1 = -lam1; mu1 <= lam1; mu1++) {
                                                        for (int mu2 = -lam2; mu2 <= lam2; mu2++) {
                                                            val2 = val1 * SA(lam1, lam1 + mu1) * SB(lam2, lam2 + mu2);

                                                            for (int mu = -lam; mu <= lam; mu++)
                                                                values(na, nb, lam + mu) +=
                                                                    val2 *
                                                                    angInts.getIntegral(alpha_x, alpha_y, alpha_z, lam,
                                                                                        mu, lam1, mu1) *
                                                                    angInts.getIntegral(beta_x, beta_y, beta_z, lam, mu,
                                                                                        lam2, mu2);
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

                    nb++;
                }
            }

            na++;
        }
    }
}

void ECPInt::compute_shell_pair(const GaussianShell &U, const GaussianShell &shellA, const GaussianShell &shellB,
                                TwoIndex<double> &values, int shiftA, int shiftB) {
    ShellPairData data;
    // Shift A and B to be relative to U
    const double *C = U.center();
    data.A[0] = shellA.center()[0] - C[0];
    data.A[1] = shellA.center()[1] - C[1];
    data.A[2] = shellA.center()[2] - C[2];
    data.B[0] = shellB.center()[0] - C[0];
    data.B[1] = shellB.center()[1] - C[1];
    data.B[2] = shellB.center()[2] - C[2];

    data.LA = shellA.am() + shiftA;
    data.LB = shellB.am() + shiftB;
    data.maxLBasis = data.LA > data.LB ? data.LA : data.LB;

    data.ncartA = (data.LA + 1) * (data.LA + 2) / 2;
    data.ncartB = (data.LB + 1) * (data.LB + 2) / 2;

    data.A2 = data.A[0] * data.A[0] + data.A[1] * data.A[1] + data.A[2] * data.A[2];
    data.Am = sqrt(data.A2);
    data.B2 = data.B[0] * data.B[0] + data.B[1] * data.B[1] + data.B[2] * data.B[2];
    data.Bm = sqrt(data.B2);
    double RAB[3] = {data.A[0] - data.B[0], data.A[1] - data.B[1], data.A[2] - data.B[2]};
    data.RAB2 = RAB[0] * RAB[0] + RAB[1] * RAB[1] + RAB[2] * RAB[2];
    data.RABm = sqrt(data.RAB2);

    // Construct coefficients
    FiveIndex<double> CA(1, data.ncartA, data.LA + 1, data.LA + 1, data.LA + 1);
    FiveIndex<double> CB(1, data.ncartB, data.LB + 1, data.LB + 1, data.LB + 1);
    makeC(CA, data.LA, data.A);
    makeC(CB, data.LB, data.B);

    // Initialize the values to zero
    values.assign(data.ncartA, data.ncartB, 0.0);

    if (U.shell_type() == ECPType1) {
        // This is a type 1 shell
        type1(U, shellA, shellB, data, CA, CB, values);
    } else if (U.shell_type() == ECPType2) {
        // This is a type 2 shell
        int l = U.am();
        ThreeIndex<double> t2vals(data.ncartA, data.ncartB, 2 * l + 1);
        t2vals.fill(0.0);
        type2(l, U, shellA, shellB, data, CA, CB, t2vals);
        for (int m = -l; m <= l; m++) {
            for (int na = 0; na < data.ncartA; na++) {
                for (int nb = 0; nb < data.ncartB; nb++) values(na, nb) += t2vals(na, nb, l + m);
            }
        }
    } else {
        throw PSIEXCEPTION("Unrecognized shell type in ECPInt::compute_shell_pair.");
    }
}

void ECPInt::compute_pair(const GaussianShell &shellA, const GaussianShell &shellB) {
    memset(buffer_, 0, shellA.ncartesian() * shellB.ncartesian() * sizeof(double));
    TwoIndex<double> tempValues;
    int ao12;
    // TODO check that bs1 and bs2 ECPs are the same
    for (int i = 0; i < bs1_->n_ecp_shell(); i++) {
        const GaussianShell &ecpshell = bs1_->ecp_shell(i);
        compute_shell_pair(ecpshell, shellA, shellB, tempValues);
        ao12 = 0;
        for (int a = 0; a < shellA.ncartesian(); a++) {
            for (int b = 0; b < shellB.ncartesian(); b++) {
                buffer_[ao12++] += tempValues(a, b);
            }
        }
    }
}

ECPSOInt::ECPSOInt(const std::shared_ptr<OneBodyAOInt> &aoint, const std::shared_ptr<IntegralFactory> &fact)
    : OneBodySOInt(aoint, fact) {
    natom_ = ob_->basis1()->molecule()->natom();
}

ECPSOInt::ECPSOInt(const std::shared_ptr<OneBodyAOInt> &aoint, const IntegralFactory *fact)
    : OneBodySOInt(aoint, fact) {
    natom_ = ob_->basis1()->molecule()->natom();
}

}  // namespace psi
