#include <libciomr/libciomr.h>

#include "mints.h"
#include <physconst.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace boost;
using namespace psi;

// Initialize overlap_recur_ to +1 basis set angular momentum
KineticInt::KineticInt(std::vector<SphericalTransform>& st, boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2, int deriv) :
    OneBodyAOInt(st, bs1, bs2, deriv), overlap_recur_(bs1->max_am()+1+deriv, bs2->max_am()+1+deriv)
{
    if (deriv > 2)
        throw std::runtime_error("KineticInt: does not support deriv over 2.");

    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = INT_NCART(maxam1);
    int maxnao2 = INT_NCART(maxam2);

    if (deriv == 1) {
        // We set chunk count for normalize_am and pure_transform
        set_chunks(6);

        maxnao1 *= 3;
        maxnao2 *= 3;
    }
    else if (deriv == 2) {
        set_chunks(6);
        maxnao1 *= 6;
        maxnao2 *= 1;
    }

    buffer_ = new double[maxnao1*maxnao2];
}

KineticInt::~KineticInt()
{
    delete[] buffer_;
}

// The engine only supports segmented basis sets
void KineticInt::compute_pair(const GaussianShell& s1, const GaussianShell& s2)
{
    int ao12;
    int am1 = s1.am();
    int am2 = s2.am();
    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();
    double A[3], B[3];
    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, s1.ncartesian() * s2.ncartesian() * sizeof(double));

    double **x = overlap_recur_.x();
    double **y = overlap_recur_.y();
    double **z = overlap_recur_.z();

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1.exp(p1);
        double c1 = s1.coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2.exp(p2);
            double c2 = s2.coef(p2);
            double gamma = a1 + a2;
            double oog = 1.0/gamma;

            double PA[3], PB[3];
            double P[3];

            P[0] = (a1*A[0] + a2*B[0])*oog;
            P[1] = (a1*A[1] + a2*B[1])*oog;
            P[2] = (a1*A[2] + a2*B[2])*oog;
            PA[0] = P[0] - A[0];
            PA[1] = P[1] - A[1];
            PA[2] = P[2] - A[2];
            PB[0] = P[0] - B[0];
            PB[1] = P[1] - B[1];
            PB[2] = P[2] - B[2];

            double over_pf = exp(-a1*a2*AB2*oog) * sqrt(M_PI*oog) * M_PI * oog * c1 * c2;

            // Do recursion
            overlap_recur_.compute(PA, PB, gamma, am1+1, am2+1);

            ao12 = 0;
            for(int ii = 0; ii <= am1; ii++) {
                int l1 = am1 - ii;
                for(int jj = 0; jj <= ii; jj++) {
                    int m1 = ii - jj;
                    int n1 = jj;
                    /*--- create all am components of sj ---*/
                    for(int kk = 0; kk <= am2; kk++) {
                        int l2 = am2 - kk;
                        for(int ll = 0; ll <= kk; ll++) {
                            int m2 = kk - ll;
                            int n2 = ll;

                            double I1, I2, I3, I4;

                            I1 = (l1 == 0 || l2 == 0) ? 0.0 : x[l1-1][l2-1] * y[m1][m2] * z[n1][n2] * over_pf;
                            I2 = x[l1+1][l2+1] * y[m1][m2] * z[n1][n2] * over_pf;
                            I3 = (l2 == 0) ? 0.0 : x[l1+1][l2-1] * y[m1][m2] * z[n1][n2] * over_pf;
                            I4 = (l1 == 0) ? 0.0 : x[l1-1][l2+1] * y[m1][m2] * z[n1][n2] * over_pf;
                            double Ix = 0.5 * l1 * l2 * I1 + 2.0 * a1 * a2 * I2 - a1 * l2 * I3 - l1 * a2 * I4;

                            I1 = (m1 == 0 || m2 == 0) ? 0.0 : x[l1][l2] * y[m1-1][m2-1] * z[n1][n2] * over_pf;
                            I2 = x[l1][l2] * y[m1+1][m2+1] * z[n1][n2] * over_pf;
                            I3 = (m2 == 0) ? 0.0 : x[l1][l2] * y[m1+1][m2-1] * z[n1][n2] * over_pf;
                            I4 = (m1 == 0) ? 0.0 : x[l1][l2] * y[m1-1][m2+1] * z[n1][n2] * over_pf;
                            double Iy = 0.5 * m1 * m2 * I1 + 2.0 * a1 * a2 * I2 - a1 * m2 * I3 - m1 * a2 * I4;

                            I1 = (n1 == 0 || n2 == 0) ? 0.0 : x[l1][l2] * y[m1][m2] * z[n1-1][n2-1] * over_pf;
                            I2 = x[l1][l2] * y[m1][m2] * z[n1+1][n2+1] * over_pf;
                            I3 = (n2 == 0) ? 0.0 : x[l1][l2] * y[m1][m2] * z[n1+1][n2-1] * over_pf;
                            I4 = (n1 == 0) ? 0.0 : x[l1][l2] * y[m1][m2] * z[n1-1][n2+1] * over_pf;
                            double Iz = 0.5 * n1 * n2 * I1 + 2.0 * a1 * a2 * I2 - a1 * n2 * I3 - n1 * a2 * I4;

                            buffer_[ao12++] += (Ix + Iy + Iz);
                        }
                    }
                }
            }
        }
    }
}

static double ke_int(double **x, double **y, double **z, double a1, int l1, int m1, int n1,
          double a2, int l2, int m2, int n2)
{
    double I1, I2, I3, I4;

    I1 = (l1 == 0 || l2 == 0) ? 0.0 : x[l1-1][l2-1] * y[m1][m2] * z[n1][n2];
    I2 = x[l1+1][l2+1] * y[m1][m2] * z[n1][n2];
    I3 = (l2 == 0) ? 0.0 : x[l1+1][l2-1] * y[m1][m2] * z[n1][n2];
    I4 = (l1 == 0) ? 0.0 : x[l1-1][l2+1] * y[m1][m2] * z[n1][n2];
    double Ix = 0.5 * l1 * l2 * I1 + 2.0 * a1 * a2 * I2 - a1 * l2 * I3 - l1 * a2 * I4;

    I1 = (m1 == 0 || m2 == 0) ? 0.0 : x[l1][l2] * y[m1-1][m2-1] * z[n1][n2];
    I2 = x[l1][l2] * y[m1+1][m2+1] * z[n1][n2];
    I3 = (m2 == 0) ? 0.0 : x[l1][l2] * y[m1+1][m2-1] * z[n1][n2];
    I4 = (m1 == 0) ? 0.0 : x[l1][l2] * y[m1-1][m2+1] * z[n1][n2];
    double Iy = 0.5 * m1 * m2 * I1 + 2.0 * a1 * a2 * I2 - a1 * m2 * I3 - m1 * a2 * I4;

    I1 = (n1 == 0 || n2 == 0) ? 0.0 : x[l1][l2] * y[m1][m2] * z[n1-1][n2-1];
    I2 = x[l1][l2] * y[m1][m2] * z[n1+1][n2+1];
    I3 = (n2 == 0) ? 0.0 : x[l1][l2] * y[m1][m2] * z[n1+1][n2-1];
    I4 = (n1 == 0) ? 0.0 : x[l1][l2] * y[m1][m2] * z[n1-1][n2+1];
    double Iz = 0.5 * n1 * n2 * I1 + 2.0 * a1 * a2 * I2 - a1 * n2 * I3 - n1 * a2 * I4;

    return (Ix + Iy + Iz);
}

void KineticInt::compute_pair_deriv1(const GaussianShell& s1, const GaussianShell& s2)
{
    int ao12;
    const int am1 = s1.am();
    const int am2 = s2.am();
    const int nprim1 = s1.nprimitive();
    const int nprim2 = s2.nprimitive();

    double A[3], B[3];
    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    // size of the length of a perturbation
    const size_t size = s1.ncartesian() * s2.ncartesian();
    const int center_i_start = 0;       // always 0
    const int center_j_start = 3*size;  // skip over x, y, z of center i

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, 6 * size * sizeof(double));

    double **x = overlap_recur_.x();
    double **y = overlap_recur_.y();
    double **z = overlap_recur_.z();

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1.exp(p1);
        double c1 = s1.coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2.exp(p2);
            double c2 = s2.coef(p2);
            double gamma = a1 + a2;
            double oog = 1.0/gamma;

            double PA[3], PB[3];
            double P[3];

            P[0] = (a1*A[0] + a2*B[0])*oog;
            P[1] = (a1*A[1] + a2*B[1])*oog;
            P[2] = (a1*A[2] + a2*B[2])*oog;
            PA[0] = P[0] - A[0];
            PA[1] = P[1] - A[1];
            PA[2] = P[2] - A[2];
            PB[0] = P[0] - B[0];
            PB[1] = P[1] - B[1];
            PB[2] = P[2] - B[2];

            double over_pf = exp(-a1*a2*AB2*oog) * sqrt(M_PI*oog) * M_PI * oog * c1 * c2;

            // Do recursion
            overlap_recur_.compute(PA, PB, gamma, am1+2, am2+2);

            ao12 = 0;
            for(int ii = 0; ii <= am1; ii++) {
                int l1 = am1 - ii;
                for(int jj = 0; jj <= ii; jj++) {
                    int m1 = ii - jj;
                    int n1 = jj;
                    /*--- create all am components of sj ---*/
                    for(int kk = 0; kk <= am2; kk++) {
                        int l2 = am2 - kk;
                        for(int ll = 0; ll <= kk; ll++) {
                            int m2 = kk - ll;
                            int n2 = ll;

                            double ix=0.0,iy=0.0,iz=0.0;
                            // x on center i
                            ix += 2.0 * a1 * ke_int(x, y, z, a1, l1+1, m1, n1, a2, l2, m2, n2) * over_pf;
                            if (l1)
                                ix -= l1 * ke_int(x, y, z, a1, l1-1, m1, n1, a2, l2, m2, n2) * over_pf;
                            // y on center i
                            iy += 2.0 * a1 * ke_int(x, y, z, a1, l1, m1+1, n1, a2, l2, m2, n2) * over_pf;
                            if (m1)
                                iy -= m1 * ke_int(x, y, z, a1, l1, m1-1, n1, a2, l2, m2, n2) * over_pf;
                            // z on center i
                            iz += 2.0 * a1 * ke_int(x, y, z, a1, l1, m1, n1+1, a2, l2, m2, n2) * over_pf;
                            if (n1)
                                iz -= n1 * ke_int(x, y, z, a1, l1, m1, n1-1, a2, l2, m2, n2) * over_pf;
                            // x on center i,j
                            buffer_[center_i_start+(0*size)+ao12] += ix;
                            buffer_[center_j_start+(0*size)+ao12] -= ix;
                            // y on center i,j
                            buffer_[center_i_start+(1*size)+ao12] += iy;
                            buffer_[center_j_start+(1*size)+ao12] -= iy;
                            // z on center i,j
                            buffer_[center_i_start+(2*size)+ao12] += iz;
                            buffer_[center_j_start+(2*size)+ao12] -= iz;

                            ao12++;
                        }
                    }
                }
            }
        }
    }
}

void KineticInt::compute_pair_deriv2(const GaussianShell&s1,
                                     const GaussianShell&s2)
{
    int ao12;
    int am1 = s1.am();
    int am2 = s2.am();
    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();
    double A[3], B[3];
    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    size_t size = s1.ncartesian() * s2.ncartesian();

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, 6 * s1.ncartesian() * s2.ncartesian() * sizeof(double));

    double **x = overlap_recur_.x();
    double **y = overlap_recur_.y();
    double **z = overlap_recur_.z();

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1.exp(p1);
        double c1 = s1.coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2.exp(p2);
            double c2 = s2.coef(p2);
            double gamma = a1 + a2;
            double oog = 1.0/gamma;

            double PA[3], PB[3];
            double P[3];

            P[0] = (a1*A[0] + a2*B[0])*oog;
            P[1] = (a1*A[1] + a2*B[1])*oog;
            P[2] = (a1*A[2] + a2*B[2])*oog;
            PA[0] = P[0] - A[0];
            PA[1] = P[1] - A[1];
            PA[2] = P[2] - A[2];
            PB[0] = P[0] - B[0];
            PB[1] = P[1] - B[1];
            PB[2] = P[2] - B[2];

            double over_pf = exp(-a1*a2*AB2*oog) * sqrt(M_PI*oog) * M_PI * oog * c1 * c2;

            // Do recursion
            overlap_recur_.compute(PA, PB, gamma, am1+2, am2+2);

            ao12 = 0;
            for(int ii = 0; ii <= am1; ii++) {
                int l1 = am1 - ii;
                for(int jj = 0; jj <= ii; jj++) {
                    int m1 = ii - jj;
                    int n1 = jj;
                        /*--- create all am components of sj ---*/
                    for(int kk = 0; kk <= am2; kk++) {
                        int l2 = am2 - kk;
                        for(int ll = 0; ll <= kk; ll++) {
                            int m2 = kk - ll;
                            int n2 = ll;

                            // T_{\mu\nu}^{a_x a_x}
                            buffer_[(0*size)+ao12] += (4.0*a1*a1*ke_int(x, y, z, a1, l1+2, m1, n1, a2, l2, m2, n2) -
                                                      2.0*a1*(2*l1+1)*ke_int(x, y, z, a1, l1, m1, n1, a2, l2, m2, n2)) * over_pf;
                            if (l1 > 1)
                                buffer_[(0*size)+ao12] += over_pf*l1*(l1-1)*ke_int(x, y, z, a1, l1-2, m1, n1, a2, l2, m2, n2);

                            // T_{\mu\nu}^{a_x a_y}
                            buffer_[(1*size)+ao12] += over_pf*4.0*a1*a1*ke_int(x, y, z, a1, l1+1, m1+1, n1, a2, l2, m2, n2);
                            if (l1)
                                buffer_[(1*size)+ao12] -= over_pf*2.0*l1*a1*ke_int(x, y, z, a1, l1-1, m1+1, n1, a2, l2, m2, n2);
                            if (m1)
                                buffer_[(1*size)+ao12] -= over_pf*2.0*m1*a1*ke_int(x, y, z, a1, l1+1, m1-1, n1, a2, l2, m2, n2);
                            if (l1 && m1)
                                buffer_[(1*size)+ao12] += over_pf*l1*m1*ke_int(x, y, z, a1, l1-1, m1-1, n1, a2, l2, m2, n2);

                            // T_{\mu\nu}^{a_x a_z}
                            buffer_[(2*size)+ao12] += over_pf*4.0*a1*a1*ke_int(x, y, z, a1, l1+1, m1, n1+1, a2, l2, m2, n2);
                            if (l1)
                                buffer_[(2*size)+ao12] -= over_pf*2.0*l1*a1*ke_int(x, y, z, a1, l1-1, m1, n1+1, a2, l2, m2, n2);
                            if (n1)
                                buffer_[(2*size)+ao12] -= over_pf*2.0*n1*a1*ke_int(x, y, z, a1, l1+1, m1, n1-1, a2, l2, m2, n2);
                            if (l1 && n1)
                                buffer_[(2*size)+ao12] += over_pf*l1*n1*ke_int(x, y, z, a1, l1-1, m1, n1-1, a2, l2, m2, n2);

                            // T_{\mu\nu}^{a_y a_y}
                            buffer_[(3*size)+ao12] += (4.0*a1*a1*over_pf*ke_int(x, y, z, a1, l1, m1+2, n1, a2, l2, m2, n2) -
                                                      2.0*a1*(2*m1+1)*ke_int(x, y, z, a1, l1, m1, n1, a2, l2, m2, n2)) * over_pf;
                            if (m1 > 1)
                                buffer_[(3*size)+ao12] += over_pf*m1*(m1-1)*ke_int(x, y, z, a1, l1, m1-2, n1, a2, l2, m2, n2);

                            // T_{\mu\nu}^{a_y a_z}
                            buffer_[(4*size)+ao12] += over_pf*4.0*a1*a1*ke_int(x, y, z, a1, l1, m1+1, n1+1, a2, l2, m2, n2);
                            if (m1)
                                buffer_[(4*size)+ao12] -= over_pf*2.0*m1*a1*ke_int(x, y, z, a1, l1, m1-1, n1+1, a2, l2, m2, n2);
                            if (n1)
                                buffer_[(4*size)+ao12] -= over_pf*2.0*n1*a1*ke_int(x, y, z, a1, l1, m1+1, n1-1, a2, l2, m2, n2);
                            if (m1 && n1)
                                buffer_[(4*size)+ao12] += over_pf*m1*n1*ke_int(x, y, z, a1, l1, m1-1, n1-1, a2, l2, m2, n2);

                            // T_{\mu\nu}^{a_z a_z}
                            buffer_[(5*size)+ao12] += (4.0*a1*a1*over_pf*ke_int(x, y, z, a1, l1, m1, n1+2, a2, l2, m2, n2) -
                                                      2.0*a1*(2*n1+1)*ke_int(x, y, z, a1, l1, m1, n1, a2, l2, m2, n2)) * over_pf;
                            if (n1 > 1)
                                buffer_[(5*size)+ao12] += over_pf*n1*(n1-1)*ke_int(x, y, z, a1, l1, m1, n1-2, a2, l2, m2, n2);

                            ao12++;
                        }
                    }
                }
            }
        }
    }
}
