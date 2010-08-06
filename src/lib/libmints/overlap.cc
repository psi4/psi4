#include <stdexcept>
#include <libciomr/libciomr.h>

#include "mints.h"
#include <physconst.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace psi;

OverlapInt::OverlapInt(std::vector<SphericalTransform>& st, shared_ptr<BasisSet> bs1, shared_ptr<BasisSet> bs2, int deriv) :
    OneBodyInt(st, bs1, bs2, deriv), overlap_recur_(bs1->max_am()+deriv, bs2->max_am()+deriv)
{
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    if (deriv > 2) {
        throw std::runtime_error("OverlapInt: does not support 3rd order derivatives and higher.");
    }
    int maxnao1 = (maxam1+1)*(maxam1+2)/2;
    int maxnao2 = (maxam2+1)*(maxam2+2)/2;
    if (deriv == 1) {
        maxnao1 *= 3 * natom_;
        maxnao2 *= 3 * natom_;
    } else {
        maxnao1 *= 9 * natom_;
        maxnao2 *= 9 * natom_;
    }
    buffer_ = new double[maxnao1*maxnao2];
}

OverlapInt::~OverlapInt()
{
    delete[] buffer_;
}

void OverlapInt::compute_shell(int sh1, int sh2)
{
    compute_pair(bs1_->shell(sh1), bs2_->shell(sh2));
}

// The engine only supports segmented basis sets
void OverlapInt::compute_pair(shared_ptr<GaussianShell> s1, shared_ptr<GaussianShell> s2)
{
    int ao12;
    int am1 = s1->am();
    int am2 = s2->am();
    int nprim1 = s1->nprimitive();
    int nprim2 = s2->nprimitive();
    double A[3], B[3];
    A[0] = s1->center()[0];
    A[1] = s1->center()[1];
    A[2] = s1->center()[2];
    B[0] = s2->center()[0];
    B[1] = s2->center()[1];
    B[2] = s2->center()[2];

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, s1->ncartesian() * s2->ncartesian() * sizeof(double));

    double **x = overlap_recur_.x();
    double **y = overlap_recur_.y();
    double **z = overlap_recur_.z();

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1->exp(p1);
        double c1 = s1->coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2->exp(p2);
            double c2 = s2->coef(p2);
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

            //printf("PA %f %f %f PB %f %f %f gamma %f am1 %d am2 %d\n", PA[0], PA[1], PA[2], PB[0], PB[1], PB[2], gamma, am1, am2);

            // Do recursion
            overlap_recur_.compute(PA, PB, gamma, am1, am2);

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

                            double x0 = x[l1][l2];
                            double y0 = y[m1][m2];
                            double z0 = z[n1][n2];

                            buffer_[ao12++] += over_pf*x0*y0*z0;
                        }
                    }
                }
            }
        }
    }

    // Integrals are done. Normalize for angular momentum
    normalize_am(s1, s2);

    // Spherical harmonic transformation
    // Wrapped up in the AO to SO transformation
}

void OverlapInt::compute_shell_deriv1(int sh1, int sh2)
{
    compute_pair_deriv1(bs1_->shell(sh1), bs2_->shell(sh2));
}

void OverlapInt::compute_pair_deriv1(shared_ptr<GaussianShell> s1, shared_ptr<GaussianShell> s2)
{
    int ao12;
    int am1 = s1->am();
    int am2 = s2->am();
    int nprim1 = s1->nprimitive();
    int nprim2 = s2->nprimitive();
    double A[3], B[3];
    A[0] = s1->center()[0];
    A[1] = s1->center()[1];
    A[2] = s1->center()[2];
    B[0] = s2->center()[0];
    B[1] = s2->center()[1];
    B[2] = s2->center()[2];

    size_t size = s1->ncartesian() * s2->ncartesian();
    int center_i = s1->ncenter()*3*size;
    int center_j = s2->ncenter()*3*size;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, 3 * natom_ * s1->ncartesian() * s2->ncartesian() * sizeof(double));

    double **x = overlap_recur_.x();
    double **y = overlap_recur_.y();
    double **z = overlap_recur_.z();

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1->exp(p1);
        double c1 = s1->coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2->exp(p2);
            double c2 = s2->coef(p2);
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

                            // x on center i
                            buffer_[center_i+(0*size)+ao12] += 2.0*a1*over_pf*x[l1+1][l2]*y[m1][m2]*z[n1][n2];
                            if (l1)
                                buffer_[center_i+(0*size)+ao12] -= over_pf*x[l1-1][l2]*y[m1][m2]*z[n1][n2];
                            // y on center i
                            buffer_[center_i+(1*size)+ao12] += 2.0*a1*over_pf*x[l1][l2]*y[m1+1][m2]*z[n1][n2];
                            if (m1)
                                buffer_[center_i+(1*size)+ao12] -= over_pf*x[l1][l2]*y[m1-1][m2]*z[n1][n2];
                            // z on center i
                            buffer_[center_i+(2*size)+ao12] += 2.0*a1*over_pf*x[l1][l2]*y[m1][m2]*z[n1+1][n2];
                            if (n1)
                                buffer_[center_i+(2*size)+ao12] -= over_pf*x[l1][l2]*y[m1][m2]*z[n1-1][n2];

                            // x on center j
                            buffer_[center_j+(0*size)+ao12] += 2.0*a2*over_pf*x[l1][l2+1]*y[m1][m2]*z[n1][n2];
                            if (l2)
                                buffer_[center_j+(0*size)+ao12] -= over_pf*x[l1][l2-1]*y[m1][m2]*z[n1][n2];
                            // y on center j
                            buffer_[center_j+(1*size)+ao12] += 2.0*a2*over_pf*x[l1][l2]*y[m1][m2+1]*z[n1][n2];
                            if (m2)
                                buffer_[center_j+(1*size)+ao12] -= over_pf*x[l1][l2]*y[m1][m2-1]*z[n1][n2];
                            // z on center j
                            buffer_[center_j+(2*size)+ao12] += 2.0*a2*over_pf*x[l1][l2]*y[m1][m2]*z[n1][n2+1];
                            if (n2)
                                buffer_[center_j+(2*size)+ao12] -= over_pf*x[l1][l2]*y[m1][m2]*z[n1][n2-1];

                            ao12++;
                        }
                    }
                }
            }
        }
    }

    // Integrals are done. Normalize for angular momentum
    normalize_am(s1, s2, 3*natom_);

    // Spherical harmonic transformation
    // Wrapped up in the AO to SO transformation
}
