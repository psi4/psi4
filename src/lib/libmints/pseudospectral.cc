#include <libciomr/libciomr.h>

#include "mints.h"

#include <physconst.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace boost;
using namespace psi;

// Initialize potential_recur_ to +1 basis set angular momentum
PseudospectralInt::PseudospectralInt(std::vector<SphericalTransform>& st, boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2, int deriv) :
    OneBodyAOInt(st, bs1, bs2, deriv), potential_recur_(bs1->max_am()+1, bs2->max_am()+1),
    potential_deriv_recur_(bs1->max_am()+2, bs2->max_am()+2)
{
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    use_omega_ = false;
    omega_ = 0.0;
    C_[0] = 0.0; C_[1] = 0.0; C_[2] = 0.0;
    // These are equivalent to INT_NCART
//    int maxnao1 = (maxam1+1)*(maxam1+2)/2;
//    int maxnao2 = (maxam2+1)*(maxam2+2)/2;
    int maxnao1 = INT_NCART(maxam1);
    int maxnao2 = INT_NCART(maxam2);

    if (deriv == 1) {
        maxnao1 *= 3;
        maxnao2 *= 3;
    }

    buffer_ = new double[maxnao1*maxnao2];
}

PseudospectralInt::~PseudospectralInt()
{
    delete[] buffer_;
}

void PseudospectralInt::compute_shell_deriv1(int sh1, int sh2)
{
    compute_pair_deriv1(bs1_->shell(sh1), bs2_->shell(sh2));
}

// The engine only supports segmented basis sets
void PseudospectralInt::compute_pair(const GaussianShell& s1,
                                     const GaussianShell& s2)
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

    int izm = 1;
    int iym = am1 + 1;
    int ixm = iym * iym;
    int jzm = 1;
    int jym = am2 + 1;
    int jxm = jym * jym;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, s1.ncartesian() * s2.ncartesian() * sizeof(double));

    double ***vi = potential_recur_.vi();

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1.exp(p1);
        double c1 = s1.coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2.exp(p2);
            double c2 = s2.coef(p2);

            // The GPT gamma of the two cartesian functions.
            // This is used for all GPT operations here, such as over_pf and P,
            // And for OS relations AFTER the modified (0|A|0)^(m) integrals are built
            double gamma0 = a1 + a2;
            double oog = 1.0/gamma0;

            // An effective gamma if range-separation is to be used.
            // This gamma is only for use in the generation of auxiliary integrals,
            // particularly the (0|A|0)^(m) auxiliary integrals as built in potential_recur_.compute().
            double gamma = gamma0;
            if (use_omega_) {
                gamma = gamma0 * omega_ * omega_ / (gamma0 + omega_ * omega_);
            }

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

            // Loop over atoms of basis set 1 (only works if bs1_ and bs2_ are on the same
            // molecule)
            double PC[3];

            // C_ is the pseudospectral grid point
            PC[0] = P[0] - C_[0];
            PC[1] = P[1] - C_[1];
            PC[2] = P[2] - C_[2];

            // Do recursion
            potential_recur_.compute_erf(PA, PB, PC, gamma0, am1, am2, gamma);

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

                            // Compute location in the recursion
                            int iind = l1 * ixm + m1 * iym + n1 * izm;
                            int jind = l2 * jxm + m2 * jym + n2 * jzm;

                            buffer_[ao12++] += vi[iind][jind][0] * over_pf;

//                            fprintf(outfile, "ao12=%d, vi[%d][%d][0] = %20.14f, over_pf = %20.14f, Z = %f\n", ao12-1, iind, jind, vi[iind][jind][0], over_pf, Z);
                        }
                    }
                }
            }
        }
    }
}

// The engine only supports segmented basis sets
void PseudospectralInt::compute_pair_deriv1(const GaussianShell& s1, const GaussianShell& s2)
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
    int center_i = s1.ncenter()*3*size;
    int center_j = s2.ncenter()*3*size;

    int izm1 = 1;
    int iym1 = am1 + 1 + 1;  // extra 1 for derivative
    int ixm1 = iym1 * iym1;
    int jzm1 = 1;
    int jym1 = am2 + 1 + 1;  // extra 1 for derivative
    int jxm1 = jym1 * jym1;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, 3 * s1.ncartesian() * s2.ncartesian() * sizeof(double));

    double ***vi = potential_deriv_recur_.vi();
    double ***vx = potential_deriv_recur_.vx();
    double ***vy = potential_deriv_recur_.vy();
    double ***vz = potential_deriv_recur_.vz();

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1.exp(p1);
        double c1 = s1.coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2.exp(p2);
            double c2 = s2.coef(p2);

            // The GPT gamma of the two cartesian functions.
            // This is used for all GPT operations here, such as over_pf and P,
            // And for OS relations AFTER the modified (0|A|0)^(m) integrals are built
            double gamma0 = a1 + a2;
            double oog = 1.0/gamma0;

            // An effective gamma if range-separation is to be used.
            // This gamma is only for use in the generation of auxiliary integrals,
            // particularly the (0|A|0)^(m) auxiliary integrals as built in potential_recur_.compute().
            double gamma = gamma0;
            if (use_omega_) {
                gamma = gamma0 * omega_ * omega_ / (gamma0 + omega_ * omega_);
            }

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

            // Loop over atoms of basis set 1 (only works if bs1_ and bs2_ are on the same
            // molecule)
            double PC[3];

            PC[0] = P[0] - C_[0];
            PC[1] = P[1] - C_[1];
            PC[2] = P[2] - C_[2];

            // Do recursion
            potential_deriv_recur_.compute_erf(PA, PB, PC, gamma0, am1+1, am2+1, gamma);

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

                            // Compute location in the recursion
                            int iind = l1 * ixm1 + m1 * iym1 + n1 * izm1;
                            int jind = l2 * jxm1 + m2 * jym1 + n2 * jzm1;

                            const double pfac = over_pf;

                            // x
                            double temp = 2.0*a1*vi[iind+ixm1][jind][0];
                            if (l1)
                                temp -= l1*vi[iind-ixm1][jind][0];
                            buffer_[center_i+(0*size)+ao12] -= temp * pfac;
                            // printf("ix temp = %f ", temp);

                            temp = 2.0*a2*vi[iind][jind+jxm1][0];
                            if (l2)
                                temp -= l2*vi[iind][jind-jxm1][0];
                            buffer_[center_j+(0*size)+ao12] -= temp * pfac;
                            // printf("jx temp = %f ", temp);

                            buffer_[3*size+ao12] -= vx[iind][jind][0] * pfac;

                            // y
                            temp = 2.0*a1*vi[iind+iym1][jind][0];
                            if (m1)
                                temp -= m1*vi[iind-iym1][jind][0];
                            buffer_[center_i+(1*size)+ao12] -= temp * pfac;
                            // printf("iy temp = %f ", temp);

                            temp = 2.0*a2*vi[iind][jind+jym1][0];
                            if (m2)
                                temp -= m2*vi[iind][jind-jym1][0];
                            buffer_[center_j+(1*size)+ao12] -= temp * pfac;
                            // printf("jy temp = %f ", temp);

                            buffer_[3*size+size+ao12] -= vy[iind][jind][0] * pfac;

                            // z
                            temp = 2.0*a1*vi[iind+izm1][jind][0];
                            if (n1)
                                temp -= n1*vi[iind-izm1][jind][0];
                            buffer_[center_i+(2*size)+ao12] -= temp * pfac;
                            // printf("iz temp = %f ", temp);

                            temp = 2.0*a2*vi[iind][jind+jzm1][0];
                            if (n2)
                                temp -= n2*vi[iind][jind-jzm1][0];
                            buffer_[center_j+(2*size)+ao12] -= temp * pfac;
                            // printf("jz temp = %f \n", temp);

                            buffer_[3*size+2*size+ao12] -= vz[iind][jind][0] * pfac;

                            ao12++;
                        }
                    }
                }
            }
        }
    }

    // Integrals are done. Normalize for angular momentum
    normalize_am(s1, s2);

    // Spherical harmonic transformation
    // Wrapped up in the AO to SO transformation (I think)
}
