#include <boost/shared_ptr.hpp>
#include <stdexcept>
#include <libciomr/libciomr.h>

#include "mints.h"
#include "angularmomentum.h"
#include <physconst.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace psi;
using namespace boost;

// Initialize overlap_recur_ to +1 basis set angular momentum, +1 on each center is sufficient
// to compute the dipole derivatives
AngularMomentumInt::AngularMomentumInt(std::vector<SphericalTransform>& spherical_transforms, boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2, int nderiv) :
    OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv), overlap_recur_(bs1->max_am()+1, bs2->max_am()+1)
{
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = (maxam1+1)*(maxam1+2)/2;
    int maxnao2 = (maxam2+1)*(maxam2+2)/2;

    // Increase buffer size to handle x, y, and z components
    if (deriv_ == 0) {
        buffer_ = new double[3*maxnao1*maxnao2];
        set_chunks(3);
    }
    else if (deriv_ == 1) {
        natom_ = bs1_->molecule()->natom();
        buffer_ = new double[3*3*natom_*maxnao1*maxnao2]; // 3 * 3 * N * maxnao1 * maxnao2
        set_chunks(3*3*natom_);
    }
}

AngularMomentumInt::~AngularMomentumInt()
{
    delete[] buffer_;
}

// The engine only supports segmented basis sets
void AngularMomentumInt::compute_pair(const GaussianShell& s1, const GaussianShell& s2)
{
    int ao12;
    int am1 = s1.am();
    int am2 = s2.am();
    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();
    double A[3], B[3];
    double Sxy1, Sxy2, Sxz1, Sxz2;
    double Syx1, Syx2, Syz1, Syz2;
    double Szx1, Szx2, Szy1, Szy2;
    double S0x1, S0x2, S0y1, S0y2, S0z1, S0z2;
    double muxy1, muxy2, muxz1, muxz2;
    double muyx1, muyx2, muyz1, muyz2;
    double muzx1, muzx2, muzy1, muzy2;

    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    int ydisp = INT_NCART(am1) * INT_NCART(am2);
    int zdisp = ydisp + INT_NCART(am1) * INT_NCART(am2);

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, 3 * INT_NCART(am1) * INT_NCART(am2) * sizeof(double));

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

                            double x00 = x[l1][l2],   y00 = y[m1][m2],   z00 = z[n1][n2];
                            double x10 = x[l1+1][l2], y10 = y[m1+1][m2], z10 = z[n1+1][n2];
                            double x01 = x[l1][l2+1], y01 = y[m1][m2+1], z01 = z[n1][n2+1];

                            Sxy1 = Sxy2 = Sxz1 = Sxz2 = 0.0;

                            //
                            // Overlaps
                            //

                            // (a+1x|b+1y)
                            Sxy1 = x10 * y01 * z00 * over_pf;
                            // (a+1x|b-1y)
                            if (m2)
                                Sxy2 = x10 * y[m1][m2-1] * z00 * over_pf;
                            // (a+1x|b+1z)
                            Sxz1 = x10 * y00 * z01 * over_pf;
                            // (a+1x|b-1z)
                            if (n2)
                                Sxz2 = x10 * y00 * z[n1][n2-1] * over_pf;

//                            fprintf(outfile, "Sxy1 %f Sxy2 %f Sxz1 %f Sxz2 %f\n", Sxy1, Sxy2, Sxz1, Sxz2);

                            Syx1 = Syx2 = Syz1 = Syz2 = 0.0;

                            // (a+1y|b+1x)
                            Syx1 = x01 * y10 * z00 * over_pf;
                            // (a+1y|b-1x)
                            if (l2)
                                Syx2 = x[l1][l2-1] * y10 * z00 * over_pf;
                            // (a+1y|b+1z)
                            Syz1 = x00 * y10 * z01 * over_pf;
                            // (a+1y|b-1z)
                            if (n2)
                                Syz2 = x00 * y10 * z[n1][n2-1] * over_pf;

//                            fprintf(outfile, "Syx1 %f Syx2 %f Syz1 %f Syz2 %f\n", Syx1, Syx2, Syz1, Syz2);

                            Szx1 = Szx2 = Szy1 = Szy2 = 0.0;

                            // (a+1z|b+1x)
                            Szx1 = x01 * y00 * z10 * over_pf;
                            // (a+1z|b-1x)
                            if (l2)
                                Szx2 = x[l1][l2-1] * y00 * z10 * over_pf;
                            // (a+1z|b+1y)
                            Szy1 = x00 * y01 * z10 * over_pf;
                            // (a+1z|b-1y)
                            if (m2)
                                Szy2 = x00 * y[m1][m2-1] * z10 * over_pf;

//                            fprintf(outfile, "Szx1 %f Szx2 %f Szy1 %f Szy2 %f\n", Szx1, Szx2, Szy1, Szy2);

                            S0x1 = S0x2 = S0y1 = S0y2 = S0z1 = S0z2 = 0.0;

                            // (a|b+1x)
                            S0x1 = x01 * y00 * z00 * over_pf;
                            // (a|b-1x)
                            if (l2)
                                S0x2 = x[l1][l2-1] * y00 * z00 * over_pf;
                            // (a|b+1y)
                            S0y1 = x00 * y01 * z00 * over_pf;
                            // (a|b-1y)
                            if (m2)
                                S0y2 = x00 * y[m1][m2-1] * z00 * over_pf;
                            // (a|b+1z)
                            S0z1 = x00 * y00 * z01 * over_pf;
                            // (a|b-1z)
                            if (n2)
                                S0z2 = x00 * y00 * z[n1][n2-1] * over_pf;

//                            fprintf(outfile, "S0x1 %f S0x2 %f S0y1 %f S0y2 %f S0z1 S0z2\n", S0x1, S0x2, S0y1, S0y2, S0z1, S0z2);

                            //
                            // Moment integrals
                            //

                            muxy1 = muxy2 = muxz1 = muxz2 = 0.0;
                            // (a|x|b+1y)
                            muxy1 = Sxy1 + (A[0] - origin_[0]) * S0y1;
                            // (a|x|b+1z)
                            muxz1 = Sxz1 + (A[0] - origin_[0]) * S0z1;
                            // (a|x|b-1y)
                            muxy2 = Sxy2 + (A[0] - origin_[0]) * S0y2;
                            // (a|x|b-1z)
                            muxz2 = Sxz2 + (A[0] - origin_[0]) * S0z2;

                            muyx1 = muyx2 = muyz1 = muyz2 = 0.0;
                            // (a|y|b+1x)
                            muyx1 = Syx1 + (A[1] - origin_[1]) * S0x1;
                            // (a|y|b+1z)
                            muyz1 = Syz1 + (A[1] - origin_[1]) * S0z1;
                            // (a|y|b-1x)
                            muyx2 = Syx2 + (A[1] - origin_[1]) * S0x2;
                            // (a|y|b-1z)
                            muyz2 = Syz2 + (A[1] - origin_[1]) * S0z2;

                            muzx1 = muzx2 = muzy1 = muzy2 = 0.0;
                            // (a|z|b+1x)
                            muzx1 = Szx1 + (A[2] - origin_[2]) * S0x1;
                            // (a|z|b+1y)
                            muzy1 = Szy1 + (A[2] - origin_[2]) * S0y1;
                            // (a|z|b+1x)
                            muzx2 = Szx2 + (A[2] - origin_[2]) * S0x2;
                            // (a|z|b+1y)
                            muzy2 = Szy2 + (A[2] - origin_[2]) * S0y2;

                            /* (a|Lx|b) = 2 a2 * (a|(y-Cy)|b+1z) - B.z * (a|(y-Cy)|b-1z)
                               - 2 a2 * (a|(z-Cz)|b+1y) + B.y * (a|(z-Cz)|b-1y) */
                            double Lx = (2.0*a2*muyz1 - n2*muyz2 - 2.0*a2*muzy1 + m2*muzy2)/* * over_pf*/;

                            /* (a|Ly|b) = 2 a2 * (a|(z-Cz)|b+1x) - B.x * (a|(z-Cz)|b-1x)
                               - 2 a2 * (a|(x-Cx)|b+1z) + B.z * (a|(x-Cx)|b-1z) */
                            double Ly = (2.0*a2*muzx1 - l2*muzx2 - 2.0*a2*muxz1 + n2*muxz2)/* * over_pf*/;

                            /* (a|Lz|b) = 2 a2 * (a|(x-Cx)|b+1y) - B.y * (a|(x-Cx)|b-1y)
                               - 2 a2 * (a|(y-Cy)|b+1x) + B.x * (a|(y-Cy)|b-1x) */
                            double Lz = (2.0*a2*muxy1 - m2*muxy2 - 2.0*a2*muyx1 + l2*muyx2)/* * over_pf*/;

                            // Electrons have a negative charge
                            buffer_[ao12]       += (Lx);
                            buffer_[ao12+ydisp] += (Ly);
                            buffer_[ao12+zdisp] += (Lz);

                            ao12++;
                        }
                    }
                }
            }
        }
    }
}

#if 0
// The engine only supports segmented basis sets
void AngularMomentumInt::compute_pair_deriv1(const GaussianShell& s1, const GaussianShell& s2)
{
    int ao12;
    int am1 = s1->am();
    int am2 = s2->am();
    int at1 = s1->ncenter();
    int at2 = s2->ncenter();
    int nprim1 = s1->nprimitive();
    int nprim2 = s2->nprimitive();
    size_t length = INT_NCART(am1) * INT_NCART(am2);
    double A[3], B[3];

    A[0] = s1->center()[0];
    A[1] = s1->center()[1];
    A[2] = s1->center()[2];
    B[0] = s2->center()[0];
    B[1] = s2->center()[1];
    B[2] = s2->center()[2];

    size_t xaxdisp = at1 * length * 9;
    size_t xaydisp = xaxdisp + length;
    size_t xazdisp = xaydisp + length;
    size_t yaxdisp = xazdisp + length;
    size_t yaydisp = yaxdisp + length;
    size_t yazdisp = yaydisp + length;
    size_t zaxdisp = yazdisp + length;
    size_t zaydisp = zaxdisp + length;
    size_t zazdisp = zaydisp + length;

    size_t xbxdisp = at2 * length * 9;
    size_t xbydisp = xbxdisp + length;
    size_t xbzdisp = xbydisp + length;
    size_t ybxdisp = xbzdisp + length;
    size_t ybydisp = ybxdisp + length;
    size_t ybzdisp = ybydisp + length;
    size_t zbxdisp = ybzdisp + length;
    size_t zbydisp = zbxdisp + length;
    size_t zbzdisp = zbydisp + length;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    // Zero out the buffer
    memset(buffer_, 0, 3 * 3 * natom_ * length * sizeof(double));

    double **x = overlap_recur_.x();
    double **y = overlap_recur_.y();
    double **z = overlap_recur_.z();
    double v1, v2, v3, v4;      // temporary value storage

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

            // Do recursion, this is sufficient information to compute dipole derivatives
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

                            // mu-x derivatives
                            v1 = v2 = v3 = v4 = 0.0;

                            //
                            // A derivatives with mu-x
                            //

                            // (a+1_x|mx|b+1_x)
                            v1 = x[l1+1][l2+1] * y[m1][m2] * z[n1][n2];
                            // (a+1_x|mx|b)
                            v2 = x[l1+1][l2]   * y[m1][m2] * z[n1][n2];
                            if (l1) {
                                // (a-1_x|mx|b+1_x)
                                v3 = x[l1-1][l2+1] * y[m1][m2] * z[n1][n2];
                                // (a-1_x|mx|b)
                                v4 = x[l1-1][l2]   * y[m1][m2] * z[n1][n2];
                            }
                            buffer_[ao12+xaxdisp] -= (2.0 * a1 * (v1 + B[0] * v2) - l1 * (v3 + B[0] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_y|mx|b+1_x)
                            v1 = x[l1][l2+1] * y[m1+1][m2] * z[n1][n2];
                            // (a+1_y|mx|b)
                            v2 = x[l1][l2]   * y[m1+1][m2] * z[n1][n2];
                            if (m1) {
                                // (a-1_y|mx|b+1_x)
                                v3 = x[l1][l2+1] * y[m1-1][m2] * z[n1][n2];
                                // (a-1_y|mx|b)
                                v4 = x[l1][l2]   * y[m1-1][m2] * z[n1][n2];
                            }
                            buffer_[ao12+xaydisp] -= (2.0 * a1 * (v1 + B[0] * v2) - m1 * (v3 + B[0] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_z|mx|b+1_x)
                            v1 = x[l1][l2+1] * y[m1][m2] * z[n1+1][n2];
                            // (a+1_z|mx|b)
                            v2 = x[l1][l2]   * y[m1][m2] * z[n1+1][n2];
                            if (n1) {
                                // (a-1_z|mx|b+1_x)
                                v3 = x[l1][l2+1] * y[m1][m2] * z[n1-1][n2];
                                // (a-1_z|mx|b)
                                v4 = x[l1][l2]   * y[m1][m2] * z[n1-1][n2];
                            }
                            buffer_[ao12+xazdisp] -= (2.0 * a1 * (v1 + B[0] * v2) - n1 * (v3 + B[0] * v4)) * over_pf;

                            //
                            // B derivatives with mu-x
                            //

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_x|mx|b+1_x)
                            v1 = x[l1+1][l2+1] * y[m1][m2] * z[n1][n2];
                            // (a|mx|b+1_x)
                            v2 = x[l1][l2+1]   * y[m1][m2] * z[n1][n2];
                            if (l2) {
                                // (a+1_x|mx|b-1_x)
                                v3 = x[l1+1][l2-1] * y[m1][m2] * z[n1][n2];
                                // (a|mx|b-1_x)
                                v4 = x[l1][l2-1]   * y[m1][m2] * z[n1][n2];
                            }
                            buffer_[ao12+xbxdisp] -= (2.0 * a2 * (v1 + A[0] * v2) - l2 * (v3 + A[0] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_x|mx|b+1_y)
                            v1 = x[l1+1][l2] * y[m1][m2+1] * z[n1][n2];
                            // (a|mx|b+1_y)
                            v2 = x[l1][l2]   * y[m1][m2+1] * z[n1][n2];
                            if (m2) {
                                // (a+1_x|mx|b-1_y)
                                v3 = x[l1+1][l2] * y[m1][m2-1] * z[n1][n2];
                                // (a|mx|b-1_y)
                                v4 = x[l1][l2]   * y[m1][m2-1] * z[n1][n2];
                            }
                            buffer_[ao12+xbydisp] -= (2.0 * a2 * (v1 + A[0] * v2) - m2 * (v3 + A[0] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_x|mx|b+1_z)
                            v1 = x[l1+1][l2] * y[m1][m2] * z[n1][n2+1];
                            // (a|mx|b+1_z)
                            v2 = x[l1][l2]   * y[m1][m2] * z[n1][n2+1];
                            if (n2) {
                                // (a+1_x|mx|b-1_z)
                                v3 = x[l1+1][l2] * y[m1][m2] * z[n1][n2-1];
                                // (a|mx|b-1_z)
                                v4 = x[l1][l2]   * y[m1][m2] * z[n1][n2-1];
                            }
                            buffer_[ao12+xbzdisp] -= (2.0 * a2 * (v1 + A[0] * v2) - n2 * (v3 + A[0] * v4)) * over_pf;

                            //
                            // A derivatives with mu-y
                            //

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_x|my|b+1_y)
                            v1 = x[l1+1][l2] * y[m1][m2+1] * z[n1][n2];
                            // (a+1_x|my|b)
                            v2 = x[l1+1][l2] * y[m1][m2] * z[n1][n2];
                            if (l1) {
                                // (a-1_x|my|b+1_y)
                                v3 = x[l1-1][l2] * y[m1][m2+1] * z[n1][n2];
                                // (a-1_x|my|b)
                                v4 = x[l1-1][l2] * y[m1][m2] * z[n1][n2];
                            }
                            buffer_[ao12+yaxdisp] -= (2.0 * a1 * (v1 + B[1] * v2) - l1 * (v3 + B[1] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_y|my|b+1_y)
                            v1 = x[l1][l2] * y[m1+1][m2+1] * z[n1][n2];
                            // (a+1_y|my|b)
                            v2 = x[l1][l2] * y[m1+1][m2] * z[n1][n2];
                            if (m1) {
                                // (a-1_y|my|b+1_y)
                                v3 = x[l1][l2] * y[m1-1][m2+1] * z[n1][n2];
                                // (a-1_y|my|b)
                                v4 = x[l1][l2] * y[m1-1][m2] * z[n1][n2];
                            }
                            buffer_[ao12+yaydisp] -= (2.0 * a1 * (v1 + B[1] * v2) - m1 * (v3 + B[1] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_z|my|b+1_y)
                            v1 = x[l1][l2] * y[m1][m2+1] * z[n1+1][n2];
                            // (a+1_z|my|b)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1+1][n2];
                            if (n1) {
                                // (a-1_z|my|b+1_y)
                                v3 = x[l1][l2] * y[m1][m2+1] * z[n1-1][n2];
                                // (a-1_z|my|b)
                                v4 = x[l1][l2] * y[m1][m2] * z[n1-1][n2];
                            }
                            buffer_[ao12+yazdisp] -= (2.0 * a1 * (v1 + B[1] * v2) - n1 * (v3 + B[1] * v4)) * over_pf;

                            //
                            // B derivatives with mu-y
                            //

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_y|my|b+1_x)
                            v1 = x[l1][l2+1] * y[m1+1][m2] * z[n1][n2];
                            // (a|my|b+1_x)
                            v2 = x[l1][l2+1] * y[m1][m2] * z[n1][n2];
                            if (l2) {
                                // (a+1_y|my|b-1_x)
                                v3 = x[l1][l2-1] * y[m1+1][m2] * z[n1][n2];
                                // (a|my|b-1_x)
                                v4 = x[l1][l2-1] * y[m1][m2] * z[n1][n2];
                            }
                            buffer_[ao12+ybxdisp] -= (2.0 * a2 * (v1 + A[1] * v2) - l2 * (v3 + A[1] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_y|my|b+1_y)
                            v1 = x[l1][l2] * y[m1+1][m2+1] * z[n1][n2];
                            // (a|my|b+1_y)
                            v2 = x[l1][l2] * y[m1][m2+1] * z[n1][n2];
                            if (m2) {
                                // (a+1_y|my|b-1_y)
                                v3 = x[l1][l2] * y[m1+1][m2-1] * z[n1][n2];
                                // (a|my|b-1_y)
                                v4 = x[l1][l2] * y[m1][m2-1] * z[n1][n2];
                            }
                            buffer_[ao12+ybydisp] -= (2.0 * a2 * (v1 + A[1] * v2) - m2 * (v3 + A[1] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_y|my|b+1_z)
                            v1 = x[l1][l2] * y[m1+1][m2] * z[n1][n2+1];
                            // (a|my|b+1_z)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1][n2+1];
                            if (n2) {
                                // (a+1_y|my|b-1_z)
                                v3 = x[l1][l2] * y[m1+1][m2] * z[n1][n2-1];
                                // (a|my|b-1_z)
                                v4 = x[l1][l2] * y[m1][m2] * z[n1][n2-1];
                            }
                            buffer_[ao12+ybzdisp] -= (2.0 * a2 * (v1 + A[1] * v2) - n2 * (v3 + A[1] * v4)) * over_pf;

                            //
                            // A derivatives with mu-z
                            //

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_x|mz|b+1_z)
                            v1 = x[l1+1][l2] * y[m1][m2] * z[n1][n2+1];
                            // (a+1_x|mz|b)
                            v2 = x[l1+1][l2] * y[m1][m2] * z[n1][n2];
                            if (l1) {
                                // (a-1_x|mz|b+1_z)
                                v3 = x[l1-1][l2] * y[m1][m2] * z[n1][n2+1];
                                // (a-1_x|mz|b)
                                v4 = x[l1-1][l2] * y[m1][m2] * z[n1][n2];
                            }
                            buffer_[ao12+zaxdisp] -= (2.0 * a1 * (v1 + B[2] * v2) - l1 * (v3 + B[2] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_y|mz|b+1_z)
                            v1 = x[l1][l2] * y[m1+1][m2] * z[n1][n2+1];
                            // (a+1_y|mz|b)
                            v2 = x[l1][l2] * y[m1+1][m2] * z[n1][n2];
                            if (m1) {
                                // (a-1_y|mz|b+1_z)
                                v3 = x[l1][l2] * y[m1-1][m2] * z[n1][n2+1];
                                // (a-1_y|mz|b)
                                v4 = x[l1][l2] * y[m1-1][m2] * z[n1][n2];
                            }
                            buffer_[ao12+zaydisp] -= (2.0 * a1 * (v1 + B[2] * v2) - m1 * (v3 + B[2] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_z|mz|b+1_z)
                            v1 = x[l1][l2] * y[m1][m2] * z[n1+1][n2+1];
                            // (a+1_z|mz|b)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1+1][n2];
                            if (n1) {
                                // (a-1_z|mz|b+1_z)
                                v3 = x[l1][l2] * y[m1][m2] * z[n1-1][n2+1];
                                // (a-1_z|mz|b)
                                v4 = x[l1][l2] * y[m1][m2] * z[n1-1][n2];
                            }
                            buffer_[ao12+zazdisp] -= (2.0 * a1 * (v1 + B[2] * v2) - n1 * (v3 + B[2] * v4)) * over_pf;

                            //
                            // B derivates with mu-z
                            //

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_z|mz|b+1_x)
                            v1 = x[l1][l2+1] * y[m1][m2] * z[n1+1][n2];
                            // (a|mz|b+1_x)
                            v2 = x[l1][l2+1] * y[m1][m2] * z[n1][n2];
                            if (l2) {
                                // (a+1_z|mz|b-1_x)
                                v3 = x[l1][l2-1] * y[m1][m2] * z[n1+1][n2];
                                // (a|mz|b-1_x)
                                v4 = x[l1][l2-1] * y[m1][m2] * z[n1][n2];
                            }
                            buffer_[ao12+zbxdisp] -= (2.0 * a2 * (v1 + A[2] * v2) - l2 * (v3 + A[2] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_z|mz|b+1_y)
                            v1 = x[l1][l2] * y[m1][m2+1] * z[n1+1][n2];
                            // (a|mz|b+1_y)
                            v2 = x[l1][l2] * y[m1][m2+1] * z[n1][n2];
                            if (m2) {
                                // (a+1_z|mz|b-1_y)
                                v3 = x[l1][l2] * y[m1][m2-1] * z[n1+1][n2];
                                // (a|mz|b-1_y)
                                v4 = x[l1][l2] * y[m1][m2-1] * z[n1][n2];
                            }
                            buffer_[ao12+zbydisp] -= (2.0 * a2 * (v1 + A[2] * v2) - m2 * (v3 + A[2] * v4)) * over_pf;

                            v1 = v2 = v3 = v4 = 0.0;
                            // (a+1_z|mz|b+1_z)
                            v1 = x[l1][l2] * y[m1][m2] * z[n1+1][n2+1];
                            // (a|mz|b+1_z)
                            v2 = x[l1][l2] * y[m1][m2] * z[n1][n2+1];
                            if (n2) {
                                // (a+1_z|mz|b-1_z)
                                v3 = x[l1][l2] * y[m1][m2] * z[n1+1][n2-1];
                                // (a|mz|b-1_z)
                                v4 = x[l1][l2] * y[m1][m2] * z[n1][n2-1];
                            }
                            buffer_[ao12+zbzdisp] -= (2.0 * a2 * (v1 + A[2] * v2) - n2 * (v3 + A[2] * v4)) * over_pf;

                            ao12++;
                        }
                    }
                }
            }
        }
    }
}
#endif
