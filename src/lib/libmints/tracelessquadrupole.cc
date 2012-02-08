#include <libciomr/libciomr.h>

#include "mints.h"
#include <physconst.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace boost;
using namespace psi;

// Initialize overlap_recur_ to +2 basis set angular momentum
TracelessQuadrupoleInt::TracelessQuadrupoleInt(std::vector<SphericalTransform>& st, boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2) :
    OneBodyAOInt(st, bs1, bs2), overlap_recur_(bs1->max_am()+2, bs2->max_am()+2)
{
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = (maxam1+1)*(maxam1+2)/2;
    int maxnao2 = (maxam2+1)*(maxam2+2)/2;

    // Increase buffer size to handle xx, xy, xz, yy, yz, zz components
    buffer_ = new double[6*maxnao1*maxnao2];
    set_chunks(6);
}

TracelessQuadrupoleInt::~TracelessQuadrupoleInt()
{
    delete[] buffer_;
}

void TracelessQuadrupoleInt::compute_pair(const GaussianShell& s1,
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

    int xydisp = INT_NCART(am1) * INT_NCART(am2);
    int xzdisp = 2 * INT_NCART(am1) * INT_NCART(am2);
    int yydisp = 3 * INT_NCART(am1) * INT_NCART(am2);
    int yzdisp = 4 * INT_NCART(am1) * INT_NCART(am2);
    int zzdisp = 5 * INT_NCART(am1) * INT_NCART(am2);

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

                            double x00 = x[l1][l2];
                            double y00 = y[m1][m2];
                            double z00 = z[n1][n2];
                            double x01 = x[l1][l2+1];
                            double y01 = y[m1][m2+1];
                            double z01 = z[n1][n2+1];
                            double x10 = x[l1+1][l2];
                            double y10 = y[m1+1][m2];
                            double z10 = z[n1+1][n2];
                            double x11 = x[l1+1][l2+1];
                            double y11 = y[m1+1][m2+1];
                            double z11 = z[n1+1][n2+1];

                            double mxx = -over_pf*(x11 + x10*(B[0] - origin_[0]) + x01*(A[0] - origin_[0]) + x00*(A[0] - origin_[0])*(B[0] - origin_[0]))*y00*z00;
                            double myy = -over_pf*(y11 + y10*(B[1] - origin_[1]) + y01*(A[1] - origin_[1]) + y00*(A[1] - origin_[1])*(B[1] - origin_[1]))*x00*z00;
                            double mzz = -over_pf*(z11 + z10*(B[2] - origin_[2]) + z01*(A[2] - origin_[2]) + z00*(A[2] - origin_[2])*(B[2] - origin_[2]))*x00*y00;
                            double mxy = -over_pf*(x01+x00*(B[0] - origin_[0]))*(y01+y00*(B[1] - origin_[1]))*z00;
                            double mxz = -over_pf*(x01+x00*(B[0] - origin_[0]))*y00*(z01+z00*(B[2] - origin_[2]));
                            double myz = -over_pf*x00*(y01+y00*(B[1] - origin_[1]))*(z01+z00*(B[2] - origin_[2]));

                            double mrr = (1.0 / 3.0) * (mxx + myy + mzz);

                            buffer_[ao12]        += (3.0/2.0) * (mxx - mrr);
                            buffer_[ao12+xydisp] += (3.0/2.0) * mxy;
                            buffer_[ao12+xzdisp] += (3.0/2.0) * mxz;
                            buffer_[ao12+yydisp] += (3.0/2.0) * (myy - mrr);
                            buffer_[ao12+yzdisp] += (3.0/2.0) * myz;
                            buffer_[ao12+zzdisp] += (3.0/2.0) * (mzz - mrr);

                            ao12++;
                        }
                    }
                }
            }
        }
    }
}
