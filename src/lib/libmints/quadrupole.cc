#include <libciomr/libciomr.h>

#include <libmints/basisset.h>
#include <libmints/gshell.h>
#include <libmints/overlap.h>
#include <libmints/quadrupole.h>
#include <physconst.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace psi;

// Initialize overlap_recur_ to +2 basis set angular momentum
QuadrupoleInt::QuadrupoleInt(IntegralFactory* integral, BasisSet* bs1, BasisSet* bs2) :
    OneBodyInt(integral, bs1, bs2), overlap_recur_(bs1->max_am()+2, bs2->max_am()+2)
{
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();
    
    int maxnao1 = (maxam1+1)*(maxam1+2)/2;
    int maxnao2 = (maxam2+1)*(maxam2+2)/2;
    
    // Increase buffer size to handle xx, xy, xz, yy, yz, zz components
    buffer_ = new double[6*maxnao1*maxnao2];
}

QuadrupoleInt::~QuadrupoleInt()
{
    delete[] buffer_;
}

void QuadrupoleInt::compute_shell(int sh1, int sh2)
{
    compute_pair(bs1_->shell(sh1), bs2_->shell(sh2));
}

// The engine only supports segmented basis sets
void QuadrupoleInt::compute_pair(GaussianShell* s1, GaussianShell* s2)
{
    int ao12;
    int am1 = s1->am(0);
    int am2 = s2->am(0);
    int nprim1 = s1->nprimitive();
    int nprim2 = s2->nprimitive();
    double A[3], B[3];
    A[0] = s1->center()[0];
    A[1] = s1->center()[1];
    A[2] = s1->center()[2];
    B[0] = s2->center()[0];
    B[1] = s2->center()[1];
    B[2] = s2->center()[2];
    
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
    
    memset(buffer_, 0, 6 * s1->ncartesian() * s2->ncartesian() * sizeof(double));
    
    double **x = overlap_recur_.x();
    double **y = overlap_recur_.y();
    double **z = overlap_recur_.z();
    
    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1->exp(p1);
        double c1 = s1->coef(0, p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2->exp(p2);
            double c2 = s2->coef(0, p2);
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
                            
                            double mxx = -over_pf*(x11 + x10*(B[0]) + x01*(A[0]) + x00*(A[0])*(B[0]))*y00*z00;
                            double myy = -over_pf*(y11 + y10*(B[1]) + y01*(A[1]) + y00*(A[1])*(B[1]))*x00*z00;
                            double mzz = -over_pf*(z11 + z10*(B[2]) + z01*(A[2]) + z00*(A[2])*(B[2]))*x00*y00;
                            double mxy = -over_pf*(x01+x00*(B[0]))*(y01+y00*(B[1]))*z00;
                            double mxz = -over_pf*(x01+x00*(B[0]))*y00*(z01+z00*(B[2]));
                            double myz = -over_pf*x00*(y01+y00*(B[1]))*(z01+z00*(B[2]));
                            
                            buffer_[ao12]        += mxx;
                            buffer_[ao12+xydisp] += mxy;
                            buffer_[ao12+xzdisp] += mxz;
                            buffer_[ao12+yydisp] += myy;
                            buffer_[ao12+yzdisp] += myz;
                            buffer_[ao12+zzdisp] += mzz;
                            
                            ao12++;
                        }
                    }
                }
            }
        }
    }
    
    // Integrals are done. Normalize for angular momentum
    normalize_am(s1, s2, 6);
}

void QuadrupoleInt::spherical_transform(GaussianShell* s1, GaussianShell* s2)
{
    do_transform(s1, s2, 6);
}

void QuadrupoleInt::compute(Matrix** result)
{
    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();
    
    for (int i=0; i<ns1; ++i) {
        for (int j=0; j<ns2; ++j) {
            // Compute the shell
            compute_shell(i, j);
            // Transform the shell to SO basis
            so_transform(result[0], i, j, 0);
            so_transform(result[1], i, j, 1);
            so_transform(result[2], i, j, 2);
            so_transform(result[3], i, j, 3);
            so_transform(result[4], i, j, 4);
            so_transform(result[5], i, j, 5);
        }
    }
}

void QuadrupoleInt::compute(SimpleMatrix** result)
{
    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();
    
    for (int i=0; i<ns1; ++i) {
        for (int j=0; j<ns2; ++j) {
            // Compute the shell
            compute_shell(i, j);
            // Transform the shell to SO basis
            so_transform(result[0], i, j, 0);
            so_transform(result[1], i, j, 1);
            so_transform(result[2], i, j, 2);
            so_transform(result[3], i, j, 3);
            so_transform(result[4], i, j, 4);
            so_transform(result[5], i, j, 5);
        }
    }
}
