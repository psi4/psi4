#include <boost/shared_ptr.hpp>
#include <stdexcept>
#include <libqt/qt.h>
#include "mints.h"

using namespace psi;

ThreeCenterOverlapInt::ThreeCenterOverlapInt(std::vector<SphericalTransform>& st,
                                             boost::shared_ptr<BasisSet> bs1,
                                             boost::shared_ptr<BasisSet> bs2,
                                             boost::shared_ptr<BasisSet> bs3)
    : overlap_recur_(bs1->max_am(), bs2->max_am(), bs3->max_am()),
      bs1_(bs1), bs2_(bs2), bs3_(bs3), st_(st)
{
    size_t size = INT_NCART(bs1->max_am()) * INT_NCART(bs2->max_am()) * INT_NCART(bs3->max_am());

    // Allocate memory for buffer_ storage
    try {
        buffer_ = new double[size];
    }
    catch (std::bad_alloc& e) {
        fprintf(stderr, "Error allocating memory for buffer_\n");
        exit(EXIT_FAILURE);
    }
    ::memset(buffer_, 0, sizeof(double)*size);

    try {
        temp_ = new double[size];
    }
    catch (std::bad_alloc& e) {
        fprintf(stderr, "Error allocating memory for temp_\n");
        exit(EXIT_FAILURE);
    }
    ::memset(temp_, 0, sizeof(double)*size);
}

ThreeCenterOverlapInt::~ThreeCenterOverlapInt()
{
    delete[] buffer_;
    delete[] temp_;
}

boost::shared_ptr<BasisSet> ThreeCenterOverlapInt::basis()
{
    return bs1_;
}

boost::shared_ptr<BasisSet> ThreeCenterOverlapInt::basis1()
{
    return bs1_;
}

boost::shared_ptr<BasisSet> ThreeCenterOverlapInt::basis2()
{
    return bs2_;
}

boost::shared_ptr<BasisSet> ThreeCenterOverlapInt::basis3()
{
    return bs3_;
}

void ThreeCenterOverlapInt::compute_shell(int sh1, int sh2, int sh3)
{
    compute_pair(bs1_->shell(sh1), bs2_->shell(sh2), bs3_->shell(sh3));
}

void ThreeCenterOverlapInt::compute_pair(const GaussianShell& sA,
                                         const GaussianShell& sB,
                                         const GaussianShell& sC)
{
    unsigned int ao123;
    int amA = sA.am();
    int amB = sB.am();
    int amC = sC.am();
    int nprimA = sA.nprimitive();
    int nprimB = sB.nprimitive();
    int nprimC = sC.nprimitive();
    double A[3], B[3], C[3], P[3], G[3], GA[3], GB[3], GC[3];
    A[0] = sA.center()[0];
    A[1] = sA.center()[1];
    A[2] = sA.center()[2];
    B[0] = sB.center()[0];
    B[1] = sB.center()[1];
    B[2] = sB.center()[2];
    C[0] = sC.center()[0];
    C[1] = sC.center()[1];
    C[2] = sC.center()[2];

    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, sA.ncartesian() * sC.ncartesian() * sB.ncartesian() * sizeof(double));

    double ***x = overlap_recur_.x();
    double ***y = overlap_recur_.y();
    double ***z = overlap_recur_.z();

    for (int pA=0; pA<nprimA; ++pA) {
        double aA = sA.exp(pA);
        double cA = sA.coef(pA);

        for (int pB=0; pB<nprimB; ++pB) {
            double aB = sB.exp(pB);
            double cB = sB.coef(pB);

            double gamma = aA + aB;
            double oog = 1.0 / gamma;

            P[0] = (aA * A[0] + aB * B[0]) * oog;
            P[1] = (aA * A[1] + aB * B[1]) * oog;
            P[2] = (aA * A[2] + aB * B[2]) * oog;

            double overlap_AB = exp(-aA*aB*AB2*oog) * sqrt(M_PI*oog) * M_PI * oog * cA * cB;

            for (int pC=0; pC<nprimC; ++pC) {
                double aC = sC.exp(pC);
                double cC = sC.coef(pC);

                double PC2 = 0.0;
                PC2 += (P[0] - C[0]) * (P[0] - C[0]);
                PC2 += (P[1] - C[1]) * (P[1] - C[1]);
                PC2 += (P[2] - C[2]) * (P[2] - C[2]);

                double gammac = gamma + aC;
                double oogc = 1.0 / (gammac);

                G[0] = (gamma * P[0] + aC * C[0]) * oogc;
                G[1] = (gamma * P[1] + aC * C[1]) * oogc;
                G[2] = (gamma * P[2] + aC * C[2]) * oogc;

                GA[0] = G[0] - A[0];
                GA[1] = G[1] - A[1];
                GA[2] = G[2] - A[2];
                GB[0] = G[0] - B[0];
                GB[1] = G[1] - B[1];
                GB[2] = G[2] - B[2];
                GC[0] = G[0] - C[0];
                GC[1] = G[1] - C[1];
                GC[2] = G[2] - C[2];

                double overlap_ACB = exp(-gamma*aC*oogc*PC2) * sqrt(gamma * oogc) * (gamma *oogc) * overlap_AB * cC;

                 // Computes (ACB) overlap
                overlap_recur_.compute(GA, GB, GC, gammac, amA, amB, amC);

                // We're going to be reordering the result of the above line.
                // The result of above B is the fast running index, but I'm going to be make it C instead.
                ao123 = 0;
                for(int ii = 0; ii <= amA; ii++) {
                    int lA = amA - ii;
                    for(int jj = 0; jj <= ii; jj++) {
                        int mA = ii - jj;
                        int nA = jj;

                        for(int mm = 0; mm <= amB; mm++) {
                            int lB = amB - mm;
                            for(int nn = 0; nn <= mm; nn++) {
                                int mB = mm - nn;
                                int nB = nn;

                                for(int kk = 0; kk <= amC; kk++) {
                                    int lC = amC - kk;
                                    for(int ll = 0; ll <= kk; ll++) {
                                        int mC = kk - ll;
                                        int nC = ll;

                                        // These are ordered (ACB) -> B fast running
                                        double x0 = x[lA][lC][lB];
                                        double y0 = y[mA][mC][mB];
                                        double z0 = z[nA][nC][nB];

                                        // But we're going to store then like (ABC) -> C fast running
                                        buffer_[ao123++] += overlap_ACB*x0*y0*z0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    normalize_am(sA, sB, sC);
    pure_transform(sA, sB, sC);
}

void ThreeCenterOverlapInt::normalize_am(const GaussianShell& sA,
                                         const GaussianShell& sB,
                                         const GaussianShell& sC)
{
    // Assume integrals are done. Normalize for angular momentum
    int amA = sA.am();
    int amB = sB.am();
    int amC = sC.am();

    size_t ao123 = 0;
    for(int ii = 0; ii <= amA; ii++) {
        int lA = amA - ii;
        for(int jj = 0; jj <= ii; jj++) {
            int mA = ii - jj;
            int nA = jj;

            double normA = GaussianShell::normalize(lA, mA, nA);

            for(int mm = 0; mm <= amB; mm++) {
                int lB = amB - mm;
                for(int nn = 0; nn <= mm; nn++) {
                    int mB = mm - nn;
                    int nB = nn;

                    double normB = GaussianShell::normalize(lB, mB, nB);

                    for(int kk = 0; kk <= amC; kk++) {
                        int lC = amC - kk;
                        for(int ll = 0; ll <= kk; ll++) {
                            int mC = kk - ll;
                            int nC = ll;

                            buffer_[ao123] *= normA * normB * GaussianShell::normalize(lC, mC, nC);
                            ao123++;
                        }
                    }
                }
            }
        }
    }
}

void ThreeCenterOverlapInt::pure_transform(const GaussianShell& s1,
                                           const GaussianShell& s2,
                                           const GaussianShell& s3)
{
    // Get the transforms from the basis set
    SphericalTransformIter trans1(st_[s1.am()]);
    SphericalTransformIter trans2(st_[s2.am()]);
    SphericalTransformIter trans3(st_[s3.am()]);

    // Get the angular momentum for each shell
    int am1 = s1.am();
    int am2 = s2.am();
    int am3 = s3.am();

    // Get number of Cartesian functions for each shell
    int nao1 = s1.ncartesian();
    int nao2 = s2.ncartesian();
    int nao3 = s3.ncartesian();

    // Get number of Basis functions for each shell
    int nso1 = s1.nfunction();
    int nso2 = s2.nfunction();
    int nso3 = s3.nfunction();

    // Get if each shell has pure functions
    bool is_pure1 = s1.is_pure();
    bool is_pure2 = s2.is_pure();
    bool is_pure3 = s3.is_pure();

    // ABC -> ABc
    if (is_pure3) {
        
        ::memset(temp_, '\0', sizeof(double) * nao1 * nao2 * nso3);

        for (trans3.first(); !trans3.is_done(); trans3.next()) {
            double *sptr = buffer_ + trans3.cartindex();
            double *tptr = temp_   + trans3.pureindex();
            double coef = trans3.coef();
            C_DAXPY(nao1 * nao2, coef, sptr, nao3, tptr, nso3);
        }

        ::memcpy((void*) buffer_, (void*) temp_, sizeof(double) * nao1 * nao2 * nso3); 
    }

    // ABc -> Abc
    if (is_pure2) {
        
        ::memset(temp_, '\0', sizeof(double) * nao1 * nso2 * nso3);

        for (trans2.first(); !trans2.is_done(); trans2.next()) {
            double coef = trans2.coef();
            double *sptr = buffer_ + trans2.cartindex() * nso3;
            double *tptr = temp_   + trans2.pureindex() * nso3;
            for (int a=0; a<nao1; ++a) {
                C_DAXPY(nso3,coef,sptr,1,tptr,1); 
                sptr += nao2 * nso3; 
                tptr += nso2 * nso3; 
            } 
        }

        ::memcpy((void*) buffer_, (void*) temp_, sizeof(double) * nao1 * nso2 * nso3); 
    }

    // Abc -> abc
    if (is_pure1) {
        
        ::memset(temp_, '\0', sizeof(double) * nso1 * nso2 * nso3);

        for (trans1.first(); !trans1.is_done(); trans1.next()) {
            double *sptr = buffer_ + trans1.cartindex() * nso2 * nso3;
            double *tptr = temp_   + trans1.pureindex() * nso2 * nso3;
            double coef = trans1.coef();
            C_DAXPY(nso2 * nso3, coef, sptr, 1, tptr, 1);
        }

        ::memcpy((void*) buffer_, (void*) temp_, sizeof(double) * nso1 * nso2 * nso3); 
    }
}
