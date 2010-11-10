#include <stdexcept>
#include <libqt/qt.h>
#include "mints.h"

using namespace psi;

ThreeCenterOverlapInt::ThreeCenterOverlapInt(std::vector<SphericalTransform>& st,
                                             boost::shared_ptr<BasisSet> bs1,
                                             boost::shared_ptr<BasisSet> bs2,
                                             boost::shared_ptr<BasisSet> bs3)
    : overlap_recur_(bs1->max_am(), bs2->max_am(), bs3->max_am()),
      bs1_(bs1), bs2_(bs2), bs3_(bs3)
{
    buffer_ = 0;
    size_t size;

    // Allocate memory for buffer_ storage
    try {
        size = INT_NCART(bs1->max_am()) * INT_NCART(bs2->max_am()) * INT_NCART(bs3->max_am());
        buffer_ = new double[size];
    }
    catch (std::bad_alloc& e) {
        fprintf(stderr, "Error allocating memory for buffer_\n");
        exit(EXIT_FAILURE);
    }
    memset(buffer_, 0, sizeof(double)*size);
}

ThreeCenterOverlapInt::~ThreeCenterOverlapInt()
{
    delete[] buffer_;
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

void ThreeCenterOverlapInt::compute_pair(shared_ptr<GaussianShell> sA, shared_ptr<GaussianShell> sC, shared_ptr<GaussianShell> sB)
{
    int ao132;
    int amA = sA->am();
    int amC = sC->am();
    int amB = sB->am();
    int nprimA = sA->nprimitive();
    int nprimC = sC->nprimitive();
    int nprimB = sB->nprimitive();
    double A[3], B[3], C[3], P[3], G[3], GA[3], GB[3], GC[3];
    A[0] = sA->center()[0];
    A[1] = sA->center()[1];
    A[2] = sA->center()[2];
    C[0] = sC->center()[0];
    C[1] = sC->center()[1];
    C[2] = sC->center()[2];
    B[0] = sB->center()[0];
    B[1] = sB->center()[1];
    B[2] = sB->center()[2];

    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, sA->ncartesian() * sC->ncartesian() * sB->ncartesian() * sizeof(double));

    double ***x = overlap_recur_.x();
    double ***y = overlap_recur_.y();
    double ***z = overlap_recur_.z();

    for (int pA=0; pA<nprimA; ++pA) {
        double aA = sA->exp(pA);
        double cA = sA->coef(pA);

        for (int pB=0; pB<nprimB; ++pB) {
            double aB = sB->exp(pB);
            double cB = sB->coef(pB);

            double gamma = aA + aB;
            double oog = 1.0 / gamma;

            P[0] = (aA * A[0] + aB * B[0]) * oog;

            double overlap_AB = exp(-aA*aB*AB2*oog) * sqrt(M_PI*oog) * M_PI * oog * cA * cB;

            for (int pC=0; pC<nprimC; ++pC) {
                double aC = sC->exp(pC);
                double cC = sC->coef(pC);

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

                overlap_recur_.compute(GA, GB, GC, gammac, amA, amC, amB);

                ao132 = 0;
                for(int ii = 0; ii <= amA; ii++) {
                    int lA = amA - ii;
                    for(int jj = 0; jj <= ii; jj++) {
                        int mA = ii - jj;
                        int nA = jj;
                        for(int kk = 0; kk <= amC; kk++) {
                            int lC = amC - kk;
                            for(int ll = 0; ll <= kk; ll++) {
                                int mC = kk - ll;
                                int nC = ll;

                                for(int mm = 0; mm <= amB; mm++) {
                                    int lB = amB - mm;
                                    for(int nn = 0; nn <= mm; nn++) {
                                        int mB = mm - nn;
                                        int nB = nn;

                                        double x0 = x[lA][lC][lB];
                                        double y0 = y[mA][mC][mB];
                                        double z0 = z[nA][nC][nB];

                                        buffer_[ao132++] += overlap_ACB*x0*y0*z0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    normalize_am(sA, sC, sB);
}

void ThreeCenterOverlapInt::normalize_am(shared_ptr<GaussianShell>& sA, shared_ptr<GaussianShell>& sC, shared_ptr<GaussianShell>& sB)
{
    // Integrals are done. Normalize for angular momentum
    int amA = sA->am();
    int amC = sC->am();
    int amB = sB->am();

    int length = INT_NCART(amA) * INT_NCART(amC) * INT_NCART(amB);

    int ao132 = 0;
    for(int ii = 0; ii <= amA; ii++) {
        int lA = amA - ii;
        for(int jj = 0; jj <= ii; jj++) {
            int mA = ii - jj;
            int nA = jj;
            for(int kk = 0; kk <= amC; kk++) {
                int lC = amC - kk;
                for(int ll = 0; ll <= kk; ll++) {
                    int mC = kk - ll;
                    int nC = ll;

                    for(int mm = 0; mm <= amB; mm++) {
                        int lB = amB - mm;
                        for(int nn = 0; nn <= mm; nn++) {
                            int mB = mm - nn;
                            int nB = nn;


                            buffer_[ao132] *= sA->normalize(lA, mA, nA) * sC->normalize(lC, mC, nC) * sB->normalize(lB, mB, nB);
                            ao132++;
                        }
                    }
                }
            }
        }
    }
}
