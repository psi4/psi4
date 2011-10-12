#include <boost/shared_ptr.hpp>
#include <stdexcept>
#include <libqt/qt.h>
#include "mints.h"

using namespace psi;

static void transform3c_1(int, SphericalTransformIter&, double*, double*, int);
static void transform3c_2(int, SphericalTransformIter&, double*, double*, int, int, int);
static void transform3c_3(int, SphericalTransformIter&, double*, double*, int, int);

static void transform3c_1(int am, SphericalTransformIter& sti, double *s, double *t, int ncb)
{
    memset(t, 0, INT_NPURE(am)*ncb*sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex() * ncb;
        double *tptr = t + sti.pureindex() * ncb;
        double coef = sti.coef();
        for (int cb=0; cb<ncb; ++cb)
            *(tptr++) += coef * *(sptr++);
    }
}

static void transform3c_2(int am, SphericalTransformIter &sti, double *s, double *t, int na, int nc, int nb)
{
    int sc = INT_NPURE(am);
    const int scb = nc * nb;
    const int tcb = sc * nb;

    memset(t, 0, na*tcb*sizeof(double));

    int interval = 0;
    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex()*nb;
        double *tptr = t + sti.pureindex()*nb;
        double coef = sti.coef();
        for (int a=0; a<na; ++a,sptr+=scb, tptr+=tcb) {
            for (int b=0; b<nb; ++b) {
                tptr[b] += coef * sptr[b];
            }
        }
        interval++;
    }
}

static void transform3c_3(int am, SphericalTransformIter &sti, double *s, double *t, int nac, int nb)
{
    const int sb = nb;
    const int tb = INT_NPURE(am);

    // Clear out target memory
    memset(t, 0, nac*tb*sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex();
        double *tptr = t + sti.pureindex();
        double coef = sti.coef();

        for (int ac=0; ac<nac; ++ac) {
            *(tptr) += coef * *(sptr);

            // skip ahead to the next ijk
            sptr += sb;
            tptr += tb;
        }
    }
}

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
    memset(buffer_, 0, sizeof(double)*size);

    try {
        target_ = new double[size];
    }
    catch (std::bad_alloc& e) {
        fprintf(stderr, "Error allocating memory for target_\n");
        exit(EXIT_FAILURE);
    }
    memset(target_, 0, sizeof(double)*size);

    try {
        tformbuf_ = new double[size];
    }
    catch (std::bad_alloc& e) {
        fprintf(stderr, "Error allocating memory for tformbuf_");
        exit(EXIT_FAILURE);
    }
    memset(tformbuf_, 0, sizeof(double)*size);
}

ThreeCenterOverlapInt::~ThreeCenterOverlapInt()
{
    delete[] buffer_;
    delete[] target_;
    delete[] tformbuf_;
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

void ThreeCenterOverlapInt::compute_pair(boost::shared_ptr<GaussianShell> sA,
                                         boost::shared_ptr<GaussianShell> sB,
                                         boost::shared_ptr<GaussianShell> sC)
{
    unsigned int ao123;
    int amA = sA->am();
    int amB = sB->am();
    int amC = sC->am();
    int nprimA = sA->nprimitive();
    int nprimB = sB->nprimitive();
    int nprimC = sC->nprimitive();
    double A[3], B[3], C[3], P[3], G[3], GA[3], GB[3], GC[3];
    A[0] = sA->center()[0];
    A[1] = sA->center()[1];
    A[2] = sA->center()[2];
    B[0] = sB->center()[0];
    B[1] = sB->center()[1];
    B[2] = sB->center()[2];
    C[0] = sC->center()[0];
    C[1] = sC->center()[1];
    C[2] = sC->center()[2];

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
            P[1] = (aA * A[1] + aB * B[1]) * oog;
            P[2] = (aA * A[2] + aB * B[2]) * oog;

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

void ThreeCenterOverlapInt::normalize_am(boost::shared_ptr<GaussianShell>& sA,
                                         boost::shared_ptr<GaussianShell>& sB,
                                         boost::shared_ptr<GaussianShell>& sC)
{
    // Assume integrals are done. Normalize for angular momentum
    int amA = sA->am();
    int amB = sB->am();
    int amC = sC->am();

    int length = INT_NCART(amA) * INT_NCART(amC) * INT_NCART(amB);

    int ao123 = 0;
    for(int ii = 0; ii <= amA; ii++) {
        int lA = amA - ii;
        for(int jj = 0; jj <= ii; jj++) {
            int mA = ii - jj;
            int nA = jj;

            double normA = sA->normalize(lA, mA, nA);

            for(int mm = 0; mm <= amB; mm++) {
                int lB = amB - mm;
                for(int nn = 0; nn <= mm; nn++) {
                    int mB = mm - nn;
                    int nB = nn;

                    double normB = sB->normalize(lB, mB, nB);

                    for(int kk = 0; kk <= amC; kk++) {
                        int lC = amC - kk;
                        for(int ll = 0; ll <= kk; ll++) {
                            int mC = kk - ll;
                            int nC = ll;

                            buffer_[ao123] *= normA * normB * sC->normalize(lC, mC, nC);
                            ao123++;
                        }
                    }
                }
            }
        }
    }
}

void ThreeCenterOverlapInt::pure_transform(boost::shared_ptr<GaussianShell>& s1,
                                           boost::shared_ptr<GaussianShell>& s2,
                                           boost::shared_ptr<GaussianShell>& s3)
{
    // Get the transforms from the basis set
    SphericalTransformIter trans1(st_[s1->am()]);
    SphericalTransformIter trans2(st_[s2->am()]);
    SphericalTransformIter trans3(st_[s3->am()]);

    // Get the angular momentum for each shell
    int am1 = s1->am();
    int am2 = s2->am();
    int am3 = s3->am();

    // Get number of Cartesian functions for each shell
    int nao1 = s1->ncartesian();
    int nao2 = s2->ncartesian();
    int nao3 = s3->ncartesian();

    int nbf1 = s1->nfunction();
    int nbf2 = s2->nfunction();
    int nbf3 = s3->nfunction();

    // Get if each shell has pure functions
    bool is_pure1 = s1->is_pure();
    bool is_pure2 = s2->is_pure();
    bool is_pure3 = s3->is_pure();

    double *source1, *target1;
    double *source2, *target2;
    double *source3, *target3;

    double *source = buffer_;
    double *target = target_;
    double *tmpbuf = tformbuf_;

    int transform_index = 4*is_pure1 + 2*is_pure2 + is_pure3;
    switch (transform_index) {
    case 0:  // no transform
        break;

    case 1:  // (a|b|cT)
        source3 = source;
        target3 = target;
        break;

    case 2:  // (a|bT|c)
        source2 = source;
        target2 = target;
        break;

    case 3:  // (a|bT|cT)
        source3 = source;
        target3 = tmpbuf;
        source2 = tmpbuf;
        target2 = target;
        break;

    case 4:  // (aT|b|c)
        source1 = source;
        target1 = target;
        break;

    case 5:  // (aT|b|cT)
        source3 = source;
        target3 = tmpbuf;
        source1 = tmpbuf;
        target1 = target;
        break;

    case 6:  // (aT|bT|c)
        source2 = source;
        target2 = tmpbuf;
        source1 = tmpbuf;
        target1 = target;
        break;

    case 7: // (aT|bT|cT)
        source3 = source;
        target3 = tmpbuf;
        source2 = tmpbuf;
        target2 = source;
        source1 = source;
        target1 = target;
        break;
    }

    size_t size = 1;
    if (is_pure3) {
        transform3c_3(am3, trans3, source3, target3, nao1*nao2, nao3);
        size *= nbf3;
    }
    if (is_pure2) {
        transform3c_2(am2, trans2, source2, target2, nao1, nao2, nbf3);
        size *= nbf2;
    }
    if (is_pure1) {
        transform3c_1(am1, trans1, source1, target1, nbf2*nbf3);
        size *= nbf1;
    }
}
