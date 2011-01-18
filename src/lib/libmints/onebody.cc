#include <stdexcept>

#include "mints.h"
#include <exception.h>

using namespace psi;

static void transform1e_1(int am, SphericalTransformIter& sti, double *s, double *t, int nl)
{
    memset(t,0,INT_NPURE(am)*nl*sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex()*nl;
        double *tptr = t + sti.pureindex()*nl;
        double coef = sti.coef();

//        fprintf(outfile, "1e_1: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(), sti.coef());

        for(int l=0; l<nl; l++) {
//            fprintf(outfile, "\ttptr = %8.5f coef = %8.5f sptr = %8.5f\n", *tptr, coef, *sptr);
            *(tptr++) += coef * *(sptr++);
        }
    }
}

static void transform1e_2(int am, SphericalTransformIter& sti, double *s, double *t, int nk, int nl)
{
    const int sl = nl;
    const int tl = INT_NPURE(am);

    memset(t,0,nk*tl*sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex();
        double *tptr = t + sti.pureindex();
        double coef = sti.coef();

//        fprintf(outfile, "1e_2: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(), sti.coef());

        for(int k=0; k<nk; k++,sptr+=sl,tptr+=tl) {
//            fprintf(outfile, "\ttptr = %8.5f coef = %8.5f sptr = %8.5f\n", *tptr, coef, *sptr);
            *(tptr) += coef * *(sptr);
        }
    }
}

OneBodyInt::OneBodyInt(std::vector<SphericalTransform> &spherical_transforms, shared_ptr<BasisSet> bs1, shared_ptr<BasisSet> bs2, int deriv)
    : bs1_(bs1), bs2_(bs2), spherical_transforms_(spherical_transforms), deriv_(deriv), nchunk_(1)
{
    buffer_ = 0;
    natom_ = bs1_->molecule()->natom();

    size_t buffsize = INT_NCART(bs1->max_am()) * INT_NCART(bs2->max_am());

    tformbuf_ = new double[buffsize];
    target_ = new double[buffsize];
}

OneBodyInt::~OneBodyInt()
{
    delete[] tformbuf_;
    delete[] target_;
}

shared_ptr<BasisSet> OneBodyInt::basis()
{
    return bs1_;
}

shared_ptr<BasisSet> OneBodyInt::basis1()
{
    return bs1_;
}

shared_ptr<BasisSet> OneBodyInt::basis2()
{
    return bs2_;
}

const double* OneBodyInt::buffer() const
{
    return buffer_;
}

bool OneBodyInt::cloneable()
{
    return false;
}

OneBodyInt* OneBodyInt::clone()
{
    throw FeatureNotImplemented("libmints", "OneBodyInt::clone()", __FILE__, __LINE__);
}

void OneBodyInt::normalize_am(const shared_ptr<GaussianShell>& s1, const shared_ptr<GaussianShell>& s2, int nchunk)
{
    // Integrals are done. Normalize for angular momentum
    int am1 = s1->am();
    int am2 = s2->am();
    int length = INT_NCART(am1) * INT_NCART(am2);

    int ao12 = 0;
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

                    for (int chunk=0; chunk<nchunk; ++chunk) {
                        buffer_[ao12+(chunk*length)] *= s1->normalize(l1, m1, n1) * s2->normalize(l2, m2, n2);
                    }
                    ao12++;
                }
            }
        }
    }
}

void OneBodyInt::pure_transform(const shared_ptr<GaussianShell>& s1,
                                const shared_ptr<GaussianShell>& s2, int chunk)
{
    double *source1, *target1;
    double *source2, *target2;
    double *source = buffer_;
    double *target = target_;
    double *tmpbuf = tformbuf_;

    int am1 = s1->am();
    int is_pure1 = s1->is_pure() && am1 > 1;
    int ncart1 = s1->ncartesian();
    int nbf1 = s1->nfunction();

    int am2 = s2->am();
    int is_pure2 = s2->is_pure() && am2 > 1;
    int ncart2 = s2->ncartesian();
    int nbf2 = s2->nfunction();

    int transform_index = 2*is_pure1 + is_pure2;
    switch(transform_index) {
    case 0:
        break;
    case 1:
        source2 = source;
        target2 = target;
        break;
    case 2:
        source1 = source;
        target1 = target;
        break;
    case 3:
        source2 = source;
        target2 = tmpbuf;
        source1 = tmpbuf;
        target1 = target;
        break;
    }

    if (is_pure2) {
        SphericalTransformIter stiter(spherical_transforms_[am2]);
        transform1e_2(am2, stiter, source2, target2, ncart1, ncart2);
    }
    if (is_pure1) {
        SphericalTransformIter stiter(spherical_transforms_[am1]);
        transform1e_1(am1, stiter, source1, target1, nbf2);
    }

    if (transform_index) {
        memcpy(buffer_, target_, sizeof(double) * ncart1 * ncart2);
    }
}

void OneBodyInt::compute_shell(int sh1, int sh2)
{
    const shared_ptr<GaussianShell> s1 = bs1_->shell(sh1);
    const shared_ptr<GaussianShell> s2 = bs2_->shell(sh2);

    // Call the child's compute_pair method, results better be in buffer_.
    compute_pair(s1, s2);

    // Normalize for angular momentum
    normalize_am(s1, s2, nchunk_);

    // Pure angular momentum (6d->5d, ...) transformation
    pure_transform(s1, s2, nchunk_);
}

void OneBodyInt::compute(shared_ptr<Matrix> result)
{
    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();

    for (int i=0; i<ns1; ++i) {
        for (int j=0; j<ns2; ++j) {
            // Compute the shell
            compute_shell(i, j);
            // Transform the shell to SO basis
//            pure_transform(result, i, j);
        }
    }
}

void OneBodyInt::compute(std::vector<shared_ptr<Matrix> > &result)
{
    throw FeatureNotImplemented("libmints", "OneBodyInt::compute(Array)", __FILE__, __LINE__);
}

void OneBodyInt::compute(std::vector<shared_ptr<SimpleMatrix> > &result)
{
    throw FeatureNotImplemented("libmints", "OneBodyInt::compute(SimpleArray)", __FILE__, __LINE__);
}

void OneBodyInt::compute_deriv1(std::vector<shared_ptr<Matrix> > &result)
{
    throw FeatureNotImplemented("libmints", "OneBodyInt::compute_deriv1(Array)", __FILE__, __LINE__);
}

void OneBodyInt::compute_deriv1(std::vector<shared_ptr<SimpleMatrix> > &result)
{
    if (deriv_ < 1)
        throw SanityCheckError("OneBodyInt::compute_deriv1(result): integral object not created to handle derivatives.", __FILE__, __LINE__);

    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();

    // Check the length of result, must be 3*natom_
    if (result.size() != 3*natom_)
        throw SanityCheckError("OneBodyInt::compute_derv1(result): result must be 3 * natom in length.", __FILE__, __LINE__);

    for (int i=0; i<ns1; ++i) {
        for (int j=0; j<ns2; ++j) {
            // Compute the shell
            compute_shell_deriv1(i, j);
            // Transform the shell to SO basis
//            for (int k=0; k<3*natom_; ++k)
//                so_transform(result[k], i, j, k);
        }
    }
}

void OneBodyInt::compute_shell_deriv1(int, int)
{
    throw FeatureNotImplemented("libmints", "OneBodyInt::compute_shell_deriv1(Array)", __FILE__, __LINE__);
}

void OneBodyInt::compute_deriv2(std::vector<shared_ptr<SimpleMatrix> > &result)
{
    if (deriv_ < 2)
        throw SanityCheckError("OneBodyInt::compute_deriv1(result): integral object not created to handle derivatives.", __FILE__, __LINE__);

    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();

    // Check the length of result, must be 9*natom_
    if (result.size() != 9*natom_)
        throw SanityCheckError("OneBodyInt::compute_derv2(result): result must be 9 * natom in length.", __FILE__, __LINE__);

    for (int i=0; i<ns1; ++i) {
        for (int j=0; j<ns2; ++i) {
            // Compute the shell
            compute_shell_deriv1(i, j);
            // Transform the shell to SO basis
//            for (int k=0; k<9*natom_; ++k)
//                so_transform(result[k], i, j, k);
        }
    }
}

void OneBodyInt::compute_shell_deriv2(int, int)
{
    throw FeatureNotImplemented("libmints", "OneBodyInt::compute_shell_deriv2(Array)", __FILE__, __LINE__);
}

