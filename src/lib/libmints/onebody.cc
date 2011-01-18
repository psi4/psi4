#include <stdexcept>

#include "mints.h"
#include <exception.h>

using namespace psi;

void do_copy1(double *source, double *target, int chunk,
    int n1, int s1, int offset1,
    int n2, int s2, int offset2)
{
    int i1, i2;

    for (i1=0; i1<n1; i1++) {
        int off = ((offset1 + i1)*s2 + offset2)*chunk;
        for (i2=0; i2<n2*chunk; i2++, off++) {
            target[off] = source[off];
        }
    }
}

void do_sparse_transform11(double *source, double *target, int chunk,
    SphericalTransformIter& trans,
    int offsetcart1,
    int offsetpure1,
    int n2, int s2, int offset2)
{
    int i2;

    for (trans.first(); !trans.is_done(); trans.next()) {
        double coef = trans.coef();
        int pure = trans.pureindex();
        int cart = trans.cartindex();
        int offtarget = ((offsetpure1 + pure)*s2 + offset2)*chunk;
        int offsource = ((offsetcart1 + cart)*s2 + offset2)*chunk;
        for (i2=0; i2<n2*chunk; i2++) {
            target[offtarget++] += coef * source[offsource++];
        }
    }
}

void do_sparse_transform12(double *source, double *target, int chunk,
    SphericalTransformIter& trans,
    int n1, int offset1,
    int s2cart, int offsetcart2,
    int s2pure, int offsetpure2)
{
    int i1, ichunk;

    for (trans.first(); !trans.is_done(); trans.next()) {
        double coef = trans.coef();
        int pure = trans.pureindex();
        int cart = trans.cartindex();
        for (i1=0; i1<n1; i1++) {
            int offtarget = ((offset1 + i1)*s2pure + offsetpure2 + pure)*chunk;
            int offsource = ((offset1 + i1)*s2cart + offsetcart2 + cart)*chunk;
            for (ichunk=0; ichunk<chunk; ichunk++) {
                target[offtarget++] += coef * source[offsource++];
            }
        }
    }
}

OneBodyInt::OneBodyInt(std::vector<SphericalTransform> &spherical_transforms, shared_ptr<BasisSet> bs1, shared_ptr<BasisSet> bs2, int deriv)
    : bs1_(bs1), bs2_(bs2), spherical_transforms_(spherical_transforms), deriv_(deriv)
{
    buffer_ = 0;
    natom_ = bs1_->molecule()->natom();
}

OneBodyInt::~OneBodyInt()
{

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

void OneBodyInt::normalize_am(shared_ptr<GaussianShell> s1, shared_ptr<GaussianShell> s2, int nchunk)
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

void OneBodyInt::do_transform(shared_ptr<GaussianShell> s1, shared_ptr<GaussianShell> s2, int chunk)
{
    int i, j;
    int ogc1, ogc2;
    int ogc1pure, ogc2pure;
    int am1, am2;
    int pure1 = s1->is_pure();
    int pure2 = s2->is_pure();
    int ncart1 = s1->ncartesian();
    int ncart2 = s2->ncartesian();
    int nfunc1 = s1->nfunction();
    int nfunc2 = s2->nfunction();
    int nfunci, nfuncj;
    double *source;

    if (!pure1 && !pure2) {
        return;
    }

    if (pure1) {
        source = new double[ncart1*ncart2*chunk];
        memcpy(source, buffer_, sizeof(double)*ncart1*ncart2*chunk);
        memset(buffer_, 0, sizeof(double)*nfunc1*ncart2*chunk);

        ogc1 = 0;
        ogc1pure = 1;
        am1 = s1->am();
        nfunci = s1->nfunction();
        ogc2 = 0;
        am2 = s2->am();
        nfuncj = s2->nfunction();

        if (s1->is_pure()) {
            SphericalTransformIter trans(spherical_transforms_[s1->am()]);
            do_sparse_transform11(source, buffer_, chunk,
                                  trans,
                                  ogc1,
                                  ogc1pure,
                                  INT_NCART(am2), ncart2, ogc2);
        }
        else {
            do_copy1(source, buffer_, chunk,
                     nfunci, nfunc1, ogc1pure,
                     INT_NCART(am2), ncart2, ogc2);
        }
        delete[] source;
    }

    if (pure2) {
        source = new double[nfunc1*ncart2*chunk];
        memcpy(source, buffer_, nfunc1*ncart2*chunk);
        memset(buffer_, 0, sizeof(double)*s1->nfunction()*s2->nfunction()*chunk);

        ogc1 = 0;
        am1 = s1->am();
        nfunci = s1->nfunction();
        ogc2 = 0;
        ogc2pure = 0;
        am2 = s2->am();
        nfuncj = s2->nfunction();

        if (s2->is_pure()) {
            SphericalTransformIter trans(spherical_transforms_[s2->am()]);
            do_sparse_transform12(source, buffer_, chunk,
                                  trans,
                                  INT_NPURE(am1), ogc1,
                                  ncart2, ogc2,
                                  s2->nfunction(), ogc2pure);
        }
        else {
            do_copy1(source, buffer_, chunk,
                     nfunci, nfunc1, ogc1,
                     nfuncj, nfunc2, ogc2pure);
        }
        delete[] source;
    }
}

void OneBodyInt::spherical_transform(shared_ptr<GaussianShell> s1, shared_ptr<GaussianShell> s2)
{
    do_transform(s1, s2, 1);
}

// Converts the AO integrals stored in the buffer to SO integrals (plus spherical transform)
// This function does a full transform. It might be of use to create two half transform
// functions.
void OneBodyInt::so_transform(shared_ptr<Matrix> result, int sh1, int sh2, int ichunk)
{
#if 0
    // Get the transforms from the basis sets
    SOTransformIter trans1(bs1_->so_transform(sh1));
    SOTransformIter trans2(bs2_->so_transform(sh2));
    int nao2 = bs2_->shell(sh2)->ncartesian();
    size_t chunkoffset = ichunk * (bs1_->shell(sh1)->ncartesian() * bs2_->shell(sh2)->ncartesian());
    double *localbuffer = buffer_+chunkoffset;

    // Assume result is zeroed out where it matters.
    for (trans1.first(); !trans1.is_done(); trans1.next()) {
        int irrep1       = trans1.irrep();
        int aofunc1      = trans1.aofunc();
        int soirrepfunc1 = trans1.sofuncirrep();
        double coef1     = trans1.coef();
        int offset1      = aofunc1 * nao2;

        for (trans2.first(); !trans2.is_done(); trans2.next()) {
            int irrep2       = trans2.irrep();
            int aofunc2      = trans2.aofunc();
            int soirrepfunc2 = trans2.sofuncirrep();
            double coef2     = trans2.coef();

            // Okay, for one-electron integrals they must be the same irrep, unless they're special
            if (irrep1 == irrep2) {
                // Compute and store
                double val = coef1 * coef2 * localbuffer[offset1 + aofunc2];

                //fprintf(outfile,"Irrep: %d, Irrepfunc 1: %d, Irrepfunc 2: %d, coef1 = %lf, coef2 = %lf, localbuffer = %lf, Val %14.10f,\n",irrep1,soirrepfunc1,soirrepfunc2,coef1,coef2,localbuffer[offset1+aofunc2], val); fflush(outfile);
                result->add(irrep1, soirrepfunc1, soirrepfunc2, val);
            }
        }
    }
#endif
}

// Converts the AO integrals stored in the buffer to SO integrals (plus spherical transform)
// This function does a full transform. It might be of use to create two half transform
// functions.
void OneBodyInt::so_transform(shared_ptr<SimpleMatrix> result, int sh1, int sh2, int ichunk)
{
#if 0
    // Get the transforms from the basis sets
    SOTransformIter trans1(bs1_->so_transform(sh1));
    SOTransformIter trans2(bs2_->so_transform(sh2));
    int nao2 = bs2_->shell(sh2)->ncartesian();
    size_t chunkoffset = ichunk * (bs1_->shell(sh1)->ncartesian() * bs2_->shell(sh2)->ncartesian());
    double *localbuffer = buffer_+chunkoffset;

    // Assume result is zeroed out where it matters.
    for (trans1.first(); !trans1.is_done(); trans1.next()) {
        int aofunc1      = trans1.aofunc();
        int sofunc1      = trans1.sofunc();
        double coef1     = trans1.coef();
        int offset1      = aofunc1 * nao2;

        for (trans2.first(); !trans2.is_done(); trans2.next()) {
            int aofunc2      = trans2.aofunc();
            int sofunc2      = trans2.sofunc();
            double coef2     = trans2.coef();

            // Compute and store
            double val = coef1 * coef2 * localbuffer[offset1 + aofunc2];

            result->add(sofunc1, sofunc2, val);
        }
    }
#endif
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
            so_transform(result, i, j);
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
            for (int k=0; k<3*natom_; ++k)
                so_transform(result[k], i, j, k);
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
            for (int k=0; k<9*natom_; ++k)
                so_transform(result[k], i, j, k);
        }
    }
}

void OneBodyInt::compute_shell_deriv2(int, int)
{
    throw FeatureNotImplemented("libmints", "OneBodyInt::compute_shell_deriv2(Array)", __FILE__, __LINE__);
}

