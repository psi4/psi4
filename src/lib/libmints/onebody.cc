#include <stdexcept>

#include "mints.h"
#include <exception.h>

using namespace psi;

//void do_copy1(double *source, double *target, int chunk,
//    int n1, int s1, int offset1,
//    int n2, int s2, int offset2)
//{
//    int i1, i2;

//    for (i1=0; i1<n1; i1++) {
//        int off = ((offset1 + i1)*s2 + offset2)*chunk;
//        for (i2=0; i2<n2*chunk; i2++, off++) {
//            target[off] = source[off];
//        }
//    }
//}

//void do_sparse_transform11(double *source, double *target, int chunk,
//    SphericalTransformIter& trans,
//    int offsetcart1,
//    int offsetpure1,
//    int n2, int s2, int offset2)
//{
//    int i2;

//    fprintf(outfile, "offsetcart1 = %d offsetpure1 = %d n2 = %d s2 = %d offset2 = %d\n",
//            offsetcart1, offsetpure1, n2, s2, offset2);

//    for (trans.first(); !trans.is_done(); trans.next()) {
//        double coef = trans.coef();
//        int pure = trans.pureindex();
//        int cart = trans.cartindex();
//        int offtarget = ((offsetpure1 + pure)*s2 + offset2)*chunk;
//        int offsource = ((offsetcart1 + cart)*s2 + offset2)*chunk;
//        for (i2=0; i2<n2*chunk; i2++) {
//            target[offtarget++] += coef * source[offsource++];
//        }
//    }
//}

//void do_sparse_transform12(double *source, double *target, int chunk,
//    SphericalTransformIter& trans,
//    int n1, int offset1,
//    int s2cart, int offsetcart2,
//    int s2pure, int offsetpure2)
//{
//    int i1, ichunk;

//    fprintf(outfile, "n1 = %d offset1 = %d s2cart = %d offsetcart2 = %d s2pure = %d offsetpure = %d\n",
//            n1, offset1, s2cart, offsetcart2, s2pure, offsetpure2);

//    for (trans.first(); !trans.is_done(); trans.next()) {
//        double coef = trans.coef();
//        int pure = trans.pureindex();
//        int cart = trans.cartindex();
//        for (i1=0; i1<n1; i1++) {
//            int offtarget = ((offset1 + i1)*s2pure + offsetpure2 + pure)*chunk;
//            int offsource = ((offset1 + i1)*s2cart + offsetcart2 + cart)*chunk;
//            for (ichunk=0; ichunk<chunk; ichunk++) {
//                target[offtarget++] += coef * source[offsource++];
//            }
//        }
//    }
//}

static void transform1e_1(int am, SphericalTransformIter& sti, double *s, double *t, int nl)
{
    memset(t,0,INT_NPURE(am)*nl*sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex()*nl;
        double *tptr = t + sti.pureindex()*nl;
        double coef = sti.coef();

        fprintf(outfile, "1e_1: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(), sti.coef());

        for(int l=0; l<nl; l++) {
            fprintf(outfile, "\ttptr = %8.5f coef = %8.5f sptr = %8.5f\n", *tptr, coef, *sptr);
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

        fprintf(outfile, "1e_2: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(), sti.coef());

        for(int k=0; k<nk; k++,sptr+=sl,tptr+=tl) {
            fprintf(outfile, "\ttptr = %8.5f coef = %8.5f sptr = %8.5f\n", *tptr, coef, *sptr);
            *(tptr) += coef * *(sptr);
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

//    for (int z=0; z < ncart1*ncart2; ++z) {
//        fprintf(outfile, "tform target: %d -> %8.5f\n", z, target_[z]);
//        fprintf(outfile, "tform buffer: %d -> %8.5f\n", z, buffer_[z]);
//    }
    if (transform_index) {
        memcpy(buffer_, target_, sizeof(double) * ncart1 * ncart2);
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

