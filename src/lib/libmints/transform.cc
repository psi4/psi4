#include "wavefunction.h"
#include "matrix.h"
#include "integral.h"

#include <psiconfig.h>
#include <psi4-dec.h>

#include <cmath>

using namespace psi;

extern void solidharmonic(int l, SimpleMatrix& coefmat);

// there ordering here is arbitrary and doesn't have to match the
// basis set ordering
static inline int ncart(int l) { return (l>=0)?((((l)+2)*((l)+1))>>1):0; }
static inline int npure(int l) { return 2*l+1; }
static inline int icart(int a, int b, int c)
{
  return (((((a+b+c+1)<<1)-a)*(a+1))>>1)-b-1;
}
static inline int ipure(int l, int m) { return m<0?2*-m:(m==0?0:2*m-1); }

void SphericalTransformComponent::init(int a, int b, int c, double coef,
                                       int cartindex, int pureindex)
{
    fprintf(outfile, "a = %d, b = %d, c = %d, coef = %f, cartindex = %d, pureindex = %d\n", a, b, c, coef, cartindex, pureindex);
    a_ = a;
    b_ = b;
    c_ = c;
    coef_ = coef;
    cartindex_ = cartindex;
    pureindex_ = pureindex;
}

SphericalTransform::SphericalTransform()
{

}

SphericalTransform::SphericalTransform(int l, int subl) : l_(l)
{
    if (subl == -1) subl_ = l;
    else subl_ = subl;

    init();
}

void SphericalTransform::init()
{
    int cartdim = INT_NCART(l_);
    SimpleMatrix coefmat(cartdim, cartdim);
    coefmat.zero();

    // Compute the solid harmonic matrix elements
    solidharmonic(l_, coefmat);

    coefmat.set_name("SphericalTransform: COEFMAT");
    coefmat.print();
    fflush(outfile);

    // Go through and grab the values.
    int pureindex = 0;
    int cartindex = 0;

    for (int i=1; i<=(l_-subl_)/2; ++i)
        pureindex += npure(subl_+2*i);

    fprintf(outfile, "in SphericalTransform::init\n");
    fprintf(outfile, "l_ = %d, subl_ = %d\n", l_, subl_);

    for (int p=0; p<npure(subl_); ++p) {
        cartindex = 0;
//        for (int ii=0; ii<=l_; ++ii) {
//            int a = l_ - ii;
//            for (int jj=0; jj<=ii; ++jj) {
//                int b = ii - jj;
//                int c = jj;

        for (int a=0; a<=l_; ++a) {
            for (int b=0; (a+b)<=l_; ++b) {
                int c = l_ - a -b;

                fprintf(outfile, "a = %d b = %d c = %d\n", a, b, c);
                int cart = icart(a, b, c);

                double coef = coefmat(cart, p+pureindex);

                if (fabs(coef) > 1.0e-16) {
                    SphericalTransformComponent component;
                    component.init(a, b, c, coef, cartindex, p);
                    components_.push_back(component);
                }
                cartindex++;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

ISphericalTransform::ISphericalTransform()
    : SphericalTransform()
{

}

ISphericalTransform::ISphericalTransform(int l, int subl)
    : SphericalTransform(l, subl)
{
    printf("in ISphericalTranform constructor\n");
    components_.clear();
    init();
}

void ISphericalTransform::init()
{
    printf("in ISphericalTransform.init\n");
    int cartdim = ncart(l_);
    SimpleMatrix coefmat(cartdim, cartdim);
//    SimpleMatrix coefmatI(cartdim, cartdim);
//    SimpleMatrix temp(cartdim, cartdim);
    coefmat.zero();

    // Compute the solid harmonic matrix elements
    solidharmonic(l_, coefmat);
//    solidharmonic(l_, coefmatI);

    coefmat.set_name("ISphericalTransform: Pre-inverted COEFMAT");
    coefmat.print();

    // Invert and transpose the coefficient matrix
    coefmat.invert();

    coefmat.print();
//    temp.gemm(false, false, 1.0, &coefmat, &coefmatI, 0.0);
//    temp.set_name("coefmat * coefmatI");
//    temp.print();

    coefmat.transpose_this();

    coefmat.set_name("ISphericalTransform: COEFMAT");
    coefmat.print();

    // Go through and grab the values.
    int pureindex = 0;
    int cartindex = 0;

    for (int i=1; i<=(l_-subl_)/2; ++i)
        pureindex += npure(subl_+2*i);

    //printf("l_ = %d, subl_ = %d\n", l_, subl_);

    for (int p=0; p<npure(subl_); ++p) {
        cartindex = 0;
        for (int ii=0; ii<=l_; ++ii) {
            int c = l_ - ii;
            for (int jj=0; jj<=ii; ++jj) {
                int b = ii - jj;
                int a = jj;
//        for (int a=0; a<=l_; ++a) {
//            for (int b=0; (a+b)<=l_; ++b) {
//                int c = l_ - a -b;

                int cart = icart(a, b, c);

                double coef = coefmat(cart, p+pureindex);

                if (fabs(coef) > 1.0e-16) {
                    SphericalTransformComponent component;
                    component.init(a, b, c, coef, cartindex, p);
                    components_.push_back(component);
                }
                cartindex++;
            }
        }
    }

    std::vector<SphericalTransformComponent>::const_iterator iter = components_.begin();
    for (; iter != components_.end(); ++iter) {
        const SphericalTransformComponent& ip = *iter;
        fprintf(outfile, "%lf %d %d %d %d\n", ip.coef(), ip.pureindex(), ip.a(), ip.b(), ip.c());
    }
}
