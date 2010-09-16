#include "wavefunction.h"
#include "matrix.h"
#include "integral.h"

#include <psiconfig.h>
#include <psi4-dec.h>

#include <cmath>

using namespace psi;

extern void solidharmonic(int l, SimpleMatrix& coefmat);

void SphericalTransformComponent::init(int a, int b, int c, double coef,
                                       int cartindex, int pureindex)
{
    //printf("a = %d, b = %d, c = %d, coef = %f, cartindex = %d, pureindex = %d\n", a, b, c, coef, cartindex, pureindex);
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

    // Go through and grab the values.
    int pureindex = 0;
    int cartindex = 0;

    for (int i=1; i<=(l_-subl_)/2; ++i)
        pureindex += INT_NPURE(subl_+2*i);

    //printf("l_ = %d, subl_ = %d\n", l_, subl_);

    for (int p=0; p<INT_NPURE(subl_); ++p) {
        cartindex = 0;
        for (int ii=0; ii<=l_; ++ii) {
            int a = l_ - ii;
            for (int jj=0; jj<=ii; ++jj) {
                int b = ii - jj;
                int c = jj;
                int cart = INT_ICART(a, b, c);

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

}

void ISphericalTransform::init()
{
    int cartdim = INT_NCART(l_);
    SimpleMatrix coefmat(cartdim, cartdim);
    coefmat.zero();

    // Compute the solid harmonic matrix elements
    solidharmonic(l_, coefmat);

    // Invert and transpose the coefficient matrix
    coefmat.invert();
    coefmat.transpose_this();

    // Go through and grab the values.
    int pureindex = 0;
    int cartindex = 0;

    for (int i=1; i<=(l_-subl_)/2; ++i)
        pureindex += INT_NPURE(subl_+2*i);

    //printf("l_ = %d, subl_ = %d\n", l_, subl_);

    for (int p=0; p<INT_NPURE(subl_); ++p) {
        cartindex = 0;
        for (int ii=0; ii<=l_; ++ii) {
            int a = l_ - ii;
            for (int jj=0; jj<=ii; ++jj) {
                int b = ii - jj;
                int c = jj;
                int cart = INT_ICART(a, b, c);

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
