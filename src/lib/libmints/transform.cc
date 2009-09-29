#include <cmath>

#include <libmints/wavefunction.h>

#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/integral.h>

#include <psiconfig.h>

#include <psi4-dec.h>

using namespace psi;

#define parity(m) ((m)%2 ? -1 : 1) // returns (-1)^m
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

double xyz2lm_coeff(int l, int m, int lx, int ly, int lz)
{
    static int use_cca_integrals_standard = (PSI_INTEGRALS_STANDARD == 1);
    int i, j, k, i_max;
    int k_min, k_max;
    int abs_m;
    int comp;
    double pfac, pfac1, sum, sum1;

    abs_m = abs(m);
    if ((lx + ly - abs(m))%2)
        return 0.0;
    else
        j = (lx + ly - abs(m))/2;

    if (j < 0)
        return 0.0;

    comp = (m >= 0) ? 1 : -1;
    i = abs_m-lx;
    if (comp != parity(i))
        return 0.0;

    pfac = sqrt(fac[2*lx]*fac[2*ly]*fac[2*lz]*fac[l-abs_m]/(fac[2*l]*fac[l]
        *fac[lx]*fac[ly]*fac[lz]*fac[l+abs_m]));
    pfac /= (1 << l);

    if (m < 0)
        pfac *= parity((i-1)/2);
    else
        pfac *= parity(i/2);

    i_max = (l-abs_m)/2;
    sum = 0.0;
    for (i=0; i<=i_max; i++) {
        pfac1 = bc[l][i]*bc[i][j];
        if (pfac1 == 0.0)
            continue;
        else
            pfac1 *= (parity(i)*fac[2*(l-i)]/fac[l-abs_m-2*i]);
        sum1 = 0.0;
        k_min = MAX((lx-abs_m)/2,0);
        k_max = MIN(j,lx/2);
        for (k=k_min; k<=k_max; k++)
            sum1 += bc[j][k]*bc[abs_m][lx-2*k]*parity(k);
        sum += pfac1*sum1;
    }

    if (use_cca_integrals_standard)
        sum *= sqrt(df[2*l]/(df[2*lx]*df[2*ly]*df[2*lz]));

    if (m == 0)
        return pfac*sum;
    else
        return M_SQRT2*pfac*sum;
}

void SphericalTransformComponent::init(int a, int b, int c, double coef, int cartindex, int pureindex)
{
    // printf("a = %d, b = %d, c = %d, coef = %f, cartindex = %d, pureindex = %d\n", a, b, c, coef, cartindex, pureindex);
    a_ = a;
    b_ = b;
    c_ = c;
    coef_ = coef;
    cartindex_ = cartindex;
    pureindex_ = pureindex;
}

SphericalTransform::SphericalTransform(int l)
{    
    int pureindex = 0;
    int cartindex = 0;

    l_ = l;
    
    // Go through all the spherical harmonic possibilities
    for (int m = -l_; m <= l_; ++m) {
        // Generate all the possible cartesian functions
        cartindex = 0;
        for(int ii = 0; ii <= l_; ii++) {
            int l1 = l_ - ii;
            for(int jj = 0; jj <= ii; jj++) {
                int m1 = ii - jj;
                int n1 = jj;
                
                // Compute the coefficient needed to convert from cartesian to spherical
                double coef = xyz2lm_coeff(l, m, l1, m1, n1);
                
                // Save the information if it is nonzero
                if (fabs(coef) > 1.0e-16) {
                    SphericalTransformComponent component;
                    component.init(l1, m1, n1, coef, cartindex, pureindex);
                    components_.push_back(component);
                }
                
                cartindex++;
            }
        }
        pureindex++;
    }
}
