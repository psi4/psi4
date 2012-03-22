#ifndef FUNCTIONAL_UTILITY
#define FUNCTIONAL_UTILITY

#include <psiconfig.h>

#ifndef HAVE_FUNC_ERF
double erf(double x);
double erfc(double x);
#endif

/// heaviside(x) = 1.0 if x >  0
//                 0.0 if x <= 0
// for matlab piecewise functions
inline double heaviside(double x)
{
    return (x > 0.0 ? 1.0 : 0.0);
}

/// dirac(...) = 0.0 for all x
// for matlab piecewise functions
inline double dirac(double x, ...)
{
    return 0.0;
}
double Ei(double x);

#endif
