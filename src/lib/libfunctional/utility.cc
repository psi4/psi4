#include <cmath>           // required for fabsl(), expl() and logl()        
#include <cfloat>          // required for LDBL_EPSILON, DBL_MAX
#include <limits>
#include <cstdio>
#include "utility.h"

// Error Function alternative
#ifndef HAVE_FUNC_ERF
double erf(double x)
{
    throw psi::PSIEXCEPTION("erf not implemented, get a C99 compiler to run this functional.");
    return 0.0; 
}
double erfc(double x)
{
    throw psi::PSIEXCEPTION("erf not implemented, get a C99 compiler to run this functional.");
    return 0.0; 
}
#endif

//                         Internally Defined Routines                        //
long double xExponential_Integral_Ei( long double x );

static long double Continued_Fraction_Ei( long double x );
static long double Power_Series_Ei( long double x );
static long double Argument_Addition_Series_Ei( long double x);


//                         Internally Defined Constants                       //
static const long double epsilon = 10.0 * LDBL_EPSILON;

static const double expei_a1 = 4.03640E0;
static const double expei_a2 = 1.15198E0;
static const double expei_a3 = 5.03627E0;
static const double expei_a4 = 4.19160E0;

////////////////////////////////////////////////////////////////////////////////
// double expei( double x )                                                   //
//                                                                            //
//  Description:                                                              //
//      exp(x) Ei(-x) for x > expei_cutoff (700.0)                            //
//  Thanks to Gus Scuseria for this one
////////////////////////////////////////////////////////////////////////////////

double expei(double x)
{
    double val;

    if (x == (std::numeric_limits<double>::infinity())) {
        val = 0.0;
    } else if (x == (-std::numeric_limits<double>::infinity())) {
        val = std::numeric_limits<double>::infinity(); 
    } else {
        val = -1.0 / x * ((x * x + expei_a1 * x + expei_a2) / (x * x + expei_a3 * x + expei_a4));
    } 

    return val;
}


////////////////////////////////////////////////////////////////////////////////
// double Ei( double x )                                                      //
//                                                                            //
//  Description:                                                              //
//     The exponential integral Ei(x) is the integral with integrand          //
//                             exp(t) / t                                     //
//     where the integral extends from -inf to x.                             //
//     Note that there is a singularity at t = 0.  Therefore for x > 0, the   //
//     integral is defined to be the Cauchy principal value:                  //
//          lim { I[-inf, -eta] exp(-t) dt / t + I[eta, x] exp(-t) dt / t }   //
//     in which the limit is taken as eta > 0 approaches 0 and I[a,b]         //
//     denotes the integral from a to b.                                      //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
//     If x = 0.0, then Ei is -inf and -DBL_MAX is returned.                  //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = _Ei( x );                                                          //
////////////////////////////////////////////////////////////////////////////////
double Ei(double x)
{
    double val;
    if (x == (-std::numeric_limits<double>::infinity())) {
        val = 0.0;
    } else if (x == (std::numeric_limits<double>::infinity())) {
        val = std::numeric_limits<double>::infinity(); 
    } else {
        val = (double) xExponential_Integral_Ei( (long double) x);
    }
    //printf("%24.16E %24.16E\n", x, val);
    return val;
}


////////////////////////////////////////////////////////////////////////////////
// long double xExponential_Integral_Ei( long double x )                      //
//                                                                            //
//  Description:                                                              //
//     The exponential integral Ei(x) is the integral with integrand          //
//                             exp(t) / t                                     //
//     where the integral extends from -inf to x.                             //
//     Note that there is a singularity at t = 0.  Therefore for x > 0, the   //
//     integral is defined to be the Cauchy principal value:                  //
//          lim { I[-inf, -eta] exp(-t) dt / t + I[eta, x] exp(-t) dt / t }   //
//     in which the limit is taken as eta > 0 approaches 0 and I[a,b]         //
//     denotes the integral from a to b.                                      //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the exponential integral Ei().         //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
//     If x = 0.0, then Ei is -inf and -DBL_MAX is returned.                  //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xExponential_Integral_Ei( x );                                     //
////////////////////////////////////////////////////////////////////////////////

long double xExponential_Integral_Ei( long double x )
{
   if ( x < -5.0L ) return Continued_Fraction_Ei(x);
   if ( x == 0.0L ) return -DBL_MAX;
   if ( x < 6.8L )  return Power_Series_Ei(x);
   if ( x < 50.0L ) return Argument_Addition_Series_Ei(x);
   return Continued_Fraction_Ei(x);
}

////////////////////////////////////////////////////////////////////////////////
// static long double Continued_Fraction_Ei( long double x )                  //
//                                                                            //
//  Description:                                                              //
//     For x < -5 or x > 50, the continued fraction representation of Ei      //
//     converges fairly rapidly.                                              //
//                                                                            //
//     The continued fraction expansion of Ei(x) is:                          //
//        Ei(x) = -exp(x) { 1/(-x+1-) 1/(-x+3-) 4/(-x+5-) 9/(-x+7-) ... }.    //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////

static long double Continued_Fraction_Ei( long double x )
{
   long double Am1 = 1.0L;
   long double A0 = 0.0L;
   long double Bm1 = 0.0L;
   long double B0 = 1.0L;
   long double a = expl(x);
   long double b = -x + 1.0L;
   long double Ap1 = b * A0 + a * Am1;
   long double Bp1 = b * B0 + a * Bm1;
   int j = 1;

   a = 1.0L;
   while ( fabsl(Ap1 * B0 - A0 * Bp1) > epsilon * fabsl(A0 * Bp1) ) {
      if ( fabsl(Bp1) > 1.0L) {
         Am1 = A0 / Bp1;
         A0 = Ap1 / Bp1;
         Bm1 = B0 / Bp1;
         B0 = 1.0L;
      } else {
         Am1 = A0;
         A0 = Ap1;
         Bm1 = B0;
         B0 = Bp1;
      }
      a = -j * j;
      b += 2.0L;
      Ap1 = b * A0 + a * Am1;
      Bp1 = b * B0 + a * Bm1;
      j += 1;
   }
   return (-Ap1 / Bp1);
}


////////////////////////////////////////////////////////////////////////////////
// static long double Power_Series_Ei( long double x )                        //
//                                                                            //
//  Description:                                                              //
//     For -5 < x < 6.8, the power series representation for                  //
//     (Ei(x) - gamma - ln|x|)/exp(x) is used, where gamma is Euler's gamma   //
//     constant.                                                              //
//     Note that for x = 0.0, Ei is -inf.  In which case -DBL_MAX is          //
//     returned.                                                              //
//                                                                            //
//     The power series expansion of (Ei(x) - gamma - ln|x|) / exp(x) is      //
//        - Sum(1 + 1/2 + ... + 1/j) (-x)^j / j!, where the Sum extends       //
//        from j = 1 to inf.                                                  //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////

static long double Power_Series_Ei( long double x )
{ 
   long double xn = -x;
   long double Sn = -x;
   long double Sm1 = 0.0L;
   long double hsum = 1.0L;
   long double g = 0.5772156649015328606065121L;
   long double y = 1.0L;
   long double factorial = 1.0L;
  
   if ( x == 0.0L ) return (long double) -DBL_MAX;
 
   while ( fabsl(Sn - Sm1) > epsilon * fabsl(Sm1) ) {
      Sm1 = Sn;
      y += 1.0L;
      xn *= (-x);
      factorial *= y;
      hsum += (1.0 / y);
      Sn += hsum * xn / factorial;
   }
   return (g + logl(fabsl(x)) - expl(x) * Sn);
}


////////////////////////////////////////////////////////////////////////////////
// static long double Argument_Addition_Series_Ei(long double x)              //
//                                                                            //
//  Description:                                                              //
//     For 6.8 < x < 50.0, the argument addition series is used to calculate  //
//     Ei.                                                                    //
//                                                                            //
//     The argument addition series for Ei(x) is:                             //
//     Ei(x+dx) = Ei(x) + exp(x) Sum j! [exp(j) expj(-dx) - 1] / x^(j+1),     //
//     where the Sum extends from j = 0 to inf, |x| > |dx| and expj(y) is     //
//     the exponential polynomial expj(y) = Sum y^k / k!, the Sum extending   //
//     from k = 0 to k = j.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the exponential integral Ei().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the exponential integral Ei evaluated at x.               //
////////////////////////////////////////////////////////////////////////////////
static long double Argument_Addition_Series_Ei(long double x)
{
   static long double ei[] = {
      1.915047433355013959531e2L,  4.403798995348382689974e2L,
      1.037878290717089587658e3L,  2.492228976241877759138e3L,
      6.071406374098611507965e3L,  1.495953266639752885229e4L,
      3.719768849068903560439e4L,  9.319251363396537129882e4L,
      2.349558524907683035782e5L,  5.955609986708370018502e5L,
      1.516637894042516884433e6L,  3.877904330597443502996e6L,
      9.950907251046844760026e6L,  2.561565266405658882048e7L,
      6.612718635548492136250e7L,  1.711446713003636684975e8L,
      4.439663698302712208698e8L,  1.154115391849182948287e9L,
      3.005950906525548689841e9L,  7.842940991898186370453e9L,
      2.049649711988081236484e10L, 5.364511859231469415605e10L,
      1.405991957584069047340e11L, 3.689732094072741970640e11L,
      9.694555759683939661662e11L, 2.550043566357786926147e12L,
      6.714640184076497558707e12L, 1.769803724411626854310e13L,
      4.669055014466159544500e13L, 1.232852079912097685431e14L,
      3.257988998672263996790e14L, 8.616388199965786544948e14L,
      2.280446200301902595341e15L, 6.039718263611241578359e15L,
      1.600664914324504111070e16L, 4.244796092136850759368e16L,
      1.126348290166966760275e17L, 2.990444718632336675058e17L,
      7.943916035704453771510e17L, 2.111342388647824195000e18L,
      5.614329680810343111535e18L, 1.493630213112993142255e19L,
      3.975442747903744836007e19L, 1.058563689713169096306e20L
   };
   int  k = (int) (x + 0.5);
   int  j = 0;
   long double xx = (long double) k;
   long double dx = x - xx;
   long double xxj = xx;
   long double edx = expl(dx);
   long double Sm = 1.0L;
   long double Sn = (edx - 1.0L) / xxj;
   long double term = DBL_MAX;
   long double factorial = 1.0L;
   long double dxj = 1.0L;

   while (fabsl(term) > epsilon * fabsl(Sn) ) {
      j++;
      factorial *= (long double) j;
      xxj *= xx;
      dxj *= (-dx);
      Sm += (dxj / factorial);
      term = ( factorial * (edx * Sm - 1.0L) ) / xxj;
      Sn += term;
   }
   
   return ei[k-7] + Sn * expl(xx); 
}

