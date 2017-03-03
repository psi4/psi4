/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*!
** \file
** \brief Simple root finding methods
** \ingroup QT
**
**      (1) Bisection method
**      (2) Newton's method
**      (3) Secant method
**
** David Sherrill
** 20 Jan 1994
**
** We know the maximum error for the bisection method, so this is used as 
** the convergence criterion.  However, we only know the general behavior
** of the error for the other two methods; we do not have an absolute value
** bounding it at each step.  Therefore use the difference between
** successive guesses as a convergence criterion in those two cases.
**
** Note also how these three routines access the mathematical function in
** question.  Rather than using an inline function (which would be faster,
** but would require recompilation each time), we pass the function to be
** solved as an argument to the three root-finding routines.  This means
** that this module could be added to a library and it would not need to
** be recompiled for each specific function to be solved.
**
*/


/* INCLUDES */
#include <cstdio>
#include <cmath>

namespace psi {

/*!
** bisect(): Finds the root of a function between two points to within a 
** given tolerance.  Iterations are limited to a given maximum, and
** a print flag specifies whether this function should print the results
** of each iteration to stdout.  Note that the values of the function at
** the two endpoints of the interval must have _different_ signs for the
** bisection method to work!  This routine checks for this initially.
**
** \param function  = pointer to function we want to examine 
**                    (must return double)
** \param low       = lower bound of interval to search for root
** \param high      = upper bound of interval
** \param tolerance = how small is the maximum allowable error
** \param maxiter   = maximum number of iterations
** \param printflag = whether or not to print results for each iteration 
**                    (1 or 0)
**
** Returns: the value of the root
** \ingroup QT
*/
double bisect(double (*function)(double), double low, double high, 
  double tolerance, int maxiter, int printflag)
{
   int iter = 1 ;
   double tmpval ;
   double center ;      /* half way between low and high, i.e. current guess */
   double fl, fc, fh ;  /* values of function at the three points */
   double errorval ; 

   /* check the input parameters */
   if (low > high) {  /* make sure low is the smaller of the 2 limits */
      tmpval = low ;
      low = high ; 
      high = tmpval ;
      }
 
   fl = function(low) ;
   fh = function(high) ;
   if (fl * fh > 0.0) {
      printf("(bisect): Product of function values is greater than 0\n") ;
      return(0.0) ;
      }

   if (printflag) {
      printf("Bisection method root finder\n") ;
      printf("%4s %12s  %12s\n", "iter", "x", "max error") ;
      }

   /* now iterate until we hit maxiter or until we get within the tolerance */
   while (iter < maxiter) {  
      center = (low + high) / 2.0 ;
      fc = function(center) ;
      errorval = (high - low) / 2.0 ;
      if (printflag) 
         printf("%4d %12.8f  %12.8f\n", iter, center, errorval) ;
      if (errorval < tolerance) return(center) ;
      if (fl * fc > 0.0) {
         low = center ;
         fl = fc ;
         }
      else {
         high = center ;
         fh = fc ;
         }
      iter++ ;
      }

   printf("(bisect): Failed to converge within %d iterations\n", iter) ;
   return(center) ;
} 


/*!
** newton(): Find the root of a function by Newton's method.  Iterations are
** limited to a maximum value.  The algorithm stops when the difference
** between successive estimates of the root is less than the specified
** tolerance.  An initial guess for the root must be given, as well as
** the function AND it's derivative.
**
** \param F         = pointer to function we want to examine 
**                    (must return double)
** \param dF        = pointer to _derivative_ of function F
** \param x         = initial guess for root
** \param tolerance = how close successive guesses must get before convergence
** \param maxiter   = maximum number of iterations
** \param printflag = whether or not to print results for each iteration 
**                    (1 or 0)
**
** Returns: the value of the root
** \ingroup QT
*/
double newton(double (*F)(double), double (*dF)(double), double x, 
  double tolerance, int maxiter, int printflag)
{
   int iter = 1 ;
   double newx ;

   if (printflag) {
      printf("Newton's method root finder\n") ;
      printf("%4s %12s\n", "iter", "x") ;
      printf("%4d %12.8f\n", iter, x) ;
      }

   while (iter < maxiter) {
      newx = x - F(x) / dF(x) ;
      iter++ ;
      if (printflag) printf("%4d %12.8f\n", iter, newx) ;
      if (fabs(newx - x) < tolerance) return (newx) ;
      x = newx ;
      }

   printf("(newton): Failed to converge within %d iterations\n", iter) ;
   return(newx) ;
}


/*!
** secant(): Find the root of a function by the Secant Method.  Iterations are
** limited to a maximum value.  The algorithm stops when the relative 
** difference between successive guesses is less than the specified
** tolerance.  An initial TWO guesses for the root must be given, as well
** as the function itself.
**
** \param F         = pointer to function we want to examine 
**                    (must return double)
** \param x0        = 1st guess for root
** \param x1        = 2nd guess for root
** \param tolerance = how close successive guesses must get before convergence
** \param maxiter   = maximum number of iterations
** \param printflag = whether or not to print results for each iteration 
**                    (1 or 0)
**
** Returns: the value of the root
** \ingroup QT
*/
double secant(double (*F)(double), double x0, double x1, double tolerance, 
  int maxiter, int printflag)
{
   int iter = 2 ;
   double x2=0.0 ;
   double fx1, fx0 ;

   if (printflag) {
      printf("Secant Method root finder\n") ;
      printf("%4s %12s\n", "iter", "x") ;
      printf("%4d %12.8f\n", 1, x0) ;
      printf("%4d %12.8f\n", 2, x1) ;
      }


   while (iter < maxiter) {
      fx1 = F(x1) ; 
      fx0 = F(x0) ;
      x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0) ;
      iter++ ;
      if (printflag) printf("%4d %12.8f\n", iter, x2) ;
      if (fabs((x2 - x1)/x1) < tolerance) return(x2) ;
      x0 = x1 ;
      x1 = x2 ;
      }

   printf("(secant): Failed to converge within %d iterations\n", iter) ;
   return(x2) ;
}

}
