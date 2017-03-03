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
  \file
  \brief Contains some probability functions
  \ingroup QT
*/

namespace psi {
	
/*!
** combinations() : Calculates the number of ways to choose k objects
**    from n objects, or "n choose k" 
**
** Parameters:
**   \param n   =  number of objects in total
**   \param k   =  number of objects taken at a time
**
** Returns:
**    number of combinations of n objects taken k at a time ("n choose k")
**    (returned as a double).
**
** \ingroup QT
*/
double combinations(int n, int k)
{
   double factorial(int) ;
   double comb ;

   if (n == k) return (1.0) ;
   else if (k > n) return(0.0) ;
   else if (k == 0) return(1.0) ; 
   comb = factorial(n) / (factorial(k) * factorial(n-k)) ;
 
   return(comb) ;
}


/*!
** factorial(): Returns n!
** 
** Parameters:
**    \param n  = number to take factorial of
**
** Returns:
**    n factorial, as a double word (since n! can get very large).
** \ingroup QT
*/
double factorial(int n)
{

   if (n == 0 || n == 1) return(1.0);
   if (n < 0) return(0.0) ;
   else {
      return ((double) n * factorial(n-1)) ;
      }
}



/*
** test combinations routines
**
#include <cstdio>

main()
{
   int i, j ;
   double factorial() ;
   double combinations() ;

   printf("Enter two numbers: ") ;
   scanf("%d %d", &i, &j) ;
   printf("i! = %.2f\n", factorial(i)) ;
   printf("j! = %.2f\n", factorial(j)) ;
   printf("i choose j = %.2f\n", combinations(i,j)) ;
}
**
*/

}
