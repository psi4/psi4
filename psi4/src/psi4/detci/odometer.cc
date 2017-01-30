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

/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/

/*
** ODOMETER.CC
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia
*/

#include "psi4/detci/odometer.h"
#include "psi4/psi4-dec.h"

#include <cstdio>

namespace psi { namespace detci {

// #define STANDALONE
#ifdef STANDALONE
#include <iostream>
#endif

Odometer::Odometer()
{
   length = 0 ;
   max = NULL ;
   min = NULL ;
   value = NULL ;
}

Odometer::~Odometer()
{
   if (length) {
      delete [] max;
      delete [] min;
      delete [] value;
      }
   length = 0;
}

void Odometer::size(unsigned n)
{
   int i ;

   length = n ;
   max = new int[n] ;
   min = new int[n] ;
   value = new int[n] ;

   for (i=0; i<n; i++) {
      max[i] = 9 ; 
      min[i] = 0 ;
      value[i] = 0 ;
      }
}

void Odometer::resize(unsigned n)
{
   int i ;

   if (length != 0) {
      delete [] max ;
      delete [] min ;
      delete [] value ;
      }
   length = n ;

   if (n) {
      max = new int[n] ;
      min = new int[n] ;
      value = new int[n] ;
      }

   for (i=0; i<n; i++) {
      max[i] = 9 ; 
      min[i] = 0 ;
      value[i] = 0 ;
      }
}


void Odometer::set_max(int m)
{
   int i ;

   for (i=0; i<length; i++)
     max[i] = m ;
}

void Odometer::set_max_lex(int m)
{
   int i ;

   for (i=0; i<length; i++)
      max[i] = m-i ;
}

void Odometer::set_max(int* m)
{
   int i ;

   for (i=0; i<length; i++)
     max[i] = m[i] ;
}

void Odometer::set_min(int m)
{
   int i ;

   for (i=0; i<length; i++)
     min[i] = m ;

}

void Odometer::set_min_lex(int m)
{
   int i,j ;

   for (i=length-1,j=0; i>=0; i--,j++)
      min[i] = m + j; 
}

void Odometer::set_min(int* m)
{
   int i ;

   for (i=0; i<length; i++)
     min[i] = m[i] ;
}


void Odometer::set_value(int m)
{
   int i ;

   for (i=0; i<length; i++)
     value[i] = m ;
}


void Odometer::get_value(int* m)
{
   int i ;

   for (i=0; i<length; i++)
     m[i] = value[i] ;
}


void Odometer::set_value(int* m)
{
   int i ;

   for (i=0; i<length; i++)
     value[i] = m[i] ;
}


void Odometer::reset()
{
   int i ;

   for (i=0; i<length; i++)
      value[i] = min[i] ;
}


void Odometer::increment()
{
   int i ;

   for (i=0; i<length; i++) {
      if (value[i] < max[i]) {
         value[i] += 1 ;
         return ;
         }
      else {
         value[i] = min[i] ;
         }
      }
}
       

void Odometer::increment_lex()
{
   int i,j ;

   for (i=0; i<length; i++) {
     if (value[i] < max[i]) {
        value[i] += 1 ;
        for (j=i-1; j>=0; j--) {
          if (value[j+1] + 1 >= min[j]) value[j] = value[j+1] + 1;
          else value[j] = min[j] ;
          }
        break ;
        }
    else {
       value[i] = min[i] ;
       }
    }
}


void Odometer::print()
{
   int i ;

   for (i=length-1; i>=0; i--) {
     outfile->Printf("%d ", value[i]);
      }
  outfile->Printf("\n");
}
 

unsigned Odometer::at_max()
{
   int i ;
   unsigned tval=1 ;

   if (length == 0) return (1) ;
   for (i=0; i<length; i++) {
      if (value[i] != max[i]) tval = 0 ;
      }
    return(tval) ;
}


unsigned Odometer::at_min()
{
   int i ;
   unsigned tval=1 ;

   if (length == 0) return (1) ;
   for (i=0; i<length; i++) {
      if (value[i] != min[i]) tval = 0 ;
      }
    return(tval) ;
}


unsigned Odometer::boundscheck()
{
   int i ;

   if (length == 0) return(1) ;

   for (i=0; i<length; i++)
      if ((max[i] < min[i]) || (min[i] > max[i])) return(0) ;

   return(1) ;
}


#ifdef STANDALONE
main()
{
Odometer od ;
int maxvals[] = {7, 6, 5} ;
int minvals[] = {3, 2, 1} ;

   od.size(3) ;
   od.print() ;
   // od.set_min(minvals) ;
   od.set_min_lex(1) ;
   od.reset() ;
   // od.set_max(maxvals) ;
   od.set_max_lex(7) ;
   while (!od.at_max()) {
      od.print(); 
      od.increment_lex() ;
      }
   od.print() ;
}
#endif

}} // namespace psi::detci
