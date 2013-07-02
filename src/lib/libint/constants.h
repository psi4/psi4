/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup INT
    \brief Enter brief description of file here 
*/

#ifndef _psi3_src_lib_libint_constants_h_
#define _psi3_src_lib_libint_constants_h_

  static int io[] = {0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153,171,190,210,231,253,276,300,325,351,378,406,435,465};
  static const char am_letter[] = "0pdfghiklmnoqrtuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  static const char *number[] = {"zero","one","two","three","four","five","six","seven","eight","nine","ten","eleven",
	  		       "twelve","thirteen","fourteen","fifteen","sixteen","seventeen","eighteen","nineteen","twenty",
                               "twentyone", "twentytwo", "twentythree", "twentyfour", "twentyfive", "twentysix", "twentyseven",
                               "twentyeight", "twentynine", "thirty"};

  /*----------------------------------------------------------------------------------
    hash(a,b) returns a number of the (a[0] a[1]) type product within a doublet.
    a contains x y and z exponents of functions on centers A and B, and b contains
    their angular momenta
   ----------------------------------------------------------------------------------*/

  static int hash(int a[2][3], int b[2])
   {
     int c[2] = {0,0};
     int i;
     
     if(b[0]){
       i=b[0]-a[0][0];
       c[0]=i+io[i]-a[0][1];
     }
     if(b[1]){
       i=b[1]-a[1][0];
       c[1]=i+io[i]-a[1][1];
     }
     
     return c[0]*io[b[1]+1]+c[1];
   }


#endif
