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
** ODOMETER(): Generalized odometer object.  Each `digit' can have its
**   own min and max values (and actually each position can be filled by an
**   integer, not just a single digit).  The low-index digits are the 
**   faster moving ones.  There is also a special lexical increment which
**   ensures that digit i is greater than digit i+1.
**
**   Methods:
**      Odometer(): Default constructor gives 0 digits and all NULL arrays.
**      Odometer(int n): Construct an odometer with n digits
**      ~Odometer(): Free's the dynamically allocated memory
**      size(int m): Set the length of a default-constructed odometer to m
**      resize(int m): Change the length of an odometer to m
**      set_max(int m): set the maximum for every digit to m
**      set_max(int *m): set the max for each digit according to array m
**      set_min(int m), set_min(int *m): similar to set_max()
**      set_value(int m), set_value(int *m): similar to set_max()
**      set_max_lex(int m): set max values for each digit such that the
**         lowest index digit has max m, the next has max m-1, etc.
**      set_min_lex(int m): set min values for each digit such that the
**         highest index digit has min m, the next has min m-1, etc.
**      get_value(int *m): copy the values array into array m
**      increment(): increment the odometer (at position 0)
**      increment_lex(): lexical index...increment but make sure that
**         value[i] > value[i+1]
**      reset(): return all digits to their min values
**      print(): print out odometer reading (high index digits first)
**      at_max(): has the odometer reached its maximum reading? (1 or 0)
**      at_min(): has the odometer reached its minimum reading? (1 or 0)
**      get_length(): returns the length (number of digits) of the odometer
**      boundscheck(): make sure max's are > min's, and vice versa
**
*/

#ifndef _psi_src_bin_detci_odometer_h
#define _psi_src_bin_detci_odometer_h

namespace psi { namespace detci {

class Odometer {

   protected:
      unsigned length ;
      int* max ;
      int* min ;
      int* value ;
       
   public:
      Odometer() ;
      Odometer(unsigned len) { size(len); }
      ~Odometer() ; 
      
      void size(unsigned s) ;
      void resize(unsigned s) ;
      void set_max(int m) ;
      void set_max_lex(int m) ;
      void set_max(int* m) ;
      void set_min(int m) ;
      void set_min_lex(int m) ;
      void set_min(int* m) ;
      void set_value(int m) ;
      void set_value(int* m) ;
      void get_value(int* m) ;
      void increment() ;
      void increment_lex() ;
      void reset() ;
      void print() ;
      unsigned at_max() ;
      unsigned at_min() ;
      unsigned get_length() 
         {return length; }
      unsigned boundscheck() ;
} ;

}} // namespace psi::detci

#endif // header guard