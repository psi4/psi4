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

#include <cstdio>
#include <cstdlib> /* was libc.h */
/* gcc 2.7.0 doesn't like #include <cstring> */
#include "psi4/detci/slaterd.h"

#include "psi4/psi4-dec.h"

namespace psi { namespace detci {

extern double get_twoel(int i, int j, int k, int l);
extern double get_onel(int i, int j);
extern int calc_orb_diff(int cnt, unsigned char *I, unsigned char *J,
   int *I_alpha_diff, int *J_alpha_diff, int *sign, int *same,
   int extended);
extern void common_orbs(int *same_alpha, int *same_beta, int cnt_alpha,
   int cnt_beta, int *common_docc, int *common_alpha_socc,
   int *common_beta_socc, int *cnt_docc, int *cnt_alpha_socc,
   int *cnt_beta_socc);


void SlaterDeterminant::set(unsigned int na, unsigned char *alpoccs,
      unsigned int nb, unsigned char *betoccs)
{
   int i;

   if (nalp_ != na) {
      if (Occs_[0] != NULL) free(Occs_[0]);
      Occs_[0] = (unsigned char *) malloc (sizeof(unsigned char) * na);
      nalp_ = na;
      }
   if (nbet_ != nb) {
      if (Occs_[1] != NULL) free(Occs_[1]);
      Occs_[1] = (unsigned char *) malloc (sizeof(unsigned char) * nb);
      nbet_ = nb;
      }

   for (i=0; i<nalp_; i++) {
      Occs_[0][i] = alpoccs[i];
      }
   for (i=0; i<nbet_; i++) {
      Occs_[1][i] = betoccs[i];
      }
}



void SlaterDeterminant::print()
{
   int i;

   outfile->Printf( "Alpha string: ");
   for (i=0; i<nalp_; i++) {
      outfile->Printf( "%3d ", Occs_[0][i]);
      }
   outfile->Printf( "\n");

   outfile->Printf( "Beta string : ");
   for (i=0; i<nbet_; i++) {
      outfile->Printf( "%3d ", Occs_[1][i]);
      }
   outfile->Printf( "\n");
}


void SlaterDeterminant::print_config()
{
   int i=0, j=0;

   while ((i < nalp_) && (j < nbet_)) {
      if (Occs_[0][i] == Occs_[1][j]) {
         outfile->Printf( "%dX ", Occs_[0][i]+1);
         i++; j++;
         }
      else if (Occs_[0][i] < Occs_[1][j]) {
         outfile->Printf( "%dA ", Occs_[0][i]+1);
         i++;
         }
      else if (Occs_[0][i] > Occs_[1][j]) {
         outfile->Printf( "%dB ", Occs_[1][j]+1);
         j++;
         }
      }

   if (i < j) {
      while (i < nalp_) {
         outfile->Printf( "%dA ", Occs_[0][i]+1);
         i++;
         }
      }
   else if (i > j) {
      while (j < nbet_) {
         outfile->Printf( "%dB ", Occs_[1][j]+1);
         j++;
         }
      }

   outfile->Printf( "\n") ;

}


SlaterDeterminant& SlaterDeterminant::operator=(const SlaterDeterminant& s)
{
   if (nalp_ != s.nalp_) {
      if (Occs_[0] != NULL) free(Occs_[0]);
      Occs_[0] = (unsigned char *) malloc (sizeof(unsigned char) * s.nalp_);
      }
   if (nbet_ != s.nbet_) {
      if (Occs_[1] != NULL) free(Occs_[1]);
      Occs_[1] = (unsigned char *) malloc (sizeof(unsigned char) * s.nbet_);
      }
   set(s.nalp_, s.Occs_[0], s.nbet_, s.Occs_[1]);
   return(*this);
}


int operator ==(SlaterDeterminant& s1, SlaterDeterminant& s2)
{
   int i;

   if (s1.nalp_ != s2.nalp_ || s1.nbet_ != s2.nbet_) return(0);

   for (i=0; i<s1.nalp_; i++) {
      if (s1.Occs_[0][i] != s2.Occs_[0][i]) return(0);
   }
   for (i=0; i<s1.nbet_; i++) {
      if (s1.Occs_[1][i] != s2.Occs_[1][i]) return(0);
   }

   return(1);
}



}} // namespace psi::detci
