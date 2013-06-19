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
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_lag(void);
void uhf_lag(void);
void sf_lag(void);
void Iij(void);
void Iab(void);
void Iai(void);
void Iia(void);

void lag(void)
{
  if(params.gradient) sf_lag();
  else {
    if(params.ref == 0) rhf_lag();
    else if(params.ref == 2) uhf_lag();
  }
}

void rhf_lag(void)
{
  dpdfile2 D; 
  dpdfile2 L; 
  dpdbuf4 I;  
  dpdbuf4 T;  

  dpd_->file2_init(&L, PSIF_CC_OEI, 0, 1, 0, "LAI");
  dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
  dpd_->dot24(&D, &I, &L, 0, 0, 2.0, 0.0); 
  dpd_->dot23(&D, &I, &L, 0, 0, -1.0, 1.0);
  dpd_->buf4_close(&I);
  dpd_->file2_close(&D);
  dpd_->buf4_init(&I, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
  dpd_->dot24(&D, &I, &L, 0, 1, 2.0, 1.0);
  dpd_->dot23(&D, &I, &L, 0, 1, -1.0, 1.0);
  dpd_->buf4_close(&I);
  dpd_->file2_close(&D);
  dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
  dpd_->buf4_init(&I, PSIF_CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk>");
  dpd_->contract442(&T, &I, &L, 2, 0, -2.0, 1.0);
  dpd_->buf4_close(&I);
  dpd_->buf4_init(&I, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
  dpd_->contract442(&I, &T, &L, 0, 0, 2.0, 1.0);
  dpd_->buf4_close(&I);
  dpd_->buf4_close(&T);
  dpd_->file2_close(&L);
}

void uhf_lag(void)
{

}

void sf_lag(void)
{
  dpdfile2 I;

  Iij();
  Iab();
  Iai();
  Iia();

   /* Multiply all I'pq components by -1/2 for compatibility with the
     final gradient expression */
    
  if(params.ref == 0) {
    dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");
    dpd_->file2_scm(&I, -0.5);
    dpd_->file2_close(&I);
    dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");
    dpd_->file2_scm(&I, -0.5);
    dpd_->file2_close(&I);
    dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 1, "I'AB");
    dpd_->file2_scm(&I, -0.5);
    dpd_->file2_close(&I);
    dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 1, "I'ab");
    dpd_->file2_scm(&I, -0.5);
    dpd_->file2_close(&I);
    dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'IA");
    dpd_->file2_scm(&I, -0.5);
    dpd_->file2_close(&I);
    dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'ia");
    dpd_->file2_scm(&I, -0.5);
    dpd_->file2_close(&I);
    dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'AI");
    dpd_->file2_scm(&I, -0.5);
    dpd_->file2_close(&I);
    dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 0, "I'ai");
    dpd_->file2_scm(&I, -0.5);
    dpd_->file2_close(&I);
  }

}

}} /* End namespaces */
