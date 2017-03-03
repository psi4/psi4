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
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void norm_HC1(int i, int C_irr) {
  double dot;
  dpdfile2 F;
  dpdbuf4 W;

  if (params.eom_ref == 0) { /* RHF */

    /** FME **/
    global_dpd_->file2_init(&F, PSIF_CC3_HC1, 0, 0, 1, "HC1 FME");
    dot = global_dpd_->file2_dot_self(&F);
    global_dpd_->file2_close(&F);
    outfile->Printf("<FME|FME> = %15.10lf\n", dot); 

    /** WAMEF **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 11, 5, 11, 5, 0, "HC1 WAmEf (Am,Ef)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("<WAmEf (Am,Ef)|WAmEf> = %15.10lf\n", dot); 

    /** WMNIE **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 0, 10, 0, 10, 0, "HC1 WMnIe (Mn,Ie)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("<WMnIe (Mn,Ie)|WMnIe> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 0, 10, 0, 10, 0, "HC1 2WMnIe - WnMIe (Mn,Ie)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("<2WMnIe - WnMIe (Mn,Ie)|2WMnIe - WnMIe> = %15.10lf\n", dot); 

    /** WMNIJ **/

/*
    dpd_buf4_init(&W, CC3_HC1, 0, 0, 0, 0, 0, 0, "HC1 WMnIj (Mn,Ij)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("<WMnIj (Mn,Ij)|WMnIj> = %15.10lf\n", dot); 
*/

    /** WMBIJ **/
/*
    dpd_buf4_init(&W, CC3_HC1, 0, 0, 10, 0, 10, 0, "HC1 WMbIj (Ij,Mb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("<WMbIj (Ij,Mb)|WMbIj> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 10, 0, 10, 0, 0, "HC1 WMbIj (Mb,Ij)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("<WMbIj (Mb,Ij)|WMbIj> = %15.10lf\n", dot);
*/

    /** WMBEJ **/

/*
    dpd_buf4_init(&W, CC3_HC1, 0, 10, 10, 10, 10, 0, "HC1 WMbEj (ME,jb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("<WMbEj (ME,jb)|WMbEj> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 10, 10, 10, 10, 0, "HC1 WMbeJ (Me,Jb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("<WMbeJ (Me,Jb)|WMbeJ> = %15.10lf\n", dot);
*/

    /** WABEI **/

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (Ie,Ab)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("<WAbEi (Ie,Ab) | WAbEi (Ie,Ab)> = %15.10lf\n", dot);

  }

  else if (params.eom_ref == 1) { /* ROHF */

    /** FME **/
    global_dpd_->file2_init(&F, PSIF_CC3_HC1, 0, 0, 1, "HC1 FME");
    dot = global_dpd_->file2_dot_self(&F);
    global_dpd_->file2_close(&F);
    outfile->Printf("\t<FME|FME> = %15.10lf\n", dot); 

    /** Fme **/
    global_dpd_->file2_init(&F, PSIF_CC3_HC1, 0, 0, 1, "HC1 Fme");
    dot = global_dpd_->file2_dot_self(&F);
    global_dpd_->file2_close(&F);
    outfile->Printf("\t<Fme|Fme> = %15.10lf\n", dot); 

    outfile->Printf("Wamef terms\n");
    /** WAMEF **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 11, 5, 11, 7, 0, "HC1 WAMEF (AM,E>F)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WAMEF (AM,E>F)|WAMEF> = %15.10lf\n", dot); 

    /** Wamef **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 11, 5, 11, 7, 0, "HC1 Wamef (am,e>f)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<Wamef (am,e>f)|Wamef> = %15.10lf\n", dot); 

    /** WAmEf **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 11, 5, 11, 5, 0, "HC1 WAmEf (Am,Ef)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WAmEf (Am,Ef)|WAmEf> = %15.10lf\n", dot); 

    /** WaMeF **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 11, 5, 11, 5, 0, "HC1 WaMeF (aM,eF)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WaMeF (aM,eF)|WaMeF> = %15.10lf\n", dot); 

    /** WMNIE **/
    outfile->Printf("Wmnie terms\n");

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 0, 10, 2, 10, 0, "HC1 WMNIE (M>N,IE)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WMNIE (M>N,IE)|WMNIE> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 0, 10, 2, 10, 0, "HC1 Wmnie (m>n,ie)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<Wmnie(m>n,ie)|Wmnie> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 0, 10, 0, 10, 0, "HC1 WMnIe (Mn,Ie)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WMnIe (Mn,Ie)|WMnIe> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 0, 10, 0, 10, 0, "HC1 WmNiE (mN,iE)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WmNiE (mN,iE)|WmNiE> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 0, 11, 2, 11, 0, "HC1 WMNIE (M>N,EI)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WMNIE(M>N,EI)|WMNIE(M>N,EI)> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 0, 11, 2, 11, 0, "HC1 Wmnie (m>n,ei)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<Wmnie(m>n,ei)|Wmnie(m>n,ei)> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 0, 11, 0, 11, 0, "HC1 WMnIe (Mn,eI)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WMnIe(Mn,eI)|WMnIe(Mn,eI)> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 0, 11, 0, 11, 0, "HC1 WmNiE (mN,Ei)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WmNiE(mN,Ei)|WmNiE(mN,Ei)> = %15.10lf\n", dot); 

    /** WMNIJ **/
/*
    dpd_buf4_init(&W, CC3_HC1, 0, 0, 0, 2, 2, 0, "HC1 WMNIJ (M>N,I>J)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMNIJ (M>N,I>J)|WMNIJ> = %15.10lf\n", dot); 

    dpd_buf4_init(&W, CC3_HC1, 0, 0, 0, 2, 2, 0, "HC1 Wmnij (m>n,i>j)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<Wmnij (m>n,i>j)|Wmnij> = %15.10lf\n", dot); 

    dpd_buf4_init(&W, CC3_HC1, 0, 0, 0, 0, 0, 0, "HC1 WMnIj (Mn,Ij)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMnIj (Mn,Ij)|WMnIj> = %15.10lf\n", dot); 
*/

    /** WMBIJ **/

/*
    dpd_buf4_init(&W, CC3_HC1, 0, 10, 0, 10, 2, 0, "HC1 WMBIJ (MB,I>J)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMBIJ (MB,I>J)|WMBIJ> = %15.10lf\n", dot);
    dpd_buf4_init(&W, CC3_HC1, 0, 0, 10, 2, 10, 0, "HC1 WMBIJ (I>J,MB)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMBIJ (I>J,MB)|WMBIJ> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 10, 0, 10, 2, 0, "HC1 Wmbij (mb,i>j)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<Wmbij (mb,i>j)|Wmbij> = %15.10lf\n", dot);
    dpd_buf4_init(&W, CC3_HC1, 0, 0, 10, 2, 10, 0, "HC1 Wmbij (i>j,mb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<Wmbij (i>j,mb)|Wmbij> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 10, 0, 10, 0, 0, "HC1 WMbIj (Mb,Ij)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMbIj (Mb,Ij)|WMbIj> = %15.10lf\n", dot);
    dpd_buf4_init(&W, CC3_HC1, 0, 0, 10, 0, 10, 0, "HC1 WMbIj (Ij,Mb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMbIj (Ij,Mb)|WMbIj> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 10, 0, 10, 0, 0, "HC1 WmBiJ (mB,iJ)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WmBiJ (mB,iJ)|WmBiJ> = %15.10lf\n", dot);
    dpd_buf4_init(&W, CC3_HC1, 0, 0, 10, 0, 10, 0, "HC1 WmBiJ (iJ,mB)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WmBiJ (iJ,mB)|WmBiJ> = %15.10lf\n", dot);
*/

    /** WMBEJ **/

/*
    dpd_buf4_init(&W, CC3_HC1, 0, 10, 10, 10, 10, 0, "HC1 WMBEJ (ME,JB)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMBEJ (ME,JB)|WMBEJ> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 10, 10, 10, 10, 0, "HC1 Wmbej (me,jb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<Wmbej (me,jb)|Wmbej> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 10, 10, 10, 10, 0, "HC1 WMbEj (ME,jb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMbEj (ME,jb)|WMbEj> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 10, 10, 10, 10, 0, "HC1 WmBeJ (me,JB)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WmBeJ (me,JB)|WmBeJ> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 10, 10, 10, 10, 0, "HC1 WMbeJ (Me,Jb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMbeJ (Me,Jb)|WMbeJ> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 10, 10, 10, 10, 0, "HC1 WmBEj (mE,jB)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WmBEj (mE,jB)|WmBEj> = %15.10lf\n", dot);
*/

    /** WABEI **/

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 10, 5, 10, 7, 0, "HC1 WABEI (IE,A>B)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WABEI (IE,A>B)|WABEI> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 10, 5, 10, 7, 0, "HC1 Wabei (ie,a>b)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<Wabei (ie,a>b)|Wabei> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 10, 5, 10, 5, 0, "HC1 WAbEi (iE,Ab)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WAbEi (iE,Ab)|WAbEi> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 10, 5, 10, 5, 0, "HC1 WaBeI (Ie,aB)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WaBeI (Ie,aB)|WaBeI> = %15.10lf\n", dot);
  }
  else { /* UHF */

    outfile->Printf("***** Norms of HC1 *****\n");
    /** FME **/
    global_dpd_->file2_init(&F, PSIF_CC3_HC1, 0, 0, 1, "HC1 FME");
    dot = global_dpd_->file2_dot_self(&F);
    global_dpd_->file2_close(&F);
    outfile->Printf("\t<FME|FME> = %15.10lf\n", dot); 

    /** Fme **/
    global_dpd_->file2_init(&F, PSIF_CC3_HC1, 0, 2, 3, "HC1 Fme");
    dot = global_dpd_->file2_dot_self(&F);
    global_dpd_->file2_close(&F);
    outfile->Printf("\t<Fme|Fme> = %15.10lf\n", dot); 

    outfile->Printf("Wamef terms\n");
    /** WAMEF **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 20, 5, 20, 7, 0, "HC1 WAMEF (MA,F>E)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WAMEF (MA,F>E)|WAMEF> = %15.10lf\n", dot); 

    /** Wamef **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 30, 15, 30, 17, 0, "HC1 Wamef (ma,f>e)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<Wamef (ma,f>e)|Wamef> = %15.10lf\n", dot); 

    /** WAmEf **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 27, 29, 27, 28, 0, "HC1 WAmEf (mA,fE)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WAmEf (mA,fE)|WAmEf> = %15.10lf\n", dot); 

    /** WaMeF **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 24, 28, 24, 28, 0, "HC1 WaMeF (Ma,Fe)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WaMeF (Ma,Fe)|WaMeF> = %15.10lf\n", dot); 

    /** WMNIE **/
    outfile->Printf("Wmnie terms\n");

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 0, 20, 2, 20, 0, "HC1 WMNIE (M>N,IE)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WMNIE (M>N,IE)|WMNIE> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 10, 30, 12, 30, 0, "HC1 Wmnie (m>n,ie)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<Wmnie (m>n,ie)|Wmnie> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 22, 24, 22, 24, 0, "HC1 WMnIe (Mn,Ie)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WMnIe (Mn,Ie)|WMnIe> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 23, 27, 23, 27, 0, "HC1 WmNiE (mN,iE)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WmNiE (mN,iE)|WmNiE> = %15.10lf\n", dot); 

    /*
    dpd_buf4_init(&W, CC3_HC1, 0, 0, 21, 2, 21, 0, "HC1 WMNIE (M>N,EI)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMNIE(M>N,EI)|WMNIE(M>N,EI)> = %15.10lf\n", dot); 

    dpd_buf4_init(&W, CC3_HC1, 0, 10, 31, 12, 31, 0, "HC1 Wmnie (m>n,ei)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<Wmnie(m>n,ei)|Wmnie(m>n,ei)> = %15.10lf\n", dot); 

    dpd_buf4_init(&W, CC3_HC1, 0, 22, 25, 22, 25, 0, "HC1 WMnIe (Mn,eI)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMnIe(Mn,eI)|WMnIe(Mn,eI)> = %15.10lf\n", dot); 

    dpd_buf4_init(&W, CC3_HC1, 0, 23, 26, 23, 26, 0, "HC1 WmNiE (mN,Ei)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WmNiE(mN,Ei)|WmNiE(mN,Ei)> = %15.10lf\n", dot); 
    */

    /** WMNIJ **/
/*
    dpd_buf4_init(&W, CC3_HC1, 0, 0, 0, 2, 2, 0, "HC1 WMNIJ (M>N,I>J)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMNIJ (M>N,I>J)|WMNIJ> = %15.10lf\n", dot); 

    dpd_buf4_init(&W, CC3_HC1, 0, 10, 10, 12, 12, 0, "HC1 Wmnij (m>n,i>j)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<Wmnij (m>n,i>j)|Wmnij> = %15.10lf\n", dot); 

    dpd_buf4_init(&W, CC3_HC1, 0, 22, 22, 22, 22, 0, "HC1 WMnIj (Mn,Ij)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMnIj (Mn,Ij)|WMnIj> = %15.10lf\n", dot); 
*/

    /** WMBIJ **/

    /*
    dpd_buf4_init(&W, CC3_HC1, 0, 20, 0, 20, 2, 0, "HC1 WMBIJ (MB,I>J)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMBIJ (MB,I>J)|WMBIJ> = %15.10lf\n", dot);
    dpd_buf4_init(&W, CC3_HC1, 0, 30, 10, 30, 12, 0, "HC1 Wmbij (mb,i>j)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<Wmbij (mb,i>j)|Wmbij> = %15.10lf\n", dot);
    dpd_buf4_init(&W, CC3_HC1, 0, 24, 22, 24, 22, 0, "HC1 WMbIj (Mb,Ij)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMbIj (Mb,Ij)|WMbIj> = %15.10lf\n", dot);
    dpd_buf4_init(&W, CC3_HC1, 0, 27, 23, 27, 23, 0, "HC1 WmBiJ (mB,iJ)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WmBiJ (mB,iJ)|WmBiJ> = %15.10lf\n", dot);
    */

/*
    dpd_buf4_init(&W, CC3_HC1, 0, 0, 20, 2, 20, 0, "HC1 WMBIJ (I>J,MB)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMBIJ (I>J,MB)|WMBIJ> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 10, 30, 12, 30, 0, "HC1 Wmbij (i>j,mb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<Wmbij (i>j,mb)|Wmbij> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 22, 24, 22, 24, 0, "HC1 WMbIj (Ij,Mb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMbIj (Ij,Mb)|WMbIj> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 23, 27, 23, 27, 0, "HC1 WmBiJ (iJ,mB)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WmBiJ (iJ,mB)|WmBiJ> = %15.10lf\n", dot);
*/

    /** WMBEJ **/
/*
    dpd_buf4_init(&W, CC3_HC1, 0, 20, 20, 20, 20, 0, "HC1 WMBEJ (ME,JB)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMBEJ (ME,JB)|WMBEJ> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 30, 30, 30, 30, 0, "HC1 Wmbej (me,jb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<Wmbej (me,jb)|Wmbej> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 24, 26, 24, 26, 0, "HC1 WMbEj (ME,jb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMbEj (ME,jb)|WMbEj> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 27, 25, 27, 25, 0, "HC1 WmBeJ (me,JB)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WmBeJ (me,JB)|WmBeJ> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 24, 24, 24, 24, 0, "HC1 WMbeJ (Me,Jb)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WMbeJ (Me,Jb)|WMbeJ> = %15.10lf\n", dot);

    dpd_buf4_init(&W, CC3_HC1, 0, 27, 27, 27, 27, 0, "HC1 WmBEj (mE,jB)");
    dot = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("\t<WmBEj (mE,jB)|WmBEj> = %15.10lf\n", dot);
*/

    /** WABEI **/

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 20, 5, 20, 7, 0, "HC1 WABEI (IE,B>A)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WABEI (IE,A>B)|WABEI> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 30, 15, 30, 17, 0, "HC1 Wabei (ie,b>a)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<Wabei (ie,a>b)|Wabei> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 27, 29, 27, 29, 0, "HC1 WAbEi (iE,bA)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WAbEi (iE,Ab)|WAbEi> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, 0, 24, 28, 24, 28, 0, "HC1 WaBeI (Ie,Ba)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("\t<WaBeI (Ie,aB)|WaBeI> = %15.10lf\n", dot);
  }

  return;
}


}} // namespace psi::cceom
