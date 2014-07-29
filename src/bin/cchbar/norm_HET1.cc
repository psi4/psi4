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
    \ingroup CCHBAR
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

void norm_HET1(void) {
  int i;
  double dot;
  dpdfile2 F;
  dpdbuf4 W;

  psi::fprintf(outfile,"Dots of (HeT1)c in names \"CC3 Wxxx\" in CC3_HET1 \n");

  if (params.ref == 0) { /* RHF */
    /** Wamef **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 11, 5, 11, 5, 0, "CC3 WAmEf (Am,Ef)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"<WAmEf (Am,Ef) | WAmEf (Am,Ef)> = %15.10lf\n", dot);

    /** Wmnie **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"<WMnIe (Mn,Ie) | WMnIe (Mn,Ie)> = %15.10lf\n", dot);

    /** Wmnij **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"<WMnIj (Mn,Ij) | WMnIj (Mn,Ij)> = %15.10lf\n", dot);

    /** Wmbij **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"<WMbIj (Ij,Mb) | WMbIj (Ij,Mb)> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WMbIj (Mb,Ij)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"<WMbIj (Mb,Ij) | WMbIj (Mb,Ij)> = %15.10lf\n", dot);

    /** Wmbej **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbEj (ME,jb)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"<WMbEj (ME,jb) | WMbEj (ME,jb)> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbeJ (Me,Jb)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"<WMbeJ (Me,Jb) | WMbeJ (Me,Jb)> = %15.10lf\n", dot);

    /** WABEI **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"<WAbEi (Ie,bA) | WAbEi (Ie,bA)> = %15.10lf\n", dot);

  }

  else if (params.ref == 1) { /* ROHF */

    psi::fprintf(outfile,"Wamef terms\n");
    /** WAMEF **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 11, 5, 11, 7, 0, "CC3 WAMEF (AM,E>F)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WAMEF|WAMEF> = %15.10lf\n", dot); 

    /** Wamef **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 11, 5, 11, 7, 0, "CC3 Wamef (am,e>f)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wamef|Wamef> = %15.10lf\n", dot); 

    /** WAmEf **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 11, 5, 11, 5, 0, "CC3 WAmEf (Am,Ef)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WAmEf|WAmEf> = %15.10lf\n", dot); 

    /** WaMeF **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 11, 5, 11, 5, 0, "CC3 WaMeF (aM,eF)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WaMeF|WaMeF> = %15.10lf\n", dot); 

    /** WMNIE **/
    psi::fprintf(outfile,"Wmnie terms\n");

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 10, 2, 10, 0, "CC3 WMNIE (M>N,IE)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMNIE|WMNIE> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 10, 2, 10, 0, "CC3 Wmnie (m>n,ie)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wmnie|Wmnie> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMnIe|WMnIe> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WmNiE (mN,iE)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WmNiE|WmNiE> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 11, 2, 11, 0, "CC3 WMNIE (M>N,EI)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMNIE(M>N,EI)|WMNIE(M>N,EI)> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 11, 2, 11, 0, "CC3 Wmnie (m>n,ei)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wmnie(m>n,ei)|Wmnie(m>n,ei)> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 11, 0, 11, 0, "CC3 WMnIe (Mn,eI)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMnIe(Mn,eI)|WMnIe(Mn,eI)> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 11, 0, 11, 0, "CC3 WmNiE (mN,Ei)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WmNiE(mN,Ei)|WmNiE(mN,Ei)> = %15.10lf\n", dot); 

    /** WMNIJ **/
    psi::fprintf(outfile,"Doing Wmnij terms.\n");

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 0, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMNIJ (M>N,IJ)|WMNIJ> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 0, 2, 2, 0, "CC3 Wmnij (m>n,i>j)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wmnij (m>n,ij)|Wmnij> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMnIj (Mn,Ij)|WMnIj> = %15.10lf\n", dot); 

    psi::fprintf(outfile,"Doing Wmbij terms.\n");
    /** WMBIJ **/

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 0, 10, 2, 0, "CC3 WMBIJ (MB,I>J)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMBIJ (MB,I>J)|WMBIJ> = %15.10lf\n", dot);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 10, 2, 10, 0, "CC3 WMBIJ (I>J,MB)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMBIJ (I>J,MB)|WMBIJ> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 0, 10, 2, 0, "CC3 Wmbij (mb,i>j)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wmbij (mb,i>j)|Wmbij> = %15.10lf\n", dot);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 10, 2, 10, 0, "CC3 Wmbij (i>j,mb)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wmbij (i>j,mb)|Wmbij> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WMbIj (Mb,Ij)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMbIj (Mb,Ij)|WMbIj> = %15.10lf\n", dot);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WMbIj (Ij,Mb)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMbIj (Ij,Mb)|WMbIj> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WmBiJ (mB,iJ)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WmBiJ (mB,iJ)|WmBiJ> = %15.10lf\n", dot);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WmBiJ (iJ,mB)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WmBiJ (iJ,mB)|WmBiJ> = %15.10lf\n", dot);

    /** WMBEJ **/
    psi::fprintf(outfile,"Doing Wmbej terms.\n");

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMBEJ (ME,JB)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMBEJ (ME,JB)|WMBEJ> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 Wmbej (me,jb)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wmbej (me,jb)|Wmbej> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbEj (ME,jb)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMbEj (ME,jb)|WMbEj> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WmBeJ (me,JB)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WmBeJ (me,JB)|WmBeJ> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbeJ (Me,Jb)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMbeJ (Me,Jb)|WMbeJ> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WmBEj (mE,jB)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WmBEj (mE,jB)|WmBEj> = %15.10lf\n", dot);

    /** WABEI **/

    psi::fprintf(outfile,"Doing Wabei terms.\n");

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 5, 10, 7, 0, "CC3 WABEI (IE,B>A)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WABEI (IE,B>A) |WABEI> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 5, 10, 7, 0, "CC3 Wabei (ie,b>a)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wabei (ie,b>a) |Wabei> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WAbEi (iE,bA) |WAbEi> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WaBeI (Ie,Ba)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WaBeI (Ie,Ba) |WaBeI> = %15.10lf\n", dot);
  }





  else { /******************** UHF */

    psi::fprintf(outfile,"Wamef terms\n");
    /** WAMEF **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 21, 5, 21, 7, 0, "CC3 WAMEF (AM,E>F)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WAMEF|WAMEF> = %15.10lf\n", dot); 

    /** Wamef **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 31, 15, 31, 17, 0, "CC3 Wamef (am,e>f)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wamef|Wamef> = %15.10lf\n", dot); 

    /** WAmEf **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 26, 28, 26, 28, 0, "CC3 WAmEf (Am,Ef)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WAmEf|WAmEf> = %15.10lf\n", dot); 

    /** WaMeF **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 25, 29, 25, 29, 0, "CC3 WaMeF (aM,eF)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WaMeF|WaMeF> = %15.10lf\n", dot); 

    /** WMNIE **/
    psi::fprintf(outfile,"Wmnie terms\n");

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMNIE (M>N,IE)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMNIE (M>N,IE)|WMNIE> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmnie (m>n,ie)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wmnie (m>n,ie)|Wmnie> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMnIe (Mn,Ie)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMnIe (Mn,Ie)|WMnIe> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmNiE (mN,iE)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WmNiE (mN,iE)|WmNiE> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 21, 2, 21, 0, "CC3 WMNIE (M>N,EI)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMNIE(M>N,EI)|WMNIE(M>N,EI)> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 31, 12, 31, 0, "CC3 Wmnie (m>n,ei)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wmnie(m>n,ei)|Wmnie(m>n,ei)> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 22, 25, 22, 25, 0, "CC3 WMnIe (Mn,eI)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMnIe(Mn,eI)|WMnIe(Mn,eI)> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 23, 26, 23, 26, 0, "CC3 WmNiE (mN,Ei)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WmNiE(mN,Ei)|WmNiE(mN,Ei)> = %15.10lf\n", dot); 

    /** WMNIJ **/
    psi::fprintf(outfile,"Doing Wmnij terms.\n");

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 0, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMNIJ (M>N,I>J)|WMNIJ> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 12, 12, 0, "CC3 Wmnij (m>n,i>j)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wmnij (m>n,i>j)|Wmnij> = %15.10lf\n", dot); 

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 22, 22, 22, 22, 0, "CC3 WMnIj (Mn,Ij)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMnIj (Mn,Ij)|WMnIj> = %15.10lf\n", dot); 

    psi::fprintf(outfile,"Doing Wmbij terms.\n");
    /** WMBIJ **/

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 20, 0, 20, 2, 0, "CC3 WMBIJ (MB,I>J)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMBIJ (MB,I>J)|WMBIJ> = %15.10lf\n", dot);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMBIJ (I>J,MB)|WMBIJ> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 30, 10, 30, 12, 0, "CC3 Wmbij (mb,i>j)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wmbij (mb,i>j)|Wmbij> = %15.10lf\n", dot);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wmbij (i>j,mb)|Wmbij> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 24, 22, 24, 22, 0, "CC3 WMbIj (Mb,Ij)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMbIj (Mb,Ij)|WMbIj> = %15.10lf\n", dot);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMbIj (Ij,Mb)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMbIj (Ij,Mb)|WMbIj> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 27, 23, 27, 23, 0, "CC3 WmBiJ (mB,iJ)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WmBiJ (mB,iJ)|WmBiJ> = %15.10lf\n", dot);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmBiJ (iJ,mB)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WmBiJ (iJ,mB)|WmBiJ> = %15.10lf\n", dot);

    /** WMBEJ **/
    psi::fprintf(outfile,"Doing Wmbej terms.\n");

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 20, 20, 20, 20, 0, "CC3 WMBEJ (ME,JB)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMBEJ (all ME,JB)|WMBEJ> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 30, 30, 30, 30, 0, "CC3 Wmbej (me,jb)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wmbej|Wmbej> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 24, 26, 24, 26, 0, "CC3 WMbEj (ME,jb)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMbEj|WMbEj> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 27, 25, 27, 25, 0, "CC3 WmBeJ (me,JB)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WmBeJ|WmBeJ> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 24, 24, 24, 24, 0, "CC3 WMbeJ (Me,Jb)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WMbeJ|WMbeJ> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 27, 27, 27, 27, 0, "CC3 WmBEj (mE,jB)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WmBEj|WmBEj> = %15.10lf\n", dot);

    /** WABEI **/

    psi::fprintf(outfile,"Doing Wabei terms.\n");

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WABEI (IE,B>A) |WABEI> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<Wabei (ie,b>a)|Wabei> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAbEi (iE,bA)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WAbEi (iE,bA)|WAbEi> = %15.10lf\n", dot);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaBeI (Ie,Ba)");
    dot = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    psi::fprintf(outfile,"\t<WAbEi (iE,Ba)|WAbEi> = %15.10lf\n", dot);
  }

  return;
}

}} // namespace psi::cchbar
