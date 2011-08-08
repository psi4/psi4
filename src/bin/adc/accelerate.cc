
/*
 *  accelerate.cc
 *  
 *
 *  Created by M.Saitow on 10/08/22.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {namespace adc {
						
    double accelerate(int iroot, int irrep)
    {
      char lbl[32];
      double value=1;
      dpdfile2 S, B;
      dpdbuf4 D, V, W, X, Y;
      
      timer_on("Extrapol.");
      sprintf(lbl, "V^(%d)[%d]", iroot, irrep);
      dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);
      
      sprintf(lbl, "Xjiab[%d]", irrep);
      dpd_buf4_init(&X, CC_TMP0,  irrep, 0, 5, 0, 5, 0, lbl);
      dpd_buf4_init(&V, CC_FINTS, 0    ,10, 5,10, 5, 0, "F <ia|bc>");
      dpd_contract424(&V, &B, &X, 1, 1, 1, 1, 0);
      dpd_buf4_close(&V);	
      sprintf(lbl, "Vijab[%d]", irrep);
      dpd_buf4_scmcopy(&X, CC_TMP1, lbl, 1.0);
      dpd_buf4_close(&X);
      
      sprintf(lbl, "Yijab[%d]", irrep);
      dpd_buf4_init(&Y, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
      dpd_buf4_init(&V, CC_EINTS,0    ,11, 0,11, 0, 0, "E <ai|jk>");
      dpd_contract424(&V, &B, &Y, 1, 0, 0, 1, 0);
      dpd_buf4_close(&V);
      sprintf(lbl, "Vijab[%d]", irrep);
      dpd_buf4_init(&V, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
      dpd_buf4_axpy(&Y, &V, -1.0);
      dpd_buf4_close(&Y);
      
      sprintf(lbl, "dWijab[%d]", irrep);
      dpd_buf4_scmcopy(&V, CC_TMP1, lbl, 2.0);
      sprintf(lbl, "Vjiab[%d]", irrep);
      dpd_buf4_sort(&V, CC_TMP1, qprs, 0, 5, lbl);
      sprintf(lbl, "Vijba[%d]", irrep);
      dpd_buf4_sort(&V, CC_TMP1, pqsr, 0, 5, lbl);
      dpd_buf4_close(&V);
      sprintf(lbl, "dWijab[%d]", irrep);
      dpd_buf4_init(&W, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
      sprintf(lbl, "Vjiab[%d]", irrep);
      dpd_buf4_init(&V, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
      sprintf(lbl, "Vjiba[%d]", irrep);
      dpd_buf4_sort(&V, CC_TMP1, pqsr, 0, 5, lbl);
      dpd_buf4_axpy(&V, &W, -1.0);
      dpd_buf4_close(&V);
      sprintf(lbl, "Vijba[%d]", irrep);
      dpd_buf4_init(&V, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
      dpd_buf4_axpy(&V, &W, -1.0);
      dpd_buf4_close(&V);
      sprintf(lbl, "Vjiba[%d]", irrep);
      dpd_buf4_init(&V, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
      dpd_buf4_axpy(&V, &W,  2.0);
      dpd_buf4_close(&V);
      
      sprintf(lbl, "Dijab[%d]",  irrep);
      dpd_buf4_init(&D,  CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
      dpd_buf4_dirprd(&D, &W);
      dpd_buf4_dirprd(&D, &W);
      dpd_buf4_close(&D);
      
      sprintf(lbl, "dS^(%d)[%d]", iroot, irrep);
      dpd_file2_init(&S, CC_OEI, irrep, 0, 1, lbl);
      
      dpd_buf4_init(&V, CC_FINTS,  0, 11, 5, 11, 5, 0, "F <ai|bc>");
      dpd_contract442(&W, &V, &S, 0, 0, 1, 0); 
      dpd_buf4_close(&V);
      dpd_buf4_init(&V, CC_EINTS,  0, 0, 10, 0, 10, 0, "E <ij|ka>");
      dpd_contract442(&V, &W, &S, 2, 2, -1, 1); 
      dpd_buf4_close(&V);
      
      value += dpd_file2_dot(&S, &B);
      
      dpd_buf4_close(&W);
      dpd_file2_close(&S);
      dpd_file2_close(&B);
      
      timer_off("Extrapol.");
      
      return value;
    }

}}


