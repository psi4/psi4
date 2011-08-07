
/*
 *  orbital_relaxation.cc
 *  
 *
 *  Created by M.Saitow on 10/07/28.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <cstring>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {namespace adc {
			
    void make_screened_interaction(int iroot, int irrep, double omega)
    {
      char lbl[32];
      dpdfile2 S, B;
      dpdbuf4 D, X, Y, V, W;
      
      if(params.spin == 0){		
	sprintf(lbl, "B^(%d)[%d]", iroot, irrep);
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
	
	dpd_file2_close(&B);
	
	sprintf(lbl, "Wijab[%d]", irrep);
	dpd_buf4_scmcopy(&V, CC_TMP1, lbl, 2.0);
	sprintf(lbl, "Vjiab[%d]", irrep);
	dpd_buf4_sort(&V, CC_TMP1, qprs, 0, 5, lbl);
	sprintf(lbl, "Vijba[%d]", irrep);
	dpd_buf4_sort(&V, CC_TMP1, pqsr, 0, 5, lbl);
	dpd_buf4_close(&V);
	sprintf(lbl, "Wijab[%d]", irrep);
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
	dpd_buf4_close(&D);
	
	sprintf(lbl, "S^(%d)[%d]", iroot, irrep);
	dpd_file2_init(&S, CC_OEI, irrep, 0, 1, lbl);
	
	dpd_buf4_init(&V, CC_FINTS,  0, 11, 5, 11, 5, 0, "F <ai|bc>");
	dpd_contract442(&W, &V, &S, 0, 0, 1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_init(&V, CC_EINTS,  0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&V, &W, &S, 2, 2, -1, 1); 
	dpd_buf4_close(&V);
	
	dpd_buf4_close(&W);
	dpd_file2_close(&S);
      }
      else{
	/*
	  sprintf(lbl, "B^(%d)[%d]", iroot, irrep);
	  dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);
	  
	  sprintf(lbl, "Yijab[%d]", irrep);
	  dpd_buf4_init(&Y, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
	  dpd_buf4_init(&V, CC_EINTS,0    ,11, 0,11, 0, 0, "E <ai|jk>");
	  dpd_contract424(&V, &B, &Y, 1, 0, 0, 1, 0);
	  dpd_buf4_close(&V);
	  
	  sprintf(lbl, "Vijab[%d]", irrep);
	  dpd_buf4_scmcopy(&Y, CC_TMP1, lbl, -1.0);
	  dpd_buf4_close(&Y);
	  
	  sprintf(lbl, "Xjiab[%d]", irrep);
	  dpd_buf4_init(&X, CC_TMP0,  irrep, 0, 5, 0, 5, 0, lbl);
	  dpd_buf4_init(&V, CC_FINTS, 0    ,10, 5,10, 5, 0, "F <ia|bc>");
	  dpd_contract424(&V, &B, &X, 1, 1, 1, 1, 0);
	  dpd_buf4_close(&V);	
	  sprintf(lbl, "Vijab[%d]", irrep);
	  dpd_buf4_sort_axpy(&X, CC_TMP1, qprs, 0, 5, lbl, 1.0);
	  dpd_buf4_close(&X);
	  
	  dpd_file2_close(&B);
	  
	  sprintf(lbl, "Vijab[%d]", irrep);
	  dpd_buf4_init(&V, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
	  
	  sprintf(lbl, "W1ijab[%d]", irrep);
	  dpd_buf4_init(&W, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
	  dpd_buf4_scm(&W, 0.0);
	  dpd_buf4_close(&W);
	  dpd_buf4_sort_axpy(&V, CC_TMP1, qprs, 0, 5, lbl,  1.5);
	  dpd_buf4_sort_axpy(&V, CC_TMP1, qpsr, 0, 5, lbl, -1.5);
	  
	  sprintf(lbl, "W2ijab[%d]", irrep);
	  dpd_buf4_scmcopy(&V, CC_TMP1, lbl, -2.0);
	  dpd_buf4_sort_axpy(&V, CC_TMP1, qprs, 0, 5, lbl,  1.0);
	  dpd_buf4_sort_axpy(&V, CC_TMP1, pqsr, 0, 5, lbl,  1.0);
	  dpd_buf4_close(&V);
	  
	  sprintf(lbl, "S^(%d)[%d]", iroot, irrep);
	  dpd_file2_init(&S, CC_OEI, irrep, 0, 1, lbl);
	  
	  sprintf(lbl, "Dijab[%d]",  irrep);
	  dpd_buf4_init(&D,  CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	  sprintf(lbl, "W1ijab[%d]", irrep);
	  dpd_buf4_init(&W, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
	  dpd_buf4_dirprd(&D, &W);
	  dpd_buf4_init(&V, CC_FINTS,  0, 11, 5, 11, 5, 0, "F <ai|bc>");
	  dpd_contract442(&W, &V, &S, 0, 0, 1, 1); // +1
	  dpd_buf4_close(&V);
	  dpd_buf4_close(&W);
	  
	  sprintf(lbl, "W2ijab[%d]", irrep);
	  dpd_buf4_init(&W, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
	  dpd_buf4_dirprd(&D, &W);
	  dpd_buf4_close(&D);
	  dpd_buf4_init(&V, CC_EINTS,  0, 0, 10, 0, 10, 0, "E <ij|ka>");
	  dpd_contract442(&V, &W, &S, 2, 2, -1, 1); // -1
	  dpd_buf4_close(&V);
	  
	  dpd_buf4_close(&W);
	  dpd_file2_close(&S);
	*/
	/* These code should be deprecated! */ 
	/*	dpdbuf4 V, U1, U2, U3, Z1, Z2, Z3;
	sprintf(lbl, "B^(%d)[%d]", iroot, irrep);
	dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);
	
	sprintf(lbl, "Xjiab[%d]", irrep);
	dpd_buf4_init(&X, CC_TMP0,  irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_init(&V, CC_FINTS, 0    ,10, 5,10, 5, 0, "F <ia|bc>");
	dpd_contract424(&V, &B, &X, 1, 1, 1, 1, 0);
	dpd_buf4_close(&V);	
	sprintf(lbl, "Xijab[%d]", irrep);
	dpd_buf4_sort(&X, CC_TMP0, qprs, 0, 5, lbl);
	dpd_buf4_close(&X);
	
	sprintf(lbl, "Yijab[%d]", irrep);
	dpd_buf4_init(&Y, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_init(&V, CC_EINTS,0    ,11, 0,11, 0, 0, "E <ai|jk>");
	dpd_contract424(&V, &B, &Y, 1, 0, 0, 1, 0);
	dpd_buf4_close(&V);
	sprintf(lbl, "Yjiab[%d]", irrep);
	dpd_buf4_sort(&Y, CC_TMP1, qprs, 0, 5, lbl);
	dpd_buf4_close(&Y);
	dpd_file2_close(&B);
	
	sprintf(lbl, "Z1ijab[%d]", irrep);
	dpd_buf4_init(&Z1, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_scm(&Z1, 0.0);
	
	sprintf(lbl, "Xjiab[%d]", irrep);
	dpd_buf4_init(&X, CC_TMP0,  irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_axpy(&X, &Z1,  1); 
	dpd_buf4_close(&X);
	sprintf(lbl, "Xijab[%d]", irrep);
	dpd_buf4_init(&X, CC_TMP0,  irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_axpy(&X, &Z1, -1); 
	dpd_buf4_close(&X);
	sprintf(lbl, "Yijab[%d]", irrep);
	dpd_buf4_init(&Y, CC_TMP1,  irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_axpy(&Y, &Z1,  1); 
	dpd_buf4_close(&Y);
	sprintf(lbl, "Yjiab[%d]", irrep);
	dpd_buf4_init(&Y, CC_TMP1,  irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_axpy(&Y, &Z1, -1); 
	dpd_buf4_close(&Y);
	
	sprintf(lbl, "U1ijab[%d]", irrep);
	dpd_buf4_init(&U1, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_scm(&U1, 0.0);
	dpd_buf4_axpy(&Z1, &U1,  1);
	sprintf(lbl, "Z1ijba[%d]", irrep);
	dpd_buf4_sort(&Z1, CC_MISC, pqsr, 0, 5, lbl);
	dpd_buf4_close(&Z1);
	dpd_buf4_init(&Z1, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_axpy(&Z1, &U1,  1);
	dpd_buf4_close(&Z1);
	dpd_buf4_close(&U1);
	
	sprintf(lbl, "Z2ijab[%d]", irrep);
	dpd_buf4_init(&Z2, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_scm(&Z2, 0.0);
	
	sprintf(lbl, "Xjiab[%d]", irrep);
	dpd_buf4_init(&X, CC_TMP0,  irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_axpy(&X, &Z2,  1); 
	dpd_buf4_close(&X);
	sprintf(lbl, "Xijab[%d]", irrep);
	dpd_buf4_init(&X, CC_TMP0,  irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_axpy(&X, &Z2, -1); 
	dpd_buf4_close(&X);
	sprintf(lbl, "Yijab[%d]", irrep);
	dpd_buf4_init(&Y, CC_TMP1,  irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_axpy(&Y, &Z2, -1); 
	dpd_buf4_close(&Y);
	sprintf(lbl, "Yjiab[%d]", irrep);
	dpd_buf4_init(&Y, CC_TMP1,  irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_axpy(&Y, &Z2,  1); 
	dpd_buf4_close(&Y);
	
	sprintf(lbl, "U2ijab[%d]", irrep);
	dpd_buf4_init(&U2, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_scm(&U2, 0.0);	
	dpd_buf4_axpy(&Z2, &U2,  1);
	sprintf(lbl, "Z2ijba[%d]", irrep);
	dpd_buf4_sort(&Z2, CC_MISC, pqsr, 0, 5, lbl);
	dpd_buf4_close(&Z2);
	dpd_buf4_init(&Z2, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_axpy(&Z2, &U2, -1);
	dpd_buf4_close(&Z2);
	dpd_buf4_close(&U2);
	
	sprintf(lbl, "Dijab[%d]",  irrep);
	dpd_buf4_init(&D,  CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	sprintf(lbl, "U1ijab[%d]", irrep);
	dpd_buf4_init(&U1, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_dirprd(&D, &U1);
	dpd_buf4_close(&U1);
	sprintf(lbl, "U2ijab[%d]", irrep);
	dpd_buf4_init(&U2, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);	
	dpd_buf4_dirprd(&D, &U2);
	dpd_buf4_close(&U2);
	dpd_buf4_close(&D);
	
	sprintf(lbl, "S^(%d)[%d]", iroot, irrep);
	dpd_file2_init(&S, CC_OEI, irrep, 0, 1, lbl);
	
	sprintf(lbl, "U1ijab[%d]", irrep);
	dpd_buf4_init(&U1, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_scm(&U1, 1.0/8.0);
	sprintf(lbl, "U1jiab[%d]", irrep);
	dpd_buf4_sort(&U1, CC_MISC, qprs, 0, 5, lbl);
	sprintf(lbl, "U1ijba[%d]", irrep);
	dpd_buf4_sort(&U1, CC_MISC, pqsr, 0, 5, lbl);
	
	dpd_buf4_init(&V, CC_FINTS,  0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_contract442(&U1, &V, &S, 1, 1, -1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_init(&V, CC_FINTS,  0, 11, 5, 11, 5, 0, "F <ai|bc>");
	dpd_contract442(&U1, &V, &S, 0, 0, -1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_init(&V, CC_EINTS,  0, 0, 11, 0, 11, 0, "E <ij|ak>");
	dpd_contract442(&V, &U1, &S, 3, 3,  1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_init(&V, CC_EINTS,  0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&V, &U1, &S, 2, 2,  1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_close(&U1);
	
	sprintf(lbl, "U1jiab[%d]", irrep); 
	dpd_buf4_init(&U1, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_init(&V, CC_FINTS,  0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_contract442(&U1, &V, &S, 1, 1, -1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_init(&V, CC_FINTS,  0, 11, 5, 11, 5, 0, "F <ai|bc>");
	dpd_contract442(&U1, &V, &S, 0, 0, -1, 1);
	dpd_buf4_close(&V);
	dpd_buf4_close(&U1);
	
	sprintf(lbl, "U1ijba[%d]", irrep); 
	dpd_buf4_init(&U1, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_init(&V, CC_EINTS,  0, 0, 11, 0, 11, 0, "E <ij|ak>");
	dpd_contract442(&V, &U1, &S, 3, 3,  1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_init(&V, CC_EINTS,  0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&V, &U1, &S, 2, 2,  1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_close(&U1);
	
	sprintf(lbl, "U2ijab[%d]", irrep);
	dpd_buf4_init(&U2, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_scm(&U2, 3.0/8.0);
	sprintf(lbl, "U2jiab[%d]", irrep);
	dpd_buf4_sort(&U2, CC_MISC, qprs, 0, 5, lbl);
	sprintf(lbl, "U2ijba[%d]", irrep);
	dpd_buf4_sort(&U2, CC_MISC, pqsr, 0, 5, lbl);
	dpd_buf4_init(&V, CC_FINTS,  0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_contract442(&U2, &V, &S, 1, 1,  1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_init(&V, CC_FINTS,  0, 11, 5, 11, 5, 0, "F <ai|bc>");
	dpd_contract442(&U2, &V, &S, 0, 0,  1, 1);
	dpd_buf4_close(&V);
	dpd_buf4_init(&V, CC_EINTS,  0, 0, 11, 0, 11, 0, "E <ij|ak>");
	dpd_contract442(&V, &U2, &S, 3, 3, -1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_init(&V, CC_EINTS,  0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&V, &U2, &S, 2, 2, -1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_close(&U2);
	
	sprintf(lbl, "U2jiab[%d]", irrep);
	dpd_buf4_init(&U2, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_init(&V, CC_FINTS,  0, 10, 5, 10, 5, 0, "F <ia|bc>");
	dpd_contract442(&U2, &V, &S, 1, 1, -1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_init(&V, CC_FINTS,  0, 11, 5, 11, 5, 0, "F <ai|bc>");
	dpd_contract442(&U2, &V, &S, 0, 0, -1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_close(&U2);
	
	sprintf(lbl, "U2ijba[%d]", irrep); 
	dpd_buf4_init(&U2, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_init(&V, CC_EINTS,  0, 0, 11, 0, 11, 0, "E <ij|ak>");
	dpd_contract442(&V, &U2, &S, 3, 3,  1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_init(&V, CC_EINTS,  0, 0, 10, 0, 10, 0, "E <ij|ka>");
	dpd_contract442(&V, &U2, &S, 2, 2,  1, 1); 
	dpd_buf4_close(&V);
	dpd_buf4_close(&U2);
	
	dpd_file2_close(&S);		
	*/
      }
    }

}}


