/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void build_X(void)
{
  dpdfile2 X, I, X2;

  if(params.ref == 0) { 

    dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");
    dpd_file2_copy(&I, CC_OEI, "XAI");
    dpd_file2_close(&I);

    dpd_file2_init(&X, CC_OEI, 0, 1, 0, "XAI");
    dpd_file2_scm(&X, -1.0);
    dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");
    dpd_file2_axpy(&I, &X, 1.0, 1);
    dpd_file2_close(&I);
    dpd_file2_close(&X);

    dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");
    dpd_file2_copy(&I, CC_OEI, "Xai");
    dpd_file2_close(&I);

    dpd_file2_init(&X, CC_OEI, 0, 1, 0, "Xai");
    dpd_file2_scm(&X, -1.0);
    dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");
    dpd_file2_axpy(&I, &X, 1.0, 1);
    dpd_file2_close(&I);
    dpd_file2_close(&X);

    /* Build spatial orbital version of X for Zvector:
       X(A,I) = 1/2[X(A,I)+X(a,i)]  */
    dpd_file2_init(&X, CC_OEI, 0, 1, 0, "XAI");
    dpd_file2_copy(&X, CC_MISC, "X(A,I)");
    dpd_file2_close(&X);
    dpd_file2_init(&X, CC_MISC, 0, 1, 0, "X(A,I)");
    dpd_file2_init(&X2, CC_OEI, 0, 1, 0, "Xai");
    dpd_file2_axpy(&X2, &X, 1.0, 0);
    dpd_file2_close(&X2);
    dpd_file2_scm(&X, 0.5);
    dpd_file2_close(&X);
  }
  else if(params.ref == 2) { 

    dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");
    dpd_file2_copy(&I, CC_OEI, "XAI");
    dpd_file2_close(&I);

    dpd_file2_init(&X, CC_OEI, 0, 1, 0, "XAI");
    dpd_file2_scm(&X, -1.0);
    dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");
    dpd_file2_axpy(&I, &X, 1.0, 1);
    dpd_file2_close(&I);
    dpd_file2_close(&X);

    dpd_file2_init(&I, CC_OEI, 0, 3, 2, "I'ai");
    dpd_file2_copy(&I, CC_OEI, "Xai");
    dpd_file2_close(&I);

    dpd_file2_init(&X, CC_OEI, 0, 3, 2, "Xai");
    dpd_file2_scm(&X, -1.0);
    dpd_file2_init(&I, CC_OEI, 0, 2, 3, "I'ia");
    dpd_file2_axpy(&I, &X, 1.0, 1);
    dpd_file2_close(&I);
    dpd_file2_close(&X);
  }

}

}} /* End namespaces */
