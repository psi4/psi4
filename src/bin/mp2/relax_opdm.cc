/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void relax_opdm(void)
{
  dpdfile2 D1, D2;

  if(params.ref == 0) {   
    dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "DAI");
    dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);

    dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "DIA");
    dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);

    dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "Dai");
    dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(a,i)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);

    dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "Dia");
    dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(a,i)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);
  }
  else if(params.ref == 2) {

    dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "DAI");
    dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);

    dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "DIA");
    dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);

    dpd_file2_init(&D1, CC_OEI, 0, 2, 3, "Dai");
    dpd_file2_init(&D2, CC_OEI, 0, 3, 2, "D(orb)(a,i)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);

    dpd_file2_init(&D1, CC_OEI, 0, 2, 3, "Dia");
    dpd_file2_init(&D2, CC_OEI, 0, 3, 2, "D(orb)(a,i)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);
  }
}

}} /* End namespaces */
