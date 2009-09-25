/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void td_setup(struct TD_Params S)
{
  dpdfile2 L1,R1;
  dpdbuf4 L2,R2;

  if((params.ref == 0) || (params.ref == 1)) {
    dpd_file2_init(&L1, CC_LAMPS, 0, 0, 1, "LIA 0 -1");
    dpd_file2_copy(&L1, CC_GLG, "LIA");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_LAMPS, 0, 0, 1, "Lia 0 -1");
    dpd_file2_copy(&L1, CC_GLG, "Lia");
    dpd_file2_close(&L1);

    dpd_buf4_init(&L2, CC_LAMPS, 0, 2, 7, 2, 7, 0, "LIJAB 0 -1");
    dpd_buf4_copy(&L2, CC_GLG, "LIJAB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, 0, 2, 7, 2, 7, 0, "Lijab 0 -1");
    dpd_buf4_copy(&L2, CC_GLG, "Lijab");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb 0 -1");
    dpd_buf4_copy(&L2, CC_GLG, "LIjAb");
    dpd_buf4_close(&L2);
  }
  else if (params.ref == 2) {
    dpd_file2_init(&L1, CC_LAMPS, 0, 0, 1, "LIA 0 -1");
    dpd_file2_copy(&L1, CC_GLG, "LIA");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_LAMPS, 0, 2, 3, "Lia 0 -1");
    dpd_file2_copy(&L1, CC_GLG, "Lia");
    dpd_file2_close(&L1);

    dpd_buf4_init(&L2, CC_LAMPS, 0, 2, 7, 2, 7, 0, "LIJAB 0 -1");
    dpd_buf4_copy(&L2, CC_GLG, "LIJAB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, 0, 12, 17, 12, 17, 0, "Lijab 0 -1");
    dpd_buf4_copy(&L2, CC_GLG, "Lijab");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, 0, 22, 28, 22, 28, 0, "LIjAb 0 -1");
    dpd_buf4_copy(&L2, CC_GLG, "LIjAb");
    dpd_buf4_close(&L2);
  }
 
  if((params.ref==0) || (params.ref==1)) {
    dpd_file2_init(&L1, CC_LAMPS, S.irrep, 0, 1, S.L1A_lbl);
    dpd_file2_copy(&L1, CC_GL, "LIA");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_LAMPS, S.irrep, 0, 1, S.L1B_lbl);
    dpd_file2_copy(&L1, CC_GL, "Lia");
    dpd_file2_close(&L1);

    dpd_buf4_init(&L2, CC_LAMPS, S.irrep, 2, 7, 2, 7, 0, S.L2AA_lbl);
    dpd_buf4_copy(&L2, CC_GL, "LIJAB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, S.irrep, 2, 7, 2, 7, 0, S.L2BB_lbl);
    dpd_buf4_copy(&L2, CC_GL, "Lijab");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, S.irrep, 0, 5, 0, 5, 0, S.L2AB_lbl);
    dpd_buf4_copy(&L2, CC_GL, "LIjAb");
    dpd_buf4_close(&L2);
  }
  else if (params.ref == 2) {
    dpd_file2_init(&L1, CC_LAMPS, S.irrep, 0, 1, S.L1A_lbl);
    dpd_file2_copy(&L1, CC_GL, "LIA");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_LAMPS, S.irrep, 2, 3, S.L1B_lbl);
    dpd_file2_copy(&L1, CC_GL, "Lia");
    dpd_file2_close(&L1);

    dpd_buf4_init(&L2, CC_LAMPS, S.irrep, 2, 7, 2, 7, 0, S.L2AA_lbl);
    dpd_buf4_copy(&L2, CC_GL, "LIJAB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, S.irrep, 12, 17, 12, 17, 0, S.L2BB_lbl);
    dpd_buf4_copy(&L2, CC_GL, "Lijab");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, S.irrep, 22, 28, 22, 28, 0, S.L2AB_lbl);
    dpd_buf4_copy(&L2, CC_GL, "LIjAb");
    dpd_buf4_close(&L2);
  }

  if((params.ref == 0) || (params.ref == 1)) {
    dpd_file2_init(&R1, CC_RAMPS, S.irrep, 0, 1, S.R1A_lbl);
    dpd_file2_copy(&R1, CC_GR, "RIA");
    dpd_file2_close(&R1);

    dpd_file2_init(&R1, CC_RAMPS, S.irrep, 0, 1, S.R1B_lbl);
    dpd_file2_copy(&R1, CC_GR, "Ria");
    dpd_file2_close(&R1);

    dpd_buf4_init(&R2, CC_RAMPS, S.irrep, 2, 7, 2, 7, 0, S.R2AA_lbl);
    dpd_buf4_copy(&R2, CC_GR, "RIJAB");
    dpd_buf4_close(&R2);

    dpd_buf4_init(&R2, CC_RAMPS, S.irrep, 2, 7, 2, 7, 0, S.R2BB_lbl);
    dpd_buf4_copy(&R2, CC_GR, "Rijab");
    dpd_buf4_close(&R2);

    dpd_buf4_init(&R2, CC_RAMPS, S.irrep, 0, 5, 0, 5, 0, S.R2AB_lbl);
    dpd_buf4_copy(&R2, CC_GR, "RIjAb");
    dpd_buf4_close(&R2);
  }
  else if (params.ref == 2) {
    dpd_file2_init(&R1, CC_RAMPS, S.irrep, 0, 1, S.R1A_lbl);
    dpd_file2_copy(&R1, CC_GR, "RIA");
    dpd_file2_close(&R1);

    dpd_file2_init(&R1, CC_RAMPS, S.irrep, 2, 3, S.R1B_lbl);
    dpd_file2_copy(&R1, CC_GR, "Ria");
    dpd_file2_close(&R1);

    dpd_buf4_init(&R2, CC_RAMPS, S.irrep, 2, 7, 2, 7, 0, S.R2AA_lbl);
    dpd_buf4_copy(&R2, CC_GR, "RIJAB");
    dpd_buf4_close(&R2);

    dpd_buf4_init(&R2, CC_RAMPS, S.irrep, 12, 17, 12, 17, 0, S.R2BB_lbl);
    dpd_buf4_copy(&R2, CC_GR, "Rijab");
    dpd_buf4_close(&R2);

    dpd_buf4_init(&R2, CC_RAMPS, S.irrep, 22, 28, 22, 28, 0, S.R2AB_lbl);
    dpd_buf4_copy(&R2, CC_GR, "RIjAb");
    dpd_buf4_close(&R2);
  }

  if((params.ref==0) || (params.ref==1)) {
    if((S.irrep == 0)) {
      dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
      dpd_file2_scm(&L1, S.R0);
      dpd_file2_close(&L1);

      dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
      dpd_file2_scm(&L1, S.R0);
      dpd_file2_close(&L1);

      dpd_buf4_init(&L2, CC_GLG, 0, 2, 7, 2, 7, 0, "LIJAB");
      dpd_buf4_scm(&L2, S.R0);
      dpd_buf4_close(&L2);

      dpd_buf4_init(&L2, CC_GLG, 0, 2, 7, 2, 7, 0, "Lijab");
      dpd_buf4_scm(&L2, S.R0);
      dpd_buf4_close(&L2);

      dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
      dpd_buf4_scm(&L2, S.R0);
      dpd_buf4_close(&L2);
    }

    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_GLG, qpsr, 0, 5, "LiJaB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_sort(&L2, CC_GLG, prqs, 10, 10, "LIAJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 2, 7, 0, "Lijab");
    dpd_buf4_sort(&L2, CC_GLG, prqs, 10, 10, "Liajb");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_GLG, prqs, 10, 10, "LIAjb");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_sort(&L2, CC_GLG, prqs, 10, 10, "LiaJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_GLG, 0, 10, 10, 10, 10, 0, "LIAjb");
    dpd_buf4_sort(&L2, CC_GLG, psrq, 10, 10, "LIbjA");
    dpd_buf4_sort(&L2, CC_GLG, rqps, 10, 10, "LjAIb");
    dpd_buf4_close(&L2);
  }
  else if(params.ref == 2) { 
    if((S.irrep == 0)) {
      dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
      dpd_file2_scm(&L1, S.R0);
      dpd_file2_close(&L1);

      dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
      dpd_file2_scm(&L1, S.R0);
      dpd_file2_close(&L1);

      dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 2, 7, 0, "LIJAB");
      dpd_buf4_scm(&L2, S.R0);
      dpd_buf4_close(&L2);

      dpd_buf4_init(&L2, CC_GLG, 0, 10, 15, 12, 17, 0, "Lijab");
      dpd_buf4_scm(&L2, S.R0);
      dpd_buf4_close(&L2);

      dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
      dpd_buf4_scm(&L2, S.R0);
      dpd_buf4_close(&L2);
    }

    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_GLG, qpsr, 23, 29, "LiJaB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_sort(&L2, CC_GLG, prqs, 20, 20, "LIAJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_GLG, 0, 10, 15, 12, 17, 0, "Lijab");
    dpd_buf4_sort(&L2, CC_GLG, prqs, 30, 30, "Liajb");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_GLG, prqs, 20, 30, "LIAjb");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_sort(&L2, CC_GLG, prqs, 30, 20, "LiaJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_GLG, 0, 20, 30, 20, 30, 0, "LIAjb");
    dpd_buf4_sort(&L2, CC_GLG, psrq, 24, 27, "LIbjA");
    dpd_buf4_sort(&L2, CC_GLG, rqps, 27, 24, "LjAIb");
    dpd_buf4_close(&L2);
  }

  if((params.ref==0) || (params.ref==1)) { /** RHF/ROHF **/
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_GL, qpsr, 0, 5, "LiJaB");
    dpd_buf4_sort(&L2, CC_GL, pqsr, 0, 5, "LIjaB");
    dpd_buf4_sort(&L2, CC_GL, qprs, 0, 5, "LiJAb");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_sort(&L2, CC_GL, prqs, 10, 10, "LIAJB");
    dpd_buf4_close(&L2);
    
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 5, 2, 7, 0, "Lijab");
    dpd_buf4_sort(&L2, CC_GL, prqs, 10, 10, "Liajb");
    dpd_buf4_close(&L2);
    
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_GL, prqs, 10, 10, "LIAjb");
    dpd_buf4_close(&L2);
    
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_sort(&L2, CC_GL, prqs, 10, 10, "LiaJB");
    dpd_buf4_close(&L2);
    
    dpd_buf4_init(&L2, CC_GL, S.irrep, 10, 10, 10, 10, 0, "LIAjb");
    dpd_buf4_sort(&L2, CC_GL, psrq, 10, 10, "LIbjA");
    dpd_buf4_sort(&L2, CC_GL, rqps, 10, 10, "LjAIb");
    dpd_buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_buf4_init(&L2, CC_GL, S.irrep, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_GL, qpsr, 23, 29, "LiJaB");
    dpd_buf4_sort(&L2, CC_GL, pqsr, 22, 29, "LIjaB");
    dpd_buf4_sort(&L2, CC_GL, qprs, 23, 28, "LiJAb");
    dpd_buf4_close(&L2);
  
    dpd_buf4_init(&L2, CC_GL, S.irrep, 0, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_sort(&L2, CC_GL, prqs, 20, 20, "LIAJB");
    dpd_buf4_close(&L2);
  
    dpd_buf4_init(&L2, CC_GL, S.irrep, 10, 15, 12, 17, 0, "Lijab");
    dpd_buf4_sort(&L2, CC_GL, prqs, 30, 30, "Liajb");
    dpd_buf4_close(&L2);
  
    dpd_buf4_init(&L2, CC_GL, S.irrep, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_GL, prqs, 20, 30, "LIAjb");
    dpd_buf4_close(&L2);
  
    dpd_buf4_init(&L2, CC_GL, S.irrep, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_sort(&L2, CC_GL, prqs, 30, 20, "LiaJB");
    dpd_buf4_close(&L2);
  
    dpd_buf4_init(&L2, CC_GL, S.irrep, 20, 30, 20, 30, 0, "LIAjb");
    dpd_buf4_sort(&L2, CC_GL, psrq, 24, 27, "LIbjA");
    dpd_buf4_sort(&L2, CC_GL, rqps, 27, 24, "LjAIb");
    dpd_buf4_close(&L2);
  }

  if((params.ref==0) || (params.ref==1)) { /** RHF/ROHF **/
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 0, 5, 0, "RIjAb");
    dpd_buf4_sort(&R2, CC_GR, qpsr, 0, 5, "RiJaB");
    dpd_buf4_sort(&R2, CC_GR, pqsr, 0, 5, "RIjaB");
    dpd_buf4_sort(&R2, CC_GR, qprs, 0, 5, "RiJAb");
    dpd_buf4_close(&R2);
    
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 2, 7, 0, "RIJAB");
    dpd_buf4_sort(&R2, CC_GR, prqs, 10, 10, "RIAJB");
    dpd_buf4_close(&R2);
    
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 2, 7, 0, "Rijab");
    dpd_buf4_sort(&R2, CC_GR, prqs, 10, 10, "Riajb");
    dpd_buf4_close(&R2);
    
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 0, 5, 0, "RIjAb");
    dpd_buf4_sort(&R2, CC_GR, prqs, 10, 10, "RIAjb");
    dpd_buf4_close(&R2);
    
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 0, 5, 0, "RiJaB");
    dpd_buf4_sort(&R2, CC_GR, prqs, 10, 10, "RiaJB");
    dpd_buf4_close(&R2);
    
    dpd_buf4_init(&R2, CC_GR, S.irrep, 10, 10, 10, 10, 0, "RIAjb");
    dpd_buf4_sort(&R2, CC_GR, psrq, 10, 10, "RIbjA");
    dpd_buf4_sort(&R2, CC_GR, rqps, 10, 10, "RjAIb");
    dpd_buf4_close(&R2);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_buf4_init(&R2, CC_GR, S.irrep, 22, 28, 22, 28, 0, "RIjAb");
    dpd_buf4_sort(&R2, CC_GR, qpsr, 23, 29, "RiJaB");
    dpd_buf4_sort(&R2, CC_GR, pqsr, 22, 29, "RIjaB");
    dpd_buf4_sort(&R2, CC_GR, qprs, 23, 28, "RiJAb");
    dpd_buf4_close(&R2);
  
    dpd_buf4_init(&R2, CC_GR, S.irrep, 0, 5, 2, 7, 0, "RIJAB");
    dpd_buf4_sort(&R2, CC_GR, prqs, 20, 20, "RIAJB");
    dpd_buf4_close(&R2);
  
    dpd_buf4_init(&R2, CC_GR, S.irrep, 10, 15, 12, 17, 0, "Rijab");
    dpd_buf4_sort(&R2, CC_GR, prqs, 30, 30, "Riajb");
    dpd_buf4_close(&R2);
  
    dpd_buf4_init(&R2, CC_GR, S.irrep, 22, 28, 22, 28, 0, "RIjAb");
    dpd_buf4_sort(&R2, CC_GR, prqs, 20, 30, "RIAjb");
    dpd_buf4_close(&R2);
  
    dpd_buf4_init(&R2, CC_GR, S.irrep, 23, 29, 23, 29, 0, "RiJaB");
    dpd_buf4_sort(&R2, CC_GR, prqs, 30, 20, "RiaJB");
    dpd_buf4_close(&R2);
  
    dpd_buf4_init(&R2, CC_GR, S.irrep, 20, 30, 20, 30, 0, "RIAjb");
    dpd_buf4_sort(&R2, CC_GR, psrq, 24, 27, "RIbjA");
    dpd_buf4_sort(&R2, CC_GR, rqps, 27, 24, "RjAIb");
    dpd_buf4_close(&R2);
  }

  return;
}

}} // namespace psi::ccdensity
