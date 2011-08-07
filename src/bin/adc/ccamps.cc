
/*
 *  ccamps.cc
 *  
 *
 *  Created by M.Saitow on 10/09/2.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <cmath>
#include <cstdlib>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace adc {

    void ccamps(void)
    {
      dpdbuf4 T2, Z;
      
      dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_buf4_sort(&T2, CC_TMP0, pqsr, 0, 5, "tIjbA");
      dpd_buf4_copy(&T2, CC_TAMPS, "2 tIjAb - tIjbA");
      dpd_buf4_close(&T2);
      dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjbA");
      dpd_buf4_scm(&T2, 2);
      dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjbA");
      dpd_buf4_axpy(&Z, &T2, -1);
      dpd_buf4_close(&Z);
      dpd_buf4_close(&T2);
    }

}}
