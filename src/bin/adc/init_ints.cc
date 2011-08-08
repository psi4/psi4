
/*
 *  init_ints.cc
 *  
 *
 *  Created by M.Saitow on 10/08/05.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {namespace adc {
		
    void init_ints(void)
    {
      dpdbuf4 V;
      
      dpd_buf4_init(&V, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
      dpd_buf4_sort(&V, CC_FINTS, qpsr, 11, 5, "F <ai|bc>");
      dpd_buf4_close(&V);
      
    }
	
}}

