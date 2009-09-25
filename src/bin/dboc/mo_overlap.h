/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/

#ifndef _psi3_bin_dboc_mooverlap_h_
#define _psi3_bin_dboc_mooverlap_h_

#include "float.h"
#include "defines.h"

namespace psi { namespace dboc {

FLOAT **eval_S_alpha(DisplacementIndex LDisp, DisplacementIndex RDisp);
FLOAT **eval_S_beta(DisplacementIndex LDisp, DisplacementIndex RDisp);

}} /* namespace psi::dboc */

#endif
