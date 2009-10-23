#ifndef _psi_src_bin_cints_DFT_calc_den_fast_h
#define _psi_src_bin_cints_DFT_calc_den_fast_h

/*! \file calc_den_fast.h
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include "data_structs.h"

namespace psi { namespace cints {
  
  struct den_info_s calc_density_fast(struct coordinates geom, int atom_num);

}}
#endif
