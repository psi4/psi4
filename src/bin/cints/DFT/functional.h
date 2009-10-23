#ifndef _psi_src_bin_cints_DFT_functional_h
#define _psi_src_bin_cints_DFT_functional_h

/*! \file functional.h
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include "data_structs.h"

namespace psi { namespace cints {
  struct fun_info_s slater_e(struct den_info_s den_info);
  struct fun_info_s slater_ed(struct den_info_s den_info);
  struct fun_info_s no_funct(struct den_info_s den_info);
  struct fun_info_s VWN5_e(struct den_info_s den_info);
  struct fun_info_s VWN5_ed(struct den_info_s den_info);
  struct fun_info_s VWN4_e(struct den_info_s den_info);
  struct fun_info_s VWN4_ed(struct den_info_s den_info);
  struct fun_info_s Becke88_e(struct den_info_s den_info);
  struct fun_info_s Becke88_ed(struct den_info_s den_info);
}}
#endif
