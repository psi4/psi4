#ifndef _psi_src_bin_cints_DFT_functional_u_h
#define _psi_src_bin_cints_DFT_functional_u_h

/*! \file functional_u.h
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include "data_structs.h"

namespace psi { namespace CINTS {
  struct fun_info_s slater_u_e(struct den_info_s den_info);
  struct fun_info_s slater_u_ed(struct den_info_s den_info);
  struct fun_info_s no_funct_u(struct den_info_s den_info);
  struct fun_info_s VWN5_u_e(struct den_info_s den_info);
  struct fun_info_s VWN5_u_ed(struct den_info_s den_info);
  struct fun_info_s VWN4_u_e(struct den_info_s den_info);
  struct fun_info_s VWN4_u_ed(struct den_info_s den_info);
};}
#endif
