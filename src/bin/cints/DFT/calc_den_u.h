#ifndef _psi_src_bin_cints_DFT_calc_den_u_h
#define _psi_src_bin_cints_DFT_calc_den_u_h

/*! \file calc_den_u.h
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include"data_structs.h"
namespace psi { namespace CINTS {
  
  struct den_info_s calc_density_u(struct coordinates geom, int atom_num);
  
}}
#endif
