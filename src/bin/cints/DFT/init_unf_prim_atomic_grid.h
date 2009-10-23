#ifndef _psi_src_bin_cints_DFT_init_unf_prim_atomic_grid_h
#define _psi_src_bin_cints_DFT_init_unf_prim_atomic_grid_h

/*! \file init_unf_prim_atomic_grid.h
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include"data_structs.h"

namespace psi { namespace cints {
  prim_atomic_grid_t init_uniform_prim_atomic_grid(int n_rpoints,int n_angpoints,int num_chunks);
}}
#endif
