#ifndef _psi_src_bin_cints_DFT_init_prun_prim_atomic_grid_h
#define _psi_src_bin_cints_DFT_init_prun_prim_atomic_grid_h

/*! \file init_prun_prim_atomic_grid.h
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include"data_structs.h"
namespace psi { namespace CINTS {
  prim_atomic_grid_t *init_pruned_prim_atomic_grid(int n_rpoints, int num_chunk, struct pruned_info_s pruned_info);
}}
#endif
