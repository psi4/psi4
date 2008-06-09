#ifndef _psi_src_bin_cints_DFT_free_grid_structs_h
#define _psi_src_bin_cints_DFT_free_grid_structs_h

/*! \file free_grid_structs.h
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include "data_structs.h"

namespace psi { namespace CINTS {

  void cleanup_sphere(leb_sphere_t sphere);
  void cleanup_prim_chunk(prim_leb_chunk_t chunk);
  void cleanup_prim_atomic_grid(prim_atomic_grid_t prim_atomic_grid);
  void cleanup_prim_atomic_grid_array(prim_atomic_grid_t *prim_array, int n);
  void cleanup_conc_chunk(struct leb_chunk_s chunk);
  void cleanup_conc_atomic_grid(struct atomic_grid_s atomic_grid);
  void cleanup_grid_type(grid_t grid);
};}
#endif
