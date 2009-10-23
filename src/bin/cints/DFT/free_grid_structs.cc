/*! \file free_grid_structs.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<libciomr/libciomr.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

/* ---- Cleanup the primitive classes ---- */

namespace psi { namespace cints {
  void cleanup_sphere(leb_sphere_t sphere){
    
    free(sphere.points);
  }

  void cleanup_prim_chunk(prim_leb_chunk_t chunk){
    
    int i;
    
    for(i=0;i<chunk.size;i++)
	cleanup_sphere(chunk.spheres[i]);
  }

  void cleanup_prim_atomic_grid(prim_atomic_grid_t prim_atomic_grid){
    int i;
    
    for(i=0;i<prim_atomic_grid.chunk_num;i++)
	cleanup_prim_chunk(prim_atomic_grid.leb_chunk[i]);
  }
  
  void cleanup_prim_atomic_grid_array(prim_atomic_grid_t *prim_array, int n){
    
    int i;
    
    for(i=0;i<n;i++)
	cleanup_prim_atomic_grid(prim_array[i]);
    
    free(prim_array);
  }

  /* ---- Cleanup the concrete classes ---- */

  /*void cleanup_conc_chunk(struct leb_chunk_s *chunk){
    
  int i;
    
    for(i=0;i<chunk.size;i++)
	cleanup_sphere(chunk.spheres[i]);
	}*/

  void cleanup_conc_atomic_grid(struct atomic_grid_s atomic_grid){
    
    int i;
    
    free(atomic_grid.leb_chunk);
    /*for(i=0;i<atomic_grid.chunk_num;i++)
	cleanup_conc_chunk(atomic_grid.leb_chunk[i]);*/
  }

  void cleanup_grid_type(grid_t grid){
    
    int i;
    
    if(grid.pruned_flag)
	cleanup_prim_atomic_grid_array(grid.prim_pruned_atomic_grids,
				       grid.pruned_info.n_param_sets);
    else
	cleanup_prim_atomic_grid(grid.prim_atomic_grid);
    
    for(i=0;i<Symmetry.num_unique_atoms;i++)
	cleanup_conc_atomic_grid(grid.atomic_grid[i]);
  }
}}
