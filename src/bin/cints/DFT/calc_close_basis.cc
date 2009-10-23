/*! \file calc_close_basis.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

#include"defines.h"
#define EXTERN
#define NIM(a,b) a < b ? a : b
#include"global.h"
#include <stdexcept>
#include"small_fns.h"
#include"bas_comp_functions.h"
#define TOL 1E-50

namespace psi { namespace cints {

  void calc_close_basis(int atom_num, int chunk_num){
    int i,j,k,l,m;
    int chunk_center;
    int shell_center;
    int shell_type;
    int am2shell;
    int max_am;
    int ndocc;
    int num_shells;
    
    double chunk_rad_in;
    double chunk_rad_out;
    double bragg;
    double dist;
    double r,rr;
    double bastmp;
    
    struct coordinates atom_point_geom;
    struct coordinates shell_geom;
    
    struct atomic_grid_s *atom_grid;
    struct leb_chunk_s *chunk;
    
    num_shells = BasisSet.num_shells;
    ndocc = MOInfo.ndocc;
    max_am = BasisSet.max_am;
    
    atom_grid = &(DFT_options.grid.atomic_grid[atom_num]);
    bragg = atom_grid->Bragg_radii;
    atom_point_geom = atom_grid->atom_center;
    chunk_center = atom_grid->atom_num;
    
    chunk = &(atom_grid->leb_chunk[chunk_num]);   
    chunk_rad_in = chunk->spheres[0].r*bragg;
    chunk_rad_out = chunk->spheres[chunk->size-1].r*bragg;
    
    bzero(DFT_options.close_shell_info.close_shells_per_am
	  ,sizeof(int)*max_am);
    
    DFT_options.close_shell_info.num_close_aos = 0;
    
    j = 0;
    l = 0;
    for(i=0;i<num_shells;i++){
      
      am2shell = BasisSet.am2shell[i];
      shell_center = BasisSet.shells[am2shell].center - 1;
      shell_geom = Molecule.centers[shell_center];
      shell_type = BasisSet.shells[am2shell].am;
      
      /* Compute the Distance between the chunk's center
	 and the shell's center */
      
      dist = distance_calc(atom_point_geom,shell_geom);
      
      
      if(dist >= chunk_rad_in){
	if(dist <= chunk_rad_out){
	  /* Doesn't matter because the atom is in the chunk */
	  /* so its basis functions will atomatically be accepted */
	  bastmp = 1.0;
	}
	else{
	  /* dist is greater than the outer chunk, so
	     use the distance between the atom and the outer
	     sphere of the chunk */ 
	  
	  r = dist-chunk_rad_out;
	  rr = r*r;
	  bastmp = calc_radial_bas(am2shell,rr,r);
	}
	/* dist is less than the inner shell */
      }
      else{
	r = chunk_rad_in-dist;
	rr = r*r;
	bastmp = calc_radial_bas(am2shell,rr,r);
      }
      
      /*fprintf(outfile,"\natom x = %10.10lf y = %10.10lf z = %10.10lf"
	,atom_point_geom.x,atom_point_geom.y,atom_point_geom.z);
	fprintf(outfile,"\nshell x = %10.10lf y = %10.10lf z = %10.10lf"
	,shell_geom.x,shell_geom.y,shell_geom.z);
	fprintf(outfile,"\nshell_center = %d chunk_center = %d dist = %10.10lf",shell_center,chunk_center,dist);
	fprintf(outfile,"\ndist = %10.10lf",dist);
	fprintf(outfile,"\nchunk_rad_in = %10.10lf chunk_rad_out = %10.10lf",chunk_rad_in,chunk_rad_out);
	fprintf(outfile, "\nr = %e rr = %e",r,rr);
	fprintf(outfile, "\nbastmp = %e",bastmp);*/
      /* ---------------------------------
	 Determine whether the basis
	 the function is close or not
	 --------------------------------*/
      
      if(fabs(bastmp) > TOL){
	DFT_options.close_shell_info.shells_close_to_chunk[j] = am2shell;
	for(k=0;k<ioff[shell_type];k++){
	  for(m=0;m<ndocc;m++){
	    DFT_options.close_shell_info.close_COCC[l][m] 
	      = Cocc[BasisSet.shells[am2shell].fao-1+k][m];
	  }
	  DFT_options.close_shell_info.aos_close_to_chunk[l] = 
	    BasisSet.shells[am2shell].fao-1+k;
	  l++;
	}
	j++;
	DFT_options.close_shell_info.close_shells_per_am[shell_type-1]++;	  	    
      }
      DFT_options.close_shell_info.num_close_shells = j;
      DFT_options.close_shell_info.num_close_aos = l;
    }
    /*print_close_shell_info(DFT_options.close_shell_info);*/
  }
}}
