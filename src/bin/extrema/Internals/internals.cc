/*###########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief Member functions for <b>internals</b> derived base class. */
/*						Joseph P. Kenny 11/29/01
  ###########################################################################*/

#define EXTERN
#include"extrema.h"

namespace psi { namespace extrema {
double **symm_matrix_invert(double **A, int dim, int print_det, int redundant);
}}

using namespace psi::extrema;

/*--------------------------------------------------------------------------*/
/*! \fn internals::mem_alloc()
    \brief Allocates memory.

    Must be called by top-level class once dimensions determined. */
/*--------------------------------------------------------------------------*/

void internals::mem_alloc() {

    coord_base::mem_alloc();

    full_geom = init_matrix(num_entries,3);
    fgrads = init_array(fnum_coords);
    fcoords = init_array(fnum_coords);
    fcoords_old = init_array(fnum_coords);
    B = (double**) malloc(fnum_coords*sizeof(double*));
    B_red = init_matrix(num_coords,3*num_entries);
    G = init_matrix(fnum_coords, fnum_coords);
    A = init_matrix(3*num_entries, fnum_coords);
    u = init_matrix(3*num_entries, 3*num_entries);
    
    return;
}

/*---------------------------------------------------------------------------*/
/*! \fn internals::grad_trans()
  \brief Performs gradient transformation from cartesians to internals. */
/*---------------------------------------------------------------------------*/

void internals :: grad_trans() {

    int i,j;
    double **temp1, **temp2, **temp3;
    
    temp1 = init_matrix(fnum_coords,fnum_coords);
    temp2 = init_matrix(fnum_coords,fnum_coords);
    temp3 = init_matrix(fnum_coords,3*num_entries);
    
    mmult(B,0,B,1,temp1,0,fnum_coords,3*num_entries,fnum_coords,0);
    temp2 = symm_matrix_invert(temp1,fnum_coords,0,1);
    mmult(temp2,0,B,0,temp3,0,fnum_coords,fnum_coords,3*num_entries,0);

    for(i=0;i<fnum_coords;++i) {
	for(j=0;j<3*num_entries;++j) {
	    fgrads[i] += temp3[i][j] * c_grads[j];
	}
    }

    free_matrix(temp1,fnum_coords);
    free_matrix(temp2,fnum_coords);
    free_matrix(temp3,fnum_coords);
    
    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn internals::back_transform(double *c_new, double *c_old)
  \brief Performs iterative back transformation from internals to cartesians.
  
  This function requires the full array of coordinate values.
  \param *c_new new internal coordinate array
  \param *c_old old internal coordinate array */
/*---------------------------------------------------------------------------*/

void internals :: back_transform(double *c_new, double *c_old ) {

  int i, j, pos;
  double conv=1.0, *dq, *dx;

  dq = init_array(fnum_coords);
  dx = init_array(3*num_entries);         

  int loop=0;
  int converged=0;
  double criteria = 1.0e-14;
  double dx_sum; 

  fprintf(outfile,"\n\n\n");
  fprintf(outfile,"  --------------------------------------");
  fprintf(outfile,"--------------------------------------\n");
  fprintf(outfile,"  Performing iterative transformation to find");
  fprintf(outfile," new cartesian coordinates\n");
  fprintf(outfile,"  --------------------------------------");
  fprintf(outfile,"--------------------------------------");
  if(print_lvl > NORMAL_PRINT) {
      fprintf(outfile,"\n  Iter    dq (internals)          dx (cartesians)");
      fprintf(outfile,  "\n  ---- ----------------------  ");
      fprintf(outfile,"----------------------");  
  }
  else fprintf(outfile,"\n");

  if(print_lvl > RIDICULOUS_PRINT) 
      print_carts(1.0);

  while((!converged) && (loop<BT_LOOP) ) {
  
      /*compute A*/
      compute_G(); 
      compute_A(); 
      
      for(i=0;i<fnum_coords;++i) {
	  dq[i] = c_new[i] - c_old[i];
      }
      if(print_lvl >= RIDICULOUS_PRINT) {
          for(i=0;i<fnum_coords;++i)
              fprintf(outfile,"\n dq[%d]=%.20lf",i+1,dq[i]);
      }


      /*compute dx = A dq */
      for(i=0;i<3*num_entries;++i) 
          dx[i]=0;
      for(i=0;i<3*num_entries;++i) {
	  for(j=0;j<fnum_coords;++j) {
	      dx[i] += A[i][j] * dq[j];
	  }
      }

      pos=0;
      dx_sum = 0.0;
      int hack=0;
      for(i=0;i<3*num_entries;++i) {
	  /* hack to keep proper orientation */
	  if(fabs(carts[i])>ANOTHER_ZERO) 
	      carts[i] += dx[i];
	  else if(fabs(dx[i])>ANOTHER_ZERO)
	      hack = 1;
          dx_sum += sqrt(dx[i]*dx[i]);
      }
      dx_sum /= (3*num_entries);
      if(hack)
	  if(print_lvl >= RIDICULOUS_PRINT)
	      fprintf(outfile,
		      "\n  WARNING: Using hack to keep proper orientation"); 

      if(print_lvl >= RIDICULOUS_PRINT) {
	  for(i=0;i<3*num_entries;++i) 
	      fprintf(outfile,"\n dx[%d]=%.20lf",i+1,dx[i]);
      }
      
      cart_to_internal(&c_old);

      compute_B();
    
      conv=0;
      for(i=0;i<fnum_coords;++i) {
	  conv += sqrt((c_new[i] - c_old[i])*
		       (c_new[i] - c_old[i]));
      }
      conv /= fnum_coords;
      if(print_lvl>NORMAL_PRINT) {
	  fprintf(outfile,"\n  %3d  %.20lf  %.20lf",loop+1,conv,dx_sum);
      }	  
      if( (conv<BT_CONV) && (dx_sum<BT_CONV) )
	  converged = 1;
      ++loop;
  }
  
  if(!converged) { 
      fprintf(outfile,
	      "\n  Check for angles near 180.0 degrees, they're bad");
      fprintf(outfile,
	      "\n  You may need to use z-matrix coordinates avoiding 180.0");
      fprintf(outfile," degree angles\n");
      punt("Back transformation to cartesians has failed");
  }
  else
      if(print_lvl > NORMAL_PRINT)
	  fprintf(outfile,"\n  Back transformation to cartesians completed\n");

  

  return;
}








