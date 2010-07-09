/*! \file
    \ingroup OPTKING
    \brief NEW_GEOM performs the back-transformation to cartesian coordinates.
It computes a new cartesian geometry from an old cartesian geometry
and a set of internal coordinate displacements
it returns 1 if a new cartesian geometry was successfully determined
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

namespace psi { //namespace optking {

bool new_geom(const cartesians &carts, simples_class &simples, const salc_set &all_salcs,
    double *dq, int print_flag, int restart_geom_file,
    char *disp_label, int disp_num, int last_disp, double *return_geom) {

  int bmat_iter_done,count,i,j,dim_carts,natom,nsalcs;
  double **A, **G, **G_inv, **B, **u, **temp_mat, *no_fx;
  double dx_sum, dq_sum, *dx, *new_x, *x, *new_q, *q, *masses, *coord;

  natom = optinfo.natom;
  dim_carts = 3*natom;

  nsalcs = all_salcs.get_num();
  dx = init_array(dim_carts);
  new_x = init_array(dim_carts);
  new_q = init_array(nsalcs);

  masses = carts.get_mass();
  u = mass_mat(masses);
  free_array(masses);
//u = unit_matrix(dim_carts);

  A = block_matrix(dim_carts, nsalcs);
  G = block_matrix(nsalcs, nsalcs);
  G_inv = block_matrix(nsalcs, nsalcs);
  temp_mat = block_matrix(dim_carts, nsalcs);

  // Compute B matrix -- Isn't this slick?
  coord = carts.get_coord();
  simples.compute(coord);
  // fix configuration for torsions, sets flag for torsions > FIX_NEAR180 or < -FIX_NEAR180 
  simples.fix_near_180(); // subsequent computes will modify torsional values
  simples.compute_s(coord);
  free_array(coord);

  B = compute_B(simples, all_salcs);
  q = compute_q(simples, all_salcs);

  if (optinfo.print_debug_backtransformation) {
    fprintf(outfile,"Current q internal coordinates\n");
    for (i=0;i<nsalcs;++i) fprintf(outfile,"%15.10lf",q[i]);
    fprintf(outfile,"\n");
  }

  for (i=0;i<nsalcs;++i) {
    q[i] += dq[i];
  }

  if (optinfo.print_debug_backtransformation) {
    fprintf(outfile,"Target q internal coordinates\n");
    for (i=0;i<nsalcs;++i) fprintf(outfile,"%15.10lf",q[i]);
    fprintf(outfile,"\n");
  }

  x = carts.get_coord();
  scalar_mult(_bohr2angstroms,x,dim_carts); // x now holds geom in Ang

  fprintf(outfile,"\n\tBack-transformation to cartesian coordinates...\n");
  fprintf(outfile,"\t------------------------------------\n");
  fprintf(outfile,"\t Iter  RMS Delta(dx)  RMS Delta(dq)\n");
  fprintf(outfile,"\t------------------------------------\n");

  // Start back transformation iterations
  bmat_iter_done = 0;
  count = 0;
  do {
    free_block(G);
    free_block(G_inv);
    G = compute_G(B, nsalcs, carts);
    G_inv = symm_matrix_invert(G, nsalcs, 0,optinfo.redundant);

    if (optinfo.print_debug_backtransformation) {
      fprintf(outfile,"G matrix\n");
      print_mat2(G,nsalcs, nsalcs,outfile);
      fprintf(outfile,"\nG matrix inverted with redundant = %d\n", optinfo.redundant);
      fprintf(outfile,"G_inv matrix\n");
      print_mat2(G_inv,nsalcs, nsalcs,outfile);
    }

    // BMAT computes G_inv only once like the following.
    // OPTKING recomputes G_inv at each iteration, which
    // is slower but gives better convergence.
    //   if (count == 0) {
    //     G_inv = symm_matrix_invert(G,nsalcs,0,optinfo.redundant);
    //   }

    // u B^t G_inv = A
    opt_mmult(B,1,G_inv,0,temp_mat,0,dim_carts, nsalcs, nsalcs,0);
    opt_mmult(u,0,temp_mat,0,A,0,dim_carts,dim_carts, nsalcs, 0);
    // A dq = dx
    opt_mmult(A,0,&dq,1,&dx,1,dim_carts, nsalcs,1,0);

    if (optinfo.print_debug_backtransformation) {
      fprintf(outfile,"dx increments\n");
      for (i=0;i<dim_carts;++i)
        fprintf(outfile,"%15.10lf\n",dx[i]);
      fprintf(outfile,"\n");
    }

    // Compute new cart coordinates in au, then B matrix
    for (i=0;i<dim_carts;++i)
      new_x[i] = (x[i] + dx[i]) / _bohr2angstroms;

    if (optinfo.print_debug_backtransformation) {
      fprintf(outfile,"new x \n");
      for (i=0;i<dim_carts;++i)
        fprintf(outfile,"%15.10lf\n", new_x[i]);
      fprintf(outfile,"\n");
    }

    simples.compute(new_x);
    // simples.print(outfile,1);
    simples.compute_s(new_x);
    free_block(B);
    B = compute_B(simples, all_salcs);

    // compute new internal coordinate values
    free_array(new_q);
    new_q = compute_q(simples, all_salcs);

    if (optinfo.print_debug_backtransformation) {
      fprintf(outfile,"Obtained q internal coordinates\n");
      for (i=0;i<nsalcs;++i) fprintf(outfile,"%15.10lf",new_q[i]);
      fprintf(outfile,"\n");
    }

    for (i=0;i< nsalcs;++i)
      dq[i] = q[i] - new_q[i];
 
    if (optinfo.print_debug_backtransformation) {
       fprintf(outfile,"New internal coordinate errors dq\n");
       for (i=0;i<nsalcs;++i)
       fprintf(outfile,"%d, %15.10lf\n",i,dq[i]);
    }

    // Test for convergence of iterations
    dx_sum = dq_sum = 0.0;
    for (i=0;i<dim_carts;++i)
      dx_sum += dx[i]*dx[i];
    dx_sum = sqrt(dx_sum / ((double) dim_carts));

    for (i=0;i<nsalcs;++i)
      dq_sum += dq[i]*dq[i];
    dq_sum = sqrt(dq_sum / ((double) nsalcs));

   // fprintf(outfile,"\n");
   // for (i=0;i<nsalcs;++i)
     // if ( (dq[i]*dq[i]) > 1E-8 )
       // fprintf(outfile,"internal coordinate %d contributing: %15.10lf \n", i, dq[i]);

    if ((dx_sum < optinfo.bt_dx_conv) && (dq_sum < optinfo.bt_dq_conv))
      bmat_iter_done = 1;
    fprintf (outfile,"\t%5d %12.1e %12.1e\n", count+1, dx_sum, dq_sum);

    // store new x in angstroms 
    for (i=0;i<dim_carts;++i)
      x[i] = new_x[i] * _bohr2angstroms;

    ++count;
  } while( (bmat_iter_done == 0) && (count < optinfo.bt_max_iter) );

  free_block(B);

  fprintf(outfile,"\t------------------------------------\n");

  if (count >= optinfo.bt_max_iter) {
    fprintf(outfile,"Could not converge new geometry in %d iterations.\n",count);
    if (optinfo.mode == MODE_OPT_STEP)
      return false; /* let opt_step try smaller steps */
    else 
      exit(2);
  }
  else {
    fprintf(outfile, "\nSuccessfully converged to displaced geometry.\n");
  }

  // take x back to bohr
  scalar_mult(1.0/_bohr2angstroms, x, dim_carts);

//print_mat(&x,1,dim_carts,outfile);

  // avoid sending slightly non-symmetric geometries into chkpt file
  /*
  j=0;
  for (i=0;i<dim_carts;++i) {
    if ( fabs(x[i]) < MIN_CART_OUT) {
      x[i] = 0.0;
      ++j;
    }
  }
  if (j>0)
    fprintf(outfile,"Setting cartesian coodinates less than %e to zero.", MIN_CART_OUT);
    */

  /* Symmetry adapt the resulting geometry */
  if (optinfo.mode == MODE_OPT_STEP)
    symmetrize_geom(x);

  //  carts.set_coord(x); can't change calling geometry - may be reused
//  no_fx = new double [3*optinfo.nallatom];

//  for (i=0;i<3*optinfo.nallatom;++i) no_fx[i] = 0.0;
//  for (i=0;i<natom;++i) {
  // make fcoord for writing to chkpt - not sure why I have to do this
  // but wt_geom doesn't seem to work if dummy atoms are present
//    no_fx[3*optinfo.to_dummy[i]+0] = x[i*3+0];
//    no_fx[3*optinfo.to_dummy[i]+1] = x[i*3+1];
//    no_fx[3*optinfo.to_dummy[i]+2] = x[i*3+2];
//  }

  // cart_temp is used to return results to output.dat and chkpt
  cartesians cart_temp;
  cart_temp.set_coord(x);
//  cart_temp.set_fcoord(no_fx);

  fprintf(outfile,"\n%s\n",disp_label);
  cart_temp.print(1,outfile,0,disp_label,0);

  for (i=0;i<dim_carts;++i)
    return_geom[i] = x[i];

  // write geometry to chkpt or geom.dat
  if (print_flag == PRINT_TO_GEOM) {
    FILE *fp_geom;
    if (restart_geom_file) {
      opt_ffile(&fp_geom, "geom.dat",0);
      fprintf(fp_geom, "geom.dat: (\n");
    }
    else
      opt_ffile(&fp_geom, "geom.dat",1);

    cart_temp.print(10,fp_geom,restart_geom_file,disp_label,disp_num);
    if(last_disp) fprintf(fp_geom,")\n");
    fclose(fp_geom);
    fflush(outfile);
  }
  else if (print_flag == PSIF_CHKPT) {
    cart_temp.print(PSIF_CHKPT,outfile,0,disp_label,disp_num);
  }

  free_array(q); free_array(new_q);
  free_array(x); free_array(dx); free_array(new_x);
  free_block(A);
  free_block(G);
  free_block(G_inv);
  free_block(u);
  free_block(temp_mat);
  return true;
}

}//} /* namespace psi::optking */
