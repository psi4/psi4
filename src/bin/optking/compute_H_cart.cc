/*! \file
    \ingroup OPTKING
    \brief This function reads in Force Constants H from
    PSIF_OPTKING (in redundant internal coodinates) does a
    BFGS update on H inverts H to form H_inv and returns H_inv.
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>

namespace psi { //namespace optking {

double **compute_H_cart(cartesians & carts, double **P) {
  double **H, **H_inv, **H_inv_new, **H_new, **temp_mat;
  int i,j,dim, n_previous;
  char buffer[MAX_LINELENGTH];

  dim = 3*carts.get_natom();
  H = block_matrix(dim,dim);

  /*** Read in force constants from PSIF_OPTKING ***/
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "Cartesian Force Constants",
      (char *) &(H[0][0]),dim*dim*sizeof(double));
  close_PSIF();
  fprintf(outfile,"\nForce Constants read from PSIF_OPTKING\n");

  // Do BFGS update on H if desired and possible
  // Current point has already been put in PSIF_OPTKING by opt_step()
  n_previous = 0;
  open_PSIF();
  if ( psio_tocscan(PSIF_OPTKING, "Num. of Previous Entries") != NULL)
    psio_read_entry(PSIF_OPTKING, "Num. of Previous Entries", (char *) &n_previous, sizeof(int));
  close_PSIF();

  if (optinfo.H_update == OPTInfo::NONE || n_previous < 2 )
    fprintf(outfile,"\nNo Hessian update performed.\n");
  else
    H_update_cart(H, carts);

  if (optinfo.print_hessian) {
    fprintf(outfile,"The Hessian (Second Derivative) Matrix\n");
    print_mat5(H,dim,dim,outfile);
  }

  if (optinfo.redundant) { // could add constraints later
    temp_mat = block_matrix(dim,dim);
    opt_mmult(P,0,H,0,temp_mat,0,dim,dim,dim,0);
    opt_mmult(temp_mat,0,P,0,H,0,dim,dim,dim,0);
    free_block(temp_mat);

    temp_mat = unit_matrix((long int) dim);
    for (i=0;i<dim;++i)
      for (j=0;j<dim;++j) {
          H[i][j] += 1000 * (temp_mat[i][j] - P[i][j]);
      }
    free_block(temp_mat);
  }

  H_inv = symm_matrix_invert(H,dim,0,1);

  /*** write new force constants H to PSIF_OPTKING ***/
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Cartesian Force Constants",
      (char *) &(H[0][0]),dim*dim*sizeof(double));
  close_PSIF();
  free_block(H);
  return H_inv;
}

/* This functions performs update of Hessian */
void H_update_cart(double **H, cartesians &carts) {
  int i,j,dim,n_previous,i_step, natom, step_start, skip;
  double xx, xg, xz, zz, *f, *f_old, *Z;
  double *dx, *dg, **X, **temp_mat, *x, *x_old, phi;
  char force_string[30], x_string[30];

  if (optinfo.H_update == OPTInfo::BFGS)
    fprintf(outfile,"\nPerforming BFGS update");
  else if (optinfo.H_update == OPTInfo::MS)
    fprintf(outfile,"\nPerforming Murtagh/Sargent update");
  else if (optinfo.H_update == OPTInfo::POWELL)
    fprintf(outfile,"\nPerforming Powell update");
  else if (optinfo.H_update == OPTInfo::BOFILL)
    fprintf(outfile,"\nPerforming Bofill update");

  dim = carts.get_natom() * 3;

  f = init_array(dim);
  f_old = init_array(dim);
  dx    = init_array(dim);
  dg    = init_array(dim);
  x = init_array(dim);
  x_old = init_array(dim);

  /*** read/compute current internals and forces from PSIF_OPTKING ***/
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "Num. of Previous Entries", (char *) &n_previous, sizeof(int));

  sprintf(force_string,"Previous Cartesian Forces %d", n_previous-1);
  psio_read_entry(PSIF_OPTKING, force_string, (char *) &(f[0]), dim * sizeof(double));

  sprintf(x_string,"Previous Cartesian Values %d", n_previous-1);
  psio_read_entry(PSIF_OPTKING, x_string, (char *) &(x[0]), dim*sizeof(double));
  close_PSIF();

  if (optinfo.H_update_use_last == 0) { /* include all available gradients */
    step_start = 0;
  }
  else {
    step_start = n_previous - optinfo.H_update_use_last - 1;
    if (step_start < 0)
      step_start = 0;
  }

  fprintf(outfile," with previous %d gradient(s).\n", n_previous-1-step_start);

  for (i_step=step_start; i_step<(n_previous-1); ++i_step) {
    /* read/compute old internals and forces from PSIF_OPTKING ***/
    sprintf(force_string,"Previous Cartesian Forces %d", i_step);
    sprintf(x_string,"Previous Cartesian Values %d", i_step);
    open_PSIF();
    psio_read_entry(PSIF_OPTKING, force_string, (char *) &(f_old[0]), dim * sizeof(double));
    psio_read_entry(PSIF_OPTKING, x_string, (char *) &(x_old[0]), dim*sizeof(double));
    close_PSIF();

    // compute delta(coordinate) and delta(gradient)
    for (i=0;i<dim;++i) {
      dx[i] = x[i] - x_old[i];
      dg[i] = (-1.0) * (f[i] - f_old[i]); // gradients -- not forces!
    }

    if (optinfo.H_update == OPTInfo::BFGS) {
      // Do BFGS update: Schlegel 1987 Ab Initio Methods in Quantum Chemistry 
      // Let a = DQT.DG and X = (I - DQ*DGT/a)
      // Then Hk = X * H_(k-1) * XT + DQ*DQT/a
      // To make formula work for Hessian (some call it "B")
      // you have to switch DQ and DG in the equation
  
      dot_array(dx,dg,dim,&xg);
      X = unit_matrix((long int) dim);
      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          X[i][j] -= (dg[i] * dx[j]) / xg ; 
  
      temp_mat = block_matrix(dim,dim);
      opt_mmult(X,0,H,0,temp_mat,0,dim,dim,dim,0);
      opt_mmult(temp_mat,0,X,1,H,0,dim,dim,dim,0);
      free_block(temp_mat);
      free_block(X);
    
      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          H[i][j] += dg[i] * dg[j] / xg ; 
    }
    else if (optinfo.H_update == OPTInfo::MS) {
      // Equations taken from Bofill article below
      Z = init_array(dim);
      opt_mmult(H,0,&dx,1,&Z,1,dim,dim,1,0);
      for (i=0;i<dim;++i)
        Z[i] = dg[i] - Z[i];

      dot_array(dx,Z,dim,&xz);

      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          H[i][j] += Z[i] * Z[j] / xz ;

      free_array(Z);
    }
    else if (optinfo.H_update == OPTInfo::POWELL) {
      // Equations taken from Bofill article below
      Z = init_array(dim);
      opt_mmult(H,0,&dx,1,&Z,1,dim,dim,1,0);
      for (i=0;i<dim;++i)
        Z[i] = dg[i] - Z[i];

      dot_array(dx,Z,dim,&xz);
      dot_array(dx,dx,dim,&xx);

      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          H[i][j] += -1.0*xz/(xx*xx)*dx[i]*dx[j] + (Z[i]*dx[j] + dx[i]*Z[j])/xx;

      free_array(Z);
    }
    else if (optinfo.H_update == OPTInfo::BOFILL) {
      /* This functions performs a Bofill update on the Hessian according to
      J. M. Bofill, J. Comp. Chem., Vol. 15, pages 1-11 (1994). */
      // Bofill = (1-phi) * MS + phi * Powell
      Z = init_array(dim);
      opt_mmult(H,0,&dx,1,&Z,1,dim,dim,1,0);
      for (i=0;i<dim;++i)
        Z[i] = dg[i] - Z[i];

      dot_array(dx,Z,dim,&xz);
      dot_array(dx,dx,dim,&xx);
      dot_array(Z,Z,dim,&zz);

      phi = 1.0 - xz*xz/(xx*zz);
      if (phi < 0.0) phi = 0.0;
      if (phi > 1.0) phi = 1.0;

      for (i=0;i<dim;++i) 
        for (j=0;j<dim;++j)
          H[i][j] += (1.0-phi) * Z[i] * Z[j] / xz ;

      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          H[i][j] += phi*(-1.0*xz/(xx*xx)*dx[i]*dx[j] + (Z[i]*dx[j] + dx[i]*Z[j])/xx);
      free_array(Z);
    }
  } //end over old steps

  free_array(f);
  free_array(f_old);
  free_array(dx);
  free_array(dg);
  free_array(x);
  free_array(x_old);
  return;
}

/*! \file fconst_init_cart()
    \ingroup OPTKING
    \brief 
  FCONST_INIT -- make sure there are _some_ force constants in PSIF_OPTKING
  1) Confirm PSIF_OPTKING has them
  2) read them from FCONST: section of input
  3) read them from fconst.dat
  4) generate empirical force constants
*/

void fconst_init_cart(cartesians &carts) {
  int i, j, dim, count, constants_in_PSIF, cnt;
  char *buffer;
  double tval, **F, **temp_mat, *x, *x2, *x3, *z;
  buffer = new char[MAX_LINELENGTH];

  dim = carts.get_natom()*3;

  open_PSIF();
  if (psio_tocscan(PSIF_OPTKING, "Cartesian Force Constants") != NULL) {
    close_PSIF();
    return;
  }
  close_PSIF();

   /* read force constants from fconst section of input */
  if (ip_exist(":FCONST",0) ) {
    ip_cwk_add(":FCONST");
    if (ip_exist("CART_FC",0)) {
      fprintf(outfile,"Reading force constants from FCONST: \n");
      temp_mat = block_matrix(dim,dim);
      ip_count("CART_FC",&i,0);
      if (i != (dim*(dim+1))/2) {
        fprintf(outfile,"fconst has wrong number of entries\n");
        exit(2);
      }
      cnt = -1;
      for (i=0;i<dim;++i) {
        for (j=0;j<=i;++j) {
          ++cnt;
          ip_data("CART_FC","%lf",&(temp_mat[i][j]), 1, cnt);
        }
      }
      for (i=0;i<dim;++i)
        for (j=0;j<=i;++j)
          temp_mat[j][i] = temp_mat[i][j];
  
      /*** write to PSIF_OPTKING ***/
      open_PSIF();
      psio_write_entry(PSIF_OPTKING, "Cartesian Force Constants",
          (char *) &(temp_mat[0][0]),dim*dim*sizeof(double));
      close_PSIF();
  
      free_block(temp_mat);
      return;
    }
  }

  temp_mat = block_matrix(dim,dim);
  fprintf(outfile,"Making diagonal Hessian guess.\n");
  for (i=0; i<dim; ++i)
    temp_mat[i][i] = 1.0;

  // fprintf(outfile,"Hessian guess\n");
  // print_mat2(temp_mat,dim,dim,outfile);

  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Cartesian Force Constants",
      (char *) &(temp_mat[0][0]),dim*dim*sizeof(double));
  close_PSIF();
  free_block(temp_mat);

  delete [] buffer;
}

}//} /* namespace psi::optking */

