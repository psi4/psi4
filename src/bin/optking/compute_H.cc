/*
    \ingroup OPTKING
    \brief This function reads in Force Constants H from
    PSIF_OPTKING (in possible redundant internal coodinates) does a
    Hessian update on H, and inverts H to form H_inv and returns H_inv
    unless an RFO step is being used in which case H is returned
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <libpsio/psio.h>
#include <libipv1/ip_lib.h>

namespace psi { namespace optking {

static void H_update(double **H, simples_class &simples, const salc_set &symm, 
  const cartesians &carts);

double **compute_H(simples_class & simples, const salc_set &symm,
    double **P, const cartesians &carts) {

  double **H, **H_inv, **H_inv_new, **H_new, **temp_mat;
  int i,j,dim, n_previous;
  char buffer[MAX_LINELENGTH];

  dim = symm.get_num();
  H = init_matrix(dim,dim);

  // Read in force constants from PSIF_OPTKING
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "Symmetric Force Constants",
      (char *) &(H[0][0]),dim*dim*sizeof(double));
  close_PSIF();
  fprintf(outfile,"\nForce Constants read from PSIF_OPTKING\n");

  // Do BFGS update on H if desired and possible
  // Current point has already been put in PSIF_OPTKING by opt_step()
  n_previous = optinfo.iteration+1;

  if (optinfo.H_update == OPTInfo::NONE || n_previous < 2 || optinfo.balked_last_time)
    fprintf(outfile,"\nNo Hessian update performed.\n");
  else
    H_update(H, simples, symm, carts);

  if (optinfo.print_hessian) {
    fprintf(outfile,"The Hessian (Second Derivative) Matrix\n");
    print_mat5(H,dim,dim,outfile);
    fprintf(outfile,"\n");
  }

  // Project redundancies out of H according to Peng JCC 1996
  // and produce invertible Hessian
  // Form H <= PHP + 1000 (1-P)
  if ((optinfo.redundant) || (optinfo.constraints_present)) {
    temp_mat = init_matrix(dim,dim);
    opt_mmult(P,0,H,0,temp_mat,0,dim,dim,dim,0);
    opt_mmult(temp_mat,0,P,0,H,0,dim,dim,dim,0);
    free_matrix(temp_mat);

    temp_mat = unit_matrix((long int) dim);
    for (i=0;i<dim;++i)
      for (j=0;j<dim;++j) {
          H[i][j] += 1000 * (temp_mat[i][j] - P[i][j]);
      }
    free_matrix(temp_mat);
  }

  /*** write new force constants H to PSIF_OPTKING ***/
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Symmetric Force Constants",
      (char *) &(H[0][0]),dim*dim*sizeof(double));
  close_PSIF();

  return H;
}

/* This functions performs update of Hessian */
void H_update(double **H, simples_class &simples, const salc_set &symm, 
    const cartesians &carts) {

  int i,j,dim,n_previous,i_step, natom, step_start, skip;
  double qq, qg, qz, zz, *q, *f, *q_old, *f_old, *Z;
  double *dq, *dg, **X, **temp_mat, *x, *x_old, phi, **H_new, max;
  char force_string[30], x_string[30];

  if (optinfo.H_update == OPTInfo::BFGS)
    fprintf(outfile,"\nPerforming BFGS update");
  else if (optinfo.H_update == OPTInfo::MS)
    fprintf(outfile,"\nPerforming Murtagh/Sargent update");
  else if (optinfo.H_update == OPTInfo::POWELL)
    fprintf(outfile,"\nPerforming Powell update");
  else if (optinfo.H_update == OPTInfo::BOFILL)
    fprintf(outfile,"\nPerforming Bofill update");

  natom = carts.get_natom();
  dim = symm.get_num();

  f = init_array(dim);
  f_old = init_array(dim);
  dq    = init_array(dim);
  dg    = init_array(dim);
  x = init_array(3*natom);
  x_old = init_array(3*natom);
  H_new = init_matrix(dim,dim);

  /*** read/compute current internals and forces from PSIF_OPTKING ***/
  n_previous = optinfo.iteration+1;

  sprintf(force_string,"Internal Forces %d", n_previous-1);
  sprintf(x_string,"Cartesian Values %d", n_previous-1);
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, force_string, (char *) &(f[0]), dim * sizeof(double));
  psio_read_entry(PSIF_OPTKING, x_string, (char *) &(x[0]), 3*natom*sizeof(double));
  close_PSIF();

  simples.compute(x);
  simples.fix_near_180(); // fix configuration for torsions
  q = compute_q(simples, symm);

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
    sprintf(force_string,"Internal Forces %d", i_step);
    sprintf(x_string,"Cartesian Values %d", i_step);
    open_PSIF();
    psio_read_entry(PSIF_OPTKING, force_string, (char *) &(f_old[0]), dim * sizeof(double));
    psio_read_entry(PSIF_OPTKING, x_string, (char *) &(x_old[0]), 3*natom*sizeof(double));
    close_PSIF();

    simples.compute(x_old);
    q_old = compute_q(simples, symm);

    /* check for scary passages through 0 */
/*
    skip=0;
    for (i=0;i<dim;++i) {
      if (q[i] * q_old[i] < 0.0)
        skip = 1;
    }
    if (skip) {
      fprintf(outfile,"Warning a coordinate has passed through 0. Skipping Hessian update.\n");
      skip = 0;
      continue;
    } */

    // compute delta(coordinate) and delta(gradient)
    for (i=0;i<dim;++i) {
      dq[i] = q[i] - q_old[i];
      dg[i] = (-1.0) * (f[i] - f_old[i]); // gradients -- not forces!
    }
    // for (i=0;i<dim;++i) fprintf(outfile,"dq[%d]: %20.15lf\n",i,dq[i]);
    // for (i=0;i<dim;++i) fprintf(outfile,"dg[%d]: %20.15lf\n",i,dg[i]);

    if (optinfo.H_update == OPTInfo::BFGS) {
      // Do BFGS update: Schlegel 1987 Ab Initio Methods in Quantum Chemistry 
      // Let a = DQT.DG and X = (I - DQ*DGT/a)
      // Then Hk = X * H_(k-1) * XT + DQ*DQT/a
      // To make formula work for Hessian (some call it "B")
      // you have to switch DQ and DG in the equation
  
      dot_array(dq,dg,dim,&qg);
      // fprintf(outfile,"dq dot dg = %20.15lf\n", qg);

      X = unit_matrix((long int) dim);
      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          X[i][j] -= (dg[i] * dq[j]) / qg ; 
      //fprintf(outfile,"Xmat\n");
      //print_mat5(X,dim,dim,outfile);
  
      temp_mat = init_matrix(dim,dim);
      opt_mmult(X,0,H,0,temp_mat,0,dim,dim,dim,0);
      opt_mmult(temp_mat,0,X,1,H_new,0,dim,dim,dim,0);
      free_matrix(temp_mat);
      free_matrix(X);
    
      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          H_new[i][j] += dg[i] * dg[j] / qg ; 
    }
    else if (optinfo.H_update == OPTInfo::MS) {
      // Equations taken from Bofill article below
      Z = init_array(dim);
      opt_mmult(H,0,&dq,1,&Z,1,dim,dim,1,0);
      for (i=0;i<dim;++i)
        Z[i] = dg[i] - Z[i];

      dot_array(dq,Z,dim,&qz);

      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          H_new[i][j] = H[i][j] + Z[i] * Z[j] / qz ;

      free_array(Z);
    }
    else if (optinfo.H_update == OPTInfo::POWELL) {
      // Equations taken from Bofill article below
      Z = init_array(dim);
      opt_mmult(H,0,&dq,1,&Z,1,dim,dim,1,0);
      for (i=0;i<dim;++i)
        Z[i] = dg[i] - Z[i];

      dot_array(dq,Z,dim,&qz);
      dot_array(dq,dq,dim,&qq);

      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          H_new[i][j] = H[i][j] -1.0*qz/(qq*qq)*dq[i]*dq[j] + (Z[i]*dq[j] + dq[i]*Z[j])/qq;

      free_array(Z);
    }
    else if (optinfo.H_update == OPTInfo::BOFILL) {
      /* This functions performs a Bofill update on the Hessian according to
      J. M. Bofill, J. Comp. Chem., Vol. 15, pages 1-11 (1994). */
      // Bofill = (1-phi) * MS + phi * Powell
      Z = init_array(dim);
      opt_mmult(H,0,&dq,1,&Z,1,dim,dim,1,0);
      for (i=0;i<dim;++i)
        Z[i] = dg[i] - Z[i];

      dot_array(dq,Z,dim,&qz);
      dot_array(dq,dq,dim,&qq);
      dot_array(Z,Z,dim,&zz);

      phi = 1.0 - qz*qz/(qq*zz);
      if (phi < 0.0) phi = 0.0;
      if (phi > 1.0) phi = 1.0;

      for (i=0;i<dim;++i) 
        for (j=0;j<dim;++j)
          H_new[i][j] = H[i][j] + (1.0-phi) * Z[i] * Z[j] / qz ;

      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          H_new[i][j] += phi*(-1.0*qz/(qq*qq)*dq[i]*dq[j] + (Z[i]*dq[j] + dq[i]*Z[j])/qq);
      free_array(Z);
    }
  } //end over old steps

  free_array(x);

  // limit allowed changes to Hessian to 0.3 or 50%
  for (i=0;i<dim;++i)
    for (j=0;j<dim;++j)
      H_new[i][j] -= H[i][j];
   //fprintf(outfile,"H_new\n");
   //print_mat5(H_new,dim,dim,outfile);

  for (i=0;i<dim;++i) {
    for (j=0;j<dim;++j) {
      max = ((0.5*fabs(H[i][j]) > 0.3) ? (0.5*fabs(H[i][j])) : 0.3);

      if (fabs(H_new[i][j]) < max)
        H[i][j] += H_new[i][j];
      else
        H[i][j] += 0.3 * H_new[i][j]/fabs(H_new[i][j]);
    }
  }


/*
  // using fragments to make hessian block diagonal
  int *a2f = simples.atom2fragment(natom);
  for (i=0; i<dim; ++i)
    for (j=0; j<dim; ++j)
      if (a2f[i] != a2f[j]) H[i][j] = 0.0;
  delete [] a2f;
*/

  // put current values of internal coordinates back in place(!)
  x = carts.get_coord();
  simples.compute(x);
  simples.fix_near_180();
  simples.compute_s(x);
  free_array(x);

  free_array(q);
  free_array(f);
  free_array(q_old);
  free_array(f_old);
  free_array(dq);
  free_array(dg);
  free_array(x_old);
  free_matrix(H_new);
  return;
}

/*! \file fconst_init()
    \ingroup OPTKING
    \brief 
  FCONST_INIT -- make sure there are _some_ force constants in PSIF_OPTKING
  1) Confirm PSIF_OPTKING has them
  2) read them from FCONST: section of input
  3) read them from fconst.dat
  4) generate empirical force constants
*/

void fconst_init(const cartesians & carts, const simples_class & simples, const salc_set &symm) {
  int i, j, dim, count, constants_in_PSIF, cnt;
  char *buffer;
  double **F, **temp_mat;
  buffer = new char[MAX_LINELENGTH];

  open_PSIF();
  if (psio_tocscan(PSIF_OPTKING, "Symmetric Force Constants") != NULL) {
    close_PSIF();
    fprintf(outfile,"\tSymmetric Force Constants found in PSIF_OPTKING\n");
    return;
  }
  close_PSIF();

   /* read force constants from fconst section of input */
  if (ip_exist(":FCONST",0) ) {
    ip_cwk_add(":FCONST");
    fprintf(outfile,"Reading force constants from FCONST: \n");
    dim = symm.get_num();
    temp_mat = init_matrix(dim,dim);
    ip_count("SYMM_FC",&i,0);
    if (i != (dim*(dim+1))/2) {
      fprintf(outfile,"fconst has wrong number of entries\n");
      exit(2);
    }
    cnt = -1;
    for (i=0;i<dim;++i) {
      for (j=0;j<=i;++j) {
        ++cnt;
        ip_data("SYMM_FC","%lf",&(temp_mat[i][j]), 1, cnt);
      }
    }
    for (i=0;i<dim;++i)
      for (j=0;j<=i;++j)
        temp_mat[j][i] = temp_mat[i][j];

    /*** write to PSIF_OPTKING ***/
    open_PSIF();
    psio_write_entry(PSIF_OPTKING, "Symmetric Force Constants",
        (char *) &(temp_mat[0][0]),dim*dim*sizeof(double));
    close_PSIF();

    free_matrix(temp_mat);
    return;
  }

  opt_ffile_noexit(&fp_fconst, "fconst.dat",2);
  if (fp_fconst == NULL) { // generate empirical Hessian
    empirical_H(simples,symm,carts);
    return;
  }
  else { // read force constants from fconst.dat
    ip_append(fp_fconst, outfile);
    if (ip_exist(":FCONST",0) ) { // file has libipv1 format
      ip_cwk_add(":FCONST");
      fprintf(outfile,"Reading force constants from FCONST: \n");
      dim = symm.get_num();
      temp_mat = init_matrix(dim,dim);
      ip_count("SYMM_FC",&i,0);
      if (i != (dim*(dim+1))/2) {
        fprintf(outfile,"fconst has wrong number of entries\n");
        exit(2);
      }
      cnt = -1;
      for (i=0;i<dim;++i) {
        for (j=0;j<=i;++j) {
          ++cnt;
          ip_data("SYMM_FC","%lf",&(temp_mat[i][j]), 1, cnt);
        }
      }
    }
    fclose(fp_fconst);

    for (i=0;i<dim;++i)
      for (j=0;j<i;++j)
        temp_mat[j][i] = temp_mat[i][j];

    open_PSIF();
    psio_write_entry(PSIF_OPTKING, "Symmetric Force Constants",
        (char *) &(temp_mat[0][0]),dim*dim*sizeof(double));
    close_PSIF();
    free_matrix(temp_mat);
  }
  delete [] buffer;
}

}} /* namespace psi::optking */

