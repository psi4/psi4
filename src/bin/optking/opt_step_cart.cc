/*! \file
    \ingroup OPTKING
    \brief OPT_STEP_CART.CC takes geometry step in cartesian coordinates
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

namespace psi { namespace optking {

int opt_step_cart(cartesians &carts) {
  int xyz, i,j,k,ii,a,b, dim, dim_carts, success,nbfgs, constraint, nf;
  double **H_inv, **temp_mat, **u, **P;
  double *f, *x, *dx, *q, tval, tval2, scale, temp;
  double **geom_A, **geom_B;
  char *disp_label, *wfn, value_string[30], force_string[30];

  dim_carts = carts.get_natom()*3;
  disp_label = new char[MAX_LINELENGTH];

  ip_string("WFN", &(wfn),0);
  fprintf(outfile,"\nCurrent %s energy before step %20.10lf\n", wfn, carts.get_energy());
  free(wfn);
  fprintf(outfile,"\nTaking geometry step number %d\n",optinfo.iteration+1);

  // masses = carts.get_mass();
  // u = mass_mat(masses);
  // free_array(masses);

  x = carts.get_coord();
  scalar_mult(_bohr2angstroms,x,dim_carts); // x holds geom in Ang
  //fprintf(outfile,"Cartesians (Ang) \n");
  //mat_print(&x, 1, 9, outfile);

  // get forces array in cartesian coordinates, f, (in aJ/Ang)
  f = carts.get_forces();
  //fprintf(outfile,"Forces in aJ/Ang\n");
  //mat_print(&f, 1, 9, outfile);

  // read/write BFGS stuff
  open_PSIF();
  nbfgs = 0;
  if (psio_tocscan(PSIF_OPTKING, "Num. of Previous Entries") != NULL)
    psio_read_entry(PSIF_OPTKING, "Num. of Previous Entries", (char *) &nbfgs, sizeof(int));

  sprintf(value_string,"Previous Cartesian Values %d", nbfgs);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) &(x[0]),
      dim_carts*sizeof(double));

  sprintf(force_string,"Previous Cartesian Forces %d", nbfgs);
  psio_write_entry(PSIF_OPTKING, force_string, (char *) &(f[0]),
      dim_carts*sizeof(double));

  ++nbfgs;
  psio_write_entry(PSIF_OPTKING, "Num. of Previous Entries", (char *) &nbfgs, sizeof(int));
  close_PSIF();

  tval = 0.0;
  tval2 = fabs(f[0]);
  for (i=0;i<dim_carts;++i) {
    tval += SQR(f[i]);
    if (fabs(f[i]) > tval2) tval2 = fabs(f[i]);
  }
  tval = tval/((double) dim_carts);
  tval = sqrt(tval);
  fprintf(outfile,"MAX force: %15.10e\n",tval2);
  fprintf(outfile,"RMS force: %15.10e\n",tval);

  if (tval2 < optinfo.conv_max_force) {
    fprintf(outfile,"\nMAX force is < %5.1e.  Optimization is complete.\n", optinfo.conv_max_force);
    ip_string("WFN", &(wfn),0);
    fprintf(outfile,"Final %s energy is %15.10lf\n", wfn, carts.get_energy());
    free(wfn);
    fprintf(stderr,"\n  OPTKING:  optimization is complete\n");
    fprintf(outfile,"The Optimized geometry in a.u.\n");
    carts.print(12,outfile,0,disp_label,0);
    fprintf(outfile,"\nThe Optimized geometry in Angstrom\n");
    carts.print(13,outfile,0,disp_label,0);
    fprintf(outfile, "\n");
    fflush(outfile);

    if (optinfo.zmat) {
      int *unique_zvars;
      unique_zvars = init_int_array(MAX_ZVARS);
      //compute_zmat(carts, unique_zvars);
      //print_zmat(outfile, unique_zvars);
      free_int_array(unique_zvars);
      fprintf(outfile,"\n");
    }
    free_array(f);
    return(PSI_RETURN_ENDLOOP);
  } /* end converged geometry */

  fconst_init_cart(carts);
  /* Take geometry step in cartesian coordinates */

  // P could contain constraints later
  P = unit_matrix((long int) dim_carts);
  H_inv = compute_H_cart(carts, P);

  /* compute cartesian step */
  dx = init_array(dim_carts);
  opt_mmult(H_inv,0,&f,1,&dx,1,dim_carts,dim_carts,1,0);

  /* scale stepsize */
  scale = 1.0;
  for (i=0;i<dim_carts;++i) {
    if (fabs(dx[i]) > optinfo.step_limit_cart)
      scale = optinfo.step_limit_cart / fabs(dx[i]);
  }
  fprintf(outfile,"\nScaling displacements by %lf\n",scale);
  scalar_mult(scale, dx, dim_carts);
  //fprintf(outfile,"Displacements dx in Ang\n");
  //mat_print(&dx, 1, 9, outfile);

  for (i=0; i<dim_carts; ++i)
    x[i] += dx[i];

  scalar_mult(1.0/_bohr2angstroms, x, dim_carts);
  symmetrize_geom(x);
  carts.set_coord(x);
  carts.print(PSIF_CHKPT,outfile,0,disp_label,0);
  fprintf(outfile,"\nNew Cartesian Geometry in a.u.\n");
  carts.print(1, outfile,0, disp_label, 0);

  free_array(f); free_array(x); free_array(dx);

  optinfo.iteration += 1;
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Iteration", (char *) &(optinfo.iteration),sizeof(int));
  close_PSIF();
  delete [] disp_label;
  return 0;
}

}} /* namespace psi::optking */

