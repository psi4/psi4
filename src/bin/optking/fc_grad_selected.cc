/*! \file
    \ingroup OPTKING
    \brief fc_grad_selected(): computes force constants from gradients for selected coordinates
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <libqt/qt.h>
#include <libpsio/psio.h>

namespace psi { //namespace optking {

// only the symmetric salcs are passed in
// compute force constants for selected coordinates

void fc_grad_selected(const cartesians &carts, simples_class & simples, const salc_set &symm) {

  int i,j,k,a,b, ii, jj,cnt;
  double **B, **G, **G_inv, *masses, **u, *geom, *forces, **force_constants;
  double energy, *energies, **displacements, cm_convert;
  double *f, **all_f_q, *f_q, *temp_arr, *temp_arr2, **all_q, tval, **geom2D;
  double **evects, *evals, **FG, tmp;
  double *micro_geom, *micro_grad, *grad;

  int ndisps ;     // # of displacements
  int ncoord ;     // ndisps/2, the # of coordinates to compute fc for
  int nsymm = symm.get_num(); // # of symmetric coordinates
  int dim_carts = 3*carts.get_natom();
  int *coord2salc; // list of absolute salc indices for coordinates to displace

  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Num. of disp.", (char *) &(ndisps), sizeof(int));

  ncoord = ndisps / 2;

  micro_grad = new double [ndisps*dim_carts];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(micro_grad[0]), ndisps*dim_carts*sizeof(double));

  micro_geom = new double [ndisps*dim_carts];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced geometries",
      (char *) &(micro_geom[0]), ndisps*dim_carts*sizeof(double));

  coord2salc = init_int_array(ncoord);
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced coords to Salc Index",
      (char *) coord2salc, ncoord * sizeof(int));

  close_PSIF();

  all_q = (double **) malloc(ndisps*sizeof(double *));
  all_f_q = (double **) malloc(ndisps*sizeof(double *));
  f = (double *) malloc(dim_carts*sizeof(double));

  // compute forces in internal coordinates, f_q = G_inv B u f
  for (i=0; i<ndisps; ++i) {
    simples.compute(&(micro_geom[i*dim_carts]));
    simples.compute_s(&(micro_geom[i*dim_carts]));

    all_q[i] = compute_q(simples, symm);
    B = compute_B(simples,symm);
    G = compute_G(B,nsymm,carts);
    G_inv = symm_matrix_invert(G,nsymm,0,optinfo.redundant);
    masses = carts.get_fmass();
    u = mass_mat(masses);

    // load up forces
    for (j=0;j<dim_carts;++j)
      f[j] = micro_grad[i*dim_carts+j] * -1.0 * _hartree2J * 1.0E18 /
        _bohr2angstroms;

    all_f_q[i] = init_array(nsymm);
    temp_arr = init_array(nsymm);
    temp_arr2 = init_array(dim_carts);

    opt_mmult(u,0,&f,1,&temp_arr2,1,dim_carts,dim_carts,1,0);
    opt_mmult(B,0,&temp_arr2,1,&temp_arr,1,nsymm,dim_carts,1,0);
    opt_mmult(G_inv,0,&temp_arr,1,&(all_f_q[i]),1,nsymm,nsymm,1,0);

    free_array(temp_arr);
    free_array(temp_arr2);
    free_block(u);
    free_block(B);
    free_block(G);
    free_block(G_inv);
  }

  free_array(f);

  if (optinfo.print_hessian) {
    fprintf(outfile,"Internal coordinate forces for each displacement\n");
    cnt = 0;
    for (i=0; i<ndisps; ++i) {
      for (j=0; j<nsymm; ++j) {
        fprintf(outfile, " %12.8lf", all_f_q[i][j]);
        ++cnt;
        if ((cnt == 8) && (j+1 != nsymm)) { fprintf(outfile,"\n"); cnt = 0; }
      }
      fprintf(outfile,"\n");
    }
  }

  fprintf(outfile,"\nApplying 3-point formula\n");
  force_constants = block_matrix(nsymm,nsymm);
  for (i=0;i<ncoord;++i) {
    for (j=0;j<nsymm;++j) {
      ii = coord2salc[i];
      // reverse sign for forces not gradients
      force_constants[ii][j] = force_constants[j][ii] =
        (all_f_q[2*i][j]-all_f_q[2*i+1][j]) / (2.0 * optinfo.disp_size);
    }
  }

  free_block(all_f_q);

  //  fprintf(outfile,"\nForce Constants\n");
  //  print_mat(force_constants, nsymm, nsymm, outfile);

  // supplement empirical or old force constants with new ones

  i=0;
  open_PSIF();
  if (psio_tocscan(PSIF_OPTKING, "Symmetric Force Constants") == NULL) i = 1;
  close_PSIF();
  if (i) {
    empirical_H(simples, symm, carts);
  }

  double **fc_old = block_matrix(nsymm,nsymm);
  open_PSIF();
  fprintf(outfile,"Reading force constants from PSIF_OPTKING\n");
  psio_read_entry(PSIF_OPTKING, "Symmetric Force Constants",
    (char *) &(fc_old[0][0]),nsymm*nsymm*sizeof(double));

  for (i=0;i<ncoord;++i) {
    for (j=0;j<nsymm;++j) {
      ii = coord2salc[i];
      fc_old[ii][j] = fc_old[j][ii] = force_constants[ii][j];
    }
  }

  fprintf(outfile,"Merging computed force constants and saving in PSIF_OPTKING.\n");
  psio_write_entry(PSIF_OPTKING, "Symmetric Force Constants",
    (char *) &(fc_old[0][0]),nsymm*nsymm*sizeof(double));

  fprintf(outfile,"\nMerged Force Constants\n");
  print_mat(fc_old, nsymm, nsymm, outfile); fflush(outfile);

  free_block(fc_old);
  free_block(force_constants);

  optinfo.disp_num = 0;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num),sizeof(int));
  close_PSIF();
}

}//} /* namespace psi::optking */

