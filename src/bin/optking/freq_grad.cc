/*! \file
    \ingroup OPTKING
    \brief freq_grad_irrep(): computes frequencies from gradients for an irrep block
           freq_grad_nosymm(): computes frequencies from gradients for all coordinates
           ignoring symmetry
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

namespace psi { namespace optking {

void freq_grad_irrep(const cartesians &carts, simples_class &simples, const salc_set &all_salcs) {

  int i,j,ii,jj,k,a,b, cnt, dim, dim_carts, ndisps,irrep;
  int nirr_salcs, nsalcs, *irrep_salcs;
  double **B, **G, **G_inv, *masses, **u, *geom, *forces, **force_constants;
  double energy, *energies, **displacements, cm_convert;
  double *f, *f_q, *temp_arr, *temp_arr2, *q, tval, **geom2D;
  double **all_f_q; // internal coordinate forces for all unique displacements
  double **full_all_f_q; // internal coordinate forces for all displacements
  double **evects, *evals, **FG;
  double *micro_geom, *micro_grad, *grad, tmp, **force_constants_symm;
  char *salc_lbl;

  dim_carts = 3*carts.get_natom();
  irrep_salcs = new int[all_salcs.get_num()]; /* irrep -> total salc list lookup */
  irrep = optinfo.irrep;
  nsalcs = all_salcs.get_num();

  /* count and identify IRREP salcs */
  nirr_salcs = 0;
  cnt = 0;
  for (i=0; i<all_salcs.get_num(); ++i) {
    salc_lbl = all_salcs.get_label(i);
    if ( strcmp(salc_lbl, syminfo.irrep_lbls[irrep]) == 0) {
      ++nirr_salcs;
      irrep_salcs[cnt++] = i;
    }
  }

  // assume all coordinates in SYMM list are symmetric
  if ( (nirr_salcs == 0) && (optinfo.mode == MODE_FREQ_GRAD_IRREP) && (optinfo.irrep == 0) ) {
    fprintf(outfile,"\tNo proper irrep labels were detected.\n");
    fprintf(outfile,"\tAssuming all coordinates in SYMM vector are symmetric.\n");
      
    salc_set symm_salcs("SYMM");
    nirr_salcs = symm_salcs.get_num();
    for (i=0; i<nirr_salcs; ++i)
      irrep_salcs[i] = i;
  }

  if (nirr_salcs == 0) { fprintf(outfile,"No coordinates of irrep %d\n.", irrep); return;  }
  fprintf(outfile,"Found %d salcs of this irrep\n",nirr_salcs);

  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Num. of disp.",
      (char *) &(ndisps), sizeof(int));

  micro_grad = new double [ndisps*3*carts.get_natom()];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(micro_grad[0]), ndisps*3*carts.get_natom()*sizeof(double));

  micro_geom = new double [ndisps*3*carts.get_natom()];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced geometries",
      (char *) &(micro_geom[0]), ndisps*3*carts.get_natom()*sizeof(double));
  close_PSIF();

  // compute forces in internal coordinates for all disps, f_q = G_inv B u f
  all_f_q = init_matrix(ndisps, nsalcs);
  f = init_array(dim_carts);
  temp_arr = init_array(nsalcs);
  temp_arr2 = init_array(dim_carts);
  masses = carts.get_fmass();
  u = mass_mat(masses);
  for (i=0; i<ndisps; ++i) {

    simples.compute(&(micro_geom[i*dim_carts]));
    simples.compute_s(&(micro_geom[i*dim_carts]));
    q = compute_q(simples, all_salcs);
    /* fprintf(outfile,"Values of internal coordinates, displacement %d\n",i);
     for (j=0; j<salcs.get_num();++j) fprintf(outfile,"%15.10lf",all_q[i][j]);
     fprintf(outfile,"\n"); */

    B = compute_B(simples,all_salcs);
    G = compute_G(B, nsalcs, carts);
    fprintf(outfile,"BuB^t ");
    G_inv = symm_matrix_invert(G, nsalcs, 1, optinfo.redundant);

    for (j=0;j<dim_carts;++j)
      f[j] = micro_grad[i*dim_carts+j] * -1.0 * _hartree2J * 1.0E18 / _bohr2angstroms;

    opt_mmult(u,0,&f,1,&temp_arr2,1,dim_carts,dim_carts,1,0);
    opt_mmult(B,0,&temp_arr2,1,&temp_arr,1, nsalcs,dim_carts,1,0);
    opt_mmult(G_inv,0,&temp_arr,1,&(all_f_q[i]),1, nsalcs, nsalcs, 1, 0);

    free_array(q);
    free_matrix(B);
    free_matrix(G);
    free_matrix(G_inv);

  }
  free_array(f);
  free_array(temp_arr);
  free_array(temp_arr2);
  free_array(masses);
  free_matrix(u);

  /* expand unique displacements to redundant displacements */
  full_all_f_q = init_matrix( 2*all_salcs.get_num(), all_salcs.get_num());
  if (irrep == 0) {
    for (i=0; i<ndisps; ++i)
      for (j=0; j<nsalcs; ++j)
        full_all_f_q[i][j] = all_f_q[i][j];
  }
  else {
    for (i=0; i<ndisps; ++i) { // loop over displacements
      for (j=0; j<nsalcs; ++j) {
        full_all_f_q[2*irrep_salcs[i]][j] = all_f_q[i][j]; // the - disp
        salc_lbl = all_salcs.get_label(j);
        if ( strcmp(salc_lbl, syminfo.irrep_lbls[irrep]) == 0)
          full_all_f_q[2*irrep_salcs[i]+1][j] = -1 * all_f_q[i][j];
        else
          full_all_f_q[2*irrep_salcs[i]+1][j] = all_f_q[i][j];
      }
    }
  }
  free_matrix(all_f_q);

/*
fprintf(outfile,"full_all_f_q in freq_grad_irrep\n");
print_mat2(&(full_all_f_q[24]), 12, all_salcs.get_num(), outfile);
*/

  // apply three point formula - to generate force constants in this irrep block
  fprintf(outfile,"Applying %d-point formula\n",optinfo.points_freq_grad_ints);
  force_constants = init_matrix(nsalcs,nsalcs);
  for (i=0;i<nirr_salcs;++i)
    for (j=0;j<nirr_salcs;++j) {
      ii = irrep_salcs[i];
      jj = irrep_salcs[j];
      force_constants[ii][jj] = force_constants[jj][ii] =
        (full_all_f_q[2*ii][jj]-full_all_f_q[2*ii+1][jj]) / (2.0 * optinfo.disp_size);
    }

  free_matrix(full_all_f_q);

  fprintf(outfile,"\n\t** Force Constants\n");
  print_mat(force_constants, nsalcs, nsalcs, outfile);
  fflush(outfile);

  if (irrep == 0) {
    force_constants_symm = init_matrix(nirr_salcs,nirr_salcs);
    for (i=0;i<nirr_salcs;++i)
      for (j=0;j<nirr_salcs;++j) {
        ii = irrep_salcs[i];
        jj = irrep_salcs[j];
        force_constants_symm[i][j] = force_constants[ii][jj];
      }
    fprintf(outfile,"\n\t ** Writing symmetric force constants to PSIF_OPTKING ** \n");
    open_PSIF();
    psio_write_entry(PSIF_OPTKING, "Symmetric Force Constants",
      (char *) &(force_constants_symm[0][0]),nirr_salcs*nirr_salcs*sizeof(double));
    close_PSIF();
    if (optinfo.print_hessian)
      print_mat(force_constants_symm, nirr_salcs, nirr_salcs, outfile);

    if (optinfo.print_fconst) { // write force constants to fconst.dat
      opt_ffile_noexit(&fp_fconst, "fconst.dat",0);
      fprintf(fp_fconst,"FCONST: (\n symm_fc = (\n");
      for (i=0;i<nirr_salcs;++i) {
        cnt = 0;
        for (j=0;j<=i;++j) {
          fprintf(fp_fconst, "%10.6f", force_constants_symm[i][j]);
          if ( (++cnt > 7) && (j!=i) ) {
            fprintf(fp_fconst,"\n");
            cnt = 0;
          }
        }
        fprintf(fp_fconst,"\n");
        cnt = 0;
      }
      fprintf(fp_fconst," )\n)\n");
      fclose(fp_fconst);
    }
    /* reset BFGS update */
    open_PSIF();
    i = 0;
    if ( psio_tocscan(PSIF_OPTKING, "Num. of Previous Entries") != NULL)
      psio_write_entry(PSIF_OPTKING, "Num. of Previous Entries", (char *) &i, sizeof(int));
    close_PSIF();
  }

  // build G = BuB^t
  B = compute_B(simples, all_salcs);
  G = compute_G(B, nsalcs, carts);
  free_matrix(B);

  // compute FG and diagonalize 
  FG = init_matrix(nsalcs, nsalcs);
  opt_mmult(force_constants,0,G,0,FG,0,
      nsalcs,nsalcs,nsalcs,0);
  free_matrix(force_constants);
  free_matrix(G);

  //fprintf(outfile,"FG Matrix\n");
  //print_mat2(FG,nsalcs,nsalcs,outfile);
  //fflush(outfile);

  evals  = init_array(nsalcs);
  G = init_matrix(nsalcs, nsalcs);
  dgeev_optking(nsalcs, FG, evals, G);
  free_matrix(FG);
  free_matrix(G);

  cm_convert = 1.0/(2.0 * _pi * _c * 100.0);
  for (i=0;i<nsalcs;++i) {
    evals[i] = evals[i] * 1.0E-18 / ( 1.0E-20 * _amu2kg );
    if(evals[i] < 0.0) evals[i] = -1.0 * cm_convert * sqrt( -evals[i] );
    else evals[i] = cm_convert * sqrt( evals[i] );
  }

  fprintf(outfile,"\n  Harmonic Vibrational Frequencies in cm^(-1) for Irrep %s\n",
      syminfo.irrep_lbls[irrep]) ;
  fprintf(outfile,"  -----------------------------------------------------------\n");
  sort_vector(evals, nsalcs); /* ascending order */
  for (i=nsalcs-1; i>=0; --i) { /* descending order */
    if(evals[i] < 0.0) 
      fprintf(outfile,"%5d       %15.1lfi\n",nirr_salcs-i,-1.0*evals[i]);
    else
      fprintf(outfile,"%5d       %15.1lf\n",nirr_salcs-i,evals[i]);
  }
  free_array(evals);

  open_PSIF();
  optinfo.disp_num = 0;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num),sizeof(int));
  close_PSIF();

  delete [] irrep_salcs;
}


void freq_grad_nosymm(const cartesians &carts, simples_class &simples,
    const salc_set &all_salcs) {

  int i,j,k,a,b, ii, cnt, dim, dim_carts, ndisps;
  int nsalcs;
  double **B, **G, **G_inv, *masses, **u, *geom, *forces, **force_constants;
  double energy, *energies, **displacements, cm_convert;
  double *f, **all_f_q, *f_q, *temp_arr, *temp_arr2, **all_q, tval, **geom2D;
  double **evects, *evals, **FG, tmp;
  double *micro_geom, *micro_grad, *grad;

  dim_carts = 3*carts.get_natom();
  nsalcs = all_salcs.get_num();

  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Num. of disp.",
      (char *) &(ndisps), sizeof(int));
  if (ndisps != 2*nsalcs) 
    punt("Error: number of displacements is incorrect.");

  micro_grad = new double [ndisps*3*carts.get_natom()];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(micro_grad[0]), ndisps*3*carts.get_natom()*sizeof(double));

  micro_geom = new double [ndisps*3*carts.get_natom()];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced geometries",
      (char *) &(micro_geom[0]), ndisps*3*carts.get_natom()*sizeof(double));

  close_PSIF();

  all_q = (double **) malloc(ndisps*sizeof(double *));
  all_f_q = (double **) malloc(ndisps*sizeof(double *));
  f = (double *) malloc(dim_carts*sizeof(double));

  // compute forces in internal coordinates, f_q = G_inv B u f
  for (i=0; i<ndisps; ++i) {
    simples.compute(&(micro_geom[i*dim_carts]));
    simples.compute_s(&(micro_geom[i*dim_carts]));

    all_q[i] = compute_q(simples, all_salcs);
    B = compute_B(simples,all_salcs);
    G = compute_G(B,nsalcs,carts);
    G_inv = symm_matrix_invert(G,nsalcs,0,optinfo.redundant);
    masses = carts.get_fmass();
    u = mass_mat(masses);

    // load up forces
    for (j=0;j<dim_carts;++j)
      f[j] = micro_grad[i*dim_carts+j] * -1.0 * _hartree2J * 1.0E18 /
        _bohr2angstroms;

    all_f_q[i] = init_array(nsalcs);
    temp_arr = init_array(nsalcs);
    temp_arr2 = init_array(dim_carts);

    opt_mmult(u,0,&f,1,&temp_arr2,1,dim_carts,dim_carts,1,0);
    opt_mmult(B,0,&temp_arr2,1,&temp_arr,1,nsalcs,dim_carts,1,0);
    opt_mmult(G_inv,0,&temp_arr,1,&(all_f_q[i]),1,nsalcs,nsalcs,1,0);

    free_array(temp_arr);
    free_array(temp_arr2);
    free_matrix(u);
    free_matrix(B);
    free_matrix(G);
    free_matrix(G_inv);
  }

  free_array(f);

  // apply three point formula
  fprintf(outfile,"Applying %d-point formula\n",optinfo.points_freq_grad_ints);
  force_constants = init_matrix(nsalcs,nsalcs);
  for (i=0;i<nsalcs;++i)
    for (j=0;j<nsalcs;++j)
      force_constants[i][j] = force_constants[j][i] =
        (all_f_q[2*i][j]-all_f_q[2*i+1][j]) / (2.0 * optinfo.disp_size);

  //print_mat2(all_f_q, num_disps, salcs.get_num(), outfile);
  free_matrix(all_f_q);

  fprintf(outfile,"\nForce Constants\n");
  print_mat(force_constants, nsalcs, nsalcs, outfile);
  fflush(outfile);
  fprintf(outfile,"\t ** Writing force constants to PSIF_OPTKING ** \n");
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Symmetric Force Constants",
    (char *) &(force_constants[0][0]),nsalcs*nsalcs*sizeof(double));
  close_PSIF();

  // build G = BuB^t
  B = init_matrix(nsalcs, 3*carts.get_natom());
  B = compute_B(simples, all_salcs);
  G = compute_G(B, nsalcs, carts);
  free_matrix(B);

  // compute FG and diagonalize 
  FG = init_matrix(nsalcs, nsalcs);
  opt_mmult(force_constants,0,G,0,FG,0,nsalcs,nsalcs,nsalcs,0);
  free_matrix(force_constants);
  free_matrix(G);

  //fprintf(outfile,"FG Matrix\n");
  //print_mat2(FG,salcs.get_num(),salcs.get_num(),outfile);
  //fflush(outfile);

  evals  = init_array(nsalcs);
  G = init_matrix(nsalcs, nsalcs);
  dgeev_optking(nsalcs, FG, evals, G);
  free_matrix(FG);
  free_matrix(G);

  cm_convert = 1.0/(2.0 * _pi * _c * 100.0);
  for (i=0;i<nsalcs;++i) {
    evals[i] = evals[i] * 1.0E-18 / ( 1.0E-20 * _amu2kg );
    if(evals[i] < 0.0) evals[i] = -1.0 * cm_convert * sqrt( -evals[i] );
    else evals[i] = cm_convert * sqrt( evals[i] );
  }

  fprintf(outfile,"\nHarmonic Vibrational Frequencies\n");
  fprintf(outfile,"  -----------------------------------------------------------\n");
  sort_vector(evals, nsalcs); /* ascending order */
  for (i=nsalcs-1; i>= 0; --i) {  /* descending order */
    if(evals[i] < 0.0) 
      fprintf(outfile,"%5d       %15.1lfi\n",nsalcs-i,-1.0*evals[i]);
    else
      fprintf(outfile,"%5d       %15.1lf\n",nsalcs-i,evals[i]);
  }
  free_array(evals);

  open_PSIF();
  optinfo.disp_num = 0;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num),sizeof(int));
  close_PSIF();

}

}} /* namespace psi::optking */

