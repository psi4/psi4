/*! \file
    \ingroup OPTKING
    \brief OPT_STEP.CC takes geometry steps using
    gradients -- optking's default operation
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <libqt/qt.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>

namespace psi { namespace optking {

inline double rfo_energy(double rfo_t, double rfo_g, double rfo_h) {
  return (rfo_t * rfo_g + 0.5 * rfo_t * rfo_t * rfo_h)/(1 + rfo_t*rfo_t)/(_hartree2J*1.0e18);
}
inline double nr_energy(double rfo_t, double rfo_g, double rfo_h) {
  return (rfo_t * rfo_g + 0.5 * rfo_t * rfo_t * rfo_h)/(_hartree2J*1.0e18);
}

static bool line_search(cartesians & carts, int num_ints, double *dq);

int opt_step(cartesians &carts, simples_class &simples, const salc_set &symm) {

  int xyz, i,j,k,ii,a,b, dim, dim_carts, nsimples, constraint, cnt;
  double **B, **G, **G2, **G_inv, **H, **H_inv, **temp_mat, **u, **P;
  double DE_error, DE, DE_old, E_old, DE_new;
  double **C, **T, **T2, **T3, R_AB, theta_A, theta_B, tau, phi_A, phi_B;
  double *temp_arr2, *temp_arr, *masses, **geom, *forces, *coord, max_force, rms_force;
  double max_disp, rms_disp;
  double *f, *f_q, *dq, *dq_to_new_geom, *q, tval, tval2, tval3, tval4, scale, temp, *djunk;
  char *disp_label, *wfn, value_string[30], force_string[30];
  bool do_line_search = false, success;

  int rfo_root, rfo_num;
  double **rfo_mat, *lambda, *work, *rfo_old_evect, *rfo_u;
  double rfo_xnorm, rfo_g, rfo_h, rfo_t, rfo_eval;
  double nr_xnorm, nr_g, nr_h, *nr_u;

  // for analytic interfragment coordinates
  int isalc, I, simple_id, sub_index, sub_index2, atom;
  double **geom_A, **geom_B, **weight_A, **weight_B;
  double inter_q[6];
  Intco_type intco_type;

  dim_carts = carts.get_natom()*3;
  dim = symm.get_num();
  djunk = new double[dim_carts];
  disp_label = new char[MAX_LINELENGTH];

  if (dim == 0) punt("No symmetric internal coordinates to optimize.\n");

  ip_string("WFN", &(wfn),0);
  fprintf(outfile,"\nCurrent %s energy before step   %20.12lf\n", wfn, carts.get_energy());
  free(wfn);
  
  open_PSIF();
  DE = 1.0e6; // for if iteration == 0
  if ((optinfo.iteration > 0) && !optinfo.balked_last_time) {
    sprintf(value_string,"Energy %d", optinfo.iteration-1);
    psio_read_entry(PSIF_OPTKING, value_string, (char *) &E_old, sizeof(double));
    fprintf(outfile,"  Previous energy                %20.12lf\n", E_old);
    DE = carts.get_energy() - E_old;
    fprintf(outfile,"  Actual Delta(E) from last step %20.12lf\n", DE);
    sprintf(value_string,"DE prediction %d", optinfo.iteration-1);
    psio_read_entry(PSIF_OPTKING, value_string, (char *) &(DE_old), sizeof(double));
    DE_error = (DE_old-DE)/fabs(DE_old);
    fprintf(outfile,"  Previous Predicted Delta(E)    %20.12lf  (%5.1lf%% error)\n", DE_old,DE_error*100);
  }
  close_PSIF();

  // searching for a minimum, supersized steps downward permitted
  if (!optinfo.ts && (DE_error < -1*optinfo.step_energy_limit) && fabs(DE_old) > optinfo.line_search_min) {
    fprintf(outfile,"\n\tInsufficient energy drop observed.\n");
    do_line_search = true;
  }

  // searching for a TS, no supersized steps permitted
  if (optinfo.ts && (fabs(DE_error) > optinfo.step_energy_limit) && fabs(DE_old) > optinfo.line_search_min) {
    fprintf(outfile,"\t\tEnergy deviated too much from projection.\n");
    do_line_search = true;
  }

  // *** Bad step - step back and do line search
  if (do_line_search) {
    dq = init_array(dim);
    if (!line_search(carts,dim,dq) ) {
      fprintf(outfile,"Unable to successfully complete line search. Balking.\n");
      // if fc_selected is used, then replace the other fc's with new empirically determined ones
      if ((optinfo.selected_fc) && (optinfo.nonselected_fc == OPTInfo::EMPIRICAL))
        empirical_H(simples, symm, carts);
      open_PSIF();
      i = 1;
      psio_write_entry(PSIF_OPTKING, "Balked last time", (char *) &i, sizeof(int));
      close_PSIF();
      exit_io();
      exit(PSI_RETURN_BALK);
    }
    f_q = init_array(dim);
    open_PSIF(); // read old forces from previous step - should be same
    sprintf(force_string,"Internal Forces %d", optinfo.iteration-1);
    psio_read_entry(PSIF_OPTKING, force_string, (char *) f_q, dim* sizeof(double));
    sprintf(value_string,"RMS force %d", optinfo.iteration-1);
    psio_read_entry(PSIF_OPTKING, value_string, (char *) &rms_force, sizeof(double));
    sprintf(value_string,"MAX force %d", optinfo.iteration-1);
    psio_read_entry(PSIF_OPTKING, value_string, (char *) &max_force, sizeof(double));
    close_PSIF();
    // update internal coordinate data to new cartesan geometry
    coord = carts.get_coord();
    simples.compute(coord);
    simples.compute_s(coord);
    q = compute_q(simples,symm); 
    free_array(coord);
  }
  else {
  fprintf(outfile,"\n ** Taking geometry step number %d **\n", optinfo.iteration+1);

  open_PSIF();
  i=1;
  psio_write_entry(PSIF_OPTKING,"Consecutive line searches", (char *) &i, sizeof(int));
  close_PSIF();

  // compute forces in internal coordinates, f_q = G_inv B u f
  q = compute_q(simples,symm);
  B = compute_B(simples,symm); //print_mat(B, dim, dim_carts,outfile);
  G = compute_G(B,dim,carts);

  fprintf(outfile,"\nBuB^t ");
  G_inv = symm_matrix_invert(G,dim,1,optinfo.redundant);

  masses = carts.get_mass();
  u = mass_mat(masses);
  free_array(masses);

  f = carts.get_forces(); // in aJ/Ang //print_mat2(&f, 1, dim_carts, outfile);

  f_q = init_array(dim);
  temp_arr = init_array(dim);
  temp_arr2 = init_array(dim_carts);

  opt_mmult(u,0,&f,1,&temp_arr2,1,dim_carts,dim_carts,1,0);
  opt_mmult(B,0,&temp_arr2,1,&temp_arr,1,dim,dim_carts,1,0);
  opt_mmult(G_inv,0,&temp_arr,1,&f_q,1,dim,dim,1,0);

  // fprintf(outfile,"internal forces (G_inv B u f) in aJ/A\n");
  // print_mat2(&f_q, 1, dim, outfile);
  // for (i=0;i<dim;++i)
  //   f_q[i] = f_q[i] * sqrt(2) * _hartree2J * 1.0E18 / _bohr2angstroms;
  // fprintf(outfile, "%d %15.10lf\n", i, f_q[i] * _hartree2J / _bohr2angstroms * 1.0E18) ;

  // test by transforming f_q back to cartesian forces and compare
  // fprintf(outfile,"computed forces in cartesian coordinates aJ/Ang\n");
  // opt_mmult(B,1,&f_q,1,&temp_arr2,1,dim_carts,dim,1,0);
  // print_mat2(&temp_arr2, 1, dim_carts, outfile);

  free_matrix(B);
  free_array(f);
  free_array(temp_arr);
  free_array(temp_arr2);
  free_matrix(u);

  // Setup projection matrix P = G * G- for redundant internals
  P = init_matrix(dim,dim);
  opt_mmult(G,0,G_inv,0,P,0,dim,dim,dim,0); 
  free_matrix(G);
  free_matrix(G_inv);
  
  // Add constraints to projection matrix
    // C has 1's on diagonal for constraints and 0's elsewhere
  if (optinfo.constraints_present) {
    C = init_matrix(dim,dim);
    double i_overlap, j_overlap;
    for (b=0; b<optinfo.nconstraints; ++b) {
      constraint = simples.index_to_id(optinfo.constraints[b]);

      for (i=0; i<dim; ++i) {
        for (j=0; j<=i; ++j) {
          i_overlap = 0.0;
          for (k=0;k<symm.get_length(i);++k) {
            a = symm.get_simple(i,k);
            if ( a != constraint ) continue;
            i_overlap += symm.get_coeff(i,k) * symm.get_prefactor(i);
          }
          j_overlap = 0.0;
          for (k=0;k<symm.get_length(j);++k) {
            a = symm.get_simple(j,k);
            if ( a != constraint ) continue;
            j_overlap += symm.get_coeff(j,k) * symm.get_prefactor(j);
          }
          C[i][j] += ( i_overlap * j_overlap );
          C[j][i] = C[i][j];
        }
      }
    }

  // P = P' - P' C (CPC)^-1 C P'
    T = init_matrix(dim,dim);
    T2 = init_matrix(dim,dim);
    opt_mmult(P,0,C,0,T,0,dim,dim,dim,0); 
    opt_mmult(C,0,T,0,T2,0,dim,dim,dim,0); 
    T3 = symm_matrix_invert(T2, dim, 0, 1);

    opt_mmult(C,0,P,0,T,0,dim,dim,dim,0); 
    opt_mmult(T3,0,T,0,T2,0,dim,dim,dim,0); 
    opt_mmult(C,0,T2,0,T3,0,dim,dim,dim,0); 
    opt_mmult(P,0,T3,0,T2,0,dim,dim,dim,0); 
    for (i=0;i<dim;++i)
      for (j=0;j<dim;++j)
        P[i][j] -= T2[i][j];
    free_matrix(T);
    free_matrix(T2);
    free_matrix(T3);
    free_matrix(C);
  }

   // Project redundancies and constraints out of forces: f_q = P f_q'
   temp_arr = init_array(dim);
   opt_mmult(P,0,&f_q,1,&temp_arr,1,dim, dim,1,0);
   for(i=0; i<dim; ++i)
     f_q[i] = temp_arr[i];
   free_array(temp_arr);

  // Compute RMS and MAX forces
  rms_force = 0.0; 
  max_force = fabs(f_q[0]);
  for (i=0;i<dim;++i) {
    rms_force += SQR(f_q[i]);
    if (fabs(f_q[i]) > max_force) max_force = fabs(f_q[i]);
  } 
  rms_force = sqrt(rms_force/((double) dim));

  // Write data to PSIF_OPTKING for later use
  open_PSIF();
  sprintf(value_string,"Internal Values %d", optinfo.iteration);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) q, dim*sizeof(double));
  sprintf(force_string,"Internal Forces %d", optinfo.iteration);
  psio_write_entry(PSIF_OPTKING, force_string, (char *) f_q, dim* sizeof(double));
  coord = carts.get_coord();
  sprintf(value_string,"Cartesian Values %d", optinfo.iteration);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) coord, dim_carts*sizeof(double));
  free_array(coord);
  tval = carts.get_energy();
  sprintf(value_string,"Energy %d", optinfo.iteration);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) &tval, sizeof(double));
  sprintf(value_string,"RMS force %d", optinfo.iteration);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) &rms_force, sizeof(double));
  sprintf(value_string,"MAX force %d", optinfo.iteration);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) &max_force, sizeof(double));
  close_PSIF();

  printf("Energy: %15.10lf MAX force: %6.2e RMS force: %6.2e\n",carts.get_energy(),max_force, rms_force);

  fconst_init(carts, simples, symm); // makes sure some force constants are in PSIF_OPTKING

  H = compute_H(simples,symm,P,carts); // get updated and projected Hessian
  free_matrix(P);

  // *** standard Newton-Raphson search ***
  if (!optinfo.rfo) { // displacements from inverted, projected Hessian, H_inv f_q = dq
    dq = init_array(dim);
    H_inv = symm_matrix_invert(H,dim,0,0);
    opt_mmult(H_inv,0,&f_q,1,&dq,1,dim,dim,1,0);
    free_matrix(H_inv);

    step_limit(simples, symm, dq);
    check_zero_angles(simples, symm, dq);

    // get norm |x| and unit vector in the step direction
    dot_array(dq, dq, dim, &tval);
    nr_xnorm = sqrt(tval);

    nr_u = init_array(dim);
    for (i=0; i<dim; ++i)
      nr_u[i] = dq[i];
    normalize(&nr_u, 1, dim);
    
    // get gradient and hessian in step direction
    dot_array(f_q, nr_u, dim, &nr_g);
    nr_g *= -1; // gradient not force
    nr_h = 0;
    for (i=0; i<dim; ++i) {
      dot_array(H[i], nr_u, dim, &tval);
      nr_h += nr_u[i] * tval;
    }

    DE_new = nr_energy(nr_xnorm, nr_g, nr_h);
    fprintf(outfile,"Projected energy change by quadratic approximation: %20.10lf\n", DE_new);
    
    open_PSIF();
    sprintf(value_string,"DE prediction %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &DE_new, sizeof(double));
    sprintf(value_string,"Unit step U %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) nr_u, dim*sizeof(double));
    sprintf(value_string,"Xnorm %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &nr_xnorm, sizeof(double));
    sprintf(value_string,"Scalar gradient %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &nr_g, sizeof(double));
    sprintf(value_string,"Scalar hesian %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &nr_h, sizeof(double));
    close_PSIF();
  }
  else { // take RFO step

    // build and diagonalize RFO matrix
    rfo_mat = init_matrix(dim+1,dim+1);
    for (i=0; i<dim; ++i) {
      rfo_mat[dim][i] = - f_q[i];
      for (j=0; j<=i; ++j)
        rfo_mat[i][j] = H[i][j];
    }

    if(optinfo.print_hessian) {
      fprintf(outfile,"RFO mat\n");
      print_mat5(rfo_mat,dim+1,dim+1,outfile);
    }

    lambda = init_array(dim+1);
    work = init_array(3*(dim+1));
    C_DSYEV('V', 'U', dim+1 , rfo_mat[0], dim+1, lambda, work, 3*(dim+1));
    free_array(work);
    for (i=0; i<dim; ++i)
      lambda[i] /= (_hartree2J * 1.0e18); //aJ -> hartree

  fprintf(outfile,"RFO eigenvalues/lambdas\n");
  print_mat2(&lambda,1,dim+1,outfile);

  // I'm not sure about the intermediate normalization that is supposed to occur
  // scaling the evects to make the last element '1' and then "scaling of the eigenvector
  // by a factor alpha used to calculate mode displacements" (see page 54) amounts to
  // the same thing as just making sure the sign of the last entry of the evect is positive
  //  During the course of an optimization some evects may appear that are bogus leads (presumably)
  // due to a bad hessian - the root following can avoid them.  also some evects have a virtual
  // 0 in the last entry and so cannot be divided by it
  for (i=0; i<dim+1; ++i) {
    if (rfo_mat[i][dim] < 0)
      for (j=0;j<dim+1;++j) rfo_mat[i][j] *= -1;
  }

  //fprintf(outfile,"RFO evects\n");
  //print_mat5(rfo_mat,dim+1,dim+1,outfile);

    // *** choose which RFO eigenvector to use
    // if not root following, then use rfo_root'th lowest eigenvalue; default is 0 (lowest)
    if (!optinfo.rfo_follow_root) {
      rfo_root = optinfo.rfo_root;
    }
    else { // do root following
      open_PSIF();
      if (psio_tocscan(PSIF_OPTKING, "Previous RFO eigenvector") != NULL) {
        rfo_old_evect = init_array(dim+1);
        psio_read_entry(PSIF_OPTKING, "Previous RFO eigenvector",
          (char *) &(rfo_old_evect[0]), (dim+1)*sizeof(double));

        tval = 0;
        for (i=0; i<dim+1; ++i) {
          dot_array(rfo_mat[i],rfo_old_evect,dim+1,&tval2);
          if (tval2 > tval) {
            tval = tval2;
            rfo_root = i;
           }
        }
        free_array(rfo_old_evect);
        fprintf(outfile,"RFO vector %d has maximal overlap with previous step\n",rfo_root+1);
        psio_write_entry(PSIF_OPTKING, "Previous RFO eigenvector",
          (char *) rfo_mat[rfo_root], (dim+1)*sizeof(double));
      }
      else {
        psio_write_entry(PSIF_OPTKING, "Previous RFO eigenvector",
          (char *) &(rfo_mat[optinfo.rfo_root][0]), (dim+1)*sizeof(double));
        fprintf(outfile,"Using initial RFO vector %d to follow.\n",optinfo.rfo_root+1);
        rfo_root = optinfo.rfo_root;
      }
      close_PSIF();
    }

  /*for (i=0; i<dim+1; ++i) {
    if ((lambda[i] < 0.0) || (i <= rfo_root))
      if (fabs(rfo_mat[i][dim]) > 1.0e-4)
        for (j=0;j<dim+1;++j)
          rfo_mat[i][j] /= rfo_mat[i][dim];
  } */

    // print out lowest energy evects
    for (i=0; i<dim+1; ++i) {
      if ((lambda[i] < 0.0) || (i <rfo_root)) {
        fprintf(outfile,"RFO eigenvalue %d: %15.10lf (or 2*%-15.10lf)\n", i+1, lambda[i],lambda[i]/2);
        fprintf(outfile,"eigenvector:\n");
        print_mat2(&(rfo_mat[i]),1,dim+1,outfile);
      }
    }

/* alternative algorithm (H-lambdaI)x + g = 0
    dq = init_array(dim);
    double **H_test, **H_inv_test;
    H_test = init_matrix(dim,dim);
    for (i=0; i<dim; ++i)
      for (j=0; j<dim; ++j)
        H_test[i][j] = H[i][j];
    for (i=0; i<dim; ++i)
      H_test[i][i] = H[i][i] - lambda[0];
    H_inv_test = symm_matrix_invert(H_test,dim,0,0);
    opt_mmult(H_inv_test,0,&f_q,1,&dq,1,dim,dim,1,0);
    fprintf(outfile,"dq solved by (H- lambda I) x + g = 0\n");
    print_mat2(&dq,1,dim,outfile);
    dot_array(f_q, dq, dim, &tval);
    fprintf(outfile,"g^T x = lambda ? = %15.10lf\n", - tval);
    free_array(dq);
    free_matrix(H_test);
    free_matrix(H_inv_test);
*/

    dq = init_array(dim);
    for (j=0; j<dim; ++j)
      dq[j] = rfo_mat[rfo_root][j];
    rfo_eval = lambda[rfo_root];
    free_array(lambda);

    // scales dq
    step_limit(simples, symm, dq);
    check_zero_angles(simples, symm, dq);

    // get norm |x| and unit vector in the step direction
    dot_array(dq, dq, dim, &tval);
    rfo_xnorm = sqrt(tval);

    rfo_u = init_array(dim);
    for (j=0; j<dim; ++j)
      rfo_u[j] = rfo_mat[rfo_root][j];
    normalize(&rfo_u, 1, dim);

    free_matrix(rfo_mat);

    // get gradient and hessian in step direction
    dot_array(f_q, rfo_u, dim, &rfo_g);
    rfo_g *= -1; // gradient not force
    rfo_h = 0;
    for (i=0; i<dim; ++i) {
      dot_array(H[i], rfo_u, dim, &tval);
      rfo_h += rfo_u[i] * tval;
    }

    // for debugging
    DE_new = rfo_energy(rfo_xnorm, rfo_g, rfo_h);
    fprintf(outfile,"DE_RFO(rfo_xnorm = %8.3e) = %20.10lf\n", rfo_xnorm, DE_new);
    //for (rfo_t=0; rfo_t<2*rfo_xnorm; rfo_t += (rfo_xnorm/10.0))
    //  fprintf(outfile,"DE_RFO(rfo_t = %8.3e) = %20.10lf\n", rfo_t, rfo_energy(rfo_t, rfo_g, rfo_h));

    fprintf(outfile,"Projected energy change by RFO approximation: %20.10lf\n", DE_new);
    open_PSIF();
    sprintf(value_string,"DE prediction %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &DE_new, sizeof(double));
    sprintf(value_string,"Unit step U %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) rfo_u, dim*sizeof(double));
    sprintf(value_string,"Xnorm %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &rfo_xnorm, sizeof(double)); 
    sprintf(value_string,"Scalar gradient %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &rfo_g, sizeof(double)); 
    sprintf(value_string,"Scalar hesian %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &rfo_h, sizeof(double));
    close_PSIF();

    free_array(rfo_u);
  } // end take RFO step
  } // take forward geometry step


  fprintf(outfile, "\n\tInternal Coordinate Update in Ang or Rad, aJ/Ang or aJ/Rad\n");
  fprintf(outfile, "\t----------------------------------------------------------\n");
  fprintf(outfile, "\t       Value         Force        Displacement  New Value\n");
  fprintf(outfile, "\t----------------------------------------------------------\n");
  for (i=0;i<dim;++i)
    fprintf(outfile,"\t%3d %13.8lf %13.8lf %13.8lf %13.8lf\n", i+1, q[i], f_q[i], dq[i], q[i]+dq[i]);
  fprintf(outfile, "\t----------------------------------------------------------\n");

  rms_disp = 0.0; 
  max_disp = fabs(dq[0]);
  for (i=0;i<dim;++i) {
    rms_disp += SQR(dq[i]);
    if (fabs(dq[i]) > max_disp) max_disp = fabs(dq[i]);
  }
  rms_disp = sqrt(rms_disp/((double) dim));

  fprintf(outfile, "\n\tConvergence Check Cycle %4d:\n", optinfo.iteration+1);
  fprintf(outfile, "\t                    Actual        Tolerance     Converged?\n");
  fprintf(outfile, "\t----------------------------------------------------------\n");
  if ( fabs(optinfo.conv_max_force) < 1.0e-15 ) fprintf(outfile, "\tMAX Force        %10.1e\n", max_force);
  else fprintf(outfile, "\tMAX Force        %10.1e %14.1e %11s\n", max_force, optinfo.conv_max_force, 
       ((max_force < optinfo.conv_max_force) ? "yes" : "no"));
  if ( fabs(optinfo.conv_rms_force) < 1.0e-15 ) fprintf(outfile, "\tRMS Force        %10.1e\n", rms_force);
  else fprintf(outfile, "\tRMS Force        %10.1e %14.1e %11s\n", rms_force, optinfo.conv_rms_force, 
       ((rms_force < optinfo.conv_rms_force) ? "yes" : "no"));
  if ( fabs(optinfo.conv_max_DE) < 1.0e-15 ) fprintf(outfile, "\tEnergy Change    %10.1e\n", fabs(DE));
  else fprintf(outfile, "\tEnergy Change    %10.1e %14.1e %11s\n", DE, optinfo.conv_max_DE, 
       ((fabs(DE) < optinfo.conv_max_DE) ? "yes" : "no"));
  if ( fabs(optinfo.conv_max_disp) < 1.0e-15 ) fprintf(outfile, "\tMAX Displacement %10.1e\n", max_disp);
  else fprintf(outfile, "\tMAX Displacement %10.1e %14.1e %11s\n", max_disp, optinfo.conv_max_disp, 
       ((max_disp < optinfo.conv_max_disp) ? "yes" : "no"));
  if ( fabs(optinfo.conv_rms_disp) < 1.0e-15 ) fprintf(outfile, "\tRMS Displacement %10.1e\n", rms_disp);
  else fprintf(outfile, "\tRMS Displacement %10.1e %14.1e %11s\n", rms_disp, optinfo.conv_rms_disp, 
       ((rms_disp < optinfo.conv_rms_disp) ? "yes" : "no"));
  fprintf(outfile, "\t----------------------------------------------------------\n");

  // mightily complex convergence check  
  if (
         ( 
              ( (optinfo.opt_conv == OPTInfo::LOOSE) || (optinfo.opt_conv == OPTInfo::NORMAL) || 
                (optinfo.opt_conv == OPTInfo::TIGHT) || (optinfo.opt_conv == OPTInfo::VERY_TIGHT) || 
                (optinfo.opt_conv == OPTInfo::BAKER) || (optinfo.opt_conv == OPTInfo::QCHEM)              )
           && ( (max_force < optinfo.conv_max_force) && ((fabs(DE) < optinfo.conv_max_DE) || (max_disp < optinfo.conv_max_disp)) )
         )
      || (   
              ( (optinfo.opt_conv == OPTInfo::G03_NORMAL) || (optinfo.opt_conv == OPTInfo::G03_TIGHT) || 
                (optinfo.opt_conv == OPTInfo::G03_VERY_TIGHT)                                             )
           && ( ((max_force < optinfo.conv_max_force) && (rms_force < optinfo.conv_rms_force) &&
                 (max_disp  < optinfo.conv_max_disp)  && (rms_disp  < optinfo.conv_rms_disp)) ||
                (rms_force * 100 < optinfo.conv_rms_force) )
         )
      || ( 
              ( (optinfo.opt_conv == OPTInfo::GENERAL)                                                    )
           && ( ((fabs(optinfo.conv_max_force) < 1.0e-15) || (max_force < optinfo.conv_max_force)) &&
                ((fabs(optinfo.conv_rms_force) < 1.0e-15) || (rms_force < optinfo.conv_rms_force)) &&
                ((fabs(optinfo.conv_max_DE)    < 1.0e-15) || (fabs(DE)  < optinfo.conv_max_DE)) &&
                ((fabs(optinfo.conv_max_disp)  < 1.0e-15) || (max_disp  < optinfo.conv_max_disp)) &&
                ((fabs(optinfo.conv_rms_disp)  < 1.0e-15) || (rms_disp  < optinfo.conv_rms_disp)) )
         )
     ) {

    fprintf(outfile,"\nOptimization is complete.\n");

    ip_string("WFN", &(wfn),0);
    fprintf(outfile,"\nFinal %s energy is %15.10lf\n", wfn, carts.get_energy());
    free(wfn);

    //optinfo.iteration += 1;
    i = 0;
    open_PSIF();
    //psio_write_entry(PSIF_OPTKING, "Iteration", (char *) &(optinfo.iteration),sizeof(int));
    psio_write_entry(PSIF_OPTKING, "Iteration", (char *) &(i),sizeof(int));
    psio_write_entry(PSIF_OPTKING, "Balked last time", (char *) &(i), sizeof(int));
    close_PSIF();

    opt_report(outfile);

    fprintf(stdout,"\n  OPTKING:  optimization is complete\n");
    fprintf(outfile,"The Optimized geometry in a.u.\n");
    carts.print(12,outfile,0,disp_label,0);
    fprintf(outfile,"\nThe Optimized geometry in Angstrom\n");
    carts.print(13,outfile,0,disp_label,0);
    fprintf(outfile, "\n");

    if (optinfo.zmat) {
      int *unique_zvars;
      unique_zvars = init_int_array(MAX_ZVARS);
      //compute_zmat(carts, unique_zvars);
      //print_zmat(outfile, unique_zvars);
      free(unique_zvars);
      fprintf(outfile,"\n");
    }
    free(f_q); free(q);
    return(PSI_RETURN_ENDLOOP);
  } // end converged geometry

  if (optinfo.analytic_interfragment) { // do interfragment steps analytically

    for (isalc=0; isalc<dim; ++isalc) {
     
      // find any interfragment coordinate sets
      simple_id = symm.get_simple(isalc,0);
      simples.locate_id(simple_id, &intco_type, &sub_index, &sub_index2);
      // only execute for single (non-salc) interfragment coordinates
      if ( (symm.get_length(isalc) != 1) || (intco_type != FRAG) ) continue;
      // only execute this code once for each fragment pair
      // assume this one time is for the distance coordinate entry
      if (sub_index2 != 0) continue;
      // don't do analytic steps unless all reference points (and all coordinates) are present
      // because we may have to fix their values in orient_fragment even if we are not optmizing
      // with them
      if ( (simples.frag[sub_index].get_A_P() != 3) || (simples.frag[sub_index].get_B_P() != 3) )
        continue;

    coord = carts.get_coord();
    simples.compute(coord);
    // fix configuration for torsions; sets flag for torsions > FIX_NEAR_180 or < -FIX_NEAR_180 
    // so that when I calculate residual below, I get the right answer
    simples.fix_near_180();
    free_array(coord);

      for (cnt=0, I=0; I<6; ++I) {
        inter_q[I] = 0;
        if ( simples.frag[sub_index].get_coord_on(I) ) {
          inter_q[I] = dq[isalc + cnt];
          dq[isalc + cnt] = 0; // so that back-transformation code if used, doesn't change q again
          ++cnt;
        }
      }

      dot_array(inter_q,inter_q,6,&tval); // check if there are any non-zero displacements
      if (fabs(tval) < 1.0e-15) continue;

      for (I=0; I<6; ++I)
        inter_q[I] += simples.frag[sub_index].get_val_A_or_rad(I);

      if (optinfo.frag_dist_rho)
        inter_q[0]  = 1.0 / inter_q[0];
      inter_q[0]  /= _bohr2angstroms;
      for (I=1; I<6; ++I) inter_q[I] *= 180.0/_pi; 

      // get geometries of fragments
      geom = carts.get_coord_2d();

      a = simples.frag[sub_index].get_A_natom();
      geom_A = init_matrix(a,3);
      for (atom=0;atom<a;++atom)
        for (xyz=0;xyz<3;++xyz)
          geom_A[atom][xyz] = geom[simples.frag[sub_index].get_A_atom(atom)][xyz];

      b = simples.frag[sub_index].get_B_natom();
      geom_B = init_matrix(b,3);
      for (atom=0;atom<b;++atom)
        for (xyz=0;xyz<3;++xyz)
          geom_B[atom][xyz] = geom[simples.frag[sub_index].get_B_atom(atom)][xyz];

      // get reference point information for fragments
      weight_A = init_matrix(3,a);
      for (cnt=0; cnt<simples.frag[sub_index].get_A_P(); ++cnt)
        for (atom=0;atom<a;++atom)
          weight_A[cnt][atom] = simples.frag[sub_index].get_A_weight(cnt,atom);

      weight_B = init_matrix(3,b);
      for (cnt=0; cnt<simples.frag[sub_index].get_B_P(); ++cnt)
        for (atom=0;atom<b;++atom)
          weight_B[cnt][atom] = simples.frag[sub_index].get_B_weight(cnt,atom);

      fprintf(outfile,"\nAnalytically doing interfragment displacements\n");
      // move fragment B and put result back into main geometry matrix
      orient_fragment(a, b, simples.frag[sub_index].get_A_P(), simples.frag[sub_index].get_B_P(),
        geom_A, geom_B, weight_A, weight_B, inter_q[0], inter_q[1], inter_q[2], inter_q[3],
        inter_q[4], inter_q[5], outfile);

      for (atom=0; atom<b; ++atom)
        for (xyz=0; xyz<3; ++xyz)
          geom[simples.frag[sub_index].get_B_atom(atom)][xyz] = geom_B[atom][xyz];

      symmetrize_geom(geom[0]);
      carts.set_coord_2d(geom);
      free_matrix(geom);

      // leave any residual error in dq for back-transformation to clean up
      // put inter_q in rad/Ang and/or 1/R
      inter_q[0]  *= _bohr2angstroms;
      if (optinfo.frag_dist_rho)
        inter_q[0]  = 1.0 / inter_q[0]; // convert back to 1/R
      for (I=1; I<6; ++I) inter_q[I] *= _pi/180.0;

      for (cnt=0, I=0; I<6; ++I) {
        coord = carts.get_coord();
        simples.compute(coord);
        if ( simples.frag[sub_index].get_coord_on(I) ) {
          dq[isalc + cnt] = simples.frag[sub_index].get_val_A_or_rad(I) - inter_q[I];
//fprintf(outfile,"val: %20.15lf\n", simples.frag[sub_index].get_val_A_or_rad(I));
//fprintf(outfile,"inter_q[%d]: %20.15lf\n", I, inter_q[I]);
          ++cnt;
        }
        free_array(coord);
      }

      free_matrix(geom_A); free_matrix(geom_B);
      free_matrix(weight_A); free_matrix(weight_B);
    }
    // sprintf(disp_label,"\nNew Geometry in au after analytic fragment orientation\n");
    // carts.print(1,outfile,0,disp_label,0);
  } // if analytic_interfragment

//for (i=0; i<dim; ++i)
//fprintf(outfile, "dq[%d]: %15.10lf\n", i, dq[i]);

  // Do displacements with iterative backtransformation
  for (i=0;i<dim;++i)
    q[i] += dq[i];

  success = false;
  // check to see if all displacements are zero
  dot_array(dq,dq,dim,&tval);
  if (tval < 1.0e-15) {
    success = true;
    carts.print(PSIF_CHKPT,outfile,0,NULL,0);
  }

  if (!success) {
    /* new_geom will overwrite dq so send a copy */
    dq_to_new_geom = init_array(dim);
    for (i=0;i<dim;++i)
      dq_to_new_geom[i] = dq[i];

    strcpy(disp_label,"New Cartesian Geometry in a.u.");
    success = new_geom(carts,simples,symm,dq_to_new_geom,32,0,disp_label,0,0,djunk);
  }

  if (!success && do_line_search) {
    fprintf(outfile,"Giving up - unable to back-transform to new cartesian coordinates.\n");
    fprintf(outfile,"Returning balk code.\n");
    open_PSIF();
    i = 1;
    psio_write_entry(PSIF_OPTKING, "Balked last time", (char *) &i, sizeof(int));
    close_PSIF();
    exit_io();
    exit(PSI_RETURN_BALK);
  }
  
  int retry = 0;
  while (!success) {
    fprintf(outfile,"\tWarning: halving displacement size and reducing back-transformation \
 convergence criteria.\n");
    optinfo.bt_dx_conv *= 5.0;
    optinfo.bt_dq_conv *= 5.0;
    for (i=0;i<dim;++i) {
      dq[i] = dq[i] / 2.0;
      dq_to_new_geom[i] = dq[i];
    }
    success = new_geom(carts,simples,symm,dq_to_new_geom,32,0,disp_label,0,0,djunk);
    ++retry;
    if (!success && retry == 4) {
      fprintf(outfile,"Giving up - unable to back-transform to new cartesian coordinates.\n");
      open_PSIF();
      i = 1;
      psio_write_entry(PSIF_OPTKING, "Balked last time", (char *) &i, sizeof(int));
      close_PSIF();
      exit_io();
      exit(PSI_RETURN_BALK);
    }
  }

  // Modify predicted energe change if displacement size was reduced
  if (retry > 0) {
    dot_array(f_q, dq, dim, &tval);
    DE_old = -tval; // g^T x
    for (i=0; i<dim; ++i) {
      dot_array(H[i], dq, dim, &tval);
      DE_old += 0.5 * dq[i] * tval; // 1/2 x^T H x
    }
    DE_old /= _hartree2J * 1.0e18;
    fprintf(outfile,"\tNew projected energy change (minus any analytically taken steps\n");
    fprintf(outfile,"\talong interfragment coordinates: %15.10lf\n", DE_old);
    open_PSIF();
    sprintf(value_string,"DE prediction %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &(DE_old), sizeof(double));
    close_PSIF();
  }

  free_array(q); free_array(dq); free_array(f_q);
  if (!do_line_search) {
    free_matrix(H);
    optinfo.iteration += 1;
    open_PSIF();
    psio_write_entry(PSIF_OPTKING, "Iteration", (char *) &(optinfo.iteration),sizeof(int));
    i = 0;
    psio_write_entry(PSIF_OPTKING, "Balked last time", (char *) &(i), sizeof(int));
    close_PSIF();
  }
  delete [] djunk;
  delete [] disp_label;
  return(PSI_RETURN_SUCCESS);
}


bool line_search(cartesians &carts, int dim, double *dq) {
  int i, dim_carts;
  double *u, xnorm, g, h, M, a, b, c, tval, tval2, t, *old_x;
  char value_string[20];
  double E, E_old, DE, DE_old;

  fprintf(outfile,"\tTrying smaller step in same direction from previous geometry.\n");

  int consecutive = 1;
  open_PSIF();
  if (psio_tocscan(PSIF_OPTKING,"Consecutive line searches") != NULL)
    psio_read_entry(PSIF_OPTKING, "Consecutive line searches", (char *) &consecutive, sizeof(int));
  close_PSIF();

  fprintf(outfile,"\tConsecutive line search %d\n", consecutive);

  if (consecutive > optinfo.max_consecutive_line_searches) {
    fprintf(outfile,"\tConsecutive line searches(%d) exceeded max_consecutive_line_searches(%d)\n",
      consecutive, optinfo.max_consecutive_line_searches);
    return false;
  }
      
  E = carts.get_energy();
  dim_carts = 3*carts.get_natom();
      
  open_PSIF();
  sprintf(value_string,"Energy %d", optinfo.iteration-1);
  psio_read_entry(PSIF_OPTKING, value_string, (char *) &E_old, sizeof(double));
  DE = E - E_old;
  fprintf(outfile,"\t\tActual DE    : %20.10lf\n", DE);

  sprintf(value_string,"DE prediction %d", optinfo.iteration-1);
  psio_read_entry(PSIF_OPTKING, value_string, (char *) &DE_old, sizeof(double));
  fprintf(outfile,"\t\tPredicted DE : %20.10lf\n", DE_old);

  u = init_array(dim);
  sprintf(value_string,"Unit step U %d", optinfo.iteration-1);
  psio_read_entry(PSIF_OPTKING, value_string, (char *) u, dim*sizeof(double));
  //fprintf(outfile,"Unit step U %d\n", optinfo.iteration-1);
  //for (i=0; i<dim; ++i)
  //  fprintf(outfile,"%15.10lf\n", u[i]);

  sprintf(value_string,"Xnorm %d", optinfo.iteration-1);
  psio_read_entry(PSIF_OPTKING, value_string, (char *) &xnorm, sizeof(double));
  fprintf(outfile,"\t\tXnorm %d     : %20.10lf\n", optinfo.iteration-1, xnorm);

  sprintf(value_string,"Scalar gradient %d", optinfo.iteration-1);
  psio_read_entry(PSIF_OPTKING, value_string, (char *) &g, sizeof(double));
  fprintf(outfile,"\t\tScalar gradient %d: %15.10lf\n", optinfo.iteration-1, g);

  sprintf(value_string,"Scalar hesian %d", optinfo.iteration-1);
  psio_read_entry(PSIF_OPTKING, value_string, (char *) &h, sizeof(double));
  fprintf(outfile,"\t\tScalar hessian %d : %15.10lf\n", optinfo.iteration-1, h);
  close_PSIF();

  if (optinfo.rfo)
    fprintf(outfile,"\t\tCheck energy      : %20.15lf\n", rfo_energy(xnorm, g, h));
  else
    fprintf(outfile,"\t\tCheck energy      : %20.15lf\n", nr_energy(xnorm, g, h));
      
  // sign of third derivative is + if energy obtained was too HIGH, and - if energy
  // if energy was too low, which we will only worry about for TS searches
  M = 6 * (DE-DE_old) / (xnorm*xnorm *xnorm) * (_hartree2J*1.0e18); //aJ/Ang^3
  fprintf(outfile,"\t\tApproximate 3rd derivative, M : %15.10lf\n",M);
      
  bool success = false;
  // reduce step size in 5% increments and test magnitude of cubic term
  for (t = xnorm; t > 0.05*xnorm; t -= 0.10*xnorm) {
    //fprintf(outfile,"delta abs(g+0.5ht) = %10.5e\n", optinfo.step_energy_limit_back*fabs(g+0.5*h*t));
    //fprintf(outfile,"1/6 Mt^2 = %10.5e\n", M*t*t/6);
    if (M*t*t/6 < optinfo.step_energy_limit_back*fabs(g+0.5*h*t)) {
      success = true;
      break;
    }
  }
  // things do not look good.  Try 10% step.
  if (!success) t = 0.10 * xnorm;

  fprintf(outfile,"\tNew scalar with which to scale unit step, t = %15.10lf\n", t);
  fprintf(outfile,"\tSetting step to previous step reduced to %.0lf%% of original.\n\n", t/xnorm*100);

  if (optinfo.rfo) DE_old = rfo_energy(t, g, h);
  else DE_old = nr_energy(t, g, h);
  fprintf(outfile,"\tNew Projected DE %d: %20.10lf\n", optinfo.iteration-1,DE_old); 

  // if we are searching for a minimum; line search has failed; and projected energy change
  //  is positive; give up!
  if (!optinfo.ts && DE_old > 0.0) {
    fprintf(outfile,"\tProjected energy change is positive for minimum seach - giving up!\n");
    return false;
  }

  ++consecutive;
  open_PSIF();
    psio_write_entry(PSIF_OPTKING, "Consecutive line searches",
       (char *) &consecutive, sizeof(int));
  close_PSIF();

  // put last geometry into cartesian object
  open_PSIF();
  old_x = init_array(dim_carts);
  sprintf(value_string,"Cartesian Values %d", optinfo.iteration-1);
  psio_read_entry(PSIF_OPTKING, value_string, (char *) old_x, dim_carts*sizeof(double));
  carts.set_coord(old_x);
  fprintf(outfile,"Setting cartesian coordinates to prevous step.\n");
  carts.print(5, outfile, 0, NULL, 0);
  free_array(old_x);

  // reduced step is dq = u * rho_t;
  for (i=0; i<dim; ++i)
    dq[i] = t * u[i]; 

  // replace old energy, gradient, hessian, unit step, unchanged
  sprintf(value_string,"DE prediction %d", optinfo.iteration-1);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) &DE_old, sizeof(double));
   
  sprintf(value_string,"Xnorm %d", optinfo.iteration-1);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) &t, sizeof(double));

  close_PSIF();
  fflush(outfile);
  return true;
}

}} /* namespace psi::optking */

