/*! \file
    \ingroup OPTKING
    \brief freq_grad_cart(): computes frequencies from gradients and cartesian disps
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
#include <libchkpt/chkpt.h>

namespace psi { //namespace optking {

void freq_grad_cart(const cartesians &carts) {
  int i,j,k,l,a,b, ii, cnt, dim, natom, xyz, cnt_eval = -1, *evals_all_irrep,op_disp;
  int xyzA, xyzB, atomA, atomB, row, col,h,nirreps,xyz_irr,I,cnt_all,match;
  int nsalcs,ncoord,start_disp, start_salc, atom, atom2, op, natom_unique;
  double **B, *masses,  **force_constants, cm_convert, k_convert;
  double *f, *f_q, *temp_arr, **evects, *evals, **force_constants_x;
  double *micro_e, **disp_grad, *grad, **grads_adapted;
  double *evals_all, **cartrep, **disp_grad_all, **normal;
  int *nsalc, *ndisp, ndisp_all, nsalc_all, **ict, *start_irr, print;
  char *line1;
  print = optinfo.print_cartesians;

  nirreps = syminfo.nirreps;
  natom = carts.get_natom();
  masses = carts.get_mass();
  ndisp = init_int_array(nirreps);
  nsalc = init_int_array(nirreps);

  chkpt_init(PSIO_OPEN_OLD);
  cartrep = chkpt_rd_cartrep();
  ict = chkpt_rd_ict();
  natom_unique = chkpt_rd_num_unique_atom();
  chkpt_close();

  /* Read in data from PSIF_OPTKING */
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Num. of disp.",
      (char *) &(ndisp_all), sizeof(int));
  psio_read_entry(PSIF_OPTKING, "OPT: Num. of coord.",
      (char *) &(nsalc_all), sizeof(int));
  psio_read_entry(PSIF_OPTKING, "OPT: Disp. per irrep",
      (char *) &(ndisp[0]), nirreps*sizeof(int));
  psio_read_entry(PSIF_OPTKING, "OPT: Coord. per irrep",
      (char *) &(nsalc[0]), nirreps*sizeof(int));

  fprintf(outfile,"\n total nsalc: %d, total ndisp: %d\n", nsalc_all, ndisp_all);
  fprintf(outfile,"nsalc per irrep: "); for (h=0; h<nirreps; ++h) fprintf(outfile,"%d ",nsalc[h]);
  fprintf(outfile,"\n");
  fprintf(outfile,"ndisp per irrep: "); for (h=0; h<nirreps; ++h) fprintf(outfile,"%d ",ndisp[h]);
  fprintf(outfile,"\n\n");

  disp_grad = block_matrix(ndisp_all,3*natom);
  if (optinfo.grad_dat) { /* read in gradients from "file11.dat" */
    opt_ffile(&fp_11, "file11.dat", 2);
    rewind(fp_11);
    line1 = new char[MAX_LINELENGTH+1];
    for (i=0;i<ndisp_all;++i) {
      /* read in 2 header lines */
      fgets(line1, MAX_LINELENGTH, fp_11);
      fgets(line1, MAX_LINELENGTH, fp_11);

      for (j=0; j<natom; ++j) {
        fgets(line1, MAX_LINELENGTH, fp_11);
      }
      /* read in xyz for N atoms */
      for (j=0; j<natom; ++j) {
        fgets(line1, MAX_LINELENGTH, fp_11);
        sscanf(line1, "%lf %lf %lf",
           &(disp_grad[i][3*j]), &(disp_grad[i][3*j+1]), &(disp_grad[i][3*j+2]) );
      }
    }
    fclose(fp_11);
  }
  else {
    psio_read_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(disp_grad[0][0]), ndisp_all*3*natom*sizeof(double));
  }
  if (print) {
    fprintf(outfile,"Gradients of displaced geometries\n");
    print_mat(disp_grad,ndisp_all,3*natom,outfile);
  }

  /**** construct full gradient list using unique ones ****/
  if (optinfo.points == 3) 
    disp_grad_all = block_matrix(2*nsalc_all, 3*natom);
  else if (optinfo.points == 5) 
    disp_grad_all = block_matrix(4*nsalc_all, 3*natom);
  /* just copy symmetric gradients */
  for (j=0; j<ndisp[0]; ++j) {
    for (I=0; I<3*natom; ++I)
      disp_grad_all[j][I] = disp_grad[j][I];
  }

  /* add gradients to list for asymmetric displacements */
  cnt_all = cnt = ndisp[0] - 1;
  for (h=1; h<nirreps; ++h) {

    if (!ndisp[h]) continue;

    /* copy computed displacement in */
    for (j=0; j<ndisp[h]; ++j) {
      ++cnt;
      ++cnt_all;
      for (I=0; I<3*natom; ++I)
        disp_grad_all[cnt_all][I] = disp_grad[cnt][I];

    /* determine which operation takes minus displacement into plus displacement */
      for (op_disp=0; op_disp<syminfo.nirreps; ++op_disp) {
        if ( syminfo.ct[h][op_disp] == -1 )
          break;
      }
      fprintf(outfile,"\tOperation that takes plus displacement %d to minus is %s.\n",
          cnt+1, syminfo.op_lbls[op_disp]);

      ++cnt_all;
      for (atom=0; atom<natom; ++atom) {
        for (xyz=0; xyz<3; ++xyz) {
          atom2 = ict[op_disp][atom] - 1;
          disp_grad_all[cnt_all][3*atom2+xyz] = disp_grad[cnt][3*atom+xyz] * cartrep[op_disp][3*xyz+xyz];
        }
      }
    }
    ndisp[h] += ndisp[h];
  }
  free_block(disp_grad);

  /* fix number of displacements to match full list of gradients */
  ndisp_all=0;
  for (h=0; h<nirreps; ++h) {
    ndisp_all += ndisp[h];
  }

  B = block_matrix(nsalc_all,3*natom);
  psio_read_entry(PSIF_OPTKING, "OPT: Adapted cartesians",
    (char *) &(B[0][0]), nsalc_all*3*natom*sizeof(double));
  if (print) {
    fprintf(outfile,"B matrix\n");
    print_mat(B,nsalc_all,3*natom,outfile);
  }

  evals_all = init_array(nsalc_all);
  evals_all_irrep = init_int_array(nsalc_all);

  /* just for printing out force constants */
  start_irr = init_int_array(nirreps);
  for (h=0; h<nirreps; ++h) {
    for (i=0;i<h;++i)
      start_irr[h] += nsalc[i];
  }

  /**** loop over irreps of salcs and vibrations ****/
  start_salc = 0;
  start_disp = 0;
  for (h=0; h<nirreps; ++h) {

    if (!nsalc[h]) continue;

    /* mass-weight gradients; g_xm = (1/sqrt(m)) * g_x */
    for (k=0; k<ndisp[h]; ++k)
      for (j=0;j<3*natom;++j)
        disp_grad_all[k+start_disp][j] /= sqrt(masses[j]);

    if (print) {
      fprintf(outfile,"Displaced gradients/sqrt(masses[j])\n");
      print_mat(&(disp_grad_all[start_disp]),ndisp[h],3*natom,outfile);
      fflush(outfile);
    }

    // compute forces in internal coordinates, f_q = G_inv B u f_x
   // In this case, B = c * masses^(1/2).  =>  G=I.
   // Thus, f_q = c * f_x / sqrt(masses)
    f_q = init_array(nsalc[h]);
    grads_adapted = block_matrix(ndisp[h],3*natom);
  
    if (print) fprintf(outfile,"Gradients recomputed from internal coordinate gradients\n");
    for (k=0; k<ndisp[h]; ++k) {

      opt_mmult(&(B[start_salc]),0,&(disp_grad_all[k+start_disp]),1,&f_q,1,nsalc[h],3*natom,1,0);

      if (print) fprintf(outfile,"ndisp[h]: %d, start_disp: %d\n", ndisp[h], start_disp);
  
      for (j=0;j<3*natom;++j)
        grads_adapted[k][j] = f_q[j];

      /* test by transforming f_q back to cartesian forces and compare */
      if (print) {
        temp_arr = init_array(3*natom);
        opt_mmult(B,1,&(grads_adapted[k]),1,&temp_arr,1,3*natom,nsalc[h],1,0);
        print_mat2(&temp_arr, 1, 3*natom, outfile);
        free_array(temp_arr);
      }
    }
    if (print) {
      fprintf(outfile,"Adapted gradients\n");
      print_mat(grads_adapted,ndisp[h],3*natom,outfile);
    }

    free_array(f_q);

    force_constants = block_matrix(nsalc[h],nsalc[h]);

    /*** Construct force constant matrix from finite differences of forces ***/
    fprintf(outfile,"\n ** Using %d-point formula.\n",optinfo.points);
    if (optinfo.points == 3) {
      for (i=0; i<nsalc[h]; ++i) {
        for (j=0; j<nsalc[h]; ++j) {
          /* Fij = fj(i+1) - fj(i-1) / (2h) */
          force_constants[i][j] = 
            (grads_adapted[2*i+1][j] - grads_adapted[2*i][j])
              / (2.0 * optinfo.disp_size);
        }
      }
    }
    else if (optinfo.points == 5) {
      for (i=0; i<nsalc[h]; ++i) {
        for (j=0; j<nsalc[h]; ++j) {
            /* fj(i-2) - fj(i+2) - 8fj(i-1) + 8fj(i+1)  / (12h) */
          force_constants[i][j] = 
            (  1.0 * grads_adapted[4*i][j]   - 1.0 * grads_adapted[4*i+1][j]
             - 8.0 * grads_adapted[4*i+2][j] + 8.0 * grads_adapted[4*i+3][j] )
             / (12.0 * optinfo.disp_size);
        }
      }
    }

    fprintf(outfile,"Force constants for irrep %s in symmetry-adapted cartesian coordinates\n",
      syminfo.clean_irrep_lbls[h]);
    print_mat(force_constants, nsalc[h], nsalc[h], outfile);
    fflush(outfile);

    start_salc += nsalc[h];
    start_disp += ndisp[h];

    dim = nsalc[h];

    /** Find eigenvalues of force constant matrix **/
    evals  = init_array(dim);
    evects = block_matrix(dim, dim);
    dgeev_optking(dim, force_constants, evals, evects);
    free_block(force_constants);
    sort(evals, evects, dim);

	normal = block_matrix(3*natom, dim);
    opt_mmult(&(B[start_irr[h]]),1,evects,0,normal,0,3*natom,dim,dim,0);
	fprintf(outfile,"\n\tNormal coordinates for irrep %s\n",syminfo.clean_irrep_lbls[h]);
    print_evects(normal, evals, 3*natom, dim, outfile);
	free_block(normal);
    free_block(evects);

    for (i=0; i<dim; ++i) {
      ++cnt_eval;
      evals_all[cnt_eval] = evals[i];
      evals_all_irrep[cnt_eval] = h;
    } 
    free_array(evals);
  }
  free_array(masses);
  free_block(disp_grad_all);

  sort_evals_all(nsalc_all,evals_all, evals_all_irrep);

  /* convert evals from H/(kg bohr^2) to J/(kg m^2) = 1/s^2 */
  /* v = 1/(2 pi c) sqrt( eval ) */
  fprintf(outfile, "\n\t Harmonic Vibrational Frequencies in cm^(-1)\n");
  fprintf(outfile,   "\t--------------------------------------------\n");
  k_convert = _hartree2J/(_bohr2m * _bohr2m * _amu2kg);
  cm_convert = 1.0/(2.0 * _pi * _c * 100.0);
  for(i=nsalc_all-1; i>-1; --i) {
    if(evals_all[i] < 0.0)
      fprintf(outfile, "\t %5s %10.3fi\n", syminfo.irrep_lbls[evals_all_irrep[i]],
          cm_convert * sqrt(-k_convert * evals_all[i]));
    else
      fprintf(outfile, "\t %5s %10.3f\n", syminfo.irrep_lbls[evals_all_irrep[i]],
          cm_convert * sqrt(k_convert * evals_all[i]));
    }
  fprintf(outfile,   "\t----------------------------\n");
  fflush(outfile);
  free_array(evals_all);
  free_int_array(evals_all_irrep);

  optinfo.disp_num = 0;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num),sizeof(int));
  close_PSIF();
  return;
}

void sort_evals_all(int nsalc_all, double *evals_all, int *evals_all_irrep) {
  int i,j,tmp_irrep,min_index;
  double min,tmp_eval;

  for (j=0; j<nsalc_all; ++j) {

    min = 1.0E6;
    for (i=j; i<nsalc_all; ++i) {
      if (evals_all[i] < min) {
        min = evals_all[i];
        min_index = i;
      }
    }
    tmp_eval = evals_all[j];
    evals_all[j]=evals_all[min_index];
    evals_all[min_index] = tmp_eval;

    tmp_irrep = evals_all_irrep[j];
    evals_all_irrep[j] = evals_all_irrep[min_index];
    evals_all_irrep[min_index] = tmp_irrep;
  }
  return;
}

}//} /* namespace psi::optking */

