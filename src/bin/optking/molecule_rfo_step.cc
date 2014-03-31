/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \file molecule.cc
    \ingroup optking
    \brief molecule class (really, molecular system class)
*/

#include "molecule.h"

#include <iostream>
#include <sstream>

#include "linear_algebra.h"
#include "print.h"
#include "atom_data.h"
#include "physconst.h"

#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif

namespace opt {

// compute change in energy according to RFO approximation
inline double DE_rfo_energy(double rfo_t, double rfo_g, double rfo_h) {
  return (rfo_t * rfo_g + 0.5 * rfo_t * rfo_t * rfo_h)/(1 + rfo_t*rfo_t);
}

// Take Rational Function Optimization step
void MOLECULE::rfo_step(void) {
  int i, j;
  int dim = g_nintco();
  int natom = g_natom();
  double tval, tval2;
  double *fq = p_Opt_data->g_forces_pointer();
  double **H = p_Opt_data->g_H_pointer();
  double *dq = p_Opt_data->g_dq_pointer();
  double **SRFO = init_matrix(dim+1,dim+1);
  double **rfo = init_matrix(dim+1,dim+1);

  // build (lower-triangle of) RFO matrix and diagonalize
//  double **rfo_mat = init_matrix(dim+1, dim+1);
  for (i=0; i<dim; ++i){
    for (j=0; j<=i; ++j){
//      rfo_mat[i][j] = H[i][j];
      rfo[i][j] = H[i][j];
    }
  }
//  matrix_copy(rfo_mat,rfo,dim+1,dim+1);
  for (i=0; i<dim; ++i) {
//    rfo_mat[dim][i] = - fq[i];
      rfo[i][dim]= rfo[dim][i] = -fq[i]; 
  }
//  if (Opt_params.print_lvl >= 3) {
//    fprintf(outfile,"RFO mat\n");
//    print_matrix(outfile, rfo_mat, dim+1, dim+1);
//  }
  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"RFO mat\n");
    print_matrix(outfile, rfo, dim+1, dim+1);
  }
  double **Hevects = matrix_return_copy(H,dim,dim); 
  double *h = init_array(dim);
  opt_symm_matrix_eig(Hevects, dim, h);
//  double *evalsRfo = init_array(dim+1);
//  opt_symm_matrix_eig(rfo_mat, dim+1, evalsRfo);
//  if (Opt_params.print_lvl >= 2) {
//    fprintf(outfile,"RFO eigenvalues/evalsRfos\n");
//    print_matrix(outfile, &(evalsRfo), 1, dim+1);
//  }
  if (Opt_params.print_lvl >= 2) {
//    fprintf(outfile,"RFO eigenvalues/evalsRfos\n");
//    print_matrix(outfile, &(h), 1, dim+1);
  }


  int rfo_root, f;
  double rfo_eval;
  double *evalsSrfo = init_array(dim+1);
  double *rfo_u;      // unit vector in step direction
  double dqtdq;
  double rfo_dqnorm;   // norm of step
  double rfo_g;         // gradient in step direction
  double rfo_h;         // hessian in step direction
  double DE_projected; // projected energy change by quadratic approximation
  double alpha; //scaling factor
  double** smat = init_matrix(dim+1,dim+1); //scaling matrix
  bool symm_rfo_step = false;
  double * rfo_old_evect = p_Opt_data->g_rfo_eigenvector_pointer();
  double lambda;
  double sum;
  double analyticDerivative;
  double trust = Opt_params.intrafragment_step_limit;
  dqtdq = 10;
  alpha = 1;

//Iterative sequence to find s
  while ( fabs( trust - sqrt(dqtdq )) > 1e-5) {

	//Creates scaling matrix
	  for(i=0; i<=dim; i++) {
		for(j=0; j<=dim; j++) {
			if(i==dim && j==dim) {
				smat[i][j]=1;
			}
			else if(i==j) {
				smat[i][j]=1/alpha;
			}
			else {
				smat[i][j]=0;
			}
		}
	  }
           if (Opt_params.print_lvl >= 3) {
              fprintf(outfile,"Scaling matrix\n");
   	      print_matrix(outfile, smat, dim+1, dim+1);
  	   }


 	 //Scales RFO matrix
	  opt_matrix_mult(smat, 0, rfo, 0, SRFO, 0, dim+1, dim+1, dim+1, 0);
          if (Opt_params.print_lvl >= 3) {
              fprintf(outfile,"Scaled matrix\n");
              print_matrix(outfile, SRFO, dim+1, dim+1);
           }
	  //Finds eigenvectors and eigenvalues
	  opt_symm_matrix_eig(SRFO, dim+1, evalsSrfo);
          if (Opt_params.print_lvl >= 3) {
              fprintf(outfile,"eigenvectors of scaled rfo\n");
              print_matrix(outfile, SRFO, dim+1,dim+1);
           }

	  // Do intermediate normalization.  
	  // RFO paper says to scale eigenvector to make the last element equal to 1.
	  // During the course of an optimization some evects may appear that are bogus leads
	  // - the root following can avoid them. 
	  for (i=0; i<dim+1; ++i) {
	    tval = SRFO[i][dim];
	    if (fabs(tval) > Opt_params.rfo_normalization_min) {
	      for (j=0;j<dim+1;++j)
	        SRFO[i][j] /= SRFO[i][dim];
	    }
	  }
	  if (Opt_params.print_lvl >= 3) {
	    fprintf(outfile,"scaled RFO eigenvectors (rows)\n");
	    print_matrix(outfile, SRFO, dim+1, dim+1);
	  }

	  // *** choose which RFO eigenvector to use
	  // if not root following, then use rfo_root'th lowest eigenvalue; default is 0 (lowest)
	  if ( (!Opt_params.rfo_follow_root) || (p_Opt_data->g_iteration() == 1)) {
	    rfo_root = Opt_params.rfo_root;
	    fprintf(outfile,"\tGoing to follow RFO solution %d.\n", rfo_root+1);
	
	    // Now test RFO eigenvector and make sure that it is totally symmetric.
	
	    while (!symm_rfo_step) {

	      symm_rfo_step = intco_combo_is_symmetric(SRFO[rfo_root], dim);

	      if (!symm_rfo_step) {
        	fprintf(outfile,"\tRejecting RFO root %d because it breaks the molecular point group.\n", rfo_root+1);
	        fprintf(outfile,"\tIf you are doing an energy minimization, there may exist a lower-energy, ");
        	fprintf(outfile,"structure with less symmetry.\n");
	        ++rfo_root;
	      }

	      if (rfo_root == dim+1) // quit in the unlikely event we've checked them all
	        break;
	    }
	  }
	  else { // do root following
	    tval = 0;
	    for (i=0; i<dim; ++i) { // dot only within H block, excluding rows and columns approaching 0
	      tval2 = array_dot(SRFO[i], rfo_old_evect,dim);
	      if (tval2 > tval) {
        	tval = tval2;
	        rfo_root = i;
	      }
	    }
	    fprintf(outfile,"RFO vector %d has maximal overlap with previous step\n", rfo_root+1);
	  }
	  p_Opt_data->set_rfo_eigenvector(SRFO[rfo_root]);

	  // print out lowest energy evects
	  if (Opt_params.print_lvl >= 3) {
	    for (i=0; i<dim+1; ++i) {
	      if ((evalsSrfo[i] < 0.0) || (i <rfo_root)) {
	        fprintf(outfile,"Scaled RFO eigenvalue %d: %15.10lf (or 2*%-15.10lf)\n", i+1, evalsSrfo[i],evalsSrfo[i]/2);
	        fprintf(outfile,"eigenvector:\n");
	        print_matrix(outfile, &(SRFO[i]), 1, dim+1);
	      }
	    }
	  }
	  free_array(evalsSrfo);

	  for (j=0; j<dim; ++j) {
	    dq[j] = SRFO[rfo_root][j]; // leave out last column
          }
 
	//DR.KING, WE DON'T REALLY KNOW WHAT THIS FROZEN FRAGMENT IS FOR AND IF IT NEEDS ALTERING
	  // Zero steps for frozen fragment  
//	  for (f=0; f<fragments.size(); ++f) {
//	    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
//	      fprintf(outfile,"\tZero'ing out displacements for frozen fragment %d\n", f+1);
//	      for (i=0; i<fragments[f]->g_nintco(); ++i)
//	        dq[ g_intco_offset(f) + i ] = 0.0;
//	    }
//	  }

//	  apply_intrafragment_step_limit(dq);
	  //check_intrafragment_zero_angles(dq);
	
	  // get norm |dq| and unit vector in the step direction
	  dqtdq = array_dot(dq, dq, dim);
	  rfo_dqnorm = sqrt( dqtdq );
	  rfo_u = init_array(dim);
	  array_copy(SRFO[rfo_root], rfo_u, dim);
	  array_normalize(rfo_u, dim+1);
	
	  // find the analytical derivative
	  lambda = array_dot(fq,dq,dim);
          if (Opt_params.print_lvl >= 2) {
            fprintf(outfile,"lambda");
            fprintf(outfile, "%10.10f",lambda);
   	    fprintf(outfile, "\n");
          }
	  sum = 0;
	  for (i=0; i<dim; i++) {
		sum = sum + ( pow(array_dot( Hevects[i], fq, dim),2) ) / ( pow(( h[i]-lambda*alpha ),3) );
	  }
	  analyticDerivative = 2*lambda / (1+alpha*dqtdq ) * sum;
	  alpha = alpha + 2*(trust*rfo_dqnorm - dqtdq) / analyticDerivative;
          if (Opt_params.print_lvl >= 2) {
            fprintf(outfile,"alpha is");
            fprintf(outfile, "%10.10f",alpha);
          }
   }
//  free_matrix(rfo_mat);
  free_matrix(SRFO);
 // get gradient and hessian in step direction
  rfo_g = -1 * array_dot(fq, rfo_u, dim);
  rfo_h = 0;
  for (i=0; i<dim; ++i)
    rfo_h += rfo_u[i] * array_dot(Hevects[i], rfo_u, dim);

  DE_projected = DE_rfo_energy(rfo_dqnorm, rfo_g, rfo_h);
  fprintf(outfile,"\tProjected energy change by RFO approximation: %20.10lf\n", DE_projected);
 free_matrix(Hevects);
/* Test step sizes
  double *x_before = g_geom_array();
  double **G = compute_G(true);
  double **G_inv = symm_matrix_inv(G, dim, 1);
  free_matrix(G);
  double *Gdq = init_array(dim);
  opt_matrix_mult(G_inv, 0, &dq, 1, &Gdq, 1, dim, dim, 1, 0);
  free_matrix(G_inv);
  double N = array_dot(dq, Gdq, dim);
  N = sqrt(N);
  free_array(Gdq);
  fprintf(outfile,"Step-size in mass-weighted internal coordinates: %20.10lf\n", N);
*/


// do displacements for each fragment separately
  for (f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      fprintf(outfile,"\tDisplacements for frozen fragment %d skipped.\n", f+1);
      continue;
    }
    fragments[f]->displace(&(dq[g_intco_offset(f)]), &(fq[g_intco_offset(f)]), g_atom_offset(f));
  }

  // do displacements for interfragment coordinates
  for (int I=0; I<interfragments.size(); ++I) {
    if (interfragments[I]->is_frozen() || Opt_params.freeze_interfragment) {
      fprintf(outfile,"\tDisplacements for frozen interfragment %d skipped.\n", I+1);
      continue;
    }
    interfragments[I]->orient_fragment( &(dq[g_interfragment_intco_offset(I)]),
                                        &(fq[g_interfragment_intco_offset(I)]) );
  }

#if defined(OPTKING_PACKAGE_QCHEM)
  // fix rotation matrix for rotations in QCHEM EFP code
  for (int I=0; I<efp_fragments.size(); ++I)
    efp_fragments[I]->displace( I, &(dq[g_efp_fragment_intco_offset(I)]) );
#endif

  symmetrize_geom(); // now symmetrize the geometry for next step
  
/* Test step sizes
  double *x_after = g_geom_array();
  double *masses = g_masses();

  double sum = 0.0;
  for (int i=0; i<g_natom(); ++i)
    for (int xyz=0; xyz<3; ++xyz)
      sum += (x_before[3*i+xyz] - x_after[3*i+xyz]) * (x_before[3*i+xyz] - x_after[3*i+xyz])
               * masses[i];

  sum = sqrt(sum);
  fprintf(outfile,"Step-size in mass-weighted cartesian coordinates [bohr (amu)^1/2] : %20.10lf\n", sum);
  free_array(x_after);
  free_array(x_before);
  free_array(masses);
*/

  // save values in step data
  p_Opt_data->save_step_info(DE_projected, rfo_u, rfo_dqnorm, rfo_g, rfo_h);

  free_array(rfo_u);
  fflush(outfile);

} // end take RFO step

}

