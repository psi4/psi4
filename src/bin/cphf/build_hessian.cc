/*! \file
    \ingroup CPHF
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

/* build_hessian(): Compute the vibrational (Cartesian) hessian using
** the skeleton second derivatives, the derivative Fock integrals, and
** the derivative overlap integrals computed by "cints --deriv2", as
** well as the nuclear-perturbation CPHF coefficients computed in
** cphf_X().  The equation for the electronic contribution to the
** Hessian is (MO basis):
**
** d E/da db = E^ab - 2 (S^ab)_ii eps_i - 2 (eta^ab)_ii eps_i
**     + 4 [ (U^b)_pj (F^a)_pj + (U^a)_pj (F^b)_pj ]
**     + 4 (U^a)_pj (U^b)_pj eps_i
**     + 4 (U^a)_pj (U^b)_ql [ 4(pj|ql) - (pq|jl) - (pl|jq) ]
**
** where E^ab is the energy expression evaluated using
** second-derivative integrals, S^ab is the overlap second derivative
** integrals, eps_p is the energy of orbital p, U^a and U^b are CPHF
** coefficients, F^a and F^b are derivative Fock integrals, and eta^ab
** is defined as
**
** (eta^ab)_pq = [(U^a)_pr (U^b)_qr + (U^a)_pr (U^b)_qr 
**                    - (S^a)_pr (S^b)_qr - (S^a)_pr (S^b)_qr ]
**
** For more details (and clearer notation) see:
**
** Y. Yamaguchi et al., "A New Dimension to Quantum Chemistry:
** Analytic Derivative Methods in Ab Initio Molecular Electronic
** Structure Theory", Oxford Press, New York, 1994.  Ch.4,
** esp. pp. 60-63.
**
** TDC, December 2001 (revised October 2002)
*/

namespace psi { namespace cphf {

void build_hessian(double ***F, double ***S, 
		   double **A, double ***U, 
		   double **hessian)
{
  int coord, coord_a, coord_b;
  int i, isym, ifirst, ilast;
  int m;
  int p, psym, pfirst, plast;
  int j, jsym, jfirst, jlast;
  int q, qsym, qfirst, qlast;
  int l, lsym, lfirst, llast;
  int pj, ql ;
  double *eta;
  FILE *file15;
  
  /* grab skeleton derivatives from cints --deriv2 **/
  psio_open(PSIF_DERINFO, PSIO_OPEN_OLD);
  psio_read_entry(PSIF_DERINFO, "Skeleton Hessian", 
		 (char *) hessian[0], natom*3*natom*3*sizeof(double));
  psio_close(PSIF_DERINFO, 1);

  if (print_lvl > 5)
  {
    fprintf(outfile,"Skeleton Hessian\n");
    print_mat(hessian, natom*3, natom*3, outfile);
  }
  
   
  eta = init_array(nmo);
  for(coord_a=0; coord_a < natom*3; coord_a++) {
    for(coord_b=0; coord_b < natom*3; coord_b++) {

      for(isym=0; isym < nirreps; isym++) {
	
	ifirst = ofirst[isym];
	ilast = olast[isym];
	
	for(i=ifirst; i <= ilast; i++) {
	
          eta[i] = 0.0;
		 
	  for(m=0; m < nmo; m++) {
	  
            eta[i] += (2.0 * U[coord_a][i][m] * U[coord_b][i][m] - 
		       2.0 * S[coord_a][i][m] * S[coord_b][i][m]);
	  }
	}
      }
      
      for(isym=0; isym < nirreps; isym++) {
	      
	ifirst = ofirst[isym];
	ilast = olast[isym];
	
	for(i=ifirst; i <= ilast; i++) {
	  
          hessian[coord_a][coord_b] -= 2.0 * eta[i] * evals[i];
	}
      }

      for(psym=0; psym < nirreps; psym++) {

	pfirst = first[psym];
	plast = last[psym];

	for(p=pfirst; p <= plast; p++) {

	  for(jsym=0; jsym < nirreps; jsym++) {
	 
            jfirst = ofirst[jsym];
	    jlast = olast[jsym];

	    for(j=jfirst; j <= jlast; j++) {

              pj = INDEX(p,j);
	      
	      hessian[coord_a][coord_b] += 4.0 * (U[coord_b][p][j] *
			                          F[coord_a][p][j] + 
						  U[coord_a][p][j] *
						  F[coord_b][p][j]);

	      hessian[coord_a][coord_b] += 4.0 * U[coord_a][p][j] *
		                                 U[coord_b][p][j] *
						 evals[p];
              for(qsym=0; qsym < nirreps; qsym++) {
	
                qfirst = first[qsym];
	        qlast = last[qsym];

	        for(q=qfirst; q <= qlast; q++) {

                  for(lsym=0; lsym < nirreps; lsym++) {
		 
                    lfirst = ofirst[lsym];
		    llast = olast[lsym];

		    for(l=lfirst; l <= llast; l++) {
		     
                      ql = INDEX(q,l);
			    
		      hessian[coord_a][coord_b] -= 4.0 * U[coord_a][p][j] * 
                                                         U[coord_b][q][l] *
                                                         A[pj][ql];
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
 
  free(eta);

  /*
  fprintf(outfile, "\n\tSCF Molecular Hessian:\n");
  print_mat(hessian, natom*3, natom*3, outfile);
  */

  /* write the hessian to file15 in the PSI2 standard format */
  ffile(&file15, "file15.dat", 0);
  fprintf(file15, "%5d%5d\n", natom, natom*6);
  for(i=0; i < natom*3; i++) {
    for(j=0; j < natom; j++) {
      fprintf(file15, "%20.10f%20.10f%20.10f\n", 
	      hessian[i][j*3], hessian[i][j*3+1], hessian[i][j*3+2]);
    }
  }
  fclose(file15);
}

}} // namespace psi::cphf
