/*! \file
    \ingroup CPHF
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi { namespace cphf {

/* mohess(): Builds the molecular orbital hessian matrix, which is
** needed for solving the CPHF equations and for Hartree-Fock wave
** function stability analysis.
**
** The MO hessian may be expressed as:
**
** A_ai,bj = delta_ab delta_ij (eps_i - eps_a) 
**            - [4 (ai|bj) - (ab|ij) - (aj|ib) ]
**
** For more details (and clearer notation) see:
**
** Y. Yamaguchi et al., "A New Dimension to Quantum Chemistry:
** Analytic Derivative Methods in Ab Initio Molecular Electronic
** Structure Theory", Oxford Press, New York, 1994.  Ch.10,
** esp. pp. 128-132.
**
** Note that, although the loop structure here involves irreps,
** symmetry is only used because of the Pitzer ordering of the MO's
** (all orbitals together in each irrep).  Symmetry is not actually
** used to streamline the calculation or storage.  This could be done
** (and *should* be) by computing offsets for symmetry blocks similar
** to the direct-product decomposition approach used in the PSI3
** coupled cluster codes.  However, when the CPHF equations are
** solved, the perturbations themselves must be symmetrized or the
** irrep structure of the CPHF coefficients will be lost.
**
** TDC, December 2001 (revised October 2002)
*/

void mohess(double **Aaibj)
{
  int asym, isym, bsym, jsym;
  int a, i, b, j;
  int afirst, alast, ifirst, ilast;
  int bfirst, blast, jfirst, jlast;
  int AI, BJ;

  for(asym=0,AI=0; asym < nirreps; asym++) {

    afirst = vfirst[asym];
    alast = vlast[asym];

    for(a=afirst; a <= alast; a++) {

      for(isym=0; isym < nirreps; isym++) {
	   
        ifirst = ofirst[isym];
	ilast = olast[isym];

	for(i=ifirst; i <= ilast; i++,AI++) {
           
	  for(bsym=0,BJ=0; bsym < nirreps; bsym++) {

	    bfirst = vfirst[bsym];
	    blast = vlast[bsym];

	    for(b=bfirst; b <= blast; b++) {
               
	      for(jsym=0; jsym < nirreps; jsym++) {

	        jfirst = ofirst[jsym];
	        jlast = olast[jsym];

	        for(j=jfirst; j <= jlast; j++,BJ++) {
                   
	          Aaibj[AI][BJ] += (a==b) * (i==j) * (evals[i] - evals[a]);
                }
              }
            }
          }
        }
      }
    }
  }
  
  if (print_lvl > 5) {
    fprintf(outfile, "MO Hessian A(ai,bj): \n");
    print_mat(Aaibj, num_ai, num_ai, outfile);
  }

  /* dump the hessian to disk */
  psio_open(PSIF_CPHF, 1);
  psio_write_entry(PSIF_CPHF, "RHF MO Hessian", (char *) Aaibj[0], 
                   num_ai*num_ai*sizeof(double));
  psio_close(PSIF_CPHF, 1);

  return;
}

}} // namespace psi::cphf
