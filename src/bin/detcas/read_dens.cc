/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
/*
** READ_DENS.C
**
** Read the one- and two-particle density matrices
**
** C. David Sherrill
** University of California, Berkeley
**
** Based on code from the CLAG program
** April 1998
*/

#include <cstdlib>
#include <cstdio>
#include <libiwl/iwl.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include "globaldefs.h"
#include "globals.h"

namespace psi { namespace detcas {

double **rdopdm(int nbf, int print_lvl, int opdm_file, int erase);
double *rdtpdm(int nbf, int print_lvl, int tpdm_file, int erase);


void read_density_matrices(void)
{

  /* read the one-particle density matrix */
  CalcInfo.opdm = rdopdm(CalcInfo.npop, Params.print_lvl, Params.opdm_file,
                         Params.opdm_erase);

  /* read the two-particle density matrix */
  CalcInfo.tpdm = rdtpdm(CalcInfo.npop, Params.print_lvl, Params.tpdm_file,
                         Params.tpdm_erase);

}



/*
** RDOPDM
**
** reads the one particle density matrix from opdm_file and returns opdm 
** as a block matrix
**
** Taken from CLAG, April 1998
** upgraded to libpsio 6/03 by CDS
*/
double **rdopdm(int nbf, int print_lvl, int opdm_file, int erase)
{

  int i, root, errcod;
  double **opdm;
  char opdm_key[80];

  psio_open(opdm_file, PSIO_OPEN_OLD);

  opdm = block_matrix(nbf, nbf);

  /* if the user hasn't specified a root, just get "the" onepdm */
  if (!ip_exist("ROOT",0)) {
    psio_read_entry(opdm_file, "MO-basis OPDM", (char *) opdm[0], 
                    nbf*nbf*sizeof(double));
  }
  else {
    root = 1;
    errcod = ip_data("ROOT","%d",&root,0);
    sprintf(opdm_key, "MO-basis OPDM Root %d", root);
    psio_read_entry(opdm_file, opdm_key, (char *) opdm[0], 
                    nbf*nbf*sizeof(double));
  }

  if (print_lvl > 3) {
    fprintf(outfile,"\nOne-Particle Density Matrix\n");
    print_mat(opdm, nbf, nbf, outfile);
    fprintf(outfile,"\n\n");
  }

  psio_close(opdm_file, erase ? 0 : 1);
  return(opdm);

}



/*
** RDTPDM
**
** reads the two particle density matrix from tpdm_file and returns tpdm 
** as an array
**
** Note that nbf is really going to be the number of _populated_ orbitals
** (i.e., subtract frozen virtuals).
**
** Taken from CLAG, adapted to symmetrize the tpdm
** C. David Sherrill
** April 1998
*/
double *rdtpdm(int nbf, int print_lvl, int tpdm_file, int erase)
{

  double *tpdm, *symm_tpdm;
  int numslots, sqnbf, ntri;
  int *ioff_lt, i;                    /* offsets for left (or right) indices */
  int p,q,r,s,smax,pq,qp,rs,sr,pqrs,qprs,pqsr,qpsr,target;
  struct iwlbuf TBuff;

  iwl_buf_init(&TBuff, tpdm_file, 0.0, 1, 1);

  sqnbf = nbf*nbf;
  numslots = (sqnbf*(sqnbf+1))/2;
  tpdm = init_array(numslots);

  /* Construct the ioff_lt array (same here as ioff_rt) : different than
   * regular ioff because there is no perm symmetry between left indices
   * or right indices.
   */
  ioff_lt = init_int_array(nbf);
  for (i=0; i<nbf; i++) {
    ioff_lt[i] = i * nbf;
  }

 iwl_buf_rd_all(&TBuff, tpdm, ioff_lt, ioff_lt, 1, ioff,
                (print_lvl>5), outfile);

  if (print_lvl > 4) {
    fprintf(outfile,"Non-symmetrized Two-Particle Density Matrix\n");
    for (p=0; p<nbf; p++) {
      for (q=0; q<nbf; q++) {
        for (r=0; r<=p; r++) {
          smax = (r==p) ? q+1 : nbf;
          for (s=0; s<smax; s++) {
            pq = p * nbf + q;
            rs = r * nbf + s;
            pqrs = INDEX(pq,rs);
            fprintf(outfile, "%2d %2d %2d %2d = %12.6lf\n", p, q, r, s, 
                    tpdm[pqrs]);
          }
        }
      }
    }
    fprintf(outfile,"\n\n");
  }

  iwl_buf_close(&TBuff, 1);
  free(ioff_lt);

  /* now that we have it in memory, let's symmetrize it       */
  /* eventually, move the symmetrization back down to DETCI   */
  /* using user-parameter to specify symmetrize or not        */

  ntri = (nbf * (nbf + 1))/2;
  symm_tpdm = init_array((ntri * (ntri + 1))/2);
  for (p=0,target=0; p<nbf; p++) {
    for (q=0; q<=p; q++) {
      for (r=0; r<=p; r++) {
        smax = (r==p) ? q+1 : r+1;
        for (s=0; s<smax; s++,target++) {

          pq = p * nbf + q;
          qp = q * nbf + p;
          rs = r * nbf + s;
          sr = s * nbf + r;
          pqrs = INDEX(pq,rs);
          qprs = INDEX(qp,rs);
          pqsr = INDEX(pq,sr);
          qpsr = INDEX(qp,sr);

          /* would be 0.25 but the formulae I used for the diag hessian 
           * seem to define the TPDM with the 1/2 back outside */
          symm_tpdm[target] = 0.5 * (tpdm[pqrs] + tpdm[qprs] +
                              tpdm[pqsr] + tpdm[qpsr]);

        }
      }
    }
  }

  free(tpdm);

  if (print_lvl > 4) {
    fprintf(outfile,"Symmetrized Two-Particle Density Matrix\n");
    for (p=0,target=0; p<nbf; p++) {
      for (q=0; q<=p; q++) {
        for (r=0; r<=p; r++) {
          smax = (r==p) ? q+1 : r+1;
          for (s=0; s<smax; s++,target++) {
            fprintf(outfile, "%2d %2d %2d %2d = %12.6lf\n", p, q, r, s, 
                    symm_tpdm[target]);
          }
        }
      }
    }
  }

  return (symm_tpdm);
}

}} // end namespace psi::detcas

