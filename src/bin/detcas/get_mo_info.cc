/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
#include <cstdlib>
#include <cstdio>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include "globaldefs.h"
#include "globals.h"

namespace psi { namespace detcas {

void pitzer_arrays(int nirreps, int *frdocc, int *fruocc, int *orbspi, 
                   int *first, int *last, int *fstact, int *lstact,
                   int *active);
double *** construct_evects(int nirreps, int *active, int *orbspi,
                            int *first, int *last, int *fstact, int *lstact,
                            int printflag);
extern void check(int a, const char *errmsg);

/*
** GET_MO_INFO
** 
** Reads PSIF_CHKPT & input.dat and gets all sorts of useful information about 
** the molecular orbitals (such as their reordering array, the docc
** array, frozen orbitals, etc.)
**
** Created by C. David Sherrill on 24 April 1998,
** based on the version in DETCI
**
*/
void get_mo_info(void)
{
   int h, i, j, k, tmp, cnt, irrep, errcod, errbad;
   int size;
   double *eig_unsrt;

   /* set these to NULL so we'll know which one(s) to free in cleanup */
   CalcInfo.mo_hess = NULL;
   CalcInfo.mo_hess_diag = NULL;

   /* information from checkpoint file */
   chkpt_init(PSIO_OPEN_OLD);
   CalcInfo.nirreps = chkpt_rd_nirreps();
   CalcInfo.nmo = chkpt_rd_nmo();
   CalcInfo.nbfso = chkpt_rd_nmo(); /* change to nbfso after conversion */
   CalcInfo.labels = chkpt_rd_irr_labs();
   CalcInfo.orbs_per_irr = chkpt_rd_orbspi();
   CalcInfo.enuc = chkpt_rd_enuc();
   CalcInfo.efzc = chkpt_rd_efzc();
   CalcInfo.docc = chkpt_rd_clsdpi();
   CalcInfo.socc = chkpt_rd_openpi();
   chkpt_close();
 
   CalcInfo.frozen_docc = init_int_array(CalcInfo.nirreps);
   CalcInfo.frozen_uocc = init_int_array(CalcInfo.nirreps);
   CalcInfo.rstr_docc = init_int_array(CalcInfo.nirreps);
   CalcInfo.rstr_uocc = init_int_array(CalcInfo.nirreps);
   CalcInfo.pitz2ci = init_int_array(CalcInfo.nmo);
   CalcInfo.ras_opi = init_int_matrix(MAX_RAS_SPACES,CalcInfo.nirreps);
      
   if (!ras_set2(CalcInfo.nirreps, CalcInfo.nmo, 1, 1,
                CalcInfo.orbs_per_irr, CalcInfo.docc, CalcInfo.socc, 
                CalcInfo.frozen_docc, CalcInfo.frozen_uocc, 
                CalcInfo.rstr_docc, CalcInfo.rstr_uocc,
                CalcInfo.ras_opi, CalcInfo.pitz2ci, 1, 0)) 
   { 
     fprintf(outfile, "Error in ras_set().  Aborting.\n");
     exit(1);
   }
   

  /* Compute maximum number of orbitals per irrep including
  ** and not including fzv
  */
  CalcInfo.max_orbs_per_irrep = 0;
  CalcInfo.max_pop_per_irrep = 0;
  for (i=0; i<CalcInfo.nirreps; i++) {
    if (CalcInfo.max_orbs_per_irrep < CalcInfo.orbs_per_irr[i])
      CalcInfo.max_orbs_per_irrep = CalcInfo.orbs_per_irr[i];
    if (CalcInfo.max_pop_per_irrep < (CalcInfo.orbs_per_irr[i] - 
                                   CalcInfo.frozen_uocc[i]))
      CalcInfo.max_pop_per_irrep = CalcInfo.orbs_per_irr[i] -
                                   CalcInfo.frozen_uocc[i];      
  }


  /* construct the "ordering" array, which maps the other direction */
  /* i.e. from a CI orbital to a Pitzer orbital                     */
  CalcInfo.ci2pitz = init_int_array(CalcInfo.nmo);
  for (i=0; i<CalcInfo.nmo; i++) {
    j = CalcInfo.pitz2ci[i];
    CalcInfo.ci2pitz[j] = i;
  }


  /* Set up an array to map absolute ci order to relative Pitzer order */
  CalcInfo.ci2relpitz = init_int_array(CalcInfo.nmo);
  for (h=0,cnt=0; h<CalcInfo.nirreps; h++) {
    for (i=0; i<CalcInfo.orbs_per_irr[h]; i++,cnt++) {
      j = CalcInfo.pitz2ci[cnt];
      CalcInfo.ci2relpitz[j] = i;
    }
  } 

  if (Params.print_lvl > 4) {
    fprintf(outfile, "\nPitzer to CI order array = \n");
    for (i=0; i<CalcInfo.nmo; i++) {
      fprintf(outfile, "%3d ", CalcInfo.pitz2ci[i]);
    }
    fprintf(outfile, "\n");
  }


  CalcInfo.nbstri = (CalcInfo.nmo * (CalcInfo.nmo + 1)) / 2 ;
  check((CalcInfo.nbstri <= IOFF_MAX), 
        "(get_mo_info): IOFF_MAX may not large enough!");

  /* transform orbsym vector to new MO order */
  CalcInfo.orbsym = init_int_array(CalcInfo.nmo);

  for (i=0,cnt=0; i<CalcInfo.nirreps; i++) {
    for (j=0; j<CalcInfo.orbs_per_irr[i]; j++,cnt++) {
      k = CalcInfo.pitz2ci[cnt];
      CalcInfo.orbsym[k] = i;
    }
  }

  CalcInfo.num_fzv_orbs = 0;  CalcInfo.num_vir_orbs = 0;
  for (i=0; i<CalcInfo.nirreps; i++) {
    CalcInfo.num_fzv_orbs += CalcInfo.frozen_uocc[i];  
    CalcInfo.num_vir_orbs += CalcInfo.rstr_uocc[i];
  }

  CalcInfo.npop = CalcInfo.nmo - CalcInfo.num_fzv_orbs -
    CalcInfo.num_vir_orbs;

  CalcInfo.num_fzc_orbs = 0;
  CalcInfo.num_cor_orbs = 0;
  for (i=0; i<CalcInfo.nirreps; i++) {
    CalcInfo.num_fzc_orbs += CalcInfo.frozen_docc[i];
  } 
  for (i=0; i<CalcInfo.nirreps; i++) {
    CalcInfo.num_cor_orbs += CalcInfo.rstr_docc[i];
  }

  /* construct the CalcInfo.ras_orbs array (may not be of any use now) */
  cnt = 0;
  CalcInfo.fzc_orbs = init_int_matrix(CalcInfo.nirreps,CalcInfo.nmo);
  CalcInfo.cor_orbs = init_int_matrix(CalcInfo.nirreps,CalcInfo.nmo);
  CalcInfo.vir_orbs = init_int_matrix(CalcInfo.nirreps,CalcInfo.nmo);
  CalcInfo.fzv_orbs = init_int_matrix(CalcInfo.nirreps,CalcInfo.nmo);

  /* FZC */
  for (irrep=0; irrep<CalcInfo.nirreps; irrep++)
    for (j=0; j<CalcInfo.frozen_docc[irrep]; j++)
      CalcInfo.fzc_orbs[irrep][j] = cnt++;

  /* COR */
  for (irrep=0; irrep<CalcInfo.nirreps; irrep++)
    for (j=0; j<CalcInfo.rstr_docc[irrep]; j++)
      CalcInfo.cor_orbs[irrep][j] = cnt++;

  /* RAS */
  CalcInfo.ras_orbs = (int ***) malloc (MAX_RAS_SPACES * sizeof(int **));
  for (i=0; i<MAX_RAS_SPACES; i++) {
    CalcInfo.ras_orbs[i] = init_int_matrix(CalcInfo.nirreps,
      CalcInfo.nmo);
    for (irrep=0; irrep<CalcInfo.nirreps; irrep++) {
      for (j=0; j<CalcInfo.ras_opi[i][irrep]; j++) {
        CalcInfo.ras_orbs[i][irrep][j] = cnt++;
      }
    }
  }

  /* VIR */
  for (irrep=0; irrep<CalcInfo.nirreps; irrep++)
    for (j=0; j<CalcInfo.rstr_uocc[irrep]; j++)
      CalcInfo.vir_orbs[irrep][j] = cnt++;

  /* FZV */
  for (irrep=0; irrep<CalcInfo.nirreps; irrep++)
    for (j=0; j<CalcInfo.frozen_uocc[irrep]; j++)
      CalcInfo.fzv_orbs[irrep][j] = cnt++;



  /* get the Pitzer arrays first, last, fstact, lstact, and active */
  CalcInfo.first = init_int_array(CalcInfo.nirreps);
  CalcInfo.last = init_int_array(CalcInfo.nirreps);
  CalcInfo.fstact = init_int_array(CalcInfo.nirreps);
  CalcInfo.lstact = init_int_array(CalcInfo.nirreps);
  CalcInfo.active = init_int_array(CalcInfo.nirreps);

  /* I think I never use this... --CDS 6/12/04
  pitzer_arrays(CalcInfo.nirreps, CalcInfo.frozen_docc, CalcInfo.frozen_uocc,
                CalcInfo.orbs_per_irr, CalcInfo.first, CalcInfo.last,
                CalcInfo.fstact, CalcInfo.lstact, CalcInfo.active);
  */

  /* allocate memory to store the MO coefficient matrix symm blocked */

  CalcInfo.mo_coeffs = (double ***) malloc(CalcInfo.nirreps * 
                                           sizeof(double **));
  for (irrep=0; irrep<CalcInfo.nirreps; irrep++) {
    i = CalcInfo.orbs_per_irr[irrep];
    if (i==0) continue;
    CalcInfo.mo_coeffs[irrep] = block_matrix(i,i);   
  }
  
  if (Params.print_lvl > 0) {
    fprintf(outfile, "ORBITALS:");
    fprintf(outfile, "\n   FROZEN_DOCC   = ");
    for (i=0; i<CalcInfo.nirreps; i++) {
      fprintf(outfile, "%2d ", CalcInfo.frozen_docc[i]);
    }
    fprintf(outfile, "\n   RESTR_DOCC    = ");
    for (i=0; i<CalcInfo.nirreps; i++) {
      fprintf(outfile, "%2d ", CalcInfo.rstr_docc[i]);
    }
    fprintf(outfile, "\n   DOCC          = ");
    for (i=0; i<CalcInfo.nirreps; i++) {
      fprintf(outfile, "%2d ", CalcInfo.docc[i]);
    }
    fprintf(outfile, "\n   SOCC          = ");
    for (i=0; i<CalcInfo.nirreps; i++) {
      fprintf(outfile, "%2d ", CalcInfo.socc[i]);
    }
    fprintf(outfile, "\n   RESTR_UOCC    = ");
    for (i=0; i<CalcInfo.nirreps; i++) {
      fprintf(outfile, "%2d ", CalcInfo.rstr_uocc[i]);
    }
    fprintf(outfile, "\n   FROZEN_UOCC   = ");
    for (i=0; i<CalcInfo.nirreps; i++) {
      fprintf(outfile, "%2d ", CalcInfo.frozen_uocc[i]);
    }

    for (i=0; i<MAX_RAS_SPACES; i++) {
      fprintf(outfile, "\n   RAS %d         = ",i+1);
      for (j=0; j<CalcInfo.nirreps; j++) {
        fprintf(outfile, "%2d ", CalcInfo.ras_opi[i][j]);
      }
    }
    fprintf(outfile, "\n");

    fprintf(outfile, "   MOL ORBS      =   %6d\n", CalcInfo.nmo);
    fprintf(outfile, "   FROZEN CORE   =   %6d      RESTR CORE   =   %6d\n",
        CalcInfo.num_fzc_orbs, CalcInfo.num_cor_orbs);
    fprintf(outfile, "\n");
  }
}



/*
** pitzer_arrays
**
** Form the first/last/active arrays for the orbitals in Pitzer order
** Based on code taken from TRANSQT
**
** C. David Sherrill
** April 1998
*/
void pitzer_arrays(int nirreps, int *frdocc, int *fruocc, int *orbspi, 
                   int *first, int *last, int *fstact, int *lstact,int *active)
{

  int h;
  int first_offset, last_offset;

  /*
   * Construct first and last index arrays: this defines the first
   * absolute orbital index (Pitzer ordering) and last absolute orbital
   * index for each irrep.  When there are no orbitals for an irrep, the
   * value is -1 for first[] and -2 for last[].  Note that there must be
   * orbitals in the first irrep (i.e. totally symmetric) for this to work.
   */
  for(h=0; h < nirreps; h++) {
    first[h] = -1;
    last[h] = -2;
  }

  first_offset = 0;
  last_offset = orbspi[0] - 1; 
  first[0] = first_offset;
  last[0] = last_offset;

  for(h=1; h < nirreps; h++) {
    first_offset += orbspi[h-1];
    last_offset += orbspi[h];
    if(orbspi[h]) {
      first[h] = first_offset;
      last[h] = last_offset;
    }
  }

  /*
   * Construct first and last active index arrays: this defines the first
   * absolute orbital index (Pitzer ordering) and last absolute orbital
   * index for each irrep, excluding frozen orbitals.  When there are no
   * orbitals for an irrep, the value is -1 for first[] and -2 for last[].
   * Note that there must be orbitals in the first irrep (i.e. totally
   * symmetric) for this to work.  
   */
  for(h=0; h < nirreps; h++) {
    fstact[h] = -1;
    lstact[h] = -2;
  }

  first_offset = frdocc[0];
  last_offset = orbspi[0] - fruocc[0] - 1; 
  fstact[0] = first_offset;
  lstact[0] = last_offset;

  for(h=1; h < nirreps; h++) {
    first_offset += orbspi[h-1]+frdocc[h]-frdocc[h-1];
    last_offset += orbspi[h] - fruocc[h] + fruocc[h-1];
    if(orbspi[h]) {
      fstact[h] = first_offset;
      lstact[h] = last_offset;
    }
  }

  /* Now define active[] such that frozen orbitals are taken into account */
  for(h=0; h < nirreps; h++) {
    active[h] = orbspi[h]-frdocc[h]-fruocc[h];
  }

}



/*
** read_cur_orbs
**
** Read in the molecular orbital matrix from PSIF_CHKPT and put them in CalcInfo
*/
void read_cur_orbs(void)
{
  int i, j, h, dim, nirreps;
  double **tmat;

  nirreps = CalcInfo.nirreps;

  chkpt_init(PSIO_OPEN_OLD);
  for (h=0; h<nirreps; h++) {
    dim = CalcInfo.orbs_per_irr[h];
    if (dim==0) continue;
    tmat = chkpt_rd_scf_irrep(h); 
    for (i=0; i<dim; i++) 
      for (j=0; j<dim; j++) 
        CalcInfo.mo_coeffs[h][i][j] = tmat[i][j];
    free_block(tmat);
  }
  chkpt_close();

}



/*
** construct_evects
**
** This function copies SCF eigenvectors from moinfo.scf_vector into
** an array of matrices.  Columns corresponding to inactive orbitals
** may be deleted if desired.
**
** Code taken from TRANSQT
** C. David Sherrill
** April 1998
**
double *** construct_evects(int nirreps, int *active, int *orbspi,
                            int *first, int *last, int *fstact, int *lstact,
                            int printflag)
{

  int h, row, col, p, q;
  double ***evects;

  evects = (double ***) malloc(nirreps * sizeof(double **));

  for (h=0; h<nirreps; h++) {
    if (active[h]) {
      evects[h] = block_matrix(orbspi[h],active[h]);
      row = -1;
      for(p=first[h]; p <= last[h]; p++) {
        row++; col = -1;
        for(q=fstact[h]; q <= lstact[h]; q++) {
          col++;
          evects[h][row][col] = CalcInfo.mo_matrix[p][q];
        }
      }

      if(printflag) {
        fprintf(outfile,"\n\tMolecular Orbitals for Irrep %s\n",
                CalcInfo.labels[h]);
        print_mat(evects[h],orbspi[h],active[h],outfile);
      }
    }

    else
      evects[h] = NULL;

  }

  return(evects);

}

void destruct_evects(int nirreps, double ***evects)
{
  int h;

  for (h=0; h<nirreps; h++)
    if (evects[h] != NULL) free_block(evects[h]);
}
*/

}} // end namespace psi::detcas

