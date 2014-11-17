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

/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libqt/qt.h>
#include <psi4-dec.h>
#include <libmints/wavefunction.h>
#include "globaldefs.h"
#include "structs.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

// void pitzer_arrays(int nirreps, int *frdocc, int *fruocc, int *orbspi, 
//                    int *first, int *last, int *fstact, int *lstact,
//                    int *active);
// double *** construct_evects(int nirreps, int *active, int *orbspi,
//                             int *first, int *last, int *fstact, int *lstact,
//                             int printflag);

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
void mcscf_get_mo_info(Options &options)
{
   int h, i, j, k, tmp, cnt, irrep, errcod, errbad;
   int size;
   double *eig_unsrt;

   /* set these to NULL so we'll know which one(s) to free in cleanup */
   MCSCF_CalcInfo.mo_hess = NULL;
   MCSCF_CalcInfo.mo_hess_diag = NULL;

   // /* information from checkpoint file */
   // chkpt_init(PSIO_OPEN_OLD);
   // CalcInfo.nirreps = chkpt_rd_nirreps();
   // CalcInfo.nmo = chkpt_rd_nmo();
   // MCSCF_CalcInfo.nso = chkpt_rd_nmo(); /* change to nbfso after conversion */
   // CalcInfo.labels = chkpt_rd_irr_labs();
   // CalcInfo.orbs_per_irr = chkpt_rd_orbspi();
   // MCSCF_CalcInfo.enuc = chkpt_rd_enuc();
   // MCSCF_CalcInfo.efzc = chkpt_rd_efzc();
   // CalcInfo.docc = chkpt_rd_clsdpi();
   // CalcInfo.socc = chkpt_rd_openpi();
   // chkpt_close();
 
   MCSCF_CalcInfo.frozen_docc = init_int_array(CalcInfo.nirreps);
   MCSCF_CalcInfo.frozen_uocc = init_int_array(CalcInfo.nirreps);

   //Process::environment.wavefunction()->frzcpi().copy_into_int_array(MCSCF_CalcInfo.frozen_docc);
   for (int h=0; h<CalcInfo.nirreps; h++) {
     MCSCF_CalcInfo.frozen_docc[h] = Process::environment.wavefunction()->frzcpi()[h];
     MCSCF_CalcInfo.frozen_uocc[h] = Process::environment.wavefunction()->frzvpi()[h];
   }
    if(options["FROZEN_DOCC"].has_changed()){
        if(options["FROZEN_DOCC"].size() != CalcInfo.nirreps)
            throw PSIEXCEPTION("FROZEN_DOCC array should be the same size as the number of irreps.");
        for(int h = 0; h < CalcInfo.nirreps; ++h)
            MCSCF_CalcInfo.frozen_docc[h] = options["FROZEN_DOCC"][h].to_integer();
    }
    if(options["FROZEN_UOCC"].has_changed()){
        if(options["FROZEN_UOCC"].size() != CalcInfo.nirreps)
            throw PSIEXCEPTION("FROZEN_UOCC array should be the same size as the number of irreps.");
        for(int h = 0; h < CalcInfo.nirreps; ++h)
            MCSCF_CalcInfo.frozen_uocc[h] = options["FROZEN_UOCC"][h].to_integer();
    }

   MCSCF_CalcInfo.rstr_docc = init_int_array(CalcInfo.nirreps);
   MCSCF_CalcInfo.rstr_uocc = init_int_array(CalcInfo.nirreps);
   CalcInfo.reorder = init_int_array(CalcInfo.nmo);
   CalcInfo.ras_opi = init_int_matrix(MAX_RAS_SPACES,CalcInfo.nirreps);
      
   if (!ras_set2(CalcInfo.nirreps, CalcInfo.nmo, 1, 1,
                CalcInfo.orbs_per_irr, CalcInfo.docc, CalcInfo.socc, 
                MCSCF_CalcInfo.frozen_docc, MCSCF_CalcInfo.frozen_uocc, 
                MCSCF_CalcInfo.rstr_docc, MCSCF_CalcInfo.rstr_uocc,
                CalcInfo.ras_opi, CalcInfo.reorder, 1, 0, options)) 
   { 
     throw PsiException("Error in ras_set().  Aborting.", __FILE__, __LINE__) ;
   }
   

  // /* Compute maximum number of orbitals per irrep including
  // ** and not including fzv
  // */
  // MCSCF_CalcInfo.max_orbs_per_irrep = 0;
  // MCSCF_CalcInfo.max_pop_per_irrep = 0;
  // for (i=0; i<CalcInfo.nirreps; i++) {
  //   if (MCSCF_CalcInfo.max_orbs_per_irrep < CalcInfo.orbs_per_irr[i])
  //     MCSCF_CalcInfo.max_orbs_per_irrep = CalcInfo.orbs_per_irr[i];
  //   if (MCSCF_CalcInfo.max_pop_per_irrep < (CalcInfo.orbs_per_irr[i] - 
  //                                  MCSCF_CalcInfo.frozen_uocc[i]))
  //     MCSCF_CalcInfo.max_pop_per_irrep = CalcInfo.orbs_per_irr[i] -
  //                                  MCSCF_CalcInfo.frozen_uocc[i];      
  // }


  // /* construct the "ordering" array, which maps the other direction */
  // /* i.e. from a CI orbital to a Pitzer orbital                     */
  // CalcInfo.order = init_int_array(CalcInfo.nmo);
  // for (i=0; i<CalcInfo.nmo; i++) {
  //   j = CalcInfo.reorder[i];
  //   CalcInfo.order[j] = i;
  // }


  /* Set up an array to map absolute ci order to relative Pitzer order */
  MCSCF_CalcInfo.ci2relpitz = init_int_array(CalcInfo.nmo);
  for (h=0,cnt=0; h<CalcInfo.nirreps; h++) {
    for (i=0; i<CalcInfo.orbs_per_irr[h]; i++,cnt++) {
      j = CalcInfo.reorder[cnt];
      MCSCF_CalcInfo.ci2relpitz[j] = i;
    }
  } 

  if (MCSCF_Parameters.print_lvl > 4) {
    outfile->Printf("\nPitzer to CI order array = \n");
    for (i=0; i<CalcInfo.nmo; i++) {
      outfile->Printf("%3d ", CalcInfo.reorder[i]);
    }
    outfile->Printf("\n");
  }


  // CalcInfo.nmotri = (CalcInfo.nmo * (CalcInfo.nmo + 1)) / 2 ;

  // if (CalcInfo.nmotri >= IOFF_MAX){
  //   throw PsiException("(get_mo_info): IOFF_MAX may not large enough!",
  //                             __FILE__, __LINE__);
  // }

  // /* transform orbsym vector to new MO order */
  // MCSCF_CalcInfo.orbsym = init_int_array(CalcInfo.nmo);

  // for (i=0,cnt=0; i<CalcInfo.nirreps; i++) {
  //   for (j=0; j<CalcInfo.orbs_per_irr[i]; j++,cnt++) {
  //     k = CalcInfo.reorder[cnt];
  //     MCSCF_CalcInfo.orbsym[k] = i;
  //   }
  // }

  MCSCF_CalcInfo.num_fzv_orbs = 0;  MCSCF_CalcInfo.num_vir_orbs = 0;
  for (i=0; i<CalcInfo.nirreps; i++) {
    MCSCF_CalcInfo.num_fzv_orbs += MCSCF_CalcInfo.frozen_uocc[i];  
    MCSCF_CalcInfo.num_vir_orbs += MCSCF_CalcInfo.rstr_uocc[i];
  }

  MCSCF_CalcInfo.npop = CalcInfo.nmo - MCSCF_CalcInfo.num_fzv_orbs -
    MCSCF_CalcInfo.num_vir_orbs;

  MCSCF_CalcInfo.num_fzc_orbs = 0;
  MCSCF_CalcInfo.num_cor_orbs = 0;
  for (i=0; i<CalcInfo.nirreps; i++) {
    MCSCF_CalcInfo.num_fzc_orbs += MCSCF_CalcInfo.frozen_docc[i];
  } 
  for (i=0; i<CalcInfo.nirreps; i++) {
    MCSCF_CalcInfo.num_cor_orbs += MCSCF_CalcInfo.rstr_docc[i];
  }

  /* construct the MCSCF_CalcInfo.ras_orbs array (may not be of any use now) */
  cnt = 0;
  MCSCF_CalcInfo.fzc_orbs = init_int_matrix(CalcInfo.nirreps,CalcInfo.nmo);
  MCSCF_CalcInfo.cor_orbs = init_int_matrix(CalcInfo.nirreps,CalcInfo.nmo);
  MCSCF_CalcInfo.vir_orbs = init_int_matrix(CalcInfo.nirreps,CalcInfo.nmo);
  MCSCF_CalcInfo.fzv_orbs = init_int_matrix(CalcInfo.nirreps,CalcInfo.nmo);

  /* FZC */
  for (irrep=0; irrep<CalcInfo.nirreps; irrep++)
    for (j=0; j<MCSCF_CalcInfo.frozen_docc[irrep]; j++)
      MCSCF_CalcInfo.fzc_orbs[irrep][j] = cnt++;

  /* COR */
  for (irrep=0; irrep<CalcInfo.nirreps; irrep++)
    for (j=0; j<MCSCF_CalcInfo.rstr_docc[irrep]; j++)
      MCSCF_CalcInfo.cor_orbs[irrep][j] = cnt++;

  /* RAS */
  MCSCF_CalcInfo.ras_orbs = (int ***) malloc (MAX_RAS_SPACES * sizeof(int **));
  for (i=0; i<MAX_RAS_SPACES; i++) {
    MCSCF_CalcInfo.ras_orbs[i] = init_int_matrix(CalcInfo.nirreps,
      CalcInfo.nmo);
    for (irrep=0; irrep<CalcInfo.nirreps; irrep++) {
      for (j=0; j<CalcInfo.ras_opi[i][irrep]; j++) {
        MCSCF_CalcInfo.ras_orbs[i][irrep][j] = cnt++;
      }
    }
  }

  /* VIR */
  for (irrep=0; irrep<CalcInfo.nirreps; irrep++)
    for (j=0; j<MCSCF_CalcInfo.rstr_uocc[irrep]; j++)
      MCSCF_CalcInfo.vir_orbs[irrep][j] = cnt++;

  /* FZV */
  for (irrep=0; irrep<CalcInfo.nirreps; irrep++)
    for (j=0; j<MCSCF_CalcInfo.frozen_uocc[irrep]; j++)
      MCSCF_CalcInfo.fzv_orbs[irrep][j] = cnt++;



  // /* get the Pitzer arrays first, last, fstact, lstact, and active */
  // MCSCF_CalcInfo.first = init_int_array(CalcInfo.nirreps);
  // MCSCF_CalcInfo.last = init_int_array(CalcInfo.nirreps);
  // MCSCF_CalcInfo.fstact = init_int_array(CalcInfo.nirreps);
  // MCSCF_CalcInfo.lstact = init_int_array(CalcInfo.nirreps);
  // MCSCF_CalcInfo.active = init_int_array(CalcInfo.nirreps);

  /* I think I never use this... --CDS 6/12/04
  // pitzer_arrays(CalcInfo.nirreps, MCSCF_CalcInfo.frozen_docc, MCSCF_CalcInfo.frozen_uocc,
  //               CalcInfo.orbs_per_irr, MCSCF_CalcInfo.first, MCSCF_CalcInfo.last,
  //               MCSCF_CalcInfo.fstact, MCSCF_CalcInfo.lstact, MCSCF_CalcInfo.active);
  */

  /* allocate memory to store the MO coefficient matrix symm blocked */

  MCSCF_CalcInfo.mo_coeffs = (double ***) malloc(CalcInfo.nirreps * 
                                           sizeof(double **));
  for (irrep=0; irrep<CalcInfo.nirreps; irrep++) {
    i = CalcInfo.orbs_per_irr[irrep];
    if (i==0) continue;
    MCSCF_CalcInfo.mo_coeffs[irrep] = block_matrix(i,i);   
  }
  
  if (MCSCF_Parameters.print_lvl > 0) {
    outfile->Printf("ORBITALS:");
    outfile->Printf("\n   FROZEN_DOCC   = ");
    for (i=0; i<CalcInfo.nirreps; i++) {
      outfile->Printf("%2d ", MCSCF_CalcInfo.frozen_docc[i]);
    }
    outfile->Printf("\n   RESTR_DOCC    = ");
    for (i=0; i<CalcInfo.nirreps; i++) {
      outfile->Printf("%2d ", MCSCF_CalcInfo.rstr_docc[i]);
    }
    outfile->Printf("\n   DOCC          = ");
    for (i=0; i<CalcInfo.nirreps; i++) {
      outfile->Printf("%2d ", CalcInfo.docc[i]);
    }
    outfile->Printf("\n   SOCC          = ");
    for (i=0; i<CalcInfo.nirreps; i++) {
      outfile->Printf("%2d ", CalcInfo.socc[i]);
    }
    outfile->Printf("\n   RESTR_UOCC    = ");
    for (i=0; i<CalcInfo.nirreps; i++) {
      outfile->Printf("%2d ", MCSCF_CalcInfo.rstr_uocc[i]);
    }
    outfile->Printf("\n   FROZEN_UOCC   = ");
    for (i=0; i<CalcInfo.nirreps; i++) {
      outfile->Printf("%2d ", MCSCF_CalcInfo.frozen_uocc[i]);
    }

    for (i=0; i<MAX_RAS_SPACES; i++) {
      outfile->Printf("\n   RAS %d         = ",i+1);
      for (j=0; j<CalcInfo.nirreps; j++) {
        outfile->Printf("%2d ", CalcInfo.ras_opi[i][j]);
      }
    }
    outfile->Printf("\n");

    outfile->Printf("   MOL ORBS      =   %6d\n", CalcInfo.nmo);
    outfile->Printf("   FROZEN CORE   =   %6d      RESTR CORE   =   %6d\n",
        MCSCF_CalcInfo.num_fzc_orbs, MCSCF_CalcInfo.num_cor_orbs);
    outfile->Printf("\n");
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
// void pitzer_arrays(int nirreps, int *frdocc, int *fruocc, int *orbspi, 
//                    int *first, int *last, int *fstact, int *lstact,int *active)
// {
// 
//   int h;
//   int first_offset, last_offset;
// 
//   /*
//    * Construct first and last index arrays: this defines the first
//    * absolute orbital index (Pitzer ordering) and last absolute orbital
//    * index for each irrep.  When there are no orbitals for an irrep, the
//    * value is -1 for first[] and -2 for last[].  Note that there must be
//    * orbitals in the first irrep (i.e. totally symmetric) for this to work.
//    */
//   for(h=0; h < nirreps; h++) {
//     first[h] = -1;
//     last[h] = -2;
//   }
// 
//   first_offset = 0;
//   last_offset = orbspi[0] - 1; 
//   first[0] = first_offset;
//   last[0] = last_offset;
// 
//   for(h=1; h < nirreps; h++) {
//     first_offset += orbspi[h-1];
//     last_offset += orbspi[h];
//     if(orbspi[h]) {
//       first[h] = first_offset;
//       last[h] = last_offset;
//     }
//   }
// 
//   /*
//    * Construct first and last active index arrays: this defines the first
//    * absolute orbital index (Pitzer ordering) and last absolute orbital
//    * index for each irrep, excluding frozen orbitals.  When there are no
//    * orbitals for an irrep, the value is -1 for first[] and -2 for last[].
//    * Note that there must be orbitals in the first irrep (i.e. totally
//    * symmetric) for this to work.  
//    */
//   for(h=0; h < nirreps; h++) {
//     fstact[h] = -1;
//     lstact[h] = -2;
//   }
// 
//   first_offset = frdocc[0];
//   last_offset = orbspi[0] - fruocc[0] - 1; 
//   fstact[0] = first_offset;
//   lstact[0] = last_offset;
// 
//   for(h=1; h < nirreps; h++) {
//     first_offset += orbspi[h-1]+frdocc[h]-frdocc[h-1];
//     last_offset += orbspi[h] - fruocc[h] + fruocc[h-1];
//     if(orbspi[h]) {
//       fstact[h] = first_offset;
//       lstact[h] = last_offset;
//     }
//   }
// 
//   /* Now define active[] such that frozen orbitals are taken into account */
//   for(h=0; h < nirreps; h++) {
//     active[h] = orbspi[h]-frdocc[h]-fruocc[h];
//   }
// 
// }



/*
** read_cur_orbs
**
** Read in the molecular orbital matrix from PSIF_CHKPT and put them in MCSCF_CalcInfo
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
        MCSCF_CalcInfo.mo_coeffs[h][i][j] = tmat[i][j];
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
          evects[h][row][col] = MCSCF_CalcInfo.mo_matrix[p][q];
        }
      }

      if(printflag) {
        outfile->Prinft("\n\tMolecular Orbitals for Irrep %s\n",
                CalcInfo.labels[h]);
        print_mat(evects[h],orbspi[h],active[h],"outfile");
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

}} // end namespace psi::detci

