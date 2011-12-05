/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libmints/wavefunction.h>
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include "structs.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

/*
** GET_MO_INFO
**
** Reads PSIF_CHKPT & input.dat and gets all sorts of useful information about
** the molecular orbitals (such as their reordering array, the docc
** array, frozen orbitals, etc.)
**
** Created by C. David Sherrill on 17 November 1994
**
** Updated
** CDS  1/18/95 to read SCF eigenvalues also (for MP2 guess vector)
** CDS  1/ 5/97 to use nifty new ras_set() function (which transqt has been
**              using for some time).
**
*/
void get_mo_info(Options &options)
{
   int i, j, k, tmp, cnt, irrep, errcod, errbad;
   int size;
   double *eig_unsrt;
   int parsed_ras1=0, parsed_ras2=0, do_ras4;
   int *rstr_docc, *rstr_uocc;

   CalcInfo.maxKlist = 0.0;
   CalcInfo.sigma_initialized = 0;

   chkpt_init(PSIO_OPEN_OLD);
   CalcInfo.nirreps = chkpt_rd_nirreps();
   CalcInfo.nso = chkpt_rd_nmo();
   CalcInfo.nmo = chkpt_rd_nmo();
   CalcInfo.iopen = chkpt_rd_iopen();
   CalcInfo.labels = chkpt_rd_irr_labs();
   CalcInfo.orbs_per_irr = chkpt_rd_orbspi();
   CalcInfo.so_per_irr = chkpt_rd_sopi();
   CalcInfo.docc = chkpt_rd_clsdpi();
   CalcInfo.socc = chkpt_rd_openpi();
   CalcInfo.enuc = chkpt_rd_enuc();
   CalcInfo.escf = chkpt_rd_escf();
   CalcInfo.efzc = chkpt_rd_efzc();
   eig_unsrt = chkpt_rd_evals();
   chkpt_close();

   if (CalcInfo.iopen && Parameters.opentype == PARM_OPENTYPE_NONE) {
      fprintf(outfile, "Warning: iopen=1,opentype=none. Making iopen=0\n");
      CalcInfo.iopen = 0;
      }
   else if (!CalcInfo.iopen && (Parameters.opentype == PARM_OPENTYPE_HIGHSPIN
      || Parameters.opentype == PARM_OPENTYPE_SINGLET)) {
      fprintf(outfile,"Warning: iopen=0,opentype!=closed. Making iopen=1\n");
      CalcInfo.iopen = 1;
      }
   if (Parameters.ref_sym >= CalcInfo.nirreps) {
      fprintf(outfile,"Warning: ref_sym >= nirreps.  Setting ref_sym=0\n");
      Parameters.ref_sym = 0;
      }

   //CalcInfo.frozen_docc = init_int_array(CalcInfo.nirreps);
   //CalcInfo.frozen_uocc = init_int_array(CalcInfo.nirreps);
   CalcInfo.frozen_docc =
     Process::environment.reference_wavefunction()->frzcpi();
   CalcInfo.frozen_uocc =
     Process::environment.reference_wavefunction()->frzvpi();

   rstr_docc = init_int_array(CalcInfo.nirreps);
   rstr_uocc = init_int_array(CalcInfo.nirreps);
   CalcInfo.explicit_core = init_int_array(CalcInfo.nirreps);
   CalcInfo.explicit_vir  = init_int_array(CalcInfo.nirreps);
   CalcInfo.reorder = init_int_array(CalcInfo.nmo);
   CalcInfo.ras_opi = init_int_matrix(4,CalcInfo.nirreps);

   if (!ras_set2(CalcInfo.nirreps, CalcInfo.nmo, 1, (Parameters.fzc) ?  1:0,
                CalcInfo.orbs_per_irr, CalcInfo.docc, CalcInfo.socc,
                CalcInfo.frozen_docc, CalcInfo.frozen_uocc,
                rstr_docc, rstr_uocc,
                CalcInfo.ras_opi, CalcInfo.reorder, 1, 0, options))
   {
     fprintf(outfile, "Error in ras_set().  Aborting.\n");
     exit(1);
   }

  /* Check if there are any restricted orbitals.  If there are, then
     for now, I need to filter the integrals, because restricted
     orbitals are presently treated as frozen -CDS, 11/22/2011
  */

   Parameters.filter_ints = 0;

   if (1) { /* for now, always treat restricted as frozen */
     for (i=0; i<CalcInfo.nirreps; i++) {
       CalcInfo.frozen_docc[i] += rstr_docc[i];
       CalcInfo.frozen_uocc[i] += rstr_uocc[i];
       if (rstr_docc[i] > 0 || rstr_uocc[i] > 0) Parameters.filter_ints = 1;
     }
   }
   else { /* for future use */
     for (i=0; i<CalcInfo.nirreps; i++) {
       CalcInfo.explicit_core[i] = rstr_docc[i];
       CalcInfo.explicit_vir[i]  = rstr_uocc[i];
     }
   }

   free(rstr_docc);  free(rstr_uocc);

   if (Parameters.dertype != "NONE" || Parameters.wfn == "DETCAS" ||
       Parameters.wfn == "CASSCF"   || Parameters.wfn == "RASSCF") {
     Parameters.filter_ints = 1;
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
   CalcInfo.order = init_int_array(CalcInfo.nmo);
   for (i=0; i<CalcInfo.nmo; i++) {
      j = CalcInfo.reorder[i];
      CalcInfo.order[j] = i;
      }


   if (Parameters.print_lvl > 4) {
      fprintf(outfile, "\nReordering array = \n");
      for (i=0; i<CalcInfo.nmo; i++) {
         fprintf(outfile, "%3d ", CalcInfo.reorder[i]);
         }
      fprintf(outfile, "\n");
      }

   CalcInfo.nmotri = (CalcInfo.nmo * (CalcInfo.nmo + 1)) / 2 ;

   /* transform orbsym vector to new MO order */
   CalcInfo.orbsym = init_int_array(CalcInfo.nmo);
   CalcInfo.scfeigval = init_array(CalcInfo.nmo);
   if(Parameters.zaptn) {
     CalcInfo.scfeigvala = init_array(CalcInfo.nmo);
     CalcInfo.scfeigvalb = init_array(CalcInfo.nmo);
     }

   for (i=0,cnt=0; i<CalcInfo.nirreps; i++) {
      for (j=0; j<CalcInfo.orbs_per_irr[i]; j++,cnt++) {
         k = CalcInfo.reorder[cnt];
         CalcInfo.orbsym[k] = i;
         }
      }

   for (i=0; i<CalcInfo.nmo; i++) {
      j = CalcInfo.reorder[i];
      CalcInfo.scfeigval[j] = eig_unsrt[i];
      if(Parameters.zaptn) {
        CalcInfo.scfeigvala[j] = eig_unsrt[i];
        CalcInfo.scfeigvalb[j] = eig_unsrt[i];
        }
      }
   free(eig_unsrt);

   /* calculate number of electrons */
   CalcInfo.num_alp = CalcInfo.num_bet = CalcInfo.spab = 0;
   if (Parameters.opentype == PARM_OPENTYPE_NONE ||
       Parameters.opentype == PARM_OPENTYPE_HIGHSPIN) {
      for (i=0; i<CalcInfo.nirreps; i++) {
         CalcInfo.num_alp += CalcInfo.docc[i] + CalcInfo.socc[i];
         CalcInfo.num_bet += CalcInfo.docc[i];
         }
      }
   else if (Parameters.opentype == PARM_OPENTYPE_SINGLET) {
      for (i=0; i<CalcInfo.nirreps; i++) { /* closed-shell part */
         CalcInfo.spab += CalcInfo.socc[i];
         CalcInfo.num_alp += CalcInfo.docc[i];
         CalcInfo.num_bet += CalcInfo.docc[i];
         }
      if (CalcInfo.spab % 2) {
         fprintf(outfile,"For opentype=singlet must have even number ");
         fprintf(outfile,"of socc electrons!\n");
         exit(1);
         }
      CalcInfo.spab /= 2;
      tmp = 0;
      for (i=0; i<CalcInfo.nirreps; i++) {
         j = CalcInfo.socc[i];
         k = 0;
         while (k < j) {
            if (tmp < CalcInfo.spab) {
               CalcInfo.num_alp++;
               tmp++;
               k++;
               }
            else {
               CalcInfo.num_bet++;
               tmp++;
               k++;
               }
            }
         }
      }
   else {
      fprintf(outfile, "(get_mo_info): Can't handle opentype = %d\n",
         Parameters.opentype);
      exit(1);
      }

   /* at this stage I've already overwritten CalcInfo.frozen_docc,
      CalcInfo.rstr_docc, etc, to be their internal DETCI meaning
      and not necessarily the user input.  Internally frozen_docc
      and frozen_uocc refer to dropped core/virt that are never
      allowed to change occupation
   */
   CalcInfo.num_fzv_orbs = 0;
   CalcInfo.num_vir_orbs = 0;
   CalcInfo.num_fzc_orbs = 0;
   CalcInfo.num_cor_orbs = 0;
   for (i=0; i<CalcInfo.nirreps; i++) {
      CalcInfo.num_fzv_orbs += CalcInfo.frozen_uocc[i];
      CalcInfo.num_vir_orbs += CalcInfo.explicit_vir[i];
      CalcInfo.num_fzc_orbs += CalcInfo.frozen_docc[i];
      CalcInfo.num_cor_orbs += CalcInfo.explicit_core[i];
   }

   /* calculate number of orbitals active in CI */
   /* maybe this changes later for cor orbs, depends on where we go w/ it */
   CalcInfo.num_ci_orbs = CalcInfo.nmo - CalcInfo.num_fzc_orbs -
     CalcInfo.num_fzv_orbs;

   if ((CalcInfo.num_ci_orbs * (CalcInfo.num_ci_orbs + 1)) / 2 > IOFF_MAX) {
      fprintf(outfile, "Error: IOFF_MAX not large enough!\n");
      exit(1);
   }

   CalcInfo.num_alp_expl = CalcInfo.num_alp - CalcInfo.num_fzc_orbs;
   CalcInfo.num_bet_expl = CalcInfo.num_bet - CalcInfo.num_fzc_orbs;

   /* construct the CalcInfo.ras_orbs array (may not be of any use now) */
   cnt = 0;
   for (i=0; i<4; i++) {
     CalcInfo.ras_orbs[i] = init_int_matrix(CalcInfo.nirreps,
       CalcInfo.num_ci_orbs);
     for (irrep=0; irrep<CalcInfo.nirreps; irrep++) {
       for (j=0; j<CalcInfo.ras_opi[i][irrep]; j++) {
         CalcInfo.ras_orbs[i][irrep][j] = cnt++;
       }
     }
   }
}

}} // namespace psi::detci

