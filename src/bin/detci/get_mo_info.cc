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
    \ingroup DETCI
    \brief This file gets information about MOs and subspaces of them
*/
#include <cstdio>
#include <cstdlib>
#include <boost/lexical_cast.hpp>
#include <libciomr/libciomr.h>
#include <libmints/wavefunction.h>
#include <libmints/mints.h>
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include "structs.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

/*
** get_mo_info()
**
** Reads the checkpoint file and the input file and gets all sorts of 
** useful information about the molecular orbitals (such as their 
** reordering array, the docc array, frozen orbitals, etc.)
**
** Created by C. David Sherrill on 17 November 1994
**
** Updated
** CDS  1/18/95 to read SCF eigenvalues also (for MP2 guess vector)
** CDS  1/ 5/97 to use nifty new ras_set() function (which transqt has been
**              using for some time).
** CDS  4/27/15 to more clearly distinguish between restricted vs frozen
**              spaces and to update to ras_set3() routine that sets all
**              the orbital subspaces
**
** options = Options object used to parse user input
*/
void get_mo_info(Options &options)
{
  int i, j, k, tmp, cnt, irrep, errcod, errbad;
  int size;
  double *eig_unsrt;
//  boost::shared_ptr<Vector> eig_unsrt;
  int parsed_ras1=0, parsed_ras2=0, do_ras4;

  CalcInfo.maxKlist = 0.0;
  CalcInfo.sigma_initialized = 0;

  boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();

  CalcInfo.nirreps = wfn->nirrep();
  CalcInfo.nso = wfn->nso();
  CalcInfo.nmo = wfn->nmo();
  CalcInfo.iopen = 0;
    for(int h=0; h < CalcInfo.nirreps; h++)
      CalcInfo.iopen += wfn->nsopi()[h];
  CalcInfo.labels = wfn->molecule()->irrep_labels();
  CalcInfo.orbs_per_irr = init_int_array(CalcInfo.nirreps);
  CalcInfo.so_per_irr = init_int_array(CalcInfo.nirreps);
  CalcInfo.docc = init_int_array(CalcInfo.nirreps);
  CalcInfo.socc = init_int_array(CalcInfo.nirreps);
  for(int h=0; h < CalcInfo.nirreps; h++) {
    CalcInfo.orbs_per_irr[h] = wfn->nmopi()[h];
    CalcInfo.so_per_irr[h] = wfn->nsopi()[h];
    CalcInfo.docc[h] = wfn->doccpi()[h];
    CalcInfo.socc[h] = wfn->soccpi()[h];
  }
  CalcInfo.enuc = wfn->molecule()->nuclear_repulsion_energy();
  if(wfn->reference_wavefunction())
      CalcInfo.escf = wfn->reference_wavefunction()->reference_energy();
  else
      CalcInfo.escf = wfn->reference_energy();
  CalcInfo.edrc = wfn->efzc();
  eig_unsrt = wfn->epsilon_a()->to_block_vector();

  if (CalcInfo.iopen && Parameters.opentype == PARM_OPENTYPE_NONE) {
    outfile->Printf( "Warning: iopen=1,opentype=none. Making iopen=0\n");
    CalcInfo.iopen = 0;
  }
  else if (!CalcInfo.iopen && (Parameters.opentype == PARM_OPENTYPE_HIGHSPIN
    || Parameters.opentype == PARM_OPENTYPE_SINGLET)) {
    outfile->Printf("Warning: iopen=0,opentype!=closed. Making iopen=1\n");
    CalcInfo.iopen = 1;
  }
  if (Parameters.ref_sym >= CalcInfo.nirreps) {
    outfile->Printf("Warning: ref_sym >= nirreps.  Setting ref_sym=0\n");
    Parameters.ref_sym = 0;
  }


  // these are all initialized to zero at this point 
  CalcInfo.frozen_docc   = init_int_array(CalcInfo.nirreps);
  CalcInfo.rstr_docc     = init_int_array(CalcInfo.nirreps);
  CalcInfo.dropped_docc  = init_int_array(CalcInfo.nirreps);

  CalcInfo.frozen_uocc   = init_int_array(CalcInfo.nirreps);
  CalcInfo.rstr_uocc     = init_int_array(CalcInfo.nirreps);
  CalcInfo.dropped_uocc  = init_int_array(CalcInfo.nirreps);

  CalcInfo.reorder = init_int_array(CalcInfo.nmo);
  CalcInfo.ras_opi = init_int_matrix(4,CalcInfo.nirreps);

  // Get any default frozen core array and pass it to ras_set3(),
  // which will use it as a default for FROZEN_DOCC (if not an MCSCF) or
  // for RESTRICTED_DOCC (if an MCSCF);
  // useful if the user wants to just say "freeze_core = true"
  int *core_guess = init_int_array(CalcInfo.nirreps);
  for (int h=0; h<CalcInfo.nirreps; h++) {
    core_guess[h] = Process::environment.wavefunction()->frzcpi()[h];
  }

  // This routine sets all orbital subspace arrays properly given 
  // some minimal starting information and an Options object
  if (!ras_set3(CalcInfo.nirreps, CalcInfo.nmo, 
                CalcInfo.orbs_per_irr, CalcInfo.docc, CalcInfo.socc,
                CalcInfo.frozen_docc, CalcInfo.frozen_uocc,
                CalcInfo.rstr_docc, CalcInfo.rstr_uocc,
                CalcInfo.ras_opi, core_guess, CalcInfo.reorder, 1, 
                (Parameters.mcscf ? true : false), options))
  {
    throw PsiException("Error in ras_set3(). Aborting.",__FILE__,__LINE__);
  }
  // We're done with the core_guess array now
  free(core_guess);

  for (int h=0; h<CalcInfo.nirreps; h++) {
    CalcInfo.dropped_docc[h] = CalcInfo.frozen_docc[h]+CalcInfo.rstr_docc[h];
    CalcInfo.dropped_uocc[h] = CalcInfo.frozen_uocc[h]+CalcInfo.rstr_uocc[h];
  }

  // Do we have any "restricted" orbitals that will need to be filtered
  // out when we read the two-electron integrals?
  Parameters.filter_ints = 0;
  for (i=0; i<CalcInfo.nirreps; i++) {
    if ((CalcInfo.rstr_docc[i] > 0) || (CalcInfo.rstr_uocc[i] > 0)) {
      Parameters.filter_ints = 1;
    }
  }
  if (Parameters.dertype != "NONE" || Parameters.mcscf == 1) {
    Parameters.filter_ints = 1;
  }


  /* This seems to not be used CDS 4/15
  CalcInfo.max_orbs_per_irrep = 0;
  CalcInfo.max_pop_per_irrep = 0;
  for (i=0; i<CalcInfo.nirreps; i++) {
    if (CalcInfo.max_orbs_per_irrep < CalcInfo.orbs_per_irr[i])
       CalcInfo.max_orbs_per_irrep = CalcInfo.orbs_per_irr[i];
    if (CalcInfo.max_pop_per_irrep < (CalcInfo.orbs_per_irr[i] -
                                   CalcInfo.frozen_uocc[i] - CalcI))
      CalcInfo.max_pop_per_irrep = CalcInfo.orbs_per_irr[i] -
                                   CalcInfo.frozen_uocc[i];
  }
  */

  // construct the "ordering" array, which maps the other direction 
  // i.e., from a CI orbital to a Pitzer orbital                     
  CalcInfo.order = init_int_array(CalcInfo.nmo);
  for (i=0; i<CalcInfo.nmo; i++) {
    j = CalcInfo.reorder[i];
    CalcInfo.order[j] = i;
  }

  if (Parameters.print_lvl > 4) {
    outfile->Printf( "\nReordering array = \n");
    for (i=0; i<CalcInfo.nmo; i++) {
      outfile->Printf( "%3d ", CalcInfo.reorder[i]);
    }
    outfile->Printf( "\n");
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

  // calculate number of electrons 
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
      throw PsiException("For opentype=singlet must have even number of socc electrons!.",__FILE__,__LINE__);
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
    std::string str = "(get_mo_info): Can't handle opentype = ";
    str += boost::lexical_cast<std::string>( Parameters.opentype) ;
    throw PsiException(str,__FILE__,__LINE__);
  }

  // Count up number of frozen and restricted core and virtuals
  // These are all redone with new naming scheme and distinction
  // between frozen and restricted as of CDS 4/15
  //
  // For convenience define CalcInfo.num_drc_orbs as the overall
  // total number of dropped (implicit) core, sum of num_fzc_orbs and
  // num_rsc_orbs.  Replaces most previous instances of CalcInfo.num_fzc_orbs.
  //
  // Likewise num_drv_orbs is the total number of dropped virtuals and
  // is the sum of num_fzv_orbs and num_rsv_orbs, and replaces most 
  // instances of the previous num_fzv_orbs.

  CalcInfo.num_expl_cor_orbs = 0; // this isn't enabled anymore, zero it

  CalcInfo.num_fzc_orbs = 0; // truly frozen core, sum(frozen_docc)
  CalcInfo.num_rsc_orbs = 0; // restricted core, sum(rstr_docc)
  CalcInfo.num_fzv_orbs = 0; // truly frozen virts, sum(frozen_uocc)
  CalcInfo.num_rsv_orbs = 0; // restricted virtuals, sum(rstr_uocc)
  CalcInfo.num_drc_orbs = 0; // dropped core  = num_fzc_orbs + num_rsc_orbs
  CalcInfo.num_drv_orbs = 0; // dropped virts = num_fzv_orbs + num_rsv_orbs

  for (i=0; i<CalcInfo.nirreps; i++) {
    CalcInfo.num_fzc_orbs += CalcInfo.frozen_docc[i];
    CalcInfo.num_rsc_orbs += CalcInfo.rstr_docc[i];
    CalcInfo.num_fzv_orbs += CalcInfo.frozen_uocc[i];
    CalcInfo.num_rsv_orbs += CalcInfo.rstr_uocc[i];
  }
  CalcInfo.num_drc_orbs = CalcInfo.num_fzc_orbs + CalcInfo.num_rsc_orbs;
  CalcInfo.num_drv_orbs = CalcInfo.num_fzv_orbs + CalcInfo.num_rsv_orbs;

  // calculate number of orbitals active in CI 
  // maybe this changes later for cor orbs, depends on where we go w/ it 
  CalcInfo.num_ci_orbs = CalcInfo.nmo - CalcInfo.num_drc_orbs -
    CalcInfo.num_drv_orbs;

  if ((CalcInfo.num_ci_orbs * (CalcInfo.num_ci_orbs + 1)) / 2 > IOFF_MAX) {
    throw PsiException("error: IOFF_MAX not large enough!",__FILE__,__LINE__);
  }

  CalcInfo.num_alp_expl = CalcInfo.num_alp - CalcInfo.num_drc_orbs;
  CalcInfo.num_bet_expl = CalcInfo.num_bet - CalcInfo.num_drc_orbs;
  CalcInfo.npop = CalcInfo.nmo - CalcInfo.num_drv_orbs;

  // construct the CalcInfo.ras_orbs array 
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

  // temporarily go back to old variables until entire switchover is done
  // X CDS 4/15
  // CalcInfo.num_fzc_orbs = CalcInfo.num_drc_orbs;
  // CalcInfo.num_fzv_orbs = CalcInfo.num_drv_orbs;
  // for (int h=0; h<CalcInfo.nirreps; h++) {
  //   CalcInfo.frozen_docc[h] = CalcInfo.dropped_docc[h];
  //   CalcInfo.frozen_uocc[h] = CalcInfo.dropped_uocc[h];
  // }

} // end get_mo_info()

}} // namespace psi::detci

