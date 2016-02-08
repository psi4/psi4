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
#include <libmints/molecule.h>
#include <libmints/vector.h>
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include <ciwave.h>
#include "structs.h"

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
** DGAS 7/1/15 Moving this into CIWavefunction so we can get rid of chkpt
** options = Options object used to parse user input
*/
void CIWavefunction::get_mo_info()
{
  int i, j, k, tmp, cnt, irrep, errcod, errbad;
  int size;
  double *eig_unsrt;
  int parsed_ras1=0, parsed_ras2=0, do_ras4;

  CalcInfo_->maxKlist = 0.0;
  CalcInfo_->sigma_initialized = 0;

  // Initial guess will overwrite some of this later.
  CalcInfo_->nirreps = reference_wavefunction_->nirrep();
  CalcInfo_->nso = reference_wavefunction_->nso();
  CalcInfo_->nmo = reference_wavefunction_->nmo();

  // DGAS: David does this make sense?
  CalcInfo_->iopen = !reference_wavefunction_->same_a_b_orbs();
  CalcInfo_->labels = reference_wavefunction_->molecule()->irrep_labels();
  CalcInfo_->orbs_per_irr = reference_wavefunction_->nmopi();
  CalcInfo_->so_per_irr = reference_wavefunction_->nsopi();
  CalcInfo_->docc = reference_wavefunction_->doccpi();
  CalcInfo_->socc = reference_wavefunction_->soccpi();
  CalcInfo_->enuc = reference_wavefunction_->molecule()->nuclear_repulsion_energy();
  CalcInfo_->escf = reference_wavefunction_->reference_energy();
  eig_unsrt = reference_wavefunction_->epsilon_a()->to_block_vector();
  CalcInfo_->edrc = 0.0;

  if (CalcInfo_->iopen && Parameters_->opentype == PARM_OPENTYPE_NONE) {
    outfile->Printf( "Warning: iopen=1,opentype=none. Making iopen=0\n");
    CalcInfo_->iopen = 0;
  }
  else if (!CalcInfo_->iopen && (Parameters_->opentype == PARM_OPENTYPE_HIGHSPIN
    || Parameters_->opentype == PARM_OPENTYPE_SINGLET)) {
    outfile->Printf("Warning: iopen=0,opentype!=closed. Making iopen=1\n");
    CalcInfo_->iopen = 1;
  }
  if (Parameters_->ref_sym >= CalcInfo_->nirreps) {
    outfile->Printf("Warning: ref_sym >= nirreps.  Setting ref_sym=0\n");
    Parameters_->ref_sym = 0;
  }


  // these are all initialized to zero at this point
  CalcInfo_->frozen_docc   = init_int_array(CalcInfo_->nirreps);
  CalcInfo_->rstr_docc     = init_int_array(CalcInfo_->nirreps);
  CalcInfo_->dropped_docc  = init_int_array(CalcInfo_->nirreps);

  CalcInfo_->ci_orbs       = init_int_array(CalcInfo_->nirreps);

  CalcInfo_->frozen_uocc   = init_int_array(CalcInfo_->nirreps);
  CalcInfo_->rstr_uocc     = init_int_array(CalcInfo_->nirreps);
  CalcInfo_->dropped_uocc  = init_int_array(CalcInfo_->nirreps);

  CalcInfo_->reorder = init_int_array(CalcInfo_->nmo);
  CalcInfo_->ras_opi = init_int_matrix(4,CalcInfo_->nirreps);

  // Get any default frozen core array and pass it to ras_set3(),
  // which will use it as a default for FROZEN_DOCC (if not an MCSCF) or
  // for RESTRICTED_DOCC (if an MCSCF);
  // useful if the user wants to just say "freeze_core = true"
  int *core_guess = init_int_array(CalcInfo_->nirreps);
  for (int h=0; h<CalcInfo_->nirreps; h++) {
    core_guess[h] = reference_wavefunction_->frzcpi()[h];
  }

  // This routine sets all orbital subspace arrays properly given
  // some minimal starting information and an Options object
  if (!ras_set3(CalcInfo_->nirreps, CalcInfo_->nmo,
                CalcInfo_->orbs_per_irr, CalcInfo_->docc, CalcInfo_->socc,
                CalcInfo_->frozen_docc, CalcInfo_->frozen_uocc,
                CalcInfo_->rstr_docc, CalcInfo_->rstr_uocc,
                CalcInfo_->ras_opi, core_guess, CalcInfo_->reorder, 1,
                (Parameters_->mcscf ? true : false), options_))
  {
    throw PsiException("Error in ras_set3(). Aborting.",__FILE__,__LINE__);
  }
  // We're done with the core_guess array now
  free(core_guess);

  for (int h=0; h<CalcInfo_->nirreps; h++) {
    CalcInfo_->dropped_docc[h] = CalcInfo_->frozen_docc[h]+CalcInfo_->rstr_docc[h];
    CalcInfo_->dropped_uocc[h] = CalcInfo_->frozen_uocc[h]+CalcInfo_->rstr_uocc[h];
  }

  // Do we have any "restricted" orbitals that will need to be filtered
  // out when we read the two-electron integrals?
  Parameters_->filter_ints = 0;
  for (i=0; i<CalcInfo_->nirreps; i++) {
    if ((CalcInfo_->rstr_docc[i] > 0) || (CalcInfo_->rstr_uocc[i] > 0)) {
      Parameters_->filter_ints = 1;
    }
  }
  if (Parameters_->dertype != "NONE" || Parameters_->mcscf == 1) {
    Parameters_->filter_ints = 1;
  }

  // construct the "ordering" array, which maps the other direction
  // i.e., from a CI orbital to a Pitzer orbital
  CalcInfo_->order = init_int_array(CalcInfo_->nmo);
  for (i=0; i<CalcInfo_->nmo; i++) {
    j = CalcInfo_->reorder[i];
    CalcInfo_->order[j] = i;
  }

  if (Parameters_->print_lvl > 4) {
    outfile->Printf( "\nReordering array = \n");
    for (i=0; i<CalcInfo_->nmo; i++) {
      outfile->Printf( "%3d ", CalcInfo_->reorder[i]);
    }
    outfile->Printf( "\n");
  }

  CalcInfo_->nmotri = (CalcInfo_->nmo * (CalcInfo_->nmo + 1)) / 2 ;

  /* transform orbsym vector to new MO order */
  CalcInfo_->orbsym = init_int_array(CalcInfo_->nmo);
  CalcInfo_->scfeigval = init_array(CalcInfo_->nmo);
  if(Parameters_->zaptn) {
    CalcInfo_->scfeigvala = init_array(CalcInfo_->nmo);
    CalcInfo_->scfeigvalb = init_array(CalcInfo_->nmo);
  }

  for (i=0,cnt=0; i<CalcInfo_->nirreps; i++) {
    for (j=0; j<CalcInfo_->orbs_per_irr[i]; j++,cnt++) {
      k = CalcInfo_->reorder[cnt];
      CalcInfo_->orbsym[k] = i;
    }
  }

  for (i=0; i<CalcInfo_->nmo; i++) {
    j = CalcInfo_->reorder[i];
    CalcInfo_->scfeigval[j] = eig_unsrt[i];
    if(Parameters_->zaptn) {
      CalcInfo_->scfeigvala[j] = eig_unsrt[i];
      CalcInfo_->scfeigvalb[j] = eig_unsrt[i];
    }
  }
  free(eig_unsrt);

  // calculate number of electrons
  CalcInfo_->num_alp = CalcInfo_->num_bet = CalcInfo_->spab = 0;
  if (Parameters_->opentype == PARM_OPENTYPE_NONE ||
      Parameters_->opentype == PARM_OPENTYPE_HIGHSPIN) {
     for (i=0; i<CalcInfo_->nirreps; i++) {
       CalcInfo_->num_alp += CalcInfo_->docc[i] + CalcInfo_->socc[i];
       CalcInfo_->num_bet += CalcInfo_->docc[i];
     }
  }
  else if (Parameters_->opentype == PARM_OPENTYPE_SINGLET) {
    for (i=0; i<CalcInfo_->nirreps; i++) { /* closed-shell part */
      CalcInfo_->spab += CalcInfo_->socc[i];
      CalcInfo_->num_alp += CalcInfo_->docc[i];
      CalcInfo_->num_bet += CalcInfo_->docc[i];
    }
    if (CalcInfo_->spab % 2) {
      throw PsiException("For opentype=singlet must have even number of socc electrons!.",__FILE__,__LINE__);
    }
    CalcInfo_->spab /= 2;
    tmp = 0;
    for (i=0; i<CalcInfo_->nirreps; i++) {
       j = CalcInfo_->socc[i];
       k = 0;
       while (k < j) {
         if (tmp < CalcInfo_->spab) {
           CalcInfo_->num_alp++;
           tmp++;
           k++;
         }
         else {
           CalcInfo_->num_bet++;
           tmp++;
           k++;
         }
       }
     }
  }
  else {
    std::string str = "(get_mo_info): Can't handle opentype = ";
    str += boost::lexical_cast<std::string>( Parameters_->opentype) ;
    throw PsiException(str,__FILE__,__LINE__);
  }

  // Count up number of frozen and restricted core and virtuals
  // These are all redone with new naming scheme and distinction
  // between frozen and restricted as of CDS 4/15
  //
  // For convenience define CalcInfo_->num_drc_orbs as the overall
  // total number of dropped (implicit) core, sum of num_fzc_orbs and
  // num_rsc_orbs.  Replaces most previous instances of CalcInfo_->num_fzc_orbs.
  //
  // Likewise num_drv_orbs is the total number of dropped virtuals and
  // is the sum of num_fzv_orbs and num_rsv_orbs, and replaces most
  // instances of the previous num_fzv_orbs.

  CalcInfo_->num_expl_cor_orbs = 0; // this isn't enabled anymore, zero it

  CalcInfo_->num_fzc_orbs = 0; // truly frozen core, sum(frozen_docc)
  CalcInfo_->num_rsc_orbs = 0; // restricted core, sum(rstr_docc)
  CalcInfo_->num_fzv_orbs = 0; // truly frozen virts, sum(frozen_uocc)
  CalcInfo_->num_rsv_orbs = 0; // restricted virtuals, sum(rstr_uocc)
  CalcInfo_->num_drc_orbs = 0; // dropped core  = num_fzc_orbs + num_rsc_orbs
  CalcInfo_->num_drv_orbs = 0; // dropped virts = num_fzv_orbs + num_rsv_orbs

  for (i=0; i<CalcInfo_->nirreps; i++) {
    CalcInfo_->num_fzc_orbs += CalcInfo_->frozen_docc[i];
    CalcInfo_->num_rsc_orbs += CalcInfo_->rstr_docc[i];
    CalcInfo_->num_fzv_orbs += CalcInfo_->frozen_uocc[i];
    CalcInfo_->num_rsv_orbs += CalcInfo_->rstr_uocc[i];
  }

  CalcInfo_->num_drc_orbs = CalcInfo_->num_fzc_orbs + CalcInfo_->num_rsc_orbs;
  CalcInfo_->num_drv_orbs = CalcInfo_->num_fzv_orbs + CalcInfo_->num_rsv_orbs;
  CalcInfo_->num_rot_orbs = CalcInfo_->nmo - CalcInfo_->num_fzc_orbs - CalcInfo_->num_fzv_orbs;
  CalcInfo_->num_ci_orbs = CalcInfo_->nmo - CalcInfo_->num_drc_orbs -
    CalcInfo_->num_drv_orbs;

  if ((CalcInfo_->num_ci_orbs * (CalcInfo_->num_ci_orbs + 1)) / 2 > IOFF_MAX) {
    throw PsiException("error: IOFF_MAX not large enough!",__FILE__,__LINE__);
  }

  CalcInfo_->num_alp_expl = CalcInfo_->num_alp - CalcInfo_->num_drc_orbs;
  CalcInfo_->num_bet_expl = CalcInfo_->num_bet - CalcInfo_->num_drc_orbs;
  CalcInfo_->npop = CalcInfo_->nmo - CalcInfo_->num_drv_orbs;

  // construct the CalcInfo_->ras_orbs array
  cnt = 0;
  for (i=0; i<4; i++) {
    CalcInfo_->ras_orbs[i] = init_int_matrix(CalcInfo_->nirreps,
      CalcInfo_->num_ci_orbs);
    for (irrep=0; irrep<CalcInfo_->nirreps; irrep++) {
      CalcInfo_->ci_orbs[irrep] += CalcInfo_->ras_opi[i][irrep];
      for (j=0; j<CalcInfo_->ras_opi[i][irrep]; j++) {
        CalcInfo_->ras_orbs[i][irrep][j] = cnt++;
      }
    }
  }

  // Construct "active reordering array"
  CalcInfo_->act_reorder = init_int_array(CalcInfo_->num_ci_orbs);
  CalcInfo_->act_order = init_int_array(CalcInfo_->num_ci_orbs);
  for (int h=0, target=0, pos=0; h<nirrep_; h++){
    target += CalcInfo_->dropped_docc[h];
    for (int i=0; i<CalcInfo_->ci_orbs[h]; i++){
      CalcInfo_->act_reorder[pos++] = CalcInfo_->reorder[target++] - CalcInfo_->num_drc_orbs;
    }
    target += CalcInfo_->dropped_uocc[h];
  }
  for (int i=0; i<CalcInfo_->num_ci_orbs; i++){
    CalcInfo_->act_order[CalcInfo_->act_reorder[i]] = i;
  }
  // for (int i=0; i<CalcInfo_->num_ci_orbs; i++) outfile->Printf("Act_reorder %d \n", CalcInfo_->act_reorder[i]);

  // Build arrays for integrals
  int ncitri = (CalcInfo_->num_ci_orbs * (CalcInfo_->num_ci_orbs + 1)) / 2 ;
  CalcInfo_->onel_ints = (double *) init_array(ncitri);
  CalcInfo_->twoel_ints = (double *) init_array(ncitri * (ncitri + 1) / 2);
  CalcInfo_->maxK = (double *) init_array(CalcInfo_->num_ci_orbs);
  CalcInfo_->gmat = init_matrix(CalcInfo_->num_ci_orbs, CalcInfo_->num_ci_orbs);
  CalcInfo_->tf_onel_ints = init_array(ncitri) ;

} // end get_mo_info()

}} // namespace psi::detci

