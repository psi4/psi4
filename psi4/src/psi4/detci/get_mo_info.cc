/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*! \file
    \ingroup DETCI
    \brief This file gets information about MOs and subspaces of them
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/vector.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "psi4/detci/ciwave.h"
#include "psi4/libmints/matrix.h"
#include "psi4/detci/structs.h"

namespace psi {
namespace detci {

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
void CIWavefunction::get_mo_info() {

    CalcInfo_->sigma_initialized = 0;

    // Initial guess will overwrite some of this later.
    CalcInfo_->nirreps = nirrep();
    CalcInfo_->nso = nso();
    CalcInfo_->nmo = nmo();

    // DGAS: David does this make sense?
    CalcInfo_->iopen = !same_a_b_orbs();
    CalcInfo_->labels = molecule()->irrep_labels();
    CalcInfo_->docc = doccpi();
    CalcInfo_->socc = soccpi();
    CalcInfo_->enuc = molecule()->nuclear_repulsion_energy();
    CalcInfo_->escf = reference_energy();
    CalcInfo_->edrc = 0.0;

    if (CalcInfo_->iopen && Parameters_->opentype == PARM_OPENTYPE_NONE) {
        outfile->Printf("Warning: iopen=1,opentype=none. Making iopen=0\n");
        CalcInfo_->iopen = 0;
    } else if (!CalcInfo_->iopen &&
               (Parameters_->opentype == PARM_OPENTYPE_HIGHSPIN ||
                Parameters_->opentype == PARM_OPENTYPE_SINGLET)) {
        outfile->Printf("Warning: iopen=0,opentype!=closed. Making iopen=1\n");
        CalcInfo_->iopen = 1;
    }
    if (Parameters_->ref_sym >= CalcInfo_->nirreps) {
        outfile->Printf("Warning: ref_sym >= nirreps.  Setting ref_sym=0\n");
        Parameters_->ref_sym = 0;
    }

    // these are all initialized to zero at this point
    CalcInfo_->reorder.resize(CalcInfo_->nmo);
    CalcInfo_->ras_opi = init_int_matrix(4, CalcInfo_->nirreps);

    CalcInfo_->frozen_docc = Dimension(CalcInfo_->nirreps, "Frozen doubly occupied orbitals");
    CalcInfo_->rstr_docc = Dimension(CalcInfo_->nirreps, "Restricted doubly occupied orbitals");

    CalcInfo_->ci_orbs = Dimension(CalcInfo_->nirreps, "Active orbitals");

    CalcInfo_->frozen_uocc = Dimension(CalcInfo_->nirreps, "Frozen virtual orbitals");
    CalcInfo_->rstr_uocc = Dimension(CalcInfo_->nirreps, "Restricted virtual orbitals");

    Dimension frzcpi = reference_wavefunction_->frzcpi();
    // This routine sets all orbital subspace arrays properly given
    // some minimal starting information and an Options object
    if (!ras_set3(CalcInfo_->nirreps, CalcInfo_->nmo, nmopi_,
                  CalcInfo_->docc, CalcInfo_->socc, CalcInfo_->frozen_docc,
                  CalcInfo_->frozen_uocc, CalcInfo_->rstr_docc,
                  CalcInfo_->rstr_uocc, CalcInfo_->ras_opi,
                  frzcpi,
                  CalcInfo_->reorder.data(), 1,
                  (Parameters_->mcscf ? true : false), options_)) {
        throw PsiException("Error in ras_set3(). Aborting.", __FILE__,
                           __LINE__);
    }

    CalcInfo_->dropped_docc = CalcInfo_->frozen_docc + CalcInfo_->rstr_docc;
    CalcInfo_->dropped_docc.set_name("Dropped occupied orbitals");
    CalcInfo_->dropped_uocc = CalcInfo_->frozen_uocc + CalcInfo_->rstr_uocc;
    CalcInfo_->dropped_uocc.set_name("Dropped virtual orbitals");

    // construct the "ordering" array, which maps the other direction
    // i.e., from a CI orbital to a Pitzer orbital
    CalcInfo_->order.resize(CalcInfo_->nmo);
    for (int i = 0; i < CalcInfo_->nmo; i++) {
        CalcInfo_->order[CalcInfo_->reorder[i]] = i;
    }

    CalcInfo_->nmotri = (CalcInfo_->nmo * (CalcInfo_->nmo + 1)) / 2;

    /* transform orbsym vector to new MO order */
    CalcInfo_->orbsym = init_int_array(CalcInfo_->nmo);
    CalcInfo_->scfeigval.resize(CalcInfo_->nmo);
    if (Parameters_->zaptn) {
        CalcInfo_->scfeigvala.resize(CalcInfo_->nmo);
        CalcInfo_->scfeigvalb.resize(CalcInfo_->nmo);
    }

    for (int i = 0, cnt = 0; i < CalcInfo_->nirreps; i++) {
        for (int j = 0; j < nmopi_[i]; j++, cnt++) {
            CalcInfo_->orbsym[CalcInfo_->reorder[cnt]] = i;
        }
    }

    for (int i = 0; i < CalcInfo_->nmo; i++) {
        int j = CalcInfo_->reorder[i];
        CalcInfo_->scfeigval[j] = epsilon_a()->get(i);
        if (Parameters_->zaptn) {
            CalcInfo_->scfeigvala[j] = epsilon_a()->get(i);
            CalcInfo_->scfeigvalb[j] = epsilon_a()->get(i);
        }
    }

    // calculate number of electrons
    CalcInfo_->num_alp = CalcInfo_->num_bet = CalcInfo_->spab = 0;
    if (Parameters_->opentype == PARM_OPENTYPE_NONE ||
        Parameters_->opentype == PARM_OPENTYPE_HIGHSPIN) {
        CalcInfo_->num_alp += CalcInfo_->docc.sum() + CalcInfo_->socc.sum();
        CalcInfo_->num_bet += CalcInfo_->docc.sum();
    } else if (Parameters_->opentype == PARM_OPENTYPE_SINGLET) {
        CalcInfo_->spab += CalcInfo_->socc.sum();
        CalcInfo_->num_alp += CalcInfo_->docc.sum();
        CalcInfo_->num_bet += CalcInfo_->docc.sum();
        if (CalcInfo_->spab % 2) {
            throw PsiException(
                "For opentype=singlet must have even number of socc "
                "electrons!.",
                __FILE__, __LINE__);
        }
        CalcInfo_->spab /= 2;
        int tmp = 0;
        for (int i = 0; i < CalcInfo_->nirreps; i++) {
            int j = CalcInfo_->socc[i];
            int k = 0;
            while (k < j) {
                if (tmp < CalcInfo_->spab) {
                    CalcInfo_->num_alp++;
                    tmp++;
                    k++;
                } else {
                    CalcInfo_->num_bet++;
                    tmp++;
                    k++;
                }
            }
        }
    } else {
        std::string str = "(get_mo_info): Can't handle opentype = ";
        str += std::to_string(Parameters_->opentype);
        throw PsiException(str, __FILE__, __LINE__);
    }

    // Count up number of frozen and restricted core and virtuals
    // These are all redone with new naming scheme and distinction
    // between frozen and restricted as of CDS 4/15
    //
    // For convenience define CalcInfo_->num_drc_orbs as the overall
    // total number of dropped (implicit) core, sum of num_fzc_orbs and
    // num_rsc_orbs.  Replaces most previous instances of
    // CalcInfo_->num_fzc_orbs.
    //
    // Likewise num_drv_orbs is the total number of dropped virtuals and
    // is the sum of num_fzv_orbs and num_rsv_orbs, and replaces most
    // instances of the previous num_fzv_orbs.

    CalcInfo_->num_expl_cor_orbs = 0;  // this isn't enabled anymore, zero it

    CalcInfo_->num_fzc_orbs = CalcInfo_->frozen_docc.sum();  // truly frozen core, sum(frozen_docc)
    CalcInfo_->num_rsc_orbs = CalcInfo_->rstr_docc.sum();  // restricted core, sum(rstr_docc)
    CalcInfo_->num_fzv_orbs = CalcInfo_->frozen_uocc.sum();  // truly frozen virts, sum(frozen_uocc)
    CalcInfo_->num_rsv_orbs = CalcInfo_->rstr_uocc.sum();  // restricted virtuals, sum(rstr_uocc)

    CalcInfo_->num_drc_orbs = CalcInfo_->num_fzc_orbs + CalcInfo_->num_rsc_orbs;
    CalcInfo_->num_drv_orbs = CalcInfo_->num_fzv_orbs + CalcInfo_->num_rsv_orbs;
    CalcInfo_->num_rot_orbs = CalcInfo_->nmo - CalcInfo_->num_fzc_orbs - CalcInfo_->num_fzv_orbs;
    CalcInfo_->num_ci_orbs = CalcInfo_->nmo - CalcInfo_->num_drc_orbs - CalcInfo_->num_drv_orbs;

    if ((CalcInfo_->num_ci_orbs * (CalcInfo_->num_ci_orbs + 1)) / 2 > IOFF_MAX) {
        throw PsiException("error: IOFF_MAX not large enough!", __FILE__,
                           __LINE__);
    }

    CalcInfo_->num_alp_expl = CalcInfo_->num_alp - CalcInfo_->num_drc_orbs;
    CalcInfo_->num_bet_expl = CalcInfo_->num_bet - CalcInfo_->num_drc_orbs;
    CalcInfo_->npop = CalcInfo_->nmo - CalcInfo_->num_drv_orbs;

    // construct the CalcInfo_->ras_orbs array
    for (int i = 0, cnt = 0; i < 4; i++) {
        CalcInfo_->ras_orbs[i] = init_int_matrix(CalcInfo_->nirreps, CalcInfo_->num_ci_orbs);
        for (int irrep = 0; irrep < CalcInfo_->nirreps; irrep++) {
            CalcInfo_->ci_orbs[irrep] += CalcInfo_->ras_opi[i][irrep];
            for (int j = 0; j < CalcInfo_->ras_opi[i][irrep]; j++) {
                CalcInfo_->ras_orbs[i][irrep][j] = cnt++;
            }
        }
    }

    // Construct "active reordering array"
    CalcInfo_->act_reorder.resize(CalcInfo_->num_ci_orbs);
    CalcInfo_->act_order.resize(CalcInfo_->num_ci_orbs);
    for (int h = 0, target = 0, pos = 0; h < nirrep_; h++) {
        target += CalcInfo_->dropped_docc[h];
        for (int i = 0; i < CalcInfo_->ci_orbs[h]; i++) {
            CalcInfo_->act_reorder[pos++] = CalcInfo_->reorder[target++] - CalcInfo_->num_drc_orbs;
        }
        target += CalcInfo_->dropped_uocc[h];
    }
    for (int i = 0; i < CalcInfo_->num_ci_orbs; i++) {
        CalcInfo_->act_order[CalcInfo_->act_reorder[i]] = i;
    }

    // Build arrays for integrals
    int ncitri = (CalcInfo_->num_ci_orbs * (CalcInfo_->num_ci_orbs + 1)) / 2;
    CalcInfo_->num_ci_tri = ncitri;
    CalcInfo_->num_ci_tri2 = (ncitri * (ncitri + 1)) / 2;
    CalcInfo_->so_onel_ints = SharedMatrix(new Matrix("SO CI One Electron Ints", nso_, nso_));
    CalcInfo_->onel_ints = SharedVector(new Vector("CI One Electron Ints", ncitri));
    CalcInfo_->twoel_ints = SharedVector(new Vector("CI Two Electron Ints", ncitri * (ncitri + 1) / 2));
    CalcInfo_->gmat = SharedVector(new Vector("CI RAS Gmat", CalcInfo_->num_ci_orbs * CalcInfo_->num_ci_orbs));
    CalcInfo_->tf_onel_ints = SharedVector(new Vector("CI TF One Electron Ints", ncitri));

}  // end get_mo_info()
}}  // namespace psi::detci
