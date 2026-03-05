/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*!
  \file
  \brief Obtain the QT orbital reordering array between Pitzer and correlated
    order
  \ingroup QT
*/

#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <cstdio>
#include <cstdlib>

namespace psi {

/*!
** reorder_qt()
**
** This function constructs a reordering array according to the
** "Quantum Trio" standard ordering, in which the orbitals are divided
** into the following sets: frozen core, then doubly occupied, then singly
** occupied, then virtuals, then deleted (frozen) virtuals.
** The reordering array takes a basis function in
** Pitzer ordering (orbitals grouped according to irrep) and gives the
** corresponding index in the Quantum Trio numbering scheme.
**
** Should give the same reordering array as in the old libread30 routines.
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia, 1995
**
** \param docc_in        = doubly occupied orbitals per irrep
** \param socc_in        = singly occupied orbitals per irrep
** \param frozen_docc_in = frozen occupied orbitals per irrep
** \param frozen_uocc_in = frozen unoccupied orbitals per irrep
** \param order          = reordering array (Pitzer->QT order)
**
** \ingroup QT
*/
PSI_API void reorder_qt(Dimension const& docc_in, Dimension const& socc_in, Dimension const& frozen_docc_in,
                        Dimension const& frozen_uocc_in, int *order, Dimension const& orbs_per_irrep) {
    int cnt = 0, irrep, point, tmpi;
    int *used, *offset;
    int *docc, *socc, *frozen_docc, *frozen_uocc;
    int *uocc;

    auto nirreps = docc_in.n();
    if (!(nirreps == socc_in.n() && nirreps == frozen_docc_in.n() &&
        nirreps == frozen_uocc_in.n() && nirreps == orbs_per_irrep.n())) {
        throw PSIEXCEPTION("Inconsistent number of irreps in Dimensions passed to reorder_qt.");
    }
    used = init_int_array(nirreps);
    offset = init_int_array(nirreps);

    docc = init_int_array(nirreps);
    socc = init_int_array(nirreps);
    frozen_docc = init_int_array(nirreps);
    frozen_uocc = init_int_array(nirreps);
    uocc = init_int_array(nirreps);

    for (irrep = 0; irrep < nirreps; irrep++) {
        docc[irrep] = docc_in[irrep];
        socc[irrep] = socc_in[irrep];
        frozen_docc[irrep] = frozen_docc_in[irrep];
        frozen_uocc[irrep] = frozen_uocc_in[irrep];
    }

    /* construct the offset array */
    offset[0] = 0;
    for (irrep = 1; irrep < nirreps; irrep++) {
        offset[irrep] = offset[irrep - 1] + orbs_per_irrep[irrep - 1];
    }

    /* construct the uocc array */
    for (irrep = 0; irrep < nirreps; irrep++) {
        tmpi = frozen_uocc[irrep] + docc[irrep] + socc[irrep];
        if (tmpi > orbs_per_irrep[irrep]) {
            outfile->Printf("(reorder_qt): orbitals don't add up for irrep %d\n", irrep);
            return;
        } else
            uocc[irrep] = orbs_per_irrep[irrep] - tmpi;
    }

    /* do the frozen core */
    for (irrep = 0; irrep < nirreps; irrep++) {
        while (frozen_docc[irrep]) {
            point = used[irrep] + offset[irrep];
            order[point] = cnt++;
            used[irrep]++;
            frozen_docc[irrep]--;
            docc[irrep]--;
        }
    }

    /* do doubly occupied orbitals */
    for (irrep = 0; irrep < nirreps; irrep++) {
        while (docc[irrep]) {
            point = used[irrep] + offset[irrep];
            order[point] = cnt++;
            used[irrep]++;
            docc[irrep]--;
        }
    }

    /* do singly-occupied orbitals */
    for (irrep = 0; irrep < nirreps; irrep++) {
        while (socc[irrep]) {
            point = used[irrep] + offset[irrep];
            order[point] = cnt++;
            used[irrep]++;
            socc[irrep]--;
        }
    }

    /* do virtual orbitals */
    for (irrep = 0; irrep < nirreps; irrep++) {
        while (uocc[irrep]) {
            point = used[irrep] + offset[irrep];
            order[point] = cnt++;
            used[irrep]++;
            uocc[irrep]--;
        }
    }

    /* do frozen uocc */
    for (irrep = 0; irrep < nirreps; irrep++) {
        while (frozen_uocc[irrep]) {
            point = used[irrep] + offset[irrep];
            order[point] = cnt++;
            used[irrep]++;
            frozen_uocc[irrep]--;
        }
    }

    /* do a final check */
    for (irrep = 0; irrep < nirreps; irrep++) {
        if (used[irrep] > orbs_per_irrep[irrep]) {
            outfile->Printf("(reorder_qt): on final check, used more orbitals");
            outfile->Printf("   than were available (%d vs %d) for irrep %d\n", used[irrep], orbs_per_irrep[irrep],
                            irrep);
        }
    }

    free(used);
    free(offset);
    free(docc);
    free(socc);
    free(frozen_docc);
    free(frozen_uocc);
    free(uocc);
}

/*!
** reorder_qt_uhf()
**
** Generalization of reorder_qt() for UHF case
**
** \param docc        = doubly occupied orbitals per irrep
** \param socc        = singly occupied orbitals per irrep
** \param frozen_docc = frozen occupied orbitals per irrep
** \param frozen_uocc = frozen unoccupied orbitals per irrep
** \param order_alpha = reordering array for alpha (Pitzer->QT order)
** \param order_beta  = reordering array for beta  (Pitzer->QT order)
**
** \ingroup QT
*/
PSI_API void reorder_qt_uhf(Dimension const& docc, Dimension const& socc, Dimension const& frozen_docc,
                            Dimension const& frozen_uocc, int *order_alpha, int *order_beta, Dimension const& orbspi) {
    int p, nmo;
    int cnt_alpha, cnt_beta, irrep, tmpi;
    int *offset, this_offset;
    int *uocc;

    auto nirreps = docc.n();
    if (!(nirreps == socc.n() && nirreps == frozen_docc.n() &&
        nirreps == frozen_uocc.n() && nirreps == orbspi.n())) {
        throw PSIEXCEPTION("Inconsistent number of irreps in Dimensions passed to reorder_qt_uhf.");
    }

    Dimension nalphapi(nirreps, "Number of alpha electrons per irrep");
    Dimension nbetapi(nirreps, "Number of beta electrons per irrep");
    for (int h = 0; h < nirreps; h++) {
        nalphapi[h] = docc[h] + socc[h];
        nbetapi[h] = docc[h];
    }

    offset = init_int_array(nirreps);

    uocc = init_int_array(nirreps);

    /* construct the offset array */
    offset[0] = 0;
    for (irrep = 1; irrep < nirreps; irrep++) {
        offset[irrep] = offset[irrep - 1] + orbspi[irrep - 1];
    }

    /* construct the uocc array */
    nmo = 0;
    for (irrep = 0; irrep < nirreps; irrep++) {
        nmo += orbspi[irrep];
        tmpi = frozen_uocc[irrep] + docc[irrep] + socc[irrep];
        if (tmpi > orbspi[irrep]) {
            outfile->Printf("(reorder_qt_uhf): orbitals don't add up for irrep %d\n", irrep);
            return;
        }

        /* construct the uocc array */
        nmo = 0;
        for (irrep = 0; irrep < nirreps; irrep++) {
            nmo += orbspi[irrep];
            tmpi = frozen_uocc[irrep] + docc[irrep] + socc[irrep];
            if (tmpi > orbspi[irrep]) {
                outfile->Printf("(reorder_qt_uhf): orbitals don't add up for irrep %d\n", irrep);
                return;
            } else
                uocc[irrep] = orbspi[irrep] - tmpi;
        }

        cnt_alpha = cnt_beta = 0;

        /* do the frozen core */
        for (irrep = 0; irrep < nirreps; irrep++) {
            this_offset = offset[irrep];
            for (p = 0; p < frozen_docc[irrep]; p++) {
                order_alpha[this_offset + p] = cnt_alpha++;
                order_beta[this_offset + p] = cnt_beta++;
            }
        }

        /* alpha occupied orbitals */
        for (irrep = 0; irrep < nirreps; irrep++) {
            this_offset = offset[irrep] + frozen_docc[irrep];
            for (p = 0; p < nalphapi[irrep] - frozen_docc[irrep]; p++) {
                order_alpha[this_offset + p] = cnt_alpha++;
            }
        }

        /* beta occupied orbitals */
        for (irrep = 0; irrep < nirreps; irrep++) {
            this_offset = offset[irrep] + frozen_docc[irrep];
            for (p = 0; p < nbetapi[irrep] - frozen_docc[irrep]; p++) {
                order_beta[this_offset + p] = cnt_beta++;
            }
        }

        /* alpha unoccupied orbitals */
        for (irrep = 0; irrep < nirreps; irrep++) {
            this_offset = offset[irrep] + nalphapi[irrep];
            for (p = 0; p < orbspi[irrep] - nalphapi[irrep] - frozen_uocc[irrep]; p++) {
                order_alpha[this_offset + p] = cnt_alpha++;
            }
        }

        /* beta unoccupied orbitals */
        for (irrep = 0; irrep < nirreps; irrep++) {
            this_offset = offset[irrep] + nbetapi[irrep];
            for (p = 0; p < orbspi[irrep] - nbetapi[irrep] - frozen_uocc[irrep]; p++) {
                order_beta[this_offset + p] = cnt_beta++;
            }
        }

        /* do the frozen uocc */
        for (irrep = 0; irrep < nirreps; irrep++) {
            this_offset = offset[irrep] + docc[irrep] + socc[irrep] + uocc[irrep];
            for (p = 0; p < frozen_uocc[irrep]; p++) {
                order_alpha[this_offset + p] = cnt_alpha++;
                order_beta[this_offset + p] = cnt_beta++;
            }
        }

        /* do a final check */
        for (irrep = 0; irrep < nirreps; irrep++) {
            if (cnt_alpha > nmo) {
                outfile->Printf("(reorder_qt_uhf): on final check, used more orbitals");
                outfile->Printf("   than were available (%d vs %d) for irrep %d\n", cnt_alpha, nmo, irrep);
            }
            if (cnt_beta > nmo) {
                outfile->Printf("(reorder_qt_uhf): on final check, used more orbitals");
                outfile->Printf("   than were available (%d vs %d) for irrep %d\n", cnt_beta, nmo, irrep);
            }
        }

        free(offset);
        free(uocc);
    }
}

}  // namespace psi
