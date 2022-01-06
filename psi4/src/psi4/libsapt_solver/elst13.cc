/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "sapt2p3.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
namespace sapt {

void SAPT2p3::elst13() {
    double e_elst130 = elst130(wBAA_, wBRR_, CHFA_, PSIF_SAPT_AMPS, "qAA Density Matrix", "qRR Density Matrix",
                               "Y3 AR Amplitudes", foccA_, noccA_, nvirA_);

    if (debug_) {
        outfile->Printf("    Elst130,r           = %18.12lf [Eh]\n", e_elst130);
    }

    double e_elst103 = elst130(wABB_, wASS_, CHFB_, PSIF_SAPT_AMPS, "qBB Density Matrix", "qSS Density Matrix",
                               "Y3 BS Amplitudes", foccB_, noccB_, nvirB_);

    if (debug_) {
        outfile->Printf("    Elst103,r           = %18.12lf [Eh]\n\n", e_elst103);
    }

    e_elst13_ = e_elst130 + e_elst103;

    if (print_) {
        outfile->Printf("    Elst13,r            = %18.12lf [Eh]\n", e_elst13_);
    }
}

double SAPT2p3::elst130(double **wBAA, double **wBRR, double **CHFA, int ampfile, const char *pAAlabel,
                        const char *pRRlabel, const char *Ylabel, size_t foccA, size_t noccA, size_t nvirA) {
    size_t aoccA = noccA - foccA;

    double **pAA = block_matrix(aoccA, aoccA);
    psio_->read_entry(ampfile, pAAlabel, (char *)pAA[0], sizeof(double) * aoccA * aoccA);

    double **pRR = block_matrix(nvirA, nvirA);
    psio_->read_entry(ampfile, pRRlabel, (char *)pRR[0], sizeof(double) * nvirA * nvirA);

    double **yAR = block_matrix(aoccA, nvirA);
    psio_->read_entry(ampfile, Ylabel, (char *)yAR[0], sizeof(double) * aoccA * nvirA);

    double e1 = 0.0, e2 = 0.0, e3 = 0.0;

    for (int a = 0; a < aoccA; a++) {
        e1 -= 4.0 * C_DDOT(aoccA, pAA[a], 1, &(wBAA[a + foccA][foccA]), 1);
    }

    e2 += 4.0 * C_DDOT(nvirA * nvirA, pRR[0], 1, wBRR[0], 1);
    e3 += 4.0 * C_DDOT(aoccA * nvirA, yAR[0], 1, CHFA[foccA], 1);

    free_block(pAA);
    free_block(pRR);
    free_block(yAR);

    if (debug_) {
        outfile->Printf("\n    Elst13_1            = %18.12lf [Eh]\n", e1);
        outfile->Printf("    Elst13_2            = %18.12lf [Eh]\n", e2);
        outfile->Printf("    Elst13_3            = %18.12lf [Eh]\n", e3);
    }

    return (e1 + e2 + e3);
}
}  // namespace sapt
}  // namespace psi
