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

#include "occwave.h"

using namespace psi;

namespace psi {
namespace occwave {

void OCCWave::idp2() {
    int dim;

    if (reference_ == "RESTRICTED") {
        // Form IDPs
        nidpA = 0;

        // V-O: I exclude symmetry broken rotations from the list of IDPs since they already have zero gradient.
        for (int h = 0; h < nirrep_; h++) {
            nidpA += virtpiA[h] * occpiA[h];
        }

        outfile->Printf("\tNumber of independent-pairs: %3d\n", nidpA);

        if (nidpA != 0) {
            wogA = new Array1d("Alpha MO grad vector", nidpA);
            wogA->zero();
        }

        // allocate memory
        idprowA = new int[nidpA];
        idpcolA = new int[nidpA];
        idpirrA = new int[nidpA];

        // initialize
        memset(idprowA, 0, sizeof(int) * nidpA);
        memset(idpcolA, 0, sizeof(int) * nidpA);
        memset(idpirrA, 0, sizeof(int) * nidpA);

        // set idpA
        dim = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (int a = 0; a < virtpiA[h]; a++) {
                for (int i = 0; i < occpiA[h]; i++) {
                    idprowA[dim] = a;
                    idpcolA[dim] = i;
                    idpirrA[dim] = h;
                    dim++;
                }
            }
        }

        if (print_ > 2) {
            for (int i = 0; i < nidpA; i++) {
                outfile->Printf("\n i, idpirrA, idprowA, idpcolA: %3d %3d %3d %3d\n", i, idpirrA[i], idprowA[i],
                                idpcolA[i]);
            }
        }

    }  // end if (reference_ == "RESTRICTED")

    else if (reference_ == "UNRESTRICTED") {
        // Form IDPs
        nidpA = 0;
        nidpB = 0;

        // V-O: I exclude symmetry broken rotations from the list of IDPs since they already have zero gradient.
        for (int h = 0; h < nirrep_; h++) {
            nidpA += virtpiA[h] * occpiA[h];
            nidpB += virtpiB[h] * occpiB[h];
        }

        outfile->Printf("\tNumber of alpha independent-pairs:%3d\n", nidpA);
        outfile->Printf("\tNumber of beta independent-pairs :%3d\n", nidpB);

        if (nidpA != 0) {
            idp_returnA = 1;
            wogA = new Array1d("Alpha MO grad vector", nidpA);
            wogA->zero();
        }

        if (nidpB != 0) {
            idp_returnB = 1;
            wogB = new Array1d("Beta MO grad vector", nidpB);
            wogB->zero();
        }

        // allocate memory
        idprowA = new int[nidpA];
        idpcolA = new int[nidpA];
        idpirrA = new int[nidpA];
        idprowB = new int[nidpB];
        idpcolB = new int[nidpB];
        idpirrB = new int[nidpB];

        // initialize
        memset(idprowA, 0, sizeof(int) * nidpA);
        memset(idpcolA, 0, sizeof(int) * nidpA);
        memset(idpirrA, 0, sizeof(int) * nidpA);
        memset(idprowB, 0, sizeof(int) * nidpB);
        memset(idpcolB, 0, sizeof(int) * nidpB);
        memset(idpirrB, 0, sizeof(int) * nidpB);

        // set idpA
        dim = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (int a = 0; a < virtpiA[h]; a++) {
                for (int i = 0; i < occpiA[h]; i++) {
                    idprowA[dim] = a;
                    idpcolA[dim] = i;
                    idpirrA[dim] = h;
                    dim++;
                }
            }
        }

        // set idpB
        dim = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (int a = 0; a < virtpiB[h]; a++) {
                for (int i = 0; i < occpiB[h]; i++) {
                    idprowB[dim] = a;
                    idpcolB[dim] = i;
                    idpirrB[dim] = h;
                    dim++;
                }
            }
        }

        if (print_ > 2) {
            for (int i = 0; i < nidpA; i++) {
                outfile->Printf("\n i, idpirrA, idprowA, idpcolA: %3d %3d %3d %3d\n", i, idpirrA[i], idprowA[i],
                                idpcolA[i]);
            }

            for (int i = 0; i < nidpB; i++) {
                outfile->Printf("\n i, idpirrB, idprowB, idpcolB: %3d %3d %3d %3d\n", i, idpirrB[i], idprowB[i],
                                idpcolB[i]);
            }
        }

    }  // end if (reference_ == "UNRESTRICTED")
}  // end of main
}
}  // End Namespaces
