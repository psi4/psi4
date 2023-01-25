/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

void OCCWave::idp() {
    int dim;

    if (reference_ == "RESTRICTED") {
        // Form IDPs
        nidpA = 0;

        auto alpha_dim = Dimension(nirrep_, "Alpha independent pairs");
        // V-O: I exclude symmetry broken rotations from the list of IDPs since they already have zero gradient.
        for (int h = 0; h < nirrep_; h++) {
            alpha_dim[h] = virtpiA[h] * occpiA[h];
            nidpA += alpha_dim[h];
        }
        idp_dimensions_ = {{SpinType::Alpha, alpha_dim}};

        outfile->Printf("\n\tNumber of independent-pairs: %3d\n", nidpA);

        if (nidpA != 0) {
            idp_returnA = 1;
            wogA = new Array1d("Alpha MO grad vector", nidpA);
            kappaA = new Array1d("Alpha orb rot params vector of current step", nidpA);
            kappa_newA = new Array1d("Alpha New orb rot params vector of current step", nidpA);
            auto kappa_barA = std::make_shared<Vector>(idp_dimensions_[SpinType::Alpha]);
            kappa_bar_ = {{SpinType::Alpha, kappa_barA}};
            wog_intA = new Array1d("Alpha Interpolated MO grad vector", nidpA);
            wogA->zero();
            kappaA->zero();

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

            idprow_ = {{SpinType::Alpha, idprowA}};
            idpcol_ = {{SpinType::Alpha, idpcolA}};
            idpirr_ = {{SpinType::Alpha, idpirrA}};

            if (print_ > 2) {
                for (int i = 0; i < nidpA; i++) {
                    outfile->Printf("\n i, idpirrA, idprowA, idpcolA: %3d %3d %3d %3d\n", i, idpirrA[i], idprowA[i],
                                    idpcolA[i]);
                }
            }
        }  // end if nidpA != 0

        else if (nidpA == 0) {
            outfile->Printf("\tThere is not any non-redundant orbital rotation pair! \n");
            tstop();
            exit(EXIT_SUCCESS);
        }

    }  // end if (reference_ == "RESTRICTED")

    else if (reference_ == "UNRESTRICTED") {
        // Form IDPs
        nidpA = 0;
        nidpB = 0;

        // V-O: I exclude symmetry broken rotations from the list of IDPs since they already have zero gradient.
        auto alpha_dim = Dimension(nirrep_, "Alpha independent pairs");
        auto beta_dim = Dimension(nirrep_, "Beta independent pairs");

        for (int h = 0; h < nirrep_; h++) {
            alpha_dim[h] = virtpiA[h] * occpiA[h];
            nidpA += alpha_dim[h];
            beta_dim[h] = virtpiB[h] * occpiB[h];
            nidpB += beta_dim[h];
        }
        idp_dimensions_ = {{SpinType::Alpha, alpha_dim}, {SpinType::Beta, beta_dim}};

        outfile->Printf("\n\tNumber of alpha independent-pairs:%3d\n", nidpA);
        outfile->Printf("\tNumber of beta independent-pairs :%3d\n", nidpB);

        if (nidpA == 0 && nidpB == 0) {
            outfile->Printf("\tThere is not any non-redundant orbital rotation pair! \n");
            tstop();
            exit(EXIT_SUCCESS);
        }

        if (nidpA != 0) {
            idp_returnA = 1;
            wogA = new Array1d("Alpha MO grad vector", nidpA);
            kappaA = new Array1d("Alpha orb rot params vector of current step", nidpA);
            kappa_newA = new Array1d("Alpha New orb rot params vector of current step", nidpA);
            wog_intA = new Array1d("Alpha Interpolated MO grad vector", nidpA);
            wogA->zero();
            kappaA->zero();
        }

        if (nidpB != 0) {
            idp_returnB = 1;
            wogB = new Array1d("Beta MO grad vector", nidpB);
            kappaB = new Array1d("Beta orb rot params vector of current step", nidpB);
            kappa_newB = new Array1d("Beta New orb rot params vector of current step", nidpB);
            wog_intB = new Array1d("Beta Interpolated MO grad vector", nidpB);
            wogB->zero();
            kappaB->zero();
        }

        auto kappa_barA = std::make_shared<Vector>(idp_dimensions_[SpinType::Alpha]);
        auto kappa_barB = std::make_shared<Vector>(idp_dimensions_[SpinType::Beta]);
        kappa_bar_ = {{SpinType::Alpha, kappa_barA}, {SpinType::Beta, kappa_barB}};

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

        idprow_ = {{SpinType::Alpha, idprowA}, {SpinType::Beta, idprowB}};
        idpcol_ = {{SpinType::Alpha, idpcolA}, {SpinType::Beta, idpcolB}};
        idpirr_ = {{SpinType::Alpha, idpirrA}, {SpinType::Beta, idpirrB}};

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
