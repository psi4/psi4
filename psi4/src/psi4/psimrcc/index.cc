/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
#include <iostream>
#include <algorithm>
#include <cstdio>

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include "index.h"

namespace psi {

namespace psimrcc {

CCIndex::CCIndex(std::shared_ptr<PSIMRCCWfn> wfn, std::string str)
    : label(str), wfn_(wfn), nelements(0), greater_than_or_equal(false), greater_than(false), ntuples(0) {
    nirreps = wfn_->nirrep();
    init();
}

CCIndex::~CCIndex() {}

void CCIndex::init() {
    // New orbital spaces must be added here
    for (size_t i = 0; i < label.size(); ++i)
        if (label[i] == 'o' || label[i] == 'a' || label[i] == 'v' || label[i] == 's' || label[i] == 'n' ||
            label[i] == 'f')
            nelements++;

    // Get the orbital spaces data pointers
    for (size_t i = 0; i < label.size(); ++i) {
        if (label[i] == 'o') {
            mospi.push_back(wfn_->moinfo()->get_occ());
            indices_to_pitzer.push_back(wfn_->moinfo()->get_occ_to_mo());
        } else if (label[i] == 'v') {
            mospi.push_back(wfn_->moinfo()->get_vir());
            indices_to_pitzer.push_back(wfn_->moinfo()->get_vir_to_mo());
        } else if (label[i] == 'a') {
            mospi.push_back(wfn_->moinfo()->get_actv());
            indices_to_pitzer.push_back(wfn_->moinfo()->get_actv_to_mo());
        } else if (label[i] == 'f') {
            mospi.push_back(wfn_->moinfo()->get_fvir());
            indices_to_pitzer.push_back(wfn_->moinfo()->get_fvir_to_mo());
        } else if (label[i] == 's') {
            mospi.push_back(wfn_->moinfo()->sopi_ref());
        } else if (label[i] == 'n') {
            mospi.push_back(wfn_->moinfo()->get_mopi());
        }
    }
    for (int i = 0; i < nelements; i++) {
        first_mos.push_back(std::vector<int>(nirreps, 0));
        dimension.push_back(0);
    }
    for (int i = 0; i < nelements; i++) {
        for (int h = 0; h < nirreps; h++) {
            first_mos[i][h] = dimension[i];
            dimension[i] += mospi[i][h];
        }
    }
    // Set the irrep of each individual element
    element_irrep = std::vector<std::vector<int>>(nelements);
    for (int i = 0; i < nelements; ++i) {
        element_irrep[i] = std::vector<int>(dimension[i], 0);
        int j_abs = 0;
        for (int h = 0; h < nirreps; ++h)
            for (int j = 0; j < mospi[i][h]; ++j) {
                element_irrep[i][j_abs] = h;
                j_abs++;
            }
    }

    switch (nelements) {
        case 0:
            make_zero_index();
            break;
        case 1:
            make_one_index();
            break;
        case 2:
            make_two_index();
            break;
        case 3:
            make_three_index();
            break;
        default: {
            throw std::logic_error("The CCIndex class cannot handle " + label +
                                   "  because there are more than three indices!\n");
        }
    }
}

void CCIndex::make_zero_index() {
    std::vector<std::vector<short>> pairs;  // The pairs ordered as a vector
    ntuples = 0;
    for (int h = 0; h < nirreps; h++) {
        first.push_back(ntuples);
        if (h == 0) {
            std::vector<short> pair;
            pairs.push_back(pair);
            ntuples++;
        }
        last.push_back(ntuples);
        tuplespi.push_back(last[h] - first[h]);
    }
    // Allocate the memory for the tuples and store them
    tuples = std::vector<std::array<short, 3>>{{0, 0, 0}};
}

void CCIndex::make_one_index() {
    // The pairs ordered as a vector
    std::vector<std::vector<short>> pairs;

    // Allocate the 1->tuple mapping array and set them to -1
    one_index_to_tuple_rel_index = std::vector<size_t>(dimension[0], 0);
    one_index_to_irrep = std::vector<int>(dimension[0], 0);

    for (size_t i = 0; i < dimension[0]; ++i) {
        one_index_to_tuple_rel_index[i] = 0;
        one_index_to_irrep[i] = -1;
    }

    ntuples = 0;
    for (int h = 0; h < nirreps; ++h) {
        first.push_back(ntuples);
        for (int p = 0; p < mospi[0][h]; ++p) {
            one_index_to_tuple_rel_index[ntuples] = p;
            one_index_to_irrep[ntuples] = h;
            std::vector<short> pair;
            pair.push_back(ntuples);
            pairs.push_back(pair);
            ntuples++;
        }
        last.push_back(ntuples);
        tuplespi.push_back(last[h] - first[h]);
    }

    tuples = std::vector<std::array<short, 3>>(ntuples, {0, 0, 0});
    for (size_t n = 0; n < pairs.size(); ++n) tuples[n][0] = pairs[n][0];
}

void CCIndex::make_two_index() {
    std::vector<std::vector<short>> pairs;  // The pairs ordered as a vector

    // Allocate the 2->tuple mapping array and set them to -1
    two_index_to_tuple_rel_index = std::vector<std::vector<size_t>>(dimension[0], std::vector<size_t>(dimension[1], 0));
    two_index_to_irrep = std::vector<std::vector<int>>(dimension[0], std::vector<int>(dimension[1], -1));

    // [X>=Y]
    if (label.find(">=") != std::string::npos) {
        greater_than_or_equal = true;
        ntuples = 0;
        for (int h = 0; h < nirreps; h++) {
            first.push_back(ntuples);
            for (int p_sym = 0; p_sym < nirreps; p_sym++) {
                int q_sym = h ^ p_sym;
                int p = first_mos[0][p_sym];
                for (int p_rel = 0; p_rel < mospi[0][p_sym]; p_rel++) {
                    int q = first_mos[1][q_sym];
                    for (int q_rel = 0; q_rel < mospi[1][q_sym]; q_rel++) {
                        if (p >= q) {
                            two_index_to_tuple_rel_index[p][q] = ntuples - first[h];
                            two_index_to_irrep[p][q] = h;
                            std::vector<short> pair;
                            pair.push_back(p);
                            pair.push_back(q);
                            pairs.push_back(pair);
                            ntuples++;
                        }
                        q++;
                    }
                    p++;
                }
            }
            last.push_back(ntuples);
            tuplespi.push_back(last[h] - first[h]);
        }
    } else if (label.find(">") != std::string::npos) {
        greater_than = true;
        ntuples = 0;
        for (int h = 0; h < nirreps; h++) {
            first.push_back(ntuples);
            for (int p_sym = 0; p_sym < nirreps; p_sym++) {
                int q_sym = h ^ p_sym;
                int p = first_mos[0][p_sym];
                for (int p_rel = 0; p_rel < mospi[0][p_sym]; p_rel++) {
                    int q = first_mos[1][q_sym];
                    for (int q_rel = 0; q_rel < mospi[1][q_sym]; q_rel++) {
                        if (p > q) {
                            two_index_to_tuple_rel_index[p][q] = ntuples - first[h];
                            two_index_to_irrep[p][q] = h;
                            std::vector<short> pair;
                            pair.push_back(p);
                            pair.push_back(q);
                            pairs.push_back(pair);
                            ntuples++;
                        }
                        q++;
                    }
                    p++;
                }
            }
            last.push_back(ntuples);
            tuplespi.push_back(last[h] - first[h]);
        }
    } else {
        ntuples = 0;
        for (int h = 0; h < nirreps; h++) {
            first.push_back(ntuples);
            for (int p_sym = 0; p_sym < nirreps; p_sym++) {
                int q_sym = h ^ p_sym;
                int p = first_mos[0][p_sym];
                for (int p_rel = 0; p_rel < mospi[0][p_sym]; p_rel++) {
                    int q = first_mos[1][q_sym];
                    for (int q_rel = 0; q_rel < mospi[1][q_sym]; q_rel++) {
                        two_index_to_tuple_rel_index[p][q] = ntuples - first[h];
                        two_index_to_irrep[p][q] = h;
                        std::vector<short> pair;
                        pair.push_back(p);
                        pair.push_back(q);
                        pairs.push_back(pair);
                        ntuples++;
                        q++;
                    }
                    p++;
                }
            }
            last.push_back(ntuples);
            tuplespi.push_back(last[h] - first[h]);
        }
    }

    // Allocate the memory for the tuples and store them
    tuples = std::vector<std::array<short, 3>>(ntuples, {0, 0, 0});
    for (size_t n = 0; n < pairs.size(); ++n) {
        tuples[n][0] = pairs[n][0];
        tuples[n][1] = pairs[n][1];
    }
}

void CCIndex::make_three_index() {
    if (label.find(">") != std::string::npos) {
        throw std::logic_error("The CCIndex class cannot handle restricted loops for triplets!\n");
    }

    std::vector<std::vector<short>> pairs;  // The pairs ordered as a vector

    // Allocate the 3->tuple mapping array and set them to -1
    three_index_to_tuple_rel_index = std::vector<std::vector<std::vector<size_t>>>(
        dimension[0], std::vector<std::vector<size_t>>(dimension[1], std::vector<size_t>(dimension[2], 0)));
    three_index_to_irrep = std::vector<std::vector<std::vector<int>>>(
        dimension[0], std::vector<std::vector<int>>(dimension[1], std::vector<int>(dimension[2], -1)));

    if (label[2] == '>' && label[4] == '>') {
        ntuples = 0;
        for (int h = 0; h < nirreps; h++) {
            first.push_back(ntuples);
            for (int p_sym = 0; p_sym < nirreps; p_sym++) {
                for (int q_sym = 0; q_sym < nirreps; q_sym++) {
                    int r_sym = h ^ p_sym ^ q_sym;
                    int p = first_mos[0][p_sym];
                    for (int p_rel = 0; p_rel < mospi[0][p_sym]; p_rel++) {
                        int q = first_mos[1][q_sym];
                        for (int q_rel = 0; q_rel < mospi[1][q_sym]; q_rel++) {
                            int r = first_mos[2][r_sym];
                            for (int r_rel = 0; r_rel < mospi[2][r_sym]; r_rel++) {
                                if ((p > q) && (q > r)) {
                                    three_index_to_tuple_rel_index[p][q][r] = ntuples - first[h];
                                    three_index_to_irrep[p][q][r] = h;
                                    std::vector<short> pair;
                                    pair.push_back(p);
                                    pair.push_back(q);
                                    pair.push_back(r);
                                    pairs.push_back(pair);
                                    ntuples++;
                                }
                                r++;
                            }
                            q++;
                        }
                        p++;
                    }
                }
            }
            last.push_back(ntuples);
            tuplespi.push_back(last[h] - first[h]);
        }
    } else {
        ntuples = 0;
        for (int h = 0; h < nirreps; h++) {
            first.push_back(ntuples);
            for (int p_sym = 0; p_sym < nirreps; p_sym++) {
                for (int q_sym = 0; q_sym < nirreps; q_sym++) {
                    int r_sym = h ^ p_sym ^ q_sym;
                    int p = first_mos[0][p_sym];
                    for (int p_rel = 0; p_rel < mospi[0][p_sym]; p_rel++) {
                        int q = first_mos[1][q_sym];
                        for (int q_rel = 0; q_rel < mospi[1][q_sym]; q_rel++) {
                            int r = first_mos[2][r_sym];
                            for (int r_rel = 0; r_rel < mospi[2][r_sym]; r_rel++) {
                                three_index_to_tuple_rel_index[p][q][r] = ntuples - first[h];
                                three_index_to_irrep[p][q][r] = h;
                                std::vector<short> pair;
                                pair.push_back(p);
                                pair.push_back(q);
                                pair.push_back(r);
                                pairs.push_back(pair);
                                ntuples++;
                                r++;
                            }
                            q++;
                        }
                        p++;
                    }
                }
            }
            last.push_back(ntuples);
            tuplespi.push_back(last[h] - first[h]);
        }
    }

    // Allocate the memory for the tuples and store them
    tuples = std::vector<std::array<short, 3>>(ntuples, {0, 0, 0});
    for (size_t n = 0; n < pairs.size(); ++n) {
        tuples[n][0] = pairs[n][0];
        tuples[n][1] = pairs[n][1];
        tuples[n][2] = pairs[n][2];
    }
}

void CCIndex::print() {
    outfile->Printf("\n\n---------------------------------");
    outfile->Printf("\n\tPair Type %s has %lu elements", label.c_str(), (size_t)ntuples);
    outfile->Printf("\n---------------------------------");
    int index = 0;
    for (int h = 0; h < nirreps; h++) {
        if (tuplespi[h] > 0) outfile->Printf("\n\t%s", wfn_->moinfo()->get_irr_lab(h).c_str());
        for (size_t tuple = 0; tuple < tuplespi[h]; ++tuple) {
            outfile->Printf("\n\t\t( ");
            for (int k = 0; k < nelements; k++) outfile->Printf("%d ", tuples[index][k]);
            outfile->Printf(")");
            index++;
        }
    }
    outfile->Printf("\n---------------------------------");
}

}  // namespace psimrcc
}  // namespace psi
