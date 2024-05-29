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

/**
 *  @file manybody.cc
 *  @ingroup (PSIMRCC)
 *  @brief The base class for all the many-body methods
 */
#include <iostream>

#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>
#include <vector>

#include "psi4/liboptions/liboptions.h"
#include "psi4/libmoinfo/libmoinfo.h"

#include "algebra_interface.h"
#include "blas.h"
#include "index.h"
#include "manybody.h"
#include "matrix.h"
#include "sort.h"

namespace psi {

namespace psimrcc {

/**
 * Allocate the effective Hamiltonian matrices and eigenvectors
 * @todo wrap the current operations in an init() function
 */
CCManyBody::CCManyBody(std::shared_ptr<PSIMRCCWfn> wfn, Options& options) : wfn_(wfn), options_(options) {
    // Allocate memory for the eigenvector and the effective Hamiltonian
    zeroth_order_eigenvector = std::vector<double>(wfn_->moinfo()->get_nrefs(), 0);
    right_eigenvector = std::vector<double>(wfn_->moinfo()->get_nrefs(), 0);
    left_eigenvector = std::vector<double>(wfn_->moinfo()->get_nrefs(), 0);
    Heff = block_matrix(wfn_->moinfo()->get_nrefs(), wfn_->moinfo()->get_nrefs());
    Heff_mrpt2 = block_matrix(wfn_->moinfo()->get_nrefs(), wfn_->moinfo()->get_nrefs());
    wfn_->free_memory_ -= 2 * sizeof(double) * wfn_->moinfo()->get_nrefs();

    huge = 1.0e100;
    norm_amps = 0.0;
    delta_t1_amps = 0.0;
    delta_t2_amps = 0.0;
}

/**
 * Deallocate the effective Hamiltonian matrices and eigenvectors
 * @todo wrap the current operations in an cleanup() function
 */
CCManyBody::~CCManyBody() {
    free_block(Heff);
    free_block(Heff_mrpt2);
    wfn_->free_memory_ += 2 * sizeof(double) * wfn_->moinfo()->get_nrefs();
    if (d3_ooo.size()) {
        // Triple denominators were allocated. They are only now de-allocated.
        auto nrefs = wfn_->moinfo()->get_nunique();
        auto ooo_indexing = wfn_->blas()->get_index("[ooo]");
        auto vvv_indexing = wfn_->blas()->get_index("[vvv]");
        for (int h = 0; h < wfn_->nirrep(); h++) {
            wfn_->free_memory_ +=
                4 * sizeof(double) * nrefs * (ooo_indexing->get_pairpi(h) + vvv_indexing->get_pairpi(h));
        }
    }
}

/**
 * Creates a CCSort object
 */
void CCManyBody::generate_integrals() {
    // CCSort reads the one and two electron integrals
    // and creates the Fock matrices
    std::make_shared<CCSort>(wfn_, out_of_core_sort);
    //   wfn_->blas()->show_storage();
    wfn_->blas()->compute_storage_strategy();
    //   wfn_->blas()->show_storage();
}

void CCManyBody::generate_triples_denominators() {
    generate_d3_ijk(d3_ooo, true, true, true);
    generate_d3_ijk(d3_ooO, true, true, false);
    generate_d3_ijk(d3_oOO, true, false, false);
    generate_d3_ijk(d3_OOO, false, false, false);
    generate_d3_abc(d3_vvv, true, true, true);
    generate_d3_abc(d3_vvV, true, true, false);
    generate_d3_abc(d3_vVV, true, false, false);
    generate_d3_abc(d3_VVV, false, false, false);
}

void CCManyBody::generate_d3_ijk(std::vector<std::vector<std::vector<double>>>& d3, bool alpha_i, bool alpha_j,
                                 bool alpha_k) {
    d3 = std::vector<std::vector<std::vector<double>>>(wfn_->moinfo()->get_nunique(),
                                                       std::vector<std::vector<double>>(wfn_->nirrep()));
    auto ooo_indexing = wfn_->blas()->get_index("[ooo]");
    for (int h = 0; h < wfn_->nirrep(); h++) {
        wfn_->free_memory_ -= sizeof(double) * wfn_->moinfo()->get_nunique() * ooo_indexing->get_pairpi(h);
    }
    // Loop over references
    for (int ref = 0; ref < wfn_->moinfo()->get_nunique(); ref++) {
        int reference = wfn_->moinfo()->get_ref_number(ref, UniqueRefs);

        // N.B. Never introduce Matrices/Vectors with O or V in the name before you compute the Fock matrix elements
        auto aocc = wfn_->moinfo()->get_aocc(reference, AllRefs);
        auto bocc = wfn_->moinfo()->get_bocc(reference, AllRefs);

        // Build the is_ arrays for reference ref
        std::vector<bool> is_aocc(wfn_->moinfo()->get_nocc(), false);
        std::vector<bool> is_bocc(wfn_->moinfo()->get_nocc(), false);
        for (size_t i = 0; i < aocc.size(); i++) is_aocc[aocc[i]] = true;
        for (size_t i = 0; i < bocc.size(); i++) is_bocc[bocc[i]] = true;

        // Read the Fock matrices
        auto f_oo_Matrix = wfn_->blas()->get_MatTmp("fock[oo]", reference, none);
        auto f_OO_Matrix = wfn_->blas()->get_MatTmp("fock[OO]", reference, none);

        CCMatrix* f_ii_Matrix;
        CCMatrix* f_jj_Matrix;
        CCMatrix* f_kk_Matrix;

        if (alpha_i)
            f_ii_Matrix = f_oo_Matrix.get_CCMatrix();
        else
            f_ii_Matrix = f_OO_Matrix.get_CCMatrix();

        if (alpha_j)
            f_jj_Matrix = f_oo_Matrix.get_CCMatrix();
        else
            f_jj_Matrix = f_OO_Matrix.get_CCMatrix();

        if (alpha_k)
            f_kk_Matrix = f_oo_Matrix.get_CCMatrix();
        else
            f_kk_Matrix = f_OO_Matrix.get_CCMatrix();

        auto& ooo_tuples = ooo_indexing->get_tuples();

        for (int h = 0; h < wfn_->nirrep(); h++) {
            size_t ooo_offset = ooo_indexing->get_first(h);
            d3[ref][h] = std::vector<double>(ooo_indexing->get_pairpi(h), 0);
            for (size_t ijk = 0; ijk < ooo_indexing->get_pairpi(h); ijk++) {
                short i = ooo_tuples[ooo_offset + ijk][0];
                short j = ooo_tuples[ooo_offset + ijk][1];
                short k = ooo_tuples[ooo_offset + ijk][2];

                bool external = true;
                if ((alpha_i && !is_aocc[i]) || (!alpha_i && !is_bocc[i])) external = false;
                if ((alpha_j && !is_aocc[j]) || (!alpha_j && !is_bocc[j])) external = false;
                if ((alpha_k && !is_aocc[k]) || (!alpha_k && !is_bocc[k])) external = false;

                if (external)
                    d3[ref][h][ijk] = f_ii_Matrix->get_two_address_element(i, i) +
                                      f_jj_Matrix->get_two_address_element(j, j) +
                                      f_kk_Matrix->get_two_address_element(k, k);
                else
                    d3[ref][h][ijk] = huge;
            }  // End loop over h,ijk
        }
    }
}

void CCManyBody::generate_d3_abc(std::vector<std::vector<std::vector<double>>>& d3, bool alpha_a, bool alpha_b,
                                 bool alpha_c) {
    d3 = std::vector<std::vector<std::vector<double>>>(wfn_->moinfo()->get_nunique(),
                                                       std::vector<std::vector<double>>(wfn_->nirrep()));
    auto vvv_indexing = wfn_->blas()->get_index("[vvv]");
    for (int h = 0; h < wfn_->nirrep(); h++) {
        wfn_->free_memory_ -= sizeof(double) * wfn_->moinfo()->get_nunique() * vvv_indexing->get_pairpi(h);
    }
    // Loop over references
    for (int ref = 0; ref < wfn_->moinfo()->get_nunique(); ref++) {
        int reference = wfn_->moinfo()->get_ref_number(ref, UniqueRefs);

        // N.B. Never introduce Matrices/Vectors with O or V in the name before you compute the Fock matrix elements
        auto avir = wfn_->moinfo()->get_avir(reference, AllRefs);
        auto bvir = wfn_->moinfo()->get_bvir(reference, AllRefs);

        // Build the is_ arrays for reference ref
        std::vector<bool> is_avir(wfn_->moinfo()->get_nvir(), false);
        std::vector<bool> is_bvir(wfn_->moinfo()->get_nvir(), false);
        for (size_t i = 0; i < avir.size(); i++) is_avir[avir[i]] = true;
        for (size_t i = 0; i < bvir.size(); i++) is_bvir[bvir[i]] = true;

        // Read the Fock matrices
        auto f_vv_Matrix = wfn_->blas()->get_MatTmp("fock[vv]", reference, none);
        auto f_VV_Matrix = wfn_->blas()->get_MatTmp("fock[VV]", reference, none);

        CCMatrix* f_aa_Matrix;
        CCMatrix* f_bb_Matrix;
        CCMatrix* f_cc_Matrix;

        if (alpha_a)
            f_aa_Matrix = f_vv_Matrix.get_CCMatrix();
        else
            f_aa_Matrix = f_VV_Matrix.get_CCMatrix();

        if (alpha_b)
            f_bb_Matrix = f_vv_Matrix.get_CCMatrix();
        else
            f_bb_Matrix = f_VV_Matrix.get_CCMatrix();

        if (alpha_c)
            f_cc_Matrix = f_vv_Matrix.get_CCMatrix();
        else
            f_cc_Matrix = f_VV_Matrix.get_CCMatrix();

        auto& vvv_tuples = vvv_indexing->get_tuples();

        for (int h = 0; h < wfn_->moinfo()->get_nirreps(); h++) {
            size_t vvv_offset = vvv_indexing->get_first(h);
            d3[ref][h] = std::vector<double>(vvv_indexing->get_pairpi(h), 0);
            for (size_t abc = 0; abc < vvv_indexing->get_pairpi(h); abc++) {
                short a = vvv_tuples[vvv_offset + abc][0];
                short b = vvv_tuples[vvv_offset + abc][1];
                short c = vvv_tuples[vvv_offset + abc][2];

                bool external = true;
                if ((alpha_a && !is_avir[a]) || (!alpha_a && !is_bvir[a])) external = false;
                if ((alpha_b && !is_avir[b]) || (!alpha_b && !is_bvir[b])) external = false;
                if ((alpha_c && !is_avir[c]) || (!alpha_c && !is_bvir[c])) external = false;

                if (external)
                    d3[ref][h][abc] = f_aa_Matrix->get_two_address_element(a, a) +
                                      f_bb_Matrix->get_two_address_element(b, b) +
                                      f_cc_Matrix->get_two_address_element(c, c);
                else
                    d3[ref][h][abc] = -huge;
            }  // End lvvp over h,ijk
        }
    }
}

/**
 * Computes the energy for each unique reference determinant
 */
void CCManyBody::compute_reference_energy() {
    Timer timer;

    // Compute the zeroth-order energy for the unique references
    for (int n = 0; n < wfn_->moinfo()->get_nunique(); n++) {
        int unique_n = wfn_->moinfo()->get_ref_number(n, UniqueRefs);
        double ref_energy = wfn_->moinfo()->get_nuc_E() + wfn_->moinfo()->get_fzcore_energy();
        // Grab reference n and the list of occupied orbitals
        auto aocc = wfn_->moinfo()->get_aocc(n, UniqueRefs);
        auto bocc = wfn_->moinfo()->get_bocc(n, UniqueRefs);

        // Read these matrices
        auto f_oo_Matrix = wfn_->blas()->get_MatTmp("fock[o][o]", unique_n, none);
        auto f_OO_Matrix = wfn_->blas()->get_MatTmp("fock[O][O]", unique_n, none);
        auto V_oooo_Matrix = wfn_->blas()->get_MatTmp("<[oo]:[oo]>", none);
        auto V_oOoO_Matrix = wfn_->blas()->get_MatTmp("<[oo]|[oo]>", none);

        for (size_t i = 0; i < aocc.size(); i++) ref_energy += f_oo_Matrix->get_two_address_element(aocc[i], aocc[i]);
        for (size_t i = 0; i < bocc.size(); i++) ref_energy += f_OO_Matrix->get_two_address_element(bocc[i], bocc[i]);

        for (size_t i = 0; i < aocc.size(); i++)
            for (size_t j = 0; j < aocc.size(); j++)
                ref_energy -= 0.5 * V_oooo_Matrix->get_four_address_element(aocc[i], aocc[j], aocc[i], aocc[j]);
        for (size_t i = 0; i < bocc.size(); i++)
            for (size_t j = 0; j < bocc.size(); j++)
                ref_energy -= 0.5 * V_oooo_Matrix->get_four_address_element(bocc[i], bocc[j], bocc[i], bocc[j]);
        for (size_t i = 0; i < aocc.size(); i++)
            for (size_t j = 0; j < bocc.size(); j++)
                ref_energy -= V_oOoO_Matrix->get_four_address_element(aocc[i], bocc[j], aocc[i], bocc[j]);
        // Write the energy to the ERef
        auto ERef_Matrix = wfn_->blas()->get_MatTmp("ERef", unique_n, none);
        ERef_Matrix->set_scalar(ref_energy);
    }
}

void CCManyBody::print_method(const char* text) {
    outfile->Printf("\n");
    outfile->Printf("\n  ==============================================================================");
    outfile->Printf("\n  %s", text);
    outfile->Printf("\n  ==============================================================================");
    outfile->Printf("\n");
}

void CCManyBody::print_eigensystem(int ndets, double** Heff, std::vector<double>& eigenvector) {
    if (ndets < 8) {
        outfile->Printf("\n\n  Heff Matrix\n");
        for (int i = 0; i < ndets; i++) {
            outfile->Printf("\n  ");
            for (int j = 0; j < ndets; j++) outfile->Printf(" %22.15f", Heff[i][j]);
        }
    }

    std::vector<std::pair<double, int>> eigenvector_index_pair;
    for (int i = 0; i < ndets; ++i) {
        eigenvector_index_pair.push_back(std::make_pair(eigenvector[i] * eigenvector[i], i));
    }
    sort(eigenvector_index_pair.begin(), eigenvector_index_pair.end(), std::greater<std::pair<double, int>>());
    int max_size_list = std::min(10, static_cast<int>(eigenvector_index_pair.size()));
    outfile->Printf("\n\n  Most important determinants in the wave function");
    outfile->Printf("\n\n  determinant  eigenvector   eigenvector^2\n");
    for (int i = 0; i < max_size_list; ++i) {
        outfile->Printf("\n  %11d   %9.6f    %9.6f  %s", eigenvector_index_pair[i].second,
                        eigenvector[eigenvector_index_pair[i].second], eigenvector_index_pair[i].first,
                        wfn_->moinfo()->get_determinant_label(eigenvector_index_pair[i].second).c_str());
    }
}

/**
 * This function computes \f$ E = \mathbf{c}^{\dagger} \mathbf{H} \mathbf{c} \f$
 * @param ndets size of the \f$ \mathbf{c} \f$ vector
 * @param H the \f$ \mathbf{H} \f$ matrix stored as a double**
 * @param c the \f$ \mathbf{c} \f$ vector stored as a double*
 * @return \f$ E \f$
 */
double CCManyBody::c_H_c(int ndets, double** H, std::vector<double>& c) {
    double energy = 0.0;
    for (int i = 0; i < ndets; i++)
        for (int j = 0; j < ndets; j++) energy += c[i] * H[i][j] * c[j];
    return (energy);
}

/**
 * This function computes the left and right eigenvalues of a generic real matrix
 * @param root selects the root for which the left-eigenvector must be saved
 * @param ndets size of the matrix
 * @param Heff the \f$ \mathbf{H}^{\mathrm{eff}} \f$ matrix stored as a double**
 * @param eigenvector the \f$ \mathbf{c} \f$ left-eigenvector stored as a double*  * @param initial a bool used to
 * enable root following. initial = true allows you to select a root while initial = false follows the root that has the
 * largest overlap with the previous eigenvector
 * @return
 */
double CCManyBody::diagonalize_Heff(int root, int ndets, double** Heff, std::vector<double>& right_eigenvector,
                                    std::vector<double>& left_eigenvector, bool initial) {
    double energy;

    int lwork = 6 * ndets * ndets;
    std::vector<double> work(lwork, 0);
    std::vector<double> real(ndets, 0);
    std::vector<double> imaginary(ndets, 0);

    auto H = block_matrix(ndets, ndets);
    auto left = block_matrix(ndets, ndets);
    auto right = block_matrix(ndets, ndets);

    for (int i = 0; i < ndets; i++)
        for (int j = 0; j < ndets; j++) H[j][i] = Heff[i][j];

    int info;

    F_DGEEV("V", "V", &ndets, &(H[0][0]), &ndets, real.data(), imaginary.data(), &(left[0][0]), &ndets, &(right[0][0]),
            &ndets, work.data(), &lwork, &info);

    sort_eigensystem(ndets, real, imaginary, left, right);

    if (initial) {
        if (ndets < 8) {
            outfile->Printf("\n\n  Heff Matrix\n");
            for (int i = 0; i < ndets; i++) {
                outfile->Printf("\n  ");
                for (int j = 0; j < ndets; j++) outfile->Printf(" %22.12f", Heff[i][j]);
            }

            outfile->Printf("\n\n  Left Matrix\n");
            for (int i = 0; i < ndets; i++) {
                outfile->Printf("\n  ");
                for (int j = 0; j < ndets; j++) outfile->Printf(" %22.12f", left[j][i]);
            }

            outfile->Printf("\n\n  Right Matrix\n");
            for (int i = 0; i < ndets; i++) {
                outfile->Printf("\n  ");
                for (int j = 0; j < ndets; j++) outfile->Printf(" %22.12f", right[j][i]);
            }

            outfile->Printf("\n\n  Real                  Imaginary\n");
            for (int i = 0; i < ndets; i++) outfile->Printf("\n  %22.12f   %22.12f", real[i], imaginary[i]);
            outfile->Printf("\n");
        } else {
            outfile->Printf("\n\n  There are too many determinants to print the eigensystem");
        }
        outfile->Printf("\n\n  The eigenvalue for root %d is %.12f (%.12f)", root, real[root], imaginary[root]);
    }

    // Select the eigenvector to follow
    if (initial) {
        for (int k = 0; k < ndets; k++) {
            zeroth_order_eigenvector[k] = right[root][k];
            right_eigenvector[k] = right[root][k];
            left_eigenvector[k] = left[root][k];
        }
        energy = real[root];
        // Eliminate the triplet solution if required
        if ((options_.get_bool("LOCK_SINGLET") == 1) && (ndets == 4)) {
            if ((std::fabs(right_eigenvector[0]) < 5.0e-2) && (std::fabs(right_eigenvector[3]) < 5.0e-2) &&
                ((right_eigenvector[1] / right_eigenvector[2]) < -0.5)) {
                outfile->Printf("\n\tSelecting root %d since original root is a triplet\n", root + 1);
                root++;
                for (int k = 0; k < ndets; k++) {
                    right_eigenvector[k] = right[root][k];
                    left_eigenvector[k] = left[root][k];
                }
                energy = real[root];
            }
        }
    } else  // find vector with maximum overlap
    {
        int select_vect = 0;
        double max_overlap = 0.0;
        double overlap = 0.0;
        for (int i = 0; i < ndets; i++) {
            overlap = 0.0;
            for (int m = 0; m < ndets; m++) overlap += zeroth_order_eigenvector[m] * right[i][m];
            overlap = sqrt(overlap * overlap);
            if (overlap > max_overlap) {
                select_vect = i;
                max_overlap = overlap;
            }
        }
        for (int m = 0; m < ndets; m++) {
            right_eigenvector[m] = right[select_vect][m];
            left_eigenvector[m] = left[select_vect][m];
        }
        energy = real[select_vect];
    }

    // Normalize the left-eigenvector to <L|R> = 1
    double lnorm = 0.0;
    for (int m = 0; m < ndets; m++) {
        lnorm += right_eigenvector[m] * left_eigenvector[m];
    }

    for (int m = 0; m < ndets; m++) {
        left_eigenvector[m] = left_eigenvector[m] / lnorm;
    }

    free_block(H);
    free_block(left);
    free_block(right);
    return (energy);
}

void CCManyBody::sort_eigensystem(int ndets, std::vector<double>& real, std::vector<double>& imaginary, double**& left,
                                  double**& right) {
    std::vector<std::pair<double, int>> pairs;
    for (int i = 0; i < ndets; i++) pairs.push_back(std::make_pair(real[i], i));
    sort(pairs.begin(), pairs.end());

    std::vector<double> tempv(ndets, 0);
    std::vector<std::vector<double>> tempm(ndets, std::vector<double>(ndets, 0));

    for (int i = 0; i < ndets; i++) tempv[i] = real[pairs[i].second];
    for (int i = 0; i < ndets; i++) real[i] = tempv[i];

    for (int i = 0; i < ndets; i++) tempv[i] = imaginary[pairs[i].second];
    for (int i = 0; i < ndets; i++) imaginary[i] = tempv[i];

    for (int i = 0; i < ndets; i++)
        for (int j = 0; j < ndets; j++) tempm[i][j] = left[pairs[i].second][j];
    for (int i = 0; i < ndets; i++)
        for (int j = 0; j < ndets; j++) left[i][j] = tempm[i][j];

    for (int i = 0; i < ndets; i++)
        for (int j = 0; j < ndets; j++) tempm[i][j] = right[pairs[i].second][j];
    for (int i = 0; i < ndets; i++)
        for (int j = 0; j < ndets; j++) right[i][j] = tempm[i][j];
}

// void CCManyBody::zero_internal_amps()
//{
//  if(options_get_bool("ZERO_INTERNAL_AMPS")){
//    // Zero internal amplitudes for unique reference i
//    for(int i=0;i<wfn_->moinfo()->get_nunique();i++){
//      int unique_i = wfn_->moinfo()->get_ref_number(i,UniqueRefs);
//      // Loop over reference j
//      for(int j=0;j<wfn_->moinfo()->get_ref_size(AllRefs);j++){
//        vector<pair<int,int> >  alpha_internal_excitation = wfn_->moinfo()->get_alpha_internal_excitation(unique_i,j);
//        vector<pair<int,int> >   beta_internal_excitation = wfn_->moinfo()->get_beta_internal_excitation(unique_i,j);
//
//        // Zero alpha-alpha single excitations
//        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==0)){
//          wfn_->blas()->get_MatTmp("t1[o][v]",unique_i,none)->set_two_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[0].second,
//                                            0.0);
//        }
//
//        // Zero beta-beta single excitations
//        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==1))
//          wfn_->blas()->get_MatTmp("t1[O][V]",unique_i,none)->set_two_address_element(
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//
//        // Zero (alpha,alpha)->(alpha,alpha) double excitations (all permutations)
//        if((alpha_internal_excitation.size()==2)&&(beta_internal_excitation.size()==0)){
//          wfn_->blas()->get_MatTmp("t2[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[0].second,
//                                            alpha_internal_excitation[1].second,
//                                            0.0);
//          wfn_->blas()->get_MatTmp("t2[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[1].second,
//                                            alpha_internal_excitation[0].second,
//                                            0.0);
//          wfn_->blas()->get_MatTmp("t2[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[0].second,
//                                            alpha_internal_excitation[1].second,
//                                            0.0);
//          wfn_->blas()->get_MatTmp("t2[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[1].second,
//                                            alpha_internal_excitation[0].second,
//                                            0.0);
//        }
//
//        // Zero (alpha,beta)->(alpha,beta) double excitations
//        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==1)){
//          wfn_->blas()->get_MatTmp("t2[oO][vV]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            beta_internal_excitation[0].first,
//                                            alpha_internal_excitation[0].second,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//        }
//
//        // Zero (beta,beta)->(beta,beta) double excitations (all permutations)
//        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==2)){
//          wfn_->blas()->get_MatTmp("t2[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[0].second,
//                                            beta_internal_excitation[1].second,
//                                            0.0);
//          wfn_->blas()->get_MatTmp("t2[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[1].second,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//          wfn_->blas()->get_MatTmp("t2[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[0].second,
//                                            beta_internal_excitation[1].second,
//                                            0.0);
//          wfn_->blas()->get_MatTmp("t2[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[1].second,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//        }
//      }
//    }
//  }else{
//    outfile->Printf("\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC.
//    Size-extensivity might be lost\n");
//  }
//}
//
//
// void CCManyBody::zero_t1_internal_amps()
//{
//  if(options_get_bool("ZERO_INTERNAL_AMPS")){
//    // Zero internal amplitudes for unique reference i
//    for(int i=0;i<wfn_->moinfo()->get_nunique();i++){
//      int unique_i = wfn_->moinfo()->get_ref_number(i,UniqueRefs);
//      // Loop over reference j
//      for(int j=0;j<wfn_->moinfo()->get_ref_size(AllRefs);j++){
//        vector<pair<int,int> >  alpha_internal_excitation = wfn_->moinfo()->get_alpha_internal_excitation(unique_i,j);
//        vector<pair<int,int> >   beta_internal_excitation = wfn_->moinfo()->get_beta_internal_excitation(unique_i,j);
//
//        // Zero alpha-alpha single excitations
//        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==0))
//          wfn_->blas()->get_MatTmp("t1[o][v]",unique_i,none)->set_two_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[0].second,
//                                            0.0);
//
//        // Zero beta-beta single excitations
//        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==1))
//          wfn_->blas()->get_MatTmp("t1[O][V]",unique_i,none)->set_two_address_element(
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//      }
//    }
//  }else{
//    outfile->Printf("\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC.
//    Size-extensivity might be lost\n");
//  }
//}
//
// void CCManyBody::zero_internal_delta_amps()
//{
//  if(options_get_bool("ZERO_INTERNAL_AMPS")){
//    // Zero internal amplitudes for unique reference i
//    for(int i=0;i<wfn_->moinfo()->get_nunique();i++){
//      int unique_i = wfn_->moinfo()->get_ref_number(i,UniqueRefs);
//      // Loop over reference j
//      for(int j=0;j<wfn_->moinfo()->get_ref_size(AllRefs);j++){
//        vector<pair<int,int> >  alpha_internal_excitation = wfn_->moinfo()->get_alpha_internal_excitation(unique_i,j);
//        vector<pair<int,int> >   beta_internal_excitation = wfn_->moinfo()->get_beta_internal_excitation(unique_i,j);
//
//        // Zero alpha-alpha single excitations
//        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==0))
//          wfn_->blas()->get_MatTmp("t1_delta[o][v]",unique_i,none)->set_two_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[0].second,
//                                            0.0);
//
//        // Zero beta-beta single excitations
//        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==1))
//          wfn_->blas()->get_MatTmp("t1_delta[O][V]",unique_i,none)->set_two_address_element(
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//
//        // Zero (alpha,alpha)->(alpha,alpha) double excitations (all permutations)
//        if((alpha_internal_excitation.size()==2)&&(beta_internal_excitation.size()==0)){
//          wfn_->blas()->get_MatTmp("t2_delta[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[0].second,
//                                            alpha_internal_excitation[1].second,
//                                            0.0);
//          wfn_->blas()->get_MatTmp("t2_delta[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[1].second,
//                                            alpha_internal_excitation[0].second,
//                                            0.0);
//          wfn_->blas()->get_MatTmp("t2_delta[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[0].second,
//                                            alpha_internal_excitation[1].second,
//                                            0.0);
//          wfn_->blas()->get_MatTmp("t2_delta[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[1].second,
//                                            alpha_internal_excitation[0].second,
//                                            0.0);
//        }
//
//        // Zero (alpha,beta)->(alpha,beta) double excitations
//        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==1)){
//          wfn_->blas()->get_MatTmp("t2_delta[oO][vV]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            beta_internal_excitation[0].first,
//                                            alpha_internal_excitation[0].second,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//        }
//
//        // Zero (beta,beta)->(beta,beta) double excitations (all permutations)
//        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==2)){
//          wfn_->blas()->get_MatTmp("t2_delta[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[0].second,
//                                            beta_internal_excitation[1].second,
//                                            0.0);
//          wfn_->blas()->get_MatTmp("t2_delta[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[1].second,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//          wfn_->blas()->get_MatTmp("t2_delta[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[0].second,
//                                            beta_internal_excitation[1].second,
//                                            0.0);
//          wfn_->blas()->get_MatTmp("t2_delta[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[1].second,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//        }
//      }
//    }
//  }else{
//    outfile->Printf("\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC.
//    Size-extensivity might be lost\n");
//  }
//}

}  // namespace psimrcc
}  // namespace psi
