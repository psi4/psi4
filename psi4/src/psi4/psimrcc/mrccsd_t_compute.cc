/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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
 *  @file ccmrcc_pert_triples.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Computes the (T) correction
 */

#include <cstdlib>

#include "psi4/liboptions/liboptions.h"
#include "psi4/libmoinfo/libmoinfo.h"

#include "blas.h"
#include "debugging.h"
#include "index_iterator.h"
#include "mrcc.h"
#include "mrccsd_t.h"
#include "special_matrices.h"

extern FILE* outfile;

namespace psi {
namespace psimrcc {
extern MOInfo* moinfo;

void MRCCSD_T::compute() {
    outfile->Printf("\n\n  Computing (T) correction using the unrestricted loop algorithm.\n");

    compute_ooo_triples();
    compute_ooO_triples();
    compute_oOO_triples();
    compute_OOO_triples();

    outfile->Printf("\n\n  Mk-MRCCSD(T) diagonal contributions to the effective Hamiltonian:\n");
    outfile->Printf("\n   Ref         E[4]              E_T[4]            E_ST[4]           E_DT[4]");
    outfile->Printf("\n  ------------------------------------------------------------------------------");
    if (nrefs < 100) {
        for (int mu = 0; mu < nrefs; ++mu) {
            outfile->Printf("\n   %2d  ", mu);
            outfile->Printf(" %17.12lf", E4_ooo[mu] + E4_ooO[mu] + E4_oOO[mu] + E4_OOO[mu]);
            outfile->Printf(" %17.12lf", E4T_ooo[mu] + E4T_ooO[mu] + E4T_oOO[mu] + E4T_OOO[mu]);
            outfile->Printf(" %17.12lf", E4ST_ooo[mu] + E4ST_ooO[mu] + E4ST_oOO[mu] + E4ST_OOO[mu]);
            outfile->Printf(" %17.12lf", E4DT_ooo[mu] + E4DT_ooO[mu] + E4DT_oOO[mu] + E4DT_OOO[mu]);
        }
    }
    outfile->Printf("\n   Tot ");
    double E4 = 0.0;
    for (int mu = 0; mu < nrefs; ++mu)
        E4 += (E4_ooo[mu] + E4_ooO[mu] + E4_oOO[mu] + E4_OOO[mu]) * h_eff->get_left_eigenvector(mu) *
              h_eff->get_right_eigenvector(mu);
    outfile->Printf(" %17.12lf", E4);
    double E4T = 0.0;
    for (int mu = 0; mu < nrefs; ++mu)
        E4T += (E4T_ooo[mu] + E4T_ooO[mu] + E4T_oOO[mu] + E4T_OOO[mu]) * h_eff->get_left_eigenvector(mu) *
               h_eff->get_right_eigenvector(mu);
    outfile->Printf(" %17.12lf", E4T);
    double E4ST = 0.0;
    for (int mu = 0; mu < nrefs; ++mu)
        E4ST += (E4ST_ooo[mu] + E4ST_ooO[mu] + E4ST_oOO[mu] + E4ST_OOO[mu]) * h_eff->get_left_eigenvector(mu) *
                h_eff->get_right_eigenvector(mu);
    outfile->Printf(" %17.12lf", E4ST);
    double E4DT = 0.0;
    for (int mu = 0; mu < nrefs; ++mu)
        E4DT += (E4DT_ooo[mu] + E4DT_ooO[mu] + E4DT_oOO[mu] + E4DT_OOO[mu]) * h_eff->get_left_eigenvector(mu) *
                h_eff->get_right_eigenvector(mu);
    outfile->Printf(" %17.12lf", E4DT);
    outfile->Printf("\n  ------------------------------------------------------------------------------");

    outfile->Printf("\n\n  Mk-MRCCSD(T) off-diagonal contributions to the effective Hamiltonian:\n");
    for (int mu = 0; mu < nrefs; ++mu) {
        outfile->Printf("\n");
        for (int nu = 0; nu < nrefs; ++nu) {
            outfile->Printf(" %17.12lf", d_h_eff[mu][nu]);
        }
    }

    if (!options_.get_bool("DIAGONALIZE_HEFF")) {
        double Heff_E = 0.0;
        for (int mu = 0; mu < nrefs; ++mu) {
            for (int nu = 0; nu < nrefs; ++nu) {
                if (mu != nu) {
                    Heff_E += h_eff->get_left_eigenvector(mu) * h_eff->get_right_eigenvector(nu) * d_h_eff[mu][nu];
                }
            }
        }
        double total = 0.0;
        if (options_.get_bool("DIAGONAL_CCSD_T")) {
            outfile->Printf("\n\n  Total     diagonal (T) correction: %17.12f", E4);
            total += E4;
        }
        if (options_.get_bool("OFFDIAGONAL_CCSD_T")) {
            outfile->Printf("\n  Total off-diagonal (T) correction: %17.12f", Heff_E);
            total += Heff_E;
        }
        outfile->Printf("\n  Total              (T) correction: %17.12f", total);
    }

    for (int mu = 0; mu < nrefs; ++mu) {
        for (int nu = 0; nu < nrefs; ++nu) {
            if (mu != nu) {
                if (options_.get_bool("OFFDIAGONAL_CCSD_T")) {  // Option to add the diagonal correction
                    h_eff->add_matrix(mu, nu, d_h_eff[mu][nu]);
                }
            } else {
                if (options_.get_bool("DIAGONAL_CCSD_T")) {  // Option to add the off-diagonal correction
                    h_eff->add_matrix(mu, nu, E4_ooo[mu] + E4_ooO[mu] + E4_oOO[mu] + E4_OOO[mu]);
                }
            }
        }
    }
    h_eff->print_matrix();
}

void MRCCSD_T::compute_ooo_triples() {
    CCIndexIterator ijk("[ooo]");

    size_t tot_cycles = 0;
    size_t tot_triplets = 0;

    for (ijk.first(); !ijk.end(); ijk.next()) {
        int i_sym = o->get_tuple_irrep(ijk.ind_abs<0>());

        size_t i_abs = o->get_tuple_abs_index(ijk.ind_abs<0>());
        size_t j_abs = o->get_tuple_abs_index(ijk.ind_abs<1>());
        size_t k_abs = o->get_tuple_abs_index(ijk.ind_abs<2>());

        size_t i_rel = o->get_tuple_rel_index(ijk.ind_abs<0>());

        size_t ij_abs = oo->get_tuple_abs_index(ijk.ind_abs<0>(), ijk.ind_abs<1>());
        size_t kj_abs = oo->get_tuple_abs_index(ijk.ind_abs<2>(), ijk.ind_abs<1>());
        size_t ik_abs = oo->get_tuple_abs_index(ijk.ind_abs<0>(), ijk.ind_abs<2>());

        int jk_sym = oo->get_tuple_irrep(ijk.ind_abs<1>(), ijk.ind_abs<2>());
        size_t jk_rel = oo->get_tuple_rel_index(ijk.ind_abs<1>(), ijk.ind_abs<2>());

        int ijk_sym = ijk.sym();

        // Compute W for all unique references (d N^7)
        for (int mu = 0; mu < nrefs; ++mu) {
            // Check if ijk belong to the occupied space of mu
            if (is_aocc[mu][i_abs] && is_aocc[mu][j_abs] && is_aocc[mu][k_abs]) {
                Z[mu][ijk_sym]->contract(T2_ij_a_b->get_block_matrix(ij_abs, mu), V_k_bc_e->get_block_matrix(k_abs),
                                         1.0, 0.0);
                Z[mu][ijk_sym]->contract(T2_ij_a_b->get_block_matrix(kj_abs, mu), V_k_bc_e->get_block_matrix(i_abs),
                                         -1.0, 1.0);
                Z[mu][ijk_sym]->contract(T2_ij_a_b->get_block_matrix(ik_abs, mu), V_k_bc_e->get_block_matrix(j_abs),
                                         -1.0, 1.0);
                Z[mu][ijk_sym]->contract(V_jk_c_m->get_block_matrix(ij_abs), T2_i_ab_j->get_block_matrix(k_abs, mu),
                                         -1.0, 1.0);
                Z[mu][ijk_sym]->contract(V_jk_c_m->get_block_matrix(kj_abs), T2_i_ab_j->get_block_matrix(i_abs, mu),
                                         1.0, 1.0);
                Z[mu][ijk_sym]->contract(V_jk_c_m->get_block_matrix(ik_abs), T2_i_ab_j->get_block_matrix(j_abs, mu),
                                         1.0, 1.0);

                W[mu][ijk_sym]->cyclical_permutation_1_2(Z[mu][ijk_sym], vvv, v, vv);
            }
        }

        for (int mu = 0; mu < nrefs; ++mu) {
            T[mu][ijk_sym]->zero();
        }

        // Compute T (d^2 N^6)
        int cycle = 0;
        double oldE = 1.0;
        double newE = 0.0;
        tot_triplets++;
        while (std::fabs(oldE - newE) > threshold) {
            tot_cycles++;
            cycle++;
            oldE = newE;
            newE = 0.0;
            // Iterate the Mk-MRCCSD(T) Equations
            for (int mu = 0; mu < nrefs; ++mu) {
                e4T[mu] = e4ST[mu] = e4DT[mu] = 0.0;
                // Check if ijk belong to the occupied space of mu
                if (is_aocc[mu][i_abs] && is_aocc[mu][j_abs] && is_aocc[mu][k_abs]) {
                    double*** F_ov_mu = F_ov[mu];
                    double*** T1_ov_mu = T1_ov[mu];
                    double*** T2_oovv_mu = T2_oovv[mu];

                    double D_ijk = e_oo[mu][i_abs] + e_oo[mu][j_abs] + e_oo[mu][k_abs];

                    // Add W
                    Z[mu][ijk_sym]->add(W[mu][ijk_sym], 0.0, 1.0);

                    // Add the coupling terms
                    for (int nu = 0; nu < nrefs; ++nu) {
                        if (nu != mu) {
                            Z[mu][ijk_sym]->add(T[nu][ijk_sym], 1.0, Mk_factor[mu][nu]);
                        }
                    }

                    // Divide by the denominator
                    std::vector<double>& e_vv_mu = e_vv[mu];
                    std::vector<bool>& is_avir_mu = is_avir[mu];

                    CCIndexIterator abc("[vvv]", ijk_sym);
                    //          abc.reset();
                    //          abc.set_irrep();
                    // Loop over abc
                    for (abc.first(); !abc.end(); abc.next()) {
                        size_t a_abs = v->get_tuple_abs_index(abc.ind_abs<0>());
                        size_t b_abs = v->get_tuple_abs_index(abc.ind_abs<1>());
                        size_t c_abs = v->get_tuple_abs_index(abc.ind_abs<2>());
                        if (is_avir_mu[a_abs] && is_avir_mu[b_abs] && is_avir_mu[c_abs]) {
                            int a_sym = v->get_tuple_irrep(abc.ind_abs<0>());
                            int bc_sym = vv->get_tuple_irrep(abc.ind_abs<1>(), abc.ind_abs<2>());
                            size_t a_rel = v->get_tuple_rel_index(abc.ind_abs<0>());
                            size_t bc_rel = vv->get_tuple_rel_index(abc.ind_abs<1>(), abc.ind_abs<2>());

                            double D_abc = e_vv_mu[a_abs] + e_vv_mu[b_abs] + e_vv_mu[c_abs];

                            // Update T
                            T[mu][ijk_sym]->set(
                                a_sym, a_rel, bc_rel,
                                Z[mu][ijk_sym]->get(a_sym, a_rel, bc_rel) / (Mk_shift[mu] + D_ijk - D_abc));

                            // Compute the energy
                            e4T[mu] += W[mu][ijk_sym]->get(a_sym, a_rel, bc_rel) *
                                       T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel) / 36.0;
                            if ((i_sym == a_sym) & (jk_sym == bc_sym)) {
                                e4ST[mu] += 0.25 * T1_ov_mu[i_sym][i_rel][a_rel] * V_oovv[jk_sym][jk_rel][bc_rel] *
                                            T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel);
                                e4DT[mu] += 0.25 * F_ov_mu[i_sym][i_rel][a_rel] * T2_oovv_mu[jk_sym][jk_rel][bc_rel] *
                                            T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel);
                            }
                        }
                    }  // End loop over abc
                    newE += std::fabs(e4T[mu]) + std::fabs(e4ST[mu]) + std::fabs(e4DT[mu]);
                }  // End loop over allowed ijk
            }      // End of iterations
        }

        // Compute the contributions to the off-diagonal elements of Heff
        for (int mu = 0; mu < nrefs; ++mu) {
            compute_ooo_contribution_to_Heff(i_abs, j_abs, k_abs, mu, T[mu][ijk_sym]);
        }

        // Add the energy contributions from ijk
        for (int mu = 0; mu < nrefs; ++mu) {
            E4T_ooo[mu] += e4T[mu];
            E4ST_ooo[mu] += e4ST[mu];
            E4DT_ooo[mu] += e4DT[mu];
        }
    }  // End loop over ijk

    for (int mu = 0; mu < nrefs; ++mu) {
        outfile->Printf("\n  E_T[4]  (aaa) = %20.15lf (%d)", E4T_ooo[mu], mu);
        outfile->Printf("\n  E_ST[4] (aaa) = %20.15lf (%d)", E4ST_ooo[mu], mu);
        outfile->Printf("\n  E_DT[4] (aaa) = %20.15lf (%d)", E4DT_ooo[mu], mu);
        E4_ooo[mu] = E4T_ooo[mu] + E4ST_ooo[mu] + E4DT_ooo[mu];
    }
}

void MRCCSD_T::compute_OOO_triples() {
    CCIndexIterator ijk("[ooo]");
    CCIndexIterator abc("[vvv]");

    for (ijk.first(); !ijk.end(); ijk.next()) {
        int i_sym = o->get_tuple_irrep(ijk.ind_abs<0>());

        size_t i_abs = o->get_tuple_abs_index(ijk.ind_abs<0>());
        size_t j_abs = o->get_tuple_abs_index(ijk.ind_abs<1>());
        size_t k_abs = o->get_tuple_abs_index(ijk.ind_abs<2>());

        size_t i_rel = o->get_tuple_rel_index(ijk.ind_abs<0>());

        size_t ij_abs = oo->get_tuple_abs_index(ijk.ind_abs<0>(), ijk.ind_abs<1>());
        size_t kj_abs = oo->get_tuple_abs_index(ijk.ind_abs<2>(), ijk.ind_abs<1>());
        size_t ik_abs = oo->get_tuple_abs_index(ijk.ind_abs<0>(), ijk.ind_abs<2>());

        int jk_sym = oo->get_tuple_irrep(ijk.ind_abs<1>(), ijk.ind_abs<2>());
        size_t jk_rel = oo->get_tuple_rel_index(ijk.ind_abs<1>(), ijk.ind_abs<2>());

        int ijk_sym = ijk.sym();

        // Compute W for all unique references (d N^7)
        for (int mu = 0; mu < nrefs; ++mu) {
            // Check if ijk belong to the occupied space of mu
            if (is_bocc[mu][i_abs] && is_bocc[mu][j_abs] && is_bocc[mu][k_abs]) {
                Z[mu][ijk_sym]->contract(T2_IJ_A_B->get_block_matrix(ij_abs, mu), V_k_bc_e->get_block_matrix(k_abs),
                                         1.0, 0.0);
                Z[mu][ijk_sym]->contract(T2_IJ_A_B->get_block_matrix(kj_abs, mu), V_k_bc_e->get_block_matrix(i_abs),
                                         -1.0, 1.0);
                Z[mu][ijk_sym]->contract(T2_IJ_A_B->get_block_matrix(ik_abs, mu), V_k_bc_e->get_block_matrix(j_abs),
                                         -1.0, 1.0);

                Z[mu][ijk_sym]->contract(V_jk_c_m->get_block_matrix(ij_abs), T2_I_AB_J->get_block_matrix(k_abs, mu),
                                         -1.0, 1.0);
                Z[mu][ijk_sym]->contract(V_jk_c_m->get_block_matrix(kj_abs), T2_I_AB_J->get_block_matrix(i_abs, mu),
                                         1.0, 1.0);
                Z[mu][ijk_sym]->contract(V_jk_c_m->get_block_matrix(ik_abs), T2_I_AB_J->get_block_matrix(j_abs, mu),
                                         1.0, 1.0);

                W[mu][ijk_sym]->cyclical_permutation_1_2(Z[mu][ijk_sym], vvv, v, vv);
            }
        }

        for (int mu = 0; mu < nrefs; ++mu) {
            T[mu][ijk_sym]->zero();
        }

        // Compute T (d^2 N^6)
        int cycle = 0;
        double oldE = 1.0;
        double newE = 0.0;
        while (std::fabs(oldE - newE) > threshold) {
            cycle++;
            oldE = newE;
            newE = 0.0;
            // Iterate the Mk-MRCCSD(T) Equations
            for (int mu = 0; mu < nrefs; ++mu) {
                e4T[mu] = e4ST[mu] = e4DT[mu] = 0.0;
                // Check if ijk belong to the occupied space of mu
                if (is_bocc[mu][i_abs] && is_bocc[mu][j_abs] && is_bocc[mu][k_abs]) {
                    double*** F_OV_mu = F_OV[mu];
                    double*** T1_OV_mu = T1_OV[mu];
                    double*** T2_OOVV_mu = T2_OOVV[mu];

                    double D_IJK = e_OO[mu][i_abs] + e_OO[mu][j_abs] + e_OO[mu][k_abs];

                    // Add W
                    Z[mu][ijk_sym]->add(W[mu][ijk_sym], 0.0, 1.0);

                    // Add the coupling terms
                    for (int nu = 0; nu < nrefs; ++nu) {
                        if (nu != mu) {
                            Z[mu][ijk_sym]->add(T[nu][ijk_sym], 1.0, Mk_factor[mu][nu]);
                        }
                    }

                    // Divide by the denominator
                    std::vector<double>& e_VV_mu = e_VV[mu];
                    std::vector<bool>& is_bvir_mu = is_bvir[mu];

                    CCIndexIterator abc("[vvv]", ijk_sym);
                    //          abc.reset();
                    //          abc.set_irrep();          // Loop over abc
                    for (abc.first(); !abc.end(); abc.next()) {
                        size_t a_abs = v->get_tuple_abs_index(abc.ind_abs<0>());
                        size_t b_abs = v->get_tuple_abs_index(abc.ind_abs<1>());
                        size_t c_abs = v->get_tuple_abs_index(abc.ind_abs<2>());
                        if (is_bvir_mu[a_abs] && is_bvir_mu[b_abs] && is_bvir_mu[c_abs]) {
                            int a_sym = v->get_tuple_irrep(abc.ind_abs<0>());
                            int bc_sym = vv->get_tuple_irrep(abc.ind_abs<1>(), abc.ind_abs<2>());
                            size_t a_rel = v->get_tuple_rel_index(abc.ind_abs<0>());
                            size_t bc_rel = vv->get_tuple_rel_index(abc.ind_abs<1>(), abc.ind_abs<2>());

                            double D_ABC = e_VV_mu[a_abs] + e_VV_mu[b_abs] + e_VV_mu[c_abs];

                            // Update T
                            T[mu][ijk_sym]->set(
                                a_sym, a_rel, bc_rel,
                                Z[mu][ijk_sym]->get(a_sym, a_rel, bc_rel) / (Mk_shift[mu] + D_IJK - D_ABC));

                            // Compute the energy
                            e4T[mu] += W[mu][ijk_sym]->get(a_sym, a_rel, bc_rel) *
                                       T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel) / 36.0;
                            if ((i_sym == a_sym) & (jk_sym == bc_sym)) {
                                e4ST[mu] += 0.25 * T1_OV_mu[i_sym][i_rel][a_rel] * V_oovv[jk_sym][jk_rel][bc_rel] *
                                            T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel);
                                e4DT[mu] += 0.25 * F_OV_mu[i_sym][i_rel][a_rel] * T2_OOVV_mu[jk_sym][jk_rel][bc_rel] *
                                            T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel);
                            }
                        }
                    }  // End loop over abc
                    newE += std::fabs(e4T[mu]) + std::fabs(e4ST[mu]) + std::fabs(e4DT[mu]);
                }  // End loop over allowed ijk
            }      // End of iterations
        }

        // Compute the contributions to the off-diagonal elements of Heff
        for (int mu = 0; mu < nrefs; ++mu) {
            compute_OOO_contribution_to_Heff(i_abs, j_abs, k_abs, mu, T[mu][ijk_sym]);
        }

        // Add the energy contributions from ijk
        for (int mu = 0; mu < nrefs; ++mu) {
            E4T_OOO[mu] += e4T[mu];
            E4ST_OOO[mu] += e4ST[mu];
            E4DT_OOO[mu] += e4DT[mu];
        }
    }

    for (int mu = 0; mu < nrefs; ++mu) {
        //    outfile->Printf("\n  E_T[4]  (bbb) = %20.15lf (%d)",E4T_OOO[mu],mu);
        //    outfile->Printf("\n  E_ST[4] (bbb) = %20.15lf (%d)",E4ST_OOO[mu],mu);
        //    outfile->Printf("\n  E_DT[4] (bbb) = %20.15lf (%d)",E4DT_OOO[mu],mu);
        E4_OOO[mu] = E4T_OOO[mu] + E4ST_OOO[mu] + E4DT_OOO[mu];
    }
}

void MRCCSD_T::compute_ooO_triples() {
    CCIndexIterator ijk("[ooo]");
    CCIndexIterator abc("[vvv]");

    for (ijk.first(); !ijk.end(); ijk.next()) {
        int i_sym = o->get_tuple_irrep(ijk.ind_abs<0>());
        int k_sym = o->get_tuple_irrep(ijk.ind_abs<2>());

        size_t i_abs = o->get_tuple_abs_index(ijk.ind_abs<0>());
        size_t j_abs = o->get_tuple_abs_index(ijk.ind_abs<1>());
        size_t k_abs = o->get_tuple_abs_index(ijk.ind_abs<2>());

        size_t i_rel = o->get_tuple_rel_index(ijk.ind_abs<0>());
        size_t k_rel = o->get_tuple_rel_index(ijk.ind_abs<2>());

        size_t ij_abs = oo->get_tuple_abs_index(ijk.ind_abs<0>(), ijk.ind_abs<1>());
        size_t ik_abs = oo->get_tuple_abs_index(ijk.ind_abs<0>(), ijk.ind_abs<2>());
        size_t jk_abs = oo->get_tuple_abs_index(ijk.ind_abs<1>(), ijk.ind_abs<2>());

        int ij_sym = oo->get_tuple_irrep(ijk.ind_abs<0>(), ijk.ind_abs<1>());
        size_t ij_rel = oo->get_tuple_rel_index(ijk.ind_abs<0>(), ijk.ind_abs<1>());
        int jk_sym = oo->get_tuple_irrep(ijk.ind_abs<1>(), ijk.ind_abs<2>());
        size_t jk_rel = oo->get_tuple_rel_index(ijk.ind_abs<1>(), ijk.ind_abs<2>());

        int ijk_sym = ijk.sym();

        // Compute W for all unique references (d N^7)
        for (int mu = 0; mu < nrefs; ++mu) {
            // Check if ijk belong to the occupied space of mu
            if (is_aocc[mu][i_abs] && is_aocc[mu][j_abs] && is_bocc[mu][k_abs]) {
                Z[mu][ijk_sym]->contract(T2_ij_a_b->get_block_matrix(ij_abs, mu), V_K_bC_e->get_block_matrix(k_abs),
                                         1.0, 0.0);

                Z[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(jk_abs, mu), V_k_bC_E->get_block_matrix(i_abs),
                                         -1.0, 1.0);
                Z[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(ik_abs, mu), V_k_bC_E->get_block_matrix(j_abs),
                                         1.0, 1.0);

                Z[mu][ijk_sym]->contract(V_jk_c_m->get_block_matrix(ij_abs), T2_J_aB_i->get_block_matrix(k_abs, mu),
                                         1.0, 1.0);

                Z[mu][ijk_sym]->contract(V_jK_c_M->get_block_matrix(jk_abs), T2_i_aB_J->get_block_matrix(i_abs, mu),
                                         1.0, 1.0);
                Z[mu][ijk_sym]->contract(V_jK_c_M->get_block_matrix(ik_abs), T2_i_aB_J->get_block_matrix(j_abs, mu),
                                         -1.0, 1.0);

                W[mu][ijk_sym]->a_b_permutation_1_2(Z[mu][ijk_sym], vvv, v, vv);

                Z[mu][ijk_sym]->contract(T2_iJ_B_a->get_block_matrix(ik_abs, mu), V_k_bc_e->get_block_matrix(j_abs),
                                         1.0, 0.0);
                Z[mu][ijk_sym]->contract(T2_iJ_B_a->get_block_matrix(jk_abs, mu), V_k_bc_e->get_block_matrix(i_abs),
                                         -1.0, 1.0);

                Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(jk_abs), T2_i_ab_j->get_block_matrix(i_abs, mu),
                                         -1.0, 1.0);
                Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(ik_abs), T2_i_ab_j->get_block_matrix(j_abs, mu),
                                         1.0, 1.0);

                W[mu][ijk_sym]->add_c_ab_permutation_1_2(Z[mu][ijk_sym], vvv, v, vv);
            }
        }

        for (int mu = 0; mu < nrefs; ++mu) {
            T[mu][ijk_sym]->zero();
        }

        // Compute T (d^2 N^6)
        int cycle = 0;
        double oldE = 1.0;
        double newE = 0.0;
        while (std::fabs(oldE - newE) > threshold) {
            cycle++;
            oldE = newE;
            newE = 0.0;
            // Iterate the Mk-MRCCSD(T) Equations
            for (int mu = 0; mu < nrefs; ++mu) {
                e4T[mu] = e4ST[mu] = e4DT[mu] = 0.0;
                // Check if ijk belong to the occupied space of mu
                if (is_aocc[mu][i_abs] && is_aocc[mu][j_abs] && is_bocc[mu][k_abs]) {
                    double*** F_ov_mu = F_ov[mu];
                    double*** F_OV_mu = F_OV[mu];
                    double*** T1_ov_mu = T1_ov[mu];
                    double*** T1_OV_mu = T1_OV[mu];
                    double*** T2_oovv_mu = T2_oovv[mu];
                    double*** T2_oOvV_mu = T2_oOvV[mu];

                    double D_ijK = e_oo[mu][i_abs] + e_oo[mu][j_abs] + e_OO[mu][k_abs];

                    // Add W
                    Z[mu][ijk_sym]->add(W[mu][ijk_sym], 0.0, 1.0);

                    // Add the coupling terms
                    for (int nu = 0; nu < nrefs; ++nu) {
                        if (nu != mu) {
                            Z[mu][ijk_sym]->add(T[nu][ijk_sym], 1.0, Mk_factor[mu][nu]);
                        }
                    }

                    // Divide by the denominator
                    std::vector<double>& e_vv_mu = e_vv[mu];
                    std::vector<double>& e_VV_mu = e_VV[mu];
                    std::vector<bool>& is_avir_mu = is_avir[mu];
                    std::vector<bool>& is_bvir_mu = is_bvir[mu];

                    CCIndexIterator abc("[vvv]", ijk_sym);
                    //          abc.reset();
                    //          abc.set_irrep();
                    // Loop over abc
                    for (abc.first(); !abc.end(); abc.next()) {
                        size_t a_abs = v->get_tuple_abs_index(abc.ind_abs<0>());
                        size_t b_abs = v->get_tuple_abs_index(abc.ind_abs<1>());
                        size_t c_abs = v->get_tuple_abs_index(abc.ind_abs<2>());
                        if (is_avir_mu[a_abs] && is_avir_mu[b_abs] && is_bvir_mu[c_abs]) {
                            int a_sym = v->get_tuple_irrep(abc.ind_abs<0>());
                            int c_sym = v->get_tuple_irrep(abc.ind_abs<2>());
                            int ab_sym = vv->get_tuple_irrep(abc.ind_abs<0>(), abc.ind_abs<1>());
                            int bc_sym = vv->get_tuple_irrep(abc.ind_abs<1>(), abc.ind_abs<2>());
                            size_t a_rel = v->get_tuple_rel_index(abc.ind_abs<0>());
                            size_t c_rel = v->get_tuple_rel_index(abc.ind_abs<2>());
                            size_t ab_rel = vv->get_tuple_rel_index(abc.ind_abs<0>(), abc.ind_abs<1>());
                            size_t bc_rel = vv->get_tuple_rel_index(abc.ind_abs<1>(), abc.ind_abs<2>());

                            double D_abC = e_vv_mu[a_abs] + e_vv_mu[b_abs] + e_VV_mu[c_abs];

                            // Update T
                            T[mu][ijk_sym]->set(
                                a_sym, a_rel, bc_rel,
                                Z[mu][ijk_sym]->get(a_sym, a_rel, bc_rel) / (Mk_shift[mu] + D_ijK - D_abC));

                            // Compute the energy
                            e4T[mu] += W[mu][ijk_sym]->get(a_sym, a_rel, bc_rel) *
                                       T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel) / 4.0;
                            if ((i_sym == a_sym) & (jk_sym == bc_sym)) {
                                e4ST[mu] += T1_ov_mu[i_sym][i_rel][a_rel] * V_oOvV[jk_sym][jk_rel][bc_rel] *
                                            T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel);
                                e4DT[mu] += F_ov_mu[i_sym][i_rel][a_rel] * T2_oOvV_mu[jk_sym][jk_rel][bc_rel] *
                                            T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel);
                            }
                            if ((k_sym == c_sym) & (ij_sym == ab_sym)) {
                                e4ST[mu] += 0.25 * T1_OV_mu[k_sym][k_rel][c_rel] * V_oovv[ij_sym][ij_rel][ab_rel] *
                                            T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel);
                                e4DT[mu] += 0.25 * F_OV_mu[k_sym][k_rel][c_rel] * T2_oovv_mu[ij_sym][ij_rel][ab_rel] *
                                            T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel);
                            }
                        }
                    }  // End loop over abc
                    newE += std::fabs(e4T[mu]) + std::fabs(e4ST[mu]) + std::fabs(e4DT[mu]);
                }  // End loop over allowed ijk
            }      // End of iterations
        }

        // Compute the contributions to the off-diagonal elements of Heff
        for (int mu = 0; mu < nrefs; ++mu) {
            compute_ooO_contribution_to_Heff(i_abs, j_abs, k_abs, mu, T[mu][ijk_sym]);
        }

        // Add the energy contributions from ijk
        for (int mu = 0; mu < nrefs; ++mu) {
            E4T_ooO[mu] += e4T[mu];
            E4ST_ooO[mu] += e4ST[mu];
            E4DT_ooO[mu] += e4DT[mu];
        }

    }  // End loop over ijk

    for (int mu = 0; mu < nrefs; ++mu) {
        outfile->Printf("\n  E_T[4]  (aab) = %20.15lf (%d)", E4T_ooO[mu], mu);
        outfile->Printf("\n  E_ST[4] (aab) = %20.15lf (%d)", E4ST_ooO[mu], mu);
        outfile->Printf("\n  E_DT[4] (aab) = %20.15lf (%d)", E4DT_ooO[mu], mu);
        E4_ooO[mu] = E4T_ooO[mu] + E4ST_ooO[mu] + E4DT_ooO[mu];
    }
}

void MRCCSD_T::compute_oOO_triples() {
    CCIndexIterator ijk("[ooo]");
    CCIndexIterator abc("[vvv]");

    for (ijk.first(); !ijk.end(); ijk.next()) {
        int i_sym = o->get_tuple_irrep(ijk.ind_abs<0>());
        int k_sym = o->get_tuple_irrep(ijk.ind_abs<2>());

        size_t i_abs = o->get_tuple_abs_index(ijk.ind_abs<0>());
        size_t j_abs = o->get_tuple_abs_index(ijk.ind_abs<1>());
        size_t k_abs = o->get_tuple_abs_index(ijk.ind_abs<2>());

        size_t i_rel = o->get_tuple_rel_index(ijk.ind_abs<0>());
        size_t k_rel = o->get_tuple_rel_index(ijk.ind_abs<2>());

        size_t ij_abs = oo->get_tuple_abs_index(ijk.ind_abs<0>(), ijk.ind_abs<1>());
        size_t ik_abs = oo->get_tuple_abs_index(ijk.ind_abs<0>(), ijk.ind_abs<2>());
        size_t jk_abs = oo->get_tuple_abs_index(ijk.ind_abs<1>(), ijk.ind_abs<2>());

        int ij_sym = oo->get_tuple_irrep(ijk.ind_abs<0>(), ijk.ind_abs<1>());
        size_t ij_rel = oo->get_tuple_rel_index(ijk.ind_abs<0>(), ijk.ind_abs<1>());
        int jk_sym = oo->get_tuple_irrep(ijk.ind_abs<1>(), ijk.ind_abs<2>());
        size_t jk_rel = oo->get_tuple_rel_index(ijk.ind_abs<1>(), ijk.ind_abs<2>());

        int ijk_sym = ijk.sym();

        // Compute W for all unique references (d N^7)
        for (int mu = 0; mu < nrefs; ++mu) {
            // Check if ijk belong to the occupied space of mu
            if (is_aocc[mu][i_abs] && is_bocc[mu][j_abs] && is_bocc[mu][k_abs]) {
                W[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(ij_abs, mu), V_k_bc_e->get_block_matrix(k_abs),
                                         1.0, 0.0);
                W[mu][ijk_sym]->contract(T2_iJ_a_B->get_block_matrix(ik_abs, mu), V_k_bc_e->get_block_matrix(j_abs),
                                         -1.0, 1.0);

                W[mu][ijk_sym]->contract(V_jK_c_M->get_block_matrix(ij_abs), T2_I_AB_J->get_block_matrix(k_abs, mu),
                                         1.0, 1.0);
                W[mu][ijk_sym]->contract(V_jK_c_M->get_block_matrix(ik_abs), T2_I_AB_J->get_block_matrix(j_abs, mu),
                                         -1.0, 1.0);

                Z[mu][ijk_sym]->contract(T2_IJ_A_B->get_block_matrix(jk_abs, mu), V_k_bC_E->get_block_matrix(i_abs),
                                         1.0, 0.0);

                Z[mu][ijk_sym]->contract(T2_iJ_B_a->get_block_matrix(ij_abs, mu), V_K_bC_e->get_block_matrix(k_abs),
                                         1.0, 1.0);
                Z[mu][ijk_sym]->contract(T2_iJ_B_a->get_block_matrix(ik_abs, mu), V_K_bC_e->get_block_matrix(j_abs),
                                         -1.0, 1.0);

                Z[mu][ijk_sym]->contract(V_jk_c_m->get_block_matrix(jk_abs), T2_i_aB_J->get_block_matrix(i_abs, mu),
                                         1.0, 1.0);

                Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(ij_abs), T2_J_aB_i->get_block_matrix(k_abs, mu),
                                         -1.0, 1.0);
                Z[mu][ijk_sym]->contract(V_jK_C_m->get_block_matrix(ik_abs), T2_J_aB_i->get_block_matrix(j_abs, mu),
                                         1.0, 1.0);

                W[mu][ijk_sym]->add_permutation_1_2(1.0, Z[mu][ijk_sym], vvv, v, vv, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0);
            }
        }

        for (int mu = 0; mu < nrefs; ++mu) {
            T[mu][ijk_sym]->zero();
        }

        // Compute T (d^2 N^6)
        int cycle = 0;
        double oldE = 1.0;
        double newE = 0.0;
        while (std::fabs(oldE - newE) > threshold) {
            cycle++;
            oldE = newE;
            newE = 0.0;
            // Iterate the Mk-MRCCSD(T) Equations
            for (int mu = 0; mu < nrefs; ++mu) {
                e4T[mu] = e4ST[mu] = e4DT[mu] = 0.0;
                // Check if ijk belong to the occupied space of mu
                if (is_aocc[mu][i_abs] && is_bocc[mu][j_abs] && is_bocc[mu][k_abs]) {
                    double*** F_ov_mu = F_ov[mu];
                    double*** F_OV_mu = F_OV[mu];
                    double*** T1_ov_mu = T1_ov[mu];
                    double*** T1_OV_mu = T1_OV[mu];
                    double*** T2_oOvV_mu = T2_oOvV[mu];
                    double*** T2_OOVV_mu = T2_OOVV[mu];

                    double D_iJK = e_oo[mu][i_abs] + e_OO[mu][j_abs] + e_OO[mu][k_abs];

                    // Add W
                    Z[mu][ijk_sym]->add(W[mu][ijk_sym], 0.0, 1.0);

                    // Add the coupling terms
                    for (int nu = 0; nu < nrefs; ++nu) {
                        if (nu != mu) {
                            Z[mu][ijk_sym]->add(T[nu][ijk_sym], 1.0, Mk_factor[mu][nu]);
                        }
                    }

                    // Divide by the denominator
                    std::vector<double>& e_vv_mu = e_vv[mu];
                    std::vector<double>& e_VV_mu = e_VV[mu];
                    std::vector<bool>& is_avir_mu = is_avir[mu];
                    std::vector<bool>& is_bvir_mu = is_bvir[mu];

                    CCIndexIterator abc("[vvv]", ijk_sym);
                    //          abc.reset();
                    //          abc.set_irrep();
                    // Loop over abc
                    for (abc.first(); !abc.end(); abc.next()) {
                        size_t a_abs = v->get_tuple_abs_index(abc.ind_abs<0>());
                        size_t b_abs = v->get_tuple_abs_index(abc.ind_abs<1>());
                        size_t c_abs = v->get_tuple_abs_index(abc.ind_abs<2>());
                        if (is_avir_mu[a_abs] && is_bvir_mu[b_abs] && is_bvir_mu[c_abs]) {
                            int a_sym = v->get_tuple_irrep(abc.ind_abs<0>());
                            int c_sym = v->get_tuple_irrep(abc.ind_abs<2>());
                            int ab_sym = vv->get_tuple_irrep(abc.ind_abs<0>(), abc.ind_abs<1>());
                            int bc_sym = vv->get_tuple_irrep(abc.ind_abs<1>(), abc.ind_abs<2>());

                            size_t a_rel = v->get_tuple_rel_index(abc.ind_abs<0>());
                            size_t c_rel = v->get_tuple_rel_index(abc.ind_abs<2>());
                            size_t ab_rel = vv->get_tuple_rel_index(abc.ind_abs<0>(), abc.ind_abs<1>());
                            size_t bc_rel = vv->get_tuple_rel_index(abc.ind_abs<1>(), abc.ind_abs<2>());

                            double D_aBC = e_vv_mu[a_abs] + e_VV_mu[b_abs] + e_VV_mu[c_abs];

                            // Update T
                            T[mu][ijk_sym]->set(
                                a_sym, a_rel, bc_rel,
                                Z[mu][ijk_sym]->get(a_sym, a_rel, bc_rel) / (Mk_shift[mu] + D_iJK - D_aBC));

                            // Compute the energy
                            e4T[mu] += W[mu][ijk_sym]->get(a_sym, a_rel, bc_rel) *
                                       T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel) / 4.0;
                            if ((i_sym == a_sym) & (jk_sym == bc_sym)) {
                                e4ST[mu] += 0.25 * T1_ov_mu[i_sym][i_rel][a_rel] * V_oovv[jk_sym][jk_rel][bc_rel] *
                                            T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel);
                                e4DT[mu] += 0.25 * F_ov_mu[i_sym][i_rel][a_rel] * T2_OOVV_mu[jk_sym][jk_rel][bc_rel] *
                                            T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel);
                            }
                            if ((k_sym == c_sym) & (ij_sym == ab_sym)) {
                                e4ST[mu] += T1_OV_mu[k_sym][k_rel][c_rel] * V_oOvV[ij_sym][ij_rel][ab_rel] *
                                            T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel);
                                e4DT[mu] += F_OV_mu[k_sym][k_rel][c_rel] * T2_oOvV_mu[ij_sym][ij_rel][ab_rel] *
                                            T[mu][ijk_sym]->get(a_sym, a_rel, bc_rel);
                            }
                        }
                    }  // End loop over abc
                    newE += std::fabs(e4T[mu]) + std::fabs(e4ST[mu]) + std::fabs(e4DT[mu]);
                }  // End loop over allowed ijk
            }      // End of iterations
        }

        // Compute the contributions to the off-diagonal elements of Heff
        for (int mu = 0; mu < nrefs; ++mu) {
            compute_oOO_contribution_to_Heff(i_abs, j_abs, k_abs, mu, T[mu][ijk_sym]);
        }

        // Add the energy contributions from ijk
        for (int mu = 0; mu < nrefs; ++mu) {
            E4T_oOO[mu] += e4T[mu];
            E4ST_oOO[mu] += e4ST[mu];
            E4DT_oOO[mu] += e4DT[mu];
        }
    }

    for (int mu = 0; mu < nrefs; ++mu) {
        outfile->Printf("\n  E_T[4]  (abb) = %20.15lf (%d)", E4T_oOO[mu], mu);
        outfile->Printf("\n  E_ST[4] (abb) = %20.15lf (%d)", E4ST_oOO[mu], mu);
        outfile->Printf("\n  E_DT[4] (abb) = %20.15lf (%d)", E4DT_oOO[mu], mu);
        E4_oOO[mu] = E4T_oOO[mu] + E4ST_oOO[mu] + E4DT_oOO[mu];
    }
}

}  // namespace psimrcc
}  // namespace psi
