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

#include <cstdlib>

#include "psi4/psi4-dec.h"
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/process.h"

#include "idmrpt2.h"
#include "matrix.h"
#include "blas.h"
#include "sort.h"
#include "updater.h"
#include "psimrcc_wfn.h"

namespace psi {
namespace psimrcc {

IDMRPT2::IDMRPT2(std::shared_ptr<PSIMRCCWfn> wfn, Options& options) : CCManyBody(wfn, options) {
    triples_type = pt2;
    updater_ = std::make_shared<MkUpdater>(wfn, options);
    add_matrices();
}

IDMRPT2::~IDMRPT2() {}

double IDMRPT2::compute_energy() {
    read_mrpt2_integrals();
    generate_denominators();
    compute_reference_energy();

    // Build H in the model space. (NB the amplitudes are zero)
    build_F_intermediates();
    build_Heff_mrpt2_diagonal();
    build_Heff_mrpt2_offdiagonal();

    for (int m = 0; m < wfn_->moinfo()->get_nrefs(); m++)
        for (int n = 0; n < wfn_->moinfo()->get_nrefs(); n++) Heff[m][n] = Heff_mrpt2[m][n];
    cas_energy = diagonalize_Heff(wfn_->moinfo()->get_root(), wfn_->moinfo()->get_nrefs(), Heff_mrpt2,
                                  zeroth_order_eigenvector, left_eigenvector, true);
    current_energy = cas_energy;
    old_energy = current_energy;

    wfn_->blas()->solve("d'1[o][v]{u}  = d1[o][v]{u}");
    wfn_->blas()->solve("d'1[O][V]{u}  = d1[O][V]{u}");
    wfn_->blas()->solve("d'2[oo][vv]{u}  = d2[oo][vv]{u}");
    wfn_->blas()->solve("d'2[oO][vV]{u}  = d2[oO][vV]{u}");
    wfn_->blas()->solve("d'2[OO][VV]{u}  = d2[OO][VV]{u}");

    // Compute the shifted denominators
    for (int n = 0; n < wfn_->moinfo()->get_nunique(); n++) {
        int m = wfn_->moinfo()->get_ref_number(n, UniqueRefs);
        std::string shift = to_string(cas_energy - Heff[m][m]);
        wfn_->blas()->solve("d'1[o][v]{" + to_string(m) + "} += " + shift);
        wfn_->blas()->solve("d'1[O][V]{" + to_string(m) + "} += " + shift);
        wfn_->blas()->solve("d'2[oo][vv]{" + to_string(m) + "} += " + shift);
        wfn_->blas()->solve("d'2[oO][vV]{" + to_string(m) + "} += " + shift);
        wfn_->blas()->solve("d'2[OO][VV]{" + to_string(m) + "} += " + shift);
    }

    print_method("  Second-order Mukherjee Multireference Perturbation Theory (Mk-MRPT2)\n    Using the DPD Library");
    outfile->Printf("\n  ------------------------------------------------------------------------------");
    outfile->Printf("\n    @PT  Cycle         Energy           Delta E   ");
    outfile->Printf("\n    @PT                 [Eh]             [Eh]  ");
    outfile->Printf("\n  ------------------------------------------------------------------------------");

    // Start MRPT cycle
    bool converged = false;
    int cycle = 0;
    while (!converged) {
        // Iterate the amps equation
        updater_->zero_internal_amps();
        build_amplitudes();
        update_amps_mkpt2();
        updater_->zero_internal_amps();
        synchronize_amps();

        // Compute the effective Hamiltonian
        build_Heff_mrpt2_diagonal();
        build_Heff_mrpt2_offdiagonal();

        // Compute the energy
        current_energy = c_H_c(wfn_->moinfo()->get_nrefs(), Heff_mrpt2, zeroth_order_eigenvector);
        delta_energy = current_energy - old_energy;
        converged = (std::fabs(delta_energy) < options_.get_double("E_CONVERGENCE"));
        outfile->Printf("\n    @PT %5d   %20.15f  %11.4e", cycle, current_energy, delta_energy);
        old_energy = current_energy;

        if (cycle > options_.get_int("MAXITER")) {
            std::ostringstream oss;
            oss << "The calculation did not converge in " << options_.get_int("MAXITER") << " cycles.\n";
            throw std::runtime_error(oss.str());
        }
        cycle++;
    }
    outfile->Printf("\n  ------------------------------------------------------------------------------");

    double second_order_energy = current_energy;

    // Compute SCS MRPT2
    build_Heff_scs_mrpt2_diagonal();
    build_Heff_mrpt2_offdiagonal();
    current_energy = c_H_c(wfn_->moinfo()->get_nrefs(), Heff_mrpt2, zeroth_order_eigenvector);

    double scs_second_order_energy = current_energy;

    // Diagonalize Heff
    build_Heff_mrpt2_diagonal();
    build_Heff_mrpt2_offdiagonal();
    current_energy = diagonalize_Heff(wfn_->moinfo()->get_root(), wfn_->moinfo()->get_nrefs(), Heff_mrpt2,
                                      right_eigenvector, left_eigenvector, false);

    double pseudo_second_order_energy = current_energy;

    // Diagonalize SCS Heff
    build_Heff_scs_mrpt2_diagonal();
    build_Heff_mrpt2_offdiagonal();
    current_energy = diagonalize_Heff(wfn_->moinfo()->get_root(), wfn_->moinfo()->get_nrefs(), Heff_mrpt2,
                                      right_eigenvector, left_eigenvector, false);

    double scs_pseudo_second_order_energy = current_energy;

    outfile->Printf("\n\n%6c* Mk-MRPT2 total energy         = %20.12f", ' ', second_order_energy);
    outfile->Printf("\n%6c* relaxed Mk-MRPT2 total energy = %20.12f\n", ' ', pseudo_second_order_energy);
    //  outfile->Printf("\n    @SCSPT@        %-5s%-10s  Energy =
    //  %-20.15f",options_.get_str("CORR_ANSATZ").c_str(),options_.get_str("CORR_WFN").c_str(),scs_second_order_energy);
    //  outfile->Printf("\n    @SCSPT-i@      %-5s%-10s  Energy =
    //  %-20.15f",options_.get_str("CORR_ANSATZ").c_str(),options_.get_str("CORR_WFN").c_str(),scs_pseudo_second_order_energy);

    std::string energy_string;
    if (options_.get_str("PT_ENERGY") == "SECOND_ORDER") {
        current_energy = second_order_energy;
        energy_string = "";
    }
    if (options_.get_str("PT_ENERGY") == "SCS_SECOND_ORDER") {
        current_energy = scs_second_order_energy;
        energy_string = "spin-component-scaled ";
    }
    if (options_.get_str("PT_ENERGY") == "PSEUDO_SECOND_ORDER") {
        current_energy = pseudo_second_order_energy;
        energy_string = "pseudo-";
    }
    if (options_.get_str("PT_ENERGY") == "SCS_PSEUDO_SECOND_ORDER") {
        current_energy = scs_pseudo_second_order_energy;
        energy_string = "spin-component-scaled pseudo-";
    }
    wfn_->set_scalar_variable("CURRENT ENERGY", current_energy);
    wfn_->set_scalar_variable("MRPT TOTAL ENERGY", current_energy);
    outfile->Printf("\n\n  Wrote %ssecond order energy to wavefunction", energy_string.c_str());

    //   print_eigensystem(wfn_->moinfo()->get_nrefs(),Heff_mrpt2,right_eigenvector);
    return wfn_->scalar_variable("CURRENT ENERGY");
}

void IDMRPT2::build_Heff_mrpt2_diagonal() {
    // Compute the diagonal elements of the effective Hamiltonian
    // using a simple UCCSD energy expression
    wfn_->blas()->solve("Eaa{u}   = t1[o][v]{u} . fock[o][v]{u}");
    wfn_->blas()->solve("Ebb{u}   = t1[O][V]{u} . fock[O][V]{u}");

    wfn_->blas()->solve("Eaaaa{u} = 1/4 t2[oo][vv]{u} . <[oo]:[vv]>");
    wfn_->blas()->solve("Eabab{u} =     t2[oO][vV]{u} . <[oo]|[vv]>");
    wfn_->blas()->solve("Ebbbb{u} = 1/4 t2[OO][VV]{u} . <[oo]:[vv]>");

    wfn_->blas()->solve("EPT2{u}  = Eaa{u} + Ebb{u} + Eaaaa{u} + Eabab{u} + Ebbbb{u} + ERef{u}");

    for (int n = 0; n < wfn_->moinfo()->get_nrefs(); n++)
        Heff_mrpt2[n][n] = wfn_->blas()->get_scalar("EPT2", wfn_->moinfo()->get_ref_number(n));
}

void IDMRPT2::build_Heff_scs_mrpt2_diagonal() {
    // Compute the diagonal elements of the effective Hamiltonian
    // using a simple UCCSD energy expression with spin-component
    // scaled energies
    wfn_->blas()->solve("Eaa{u}   = t1[o][v]{u} . fock[o][v]{u}");
    wfn_->blas()->solve("Ebb{u}   = t1[O][V]{u} . fock[O][V]{u}");

    wfn_->blas()->solve("Eaaaa{u} = 1/4 t2[oo][vv]{u} . <[oo]:[vv]>");
    wfn_->blas()->solve("Eabab{u} =     t2[oO][vV]{u} . <[oo]|[vv]>");
    wfn_->blas()->solve("Ebbbb{u} = 1/4 t2[OO][VV]{u} . <[oo]:[vv]>");

    wfn_->blas()->solve("EPT2{u}  = Eaa{u} + Ebb{u} + 1/3 Eaaaa{u} + 6/5 Eabab{u} + 1/3 Ebbbb{u} + ERef{u}");

    for (int n = 0; n < wfn_->moinfo()->get_nrefs(); n++)
        Heff_mrpt2[n][n] = wfn_->blas()->get_scalar("EPT2", wfn_->moinfo()->get_ref_number(n));
}

void IDMRPT2::synchronize_amps() {
    START_TIMER("Synchronizing the Amplitudes");

    wfn_->blas()->solve("t1[ov]{u}     = #12# t1[o][v]{u}");
    wfn_->blas()->solve("t1[OV]{u}     = #12# t1[O][V]{u}");

    wfn_->blas()->reduce_spaces("t1_ov[a][v]{u}", "t1[o][v]{u}");
    wfn_->blas()->reduce_spaces("t1_OV[A][V]{u}", "t1[O][V]{u}");

    wfn_->blas()->reduce_spaces("t1_ov[o][a]{u}", "t1[o][v]{u}");
    wfn_->blas()->reduce_spaces("t1_OV[O][A]{u}", "t1[O][V]{u}");

    wfn_->blas()->solve("t2[o][ovv]{u} = #1234# t2[oo][vv]{u}");
    wfn_->blas()->solve("t2[o][OvV]{u} = #1234# t2[oO][vV]{u}");
    wfn_->blas()->solve("t2[O][oVv]{u} = #2143# t2[oO][vV]{u}");
    wfn_->blas()->solve("t2[O][OVV]{u} = #1234# t2[OO][VV]{u}");

    wfn_->blas()->solve("t2[v][voo]{u} = #3412# t2[oo][vv]{u}");
    wfn_->blas()->solve("t2[v][VoO]{u} = #3412# t2[oO][vV]{u}");
    wfn_->blas()->solve("t2[V][vOo]{u} = #4321# t2[oO][vV]{u}");
    wfn_->blas()->solve("t2[V][VOO]{u} = #3412# t2[OO][VV]{u}");

    wfn_->blas()->reduce_spaces("t2_oovv[o][aaa]{u}", "t2[o][ovv]{u}");
    wfn_->blas()->reduce_spaces("t2_OoVv[O][aAa]{u}", "t2[O][oVv]{u}");
    wfn_->blas()->reduce_spaces("t2_oOvV[o][AaA]{u}", "t2[o][OvV]{u}");
    wfn_->blas()->reduce_spaces("t2_OOVV[O][AAA]{o}", "t2[O][OVV]{o}");

    wfn_->blas()->reduce_spaces("t2_oovv[oo][aa]{u}", "t2[oo][vv]{u}");
    wfn_->blas()->reduce_spaces("t2_oOvV[oO][aA]{u}", "t2[oO][vV]{u}");
    wfn_->blas()->reduce_spaces("t2_OOVV[OO][AA]{o}", "t2[OO][VV]{o}");

    wfn_->blas()->reduce_spaces("t2_oovv[a][ovv]{u}", "t2[o][ovv]{u}");
    wfn_->blas()->reduce_spaces("t2_oOvV[a][OvV]{u}", "t2[o][OvV]{u}");
    wfn_->blas()->reduce_spaces("t2_vvoo[a][voo]{u}", "t2[v][voo]{u}");
    wfn_->blas()->reduce_spaces("t2_vVoO[a][VoO]{u}", "t2[v][VoO]{u}");
    wfn_->blas()->reduce_spaces("t2_OOVV[A][OVV]{o}", "t2[O][OVV]{o}");
    wfn_->blas()->reduce_spaces("t2_OoVv[A][oVv]{o}", "t2[O][oVv]{o}");
    wfn_->blas()->reduce_spaces("t2_VVOO[A][VOO]{o}", "t2[V][VOO]{o}");
    wfn_->blas()->reduce_spaces("t2_VvOo[A][vOo]{o}", "t2[V][vOo]{o}");

    wfn_->blas()->reduce_spaces("t2_oovv[aa][vv]{u}", "t2[oo][vv]{u}");
    wfn_->blas()->reduce_spaces("t2_oOvV[aA][vV]{u}", "t2[oO][vV]{u}");
    wfn_->blas()->reduce_spaces("t2_OOVV[AA][VV]{o}", "t2[OO][VV]{o}");

    wfn_->blas()->reduce_spaces("t2_vvoo[v][aaa]{u}", "t2[v][voo]{u}");
    wfn_->blas()->reduce_spaces("t2_VvOo[V][aAa]{u}", "t2[V][vOo]{u}");
    wfn_->blas()->reduce_spaces("t2_vVoO[v][AaA]{u}", "t2[v][VoO]{u}");
    wfn_->blas()->reduce_spaces("t2_VVOO[V][AAA]{o}", "t2[V][VOO]{o}");

    wfn_->blas()->reduce_spaces("t2_oovv[ao][av]{u}", "t2[oo][vv]{u}");
    wfn_->blas()->solve("t2_ovov[aa][ov]{u} = #1324# t2_oovv[ao][av]{u}");

    wfn_->blas()->reduce_spaces("t2_oOvV[oA][vA]{u}", "t2[oO][vV]{u}");
    wfn_->blas()->solve("t2_ovOV[ov][AA]{u} = #1324# t2_oOvV[oA][vA]{u}");

    wfn_->blas()->reduce_spaces("t2_oOvV[aO][aV]{u}", "t2[oO][vV]{u}");
    wfn_->blas()->solve("t2_ovOV[aa][OV]{u} = #1324# t2_oOvV[aO][aV]{u}");

    wfn_->blas()->reduce_spaces("t2_OOVV[AO][AV]{u}", "t2[OO][VV]{u}");
    wfn_->blas()->solve("t2_OVOV[AA][OV]{u} = #1324# t2_OOVV[AO][AV]{u}");

    wfn_->blas()->reduce_spaces("t2_oOvV[aO][vA]{u}", "t2[oO][vV]{u}");
    wfn_->blas()->solve("t2_oVOv[aA][Ov]{u} = #1342# t2_oOvV[aO][vA]{u}");

    wfn_->blas()->reduce_spaces("t2_oOvV[oA][aV]{u}", "t2[oO][vV]{u}");
    wfn_->blas()->solve("t2_oVOv[oV][Aa]{u} = #1342# t2_oOvV[oA][aV]{u}");

    END_TIMER("Synchronizing the Amplitudes");
}

void IDMRPT2::build_amplitudes() {
    // These are required by the t1 amplitude equations
    build_t1_ia_amplitudes();
    build_t1_IA_amplitudes();
    build_t2_iJaB_amplitudes();
    build_t2_ijab_amplitudes();
    build_t2_IJAB_amplitudes();
}

void IDMRPT2::update_amps_mkpt2() {
    for (int i = 0; i < wfn_->moinfo()->get_nunique(); i++) {
        int unique_i = wfn_->moinfo()->get_ref_number(i, UniqueRefs);
        std::string i_str = to_string(unique_i);
        // Form the coupling terms
        for (int j = 0; j < wfn_->moinfo()->get_nrefs(); j++) {
            int unique_j = wfn_->moinfo()->get_ref_number(j);
            std::string j_str = to_string(unique_j);
            double term = zeroth_order_eigenvector[j] / zeroth_order_eigenvector[unique_i];
            if (std::fabs(term) > 1.0e5) term = 0.0;
            wfn_->blas()->set_scalar("factor_mk", unique_j, Heff[unique_i][j] * term);
            if (unique_i != j) {
                if (j == unique_j) {
                    wfn_->blas()->solve("t1_eqns[o][v]{" + i_str + "} += factor_mk{" + j_str + "} t1[o][v]{" + j_str +
                                        "}");
                    wfn_->blas()->solve("t1_eqns[O][V]{" + i_str + "} += factor_mk{" + j_str + "} t1[O][V]{" + j_str +
                                        "}");
                } else {
                    wfn_->blas()->solve("t1_eqns[o][v]{" + i_str + "} += factor_mk{" + j_str + "} t1[O][V]{" + j_str +
                                        "}");
                    wfn_->blas()->solve("t1_eqns[O][V]{" + i_str + "} += factor_mk{" + j_str + "} t1[o][v]{" + j_str +
                                        "}");
                }
            }
        }
        // Update t1 for reference i
        wfn_->blas()->solve("t1[o][v]{" + i_str + "} = t1_eqns[o][v]{" + i_str + "} / d'1[o][v]{" + i_str + "}");
        wfn_->blas()->solve("t1[O][V]{" + i_str + "} = t1_eqns[O][V]{" + i_str + "} / d'1[O][V]{" + i_str + "}");
        updater_->zero_internal_amps();

        // Add the contribution from the other references
        for (int j = 0; j < wfn_->moinfo()->get_nrefs(); j++) {
            int unique_j = wfn_->moinfo()->get_ref_number(j);
            std::string j_str = to_string(unique_j);
            double term = zeroth_order_eigenvector[j] / zeroth_order_eigenvector[unique_i];
            if (std::fabs(term) > 1.0e5) term = 0.0;
            wfn_->blas()->set_scalar("factor_mk", unique_j, Heff[unique_i][j] * term);
            if (unique_i != j) {
                if (j == unique_j) {
                    // aaaa case
                    // + t_ij^ab(nu/mu)
                    wfn_->blas()->solve("t2_eqns[oo][vv]{" + i_str + "} += factor_mk{" + j_str + "} t2[oo][vv]{" +
                                        j_str + "}");

                    // abab case
                    // + t_ij^ab(nu/mu)
                    wfn_->blas()->solve("t2_eqns[oO][vV]{" + i_str + "} += factor_mk{" + j_str + "} t2[oO][vV]{" +
                                        j_str + "}");

                    // bbbb case
                    // + t_ij^ab(nu/mu)
                    wfn_->blas()->solve("t2_eqns[OO][VV]{" + i_str + "} += factor_mk{" + j_str + "} t2[OO][VV]{" +
                                        j_str + "}");
                } else {
                    // aaaa case
                    // + t_ij^ab(nu/mu)
                    wfn_->blas()->solve("t2_eqns[oo][vv]{" + i_str + "} += factor_mk{" + j_str + "} t2[OO][VV]{" +
                                        j_str + "}");

                    // abab case
                    // + t_ij^ab(nu/mu)
                    wfn_->blas()->solve("t2_eqns[oO][vV]{" + i_str + "} += #2143# factor_mk{" + j_str +
                                        "} t2[oO][vV]{" + j_str + "}");

                    // bbbb case
                    // + t_ij^ab(nu/mu)
                    wfn_->blas()->solve("t2_eqns[OO][VV]{" + i_str + "} += factor_mk{" + j_str + "} t2[oo][vv]{" +
                                        j_str + "}");
                }
            }
        }

        wfn_->blas()->solve("t2[oo][vv]{" + i_str + "} = t2_eqns[oo][vv]{" + i_str + "} / d'2[oo][vv]{" + i_str + "}");
        wfn_->blas()->solve("t2[oO][vV]{" + i_str + "} = t2_eqns[oO][vV]{" + i_str + "} / d'2[oO][vV]{" + i_str + "}");
        wfn_->blas()->solve("t2[OO][VV]{" + i_str + "} = t2_eqns[OO][VV]{" + i_str + "} / d'2[OO][VV]{" + i_str + "}");
    }
}

void IDMRPT2::read_mrpt2_integrals() {
    START_TIMER("Reading the MRPT2 integrals");

    // CCSort reads the one and two electron integrals
    // and creates the Fock matrices
    std::make_shared<CCSort>(wfn_, mrpt2_sort);

    END_TIMER("Reading the MRPT2 integrals");
}

}  // namespace psimrcc
}  // namespace psi
