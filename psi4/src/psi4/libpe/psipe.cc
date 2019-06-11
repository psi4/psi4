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
#ifdef USING_cppe

#include <iomanip>

#include "psi4/psi4-dec.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/multipolepotential.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/potentialint.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <cppe/core/math.hh>

#include "psipe.h"

// TODO: implement output stream option in CPPE for printing

namespace psi {

libcppe::Molecule make_molecule(std::shared_ptr<Molecule> molecule) {
    libcppe::Molecule mol;
    Matrix geom = molecule->geometry();
    for (int i = 0; i < molecule->natom(); ++i) {
        libcppe::Atom a(static_cast<int>(molecule->Z(i)), geom.get(i, 0), geom.get(i, 1), geom.get(i, 2));
        mol.push_back(a);
    }
    return mol;
}

namespace {
std::vector<double> mult_coeffs_k(int max_k = 0) {
    std::vector<double> mult_coeffs;
    for (int k = 0; k <= max_k; ++k) {
        auto prefs = libcppe::prefactors(k);
        mult_coeffs.insert(std::end(mult_coeffs), std::begin(prefs), std::end(prefs));
    }
    return mult_coeffs;
}
std::vector<double> mult_coeffs = mult_coeffs_k(2);

int have_moment(int k, int max_k) { return k <= max_k; }
}  // unnamed namespace

PeState::PeState(libcppe::PeOptions options, std::shared_ptr<BasisSet> basisset)
    : cppe_state_(libcppe::CppeState(options, make_molecule(basisset->molecule()), *(outfile->stream()))),
      int_helper_(PeIntegralHelper(basisset)) {
    potentials_ = cppe_state_.get_potentials();
    cppe_state_.calculate_static_energies_and_fields();
    nbf_ = basisset->nbf();

    V_es_ = std::make_shared<Matrix>("V_es", nbf_, nbf_);
    for (auto& p : potentials_) {
        Vector3 site(p.m_x, p.m_y, p.m_z);
        std::vector<double> moments;
        int max_k = 0;
        for (auto& m : p.get_multipoles()) {
            if (m.m_k > max_k) max_k = m.m_k;
            m.remove_trace();
            if (m.m_k > 2) throw std::runtime_error("PE is only implemented up to quadrupoles.");
            for (auto d : m.get_values()) moments.push_back(d);
        }
        V_es_->add(int_helper_.compute_multipole_potential_integrals(site, max_k, moments));
    }
}

std::pair<double, SharedMatrix> PeState::compute_pe_contribution(const SharedMatrix& Dm, CalcType type) {
    cppe_state_.get_energies().set("Electrostatic/Electronic", Dm->vector_dot(V_es_));

    size_t n_sitecoords = 3 * cppe_state_.get_polarizable_site_number();

    // induction operator
    auto V_ind = std::make_shared<Matrix>("V_ind", nbf_, nbf_);
    if (n_sitecoords) {
        int current_polsite = 0;
        Eigen::VectorXd elec_fields(n_sitecoords);
        for (auto& p : potentials_) {
            if (!p.is_polarizable()) continue;
            Vector3 site(p.m_x, p.m_y, p.m_z);
            Vector elec_fields_s = int_helper_.compute_field(site, Dm);
            elec_fields[current_polsite * 3] = elec_fields_s.get(0);
            elec_fields[current_polsite * 3 + 1] = elec_fields_s.get(1);
            elec_fields[current_polsite * 3 + 2] = elec_fields_s.get(2);
            current_polsite += 1;
        }
        cppe_state_.update_induced_moments(elec_fields, type == PeState::CalcType::electronic_only);
        Eigen::VectorXd induced_moments = cppe_state_.get_induced_moments_vec();

        current_polsite = 0;
        for (auto& p : potentials_) {
            if (!p.is_polarizable()) continue;
            Vector3 site(p.m_x, p.m_y, p.m_z);
            SharedMatrix V_ind_s =
                int_helper_.compute_field_integrals(site, induced_moments.segment(3 * current_polsite, 3));
            V_ind->add(V_ind_s);
            current_polsite += 1;
        }
    }

    if (type == PeState::CalcType::total) {
        V_ind->add(V_es_);
    }

    double energy_result =
        (type == PeState::CalcType::total ? cppe_state_.get_energies().get_total_energy()
                                          : cppe_state_.get_energies().get("Polarization/Electronic"));
    return std::pair<double, SharedMatrix>(energy_result, V_ind);
}

void PeState::print_energy_summary() { cppe_state_.print_summary(); }

SharedMatrix PeIntegralHelper::compute_multipole_potential_integrals(Vector3 site, int order,
                                                                     std::vector<double>& moments) {
    if (order > 2) {
        throw std::runtime_error(
            "Can only compute multipole potential through "
            "second order.");
    }
    auto integrals = std::make_shared<IntegralFactory>(basisset_, basisset_, basisset_, basisset_);
    auto multipole_integrals = integrals->ao_multipole_potential(order);

    std::vector<SharedMatrix> mult;
    int total_components = 0;
    for (int i = 0; i <= order; ++i) {
        total_components += libcppe::multipole_components(i);
    }

    mult.push_back(std::make_shared<Matrix>("AO Charge Potential 0", basisset_->nbf(), basisset_->nbf()));
    if (have_moment(1, order)) {
        mult.push_back(std::make_shared<Matrix>("AO Dipole Potential X", basisset_->nbf(), basisset_->nbf()));
        mult.push_back(std::make_shared<Matrix>("AO Dipole Potential Y", basisset_->nbf(), basisset_->nbf()));
        mult.push_back(std::make_shared<Matrix>("AO Dipole Potential Z", basisset_->nbf(), basisset_->nbf()));
    }
    if (have_moment(2, order)) {
        mult.push_back(std::make_shared<Matrix>("AO Quadrupole Potential XX", basisset_->nbf(), basisset_->nbf()));
        mult.push_back(std::make_shared<Matrix>("AO Quadrupole Potential XY", basisset_->nbf(), basisset_->nbf()));
        mult.push_back(std::make_shared<Matrix>("AO Quadrupole Potential XZ", basisset_->nbf(), basisset_->nbf()));
        mult.push_back(std::make_shared<Matrix>("AO Quadrupole Potential YY", basisset_->nbf(), basisset_->nbf()));
        mult.push_back(std::make_shared<Matrix>("AO Quadrupole Potential YZ", basisset_->nbf(), basisset_->nbf()));
        mult.push_back(std::make_shared<Matrix>("AO Quadrupole Potential ZZ", basisset_->nbf(), basisset_->nbf()));
    }
    multipole_integrals->set_origin(site);
    multipole_integrals->compute(mult);

    auto res = std::make_shared<Matrix>("result", basisset_->nbf(), basisset_->nbf());
    for (int l = 0; l < total_components; ++l) {
        mult[l]->scale(mult_coeffs[l] * moments[l]);
        res->add(mult[l]);
    }
    return res;
}

SharedMatrix PeIntegralHelper::compute_field_integrals(Vector3 site, Eigen::VectorXd moment) {
    std::vector<SharedMatrix> fields;
    fields.push_back(std::make_shared<Matrix>("AO Field X", basisset_->nbf(), basisset_->nbf()));
    fields.push_back(std::make_shared<Matrix>("AO Field Y", basisset_->nbf(), basisset_->nbf()));
    fields.push_back(std::make_shared<Matrix>("AO Field Z", basisset_->nbf(), basisset_->nbf()));

    auto integrals = std::make_shared<IntegralFactory>(basisset_, basisset_, basisset_, basisset_);
    auto field_integrals_ = integrals->electric_field();
    field_integrals_->set_origin(site);
    field_integrals_->compute(fields);

    for (int i = 0; i < 3; ++i) {
        fields[i]->scale(-1.0 * moment[i]);
    }

    fields[0]->add(fields[1]);
    fields[0]->add(fields[2]);

    return fields[0];
}

Vector PeIntegralHelper::compute_field(Vector3 site, const SharedMatrix& D) {
    std::vector<SharedMatrix> fields;
    fields.push_back(std::make_shared<Matrix>("AO Field X", basisset_->nbf(), basisset_->nbf()));
    fields.push_back(std::make_shared<Matrix>("AO Field Y", basisset_->nbf(), basisset_->nbf()));
    fields.push_back(std::make_shared<Matrix>("AO Field Z", basisset_->nbf(), basisset_->nbf()));

    auto integrals = std::make_shared<IntegralFactory>(basisset_, basisset_, basisset_, basisset_);
    auto field_integrals_ = integrals->electric_field();
    field_integrals_->set_origin(site);
    field_integrals_->compute(fields);

    Vector efields("efields", 3);
    for (int i = 0; i < 3; ++i) {
        efields[i] = D->vector_dot(fields[i]);
    }
    return efields;
}

}  // namespace psi

#endif
