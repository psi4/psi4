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

#include <armadillo>

#include "psi4/psi4-dec.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/efpmultipolepotential.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/potentialint.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"

// #include <cppe/libcppe.hh>
#include <cppe/utils/potfile_reader.hh>
#include <cppe/utils/pot_manipulation.hh>
#include <cppe/core/multipole_expansion.hh>

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

// TODO: generate these factors!
namespace {
std::vector<double> mult_coeffs = {
    1.0,                       // q
    1.0, 1.0, 1.0,             // mu
    0.5, 0.5, 0.5, 1., 1., 1.  // theta
};

std::vector<int> from_alphabetical_to_psi4 = {
    0, 1, 2, 3, 4, 7, 9, 5, 6, 8  // 0, 3, 5, 1, 2, 4   += 4
};
}  // unnamed namespace

PeState::PeState(libcppe::PeOptions options, std::shared_ptr<BasisSet> basisset)
    : basisset_(basisset),
      cppe_state_(libcppe::CppeState(options, make_molecule(basisset_->molecule()), *(outfile->stream()))),
      int_helper_(PeIntegralHelper(basisset)) {
    potentials_ = cppe_state_.get_potentials();
    cppe_state_.calculate_static_energies_and_fields();
}

std::pair<double, SharedMatrix> PeState::compute_pe_contribution(const SharedMatrix& Dm, CalcType type,
                                                                 bool subtract_scf_density) {
    // build the electrostatics operator only once!
    if (!iteration && (type == PeState::CalcType::total)) {
        V_es_ = std::make_shared<Matrix>("V_es", basisset_->nbf(), basisset_->nbf());
        for (auto& p : potentials_) {
            Vector3 site(p.m_x, p.m_y, p.m_z);
            std::vector<double> moments;
            for (auto& m : p.get_multipoles()) {
                m.remove_trace();
                if (m.m_k > 2) throw std::runtime_error("PE is only implemented up to quadrupoles.");
                for (auto d : m.get_values()) moments.push_back(d);
            }
            V_es_->add(int_helper_.compute_multipole_potential_integrals(site, 2, moments));
        }
        cppe_state_.set_es_operator(arma::mat(V_es_->get_pointer(), V_es_->nrow(), V_es_->ncol(), true));
        cppe_state_.es_operator_copy().save("V_es_psi4.txt", arma::raw_ascii);
    }

    SharedMatrix D = Dm;
    if (subtract_scf_density && type == PeState::CalcType::electronic_only) {
        D->subtract(D_scf_);
    }
    cppe_state_.update_energies(arma::mat(D->get_pointer(), D->nrow(), D->ncol(), true));

    size_t n_sitecoords = 3 * cppe_state_.get_polarizable_site_number();

    // induction operator
    auto V_ind = std::make_shared<Matrix>("V_ind", basisset_->nbf(), basisset_->nbf());
    if (n_sitecoords) {
        int current_polsite = 0;
        arma::vec elec_fields(n_sitecoords, arma::fill::zeros);
        for (auto& p : potentials_) {
            if (!p.is_polarizable()) continue;
            Vector3 site(p.m_x, p.m_y, p.m_z);
            Vector elec_fields_s = int_helper_.compute_field(site, D);
            elec_fields[current_polsite * 3] = elec_fields_s.get(0);
            elec_fields[current_polsite * 3 + 1] = elec_fields_s.get(1);
            elec_fields[current_polsite * 3 + 2] = elec_fields_s.get(2);
            current_polsite += 1;
        }
        // std::cout << "elec fields : " << n_sitecoords << std::endl << elec_fields << std::endl;
        cppe_state_.update_induced_moments(elec_fields, iteration, type == PeState::CalcType::electronic_only);
        arma::vec induced_moments = cppe_state_.get_induced_moments();

        current_polsite = 0;
        for (auto& p : potentials_) {
            if (!p.is_polarizable()) continue;
            Vector3 site(p.m_x, p.m_y, p.m_z);
            SharedMatrix V_ind_s = int_helper_.compute_field_integrals(
                site, induced_moments.subvec(3 * current_polsite, 3 * current_polsite + 2));
            V_ind->add(V_ind_s);
            current_polsite += 1;
        }
    }

    if (type == PeState::CalcType::total) {
        V_ind->add(V_es_);
        D_scf_ = D;
    }
    // TODO: only for debugging
    // cppe_state_.print_summary();

    if (type == PeState::CalcType::total) iteration++;

    double energy_result =
        (type == PeState::CalcType::total ? cppe_state_.get_current_energies().get_total_energy()
                                          : cppe_state_.get_current_energies().get("Polarization/Electronic"));
    return std::pair<double, SharedMatrix>(energy_result, V_ind);
}

void PeState::print_energy_summary() { cppe_state_.print_summary(); }

SharedMatrix PeIntegralHelper::compute_multipole_potential_integrals(Vector3 site, int order,
                                                                     std::vector<double>& moments) {
    auto integrals = std::make_shared<IntegralFactory>(basisset_, basisset_, basisset_, basisset_);
    std::shared_ptr<OneBodyAOInt> multipole_integrals =
        static_cast<std::shared_ptr<OneBodyAOInt>>(integrals->ao_efp_multipole_potential());

    // TODO: only compute up to the needed order!
    std::vector<SharedMatrix> mult;
    mult.push_back(std::make_shared<Matrix>("AO EFP Charge 0", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Dipole X", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Dipole Y", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Dipole Z", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Quadrupole XX", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Quadrupole YY", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Quadrupole ZZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Quadrupole XY", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Quadrupole XZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Quadrupole YZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole XXX", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole YYY", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole ZZZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole XXY", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole XXZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole XYY", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole YYZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole XZZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole YZZ", basisset_->nbf(), basisset_->nbf()));
    mult.push_back(std::make_shared<Matrix>("AO EFP Octupole XYZ", basisset_->nbf(), basisset_->nbf()));
    multipole_integrals->set_origin(site);
    multipole_integrals->compute(mult);

    auto res = std::make_shared<Matrix>("result", basisset_->nbf(), basisset_->nbf());

    int total_components = 0;
    for (int i = 0; i <= order; ++i) {
        total_components += libcppe::multipole_components(i);
    }

    for (int l = 0; l < total_components; ++l) {
        // std::cout << "---" << std::setprecision(12) << std::endl;
        // std::cout << "L: " << l << std::endl;
        // std::cout << "coeff = " << mult_coeffs[l] << " moment = " << moments[from_alphabetical_to_psi4[l]] <<
        // std::endl;
        mult[l]->scale(-1.0 * mult_coeffs[l] * moments[from_alphabetical_to_psi4[l]]);
        res->add(mult[l]);
    }

    return res;
}

SharedMatrix PeIntegralHelper::compute_field_integrals(Vector3 site, arma::vec moment) {
    std::vector<SharedMatrix> fields;
    fields.push_back(std::make_shared<Matrix>("AO Field X", basisset_->nbf(), basisset_->nbf()));
    fields.push_back(std::make_shared<Matrix>("AO Field Y", basisset_->nbf(), basisset_->nbf()));
    fields.push_back(std::make_shared<Matrix>("AO Field Z", basisset_->nbf(), basisset_->nbf()));

    auto integrals = std::make_shared<IntegralFactory>(basisset_, basisset_, basisset_, basisset_);
    std::shared_ptr<OneBodyAOInt> field_integrals_ =
        static_cast<std::shared_ptr<OneBodyAOInt>>(integrals->electric_field());
    field_integrals_->set_origin(site);
    field_integrals_->compute(fields);

    fields[0]->scale(-1.0 * moment[0]);
    fields[1]->scale(-1.0 * moment[1]);
    fields[2]->scale(-1.0 * moment[2]);

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
    std::shared_ptr<OneBodyAOInt> field_integrals_ =
        static_cast<std::shared_ptr<OneBodyAOInt>>(integrals->electric_field());
    field_integrals_->set_origin(site);
    field_integrals_->compute(fields);

    Vector efields("efields", 3);
    efields.set(0, D->vector_dot(fields[0]));
    efields.set(1, D->vector_dot(fields[1]));
    efields.set(2, D->vector_dot(fields[2]));
    return efields;
}

}  // namespace psi

#endif