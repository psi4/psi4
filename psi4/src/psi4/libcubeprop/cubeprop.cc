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

#include "psi4/psi4-dec.h"

#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/process.h"

#include "cubeprop.h"
#include "csg.h"

#include <tuple>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

CubeProperties::CubeProperties(SharedWavefunction wfn) : options_(Process::environment.options) {
    basisset_ = wfn->basisset();
    if (wfn->basisset_exists("DF_BASIS_SCF")) auxiliary_ = wfn->get_basisset("DF_BASIS_SCF");

    Ca_ = wfn->Ca_subset("AO", "ALL");
    Da_ = wfn->Da_subset("AO");

    if (wfn->same_a_b_orbs()) {
        Cb_ = Ca_;
    } else {
        Cb_ = wfn->Cb_subset("AO", "ALL");
    }

    if (wfn->same_a_b_dens()) {
        Db_ = Da_;
    } else {
        Db_ = wfn->Db_subset("AO");
    }

    int nirrep = wfn->nirrep();
    Dimension nmopi = wfn->nmopi();
    nalpha_ = 0;
    // Gather orbital information
    for (int h = 0; h < nirrep; h++) {
        nalpha_ += wfn->nalphapi()[h];
        for (int i = 0; i < (int)nmopi[h]; i++) {
            info_a_.push_back(std::tuple<double, int, int>(wfn->epsilon_a()->get(h, i), i, h));
        }
    }
    std::sort(info_a_.begin(), info_a_.end(), std::less<std::tuple<double, int, int> >());  // Sort as in wfn
    nbeta_ = 0;
    for (int h = 0; h < nirrep; h++) {
        nbeta_ += wfn->nbetapi()[h];
        for (int i = 0; i < (int)nmopi[h]; i++) {
            info_b_.push_back(std::tuple<double, int, int>(wfn->epsilon_b()->get(h, i), i, h));
        }
    }
    std::sort(info_b_.begin(), info_b_.end(), std::less<std::tuple<double, int, int> >());  // Sort as in wfn

    common_init();
}
CubeProperties::~CubeProperties() {}
void CubeProperties::common_init() {
    grid_ = std::make_shared<CubicScalarGrid>(basisset_, options_);
    grid_->set_filepath(options_.get_str("CUBEPROP_FILEPATH"));
    grid_->set_auxiliary_basis(auxiliary_);
}
void CubeProperties::print_header() {
    outfile->Printf("  ==> One Electron Grid Properties (v2.0) <==\n\n");
    grid_->print_header();
}
void CubeProperties::raw_compute_properties() {
    print_header();

    for (size_t ind = 0; ind < options_["CUBEPROP_TASKS"].size(); ind++) {
        std::string task = options_["CUBEPROP_TASKS"][ind].to_string();

        if (task == "DENSITY") {
            std::shared_ptr<Matrix> Dt(Da_->clone());
            std::shared_ptr<Matrix> Ds(Da_->clone());
            Dt->copy(Da_);
            Ds->copy(Da_);
            Dt->add(Db_);
            Ds->subtract(Db_);
            compute_density(Dt, "Dt");
            compute_density(Ds, "Ds");
            compute_density(Da_, "Da");
            compute_density(Db_, "Db");
        } else if (task == "ESP") {
            std::shared_ptr<Matrix> Dt(Da_->clone());
            Dt->copy(Da_);
            Dt->add(Db_);
            compute_esp(Dt);
        } else if (task == "ORBITALS") {
            std::vector<int> indsa0;
            std::vector<int> indsb0;

            if (options_["CUBEPROP_ORBITALS"].size() == 0) {
                for (int ind = 0; ind < Ca_->colspi()[0]; ind++) {
                    indsa0.push_back(ind);
                }
                for (int ind = 0; ind < Cb_->colspi()[0]; ind++) {
                    indsb0.push_back(ind);
                }
            } else {
                for (size_t ind = 0; ind < options_["CUBEPROP_ORBITALS"].size(); ind++) {
                    int val = options_["CUBEPROP_ORBITALS"][ind].to_integer();
                    if (val > 0) {
                        indsa0.push_back(std::abs(val) - 1);
                    } else {
                        indsb0.push_back(std::abs(val) - 1);
                    }
                }
            }
            std::vector<std::string> labelsa;
            std::vector<std::string> labelsb;
            CharacterTable ct = basisset_->molecule()->point_group()->char_table();
            for (size_t ind = 0; ind < indsa0.size(); ++ind) {
                int i = std::get<1>(info_a_[indsa0[ind]]);
                int h = std::get<2>(info_a_[indsa0[ind]]);
                labelsa.push_back(std::to_string(i + 1) + "-" + ct.gamma(h).symbol());
            }
            for (size_t ind = 0; ind < indsb0.size(); ++ind) {
                int i = std::get<1>(info_b_[indsb0[ind]]);
                int h = std::get<2>(info_b_[indsb0[ind]]);
                labelsb.push_back(std::to_string(i + 1) + "-" + ct.gamma(h).symbol());
            }
            if (indsa0.size()) compute_orbitals(Ca_, indsa0, labelsa, "Psi_a");
            if (indsb0.size()) compute_orbitals(Cb_, indsb0, labelsb, "Psi_b");
        } else if (task == "FRONTIER_ORBITALS") {
            std::vector<int> indsa0;
            std::vector<int> indsb0;
            std::vector<std::string> labelsa;
            std::vector<std::string> labelsb;
            CharacterTable ct = basisset_->molecule()->point_group()->char_table();
            if (nalpha_ == nbeta_) {
                indsa0.push_back(nalpha_-1);
                labelsa.push_back(std::to_string(std::get<1>(info_a_[nalpha_-1]) + 1) + "-" +
                                  ct.gamma(std::get<2>(info_a_[nalpha_-1])).symbol() + "_HOMO");
                indsa0.push_back(nalpha_);
                labelsa.push_back(std::to_string(std::get<1>(info_a_[nalpha_]) + 1) + "-" +
                                  ct.gamma(std::get<2>(info_a_[nalpha_])).symbol() + "_LUMO");
            } else {
                int orb_index = nalpha_;
                indsa0.push_back(orb_index);
                labelsa.push_back(std::to_string(std::get<1>(info_a_[orb_index]) + 1) + "-" +
                                  ct.gamma(std::get<2>(info_a_[orb_index])).symbol() + "_LVMO");
                indsb0.push_back(orb_index);
                labelsb.push_back(std::to_string(std::get<1>(info_b_[orb_index]) + 1) + "-" +
                                  ct.gamma(std::get<2>(info_b_[orb_index])).symbol() + "_LVMO");
                for (int i = 1; i <= nalpha_ - nbeta_; i++) {
                    orb_index = nalpha_ - i;
                    indsa0.push_back(orb_index);
                    labelsa.push_back(std::to_string(std::get<1>(info_a_[orb_index]) + 1) + "-" +
                                  ct.gamma(std::get<2>(info_a_[orb_index])).symbol() + "_SOMO");
                }
                orb_index = nbeta_ - 1;
                indsa0.push_back(orb_index);
                labelsa.push_back(std::to_string(std::get<1>(info_a_[orb_index]) + 1) + "-" +
                                  ct.gamma(std::get<2>(info_a_[orb_index])).symbol() + "_DOMO");
                orb_index = nbeta_ - 1;
                indsb0.push_back(orb_index);
                labelsb.push_back(std::to_string(std::get<1>(info_b_[orb_index]) + 1) + "-" +
                                  ct.gamma(std::get<2>(info_b_[orb_index])).symbol() + "_DOMO");
            }
            if (indsa0.size()) compute_orbitals(Ca_, indsa0, labelsa, "Psi_a");
            if (indsb0.size()) compute_orbitals(Cb_, indsb0, labelsb, "Psi_b");
        } else if (task == "DUAL_DESCRIPTOR") {
            // Calculates the dual descriptor from frontier molecular orbitals.
            // The dual descriptor is a good measure of electro-/nucleophilicity:
            // f^2(r) = rho_lumo(r) - rho_homo(r)
            // See 10.1021/jp046577a and 10.1007/s10910-014-0437-7
            std::vector<int> indsa0;
            if (nalpha_ != nbeta_) {
                throw PSIEXCEPTION(task + "is not implemented for open-shell systems");
            } else {
                indsa0.push_back(nalpha_);
                indsa0.push_back(nalpha_-1);
                std::stringstream ss;
                ss << "DUAL_" << (nalpha_+1) << "_LUMO-" << nalpha_ << "_HOMO";
                compute_difference(Ca_, indsa0, ss.str(), true);
            }
        } else if (task == "BASIS_FUNCTIONS") {
            std::vector<int> inds0;
            if (options_["CUBEPROP_BASIS_FUNCTIONS"].size() == 0) {
                for (int ind = 0; ind < basisset_->nbf(); ind++) {
                    inds0.push_back(ind);
                }
            } else {
                for (size_t ind = 0; ind < options_["CUBEPROP_BASIS_FUNCTIONS"].size(); ind++) {
                    inds0.push_back(options_["CUBEPROP_BASIS_FUNCTIONS"][ind].to_integer() - 1);
                }
            }
            compute_basis_functions(inds0, "Phi");
        } else if (task == "LOL") {
            compute_LOL(Da_, "LOLa");
            compute_LOL(Db_, "LOLb");
        } else if (task == "ELF") {
            compute_ELF(Da_, "ELFa");
            compute_ELF(Db_, "ELFb");
        } else {
            throw PSIEXCEPTION(task + "is an unrecognized PROPERTY_TASKS value");
        }
    }
}
void CubeProperties::compute_density(std::shared_ptr<Matrix> D, const std::string& key) {
    grid_->compute_density(D, key);
}
void CubeProperties::compute_esp(std::shared_ptr<Matrix> Dt, const std::vector<double>& w) {
    grid_->compute_density(Dt, "Dt");
    grid_->compute_esp(Dt, w, "ESP");
}
void CubeProperties::compute_orbitals(std::shared_ptr<Matrix> C, const std::vector<int>& indices,
                                      const std::vector<std::string>& labels, const std::string& key) {
    grid_->compute_orbitals(C, indices, labels, key);
}
void CubeProperties::compute_difference(std::shared_ptr<Matrix> C, const std::vector<int>& indices,
                                      const std::string& label, bool square) {
    grid_->compute_difference(C, indices, label, square);
}
void CubeProperties::compute_basis_functions(const std::vector<int>& indices, const std::string& key) {
    grid_->compute_basis_functions(indices, key);
}
void CubeProperties::compute_LOL(std::shared_ptr<Matrix> D, const std::string& key) { grid_->compute_LOL(D, key); }
void CubeProperties::compute_ELF(std::shared_ptr<Matrix> D, const std::string& key) { grid_->compute_ELF(D, key); }
}  // namespace psi
