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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"
#include "psi4/physconst.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/vector3.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/quadrupole.h"
#include "psi4/libmints/multipolesymmetry.h"
#include "psi4/libmints/shellrotation.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/electricfield.h"
#include "psi4/libmints/electrostatic.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/multipoles.h"
#include "psi4/libmints/dipole.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <utility>
#include <fstream>
#include <regex>
#include <tuple>
#include <functional>
#include <memory>
#include <map>
#include <vector>

namespace psi {

Prop::Prop(std::shared_ptr<Wavefunction> wfn) : wfn_(wfn) {
    if (wfn_.get() == nullptr) throw PSIEXCEPTION("Prop: Wavefunction is null");
    set_wavefunction(wfn_);
}
Prop::~Prop() {}
void Prop::common_init() { set_wavefunction(wfn_); }
void Prop::set_wavefunction(std::shared_ptr<Wavefunction> wfn) {
    wfn_ = wfn;

    basisset_ = wfn_->basisset();
    same_orbs_ = wfn_->same_a_b_orbs();
    same_dens_ = wfn_->same_a_b_dens();

    integral_ = std::make_shared<IntegralFactory>(basisset_, basisset_, basisset_, basisset_);

    auto pet = std::make_shared<PetiteList>(basisset_, integral_);
    AO2USO_ = pet->aotoso();
    factory_ = wfn_->matrix_factory();

    epsilon_a_ = wfn_->epsilon_a();
    Ca_so_ = wfn_->Ca();
    Da_so_ = wfn_->Da();

    if (same_dens_) {
        Db_so_ = Da_so_;
    } else {
        Db_so_ = wfn_->Db();
    }

    if (same_orbs_) {
        epsilon_b_ = epsilon_a_;
        Cb_so_ = Ca_so_;
    } else {
        epsilon_b_ = wfn_->epsilon_b();
        Cb_so_ = wfn_->Cb();
    }
}
void Prop::set_restricted(bool restricted) {
    if (restricted == same_orbs_) return;

    same_orbs_ = restricted;

    epsilon_a_ = wfn_->epsilon_a();
    Ca_so_ = wfn_->Ca();
    Da_so_ = wfn_->Da();

    if (same_dens_) {
        Db_so_ = Da_so_;
    } else {
        Db_so_ = wfn_->Db();
    }

    if (same_orbs_) {
        epsilon_b_ = epsilon_a_;
        Cb_so_ = Ca_so_;
    } else {
        epsilon_b_ = wfn_->epsilon_b();
        Cb_so_ = wfn_->Cb();
    }
}
void Prop::set_epsilon_a(SharedVector epsilon_a) {
    epsilon_a_ = epsilon_a;
    if (same_orbs_) {
        epsilon_b_ = epsilon_a_;
    }
}
void Prop::set_epsilon_b(SharedVector epsilon_b) {
    if (same_orbs_) throw PSIEXCEPTION("Wavefunction is restricted, setting epsilon_b makes no sense");
    epsilon_b_ = epsilon_b;
}
void Prop::set_Ca(SharedMatrix C) {
    Ca_so_ = C;
    if (same_orbs_) {
        Cb_so_ = Ca_so_;
    }
}
void Prop::set_Cb(SharedMatrix C) {
    if (same_orbs_) throw PSIEXCEPTION("Wavefunction is restricted, setting Cb makes no sense");

    Cb_so_ = C;
}
void Prop::set_Da_ao(SharedMatrix D, int symm) {
    Da_so_ = std::make_shared<Matrix>("Da_so", Ca_so_->rowspi(), Ca_so_->rowspi(), symm);
    std::vector<double> temp(AO2USO_->max_ncol() * AO2USO_->max_nrow());
    double* temp_ptr = temp.data();
    for (int h = 0; h < AO2USO_->nirrep(); ++h) {
        int nao = AO2USO_->rowspi()[0];
        int nsol = AO2USO_->colspi()[h];
        int nsor = AO2USO_->colspi()[h ^ symm];

        if (!nsol || !nsor) continue;

        double** Ulp = AO2USO_->pointer(h);
        double** Urp = AO2USO_->pointer(h ^ symm);
        double** DAOp = D->pointer();
        double** DSOp = Da_so_->pointer(h);
        C_DGEMM('N', 'N', nao, nsor, nao, 1.0, DAOp[0], nao, Urp[0], nsor, 0.0, temp_ptr, nsor);
        C_DGEMM('T', 'N', nsol, nsor, nao, 1.0, Ulp[0], nsol, temp_ptr, nsor, 0.0, DSOp[0], nsor);
    }

    if (same_dens_) {
        Db_so_ = Da_so_;
    }
}
void Prop::set_Db_ao(SharedMatrix D, int symm) {
    if (same_dens_) throw PSIEXCEPTION("Wavefunction is restricted, setting Db makes no sense");

    Db_so_ = std::make_shared<Matrix>("Db_so", Cb_so_->rowspi(), Cb_so_->rowspi(), symm);

    std::vector<double> temp(AO2USO_->max_ncol() * AO2USO_->max_nrow());
    double* temp_ptr = temp.data();
    for (int h = 0; h < AO2USO_->nirrep(); ++h) {
        int nao = AO2USO_->rowspi()[0];
        int nsol = AO2USO_->colspi()[h];
        int nsor = AO2USO_->colspi()[h ^ symm];

        if (!nsol || !nsor) continue;

        double** Ulp = AO2USO_->pointer(h);
        double** Urp = AO2USO_->pointer(h ^ symm);
        double** DAOp = D->pointer();
        double** DSOp = Db_so_->pointer(h);
        C_DGEMM('N', 'N', nao, nsor, nao, 1.0, DAOp[0], nao, Urp[0], nsor, 0.0, temp_ptr, nsor);
        C_DGEMM('T', 'N', nsol, nsor, nao, 1.0, Ulp[0], nsol, temp_ptr, nsor, 0.0, DSOp[0], nsor);
    }
}
void Prop::set_Da_so(SharedMatrix D) {
    Da_so_ = D;
    if (same_dens_) {
        Db_so_ = Da_so_;
    }
}
void Prop::set_Db_so(SharedMatrix D) {
    if (same_dens_) throw PSIEXCEPTION("Wavefunction is restricted, setting Db makes no sense");

    Db_so_ = D;
}
void Prop::set_Da_mo(SharedMatrix D) {
    Da_so_ = std::make_shared<Matrix>("Da_so", Ca_so_->rowspi(), Ca_so_->rowspi(), D->symmetry());

    int symm = D->symmetry();
    int nirrep = D->nirrep();

    std::vector<double> temp(Ca_so_->max_ncol() * Ca_so_->max_nrow());
    double* temp_ptr = temp.data();
    for (int h = 0; h < nirrep; h++) {
        int nmol = Ca_so_->colspi()[h];
        int nmor = Ca_so_->colspi()[h ^ symm];
        int nsol = Ca_so_->rowspi()[h];
        int nsor = Ca_so_->rowspi()[h ^ symm];
        if (!nmol || !nmor || !nsol || !nsor) continue;
        double** Clp = Ca_so_->pointer(h);
        double** Crp = Ca_so_->pointer(h ^ symm);
        double** Dmop = D->pointer(h ^ symm);
        double** Dsop = Da_so_->pointer(h ^ symm);
        C_DGEMM('N', 'T', nmol, nsor, nmor, 1.0, Dmop[0], nmor, Crp[0], nmor, 0.0, temp_ptr, nsor);
        C_DGEMM('N', 'N', nsol, nsor, nmol, 1.0, Clp[0], nmol, temp_ptr, nsor, 0.0, Dsop[0], nsor);
    }

    if (same_dens_) {
        Db_so_ = Da_so_;
    }
}
void Prop::set_Db_mo(SharedMatrix D) {
    if (same_dens_) throw PSIEXCEPTION("Wavefunction is restricted, setting Db makes no sense");

    Db_so_ = std::make_shared<Matrix>("Db_so", Cb_so_->rowspi(), Cb_so_->rowspi(), D->symmetry());

    int symm = D->symmetry();
    int nirrep = D->nirrep();

    std::vector<double> temp(Cb_so_->max_ncol() * Cb_so_->max_nrow());
    double* temp_ptr = temp.data();
    for (int h = 0; h < nirrep; h++) {
        int nmol = Cb_so_->colspi()[h];
        int nmor = Cb_so_->colspi()[h ^ symm];
        int nsol = Cb_so_->rowspi()[h];
        int nsor = Cb_so_->rowspi()[h ^ symm];
        if (!nmol || !nmor || !nsol || !nsor) continue;
        double** Clp = Cb_so_->pointer(h);
        double** Crp = Cb_so_->pointer(h ^ symm);
        double** Dmop = D->pointer(h ^ symm);
        double** Dsop = Db_so_->pointer(h ^ symm);
        C_DGEMM('N', 'T', nmol, nsor, nmor, 1.0, Dmop[0], nmor, Crp[0], nmor, 0.0, temp_ptr, nsor);
        C_DGEMM('N', 'N', nsol, nsor, nmol, 1.0, Clp[0], nmol, temp_ptr, nsor, 0.0, Dsop[0], nsor);
    }
}
void OEProp::add(const std::string& prop) { tasks_.insert(prop); }
void OEProp::add(std::vector<std::string> props) { tasks_.insert(props.begin(), props.end()); }
void OEProp::clear() { tasks_.clear(); }
SharedVector Prop::epsilon_a() { return std::make_shared<Vector>(std::move(epsilon_a_->clone())); }
SharedVector Prop::epsilon_b() { return std::make_shared<Vector>(std::move(epsilon_b_->clone())); }
SharedMatrix Prop::Da_ao() {
    std::vector<double> temp(AO2USO_->max_ncol() * AO2USO_->max_nrow());
    double* temp_ptr = temp.data();
    auto D = std::make_shared<Matrix>("Da (AO basis)", basisset_->nbf(), basisset_->nbf());
    int symm = Da_so_->symmetry();
    for (int h = 0; h < AO2USO_->nirrep(); ++h) {
        int nao = AO2USO_->rowspi()[0];
        int nsol = AO2USO_->colspi()[h];
        int nsor = AO2USO_->colspi()[h ^ symm];
        if (!nsol || !nsor) continue;
        double** Ulp = AO2USO_->pointer(h);
        double** Urp = AO2USO_->pointer(h ^ symm);
        double** DSOp = Da_so_->pointer(h ^ symm);
        double** DAOp = D->pointer();
        C_DGEMM('N', 'T', nsol, nao, nsor, 1.0, DSOp[0], nsor, Urp[0], nsor, 0.0, temp_ptr, nao);
        C_DGEMM('N', 'N', nao, nao, nsol, 1.0, Ulp[0], nsol, temp_ptr, nao, 1.0, DAOp[0], nao);
    }
    return D;
}
SharedMatrix Prop::Db_ao() {
    if (same_dens_) throw PSIEXCEPTION("Wavefunction is restricted, asking for Db makes no sense");

    std::vector<double> temp(AO2USO_->max_ncol() * AO2USO_->max_nrow());
    double* temp_ptr = temp.data();
    auto D = std::make_shared<Matrix>("Db (AO basis)", basisset_->nbf(), basisset_->nbf());
    int symm = Db_so_->symmetry();
    for (int h = 0; h < AO2USO_->nirrep(); ++h) {
        int nao = AO2USO_->rowspi()[0];
        int nsol = AO2USO_->colspi()[h];
        int nsor = AO2USO_->colspi()[h ^ symm];
        if (!nsol || !nsor) continue;
        double** Ulp = AO2USO_->pointer(h);
        double** Urp = AO2USO_->pointer(h ^ symm);
        double** DSOp = Db_so_->pointer(h ^ symm);
        double** DAOp = D->pointer();
        C_DGEMM('N', 'T', nsol, nao, nsor, 1.0, DSOp[0], nsor, Urp[0], nsor, 0.0, temp_ptr, nao);
        C_DGEMM('N', 'N', nao, nao, nsol, 1.0, Ulp[0], nsol, temp_ptr, nao, 1.0, DAOp[0], nao);
    }
    return D;
}
SharedMatrix Prop::Ca_so() { return SharedMatrix(Ca_so_->clone()); }
SharedMatrix Prop::Cb_so() { return SharedMatrix(Cb_so_->clone()); }
SharedMatrix Prop::Ca_ao() { return wfn_->Ca_subset("AO"); }
SharedMatrix Prop::Cb_ao() { return wfn_->Cb_subset("AO"); }
SharedMatrix Prop::Da_so() { return SharedMatrix(Da_so_->clone()); }
SharedMatrix Prop::Db_so() { return SharedMatrix(Db_so_->clone()); }
SharedMatrix Prop::Da_mo() {
    auto D = std::make_shared<Matrix>("Da_mo", Ca_so_->colspi(), Ca_so_->colspi(), Da_so_->symmetry());

    int symm = D->symmetry();
    int nirrep = D->nirrep();

    SharedMatrix S = overlap_so();

    std::vector<double> SC(Ca_so_->max_ncol() * Ca_so_->max_nrow());
    std::vector<double> temp(Ca_so_->max_ncol() * Ca_so_->max_nrow());
    double* SC_ptr = SC.data();
    double* temp_ptr = temp.data();
    for (int h = 0; h < nirrep; h++) {
        int nmol = Ca_so_->colspi()[h];
        int nmor = Ca_so_->colspi()[h ^ symm];
        int nsol = Ca_so_->rowspi()[h];
        int nsor = Ca_so_->rowspi()[h ^ symm];
        if (!nmol || !nmor || !nsol || !nsor) continue;
        double** Slp = S->pointer(h);
        double** Srp = S->pointer(h ^ symm);
        double** Clp = Ca_so_->pointer(h);
        double** Crp = Ca_so_->pointer(h ^ symm);
        double** Dmop = D->pointer(h);
        double** Dsop = Da_so_->pointer(h);

        C_DGEMM('N', 'N', nsor, nmor, nsor, 1.0, Srp[0], nsor, Crp[0], nmor, 0.0, SC_ptr, nmor);
        C_DGEMM('N', 'N', nsol, nmor, nsor, 1.0, Dsop[0], nsor, SC_ptr, nmor, 0.0, temp_ptr, nmor);
        C_DGEMM('N', 'N', nsol, nmol, nsol, 1.0, Slp[0], nsol, Clp[0], nmol, 0.0, SC_ptr, nmol);
        C_DGEMM('T', 'N', nmol, nmor, nsol, 1.0, SC_ptr, nmol, temp_ptr, nmor, 0.0, Dmop[0], nmor);
    }
    return D;
}
SharedMatrix Prop::Db_mo() {
    if (same_dens_) throw PSIEXCEPTION("Wavefunction is restricted, asking for Db makes no sense");

    auto D = std::make_shared<Matrix>("Db_mo", Cb_so_->colspi(), Cb_so_->colspi(), Db_so_->symmetry());

    int symm = D->symmetry();
    int nirrep = D->nirrep();

    SharedMatrix S = overlap_so();

    std::vector<double> SC(Cb_so_->max_ncol() * Cb_so_->max_nrow());
    std::vector<double> temp(Cb_so_->max_ncol() * Cb_so_->max_nrow());

    double* SC_ptr = SC.data();
    double* temp_ptr = temp.data();
    for (int h = 0; h < nirrep; h++) {
        int nmol = Cb_so_->colspi()[h];
        int nmor = Cb_so_->colspi()[h ^ symm];
        int nsol = Cb_so_->rowspi()[h];
        int nsor = Cb_so_->rowspi()[h ^ symm];
        if (!nmol || !nmor || !nsol || !nsor) continue;
        double** Slp = S->pointer(h);
        double** Srp = S->pointer(h ^ symm);
        double** Clp = Cb_so_->pointer(h);
        double** Crp = Cb_so_->pointer(h ^ symm);
        double** Dmop = D->pointer(h);
        double** Dsop = Db_so_->pointer(h);

        C_DGEMM('N', 'N', nsor, nmor, nsor, 1.0, Srp[0], nsor, Crp[0], nmor, 0.0, SC_ptr, nmor);
        C_DGEMM('N', 'N', nsol, nmor, nsor, 1.0, Dsop[0], nsor, SC_ptr, nmor, 0.0, temp_ptr, nmor);
        C_DGEMM('N', 'N', nsol, nmol, nsol, 1.0, Slp[0], nsol, Clp[0], nmol, 0.0, SC_ptr, nmol);
        C_DGEMM('T', 'N', nmol, nmor, nsol, 1.0, SC_ptr, nmol, temp_ptr, nmor, 0.0, Dmop[0], nmor);
    }
    return D;
}
SharedMatrix Prop::Dt_so(bool total) {
    SharedMatrix Da = Da_so();
    SharedMatrix D(Da->clone());
    if (same_dens_) {
        D->set_name((total ? "Dt_so" : "Ds_so"));
        D->scale((total ? 2.0 : 0.0));
    } else {
        D->set_name((total ? "Dt_so" : "Ds_so"));
        SharedMatrix Db = Db_so();
        if (total)
            D->add(Db);
        else
            D->subtract(Db);
    }
    return D;
}
SharedMatrix Prop::Dt_mo(bool total) {
    SharedMatrix Da = Da_mo();
    if (same_dens_) {
        Da->set_name((total ? "Dt_mo" : "Ds_mo"));
        Da->scale((total ? 2.0 : 0.0));
    } else {
        Da->set_name((total ? "Dt_mo" : "Ds_mo"));
        SharedMatrix Db = Db_mo();
        if (total)
            Da->add(Db);
        else
            Da->subtract(Db);
    }
    return Da;
}
std::pair<SharedMatrix, SharedVector> Prop::Na_mo() {
    SharedMatrix D = Da_mo();
    auto C = std::make_shared<Matrix>("Na_mo", D->nirrep(), D->rowspi(), D->rowspi());
    auto O = std::make_shared<Vector>("Alpha Occupation", D->rowspi());

    D->diagonalize(C, O, descending);

    return make_pair(C, O);
}
std::pair<SharedMatrix, SharedVector> Prop::Nb_mo() {
    if (same_dens_) throw PSIEXCEPTION("Wavefunction is restricted, asking for Nb makes no sense");

    SharedMatrix D = Db_mo();
    auto C = std::make_shared<Matrix>("Nb_mo", D->nirrep(), D->rowspi(), D->rowspi());
    auto O = std::make_shared<Vector>("Beta Occupation", D->rowspi());

    D->diagonalize(C, O, descending);

    return make_pair(C, O);
}
std::pair<SharedMatrix, SharedVector> Prop::Na_so() {
    std::pair<SharedMatrix, std::shared_ptr<Vector>> pair = Na_mo();
    SharedMatrix N = pair.first;
    std::shared_ptr<Vector> O = pair.second;

    auto N2 = std::make_shared<Matrix>("Na_so", Ca_so_->nirrep(), Ca_so_->rowspi(), Ca_so_->colspi());

    for (int h = 0; h < N->nirrep(); h++) {
        int nmo = Ca_so_->colspi()[h];
        int nso = Ca_so_->rowspi()[h];

        if (!nmo || !nso) continue;

        double** Np = N->pointer(h);
        double** Cp = Ca_so_->pointer(h);
        double** N2p = N2->pointer(h);

        C_DGEMM('N', 'N', nso, nmo, nmo, 1.0, Cp[0], nmo, Np[0], nmo, 0.0, N2p[0], nmo);
    }
    return make_pair(N2, O);
}
std::pair<SharedMatrix, SharedVector> Prop::Nb_so() {
    if (same_dens_) throw PSIEXCEPTION("Wavefunction is restricted, asking for Nb makes no sense");

    std::pair<SharedMatrix, std::shared_ptr<Vector>> pair = Nb_mo();
    SharedMatrix N = pair.first;
    std::shared_ptr<Vector> O = pair.second;

    auto N2 = std::make_shared<Matrix>("Nb_so", Cb_so_->nirrep(), Cb_so_->rowspi(), Cb_so_->colspi());

    for (int h = 0; h < N->nirrep(); h++) {
        int nmo = Cb_so_->colspi()[h];
        int nso = Cb_so_->rowspi()[h];

        if (!nmo || !nso) continue;

        double** Np = N->pointer(h);
        double** Cp = Cb_so_->pointer(h);
        double** N2p = N2->pointer(h);

        C_DGEMM('N', 'N', nso, nmo, nmo, 1.0, Cp[0], nmo, Np[0], nmo, 0.0, N2p[0], nmo);
    }
    return make_pair(N2, O);
}
std::pair<SharedMatrix, SharedVector> Prop::Na_ao() {
    std::pair<SharedMatrix, std::shared_ptr<Vector>> pair = Na_so();
    SharedMatrix N = pair.first;
    std::shared_ptr<Vector> O = pair.second;

    auto N2 = std::make_shared<Matrix>("Na_ao", Ca_so_->nrow(), Ca_so_->ncol());
    auto N3 = std::make_shared<Matrix>("Na_ao", Ca_so_->nrow(), Ca_so_->ncol());
    auto O2 = std::make_shared<Vector>("Alpha Occupation", Ca_so_->ncol());

    int offset = 0;
    std::vector<std::pair<double, int>> index;
    for (int h = 0; h < Ca_so_->nirrep(); h++) {
        int ncol = Ca_so_->ncol();
        int nmo = Ca_so_->colspi()[h];
        int nso = AO2USO_->colspi()[h];
        int nao = AO2USO_->rowspi()[h];

        if (!nmo || !nso || !nao) continue;

        for (int i = 0; i < nmo; i++) {
            index.push_back(std::make_pair(O->get(h, i), i + offset));
        }

        double** Np = N->pointer(h);
        double** Up = AO2USO_->pointer(h);
        double** N2p = N2->pointer(h);

        C_DGEMM('N', 'N', nao, nmo, nso, 1.0, Up[0], nso, Np[0], nmo, 0.0, &N2p[0][offset], ncol);

        offset += nmo;
    }

    std::sort(index.begin(), index.end(), std::greater<std::pair<double, int>>());

    int nmo = N2->colspi()[0];
    int nao = N2->rowspi()[0];

    for (int i = 0; i < nmo; i++) {
        double occ = index[i].first;
        int ind = index[i].second;
        O2->set(i, occ);

        C_DCOPY(nao, &(N2->pointer()[0][ind]), nmo, &(N3->pointer()[0][i]), nmo);
    }

    return make_pair(N3, O2);
}
std::pair<SharedMatrix, SharedVector> Prop::Nb_ao() {
    if (same_dens_) throw PSIEXCEPTION("Wavefunction is restricted, asking for Nb makes no sense");

    std::pair<SharedMatrix, std::shared_ptr<Vector>> pair = Nb_so();
    SharedMatrix N = pair.first;
    std::shared_ptr<Vector> O = pair.second;

    auto N2 = std::make_shared<Matrix>("Nb_ao", Cb_so_->nrow(), Cb_so_->ncol());
    auto N3 = std::make_shared<Matrix>("Nb_ao", Cb_so_->nrow(), Cb_so_->ncol());
    auto O2 = std::make_shared<Vector>("Beta Occupation", Cb_so_->ncol());

    int offset = 0;
    std::vector<std::pair<double, int>> index;
    for (int h = 0; h < Cb_so_->nirrep(); h++) {
        int ncol = Cb_so_->ncol();
        int nmo = Cb_so_->colspi()[h];
        int nso = AO2USO_->colspi()[h];
        int nao = AO2USO_->rowspi()[h];

        if (!nmo || !nso || !nao) continue;

        for (int i = 0; i < nmo; i++) {
            index.push_back(std::make_pair(O->get(h, i), i + offset));
        }

        double** Np = N->pointer(h);
        double** Up = AO2USO_->pointer(h);
        double** N2p = N2->pointer(h);

        C_DGEMM('N', 'N', nao, nmo, nso, 1.0, Up[0], nso, Np[0], nmo, 0.0, &N2p[0][offset], ncol);

        offset += nmo;
    }

    std::sort(index.begin(), index.end(), std::greater<std::pair<double, int>>());

    int nmo = N2->colspi()[0];
    int nao = N2->rowspi()[0];

    for (int i = 0; i < nmo; i++) {
        double occ = index[i].first;
        int ind = index[i].second;
        O2->set(i, occ);

        C_DCOPY(nao, &(N2->pointer()[0][ind]), nmo, &(N3->pointer()[0][i]), nmo);
    }

    return make_pair(N3, O2);
}
std::pair<SharedMatrix, SharedVector> Prop::Nt_mo() {
    SharedMatrix D = Dt_mo();
    auto C = std::make_shared<Matrix>("Nt_mo", D->nirrep(), D->rowspi(), D->rowspi());
    auto O = std::make_shared<Vector>("Total Occupation", D->rowspi());

    D->diagonalize(C, O, descending);

    return make_pair(C, O);
}
std::pair<SharedMatrix, SharedVector> Prop::Nt_so() {
    std::pair<SharedMatrix, std::shared_ptr<Vector>> pair = Nt_mo();
    SharedMatrix N = pair.first;
    std::shared_ptr<Vector> O = pair.second;

    auto N2 = std::make_shared<Matrix>("Nt_so", Cb_so_->nirrep(), Cb_so_->rowspi(), Cb_so_->colspi());

    for (int h = 0; h < N->nirrep(); h++) {
        int nmo = Cb_so_->colspi()[h];
        int nso = Cb_so_->rowspi()[h];

        if (!nmo || !nso) continue;

        double** Np = N->pointer(h);
        double** Cp = Cb_so_->pointer(h);
        double** N2p = N2->pointer(h);

        C_DGEMM('N', 'N', nso, nmo, nmo, 1.0, Cp[0], nmo, Np[0], nmo, 0.0, N2p[0], nmo);
    }
    return make_pair(N2, O);
}
std::pair<SharedMatrix, SharedVector> Prop::Nt_ao() {
    std::pair<SharedMatrix, std::shared_ptr<Vector>> pair = Nt_so();
    SharedMatrix N = pair.first;
    std::shared_ptr<Vector> O = pair.second;

    auto N2 = std::make_shared<Matrix>("Nt_ao", Cb_so_->nrow(), Cb_so_->ncol());
    auto N3 = std::make_shared<Matrix>("Nt_ao", Cb_so_->nrow(), Cb_so_->ncol());
    auto O2 = std::make_shared<Vector>("Total Occupation", Cb_so_->ncol());

    int offset = 0;
    std::vector<std::pair<double, int>> index;
    for (int h = 0; h < Cb_so_->nirrep(); h++) {
        int ncol = Cb_so_->ncol();
        int nmo = Cb_so_->colspi()[h];
        int nso = AO2USO_->colspi()[h];
        int nao = AO2USO_->rowspi()[h];

        if (!nmo || !nso || !nao) continue;

        for (int i = 0; i < nmo; i++) {
            index.push_back(std::make_pair(O->get(h, i), i + offset));
        }

        double** Np = N->pointer(h);
        double** Up = AO2USO_->pointer(h);
        double** N2p = N2->pointer(h);

        C_DGEMM('N', 'N', nao, nmo, nso, 1.0, Up[0], nso, Np[0], nmo, 0.0, &N2p[0][offset], ncol);

        offset += nmo;
    }

    std::sort(index.begin(), index.end(), std::greater<std::pair<double, int>>());

    int nmo = N2->colspi()[0];
    int nao = N2->rowspi()[0];

    for (int i = 0; i < nmo; i++) {
        double occ = index[i].first;
        int ind = index[i].second;
        O2->set(i, occ);

        C_DCOPY(nao, &(N2->pointer()[0][ind]), nmo, &(N3->pointer()[0][i]), nmo);
    }

    return make_pair(N3, O2);
}
SharedMatrix Prop::overlap_so() {
    auto S = std::make_shared<Matrix>("S", Ca_so_->rowspi(), Ca_so_->rowspi());
    std::shared_ptr<OneBodySOInt> Sint(integral_->so_overlap());
    Sint->compute(S);
    return S;
}

std::string Prop::Da_name() const { return Da_so_->name(); }
std::string Prop::Db_name() const { return Db_so_->name(); }

Vector3 OEProp::compute_center(const double* property) const {
    std::shared_ptr<Molecule> mol = wfn_->molecule();
    int natoms = mol->natom();
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double sum = 0.0;
    for (int atom = 0; atom < natoms; ++atom) {
        Vector3 xyz = mol->xyz(atom);
        double prop = property[atom];
        x += xyz[0] * prop;
        y += xyz[1] * prop;
        z += xyz[2] * prop;
        sum += prop;
    }
    x /= sum;
    y /= sum;
    z /= sum;
    return Vector3(x, y, z);
}

OEProp::OEProp(std::shared_ptr<Wavefunction> wfn)
    : wfn_(wfn), mpc_(wfn, get_origin_from_environment()), pac_(wfn), epc_(wfn) {
    if (wfn_.get() == nullptr) throw PSIEXCEPTION("Prop: Wavefunction is null");
    common_init();
}
OEProp::~OEProp() {}

void OEProp::common_init() {
    Options& options = Process::environment.options;
    print_ = options.get_int("PRINT");

    // Determine number of NOONs to print; default is 3
    if (options.get_str("PRINT_NOONS") == "ALL")
        max_noon_ = wfn_->nmo();
    else
        max_noon_ = to_integer(options.get_str("PRINT_NOONS"));
}

MultipolePropCalc::MultipolePropCalc(std::shared_ptr<Wavefunction> wfn, Vector3 const& origin)
    : Prop(wfn), origin_(origin) {
    std::shared_ptr<Molecule> mol = basisset_->molecule();

    /*
     * Check the symmetry of the origin; if it's off-axis we can't use symmetry for multipoles anymore
     */
    auto ct = mol->point_group()->char_table();
    int nirrep = ct.nirrep();

    origin_preserves_symmetry_ = true;
    for (int irrep = 1; irrep < nirrep; ++irrep) {
        IrreducibleRepresentation gamma = ct.gamma(irrep);
        double t[] = {0.0, 0.0, 0.0};

        // Apply the projection
        for (int G = 0; G < nirrep; ++G) {
            SymmetryOperation so = ct.symm_operation(G);
            ShellRotation rr(1, so, integral_.get(), false);

            // rr(xyz, xyz) tells us how the orbitals transform in this
            // symmetry operation, then we multiply by the character in
            // the irrep
            for (int xyz = 0; xyz < 3; ++xyz) t[xyz] += origin_[xyz] * rr(xyz, xyz) * gamma.character(G) / nirrep;
        }

        for (int xyz = 0; xyz < 3; ++xyz) {
            if (std::fabs(t[xyz]) > 1.0E-8) {
                outfile->Printf("The origin chosen breaks symmetry; multipoles will be computed without symmetry.\n");
                origin_preserves_symmetry_ = false;
            }
        }
    }
}

Vector3 OEProp::get_origin_from_environment() const {
    // This function gets called early in the constructor of OEProp.
    // Take care not to use or initialize members, which are only initialized later.
    // The only member used here is basisset, which is initialized in the base class.
    // See if the user specified the origin
    Options& options = Process::environment.options;
    Vector3 origin(0.0, 0.0, 0.0);

    std::shared_ptr<Molecule> mol = wfn_->molecule();
    int natoms = mol->natom();
    if (options["PROPERTIES_ORIGIN"].has_changed()) {
        int size = options["PROPERTIES_ORIGIN"].size();

        if (size == 1) {
            std::vector<double> property(natoms);
            std::string str = options["PROPERTIES_ORIGIN"][0].to_string();
            if (str == "COM") {
                for (int atom = 0; atom < natoms; ++atom) property[atom] = mol->mass(atom);
            } else if (str == "NUCLEAR_CHARGE") {
                for (int atom = 0; atom < natoms; ++atom) property[atom] = mol->charge(atom);
            } else {
                throw PSIEXCEPTION("Invalid specification of PROPERTIES_ORIGIN.  Please consult the manual.");
            }
            origin = compute_center(property.data());
        } else if (size == 3) {
            double x = options["PROPERTIES_ORIGIN"][0].to_double();
            double y = options["PROPERTIES_ORIGIN"][1].to_double();
            double z = options["PROPERTIES_ORIGIN"][2].to_double();
            bool convert = mol->units() == Molecule::Angstrom;
            if (convert) {
                x /= pc_bohr2angstroms;
                y /= pc_bohr2angstroms;
                z /= pc_bohr2angstroms;
            }
            origin = Vector3(x, y, z);
        } else {
            throw PSIEXCEPTION("Invalid specification of PROPERTIES_ORIGIN.  Please consult the manual.");
        }
    }
    outfile->Printf("\n\nProperties will be evaluated at %10.6f, %10.6f, %10.6f [a0]\n", origin[0], origin[1],
                    origin[2]);

    return origin;
}

void OEProp::print_header() {
    outfile->Printf("\n OEPROP: One-electron properties/analyses.\n");
    outfile->Printf("  by Rob Parrish and Justin Turney.\n");
    outfile->Printf("  built on LIBMINTS.\n\n");
}

template <class T>
bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&)) {
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
}

void OEProp::compute() {
    if (title_ == "") {
        outfile->Printf("OEProp: No title given, name of density matrix used for the following properties is '%s'\n",
                        mpc_.Da_name().c_str());
    } else {
        outfile->Printf("\nProperties computed using the %s density matrix\n\n", title_.c_str());
    }

    // Search for multipole strings, which are handled separately
    std::set<std::string>::const_iterator iter = tasks_.begin();
    std::regex mpoles("^MULTIPOLE(?:S)?\\s*\\((\\d+)\\)$");
    std::smatch matches;
    for (; iter != tasks_.end(); ++iter) {
        std::string str = *iter;
        if (std::regex_match(str, matches, mpoles)) {
            int order;
            if (!from_string<int>(order, matches[1], std::dec))
                throw PSIEXCEPTION("Problem determining multipole order!  Specify, e.g., MULTIPOLE(5)");
            compute_multipoles(order);
        }
    }
    // print_header();  // Not by default, happens too often -CDS
    if (tasks_.count("ESP_AT_NUCLEI")) compute_esp_at_nuclei();
    if (tasks_.count("QUADRUPOLE")) compute_multipoles(2, false);
    else if (tasks_.count("DIPOLE")) compute_multipoles(1, false);
    if (tasks_.count("TRANSITION_QUADRUPOLE")) compute_multipoles(2, true);
    else if (tasks_.count("TRANSITION_DIPOLE")) compute_multipoles(1, true);
    if (tasks_.count("MO_EXTENTS")) compute_mo_extents();
    if (tasks_.count("MULLIKEN_CHARGES")) compute_mulliken_charges();
    if (tasks_.count("LOWDIN_CHARGES")) compute_lowdin_charges();
    if (tasks_.count("MBIS_VOLUME_RATIOS")) compute_mbis_multipoles(true);
    else if (tasks_.count("MBIS_CHARGES")) compute_mbis_multipoles(false);
    if (tasks_.count("MAYER_INDICES")) compute_mayer_indices();
    if (tasks_.count("WIBERG_LOWDIN_INDICES")) compute_wiberg_lowdin_indices();
    if (tasks_.count("NO_OCCUPATIONS")) compute_no_occupations();
    if (tasks_.count("GRID_FIELD")) compute_field_over_grid();
    if (tasks_.count("GRID_ESP")) compute_esp_over_grid();
}

void OEProp::compute_multipoles(int order, bool transition) {
    auto mpoles = mpc_.compute_multipoles(order, transition, true, print_ > 4);

    for (int l = 1; l <= order; ++l) {
        int component = 0;
        int ncomponents = (l + 1) * (l + 2) / 2;
        auto multipole_array = std::make_shared<Matrix>(1, ncomponents);


        std::string mtdname;
        if (l == 1) {
            mtdname = "DIPOLE";
        } else if (l == 2) {
            mtdname = "QUADRUPOLE";
        } else if (l == 3) {
            mtdname = "OCTUPOLE";
        } else if (l == 4) {
            mtdname = "HEXADECAPOLE";
        } else {
            int n = (1 << l);
            mtdname = std::to_string(n) + "-POLE";
        }

        for (auto it = mpoles->begin(); it != mpoles->end(); ++it) {
            std::string name;
            double total_mpole = 0.0;
            int order_mpole;
            // unpack the multipole, which is: name, nuc, elec, total, order, ignore nuc and elec:
            std::tie(name, std::ignore, std::ignore, total_mpole, order_mpole) = *it;
            /*- Process::environment.globals["mtd DIPOLE"] -*/
            /*- Process::environment.globals["mtd QUADRUPOLE"] -*/
            /*- Process::environment.globals["mtd OCTUPOLE"] -*/
            /*- Process::environment.globals["mtd HEXADECAPOLE"] -*/
            /*- Process::environment.globals["mtd 32-POLE"] -*/
            /*- Process::environment.globals["mtd 64-POLE"] -*/
            /*- Process::environment.globals["mtd 128-POLE"] -*/

            if (order_mpole == l) {
                multipole_array->set(0, component, total_mpole);
                ++component;
            }
        }

        std::unordered_set<std::string> names;
        if (names_.empty()) {
            names = {"{}"};
        } else {
            names = names_;
        }

        for (auto name: names) {
            // TODO: Use fmt strings when Psi uses C++20.
            name.replace(name.find("{}"), sizeof("{}") - 1, mtdname);
            Process::environment.arrays[name] = multipole_array;
            wfn_->set_array_variable(name, multipole_array);
        }
    }
}

MultipolePropCalc::MultipoleOutputType MultipolePropCalc::compute_multipoles(int order, bool transition,
                                                                             bool print_output, bool verbose) {
    auto mot = std::make_shared<MultipolePropCalc::MultipoleOutputTypeBase>();
    auto mol = basisset_->molecule();

    SharedMatrix Da;
    SharedMatrix Db;

    std::vector<SharedMatrix> mp_ints;
    MultipoleSymmetry mpsymm(order, mol, integral_, factory_);

    if (origin_preserves_symmetry_) {
        // We can use symmetry here
        std::shared_ptr<OneBodySOInt> sompOBI(integral_->so_multipoles(order));
        sompOBI->ob()->set_origin(origin_);
        mp_ints = mpsymm.create_matrices("", false);
        sompOBI->compute(mp_ints);

        if (same_dens_) {
            Da = Da_so_;
            Db = Da;
        } else {
            Da = Da_so_;
            Db = Db_so_;
        }
    } else {
        // No symmetry
        mp_ints = mpsymm.create_matrices("", true);
        std::shared_ptr<OneBodyAOInt> aompOBI(integral_->ao_multipoles(order));
        aompOBI->set_origin(origin_);
        aompOBI->compute(mp_ints);

        if (same_dens_) {
            Da = wfn_->matrix_subset_helper(Da_so_, Ca_so_, "AO", "D");
            Db = Da;
        } else {
            Da = wfn_->matrix_subset_helper(Da_so_, Ca_so_, "AO", "D alpha");
            Db = wfn_->matrix_subset_helper(Db_so_, Cb_so_, "AO", "D beta");
        }
    }

    if (verbose) {
        std::vector<SharedMatrix>::iterator iter;
        for (iter = mp_ints.begin(); iter != mp_ints.end(); ++iter) iter->get()->print();
    }

    auto nuclear_contributions = MultipoleInt::nuclear_contribution(mol, order, origin_);

    // => Print and populate the output tuple <=
    if (print_output) {
        outfile->Printf("\n%s Multipole Moments:\n", transition ? "Transition" : "");
        outfile->Printf("\n ------------------------------------------------------------------------------------\n");
        outfile->Printf("     Multipole            Electronic (a.u.)      Nuclear  (a.u.)        Total (a.u.)\n");
        outfile->Printf(" ------------------------------------------------------------------------------------\n\n");
    }
    double convfac = pc_dipmom_au2debye;
    int address = 0; // The index of the desired component in the mp_ints object, which covers _all_ orders.
    for (int l = 1; l <= order; ++l) {
        int ncomponents = (l + 1) * (l + 2) / 2;
        std::stringstream ss, tt;
        if (l > 1) tt << "^" << l;
        if (l > 1) ss << " ang";
        if (l > 2) ss << "^" << l - 1;
        std::string exp = ss.str();
        if (print_output) {
            outfile->Printf(" L = %d.  Multiply by %.10f to convert [e a0%s] to [Debye%s]\n", l, convfac, tt.str().c_str(), exp.c_str());
        }
        for (int component = 0; component < ncomponents; ++component) {
            auto mpmat = mp_ints[address];
            auto name = mpmat->name();
            double nuc = transition ? 0.0 : nuclear_contributions->get(address);
            double elec = Da->vector_dot(mpmat) + Db->vector_dot(mpmat);
            double tot = nuc + elec;
            if (print_output) {
                outfile->Printf(" %-20s: %18.7f   %18.7f   %18.7f\n", name.c_str(), elec, nuc, tot);
            }
            auto upper_name = to_upper_copy(name);
            mot->emplace_back(upper_name, nuc, elec, tot, l);
            ++address;
        }
        if (print_output) {
            // For historical reasons, we print dipole magnitudes now.
            // Printing magnitudes of higher-order multipole tensors would be nice, but we don't have the C-side machinery
            // to convert our vector of unique multipole components into a container of all components.
            if (l == 1) {
                double tot = sqrt(pow(std::get<3>((*mot)[0]), 2) + pow(std::get<3>((*mot)[1]), 2) + pow(std::get<3>((*mot)[2]), 2));
                outfile->Printf(" %-20s: %18s   %18s   %18.7f\n", "Magnitude", "", "", tot);
            }
            // For historical reasons, we print traceless quadrupole magnitudes now.
            // This can be generalized to any l < 1, but that's hard.
            // Computing integrals is also hard, which is why we don't recompute traceless integrals.
            if (l == 2) {
                double tr_n = (std::get<1>((*mot)[3]) + std::get<1>((*mot)[6]) + std::get<1>((*mot)[8])) / 3.0;
                double tr_e = (std::get<2>((*mot)[3]) + std::get<2>((*mot)[6]) + std::get<2>((*mot)[8])) / 3.0;
                double tr_t = (std::get<3>((*mot)[3]) + std::get<3>((*mot)[6]) + std::get<3>((*mot)[8])) / 3.0;
                outfile->Printf(" %-20s: %18.7f   %18.7f   %18.7f\n", "Traceless XX", std::get<2>((*mot)[3]) - tr_e, std::get<1>((*mot)[3]) - tr_n, std::get<3>((*mot)[3]) - tr_t);
                outfile->Printf(" %-20s: %18.7f   %18.7f   %18.7f\n", "Traceless YY", std::get<2>((*mot)[6]) - tr_e, std::get<1>((*mot)[6]) - tr_n, std::get<3>((*mot)[6]) - tr_t);
                outfile->Printf(" %-20s: %18.7f   %18.7f   %18.7f\n", "Traceless ZZ", std::get<2>((*mot)[8]) - tr_e, std::get<1>((*mot)[8]) - tr_n, std::get<3>((*mot)[8]) - tr_t);
            }
            outfile->Printf("\n");
        }
        convfac *= pc_bohr2angstroms;
    }
    if (print_output) {
        outfile->Printf(" ------------------------------------------------------------------------------------\n");
    }

    return mot;
}

/**
 * @brief The GridIterator class:  A class to iterate over a user-provided grid for
 *                                 computing electrostatic properties.
 */
class PSI_API GridIterator {
    std::ifstream gridfile_;
    Vector3 gridpoints_;

   public:
    const Vector3& gridpoints() const { return gridpoints_; }
    GridIterator(const std::string& filename) {
        gridfile_.open(filename.c_str());
        if (!gridfile_) throw PSIEXCEPTION("Unable to open the grid.dat file.");
    }
    void first() {
        if (!gridfile_) throw PSIEXCEPTION("File not initialized in Griditer::first.");
        gridfile_.clear();
        gridfile_.seekg(0, std::ios::beg);
        next();
    }
    void next() {
        if (!gridfile_) throw PSIEXCEPTION("Griditer::next called before file stream was initialized.");
        if (!(gridfile_ >> gridpoints_[0])) {
            if (gridfile_.eof()) {
                // Hitting the end of file is OK in this situation.
                return;
            } else {
                throw PSIEXCEPTION("Problem reading x gridpoint from the grid file.");
            }
        }
        if (!(gridfile_ >> gridpoints_[1])) throw PSIEXCEPTION("Problem reading y gridpoint from the grid file.");
        if (!(gridfile_ >> gridpoints_[2])) throw PSIEXCEPTION("Problem reading z gridpoint from the grid file.");
    }
    bool last() const { return gridfile_.eof(); }
    ~GridIterator() { gridfile_.close(); }
};

ESPPropCalc::ESPPropCalc(std::shared_ptr<Wavefunction> wfn) : Prop(wfn) {}

ESPPropCalc::~ESPPropCalc() {}

void OEProp::compute_esp_over_grid() { epc_.compute_esp_over_grid(true); }

void ESPPropCalc::compute_esp_over_grid(bool print_output) {
    auto mol = basisset_->molecule();

    std::shared_ptr<ElectrostaticInt> epot(dynamic_cast<ElectrostaticInt*>(integral_->electrostatic().release()));

    if (print_output) {
        outfile->Printf("\n Electrostatic potential to be computed on the grid and written to grid_esp.dat\n");
    }
    SharedMatrix Dtot = wfn_->matrix_subset_helper(Da_so_, Ca_so_, "AO", "D");
    if (same_dens_) {
        Dtot->scale(2.0);
    } else {
        Dtot->add(wfn_->matrix_subset_helper(Db_so_, Cb_so_, "AO", "D beta"));
    }

    int nbf = basisset_->nbf();
    auto ints = std::make_shared<Matrix>("Ex integrals", nbf, nbf);

    Vvals_.clear();
    FILE* gridout = fopen("grid_esp.dat", "w");
    if (!gridout) throw PSIEXCEPTION("Unable to write to grid_esp.dat");
    GridIterator griditer("grid.dat");
    for (griditer.first(); !griditer.last(); griditer.next()) {
        Vector3 origin(griditer.gridpoints());
        if (mol->units() == Molecule::Angstrom) origin /= pc_bohr2angstroms;
        ints->zero();
        epot->compute(ints, origin);
        double Velec = Dtot->vector_dot(ints);
        double Vnuc = 0.0;
        int natom = mol->natom();
        for (int i = 0; i < natom; i++) {
            Vector3 dR = origin - mol->xyz(i);
            double r = dR.norm();
            if (r > 1.0E-8) Vnuc += mol->Z(i) / r;
        }
        Vvals_.push_back(Velec + Vnuc);
        fprintf(gridout, "%16.10f\n", Velec + Vnuc);
    }
    fclose(gridout);
}

SharedVector ESPPropCalc::compute_esp_over_grid_in_memory(SharedMatrix input_grid) const {
    // We only want a plain matrix to work with here:
    if (input_grid->nirrep() != 1) {
        throw PSIEXCEPTION("ESPPropCalc only allows \"plain\" input matrices with, i.e. nirrep == 1.");
    }
    if (input_grid->coldim() != 3) {
        throw PSIEXCEPTION("ESPPropCalc only allows \"plain\" input matrices with a dimension of N (rows) x 3 (cols)");
    }

    int number_of_grid_points = input_grid->rowdim();
    SharedVector output = std::make_shared<Vector>(number_of_grid_points);

    std::shared_ptr<Molecule> mol = basisset_->molecule();

    SharedMatrix Dtot = wfn_->matrix_subset_helper(Da_so_, Ca_so_, "AO", "D");
    if (same_dens_) {
        Dtot->scale(2.0);
    } else {
        Dtot->add(wfn_->matrix_subset_helper(Db_so_, Cb_so_, "AO", "D beta"));
    }

    int const nbf = basisset_->nbf();

    bool convert = mol->units() == Molecule::Angstrom;

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    std::vector<std::shared_ptr<Matrix>> VtempT;
    std::vector<std::shared_ptr<ElectrostaticInt>> VintT;

    for (int thread = 0; thread < nthreads; thread++) {
        VtempT.push_back(std::make_shared<Matrix>("ints", nbf, nbf));
        VintT.push_back(std::shared_ptr<ElectrostaticInt>(static_cast<ElectrostaticInt*>(integral_->electrostatic().release())));
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < number_of_grid_points; ++i) {
        Vector3 origin(input_grid->get(i, 0), input_grid->get(i, 1), input_grid->get(i, 2));
        if (convert) origin /= pc_bohr2angstroms;
         
        // Thread info
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        // => Electronic part <= //
        VtempT[thread]->zero();
        VintT[thread]->compute(VtempT[thread],origin);
        
        double Velec = Dtot->vector_dot(VtempT[thread]);

        // => Nuclear part <= //
        double Vnuc = 0.0;
        int natom = mol->natom();
        for (int iat = 0; iat < natom; iat++) {
            Vector3 dR = origin - mol->xyz(iat);
            double r = dR.norm();
            if (r > 1.0E-8) Vnuc += mol->Z(iat) / r;
        }

        (*output)[i] = Velec + Vnuc;
    }
    return output;
}

void OEProp::compute_field_over_grid() { epc_.compute_field_over_grid(true); }

void ESPPropCalc::compute_field_over_grid(bool print_output) {
    std::shared_ptr<Molecule> mol = basisset_->molecule();

    std::shared_ptr<ElectrostaticInt> epot(dynamic_cast<ElectrostaticInt*>(integral_->electrostatic().release()));

    if (print_output) {
        outfile->Printf("\n Field computed on the grid and written to grid_field.dat\n");
    }

    SharedMatrix Dtot = wfn_->matrix_subset_helper(Da_so_, Ca_so_, "AO", "D");
    if (same_dens_) {
        Dtot->scale(2.0);
    } else {
        Dtot->add(wfn_->matrix_subset_helper(Db_so_, Cb_so_, "AO", "D beta"));
    }

    std::shared_ptr<ElectricFieldInt> field_ints(dynamic_cast<ElectricFieldInt*>(wfn_->integral()->electric_field().release()));

    int nbf = basisset_->nbf();
    std::vector<SharedMatrix> intmats;
    intmats.push_back(std::make_shared<Matrix>("Ex integrals", nbf, nbf));
    intmats.push_back(std::make_shared<Matrix>("Ey integrals", nbf, nbf));
    intmats.push_back(std::make_shared<Matrix>("Ez integrals", nbf, nbf));

    Exvals_.clear();
    Eyvals_.clear();
    Ezvals_.clear();

    FILE* gridout = fopen("grid_field.dat", "w");
    if (!gridout) throw PSIEXCEPTION("Unable to write to grid_field.dat");
    GridIterator griditer("grid.dat");
    for (griditer.first(); !griditer.last(); griditer.next()) {
        Vector3 origin(griditer.gridpoints());
        if (mol->units() == Molecule::Angstrom) origin /= pc_bohr2angstroms;
        field_ints->set_origin(origin);
        for (int m = 0; m < 3; ++m) intmats[m]->zero();
        field_ints->compute(intmats);
        double Ex = Dtot->vector_dot(intmats[0]);
        double Ey = Dtot->vector_dot(intmats[1]);
        double Ez = Dtot->vector_dot(intmats[2]);
        Vector3 nuc = field_ints->nuclear_contribution(origin, mol);
        Exvals_.push_back(Ex + nuc[0]);
        Eyvals_.push_back(Ey + nuc[1]);
        Ezvals_.push_back(Ez + nuc[2]);
        fprintf(gridout, "%16.10f %16.10f %16.10f\n", Ex + nuc[0], Ey + nuc[1], Ez + nuc[2]);
    }
    fclose(gridout);
}

SharedMatrix ESPPropCalc::compute_field_over_grid_in_memory(SharedMatrix input_grid) const {
    // We only want a plain matrix to work with here:
    if (input_grid->nirrep() != 1) {
        throw PSIEXCEPTION("ESPPropCalc only allows \"plain\" input matrices with, i.e. nirrep == 1.");
    }
    if (input_grid->coldim() != 3) {
        throw PSIEXCEPTION("ESPPropCalc only allows \"plain\" input matrices with a dimension of N (rows) x 3 (cols)");
    }

    std::shared_ptr<Molecule> mol = basisset_->molecule();

    SharedMatrix Dtot = wfn_->Da_subset("AO");
    if (same_dens_) {
        Dtot->scale(2.0);
    } else {
        Dtot->add(wfn_->Db_subset("AO"));
    }

    std::shared_ptr<ElectricFieldInt> field_ints(dynamic_cast<ElectricFieldInt*>(wfn_->integral()->electric_field().release()));

    // Scale the coordinates if needed
    auto coords = input_grid;

    if (mol->units() == Molecule::Angstrom) {
        coords = input_grid->clone();
        coords->scale(1.0 / pc_bohr2angstroms);
    }

    // Compute the electric field at all grid points.
    int number_of_grid_points = coords->rowdim();

    SharedMatrix efield = std::make_shared<Matrix>("efield", number_of_grid_points, 3);

    auto field_functor = ContractOverDensityFieldFunctor(efield, Dtot);
    field_ints->compute_with_functor(field_functor, coords);

    // Add the nuclear contribution.
    for (int i = 0; i < number_of_grid_points; ++i) {
        Vector3 origin(coords->get(i, 0), coords->get(i, 1), coords->get(i, 2));
        Vector3 nuc = field_ints->nuclear_contribution(origin, mol);
        efield->set(i, 0, efield->get(i, 0) + nuc[0]);
        efield->set(i, 1, efield->get(i, 1) + nuc[1]);
        efield->set(i, 2, efield->get(i, 2) + nuc[2]);
    }
    return efield;
}

void OEProp::compute_esp_at_nuclei() {
    std::shared_ptr<std::vector<double>> nesps = epc_.compute_esp_at_nuclei(true, print_ > 2);
    for (int atom1 = 0; atom1 < nesps->size(); ++atom1) {
        std::stringstream s;
        s << "ESP AT CENTER " << atom1 + 1;
        /*- Process::environment.globals["ESP AT CENTER n"] -*/
        Process::environment.globals[s.str()] = (*nesps)[atom1];
        wfn_->set_scalar_variable(s.str(), (*nesps)[atom1]);
    }
    wfn_->set_esp_at_nuclei(nesps);
}

std::shared_ptr<std::vector<double>> ESPPropCalc::compute_esp_at_nuclei(bool print_output, bool verbose) {
    std::shared_ptr<Molecule> mol = basisset_->molecule();

    auto nesps = std::make_shared<std::vector<double>>(mol->natom());
    std::shared_ptr<ElectrostaticInt> epot(dynamic_cast<ElectrostaticInt*>(integral_->electrostatic().release()));

    int nbf = basisset_->nbf();
    int natoms = mol->natom();

    SharedMatrix Dtot = wfn_->matrix_subset_helper(Da_so_, Ca_so_, "AO", "D");
    if (same_dens_) {
        Dtot->scale(2.0);
    } else {
        Dtot->add(wfn_->matrix_subset_helper(Db_so_, Cb_so_, "AO", "D beta"));
    }

    Matrix dist = mol->distance_matrix();
    if (print_output) {
        outfile->Printf("\n Electrostatic potentials at the nuclear coordinates:\n");
        outfile->Printf(" ---------------------------------------------\n");
        outfile->Printf("   Center     Electrostatic Potential (a.u.)\n");
        outfile->Printf(" ---------------------------------------------\n");
    }
    for (int atom1 = 0; atom1 < natoms; ++atom1) {
        std::stringstream s;
        s << "ESP AT CENTER " << atom1 + 1;
        auto ints = std::make_shared<Matrix>(s.str(), nbf, nbf);
        epot->compute(ints, mol->xyz(atom1));
        if (verbose) {
            ints->print();
        }
        double elec = Dtot->vector_dot(ints);
        double nuc = 0.0;
        for (int atom2 = 0; atom2 < natoms; ++atom2) {
            if (atom1 == atom2) continue;
            nuc += mol->Z(atom2) / dist[0][atom1][atom2];
        }
        if (print_output) {
            outfile->Printf("  %3d %2s           %16.12f\n", atom1 + 1, mol->label(atom1).c_str(), nuc + elec);
        }
        (*nesps)[atom1] = nuc + elec;
    }
    if (print_output) {
        outfile->Printf(" ---------------------------------------------\n");
    }
    return nesps;
}

void OEProp::compute_mo_extents() {
    std::vector<SharedVector> mo_es = mpc_.compute_mo_extents(true);
    wfn_->set_mo_extents(mo_es);
}

std::vector<SharedVector> MultipolePropCalc::compute_mo_extents(bool print_output) {
    SharedMatrix Ca = Ca_ao();

    std::vector<SharedVector> mo_es;
    mo_es.push_back(std::make_shared<Vector>("<x^2>", basisset_->nbf()));
    mo_es.push_back(std::make_shared<Vector>("<y^2>", basisset_->nbf()));
    mo_es.push_back(std::make_shared<Vector>("<z^2>", basisset_->nbf()));
    mo_es.push_back(std::make_shared<Vector>("<r^2>", basisset_->nbf()));

    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> ao_Qpole;
    ao_Qpole.push_back(std::make_shared<Matrix>("Quadrupole XX", basisset_->nbf(), basisset_->nbf()));
    ao_Qpole.push_back(std::make_shared<Matrix>("Quadrupole XY", basisset_->nbf(), basisset_->nbf()));
    ao_Qpole.push_back(std::make_shared<Matrix>("Quadrupole XZ", basisset_->nbf(), basisset_->nbf()));
    ao_Qpole.push_back(std::make_shared<Matrix>("Quadrupole YY", basisset_->nbf(), basisset_->nbf()));
    ao_Qpole.push_back(std::make_shared<Matrix>("Quadrupole YZ", basisset_->nbf(), basisset_->nbf()));
    ao_Qpole.push_back(std::make_shared<Matrix>("Quadrupole ZZ", basisset_->nbf(), basisset_->nbf()));

    // Form the one-electron integral objects from the integral factory
    std::shared_ptr<OneBodyAOInt> aoqOBI(integral_->ao_quadrupole());

    // Compute multipole moment integrals
    aoqOBI->set_origin(origin_);
    aoqOBI->compute(ao_Qpole);
    aoqOBI.reset();

    std::vector<SharedVector> quadrupole;
    quadrupole.push_back(std::make_shared<Vector>("Orbital Quadrupole XX", Ca->ncol()));
    quadrupole.push_back(std::make_shared<Vector>("Orbital Quadrupole YY", Ca->ncol()));
    quadrupole.push_back(std::make_shared<Vector>("Orbital Quadrupole ZZ", Ca->ncol()));

    if (same_orbs_) {
        // Quadrupoles
        for (int i = 0; i < Ca->ncol(); ++i) {
            double sumx = 0.0, sumy = 0.0, sumz = 0.0;
            for (int k = 0; k < Ca->nrow(); ++k) {
                for (int l = 0; l < Ca->nrow(); ++l) {
                    double tmp = Ca->get(0, k, i) * Ca->get(0, l, i);
                    sumx += ao_Qpole[0]->get(0, k, l) * tmp;
                    sumy += ao_Qpole[3]->get(0, k, l) * tmp;
                    sumz += ao_Qpole[5]->get(0, k, l) * tmp;
                }
            }

            quadrupole[0]->set(0, i, std::fabs(sumx));
            quadrupole[1]->set(0, i, std::fabs(sumy));
            quadrupole[2]->set(0, i, std::fabs(sumz));
        }

        std::vector<std::string> labels = basisset_->molecule()->irrep_labels();
        std::vector<std::tuple<double, int, int>> metric;
        for (int h = 0; h < epsilon_a_->nirrep(); h++) {
            for (int i = 0; i < epsilon_a_->dimpi()[h]; i++) {
                metric.push_back(std::tuple<double, int, int>(epsilon_a_->get(h, i), i, h));
            }
        }
        std::sort(metric.begin(), metric.end());
        if (print_output) {
            outfile->Printf("\n  Orbital extents (a.u.):\n");
            outfile->Printf("    %10s%15s%15s%15s%15s\n", "MO", "<x^2>", "<y^2>", "<z^2>", "<r^2>");
        }

        for (int i = 0; i < Ca->ncol(); i++) {
            int n = std::get<1>(metric[i]);
            int h = std::get<2>(metric[i]);

            double xx = quadrupole[0]->get(0, i), yy = quadrupole[1]->get(0, i), zz = quadrupole[2]->get(0, i);
            if (print_output) {
                outfile->Printf("    %4d%3s%3d%15.10f%15.10f%15.10f%15.10f\n", i, labels[h].c_str(), n,
                                std::fabs(quadrupole[0]->get(0, i)), std::fabs(quadrupole[1]->get(0, i)),
                                std::fabs(quadrupole[2]->get(0, i)), std::fabs(xx + yy + zz));
            }
            mo_es[0]->set(0, i, quadrupole[0]->get(0, i));
            mo_es[1]->set(0, i, quadrupole[1]->get(0, i));
            mo_es[2]->set(0, i, quadrupole[2]->get(0, i));
            mo_es[3]->set(0, i, fabs(xx + yy + zz));
        }
        if (print_output) {
            outfile->Printf("\n");
        }

    } else {
        // TODO: Both alpha and beta orbitals are reported separately
        // This helps identify symmetry breaking
    }
    return mo_es;
}

typedef PopulationAnalysisCalc PAC;
PopulationAnalysisCalc::PopulationAnalysisCalc(std::shared_ptr<Wavefunction> wfn) : Prop(wfn) {
    // No internal state. num_noon is now an argument.
}

PopulationAnalysisCalc::~PopulationAnalysisCalc() {}

void OEProp::compute_mulliken_charges() {
    PAC::SharedStdVector Qa, Qb, apcs;
    std::tie(Qa, Qb, apcs) = pac_.compute_mulliken_charges(true);
    wfn_->set_atomic_point_charges(apcs);

    auto vec_apcs = std::make_shared<Matrix>("Mulliken Charges: (a.u.)", 1, apcs->size());
    for (size_t i = 0; i < apcs->size(); i++) {
        vec_apcs->set(0, i, (*apcs)[i]);
    }
    wfn_->set_array_variable("MULLIKEN CHARGES", vec_apcs);
}

std::tuple<PAC::SharedStdVector, PAC::SharedStdVector, PAC::SharedStdVector>
PopulationAnalysisCalc::compute_mulliken_charges(bool print_output) {
    if (print_output) {
        outfile->Printf("  Mulliken Charges: (a.u.)\n");
    }
    std::shared_ptr<Molecule> mol = basisset_->molecule();

    auto Qa = std::make_shared<std::vector<double>>(mol->natom());
    auto Qb = std::make_shared<std::vector<double>>(mol->natom());

    auto apcs = std::make_shared<std::vector<double>>(mol->natom());
    std::vector<double> PSa(basisset_->nbf());
    double suma = 0.0;

    std::vector<double> PSb(basisset_->nbf());
    double sumb = 0.0;

    SharedMatrix Da;
    SharedMatrix Db;

    //    Get the Density Matrices for alpha and beta spins
    if (same_dens_) {
        Da = wfn_->matrix_subset_helper(Da_so_, Ca_so_, "AO", "D");
        Db = Da;
    } else {
        Da = wfn_->matrix_subset_helper(Da_so_, Ca_so_, "AO", "D alpha");
        Db = wfn_->matrix_subset_helper(Db_so_, Cb_so_, "AO", "D beta");
    }

    //    Compute the overlap matrix

    std::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    auto S = std::make_shared<Matrix>("S", basisset_->nbf(), basisset_->nbf());
    overlap->compute(S);

    //    Form the idempotent D*S matrix

    auto PSam = std::make_shared<Matrix>("PSa", basisset_->nbf(), basisset_->nbf());
    PSam->gemm(false, false, 1.0, Da, S, 0.0);
    auto PSbm = std::make_shared<Matrix>("PSb", basisset_->nbf(), basisset_->nbf());
    PSbm->gemm(false, false, 1.0, Db, S, 0.0);

    //     Accumulate registers

    for (int mu = 0; mu < basisset_->nbf(); mu++) {
        PSa[mu] = PSam->get(0, mu, mu);
        PSb[mu] = PSbm->get(0, mu, mu);

        int shell = basisset_->function_to_shell(mu);
        int A = basisset_->shell_to_center(shell);

        (*Qa)[A] += PSa[mu];
        (*Qb)[A] += PSb[mu];

        suma += PSa[mu];
        sumb += PSb[mu];
    }

    //    Print out the Mulliken populations and charges
    if (print_output) {
        outfile->Printf("   Center  Symbol    Alpha    Beta     Spin     Total\n");
    }
    double nuc = 0.0;
    for (int A = 0; A < mol->natom(); A++) {
        double Qs = (*Qa)[A] - (*Qb)[A];
        double Qt = mol->Z(A) - ((*Qa)[A] + (*Qb)[A]);
        (*apcs)[A] = Qt;
        if (print_output) {
            outfile->Printf("   %5d    %2s    %8.5f %8.5f %8.5f %8.5f\n", A + 1, mol->label(A).c_str(), (*Qa)[A],
                            (*Qb)[A], Qs, Qt);
        }
        nuc += (double)mol->Z(A);
    }

    if (print_output) {
        outfile->Printf("\n   Total alpha = %8.5f, Total beta = %8.5f, Total charge = %8.5f\n", suma, sumb,
                        nuc - suma - sumb);
    }

    if (print_output) outfile->Printf("\n");

    return std::make_tuple(Qa, Qb, apcs);
}
void OEProp::compute_lowdin_charges() {
    PAC::SharedStdVector Qa, Qb, apcs;
    std::tie(Qa, Qb, apcs) = pac_.compute_lowdin_charges(true);
    wfn_->set_atomic_point_charges(apcs);

    auto vec_apcs = std::make_shared<Matrix>("Lowdin Charges: (a.u.)", 1, apcs->size());
    for (size_t i = 0; i < apcs->size(); i++) {
        vec_apcs->set(0, i, (*apcs)[i]);
    }
    wfn_->set_array_variable("LOWDIN CHARGES", vec_apcs);
}

std::tuple<PAC::SharedStdVector, PAC::SharedStdVector, PAC::SharedStdVector>
PopulationAnalysisCalc::compute_lowdin_charges(bool print_output) {
    if (print_output) {
        outfile->Printf("  Lowdin Charges: (a.u.)\n");
    }
    std::shared_ptr<Molecule> mol = basisset_->molecule();

    auto Qa = std::make_shared<std::vector<double>>(mol->natom());
    auto Qb = std::make_shared<std::vector<double>>(mol->natom());

    auto apcs = std::make_shared<std::vector<double>>(mol->natom());

    double suma = 0.0;

    double sumb = 0.0;

    SharedMatrix Da;
    SharedMatrix Db;
    auto evecs = std::make_shared<Matrix>("Eigenvectors of S matrix", basisset_->nbf(), basisset_->nbf());
    auto temp = std::make_shared<Matrix>("Temporary matrix", basisset_->nbf(), basisset_->nbf());
    auto SDSa = std::make_shared<Matrix>("S_12 * D * S_12 alpha matrix", basisset_->nbf(), basisset_->nbf());
    auto SDSb = std::make_shared<Matrix>("S_12 * D * S_12 beta matrix", basisset_->nbf(), basisset_->nbf());
    auto evals = std::make_shared<Vector>(basisset_->nbf());

    //    Get the Density Matrices for alpha and beta spins
    if (same_dens_) {
        Da = wfn_->matrix_subset_helper(Da_so_, Ca_so_, "AO", "D");
        Db = Da;
    } else {
        Da = wfn_->matrix_subset_helper(Da_so_, Ca_so_, "AO", "D alpha");
        Db = wfn_->matrix_subset_helper(Db_so_, Cb_so_, "AO", "D beta");
    }

    //    Compute the overlap matrix

    std::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    auto S = std::make_shared<Matrix>("S", basisset_->nbf(), basisset_->nbf());
    overlap->compute(S);

    //    Form the S^(1/2) matrix

    S->power(1.0 / 2.0);

    //    Compute the S^(1/2)*D*S^(1/2) matrix

    temp->gemm(false, false, 1.0, Da, S, 0.0);
    SDSa->gemm(false, false, 1.0, S, temp, 0.0);
    temp->gemm(false, false, 1.0, Db, S, 0.0);
    SDSb->gemm(false, false, 1.0, S, temp, 0.0);

    //    Accumulate AO populations for each atom

    for (int mu = 0; mu < basisset_->nbf(); mu++) {
        int shell = basisset_->function_to_shell(mu);
        int A = basisset_->shell_to_center(shell);

        (*Qa)[A] += SDSa->get(0, mu, mu);
        (*Qb)[A] += SDSb->get(0, mu, mu);

        suma += SDSa->get(0, mu, mu);
        sumb += SDSb->get(0, mu, mu);
    }

    //    Print out the populations and charges

    if (print_output) {
        outfile->Printf("   Center  Symbol    Alpha    Beta     Spin     Total\n");
    }
    double nuc = 0.0;
    for (int A = 0; A < mol->natom(); A++) {
        double Qs = (*Qa)[A] - (*Qb)[A];
        double Qt = mol->Z(A) - ((*Qa)[A] + (*Qb)[A]);
        (*apcs)[A] = Qt;
        if (print_output) {
            outfile->Printf("   %5d    %2s    %8.5f %8.5f %8.5f %8.5f\n", A + 1, mol->label(A).c_str(), (*Qa)[A],
                            (*Qb)[A], Qs, Qt);
        }
        nuc += (double)mol->Z(A);
    }
    if (print_output) {
        outfile->Printf("\n  Total alpha = %8.5f, Total beta = %8.5f, Total charge = %8.5f\n", suma, sumb,
                        nuc - suma - sumb);
    }

    return std::make_tuple(Qa, Qb, apcs);
}

// See PopulationAnalysisCalc::compute_mbis_multipoles
void OEProp::compute_mbis_multipoles(bool free_atom_volumes) {
    SharedMatrix mpole, dpole, qpole, opole;
    std::tie(mpole, dpole, qpole, opole) = pac_.compute_mbis_multipoles(free_atom_volumes, true);

    wfn_->set_array_variable("MBIS CHARGES", mpole);
    wfn_->set_array_variable("MBIS DIPOLES", dpole);
    wfn_->set_array_variable("MBIS QUADRUPOLES", qpole);
    wfn_->set_array_variable("MBIS OCTUPOLES", opole);
}

/// Helper Methods for MBIS (JCTC, 2016, p. 3894-3912, Verstraelen et al.)

// Proatomic density of a specific shell of an atom  (Equation 7 in Verstraelen et al.)
double rho_ai_0(double n, double sigma, double distance) {
    return n * exp(-distance / sigma) / (pow(sigma, 3) * 8 * M_PI);
}

/* Updated parameters for initial MBIS params Nai and 1/Sai
   (from https://github.com/theochem/denspart/blob/740449c375f7e1c7b14cc8a978a8d0b1ed4617e3/denspart/mbis.py#L82 by Toon
   Verstraelen) */
// Attention: MBIS is not supported for elements with atomic number greater than 36
std::tuple<double, double> get_mbis_params(int atomic_num, int shell_num) {
    std::map<int, std::vector<std::tuple<double, double>>> params;
    params[1] = {std::make_tuple(1.00000, 1.76216)};
    params[2] = {std::make_tuple(2.00000, 3.11975)};

    params[3] = {std::make_tuple(1.86359, 5.56763), std::make_tuple(1.13641, 0.80520)};
    params[4] = {std::make_tuple(1.75663, 8.15111), std::make_tuple(2.24337, 1.22219)};
    params[5] = {std::make_tuple(1.73486, 10.46135), std::make_tuple(3.26514, 1.51797)};
    params[6] = {std::make_tuple(1.70730, 12.79758), std::make_tuple(4.29270, 1.85580)};
    params[7] = {std::make_tuple(1.68283, 15.13096), std::make_tuple(5.31717, 2.19942)};
    params[8] = {std::make_tuple(1.66122, 17.46129), std::make_tuple(6.33878, 2.54326)};
    params[9] = {std::make_tuple(1.64171, 19.78991), std::make_tuple(7.35829, 2.88601)};
    params[10] = {std::make_tuple(1.62380, 22.11938), std::make_tuple(8.37620, 3.22746)};

    params[11] = {std::make_tuple(1.48140, 25.82522), std::make_tuple(8.28761, 4.02120),
                  std::make_tuple(1.23098, 0.80897)};
    params[12] = {std::make_tuple(1.39674, 29.19802), std::make_tuple(8.10904, 4.76791),
                  std::make_tuple(2.49422, 1.08302)};
    params[13] = {std::make_tuple(1.34503, 32.33363), std::make_tuple(8.12124, 5.42812),
                  std::make_tuple(3.53372, 1.15994)};
    params[14] = {std::make_tuple(1.28865, 35.65432), std::make_tuple(7.98931, 6.17545),
                  std::make_tuple(4.72204, 1.33797)};
    params[15] = {std::make_tuple(1.23890, 39.00531), std::make_tuple(7.83125, 6.95265),
                  std::make_tuple(5.92985, 1.52690)};
    params[16] = {std::make_tuple(1.19478, 42.38177), std::make_tuple(7.66565, 7.75584),
                  std::make_tuple(7.13957, 1.71687)};
    params[17] = {std::make_tuple(1.15482, 45.79189), std::make_tuple(7.50031, 8.58542),
                  std::make_tuple(8.34487, 1.90546)};
    params[18] = {std::make_tuple(1.11803, 49.24317), std::make_tuple(7.33917, 9.44200),
                  std::make_tuple(9.54280, 2.09210)};

    params[19] = {std::make_tuple(1.09120, 52.59376), std::make_tuple(7.15086, 10.29851),
                  std::make_tuple(9.57061, 2.42121), std::make_tuple(1.18733, 0.67314)};
    params[20] = {std::make_tuple(1.07196, 55.86008), std::make_tuple(7.01185, 11.11887),
                  std::make_tuple(9.29555, 2.76621), std::make_tuple(2.62063, 0.88692)};
    params[21] = {std::make_tuple(1.05870, 59.04659), std::make_tuple(6.96404, 11.86718),
                  std::make_tuple(9.97866, 2.93024), std::make_tuple(2.99860, 0.98040)};
    params[22] = {std::make_tuple(1.04755, 62.22091), std::make_tuple(6.90438, 12.62229),
                  std::make_tuple(10.84355, 3.08264), std::make_tuple(3.20452, 1.05403)};
    params[23] = {std::make_tuple(1.03828, 65.38117), std::make_tuple(6.83516, 13.38417),
                  std::make_tuple(11.79532, 3.23508), std::make_tuple(3.33124, 1.11609)};
    params[24] = {std::make_tuple(1.03069, 68.52633), std::make_tuple(6.75998, 14.15132),
                  std::make_tuple(12.79256, 3.38991), std::make_tuple(3.41677, 1.17116)};
    params[25] = {std::make_tuple(1.02450, 71.65908), std::make_tuple(6.68141, 14.92337),
                  std::make_tuple(13.81149, 3.54730), std::make_tuple(3.48260, 1.22220)};
    params[26] = {std::make_tuple(1.01960, 74.77846), std::make_tuple(6.60101, 15.69935),
                  std::make_tuple(14.84330, 3.70685), std::make_tuple(3.53609, 1.27026)};
    params[27] = {std::make_tuple(1.01575, 77.88779), std::make_tuple(6.51976, 16.47941),
                  std::make_tuple(15.88061, 3.86829), std::make_tuple(3.58388, 1.31647)};
    params[28] = {std::make_tuple(1.01282, 80.98814), std::make_tuple(6.43837, 17.26336),
                  std::make_tuple(16.92012, 4.03115), std::make_tuple(3.62869, 1.36133)};
    params[29] = {std::make_tuple(1.01839, 83.81831), std::make_tuple(6.47823, 17.85149),
                  std::make_tuple(18.65720, 4.05312), std::make_tuple(2.84618, 1.37570)};
    params[30] = {std::make_tuple(1.00931, 87.16777), std::make_tuple(6.27682, 18.84319),
                  std::make_tuple(18.99747, 4.35989), std::make_tuple(3.71640, 1.44857)};
    params[31] = {std::make_tuple(1.00600, 90.34057), std::make_tuple(6.16315, 19.71091),
                  std::make_tuple(19.81836, 4.57852), std::make_tuple(4.01249, 1.29122)};
    params[32] = {std::make_tuple(0.99467, 93.80965), std::make_tuple(5.91408, 20.85993),
                  std::make_tuple(19.89501, 4.95158), std::make_tuple(5.19624, 1.39361)};
    params[33] = {std::make_tuple(0.98548, 97.22822), std::make_tuple(5.68319, 22.01684),
                  std::make_tuple(19.83497, 5.33969), std::make_tuple(6.49637, 1.51963)};
    params[34] = {std::make_tuple(0.97822, 100.60094), std::make_tuple(5.47209, 23.17528),
                  std::make_tuple(19.68845, 5.73803), std::make_tuple(7.86124, 1.65366)};
    params[35] = {std::make_tuple(0.97231, 103.94730), std::make_tuple(5.27765, 24.33975),
                  std::make_tuple(19.48822, 6.14586), std::make_tuple(9.26182, 1.78869)};
    params[36] = {std::make_tuple(0.96735, 107.28121), std::make_tuple(5.09646, 25.51581),
                  std::make_tuple(19.25332, 6.56380), std::make_tuple(10.68288, 1.92256)};

    return params[atomic_num][shell_num];
}

// A Helper Method That Calculates Radial Moments using Atomic Electron Densities derived from Charge Partitioning
std::vector<SharedMatrix> compute_radial_moments(const std::shared_ptr<DFTGrid>& grid,
                                                 const std::vector<double>& rho_a_points,
                                                 const std::vector<double>& distances, int num_atoms) {
    Options& options = Process::environment.options;
    const int max_power = std::max(4, options.get_int("MAX_RADIAL_MOMENT"));

    std::vector<SharedMatrix> rmoms;

    for (int n = 2; n <= max_power; n++) {
        std::stringstream sstream;
        sstream << "ATOMIC RADIAL MOMENTS <R^" << n << ">";
        auto mat_name = sstream.str();

        rmoms.push_back(std::make_shared<Matrix>(mat_name, num_atoms, 1));
    }

    auto blocks = grid->blocks();
    size_t total_points = grid->npoints();

#pragma omp parallel for
    for (int a = 0; a < num_atoms; a++) {
        for (int n = 2; n <= max_power; n++) {
            size_t running_points = 0;
            double val = 0;
            for (int b = 0; b < blocks.size(); b++) {
                auto block = blocks[b];
                SharedVector rho_block;
                size_t num_points = block->npoints();

                double* x = block->x();
                double* y = block->y();
                double* z = block->z();
                double* w = block->w();

                for (size_t p = running_points; p < running_points + num_points; p++) {
                    auto ap = a * total_points + p;
                    val += w[p - running_points] * rho_a_points[ap] * pow(distances[ap], n);
                }
                running_points += num_points;
            }
            rmoms[n - 2]->set(a, 0, val);
        }
    }

    return rmoms;
}

// Minimal Basis Iterative Stockholder (JCTC, 2016, p. 3894-3912, Verstraelen et al.)
std::tuple<SharedMatrix, SharedMatrix, SharedMatrix, SharedMatrix> PopulationAnalysisCalc::compute_mbis_multipoles(
    bool free_atom_volumes, bool print_output) {
    if (print_output) outfile->Printf("  ==> Computing MBIS Charges <==\n\n");
    timer_on("MBIS");

    // => Setup 1RDM on DFTGrid <= //

    Options& options = Process::environment.options;
    const int max_iter = options.get_int("MBIS_MAXITER");
    const double conv = options.get_double("MBIS_D_CONVERGENCE");
    const int debug = options.get_int("DEBUG");
    std::shared_ptr<Molecule> mol = basisset_->molecule();

    // MBIS grid options
    std::map<std::string, int> mbis_grid_options_int;
    std::map<std::string, std::string> mbis_grid_options_str;

    mbis_grid_options_int["DFT_RADIAL_POINTS"] = options.get_int("MBIS_RADIAL_POINTS");
    mbis_grid_options_int["DFT_SPHERICAL_POINTS"] = options.get_int("MBIS_SPHERICAL_POINTS");
    mbis_grid_options_str["DFT_PRUNING_SCHEME"] = options.get_str("MBIS_PRUNING_SCHEME");

    auto grid = std::make_shared<DFTGrid>(mol, basisset_, mbis_grid_options_int, mbis_grid_options_str, options);

    if (print_output && debug >= 1) grid->print();

    int nbf = basisset_->nbf();
    int num_atoms = mol->natom();
    size_t total_points = grid->npoints();

    SharedMatrix Da = wfn_->Da_subset("AO");
    SharedMatrix Db;
    auto point_func = (std::shared_ptr<PointFunctions>)(std::make_shared<RKSFunctions>(basisset_, total_points, nbf));

    if (same_dens_) {
        Db = Da->clone();
        point_func->set_pointers(Da);
    } else {
        Db = wfn_->Db_subset("AO");
        point_func = (std::shared_ptr<PointFunctions>)(std::make_shared<UKSFunctions>(basisset_, total_points, nbf));
        point_func->set_pointers(Da, Db);
    }

    std::vector<std::shared_ptr<BlockOPoints>> blocks = grid->blocks();
    size_t running_points = 0;

    // Coordinates, weights, and rho (molecular electron density) at each grid point
    std::vector<double> x_points(total_points, 0.0);
    std::vector<double> y_points(total_points, 0.0);
    std::vector<double> z_points(total_points, 0.0);
    std::vector<double> weights(total_points, 0.0);
    std::vector<double> rho(total_points, 0.0);

    for (int b = 0; b < blocks.size(); b++) {
        std::shared_ptr<BlockOPoints> block = blocks[b];
        SharedVector rho_block;
        size_t num_points = block->npoints();

        point_func->set_max_functions(block->local_nbf());
        point_func->set_max_points(num_points);
        point_func->compute_points(block);
        rho_block = point_func->point_values()["RHO_A"];

        if (!same_dens_) {
            rho_block->add(*point_func->point_values()["RHO_B"]);
        }

        auto x = block->x();
        auto y = block->y();
        auto z = block->z();
        auto w = block->w();

        for (size_t p = 0; p < num_points; p++) {
            x_points[running_points + p] = x[p];
            y_points[running_points + p] = y[p];
            z_points[running_points + p] = z[p];
            weights[running_points + p] = w[p];
            rho[running_points + p] = rho_block->get(p);
        }

        running_points += num_points;
    }

    // Electron count via numerical interagration
    double grid_electrons = 0.0;
    for (int p = 0; p < total_points; p++) {
        grid_electrons += weights[p] * rho[p];
    }

    // Actual electron count
    double mol_electrons = -1.0 * mol->molecular_charge();
    for (int a = 0; a < num_atoms; a++) {
        mol_electrons += mol->Z(a);
    }

    double electron_diff = grid_electrons - mol_electrons;

    if (print_output)
        outfile->Printf("  Electron Count from Grid (Expected Number): %8.5f (%8.5f)\n", grid_electrons, mol_electrons);
    if (print_output) outfile->Printf("  Difference: %8.5f\n\n", electron_diff);

    if (fabs(electron_diff) > 0.1) {
        throw PsiException("Number of electrons calculated using the grid (" +
                               std::to_string(round(100000.0 * grid_electrons) / 100000.0) +
                               ") does not match the number of electrons in the molecule. "
                               "Try increasing the number of radial or spherical points (mbis_radial_points and "
                               "mbis_spherical_points options).",
                           __FILE__, __LINE__);
    } else if (fabs(electron_diff) > 0.001) {
        outfile->Printf(
            "WARNING: The number of electrons calculated using the grid (%8.5f) differs from the number of electrons "
            "in the molecule by more than 0.001. "
            "Try increasing the number of radial or spherical points (mbis_radial_points and mbis_spherical_points "
            "options).\n\n",
            grid_electrons);
    }

    // Distances and displacements between grid points and nuclei
    std::vector<double> distances(num_atoms * total_points, 0.0);
    std::vector<std::vector<double>> disps(3, std::vector<double>(num_atoms * total_points, 0.0));

#pragma omp parallel for
    for (int atom = 0; atom < num_atoms; atom++) {
        for (size_t point = 0; point < total_points; point++) {
            Vector3 dr = Vector3(x_points[point], y_points[point], z_points[point]) - mol->xyz(atom);
            distances[atom * total_points + point] = dr.norm();
            disps[0][atom * total_points + point] = dr[0];
            disps[1][atom * total_points + point] = dr[1];
            disps[2][atom * total_points + point] = dr[2];
        }
    }

    // => Setup Proatom Basis Functions <= //

    // mA is the number of shells in an isolated atom
    std::vector<int> mA(num_atoms);
    for (int atom = 0; atom < num_atoms; atom++) {
        int atomic_num = static_cast<int>(mol->Z(atom));
        if (atomic_num <= 2) {
            mA[atom] = 1;
        } else if (atomic_num <= 10) {
            mA[atom] = 2;
        } else if (atomic_num <= 18) {
            mA[atom] = 3;
        } else if (atomic_num <= 36) {
            mA[atom] = 4;
        } else {
            throw PsiException("Atomic Number " + std::to_string(atomic_num) + " unsupported by MBIS", __FILE__,
                               __LINE__);
        }
    }

    // Atomic shell populations (N) and widths (S) (up to 4 shells per atom), from equations 18 and 19 in Verstraelen et
    // al.
    std::vector<std::vector<double>> Nai(num_atoms, std::vector<double>(4, 0.0));
    std::vector<std::vector<double>> Sai(num_atoms, std::vector<double>(4, 0.0));

    // Next iteration populations and widths
    std::vector<std::vector<double>> Nai_next(num_atoms, std::vector<double>(4, 0.0));
    std::vector<std::vector<double>> Sai_next(num_atoms, std::vector<double>(4, 0.0));

    // Population and width guesses, from get_mbis_params function
    for (int atom = 0; atom < num_atoms; atom++) {
        int atomic_num = static_cast<int>(mol->Z(atom));

        for (int m = 0; m < mA[atom]; m++) {
            std::tie(Nai[atom][m], Sai[atom][m]) = get_mbis_params(atomic_num, m);
            Sai[atom][m] = 1.0 / Sai[atom][m];

            if (print_output && debug >= 1) {
                outfile->Printf("  INITIAL ATOM %d, SHELL %d, POP %8.5f, WIDTH %8.5f\n", atom + 1, m + 1, Nai[atom][m],
                                Sai[atom][m]);
            }

            Nai_next[atom][m] = Nai[atom][m];
            Sai_next[atom][m] = Sai[atom][m];
        }
    }

    // Promolecular and proatomic densities
    std::vector<double> rho_0_points(total_points, 0.0);
    std::vector<double> rho_a_0_points(num_atoms * total_points, 0.0);

    // Next iteration densities
    std::vector<double> rho_0_points_next(total_points, 0.0);
    std::vector<double> rho_a_0_points_next(num_atoms * total_points, 0.0);

// Calculate initial proatom and promolecule density at all points
#pragma omp parallel for
    for (size_t point = 0; point < total_points; point++) {
        rho_0_points[point] = 0.0;
        for (int atom = 0; atom < num_atoms; atom++) {
            rho_a_0_points[atom * total_points + point] = 0.0;
            for (int m = 0; m < mA[atom]; m++) {
                rho_a_0_points[atom * total_points + point] +=
                    rho_ai_0(Nai[atom][m], Sai[atom][m], distances[atom * total_points + point]);
            }
            rho_0_points[point] += rho_a_0_points[atom * total_points + point];
        }
    }

    // => Main Stockholder Loop <= //

    int iter = 1;
    bool is_converged = false;
    double delta_rho_max_0;

    if (print_output && debug >= 1) outfile->Printf("                     Delta D\n");
    while (iter < max_iter) {
// Self-consistent update of population and density
#pragma omp parallel for
        for (int atom = 0; atom < num_atoms; atom++) {
            for (int m = 0; m < mA[atom]; m++) {
                double sum_n = 0.0;
                double sum_s = 0.0;

                for (int point = 0; point < total_points; point++) {
                    double rho_ai_0_point =
                        rho_ai_0(Nai[atom][m], Sai[atom][m], distances[atom * total_points + point]);
                    sum_n += weights[point] * rho[point] * rho_ai_0_point / rho_0_points[point];
                    sum_s += weights[point] * distances[atom * total_points + point] * rho[point] * rho_ai_0_point /
                             rho_0_points[point];
                }

                Nai_next[atom][m] = sum_n;
                Sai_next[atom][m] = sum_s / (3 * Nai_next[atom][m]);
            }
        }

        std::fill(rho_0_points_next.begin(), rho_0_points_next.end(), 0.0);
        std::fill(rho_a_0_points_next.begin(), rho_a_0_points_next.end(), 0.0);

#pragma omp parallel for
        for (size_t point = 0; point < total_points; point++) {
            for (int atom = 0; atom < num_atoms; atom++) {
                for (int m = 0; m < mA[atom]; m++) {
                    rho_a_0_points_next[atom * total_points + point] +=
                        rho_ai_0(Nai_next[atom][m], Sai_next[atom][m], distances[atom * total_points + point]);
                }
                rho_0_points_next[point] += rho_a_0_points_next[atom * total_points + point];
            }
        }

        // Convergence check (Equation 20 in Verstraelen et al.) and update of pro-densities
        std::vector<double> delta_rho_atoms_0(num_atoms, 0.0);

#pragma omp parallel for
        for (int atom = 0; atom < num_atoms; atom++) {
            double delta;
            for (size_t point = 0; point < total_points; point++) {
                delta = rho_a_0_points_next[atom * total_points + point] - rho_a_0_points[atom * total_points + point];
                delta_rho_atoms_0[atom] += weights[point] * delta * delta;
            }
            delta_rho_atoms_0[atom] = sqrt(delta_rho_atoms_0[atom]);
        }

        delta_rho_max_0 = 0.0;
        for (int atom = 0; atom < num_atoms; atom++) {
            if (delta_rho_atoms_0[atom] > delta_rho_max_0) delta_rho_max_0 = delta_rho_atoms_0[atom];
        }

        // Update populations, widths, and densities
        Nai = Nai_next;
        Sai = Sai_next;
        rho_0_points = rho_0_points_next;
        rho_a_0_points = rho_a_0_points_next;

        if (print_output && debug >= 1) outfile->Printf("   @MBIS iter %3d:  %.3e\n", iter, delta_rho_max_0);

        if (delta_rho_max_0 < conv) {
            if (print_output && debug >= 1) outfile->Printf("  MBIS Atomic Density Converged\n\n");
            is_converged = true;
            break;
        }

        iter += 1;
    }

    if (!is_converged) throw ConvergenceError<int>("MBIS", max_iter, conv, delta_rho_max_0, __FILE__, __LINE__);

    // Final population and width
    if (print_output && debug >= 1) {
        for (int atom = 0; atom < num_atoms; atom++) {
            for (int m = 0; m < mA[atom]; m++) {
                outfile->Printf("  FINAL ATOM %d, SHELL %d, POP %8.5f, WIDTH %8.5f\n", atom + 1, m + 1, Nai[atom][m],
                                Sai[atom][m]);
            }
        }
        outfile->Printf("\n");
    }    

    // => Post-Processing <= //

    // Atomic density, as defined in Equation 5 of Verstraelen et al.
    std::vector<double> rho_a(num_atoms * total_points, 0.0);

#pragma omp parallel for
    for (int atom = 0; atom < num_atoms; atom++) {
        for (size_t point = 0; point < total_points; point++) {
            rho_a[atom * total_points + point] =
                rho[point] * rho_a_0_points[atom * total_points + point] / rho_0_points[point];
        }
    }

    // Kronecker Delta
    std::vector<std::vector<double>> k_delta = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    // Non-redundant cartesian quadrupole components
    std::vector<std::vector<int>> qpole_inds = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};

    // Non-redundant cartesian octupole components
    std::vector<std::vector<int>> opole_inds = {{0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {0, 1, 1}, {0, 1, 2},
                                                {0, 2, 2}, {1, 1, 1}, {1, 1, 2}, {1, 2, 2}, {2, 2, 2}};

    // Non-redundant cartesian multipoles, as given in the convention in the HORTON software
    /* Commented out lines ending in _stone represent multipole conventions given in
       "The Theory of Intermolecular Forces (2nd Edition)" by Anthony Stone (A possible future addition) */
    auto mpole = std::make_shared<Matrix>("MBIS Charges: (a.u.)", num_atoms, 1);
    auto dpole = std::make_shared<Matrix>("MBIS Dipoles: (a.u.)", num_atoms, 3);
    auto qpole = std::make_shared<Matrix>("MBIS Quadrupoles: (a.u.)", num_atoms, 6);
    auto opole = std::make_shared<Matrix>("MBIS Octupoles: (a.u.)", num_atoms, 10);
// auto qpole_stone = std::make_shared<Matrix>("MBIS Quadrupoles: (a.u.)", num_atoms, 6);
// auto opole_stone = std::make_shared<Matrix>("MBIS Octupoles: (a.u.)", num_atoms, 10);
// Calculate atomic multipoles
#pragma omp parallel for
    for (int a = 0; a < num_atoms; a++) {
        mpole->add(a, 0, mol->Z(a));
        for (size_t p = 0; p < total_points; p++) {
            auto ap = a * total_points + p;

            // Atomic monopole
            mpole->add(a, 0, weights[p] * -rho_a[ap]);

            // Atomic dipole
            for (int i = 0; i < 3; i++) {
                dpole->add(a, i, weights[p] * -rho_a[ap] * disps[i][ap]);
            }

            // Atomic quadrupole
            for (int q = 0; q < 6; q++) {
                int i = qpole_inds[q][0], j = qpole_inds[q][1];
                qpole->add(a, q, weights[p] * -rho_a[ap] * (disps[i][ap] * disps[j][ap]));
                // qpole_stone->add(a, q, weights[p] * rho_a[ap] * (1.5 * disps[i][ap] * disps[j][ap] - 0.5 *
                // pow(distances[ap], 2) * k_delta[i][j]));
            }

            // Atomic octupole
            for (int o = 0; o < 10; o++) {
                int i = opole_inds[o][0], j = opole_inds[o][1], k = opole_inds[o][2];
                opole->add(a, o, weights[p] * -rho_a[ap] * (disps[i][ap] * disps[j][ap] * disps[k][ap]));
                // opole_stone->add(a, o, weights[p] * rho_a[ap]
                //    * (2.5 * disps[i][ap] * disps[j][ap] * disps[k][ap] - 0.5 * pow(distances[ap], 2)
                //        * (k_delta[j][k] * disps[i][ap] + k_delta[i][k] * disps[j][ap] + k_delta[i][j] *
                //        disps[k][ap])));
            }
        }
    }

    if (print_output) {
        outfile->Printf("  MBIS Charges: (a.u.)\n");
        outfile->Printf("   Center  Symbol  Z      Pop.       Charge\n");

        for (int a = 0; a < num_atoms; a++) {
            outfile->Printf("  %5d      %2s %4d   %9.6f   %9.6f\n", a + 1, mol->label(a).c_str(),
                            static_cast<int>(mol->Z(a)), mol->Z(a) - mpole->get(a, 0), mpole->get(a, 0));
        }

        outfile->Printf("\n  MBIS Dipoles: [e a0]\n");
        outfile->Printf("   Center  Symbol  Z        X           Y           Z\n");

        for (int a = 0; a < num_atoms; a++) {
            outfile->Printf("  %5d      %2s %4d   %9.6f   %9.6f   %9.6f\n", a + 1, mol->label(a).c_str(),
                            static_cast<int>(mol->Z(a)), dpole->get(a, 0), dpole->get(a, 1), dpole->get(a, 2));
        }

        outfile->Printf("\n  MBIS Quadrupoles: [e a0^2]\n");
        outfile->Printf("   Center  Symbol  Z      XX        XY        XZ        YY        YZ        ZZ\n");

        for (int a = 0; a < num_atoms; a++) {
            outfile->Printf("  %5d      %2s %4d   %7.4f   %7.4f   %7.4f   %7.4f   %7.4f   %7.4f\n", a + 1,
                            mol->label(a).c_str(), static_cast<int>(mol->Z(a)), qpole->get(a, 0), qpole->get(a, 1),
                            qpole->get(a, 2), qpole->get(a, 3), qpole->get(a, 4), qpole->get(a, 5));
        }

        outfile->Printf("\n  MBIS Octupoles: [e a0^3]\n");
        outfile->Printf(
            "   Center  Symbol  Z      XXX       XXY       XXZ       XYY       XYZ       XZZ       YYY       YYZ       "
            "YZZ       ZZZ\n");

        for (int a = 0; a < num_atoms; a++) {
            outfile->Printf(
                "  %5d      %2s %4d   %7.4f   %7.4f   %7.4f   %7.4f   %7.4f   %7.4f   %7.4f   %7.4f   %7.4f   %7.4f\n",
                a + 1, mol->label(a).c_str(), static_cast<int>(mol->Z(a)), opole->get(a, 0), opole->get(a, 1),
                opole->get(a, 2), opole->get(a, 3), opole->get(a, 4), opole->get(a, 5), opole->get(a, 6),
                opole->get(a, 7), opole->get(a, 8), opole->get(a, 9));
        }
    }

    const int max_power = options.get_int("MAX_RADIAL_MOMENT");
    auto rmoms = compute_radial_moments(grid, rho_a, distances, num_atoms);

    for (int n = 2; n <= max_power; n++) {
        std::stringstream sstream;
        sstream << "MBIS RADIAL MOMENTS <R^" << n << ">";

        auto var_name = sstream.str();
        wfn_->set_array_variable(var_name, rmoms[n - 2]);
    }

    auto valence_widths = std::make_shared<Matrix>("MBIS Valence Widths", num_atoms, 1);
    for (int atom = 0; atom < num_atoms; atom++) {
        valence_widths->set(atom, 0, Sai[atom][mA[atom] - 1]);
    }
    wfn_->set_array_variable("MBIS VALENCE WIDTHS", valence_widths);

    // Compute the volume widths, only for molecules
    bool free_atom = (num_atoms == 1);
    auto volume_ratios = valence_widths->clone();
    if (free_atom_volumes) {
        volume_ratios->zero();
        if (free_atom == false) {
            for (int a = 0; a < num_atoms; ++a) {
                double free_atom = wfn_->scalar_variable("MBIS FREE ATOM " + mol->label(a) + " VOLUME");
                double vr = rmoms[1]->get(a, 0) / free_atom;
                volume_ratios->set(a, 0, vr);
            }
            wfn_->set_array_variable("MBIS VOLUME RATIOS", volume_ratios);
        }
    }

    if (print_output) {
        outfile->Printf("\n  MBIS Radial Moments:");
        outfile->Printf("\n   Center  Symbol  Z      ");
        for (int n = 2; n <= max_power; n++) {
            outfile->Printf("[a0^%d]      ", n);
        }
        for (int a = 0; a < num_atoms; a++) {
            outfile->Printf("\n  %5d      %2s %4d", a + 1, mol->label(a).c_str(), static_cast<int>(mol->Z(a)));
            for (int n = 2; n <= max_power; n++) {
                outfile->Printf("   %9.6f", rmoms[n - 2]->get(a, 0));
            }
        }

        outfile->Printf("\n\n  MBIS Valence Widths: [a0]\n");
        outfile->Printf("   Center  Symbol  Z     Width\n");

        for (int a = 0; a < num_atoms; a++) {
            outfile->Printf("  %5d      %2s %4d   %9.6f\n", a + 1, mol->label(a).c_str(), static_cast<int>(mol->Z(a)),
                            valence_widths->get(a, 0));
        }
        if (free_atom_volumes) {
            if (free_atom == false) {
                outfile->Printf("\n\n  MBIS Volume Ratios: \n");
                outfile->Printf("   Center  Symbol  Z     \n");
                for (int a = 0; a < num_atoms; a++) {
                    outfile->Printf("  %5d      %2s %4d   %9.6f\n", a + 1, mol->label(a).c_str(),
                                    static_cast<int>(mol->Z(a)), volume_ratios->get(a, 0));
                }
            }
        }
    }

    timer_off("MBIS");

    return std::make_tuple(mpole, dpole, qpole, opole);
}

void OEProp::compute_mayer_indices() {
    SharedMatrix MBI_total, MBI_alpha, MBI_beta;
    SharedVector MBI_valence;
    std::tie(MBI_total, MBI_alpha, MBI_beta, MBI_valence) = pac_.compute_mayer_indices(true);

    wfn_->set_array_variable("MAYER INDICES", MBI_total);
}

std::tuple<SharedMatrix, SharedMatrix, SharedMatrix, SharedVector> PopulationAnalysisCalc::compute_mayer_indices(
    bool print_output) {
    if (print_output) {
        outfile->Printf("\n\n  Mayer Bond Indices:\n\n");
    }

    std::shared_ptr<Molecule> mol = basisset_->molecule();

    SharedMatrix bis;

    int nbf = basisset_->nbf();

    SharedMatrix Da;  // Density matrix for alpha spin
    SharedMatrix Db;  // Density matrix for beta spin

    auto DSa = std::make_shared<Matrix>("D * S alpha matrix", nbf, nbf);
    auto DSb = std::make_shared<Matrix>("D * S beta matrix", nbf, nbf);

    //    Get the Density Matrices for alpha and beta spins
    if (same_dens_) {
        Da = wfn_->matrix_subset_helper(Da_so_, Ca_so_, "AO", "D");
        Db = Da;
    } else {
        Da = wfn_->matrix_subset_helper(Da_so_, Ca_so_, "AO", "D alpha");
        Db = wfn_->matrix_subset_helper(Db_so_, Cb_so_, "AO", "D beta");
    }

    //    Compute the overlap matrix

    std::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    auto S = std::make_shared<Matrix>("S matrix", nbf, nbf);
    overlap->compute(S);

    //    Form the idempotent D*S matrix

    DSa->gemm(false, false, 1.0, Da, S, 0.0);
    DSb->gemm(false, false, 1.0, Db, S, 0.0);

    //     Compute Mayer bond indices

    int natom = mol->natom();

    auto MBI_total = std::make_shared<Matrix>(natom, natom);
    SharedMatrix MBI_alpha;
    SharedMatrix MBI_beta;

    if (!same_dens_) {
        MBI_alpha = std::make_shared<Matrix>(natom, natom);
        MBI_beta = std::make_shared<Matrix>(natom, natom);
    }

    for (int mu = 0; mu < nbf; mu++) {
        for (int nu = 0; nu < mu; nu++) {
            int shell_mu = basisset_->function_to_shell(mu);
            int shell_nu = basisset_->function_to_shell(nu);
            int atom_mu = basisset_->shell_to_center(shell_mu);
            int atom_nu = basisset_->shell_to_center(shell_nu);
            if (atom_mu == atom_nu) continue;

            double alpha = DSa->get(0, mu, nu) * DSa->get(0, nu, mu);
            double beta = DSb->get(0, mu, nu) * DSb->get(0, nu, mu);
            MBI_total->add(0, atom_mu, atom_nu, 2 * (alpha + beta));
            MBI_total->add(0, atom_nu, atom_mu, 2 * (alpha + beta));

            if (!same_dens_) {
                MBI_alpha->add(0, atom_mu, atom_nu, 2 * (alpha));
                MBI_alpha->add(0, atom_nu, atom_mu, 2 * (alpha));
                MBI_beta->add(0, atom_mu, atom_nu, 2 * (beta));
                MBI_beta->add(0, atom_nu, atom_mu, 2 * (beta));
            }
        }
    }

    //    Compute valences

    auto MBI_valence = std::make_shared<Vector>(natom);

    for (int iat = 0; iat < natom; iat++) {
        for (int jat = 0; jat < natom; jat++) {
            double valence = MBI_valence->get(0, iat);
            MBI_valence->set(0, iat, valence + MBI_total->get(0, iat, jat));
        }
    }

    //    Print out the bond index matrix and valences

    //    Note: The computed Mayer bond indices (MBI) will be different from the MBI values computed using
    //    some other program packages for the unrestricted case. The reason is that these programs left out
    //    the spin density contribution in the MBI equation for the unrestricted wavefunctions. As the result,
    //    the MBI value will be underestimated. For example, the MBI value for the H-H bond of H2+
    //    is calculated to be 0.25 using the NBO program. The equation coded above gives the correct value of 0.5.
    //    For reference, see IJQC 29 (1986) P. 73 and IJQC 29 (1986) P. 477.

    //    A nicer output is needed ...

    if (print_output) {
        if (same_dens_) {
            MBI_total->print();
            outfile->Printf("  Atomic Valences: \n");
            MBI_valence->print();
        } else {
            outfile->Printf("  Total Bond Index: \n");
            MBI_total->print();
            outfile->Printf("  Alpha Contribution: \n");
            MBI_alpha->print();
            outfile->Printf("  Beta Contribution: \n");
            MBI_beta->print();
            outfile->Printf("  Atomic Valences: \n");
            MBI_valence->print();
        }
    }

    return std::make_tuple(MBI_total, MBI_alpha, MBI_beta, MBI_valence);
}
void OEProp::compute_wiberg_lowdin_indices() {
    SharedMatrix WBI_total, WBI_alpha, WBI_beta;
    SharedVector WBI_valence;
    std::tie(WBI_total, WBI_alpha, WBI_beta, WBI_valence) = pac_.compute_wiberg_lowdin_indices(true);
    wfn_->set_array_variable("WIBERG LOWDIN INDICES", WBI_total);
}

std::tuple<SharedMatrix, SharedMatrix, SharedMatrix, SharedVector>
PopulationAnalysisCalc::compute_wiberg_lowdin_indices(bool print_output) {
    if (print_output) {
        outfile->Printf("\n\n  Wiberg Bond Indices using Orthogonal Lowdin Orbitals:\n\n");
    }

    //    We may wanna get rid of these if we have NAOs...

    std::shared_ptr<Molecule> mol = basisset_->molecule();

    int nbf = basisset_->nbf();

    SharedMatrix Da;
    SharedMatrix Db;
    auto evecs = std::make_shared<Matrix>("Eigenvectors of S matrix", nbf, nbf);
    auto temp = std::make_shared<Matrix>("Temporary matrix", nbf, nbf);
    auto SDSa = std::make_shared<Matrix>("S_12 * D * S_12 alpha matrix", nbf, nbf);
    auto SDSb = std::make_shared<Matrix>("S_12 * D * S_12 beta matrix", nbf, nbf);
    auto evals = std::make_shared<Vector>(nbf);

    //    Get the Density Matrices for alpha and beta spins
    if (same_dens_) {
        Da = wfn_->matrix_subset_helper(Da_so_, Ca_so_, "AO", "D");
        Db = Da;
    } else {
        Da = wfn_->matrix_subset_helper(Da_so_, Ca_so_, "AO", "D alpha");
        Db = wfn_->matrix_subset_helper(Db_so_, Cb_so_, "AO", "D beta");
    }

    //    Compute the overlap matrix

    std::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    auto S = std::make_shared<Matrix>("S", basisset_->nbf(), basisset_->nbf());
    overlap->compute(S);

    //    Form the S^(1/2) matrix

    S->diagonalize(evecs, evals);
    S->zero();
    for (int p = 0; p < basisset_->nbf(); ++p) S->set(0, p, p, sqrt(evals->get(0, p)));
    S->back_transform(evecs);

    //    Compute the S^(1/2)*D*S^(1/2) matrix

    temp->gemm(false, false, 1.0, Da, S, 0.0);
    SDSa->gemm(false, false, 1.0, S, temp, 0.0);
    temp->gemm(false, false, 1.0, Db, S, 0.0);
    SDSb->gemm(false, false, 1.0, S, temp, 0.0);

    //    Compute Wiberg bond indices

    int natom = mol->natom();

    auto WBI_total = std::make_shared<Matrix>(natom, natom);
    SharedMatrix WBI_alpha;
    SharedMatrix WBI_beta;

    if (!same_dens_) {
        WBI_alpha = std::make_shared<Matrix>(natom, natom);
        WBI_beta = std::make_shared<Matrix>(natom, natom);
    }

    for (int mu = 0; mu < nbf; mu++) {
        for (int nu = 0; nu < mu; nu++) {
            int shell_mu = basisset_->function_to_shell(mu);
            int shell_nu = basisset_->function_to_shell(nu);
            int atom_mu = basisset_->shell_to_center(shell_mu);
            int atom_nu = basisset_->shell_to_center(shell_nu);
            if (atom_mu == atom_nu) continue;

            double alpha = SDSa->get(0, mu, nu) * SDSa->get(0, nu, mu);
            double beta = SDSb->get(0, mu, nu) * SDSb->get(0, nu, mu);
            WBI_total->add(0, atom_mu, atom_nu, 2 * (alpha + beta));
            WBI_total->add(0, atom_nu, atom_mu, 2 * (alpha + beta));

            if (!same_dens_) {
                WBI_alpha->add(0, atom_mu, atom_nu, 2 * (alpha));
                WBI_alpha->add(0, atom_nu, atom_mu, 2 * (alpha));
                WBI_beta->add(0, atom_mu, atom_nu, 2 * (beta));
                WBI_beta->add(0, atom_nu, atom_mu, 2 * (beta));
            }
        }
    }

    //    Compute valences

    auto WBI_valence = std::make_shared<Vector>(natom);

    for (int iat = 0; iat < natom; iat++) {
        for (int jat = 0; jat < natom; jat++) {
            double valence = WBI_valence->get(0, iat);
            WBI_valence->set(0, iat, valence + WBI_total->get(0, iat, jat));
        }
    }

    //    Print out the bond index matrix
    //    A nicer output is needed ...
    if (print_output) {
        if (same_dens_) {
            WBI_total->print();
            outfile->Printf("  Atomic Valences: \n");
            WBI_valence->print();
        } else {
            outfile->Printf("  Total Bond Index: \n");
            WBI_total->print();
            outfile->Printf("  Alpha Contribution: \n");
            WBI_alpha->print();
            outfile->Printf("  Beta Contribution: \n");
            WBI_beta->print();
            outfile->Printf("  Atomic Valences: \n");
            WBI_valence->print();
        }
    }
    return std::make_tuple(WBI_total, WBI_alpha, WBI_beta, WBI_valence);
}

void OEProp::compute_no_occupations() {
    auto metrics = pac_.compute_no_occupations(max_noon_, true);
    wfn_->set_no_occupations(*metrics);
}
std::shared_ptr<std::vector<std::vector<std::tuple<double, int, int>>>> PopulationAnalysisCalc::compute_no_occupations(
    int max_noon, bool print_output) {
    std::shared_ptr<std::vector<std::vector<std::tuple<double, int, int>>>> metrics =
        std::make_shared<std::vector<std::vector<std::tuple<double, int, int>>>>();
    std::vector<std::string> labels = basisset_->molecule()->irrep_labels();

    if (print_output) {
        outfile->Printf("  Natural Orbital Occupations:\n\n");
    }

    // Terminally, it will be [metric_a , metric_b, metric] or [metric] depending on same_dens

    if (!same_dens_) {
        SharedVector Oa;
        SharedVector Ob;
        if (same_dens_) {
            std::pair<SharedMatrix, SharedVector> vals = Na_mo();
            Oa = vals.second;
            Ob = vals.second;
        } else {
            std::pair<SharedMatrix, SharedVector> vals = Na_mo();
            Oa = vals.second;
            std::pair<SharedMatrix, SharedVector> vals2 = Nb_mo();
            Ob = vals2.second;
        }
        std::vector<std::tuple<double, int, int>> metric_a;
        for (int h = 0; h < Oa->nirrep(); h++) {
            for (int i = 0; i < Oa->dimpi()[h]; i++) {
                metric_a.push_back(std::tuple<double, int, int>(Oa->get(h, i), i, h));
            }
        }

        metrics->push_back(metric_a);

        std::sort(metric_a.begin(), metric_a.end(), std::greater<std::tuple<double, int, int>>());
        int offset_a = wfn_->nalpha();
        int start_occ_a = offset_a - max_noon;
        start_occ_a = (start_occ_a < 0 ? 0 : start_occ_a);
        int stop_vir_a = offset_a + max_noon + 1;
        stop_vir_a = (int)((size_t)stop_vir_a >= metric_a.size() ? metric_a.size() : stop_vir_a);

        if (print_output) {
            outfile->Printf("  Alpha Occupations:\n");

            for (int index = start_occ_a; index < stop_vir_a; index++) {
                if (index < offset_a) {
                    outfile->Printf("  HONO-%-2d: %4d%3s %8.3f\n", offset_a - index - 1,
                                    std::get<1>(metric_a[index]) + 1, labels[std::get<2>(metric_a[index])].c_str(),
                                    std::get<0>(metric_a[index]));
                } else {
                    outfile->Printf("  LUNO+%-2d: %4d%3s %8.3f\n", index - offset_a, std::get<1>(metric_a[index]) + 1,
                                    labels[std::get<2>(metric_a[index])].c_str(), std::get<0>(metric_a[index]));
                }
            }
            outfile->Printf("\n");
        }

        std::vector<std::tuple<double, int, int>> metric_b;
        for (int h = 0; h < Ob->nirrep(); h++) {
            for (int i = 0; i < Ob->dimpi()[h]; i++) {
                metric_b.push_back(std::tuple<double, int, int>(Ob->get(h, i), i, h));
            }
        }

        metrics->push_back(metric_b);

        std::sort(metric_b.begin(), metric_b.end(), std::greater<std::tuple<double, int, int>>());

        int offset_b = wfn_->nbeta();
        int start_occ_b = offset_b - max_noon;
        start_occ_b = (start_occ_b < 0 ? 0 : start_occ_b);
        int stop_vir_b = offset_b + max_noon + 1;
        stop_vir_b = (int)((size_t)stop_vir_b >= metric_b.size() ? metric_b.size() : stop_vir_b);

        if (print_output) {
            outfile->Printf("  Beta Occupations:\n");
            for (int index = start_occ_b; index < stop_vir_b; index++) {
                if (index < offset_b) {
                    outfile->Printf("  HONO-%-2d: %4d%3s %8.3f\n", offset_b - index - 1,
                                    std::get<1>(metric_b[index]) + 1, labels[std::get<2>(metric_b[index])].c_str(),
                                    std::get<0>(metric_b[index]));
                } else {
                    outfile->Printf("  LUNO+%-2d: %4d%3s %8.3f\n", index - offset_b, std::get<1>(metric_b[index]) + 1,
                                    labels[std::get<2>(metric_b[index])].c_str(), std::get<0>(metric_b[index]));
                }
            }
            outfile->Printf("\n");
        }
    }

    std::pair<SharedMatrix, SharedVector> vals = Nt_mo();
    SharedVector Ot = vals.second;

    std::vector<std::tuple<double, int, int>> metric;
    for (int h = 0; h < Ot->nirrep(); h++) {
        for (int i = 0; i < Ot->dimpi()[h]; i++) {
            metric.push_back(std::tuple<double, int, int>(Ot->get(h, i), i, h));
        }
    }

    metrics->push_back(metric);

    std::sort(metric.begin(), metric.end(), std::greater<std::tuple<double, int, int>>());

    int offset = wfn_->nbeta();
    int start_occ = offset - max_noon;
    start_occ = (start_occ < 0 ? 0 : start_occ);
    int stop_vir = offset + max_noon + 1;
    stop_vir = (int)((size_t)stop_vir >= metric.size() ? metric.size() : stop_vir);

    if (print_output) {
        outfile->Printf("  Total Occupations:\n");
        for (int index = start_occ; index < stop_vir; index++) {
            if (index < offset) {
                outfile->Printf("  HONO-%-2d: %4d%3s %8.3f\n", offset - index - 1, std::get<1>(metric[index]) + 1,
                                labels[std::get<2>(metric[index])].c_str(), std::get<0>(metric[index]));
            } else {
                outfile->Printf("  LUNO+%-2d: %4d%3s %8.3f\n", index - offset, std::get<1>(metric[index]) + 1,
                                labels[std::get<2>(metric[index])].c_str(), std::get<0>(metric[index]));
            }
        }
        outfile->Printf("\n");
    }
    return metrics;
    // for(int h = 0; h < epsilon_a_->nirrep(); h++) free(labels[h]); free(labels);
}

void OEProp::set_wavefunction(std::shared_ptr<Wavefunction> wfn) {
    mpc_.set_wavefunction(wfn);
    pac_.set_wavefunction(wfn);
    epc_.set_wavefunction(wfn);
}

void OEProp::set_restricted(bool restricted) {
    mpc_.set_restricted(restricted);
    pac_.set_restricted(restricted);
    epc_.set_restricted(restricted);
}

void OEProp::set_epsilon_a(SharedVector epsilon_a) {
    mpc_.set_epsilon_a(epsilon_a);
    pac_.set_epsilon_a(epsilon_a);
    epc_.set_epsilon_a(epsilon_a);
}

void OEProp::set_epsilon_b(SharedVector epsilon_b) {
    mpc_.set_epsilon_b(epsilon_b);
    pac_.set_epsilon_b(epsilon_b);
    epc_.set_epsilon_b(epsilon_b);
}

void OEProp::set_Ca(SharedMatrix Ca) {
    mpc_.set_Ca(Ca);
    pac_.set_Ca(Ca);
    epc_.set_Ca(Ca);
}

void OEProp::set_Cb(SharedMatrix Cb) {
    mpc_.set_Cb(Cb);
    pac_.set_Cb(Cb);
    epc_.set_Cb(Cb);
}

void OEProp::set_Da_ao(SharedMatrix Da, int symmetry) {
    mpc_.set_Da_ao(Da, symmetry);
    pac_.set_Da_ao(Da, symmetry);
    epc_.set_Da_ao(Da, symmetry);
}

void OEProp::set_Db_ao(SharedMatrix Db, int symmetry) {
    mpc_.set_Db_ao(Db, symmetry);
    pac_.set_Db_ao(Db, symmetry);
    epc_.set_Db_ao(Db, symmetry);
}

void OEProp::set_Da_so(SharedMatrix Da) {
    mpc_.set_Da_so(Da);
    pac_.set_Da_so(Da);
    epc_.set_Da_so(Da);
}

void OEProp::set_Db_so(SharedMatrix Db) {
    mpc_.set_Db_so(Db);
    pac_.set_Db_so(Db);
    epc_.set_Db_so(Db);
}

void OEProp::set_Da_mo(SharedMatrix Da) {
    mpc_.set_Da_mo(Da);
    pac_.set_Da_mo(Da);
    epc_.set_Da_mo(Da);
}

void OEProp::set_Db_mo(SharedMatrix Db) {
    mpc_.set_Db_mo(Db);
    pac_.set_Db_mo(Db);
    epc_.set_Db_mo(Db);
}

}  // namespace psi
