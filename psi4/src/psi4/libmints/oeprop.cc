/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include <Eigen/Core>

#include <limits>

#ifdef USING_LAPACK_MKL
#include <mkl.h>
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
    if ((tasks_.count("LOWDIN_CHARGES")) || (tasks_.count("LOWDIN_SPINS"))) compute_lowdin_charges();
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
    for (size_t atom1 = 0; atom1 < nesps->size(); ++atom1) {
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
    PAC::SharedStdVector Qa, Qb, apcs, asps;
    std::tie(Qa, Qb, apcs, asps) = pac_.compute_lowdin_charges(true);
    wfn_->set_atomic_point_charges(apcs);

    auto vec_apcs = std::make_shared<Matrix>("Lowdin Charges: (a.u.)", 1, apcs->size());
    for (size_t i = 0; i < apcs->size(); i++) {
        vec_apcs->set(0, i, (*apcs)[i]);
    }
    wfn_->set_array_variable("LOWDIN CHARGES", vec_apcs);

    auto vec_asps = std::make_shared<Matrix>("Lowdin Spins: (a.u.)", 1, apcs->size());
    for (size_t i = 0; i < asps->size(); i++) {
        vec_asps->set(0, i, (*asps)[i]);
    }
    wfn_->set_array_variable("LOWDIN SPINS", vec_asps);
}

std::tuple<PAC::SharedStdVector, PAC::SharedStdVector, PAC::SharedStdVector, PAC::SharedStdVector>
PopulationAnalysisCalc::compute_lowdin_charges(bool print_output) {
    if (print_output) {
        outfile->Printf("  Lowdin Charges: (a.u.)\n");
    }
    std::shared_ptr<Molecule> mol = basisset_->molecule();

    auto Qa = std::make_shared<std::vector<double>>(mol->natom());
    auto Qb = std::make_shared<std::vector<double>>(mol->natom());

    auto apcs = std::make_shared<std::vector<double>>(mol->natom());
    auto asps = std::make_shared<std::vector<double>>(mol->natom());

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
        (*asps)[A] = Qs;
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

    return std::make_tuple(Qa, Qb, apcs, asps);
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

/* Initial MBIS proatom parameters (N, 1/sigma) per shell, sourced from denspart:
   https://github.com/theochem/denspart/blob/main/src/denspart/mbis.py by Toon Verstraelen.
   Covers Z=1-118 (H through Og). Shell counts vary due to merging of near-degenerate
   shells during the free-atom optimization; the maximum is 7 shells (Z=87, 88, 115-118). */
const std::vector<std::tuple<double, double>>& get_mbis_params(int atomic_num) {
    static const std::map<int, std::vector<std::tuple<double, double>>> params = {
        // Period 1
        {1,  {{1.00000, 1.76216}}},
        {2,  {{2.00000, 3.11975}}},
        // Period 2
        {3,  {{1.86359, 5.56763}, {1.13641, 0.80520}}},
        {4,  {{1.75663, 8.15111}, {2.24337, 1.22219}}},
        {5,  {{1.73486, 10.46135}, {3.26514, 1.51797}}},
        {6,  {{1.70730, 12.79758}, {4.29270, 1.85580}}},
        {7,  {{1.68283, 15.13096}, {5.31717, 2.19942}}},
        {8,  {{1.66122, 17.46129}, {6.33878, 2.54326}}},
        {9,  {{1.64171, 19.78991}, {7.35829, 2.88601}}},
        {10, {{1.62380, 22.11938}, {8.37620, 3.22746}}},
        // Period 3
        {11, {{1.48140, 25.82522}, {8.28761, 4.02120}, {1.23098, 0.80897}}},
        {12, {{1.39674, 29.19802}, {8.10904, 4.76791}, {2.49422, 1.08302}}},
        {13, {{1.34503, 32.33363}, {8.12124, 5.42812}, {3.53372, 1.15994}}},
        {14, {{1.28865, 35.65432}, {7.98931, 6.17545}, {4.72204, 1.33797}}},
        {15, {{1.23890, 39.00531}, {7.83125, 6.95265}, {5.92985, 1.52690}}},
        {16, {{1.19478, 42.38177}, {7.66565, 7.75584}, {7.13957, 1.71687}}},
        {17, {{1.15482, 45.79189}, {7.50031, 8.58542}, {8.34487, 1.90546}}},
        {18, {{1.11803, 49.24317}, {7.33917, 9.44200}, {9.54280, 2.09210}}},
        // Period 4
        {19, {{1.09120, 52.59376}, {7.15086, 10.29851}, {9.57061, 2.42121}, {1.18733, 0.67314}}},
        {20, {{1.07196, 55.86008}, {7.01185, 11.11887}, {9.29555, 2.76621}, {2.62063, 0.88692}}},
        {21, {{1.05870, 59.04659}, {6.96404, 11.86718}, {9.97866, 2.93024}, {2.99860, 0.98040}}},
        {22, {{1.04755, 62.22091}, {6.90438, 12.62229}, {10.84355, 3.08264}, {3.20452, 1.05403}}},
        {23, {{1.03828, 65.38117}, {6.83516, 13.38417}, {11.79532, 3.23508}, {3.33124, 1.11609}}},
        {24, {{1.03069, 68.52633}, {6.75998, 14.15132}, {12.79256, 3.38991}, {3.41677, 1.17116}}},
        {25, {{1.02450, 71.65908}, {6.68141, 14.92337}, {13.81149, 3.54730}, {3.48260, 1.22220}}},
        {26, {{1.01960, 74.77846}, {6.60101, 15.69935}, {14.84330, 3.70685}, {3.53609, 1.27026}}},
        {27, {{1.01575, 77.88779}, {6.51976, 16.47941}, {15.88061, 3.86829}, {3.58388, 1.31647}}},
        {28, {{1.01282, 80.98814}, {6.43837, 17.26336}, {16.92012, 4.03115}, {3.62869, 1.36133}}},
        {29, {{1.01839, 83.81831}, {6.47823, 17.85149}, {18.65720, 4.05312}, {2.84618, 1.37570}}},
        {30, {{1.00931, 87.16777}, {6.27682, 18.84319}, {18.99747, 4.35989}, {3.71640, 1.44857}}},
        {31, {{1.00600, 90.34057}, {6.16315, 19.71091}, {19.81836, 4.57852}, {4.01249, 1.29122}}},
        {32, {{0.99467, 93.80965}, {5.91408, 20.85993}, {19.89501, 4.95158}, {5.19624, 1.39361}}},
        {33, {{0.98548, 97.22822}, {5.68319, 22.01684}, {19.83497, 5.33969}, {6.49637, 1.51963}}},
        {34, {{0.97822, 100.60094}, {5.47209, 23.17528}, {19.68845, 5.73803}, {7.86124, 1.65366}}},
        {35, {{0.97231, 103.94730}, {5.27765, 24.33975}, {19.48822, 6.14586}, {9.26182, 1.78869}}},
        {36, {{0.96735, 107.28121}, {5.09646, 25.51581}, {19.25332, 6.56380}, {10.68288, 1.92256}}},
        // Period 5
        {37, {{0.96706, 110.48309}, {4.99899, 26.50212}, {18.99122, 6.92726}, {10.76759, 2.18101}, {1.27514, 0.66954}}},
        {38, {{0.96801, 113.67680}, {4.93897, 27.41353}, {18.80330, 7.25498}, {10.35813, 2.46293}, {2.93159, 0.86625}}},
        {39, {{0.96684, 116.97488}, {4.85128, 28.43242}, {18.94517, 7.56989}, {10.53663, 2.55617}, {3.70008, 0.97332}}},
        {40, {{0.96535, 120.30371}, {4.74742, 29.51436}, {19.08019, 7.90688}, {11.05560, 2.60522}, {4.15144, 1.06342}}},
        {41, {{0.96228, 123.69423}, {4.59559, 30.76696}, {19.29222, 8.28677}, {13.04646, 2.49390}, {3.10346, 1.06692}}},
        {42, {{0.96141, 127.02899}, {4.47623, 31.93860}, {19.27950, 8.67584}, {14.22293, 2.54514}, {3.05992, 1.12354}}},
        {43, {{0.96101, 130.36147}, {4.35465, 33.14068}, {19.21943, 9.08439}, {15.54548, 2.60482}, {2.91943, 1.16411}}},
        {44, {{0.96111, 133.68940}, {4.23256, 34.36928}, {19.12341, 9.51019}, {16.94562, 2.67378}, {2.73729, 1.19329}}},
        {45, {{0.96172, 137.01241}, {4.11107, 35.62210}, {19.00142, 9.95142}, {18.37578, 2.75113}, {2.55000, 1.21476}}},
        {46, {{0.96222, 140.34382}, {3.97975, 36.95092}, {18.97452, 10.39657}, {20.49458, 2.75967}, {1.58893, 1.30000}}},
        {47, {{0.96441, 143.64600}, {3.87219, 38.19648}, {18.71042, 10.87495}, {21.22900, 2.92470}, {2.22399, 1.24459}}},
        {48, {{0.96721, 146.96006}, {3.78203, 39.39447}, {18.45110, 11.34436}, {21.44461, 3.11456}, {3.35504, 1.32681}}},
        {49, {{0.96990, 150.29845}, {3.69397, 40.61549}, {18.22732, 11.81657}, {22.45558, 3.25901}, {3.65322, 1.16914}}},
        {50, {{0.97329, 153.61921}, {3.60552, 41.84383}, {17.87583, 12.33351}, {22.43949, 3.49629}, {5.10586, 1.27096}}},
        {51, {{0.97679, 156.98207}, {3.53402, 43.02243}, {17.55869, 12.83309}, {22.26210, 3.73890}, {6.66839, 1.37644}}},
        {52, {{0.98027, 160.39249}, {3.47671, 44.15853}, {17.26668, 13.31791}, {21.90419, 3.99243}, {8.37215, 1.49131}}},
        {53, {{0.98368, 163.85359}, {3.43054, 45.26231}, {16.99206, 13.79195}, {21.42699, 4.25770}, {10.16674, 1.60670}}},
        {54, {{0.98698, 167.36704}, {3.39333, 46.34048}, {16.72687, 14.25880}, {20.87117, 4.53670}, {12.02165, 1.72006}}},
        // Period 6 (some shells merged; shell count varies)
        {55, {{0.98855, 171.08380}, {3.39137, 47.32865}, {16.81286, 14.56630}, {20.37387, 4.67981}, {12.00664, 1.93716}, {1.42671, 0.64967}}},
        {56, {{0.99731, 175.01327}, {3.94901, 45.62379}, {19.47065, 13.32231}, {25.53046, 3.58081}, {6.05257, 1.01688}}},
        {57, {{0.99490, 178.92585}, {3.78923, 47.36151}, {19.03694, 13.96916}, {25.48253, 3.84663}, {7.69641, 1.12920}}},
        {58, {{0.99773, 182.69729}, {3.79809, 48.23836}, {19.05383, 14.26200}, {26.50828, 3.89345}, {7.64207, 1.14501}}},
        {59, {{1.00049, 186.51613}, {3.80487, 49.13072}, {19.03901, 14.56672}, {27.58038, 3.94881}, {7.57525, 1.15879}}},
        {60, {{1.00319, 190.38490}, {3.81025, 50.03637}, {18.99897, 14.88121}, {28.67528, 4.01194}, {7.51231, 1.17175}}},
        {61, {{0.99995, 194.01705}, {3.24882, 54.04025}, {16.47119, 16.98218}, {23.00328, 5.32410}, {12.91959, 2.51531}, {4.35716, 0.98405}}},
        {62, {{1.00209, 197.99687}, {3.22761, 55.17619}, {16.34706, 17.42152}, {23.80125, 5.42023}, {13.27747, 2.57444}, {4.34451, 0.99393}}},
        {63, {{1.00429, 202.02953}, {3.20862, 56.30909}, {16.21652, 17.86527}, {24.63902, 5.51636}, {13.59806, 2.63146}, {4.33349, 1.00348}}},
        {64, {{1.00654, 206.11731}, {3.19138, 57.44219}, {16.08281, 18.31332}, {25.49620, 5.61277}, {13.89990, 2.68701}, {4.32318, 1.01283}}},
        {65, {{1.00883, 210.26214}, {3.17612, 58.57320}, {15.94486, 18.76554}, {26.38056, 5.70951}, {14.17647, 2.74080}, {4.31316, 1.02181}}},
        {66, {{1.01117, 214.46688}, {3.16274, 59.70253}, {15.80365, 19.22175}, {27.28732, 5.80665}, {14.43201, 2.79301}, {4.30312, 1.03043}}},
        {67, {{1.01354, 218.73341}, {3.15115, 60.83054}, {15.65984, 19.68187}, {28.21315, 5.90421}, {14.67005, 2.84364}, {4.29227, 1.03860}}},
        {68, {{1.01595, 223.06415}, {3.14128, 61.95766}, {15.51409, 20.14581}, {29.15557, 6.00216}, {14.89317, 2.89270}, {4.27995, 1.04625}}},
        {69, {{1.01840, 227.46136}, {3.13302, 63.08447}, {15.36688, 20.61357}, {30.11252, 6.10047}, {15.10393, 2.94016}, {4.26525, 1.05328}}},
        {70, {{1.02088, 231.92848}, {3.12631, 64.21152}, {15.21913, 21.08480}, {31.08374, 6.19886}, {15.30077, 2.98599}, {4.24917, 1.05982}}},
        {71, {{1.02340, 236.46708}, {3.12101, 65.33964}, {15.07067, 21.55992}, {32.06704, 6.29753}, {15.48957, 3.02994}, {4.22831, 1.06533}}},
        {72, {{1.02614, 241.05912}, {3.11380, 66.48612}, {14.91939, 22.05124}, {33.32537, 6.39717}, {14.51375, 3.06086}, {5.10155, 1.16902}}},
        {73, {{1.02929, 245.67915}, {3.10090, 67.66413}, {14.72483, 22.58997}, {34.56189, 6.52150}, {13.79943, 3.03942}, {5.78367, 1.25985}}},
        {74, {{1.03282, 250.33518}, {3.08458, 68.85920}, {14.49437, 23.16956}, {35.56831, 6.67658}, {13.59733, 2.99108}, {6.22259, 1.33543}}},
        {75, {{1.03659, 255.04631}, {3.06638, 70.07019}, {14.24505, 23.77881}, {36.35311, 6.85375}, {13.82288, 2.93796}, {6.47600, 1.39947}}},
        {76, {{1.04065, 259.80929}, {3.04708, 71.28687}, {13.97524, 24.41922}, {36.87935, 7.05641}, {14.58843, 2.89227}, {6.46926, 1.45041}}},
        {77, {{1.04497, 264.63166}, {3.02753, 72.50494}, {13.69135, 25.08677}, {37.19139, 7.28053}, {15.80637, 2.86375}, {6.23840, 1.48910}}},
        {78, {{1.04951, 269.52094}, {3.00836, 73.72166}, {13.39860, 25.77813}, {37.33902, 7.52233}, {17.34748, 2.85595}, {5.85704, 1.51755}}},
        {79, {{1.05572, 274.27910}, {2.96900, 75.00010}, {12.96851, 26.66074}, {37.27090, 7.84449}, {20.99429, 2.78498}, {3.74158, 1.47920}}},
        {80, {{1.05910, 279.52863}, {2.97260, 76.14674}, {12.80192, 27.22274}, {37.30541, 8.04711}, {20.91497, 2.89600}, {4.94600, 1.55349}}},
        {81, {{1.06254, 284.86830}, {2.97535, 77.31377}, {12.63238, 27.79957}, {37.40661, 8.25124}, {23.16204, 2.90463}, {3.76108, 1.25293}}},
        {82, {{1.06627, 290.28925}, {2.98032, 78.45692}, {12.42797, 28.40635}, {36.95686, 8.51126}, {23.07427, 3.10731}, {5.49431, 1.35142}}},
        {83, {{1.06952, 295.89089}, {2.99110, 79.60411}, {12.26809, 28.96419}, {36.60258, 8.74903}, {23.18839, 3.28577}, {6.88031, 1.38803}}},
        {84, {{1.07247, 301.65782}, {3.00695, 80.74997}, {12.13336, 29.48765}, {36.12686, 8.98696}, {22.79251, 3.51606}, {8.86784, 1.48179}}},
        {85, {{1.07504, 307.60992}, {3.02761, 81.90658}, {12.03635, 29.96333}, {35.62798, 9.21085}, {22.13626, 3.77443}, {11.09675, 1.58575}}},
        {86, {{1.07727, 313.74785}, {3.05219, 83.07870}, {11.97392, 30.39458}, {35.09317, 9.42149}, {21.34497, 4.06298}, {13.45848, 1.68982}}},
        // Period 7 (some shells merged; shell count varies)
        {87,  {{1.07866, 320.13828}, {3.07367, 84.37488}, {12.01829, 30.74418}, {35.74839, 9.48671}, {19.04804, 4.19650}, {14.27307, 1.97507}, {1.75989, 0.71651}}},
        {88,  {{1.08004, 326.68617}, {3.09096, 85.73578}, {12.07540, 31.10148}, {37.02057, 9.50368}, {15.73782, 4.21159}, {14.71155, 2.36783}, {4.28366, 0.89370}}},
        {89,  {{1.07846, 333.89969}, {3.17050, 86.83801}, {12.57061, 30.80584}, {40.47386, 9.16822}, {24.82719, 3.11081}, {6.87937, 1.04676}}},
        {90,  {{1.08208, 340.45553}, {3.16504, 88.26458}, {12.37624, 31.47596}, {40.08726, 9.42733}, {24.56010, 3.31023}, {8.72929, 1.15770}}},
        {91,  {{1.08553, 347.22535}, {3.17119, 89.63206}, {12.22718, 32.06968}, {40.31735, 9.62313}, {25.55323, 3.30511}, {8.64552, 1.17701}}},
        {92,  {{1.08939, 354.11371}, {3.17670, 90.99571}, {12.03566, 32.71878}, {40.75159, 9.82514}, {28.06165, 3.18256}, {6.88501, 1.11150}}},
        {93,  {{1.09355, 361.13749}, {3.17538, 92.42691}, {11.81000, 33.44525}, {40.69294, 10.08473}, {29.47989, 3.22145}, {6.74825, 1.12213}}},
        {94,  {{1.09794, 368.32500}, {3.17404, 93.87067}, {11.57027, 34.20518}, {40.57749, 10.36024}, {30.96551, 3.26929}, {6.61475, 1.13111}}},
        {95,  {{1.10254, 375.69056}, {3.17360, 95.32199}, {11.32206, 34.99094}, {40.42154, 10.64838}, {32.48783, 3.32401}, {6.49243, 1.13917}}},
        {96,  {{1.10729, 383.24958}, {3.17377, 96.78830}, {11.07137, 35.79948}, {40.25985, 10.94334}, {33.99604, 3.38149}, {6.39168, 1.14787}}},
        {97,  {{1.11222, 391.01298}, {3.17544, 98.26135}, {10.81808, 36.62729}, {40.07925, 11.24693}, {35.51639, 3.44250}, {6.29863, 1.15589}}},
        {98,  {{1.11730, 398.99467}, {3.17879, 99.74179}, {10.56415, 37.47189}, {39.88651, 11.55781}, {37.04014, 3.50607}, {6.21310, 1.16337}}},
        {99,  {{1.12255, 407.20909}, {3.18387, 101.23136}, {10.31099, 38.33175}, {39.68650, 11.87502}, {38.56155, 3.57146}, {6.13453, 1.17038}}},
        {100, {{1.12794, 415.67191}, {3.19069, 102.73254}, {10.05989, 39.20531}, {39.48362, 12.19755}, {40.07701, 3.63798}, {6.06085, 1.17681}}},
        {101, {{1.13348, 424.39937}, {3.19917, 104.24845}, {9.81173, 40.09171}, {39.28078, 12.52472}, {41.58459, 3.70516}, {5.99026, 1.18259}}},
        {102, {{1.13915, 433.40929}, {3.20918, 105.78307}, {9.56720, 40.99055}, {39.08062, 12.85589}, {43.08210, 3.77261}, {5.92175, 1.18763}}},
        {103, {{1.14496, 442.72080}, {3.22058, 107.34061}, {9.32700, 41.90122}, {38.88515, 13.19043}, {44.56942, 3.83997}, {5.85289, 1.19177}}},
        {104, {{1.15135, 452.26605}, {3.22875, 108.92084}, {9.04304, 42.93588}, {38.34711, 13.60309}, {45.02618, 4.00670}, {7.20356, 1.29377}}},
        {105, {{1.15720, 462.29654}, {3.24386, 110.55376}, {8.82683, 43.83101}, {37.89762, 13.97305}, {45.25449, 4.16409}, {8.62000, 1.39805}}},
        {106, {{1.16267, 472.80654}, {3.26226, 112.26116}, {8.65940, 44.62894}, {37.51680, 14.31134}, {45.31035, 4.31472}, {10.08853, 1.50038}}},
        {107, {{1.16812, 483.75845}, {3.28113, 114.03669}, {8.50793, 45.40411}, {37.13867, 14.64637}, {45.28879, 4.46883}, {11.61536, 1.59154}}},
        {108, {{1.17350, 495.19373}, {3.30025, 115.89209}, {8.37792, 46.14357}, {36.77370, 14.97237}, {45.11349, 4.62719}, {13.26114, 1.68611}}},
        {109, {{1.17885, 507.14435}, {3.31904, 117.83528}, {8.26742, 46.85201}, {36.42110, 15.29007}, {44.80854, 4.78963}, {15.00504, 1.78141}}},
        {110, {{1.18421, 519.64313}, {3.33703, 119.87184}, {8.17332, 47.53691}, {36.07510, 15.60183}, {44.38981, 4.95728}, {16.84052, 1.87653}}},
        {111, {{1.18963, 532.72803}, {3.35389, 122.00711}, {8.09284, 48.20504}, {35.72896, 15.91017}, {43.86996, 5.13167}, {18.76472, 1.97119}}},
        {112, {{1.19513, 546.44374}, {3.36937, 124.24696}, {8.02369, 48.86232}, {35.37591, 16.21753}, {43.26106, 5.31444}, {20.77485, 2.06523}}},
        {113, {{1.19831, 561.46167}, {3.38080, 126.93676}, {8.19964, 49.01826}, {36.01255, 16.19608}, {44.61015, 5.20890}, {19.59855, 1.98048}}},
        {114, {{1.20363, 576.73694}, {3.38927, 129.51002}, {8.18834, 49.61420}, {35.99333, 16.41487}, {45.19984, 5.27324}, {20.02560, 1.97092}}},
        {115, {{1.21246, 591.92054}, {3.39659, 131.76192}, {7.83915, 50.99178}, {34.63427, 17.11127}, {42.94186, 5.74241}, {22.77382, 2.32051}, {2.20185, 1.02716}}},
        {116, {{1.21804, 608.90539}, {3.40088, 134.61027}, {7.83224, 51.62782}, {34.54736, 17.35046}, {42.66989, 5.85605}, {22.24232, 2.46122}, {4.08928, 1.16490}}},
        {117, {{1.22355, 626.95975}, {3.40140, 137.66654}, {7.84820, 52.25136}, {34.58793, 17.55428}, {42.79420, 5.92328}, {20.78333, 2.58696}, {6.36139, 1.30054}}},
        {118, {{1.22909, 646.17926}, {3.39867, 140.93076}, {7.87683, 52.88265}, {34.70741, 17.73814}, {43.32835, 5.95447}, {18.42228, 2.70154}, {9.03736, 1.42977}}},
    };
    auto it = params.find(atomic_num);
    if (it == params.end()) {
        throw PsiException("Atomic Number " + std::to_string(atomic_num) + " unsupported by MBIS", __FILE__, __LINE__);
    }
    return it->second;
}


/* Batched exp() for the MBIS stockholder sweep.
 *
 * The sweep spends nearly all of its time in exp(-r/sigma). A plain loop over std::exp runs one
 * element at a time however wide the machine is, because glibc only vectorises exp through libmvec
 * under -ffast-math, which Psi4 does not build with. Eigen ships its own vectorised exp, so this
 * gets SIMD throughput from a dependency Psi4 already requires rather than from a hand-rolled
 * range reduction that would have to restate the IEEE-754 contract itself.
 *
 * Measured on 96-element batches: 526 Mexp/s against libm's 126 with -march=native, and 138
 * against 130 without it. That second column is the reason this is a library call -- a hand-rolled
 * polynomial kernel is about 2x *slower* than libm on a build that lost its -march flag, and
 * nothing warns you when that happens.
 */
static inline void mbis_exp_batch(const double* __restrict arg, double* __restrict out, int n) {
    Eigen::Map<Eigen::ArrayXd>(out, n) = Eigen::Map<const Eigen::ArrayXd>(arg, n).exp();
}

/* Anderson acceleration for the MBIS stockholder fixed point.
 *
 * The stockholder update is a plain Picard iteration x_{k+1} = g(x_k) on the shell populations and
 * widths, and it converges linearly: 60-80 iterations to 1e-8 on the systems profiled, essentially
 * independent of how fast a single iteration is made. Anderson mixing reuses the last few residuals
 * to extrapolate, which typically removes most of those iterations for a few hundred flops.
 *
 * Anderson is not unconditionally convergent, though, and on an ionic system (NaCl was the case
 * that caught this) an unsafeguarded version oscillates indefinitely and eventually drives a shell
 * population to zero, which turns the whole result into NaN. Four guards, in order of how often
 * they matter:
 *
 *  - N and sigma are strictly positive, and a negative sigma would make exp(-r/sigma) diverge. The
 *    mixing is therefore done on log(N), log(sigma), so every combination is positive by
 *    construction, whatever coefficients the least-squares returns.
 *  - The residual least-squares problem is Tikhonov-regularised. Nearly-parallel residual histories
 *    are the usual source of the huge coefficients that produce a wild extrapolation.
 *  - The extrapolation is rejected outright if the mixing coefficients grow past `kMaxCoefSum`,
 *    if the least-squares is singular, or if the result is not finite. Acceleration is never
 *    allowed to be worse than doing nothing.
 *  - The caller watches the residual it already computes and calls `reset()` when an accelerated
 *    step made it worse, which falls back to plain Picard until the history rebuilds.
 */
class MBISAnderson {
    static constexpr double kRidge = 1.0e-12;    // relative Tikhonov penalty on the mixing coefficients
    static constexpr double kMaxCoefSum = 20.0;  // reject extrapolations that reach this far out

    size_t n_;      // length of the parameter vector
    size_t depth_;  // history length
    // Iterates and residuals, one column each. A ring buffer over fixed columns rather than a
    // deque of vectors: dropping the oldest entry is an index update instead of a memmove of the
    // whole history, and the least-squares matrix below can be filled with block writes instead of
    // hand-rolled column-major arithmetic.
    Eigen::MatrixXd xs_, fs_;
    size_t count_ = 0;  // how many columns are live, capped at depth_
    size_t head_ = 0;   // column holding the oldest live entry

   public:
    MBISAnderson(size_t n, size_t depth)
        : n_(n),
          depth_(depth),
          xs_(Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(depth))),
          fs_(Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(depth))) {}

    /// Column of the history holding entry i, counting 0 as the oldest live entry.
    Eigen::Index slot(size_t i) const { return static_cast<Eigen::Index>((head_ + i) % depth_); }

    /// Given the current iterate x and the Picard image g(x), return the next iterate.
    std::vector<double> next(const std::vector<double>& x, const std::vector<double>& g) {
        const Eigen::Map<const Eigen::VectorXd> xv(x.data(), static_cast<Eigen::Index>(n_));
        const Eigen::Map<const Eigen::VectorXd> gv(g.data(), static_cast<Eigen::Index>(n_));

        // Store the images; the mix is over g's, Anderson type-II.
        const Eigen::Index write = (count_ < depth_) ? slot(count_) : head_;
        xs_.col(write) = gv;
        fs_.col(write) = gv - xv;
        if (count_ < depth_) {
            count_++;
        } else {
            head_ = (head_ + 1) % depth_;  // overwrote the oldest; it is now the newest
        }

        const size_t m = count_;
        if (m < 2) return g;  // nothing to extrapolate from yet

        // Minimise ||sum_i c_i f_i|| subject to sum_i c_i = 1. Eliminate the last coefficient,
        // c_last = 1 - sum_i alpha_i, and solve the resulting least-squares problem directly. This
        // avoids squaring the condition number in normal equations and delegates rank decisions to
        // LAPACK's pivoted QR implementation.
        const int ncols = static_cast<int>(m - 1);
        const int nrows = static_cast<int>(n_) + static_cast<int>(m);
        const int ldb = std::max(nrows, ncols);
        double maxdiag = 0.0;
        for (size_t i = 0; i < m; i++) maxdiag = std::max(maxdiag, fs_.col(slot(i)).squaredNorm());
        if (!(maxdiag > 0.0) || !std::isfinite(maxdiag)) return g;

        // Eigen is column-major, which is what DGELSY wants, so A can be handed over directly.
        // Rows [0, n_) hold f_i - f_last. The next ncols rows penalise each alpha_i, and the final
        // row penalises c_last = 1 - sum_i alpha_i, so the ridge acts on the whole coefficient
        // vector rather than only the part that survived the elimination.
        const Eigen::Index nrows_e = nrows, ncols_e = ncols, n_e = static_cast<Eigen::Index>(n_);
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(nrows_e, ncols_e);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(ldb);

        const auto f_last = fs_.col(slot(m - 1));
        rhs.head(n_e) = -f_last;
        for (Eigen::Index i = 0; i < ncols_e; i++) {
            A.col(i).head(n_e) = fs_.col(slot(static_cast<size_t>(i))) - f_last;
        }
        const double sqrt_ridge = std::sqrt(kRidge * maxdiag);
        A.block(n_e, 0, ncols_e, ncols_e).diagonal().setConstant(sqrt_ridge);
        A.row(n_e + ncols_e).setConstant(sqrt_ridge);
        rhs(n_e + ncols_e) = sqrt_ridge;

        // The problem is tall and extremely thin -- about 550 x 5 -- so threading it buys nothing
        // and costs reproducibility: a threaded QR picks its reduction order at run time, which
        // perturbs the mixing coefficients between otherwise identical runs. Because the iteration
        // stops on a threshold, that perturbation lands in the converged answer (measured: 5.8e-13
        // on the charges of a 109-atom case, against 2.7e-15 with acceleration switched off).
        // Pinning the solve to one thread makes it deterministic and is not slower at this size.
#ifdef USING_LAPACK_MKL
        const int mkl_threads_on_entry = mkl_get_max_threads();
        mkl_set_num_threads(1);
#endif
        // The ridge block is itself full column rank, so the stacked matrix always is too. That
        // makes rank deficiency unreachable and DGELSY's rcond inert; conditioning is handled by
        // kRidge, and a wild extrapolation by the coefficient-sum test below.
        std::vector<int> jpvt(ncols, 0);
        int rank = 0;
        double work_query = 0.0;
        int info = C_DGELSY(nrows, ncols, 1, A.data(), nrows, rhs.data(), ldb, jpvt.data(), 0.0, &rank,
                            &work_query, -1);
        const int lwork = std::max(1, static_cast<int>(work_query));
        std::vector<double> work(lwork);
        if (info == 0 && std::isfinite(work_query)) {
            std::fill(jpvt.begin(), jpvt.end(), 0);
            info = C_DGELSY(nrows, ncols, 1, A.data(), nrows, rhs.data(), ldb, jpvt.data(), 0.0, &rank, work.data(),
                            lwork);
        }
#ifdef USING_LAPACK_MKL
        mkl_set_num_threads(mkl_threads_on_entry);
#endif
        if (info != 0) return g;

        std::vector<double> c(m);
        const double alpha_sum = rhs.head(ncols_e).sum();
        for (int i = 0; i < ncols; i++) c[i] = rhs(i);
        c.back() = 1.0 - alpha_sum;

        // sum_i c_i == 1 by construction, so sum_i |c_i| measures how far outside the history the
        // extrapolation reaches. Large values are the signature of a step about to go wild.
        double coefsum = 0.0;
        for (size_t i = 0; i < m; i++) coefsum += std::fabs(c[i]);
        if (!std::isfinite(coefsum) || coefsum > kMaxCoefSum) return g;

        std::vector<double> out(n_);
        Eigen::Map<Eigen::VectorXd> outv(out.data(), n_e);
        outv.setZero();
        for (size_t i = 0; i < m; i++) outv += c[i] * xs_.col(slot(i));

        if (!outv.allFinite()) return g;  // reject a bad extrapolation wholesale
        return out;
    }

    void reset() {
        count_ = 0;
        head_ = 0;
    }
};

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

    timer_on("MBIS: grid build");
    auto grid = std::make_shared<DFTGrid>(mol, basisset_, mbis_grid_options_int, mbis_grid_options_str, options);
    timer_off("MBIS: grid build");

    if (print_output && debug >= 1) grid->print();

    int nbf = basisset_->nbf();
    int num_atoms = mol->natom();
    size_t total_points = grid->npoints();

    // The pro-atom density array below is one double per (grid point, atom), and it is the one
    // allocation here large enough to fail. Both factors grow with the molecule -- the grid has
    // more points and each point carries more atoms -- so it is quadratic in system size, and the
    // failure would otherwise arrive as a bare std::bad_alloc from deep inside a property
    // calculation that has already paid for an SCF. Check it against the memory the user actually
    // gave Psi4, before the grid density evaluation rather than after, and say what to do about it.
    const size_t mbis_doubles = (static_cast<size_t>(num_atoms) + 7) * total_points;
    const double mbis_gib = static_cast<double>(mbis_doubles) * sizeof(double) / (1024.0 * 1024.0 * 1024.0);
    const double avail_gib = static_cast<double>(Process::environment.get_memory()) / (1024.0 * 1024.0 * 1024.0);
    if (print_output) outfile->Printf("  MBIS Memory: needs %.3f GiB; user supplied %.3f GiB.\n\n", mbis_gib, avail_gib);
    if (mbis_gib > avail_gib) {
        throw PSIEXCEPTION(
            "MBIS needs " + std::to_string(mbis_gib) + " GiB for its grid arrays but only " +
            std::to_string(avail_gib) +
            " GiB is available. The cost is (natom + 7) doubles per grid point, so either raise the "
            "memory given to Psi4, or shrink the grid with the mbis_radial_points and "
            "mbis_spherical_points options.");
    }

    size_t max_points = 0;
    size_t max_nbf = 0;

    std::vector<std::shared_ptr<BlockOPoints>> blocks = grid->blocks();
    for (size_t b = 0; b < blocks.size(); b++) {
        max_points = std::max(max_points, blocks[b]->npoints());
        max_nbf = std::max(max_nbf, blocks[b]->local_nbf());
    }
    timer_on("MBIS: density subset");
    SharedMatrix Da = wfn_->Da_subset("AO");
    SharedMatrix Db;
    auto point_func = (std::shared_ptr<PointFunctions>)(std::make_shared<RKSFunctions>(basisset_, max_points, max_nbf));

    if (same_dens_) {
        Db = Da->clone();
        point_func->set_pointers(Da);
    } else {
        Db = wfn_->Db_subset("AO");
        point_func = (std::shared_ptr<PointFunctions>)(std::make_shared<UKSFunctions>(basisset_, max_points, max_nbf));
        point_func->set_pointers(Da, Db);
    }
    timer_off("MBIS: density subset");

    size_t running_points = 0;

    // Block extents, used for the distance screening in the stockholder sweep below.
    const size_t n_blocks = blocks.size();
    std::vector<size_t> block_offset(n_blocks, 0), block_size(n_blocks, 0);
    std::vector<double> block_cx(n_blocks), block_cy(n_blocks), block_cz(n_blocks), block_R(n_blocks);

    // Coordinates, weights, and rho (molecular electron density) at each grid point
    std::vector<double> x_points(total_points, 0.0);
    std::vector<double> y_points(total_points, 0.0);
    std::vector<double> z_points(total_points, 0.0);
    std::vector<double> weights(total_points, 0.0);
    std::vector<double> rho(total_points, 0.0);

    timer_on("MBIS: grid density");
    for (size_t b = 0; b < blocks.size(); b++) {
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

        block_offset[b] = running_points;
        block_size[b] = num_points;
        const Vector3 block_center = block->center();
        block_cx[b] = block_center[0];
        block_cy[b] = block_center[1];
        block_cz[b] = block_center[2];
        block_R[b] = block->radius();
        running_points += num_points;
    }
    timer_off("MBIS: grid density");

    // Electron count via numerical interagration
    double grid_electrons = 0.0;
    for (size_t p = 0; p < total_points; p++) {
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

    // Nuclear coordinates, kept flat so the hot loop below does not touch Molecule.
    std::vector<double> atom_x(num_atoms), atom_y(num_atoms), atom_z(num_atoms);
    for (int atom = 0; atom < num_atoms; atom++) {
        Vector3 R = mol->xyz(atom);
        atom_x[atom] = R[0];
        atom_y[atom] = R[1];
        atom_z[atom] = R[2];
    }

    // => Setup Proatom Basis Functions <= //

    // mA is the number of proatom shells per atom, read directly from the parameter table.
    std::vector<int> mA(num_atoms);
    for (int atom = 0; atom < num_atoms; atom++) {

        // MBIS needs an all-electron density, so an atom whose nuclear charge does not match its
        // element is out of scope. Two quite different things land here, and they need different
        // messages: an ECP replaces some of the electrons (0 < Z < true Z), whereas a ghost atom
        // has no nucleus and no electrons at all (Z == 0), and telling someone doing a
        // counterpoise correction that their ghost is an ECP sends them somewhere useless.
        int n_valence_electrons = static_cast<int>(mol->Z(atom));
        int true_atomic_num = static_cast<int>(mol->true_atomic_number(atom));
        if (n_valence_electrons == 0) {
            throw PSIEXCEPTION("MBIS does not support ghost atoms. Atom " + std::to_string(atom + 1) + " (" +
                               mol->symbol(atom) +
                               ") carries basis functions but no electrons, so it has no pro-atom to fit. Run MBIS on "
                               "the molecule without ghosts.");
        }
        if (n_valence_electrons != true_atomic_num) {
            throw PSIEXCEPTION("MBIS incompatible with ECP. ECP detected on atom " + std::to_string(atom + 1) + " (" +
                               mol->symbol(atom) + "). Use all-electron basis or reconstruct density with denspart.");
        }

        int atomic_num = static_cast<int>(mol->Z(atom));
        mA[atom] = static_cast<int>(get_mbis_params(atomic_num).size());
    }

    // Shells are flattened to a single running index so the inner loop is a flat sweep.
    // shell_off[a] is where atom a's shells begin; total_shells = sum_a mA[a].
    std::vector<int> shell_off(num_atoms + 1, 0);
    for (int atom = 0; atom < num_atoms; atom++) shell_off[atom + 1] = shell_off[atom] + mA[atom];
    const int total_shells = shell_off[num_atoms];

    // Atomic shell populations (N) and widths (S), from equations 18 and 19 in Verstraelen et al.
    std::vector<double> Nf(total_shells, 0.0), Sf(total_shells, 0.0);
    std::vector<double> Nf_next(total_shells, 0.0), Sf_next(total_shells, 0.0);

    for (int atom = 0; atom < num_atoms; atom++) {
        auto shells = get_mbis_params(static_cast<int>(mol->Z(atom)));
        for (int m = 0; m < mA[atom]; m++) {
            Nf[shell_off[atom] + m] = std::get<0>(shells[m]);
            Sf[shell_off[atom] + m] = 1.0 / std::get<1>(shells[m]);
            if (print_output && debug >= 1) {
                outfile->Printf("  INITIAL ATOM %d, SHELL %d, POP %8.5f, WIDTH %8.5f\n", atom + 1, m + 1,
                                Nf[shell_off[atom] + m], Sf[shell_off[atom] + m]);
            }
        }
    }

    // Pro-atom densities, point-major (point * num_atoms + atom): the fused sweep below writes all
    // atoms of one point together, so this layout keeps those writes contiguous. Only the previous
    // iteration's values are needed (for the convergence test), so this is the single large array
    // MBIS now holds -- distances, displacements and the per-atom density are all recomputed.
    // (Its size was checked against the available memory above, before the grid density pass.)
    timer_on("MBIS: alloc");
    std::vector<double> rho_a_0_points(static_cast<size_t>(num_atoms) * total_points, 0.0);
    timer_off("MBIS: alloc");
    std::vector<double> rho_0_points(total_points, 0.0);

    // rho_ai_0 = N * exp(-r/S) / (8 pi S^3). Hoisting N/(8 pi S^3) and -1/S out of the point loop
    // removes a pow() and a division from every one of the natom * mA * npoints evaluations.
    std::vector<double> shell_coef(total_shells), shell_rate(total_shells);
    auto refresh_shell_constants = [&](const std::vector<double>& N, const std::vector<double>& S) {
        for (int s = 0; s < total_shells; s++) {
            shell_coef[s] = N[s] / (8.0 * M_PI * S[s] * S[s] * S[s]);
            shell_rate[s] = -1.0 / S[s];
        }
    };

    // Distance screening. The pro-atom density falls off as exp(-r/sigma) with sigma ~ 0.3-1 bohr,
    // so beyond a few tens of bohr an atom's contribution to a grid point is far below any
    // meaningful threshold -- yet the unscreened loop still visits every (atom, point) pair, which
    // is what makes MBIS quadratic in system size. For each block we keep only the atoms whose
    // largest possible contribution anywhere in that block exceeds MBIS_SCREENING_THRESHOLD.
    //
    // This is bounded, not heuristic: a discarded atom contributes less than the threshold to
    // rho_0, and its own population/width sums pick up an error below (w * rho * threshold / rho_0)
    // at that point. Where rho_0 is itself near the threshold the molecular density rho is
    // negligible too, so the absolute error stays bounded by the threshold either way.
    const double screen_eps = options.get_double("MBIS_SCREENING_THRESHOLD");
    const bool do_screen = screen_eps > 0.0;
    std::vector<int> block_atoms;  // flat concatenation of per-block atom lists, ascending within a block
    std::vector<size_t> block_atoms_off(n_blocks + 1, 0);
    std::vector<double> atom_cut(num_atoms, 0.0);

    // rho_a_0_points is kept valid for *every* (point, atom), with a hard zero wherever the atom is
    // screened out, so the post-processing below needs no knowledge of the screening. Maintaining
    // that by memset-ing the array each sweep would cost more traffic than the screened sweep
    // itself, so instead the zeros are written only when an atom actually drops out of a block's
    // list -- which happens in the first iteration or two and then essentially never. block_drop
    // carries those (block, atom) pairs from rebuild_screening to the next sweep.
    std::vector<int> block_drop;
    std::vector<size_t> block_drop_off(n_blocks + 1, 0);
    std::vector<int> prev_atoms;
    std::vector<size_t> prev_atoms_off(n_blocks + 1, 0);
    bool first_screening = true;

    auto rebuild_screening = [&]() {
        for (int atom = 0; atom < num_atoms; atom++) {
            double cut = std::numeric_limits<double>::infinity();
            if (do_screen) {
                cut = 0.0;
                const double shell_eps = screen_eps / (shell_off[atom + 1] - shell_off[atom]);
                for (int s = shell_off[atom]; s < shell_off[atom + 1]; s++) {
                    // Requiring every shell to be below eps / nshell guarantees that their sum,
                    // the pro-atom density documented by MBIS_SCREENING_THRESHOLD, is below eps.
                    if (shell_coef[s] > shell_eps) {
                        const double S = -1.0 / shell_rate[s];
                        cut = std::max(cut, S * std::log(shell_coef[s] / shell_eps));
                    }
                }
            }
            atom_cut[atom] = cut;
        }

        prev_atoms.swap(block_atoms);
        prev_atoms_off.swap(block_atoms_off);
        block_atoms.clear();
        block_drop.clear();
        for (size_t b = 0; b < n_blocks; b++) {
            block_atoms_off[b] = block_atoms.size();
            block_drop_off[b] = block_drop.size();
            for (int atom = 0; atom < num_atoms; atom++) {
                const double dx = block_cx[b] - atom_x[atom], dy = block_cy[b] - atom_y[atom],
                             dz = block_cz[b] - atom_z[atom];
                const double d = std::sqrt(dx * dx + dy * dy + dz * dz) - block_R[b];
                if (d < atom_cut[atom]) block_atoms.push_back(atom);
            }
            // Atoms present last time but not now: their stale densities must be zeroed. Both lists
            // are ascending, so this is a merge. The sweep accounts for the small transition from
            // the stale value to zero in the convergence delta.
            if (!first_screening) {
                size_t i = prev_atoms_off[b], j = block_atoms_off[b];
                const size_t iend = prev_atoms_off[b + 1], jend = block_atoms.size();
                while (i < iend) {
                    while (j < jend && block_atoms[j] < prev_atoms[i]) j++;
                    if (j == jend || block_atoms[j] != prev_atoms[i]) block_drop.push_back(prev_atoms[i]);
                    i++;
                }
            }
        }
        block_atoms_off[n_blocks] = block_atoms.size();
        block_drop_off[n_blocks] = block_drop.size();
        first_screening = false;
    };

    // One sweep over the grid at fixed parameters. It simultaneously
    //   (a) rebuilds the pro-atom/promolecule densities for these parameters,
    //   (b) accumulates the convergence delta against the previous sweep's densities, and
    //   (c) accumulates the population/width sums that define the *next* parameters.
    // The original code did (a)+(b) and (c) in two separate sweeps one iteration apart, evaluating
    // the identical set of exponentials twice; fusing them halves the exponential count.
    auto sweep = [&](std::vector<double>& sum_n, std::vector<double>& sum_s,
                     std::vector<double>& delta_atom, bool accumulate_delta) {
        std::fill(sum_n.begin(), sum_n.end(), 0.0);
        std::fill(sum_s.begin(), sum_s.end(), 0.0);
        std::fill(delta_atom.begin(), delta_atom.end(), 0.0);

        // Fixed chunks preserve dynamic load balancing without making the summation order depend on
        // which thread finishes first. The chunk count also caps the extra reduction storage.
        constexpr size_t kMaxReductionChunks = 64;
        const size_t n_chunks = std::min(n_blocks, kMaxReductionChunks);
        std::vector<double> chunk_n(n_chunks * total_shells, 0.0), chunk_s(n_chunks * total_shells, 0.0);
        std::vector<double> chunk_d(n_chunks * num_atoms, 0.0);

#pragma omp parallel
        {
            std::vector<double> atom_r(num_atoms);
            // Per-block flattened shell tables, so the exponentials of one grid point form a
            // single contiguous batch instead of num_atoms batches of one to five.
            std::vector<double> bcoef, brate, barg, bval;
            std::vector<int> bglob, bstart;
            // Accumulate into thread-local buffers and flush into the chunk slot at the end. The
            // arithmetic order is still fixed by the chunk, so the result stays deterministic, but
            // the hot accumulator stays in this thread's cache instead of living in a shared array.
            std::vector<double> tls_n(total_shells), tls_s(total_shells), tls_d(num_atoms);

            // Chunks vary in surviving-atom count, so hand them out dynamically. Each chunk covers a
            // fixed contiguous range of blocks and owns its reduction buffers.
#pragma omp for schedule(dynamic, 1)
            for (size_t chunk = 0; chunk < n_chunks; chunk++) {
                double* loc_n = tls_n.data();
                double* loc_s = tls_s.data();
                double* loc_d = tls_d.data();
                std::fill(tls_n.begin(), tls_n.end(), 0.0);
                std::fill(tls_s.begin(), tls_s.end(), 0.0);
                std::fill(tls_d.begin(), tls_d.end(), 0.0);
                const size_t block_begin = n_blocks * chunk / n_chunks;
                const size_t block_end = n_blocks * (chunk + 1) / n_chunks;
                for (size_t b = block_begin; b < block_end; b++) {
                    const size_t off = block_offset[b], np = block_size[b];
                    const int* atoms = block_atoms.data() + block_atoms_off[b];
                    const int na = static_cast<int>(block_atoms_off[b + 1] - block_atoms_off[b]);

                    // Atoms that left this block's list since the last sweep still hold that sweep's
                    // densities; clear them so the (point, atom) entries stay uniformly valid.
                    for (size_t k = block_drop_off[b]; k < block_drop_off[b + 1]; k++) {
                        const int atom = block_drop[k];
                        for (size_t p = off; p < off + np; p++) {
                            double& stale_density = rho_a_0_points[p * num_atoms + atom];
                            if (accumulate_delta) loc_d[atom] += weights[p] * stale_density * stale_density;
                            stale_density = 0.0;
                        }
                    }
                    if (na == 0) continue;

                    bstart.assign(na + 1, 0);
                    bcoef.clear();
                    brate.clear();
                    bglob.clear();
                    for (int ia = 0; ia < na; ia++) {
                        const int atom = atoms[ia];
                        for (int s = shell_off[atom]; s < shell_off[atom + 1]; s++) {
                            bcoef.push_back(shell_coef[s]);
                            brate.push_back(shell_rate[s]);
                            bglob.push_back(s);
                        }
                        bstart[ia + 1] = static_cast<int>(bcoef.size());
                    }
                    const int ns = static_cast<int>(bcoef.size());
                    barg.resize(ns);
                    bval.resize(ns);

                    for (size_t point = off; point < off + np; point++) {
                        const double xp = x_points[point], yp = y_points[point], zp = z_points[point];

                        for (int ia = 0; ia < na; ia++) {
                            const int atom = atoms[ia];
                            const double dx = xp - atom_x[atom], dy = yp - atom_y[atom], dz = zp - atom_z[atom];
                            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                            atom_r[atom] = r;
                            for (int k = bstart[ia]; k < bstart[ia + 1]; k++) barg[k] = brate[k] * r;
                        }

                        mbis_exp_batch(barg.data(), bval.data(), ns);

                        double rho_0 = 0.0;
#pragma omp simd reduction(+ : rho_0)
                        for (int k = 0; k < ns; k++) {
                            bval[k] *= bcoef[k];
                            rho_0 += bval[k];
                        }

                        const double w = weights[point];
                        const double scale = w * rho[point] / rho_0;
                        double* rho_a_0_here = &rho_a_0_points[point * num_atoms];

                        for (int ia = 0; ia < na; ia++) {
                            const int atom = atoms[ia];
                            const double r = atom_r[atom];
                            double rho_a_0 = 0.0;
                            for (int k = bstart[ia]; k < bstart[ia + 1]; k++) {
                                const double v = bval[k];
                                loc_n[bglob[k]] += scale * v;
                                loc_s[bglob[k]] += scale * r * v;
                                rho_a_0 += v;
                            }
                            if (accumulate_delta) {
                                const double d = rho_a_0 - rho_a_0_here[atom];
                                loc_d[atom] += w * d * d;
                            }
                            rho_a_0_here[atom] = rho_a_0;
                        }
                        rho_0_points[point] = rho_0;
                    }
                }
                std::copy(tls_n.begin(), tls_n.end(), chunk_n.begin() + chunk * total_shells);
                std::copy(tls_s.begin(), tls_s.end(), chunk_s.begin() + chunk * total_shells);
                std::copy(tls_d.begin(), tls_d.end(), chunk_d.begin() + chunk * num_atoms);
            }
        }

        for (size_t chunk = 0; chunk < n_chunks; chunk++) {
            const double* loc_n = chunk_n.data() + chunk * total_shells;
            const double* loc_s = chunk_s.data() + chunk * total_shells;
            const double* loc_d = chunk_d.data() + chunk * num_atoms;
            for (int s = 0; s < total_shells; s++) {
                sum_n[s] += loc_n[s];
                sum_s[s] += loc_s[s];
            }
            for (int atom = 0; atom < num_atoms; atom++) delta_atom[atom] += loc_d[atom];
        }

        for (int s = 0; s < total_shells; s++) {
            Nf_next[s] = sum_n[s];
            Sf_next[s] = sum_s[s] / (3.0 * Nf_next[s]);
            // A shell that has lost all of its population makes the width 0/0. There is no way to
            // continue from that, and silently carrying the NaN through to the printed charges is
            // the worst possible outcome, so say what happened.
            if (!(Nf_next[s] > 0.0) || !std::isfinite(Sf_next[s])) {
                throw PSIEXCEPTION(
                    "MBIS: the stockholder iteration reached a shell with zero population, which "
                    "leaves its width undefined. This usually means the pro-atom parameters have "
                    "diverged; try MBIS_ANDERSON false, or a finer MBIS grid.");
            }
        }
    };

    std::vector<double> sum_n(total_shells), sum_s(total_shells), delta_rho_atoms_0(num_atoms);

    // Priming sweep: builds the pro-densities for the initial guess and produces the first update.
    // Its delta is meaningless (there is no previous density), so it is not accumulated.
    timer_on("MBIS: stockholder");
    refresh_shell_constants(Nf, Sf);
    rebuild_screening();
    sweep(sum_n, sum_s, delta_rho_atoms_0, false);
    Nf = Nf_next;
    Sf = Sf_next;

    // => Main Stockholder Loop <= //

    int iter = 1;
    bool is_converged = false;
    double delta_rho_max_0 = 0.0;

    // Anderson mixing is a pure convergence accelerator: switching it off must give the same
    // fixed point, only more slowly. Keep that escape hatch genuinely reachable.
    const bool anderson_enabled = options.get_bool("MBIS_ANDERSON");
    bool accelerate = anderson_enabled;
    MBISAnderson anderson(2 * static_cast<size_t>(total_shells), 6);
    std::vector<double> log_x(2 * total_shells), log_g(2 * total_shells);

    // Anderson steps are not monotone even when they are working, so a single uptick in the
    // residual is normal. A large one is not: it means the extrapolation is fighting the
    // iteration rather than helping it, and on ionic systems it can keep that up indefinitely.
    // Dropping the history at that point falls back to plain Picard until it rebuilds.
    constexpr double kResidualBlowup = 2.0;
    double prev_delta = std::numeric_limits<double>::infinity();
    bool last_step_accelerated = false;

    if (print_output && debug >= 1) outfile->Printf("                     Delta D\n");
    while (iter < max_iter) {
        refresh_shell_constants(Nf, Sf);
        rebuild_screening();
        sweep(sum_n, sum_s, delta_rho_atoms_0, true);

        delta_rho_max_0 = 0.0;
        for (int atom = 0; atom < num_atoms; atom++) {
            const double d = std::sqrt(delta_rho_atoms_0[atom]);
            if (d > delta_rho_max_0) delta_rho_max_0 = d;
        }

        if (print_output && debug >= 1) outfile->Printf("   @MBIS iter %3d:  %.3e\n", iter, delta_rho_max_0);

        if (last_step_accelerated && delta_rho_max_0 > kResidualBlowup * prev_delta) {
            anderson.reset();
        }
        prev_delta = delta_rho_max_0;

        // The convergence test measures the density change between successive iterates. For a plain
        // Picard step that change *is* the fixed-point residual, but an Anderson step is an
        // extrapolation, so two accelerated iterates can be close together without being at the
        // fixed point -- accepting that would silently loosen MBIS_D_CONVERGENCE (measured: errors
        // grew from 1e-11 to 1.6e-5). So the first time the accelerated iteration looks converged,
        // acceleration is switched off and the test must be passed again on an unaccelerated step,
        // where step and residual coincide. That costs one or two extra sweeps and restores the
        // original convergence semantics exactly.
        if (delta_rho_max_0 < conv) {
            if (accelerate) {
                accelerate = false;
                anderson.reset();
            } else {
                if (print_output && debug >= 1) outfile->Printf("  MBIS Atomic Density Converged\n\n");
                is_converged = true;
                break;
            }
        } else if (anderson_enabled && !accelerate && delta_rho_max_0 > 100.0 * conv) {
            // The switch-off was premature; the iterate was not as close as the step suggested.
            accelerate = true;
        }

        // rho_a_0_points/rho_0_points now describe Nf/Sf; step the parameters forward, with the
        // Picard step handed to Anderson (log space, so positivity of N and sigma is automatic).
        if (accelerate) {
            for (int s = 0; s < total_shells; s++) {
                log_x[s] = std::log(Nf[s]);
                log_x[total_shells + s] = std::log(Sf[s]);
                log_g[s] = std::log(Nf_next[s]);
                log_g[total_shells + s] = std::log(Sf_next[s]);
            }
            const std::vector<double> log_new = anderson.next(log_x, log_g);
            for (int s = 0; s < total_shells; s++) {
                Nf[s] = std::exp(log_new[s]);
                Sf[s] = std::exp(log_new[total_shells + s]);
            }
            last_step_accelerated = true;
        } else {
            Nf = Nf_next;
            Sf = Sf_next;
            last_step_accelerated = false;
        }
        iter += 1;
    }

    timer_off("MBIS: stockholder");

    // Repack into the per-atom form the printing and post-processing below expect.
    std::vector<std::vector<double>> Nai(num_atoms), Sai(num_atoms);
    for (int atom = 0; atom < num_atoms; atom++) {
        Nai[atom].assign(Nf.begin() + shell_off[atom], Nf.begin() + shell_off[atom + 1]);
        Sai[atom].assign(Sf.begin() + shell_off[atom], Sf.begin() + shell_off[atom + 1]);
    }

    if (!is_converged) throw ConvergenceError<int>("MBIS", max_iter, conv, delta_rho_max_0, __FILE__, __LINE__);

    // => Post-Processing <= //

    // Multipoles and radial moments are integrals of the same atomic density (Equation 5 in
    // Verstraelen et al.),
    //     rho_a(p) = rho(p) * rho_a_0(p) / rho_0(p),
    // against different polynomials in the point-to-nucleus displacement, so they are accumulated
    // together in one screened pass over the grid. The previous code made num_atoms separate
    // passes for the multipoles and num_atoms * (max_power - 1) more for the radial moments, each
    // one streaming the whole (point, atom) pro-density array -- 12.9 GB of traffic at 48 atoms --
    // and calling pow() once per point. Displacements are recomputed rather than stored: at 109
    // atoms they would be another ~4 GB.
    //
    // Screened-out atoms hold an exact zero in rho_a_0_points, so restricting the inner loop to
    // each block's atom list drops only terms that are identically zero.

    // Non-redundant cartesian quadrupole and octupole components
    static const int qpole_inds[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
    static const int opole_inds[10][3] = {{0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {0, 1, 1}, {0, 1, 2},
                                          {0, 2, 2}, {1, 1, 1}, {1, 1, 2}, {1, 2, 2}, {2, 2, 2}};

    // Non-redundant cartesian multipoles, as given in the convention in the HORTON software
    /* Multipole conventions given in "The Theory of Intermolecular Forces (2nd Edition)" by
       Anthony Stone are a possible future addition. */
    auto mpole = std::make_shared<Matrix>("MBIS Charges: (a.u.)", num_atoms, 1);
    auto dpole = std::make_shared<Matrix>("MBIS Dipoles: (a.u.)", num_atoms, 3);
    auto qpole = std::make_shared<Matrix>("MBIS Quadrupoles: (a.u.)", num_atoms, 6);
    auto opole = std::make_shared<Matrix>("MBIS Octupoles: (a.u.)", num_atoms, 10);

    timer_on("MBIS: multipoles");

    // <r^n> for n = 2 .. rmom_top. MAX_RADIAL_MOMENT can ask for more, but <r^3> (the volume used
    // for the free-atom volume ratios) must always be available.
    const int rmom_top = std::max(4, options.get_int("MAX_RADIAL_MOMENT"));
    const int n_rmom = rmom_top - 1;
    const int n_acc = 20 + n_rmom;  // 1 monopole + 3 dipole + 6 quadrupole + 10 octupole + moments

    std::vector<double> acc(static_cast<size_t>(num_atoms) * n_acc, 0.0);
    constexpr size_t kMaxReductionChunks = 64;
    const size_t n_chunks = std::min(n_blocks, kMaxReductionChunks);
    std::vector<double> chunk_acc(n_chunks * acc.size(), 0.0);

#pragma omp parallel
    {
#pragma omp for schedule(dynamic, 1)
        for (size_t chunk = 0; chunk < n_chunks; chunk++) {
            double* loc = chunk_acc.data() + chunk * acc.size();
            const size_t block_begin = n_blocks * chunk / n_chunks;
            const size_t block_end = n_blocks * (chunk + 1) / n_chunks;
            for (size_t b = block_begin; b < block_end; b++) {
                const size_t off = block_offset[b], np = block_size[b];
                const int* atoms = block_atoms.data() + block_atoms_off[b];
                const int na = static_cast<int>(block_atoms_off[b + 1] - block_atoms_off[b]);
                if (na == 0) continue;  // no atom reaches this block; rho_0 was never formed there

                for (size_t p = off; p < off + np; p++) {
                    const double xp = x_points[p], yp = y_points[p], zp = z_points[p];
                    const double wrho = weights[p] * rho[p] / rho_0_points[p];
                    const double* rho_a_0_here = &rho_a_0_points[p * num_atoms];

                    for (int ia = 0; ia < na; ia++) {
                        const int a = atoms[ia];
                        const double d[3] = {xp - atom_x[a], yp - atom_y[a], zp - atom_z[a]};
                        // c = -w * rho_a; the multipoles carry the electron's negative charge, the
                        // radial moments do not.
                        const double c = -wrho * rho_a_0_here[a];
                        double* la = &loc[static_cast<size_t>(a) * n_acc];

                        la[0] += c;
                        for (int i = 0; i < 3; i++) la[1 + i] += c * d[i];
                        for (int q = 0; q < 6; q++) la[4 + q] += c * d[qpole_inds[q][0]] * d[qpole_inds[q][1]];
                        for (int o = 0; o < 10; o++)
                            la[10 + o] += c * d[opole_inds[o][0]] * d[opole_inds[o][1]] * d[opole_inds[o][2]];

                        const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                        const double r = std::sqrt(r2);
                        double rn = r2;  // n = 2
                        la[20] -= c * rn;
                        for (int n = 3; n <= rmom_top; n++) {
                            rn *= r;
                            la[20 + n - 2] -= c * rn;
                        }
                    }
                }
            }
        }
    }
    for (size_t chunk = 0; chunk < n_chunks; chunk++) {
        const double* loc = chunk_acc.data() + chunk * acc.size();
        for (size_t i = 0; i < acc.size(); i++) acc[i] += loc[i];
    }

    std::vector<SharedMatrix> rmoms;
    for (int n = 2; n <= rmom_top; n++) {
        rmoms.push_back(std::make_shared<Matrix>("ATOMIC RADIAL MOMENTS <R^" + std::to_string(n) + ">", num_atoms, 1));
    }
    for (int a = 0; a < num_atoms; a++) {
        const double* la = &acc[static_cast<size_t>(a) * n_acc];
        mpole->set(a, 0, mol->Z(a) + la[0]);
        for (int i = 0; i < 3; i++) dpole->set(a, i, la[1 + i]);
        for (int q = 0; q < 6; q++) qpole->set(a, q, la[4 + q]);
        for (int o = 0; o < 10; o++) opole->set(a, o, la[10 + o]);
        for (int n = 2; n <= rmom_top; n++) rmoms[n - 2]->set(a, 0, la[20 + n - 2]);
    }

    timer_off("MBIS: multipoles");

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

    for (int n = 2; n <= max_power; n++) {
        std::stringstream sstream;
        sstream << "MBIS RADIAL MOMENTS <R^" << n << ">";

        auto var_name = sstream.str();
        wfn_->set_array_variable(var_name, rmoms[n - 2]);
    }

    auto valence_widths = std::make_shared<Matrix>("MBIS Valence Widths", num_atoms, 1);
    auto valence_charges = std::make_shared<Matrix>("MBIS Valence Charges", num_atoms, 1);
    for (int atom = 0; atom < num_atoms; atom++) {
        valence_widths->set(atom, 0, Sai[atom][mA[atom] - 1]);
        valence_charges->set(atom, 0, -Nai[atom][mA[atom] - 1]);
    }
    wfn_->set_array_variable("MBIS VALENCE WIDTHS", valence_widths);
    wfn_->set_array_variable("MBIS VALENCE CHARGES", valence_charges);

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

        outfile->Printf("\n\n  MBIS Valence Charges: (a.u.)\n");
        outfile->Printf("   Center  Symbol  Z     Charge\n");

        for (int a = 0; a < num_atoms; a++) {
            outfile->Printf("  %5d      %2s %4d   %9.6f\n", a + 1, mol->label(a).c_str(), static_cast<int>(mol->Z(a)),
                            valence_charges->get(a, 0));
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
