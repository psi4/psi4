/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

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
#include "psi4/libmints/wavefunction.h"
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

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <utility>
#include <fstream>
#include <regex>
#include <tuple>

namespace psi {

Prop::Prop(std::shared_ptr<Wavefunction> wfn) : wfn_(wfn)
{
    if (wfn_.get() == NULL)
        throw PSIEXCEPTION("Prop: Wavefunction is null");
    set_wavefunction(wfn_);
}
Prop::~Prop()
{
}
void Prop::common_init()
{
    title_ = "";
    print_ = 1;
    debug_ = 0;
    tasks_.clear();
    set_wavefunction(wfn_);
}
void Prop::set_wavefunction(std::shared_ptr<Wavefunction> wfn)
{
    wfn_ = wfn;

    basisset_ = wfn_->basisset();
    same_orbs_ = wfn_->same_a_b_orbs();
    same_dens_ = wfn_->same_a_b_dens();

    integral_ = std::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));

    std::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral_));
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
void Prop::set_restricted(bool restricted)
{
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
void Prop::set_epsilon_a(SharedVector epsilon_a)
{
    epsilon_a_ = epsilon_a;
    if (same_orbs_) {
        epsilon_b_ = epsilon_a_;
    }
}
void Prop::set_epsilon_b(SharedVector epsilon_b)
{
    if (same_orbs_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting epsilon_b makes no sense");
    epsilon_b_ = epsilon_b;
}
void Prop::set_Ca(SharedMatrix C)
{
    Ca_so_ = C;
    if (same_orbs_) {
        Ca_so_ = Ca_so_;
    }
}
void Prop::set_Cb(SharedMatrix C)
{
    if (same_orbs_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting Cb makes no sense");

    Cb_so_ = C;
}
void Prop::set_Da_ao(SharedMatrix D, int symm)
{
    Da_so_ = SharedMatrix(new Matrix("Da_so", Ca_so_->rowspi(),Ca_so_->rowspi(),symm));

    double* temp = new double[AO2USO_->max_ncol() * AO2USO_->max_nrow()];
    for (int h = 0; h < AO2USO_->nirrep(); ++h) {
        int nao = AO2USO_->rowspi()[0];
        int nsol = AO2USO_->colspi()[h];
        int nsor = AO2USO_->colspi()[h^symm];

        if (!nsol || !nsor) continue;

        double** Ulp = AO2USO_->pointer(h);
        double** Urp = AO2USO_->pointer(h^symm);
        double** DAOp = D->pointer();
        double** DSOp = Da_so_->pointer(h);
        C_DGEMM('N','N',nao,nsor,nao,1.0,DAOp[0],nao,Urp[0],nsor,0.0,temp,nsor);
        C_DGEMM('T','N',nsol,nsor,nao,1.0,Ulp[0],nsol,temp,nsor,0.0,DSOp[0],nsor);
    }
    delete[] temp;

    if (same_dens_) {
        Db_so_ = Da_so_;
    }
}
void Prop::set_Db_ao(SharedMatrix D, int symm)
{
    if (same_dens_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting Db makes no sense");

    Db_so_ = SharedMatrix(new Matrix("Db_so", Cb_so_->rowspi(),Cb_so_->rowspi(),symm));

    double* temp = new double[AO2USO_->max_ncol() * AO2USO_->max_nrow()];
    for (int h = 0; h < AO2USO_->nirrep(); ++h) {
        int nao = AO2USO_->rowspi()[0];
        int nsol = AO2USO_->colspi()[h];
        int nsor = AO2USO_->colspi()[h^symm];

        if (!nsol || !nsor) continue;

        double** Ulp = AO2USO_->pointer(h);
        double** Urp = AO2USO_->pointer(h^symm);
        double** DAOp = D->pointer();
        double** DSOp = Db_so_->pointer(h);
        C_DGEMM('N','N',nao,nsor,nao,1.0,DAOp[0],nao,Urp[0],nsor,0.0,temp,nsor);
        C_DGEMM('T','N',nsol,nsor,nao,1.0,Ulp[0],nsol,temp,nsor,0.0,DSOp[0],nsor);
    }
    delete[] temp;
}
void Prop::set_Da_so(SharedMatrix D)
{
    Da_so_ = D;
    if (same_dens_) {
        Db_so_ = Da_so_;
    }
}
void Prop::set_Db_so(SharedMatrix D)
{
    if (same_dens_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting Db makes no sense");

    Db_so_ = D;
}
void Prop::set_Da_mo(SharedMatrix D)
{
    Da_so_ = SharedMatrix(new Matrix("Da_so", Ca_so_->rowspi(), Ca_so_->rowspi(), D->symmetry()));

    int symm = D->symmetry();
    int nirrep = D->nirrep();

    double* temp = new double[Ca_so_->max_ncol() * Ca_so_->max_nrow()];
    for (int h = 0; h < nirrep; h++) {
        int nmol = Ca_so_->colspi()[h];
        int nmor = Ca_so_->colspi()[h^symm];
        int nsol = Ca_so_->rowspi()[h];
        int nsor = Ca_so_->rowspi()[h^symm];
        if (!nmol || !nmor || !nsol || !nsor) continue;
        double** Clp = Ca_so_->pointer(h);
        double** Crp = Ca_so_->pointer(h^symm);
        double** Dmop = D->pointer(h^symm);
        double** Dsop = Da_so_->pointer(h^symm);
        C_DGEMM('N','T',nmol,nsor,nmor,1.0,Dmop[0],nmor,Crp[0],nmor,0.0,temp,nsor);
        C_DGEMM('N','N',nsol,nsor,nmol,1.0,Clp[0],nmol,temp,nsor,0.0,Dsop[0],nsor);
    }
    delete[] temp;

    if (same_dens_) {
        Db_so_ = Da_so_;
    }
}
void Prop::set_Db_mo(SharedMatrix D)
{
    if (same_dens_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting Db makes no sense");

    Db_so_ = SharedMatrix(new Matrix("Db_so", Cb_so_->rowspi(), Cb_so_->rowspi(), D->symmetry()));

    int symm = D->symmetry();
    int nirrep = D->nirrep();

    double* temp = new double[Cb_so_->max_ncol() * Cb_so_->max_nrow()];
    for (int h = 0; h < nirrep; h++) {
        int nmol = Cb_so_->colspi()[h];
        int nmor = Cb_so_->colspi()[h^symm];
        int nsol = Cb_so_->rowspi()[h];
        int nsor = Cb_so_->rowspi()[h^symm];
        if (!nmol || !nmor || !nsol || !nsor) continue;
        double** Clp = Cb_so_->pointer(h);
        double** Crp = Cb_so_->pointer(h^symm);
        double** Dmop = D->pointer(h^symm);
        double** Dsop = Db_so_->pointer(h^symm);
        C_DGEMM('N','T',nmol,nsor,nmor,1.0,Dmop[0],nmor,Crp[0],nmor,0.0,temp,nsor);
        C_DGEMM('N','N',nsol,nsor,nmol,1.0,Clp[0],nmol,temp,nsor,0.0,Dsop[0],nsor);
    }
    delete[] temp;
}
void Prop::add(const std::string& prop)
{
    tasks_.insert(prop);
}
void Prop::add(std::vector<std::string> props)
{
    for (int i = 0; i < (int)props.size(); i++) {
        tasks_.insert(props[i]);
    }
}
void Prop::clear()
{
    tasks_.clear();
}
SharedVector Prop::epsilon_a()
{
    return SharedVector(epsilon_a_->clone());
}
SharedVector Prop::epsilon_b()
{
    return SharedVector(epsilon_b_->clone());
}
SharedMatrix Prop::Da_ao()
{
    double* temp = new double[AO2USO_->max_ncol() * AO2USO_->max_nrow()];
    SharedMatrix D = SharedMatrix(new Matrix("Da (AO basis)", basisset_->nbf(), basisset_->nbf()));
    int symm = Da_so_->symmetry();
    for (int h = 0; h < AO2USO_->nirrep(); ++h) {
        int nao = AO2USO_->rowspi()[0];
        int nsol = AO2USO_->colspi()[h];
        int nsor = AO2USO_->colspi()[h^symm];
        if (!nsol || !nsor) continue;
        double** Ulp = AO2USO_->pointer(h);
        double** Urp = AO2USO_->pointer(h^symm);
        double** DSOp = Da_so_->pointer(h^symm);
        double** DAOp = D->pointer();
        C_DGEMM('N','T',nsol,nao,nsor,1.0,DSOp[0],nsor,Urp[0],nsor,0.0,temp,nao);
        C_DGEMM('N','N',nao,nao,nsol,1.0,Ulp[0],nsol,temp,nao,1.0,DAOp[0],nao);
    }
    delete[] temp;
    return D;
}
SharedMatrix Prop::Db_ao()
{
    if (same_dens_)
        throw PSIEXCEPTION("Wavefunction is restricted, asking for Db makes no sense");

    double* temp = new double[AO2USO_->max_ncol() * AO2USO_->max_nrow()];
    SharedMatrix D = SharedMatrix(new Matrix("Db (AO basis)", basisset_->nbf(), basisset_->nbf()));
    int symm = Db_so_->symmetry();
    for (int h = 0; h < AO2USO_->nirrep(); ++h) {
        int nao = AO2USO_->rowspi()[0];
        int nsol = AO2USO_->colspi()[h];
        int nsor = AO2USO_->colspi()[h^symm];
        if (!nsol || !nsor) continue;
        double** Ulp = AO2USO_->pointer(h);
        double** Urp = AO2USO_->pointer(h^symm);
        double** DSOp = Db_so_->pointer(h^symm);
        double** DAOp = D->pointer();
        C_DGEMM('N','T',nsol,nao,nsor,1.0,DSOp[0],nsor,Urp[0],nsor,0.0,temp,nao);
        C_DGEMM('N','N',nao,nao,nsol,1.0,Ulp[0],nsol,temp,nao,1.0,DAOp[0],nao);
    }
    delete[] temp;
    return D;
}
SharedMatrix Prop::Ca_so()
{
    return SharedMatrix(Ca_so_->clone());
}
SharedMatrix Prop::Cb_so()
{
    return SharedMatrix(Cb_so_->clone());
}
SharedMatrix Prop::Ca_ao()
{
    return wfn_->Ca_subset("AO");
}
SharedMatrix Prop::Cb_ao()
{
    return wfn_->Cb_subset("AO");
}
SharedMatrix Prop::Da_so()
{
    return SharedMatrix(Da_so_->clone());
}
SharedMatrix Prop::Db_so()
{
    return SharedMatrix(Db_so_->clone());
}
SharedMatrix Prop::Da_mo()
{
    SharedMatrix D(new Matrix("Da_mo", Ca_so_->colspi(), Ca_so_->colspi(), Da_so_->symmetry()));

    int symm = D->symmetry();
    int nirrep = D->nirrep();

    SharedMatrix S = overlap_so();

    double* SC = new double[Ca_so_->max_ncol() * Ca_so_->max_nrow()];
    double* temp = new double[Ca_so_->max_ncol() * Ca_so_->max_nrow()];
    for (int h = 0; h < nirrep; h++) {
        int nmol = Ca_so_->colspi()[h];
        int nmor = Ca_so_->colspi()[h^symm];
        int nsol = Ca_so_->rowspi()[h];
        int nsor = Ca_so_->rowspi()[h^symm];
        if (!nmol || !nmor || !nsol || !nsor) continue;
        double** Slp = S->pointer(h);
        double** Srp = S->pointer(h^symm);
        double** Clp = Ca_so_->pointer(h);
        double** Crp = Ca_so_->pointer(h^symm);
        double** Dmop = D->pointer(h);
        double** Dsop = Da_so_->pointer(h);

        C_DGEMM('N','N',nsor,nmor,nsor,1.0,Srp[0],nsor,Crp[0],nmor,0.0,SC,nmor);
        C_DGEMM('N','N',nsol,nmor,nsor,1.0,Dsop[0],nsor,SC,nmor,0.0,temp,nmor);
        C_DGEMM('N','N',nsol,nmol,nsol,1.0,Slp[0],nsol,Clp[0],nmol,0.0,SC,nmol);
        C_DGEMM('T','N',nmol,nmor,nsol,1.0,SC,nmol,temp,nmor,0.0,Dmop[0],nmor);
    }
    delete[] temp;
    delete[] SC;
    return D;
}
SharedMatrix Prop::Db_mo()
{
    if (same_dens_)
        throw PSIEXCEPTION("Wavefunction is restricted, asking for Db makes no sense");

    SharedMatrix D(new Matrix("Db_mo", Cb_so_->colspi(), Cb_so_->colspi(), Db_so_->symmetry()));

    int symm = D->symmetry();
    int nirrep = D->nirrep();

    SharedMatrix S = overlap_so();

    double* SC = new double[Cb_so_->max_ncol() * Cb_so_->max_nrow()];
    double* temp = new double[Cb_so_->max_ncol() * Cb_so_->max_nrow()];
    for (int h = 0; h < nirrep; h++) {
        int nmol = Cb_so_->colspi()[h];
        int nmor = Cb_so_->colspi()[h^symm];
        int nsol = Cb_so_->rowspi()[h];
        int nsor = Cb_so_->rowspi()[h^symm];
        if (!nmol || !nmor || !nsol || !nsor) continue;
        double** Slp = S->pointer(h);
        double** Srp = S->pointer(h^symm);
        double** Clp = Cb_so_->pointer(h);
        double** Crp = Cb_so_->pointer(h^symm);
        double** Dmop = D->pointer(h);
        double** Dsop = Db_so_->pointer(h);

        C_DGEMM('N','N',nsor,nmor,nsor,1.0,Srp[0],nsor,Crp[0],nmor,0.0,SC,nmor);
        C_DGEMM('N','N',nsol,nmor,nsor,1.0,Dsop[0],nsor,SC,nmor,0.0,temp,nmor);
        C_DGEMM('N','N',nsol,nmol,nsol,1.0,Slp[0],nsol,Clp[0],nmol,0.0,SC,nmol);
        C_DGEMM('T','N',nmol,nmor,nsol,1.0,SC,nmol,temp,nmor,0.0,Dmop[0],nmor);
    }
    delete[] temp;
    delete[] SC;
    return D;
}
SharedMatrix Prop::Dt_so(bool total)
{
    SharedMatrix Da = Da_so();
    SharedMatrix D(Da->clone());
    if (same_dens_) {
        D->set_name((total ? "Dt_so" : "Ds_so"));
        D->scale((total ? 2.0 : 0.0));
    } else {
        D->set_name((total ? "Dt_so" : "Ds_so"));
        SharedMatrix Db = Db_so();
        if (total) D->add(Db);
        else D->subtract(Db);
    }
    return D;
}
SharedMatrix Prop::Dt_mo(bool total)
{
    SharedMatrix Da = Da_mo();
    if (same_dens_) {
        Da->set_name((total ? "Dt_mo" : "Ds_mo"));
        Da->scale((total ? 2.0 : 0.0));
    } else {
        Da->set_name((total ? "Dt_mo" : "Ds_mo"));
        SharedMatrix Db = Db_mo();
        if (total) Da->add(Db);
        else Da->subtract(Db);
    }
    return Da;
}
std::pair<SharedMatrix, SharedVector> Prop::Na_mo()
{
    SharedMatrix D = Da_mo();
    SharedMatrix C(new Matrix("Na_mo", D->nirrep(), D->rowspi(), D->rowspi()));
    std::shared_ptr<Vector> O(new Vector("Alpha Occupation", D->rowspi()));

    D->diagonalize(C,O,descending);

    return make_pair(C,O);
}
std::pair<SharedMatrix, SharedVector> Prop::Nb_mo()
{
    if (same_dens_)
        throw PSIEXCEPTION("Wavefunction is restricted, asking for Nb makes no sense");

    SharedMatrix D = Db_mo();
    SharedMatrix C(new Matrix("Nb_mo", D->nirrep(), D->rowspi(), D->rowspi()));
    std::shared_ptr<Vector> O(new Vector("Beta Occupation", D->rowspi()));

    D->diagonalize(C,O,descending);

    return make_pair(C,O);
}
std::pair<SharedMatrix, SharedVector> Prop::Na_so()
{
    std::pair<SharedMatrix, std::shared_ptr<Vector> > pair = Na_mo();
    SharedMatrix N = pair.first;
    std::shared_ptr<Vector> O = pair.second;

    SharedMatrix N2(new Matrix("Na_so", Ca_so_->nirrep(), Ca_so_->rowspi(), Ca_so_->colspi()));

    for (int h = 0; h < N->nirrep(); h++) {

        int nmo = Ca_so_->colspi()[h];
        int nso = Ca_so_->rowspi()[h];

        if (!nmo || !nso) continue;

        double** Np = N->pointer(h);
        double** Cp = Ca_so_->pointer(h);
        double** N2p = N2->pointer(h);

        C_DGEMM('N','N',nso,nmo,nmo,1.0,Cp[0],nmo,Np[0],nmo,0.0,N2p[0],nmo);
    }
    return make_pair(N2,O);
}
std::pair<SharedMatrix, SharedVector> Prop::Nb_so()
{
    if (same_dens_)
        throw PSIEXCEPTION("Wavefunction is restricted, asking for Nb makes no sense");

    std::pair<SharedMatrix, std::shared_ptr<Vector> > pair = Nb_mo();
    SharedMatrix N = pair.first;
    std::shared_ptr<Vector> O = pair.second;

    SharedMatrix N2(new Matrix("Nb_so", Cb_so_->nirrep(), Cb_so_->rowspi(), Cb_so_->colspi()));

    for (int h = 0; h < N->nirrep(); h++) {

        int nmo = Cb_so_->colspi()[h];
        int nso = Cb_so_->rowspi()[h];

        if (!nmo || !nso) continue;

        double** Np = N->pointer(h);
        double** Cp = Cb_so_->pointer(h);
        double** N2p = N2->pointer(h);

        C_DGEMM('N','N',nso,nmo,nmo,1.0,Cp[0],nmo,Np[0],nmo,0.0,N2p[0],nmo);
    }
    return make_pair(N2,O);
}
std::pair<SharedMatrix, SharedVector> Prop::Na_ao()
{
    std::pair<SharedMatrix, std::shared_ptr<Vector> > pair = Na_so();
    SharedMatrix N = pair.first;
    std::shared_ptr<Vector> O = pair.second;

    SharedMatrix N2(new Matrix("Na_ao", Ca_so_->nrow(), Ca_so_->ncol()));
    SharedMatrix N3(new Matrix("Na_ao", Ca_so_->nrow(), Ca_so_->ncol()));
    std::shared_ptr<Vector> O2(new Vector("Alpha Occupation", Ca_so_->ncol()));

    int offset = 0;
    std::vector<std::pair<double,int> > index;
    for (int h = 0; h < Ca_so_->nirrep(); h++) {

        int ncol = Ca_so_->ncol();
        int nmo = Ca_so_->colspi()[h];
        int nso = AO2USO_->colspi()[h];
        int nao = AO2USO_->rowspi()[h];

        if (!nmo || !nso || !nao) continue;

        for (int i = 0; i < nmo; i++) {
            index.push_back(std::make_pair(O->get(h,i),i+offset));
        }

        double** Np = N->pointer(h);
        double** Up = AO2USO_->pointer(h);
        double** N2p = N2->pointer(h);

        C_DGEMM('N','N',nao,nmo,nso,1.0,Up[0],nso,Np[0],nmo,0.0,&N2p[0][offset],ncol);

        offset += nmo;
    }

    std::sort(index.begin(), index.end(), std::greater<std::pair<double,int> >());

    int nmo = N2->colspi()[0];
    int nao = N2->rowspi()[0];

    for (int i = 0; i < nmo; i++) {
        double occ = index[i].first;
        int ind    = index[i].second;
        O2->set(i,occ);

        C_DCOPY(nao, &(N2->pointer()[0][ind]), nmo, &(N3->pointer()[0][i]), nmo);
    }

    return make_pair(N3,O2);
}
std::pair<SharedMatrix, SharedVector> Prop::Nb_ao()
{
    if (same_dens_)
        throw PSIEXCEPTION("Wavefunction is restricted, asking for Nb makes no sense");

    std::pair<SharedMatrix, std::shared_ptr<Vector> > pair = Nb_so();
    SharedMatrix N = pair.first;
    std::shared_ptr<Vector> O = pair.second;

    SharedMatrix N2(new Matrix("Nb_ao", Cb_so_->nrow(), Cb_so_->ncol()));
    SharedMatrix N3(new Matrix("Nb_ao", Cb_so_->nrow(), Cb_so_->ncol()));
    std::shared_ptr<Vector> O2(new Vector("Beta Occupation", Cb_so_->ncol()));

    int offset = 0;
    std::vector<std::pair<double,int> > index;
    for (int h = 0; h < Cb_so_->nirrep(); h++) {

        int ncol = Cb_so_->ncol();
        int nmo = Cb_so_->colspi()[h];
        int nso = AO2USO_->colspi()[h];
        int nao = AO2USO_->rowspi()[h];

        if (!nmo || !nso || !nao) continue;

        for (int i = 0; i < nmo; i++) {
            index.push_back(std::make_pair(O->get(h,i),i+offset));
        }

        double** Np = N->pointer(h);
        double** Up = AO2USO_->pointer(h);
        double** N2p = N2->pointer(h);

        C_DGEMM('N','N',nao,nmo,nso,1.0,Up[0],nso,Np[0],nmo,0.0,&N2p[0][offset],ncol);

        offset += nmo;
    }

    std::sort(index.begin(), index.end(), std::greater<std::pair<double,int> >());

    int nmo = N2->colspi()[0];
    int nao = N2->rowspi()[0];

    for (int i = 0; i < nmo; i++) {
        double occ = index[i].first;
        int ind    = index[i].second;
        O2->set(i,occ);

        C_DCOPY(nao, &(N2->pointer()[0][ind]), nmo, &(N3->pointer()[0][i]), nmo);
    }

    return make_pair(N3,O2);
}
std::pair<SharedMatrix, SharedVector> Prop::Nt_mo()
{
    SharedMatrix D = Dt_mo();
    SharedMatrix C(new Matrix("Nt_mo", D->nirrep(), D->rowspi(), D->rowspi()));
    std::shared_ptr<Vector> O(new Vector("Total Occupation", D->rowspi()));

    D->diagonalize(C,O,descending);

    return make_pair(C,O);
}
std::pair<SharedMatrix, SharedVector> Prop::Nt_so()
{
    std::pair<SharedMatrix, std::shared_ptr<Vector> > pair = Nt_mo();
    SharedMatrix N = pair.first;
    std::shared_ptr<Vector> O = pair.second;

    SharedMatrix N2(new Matrix("Nt_so", Cb_so_->nirrep(), Cb_so_->rowspi(), Cb_so_->colspi()));

    for (int h = 0; h < N->nirrep(); h++) {

        int nmo = Cb_so_->colspi()[h];
        int nso = Cb_so_->rowspi()[h];

        if (!nmo || !nso) continue;

        double** Np = N->pointer(h);
        double** Cp = Cb_so_->pointer(h);
        double** N2p = N2->pointer(h);

        C_DGEMM('N','N',nso,nmo,nmo,1.0,Cp[0],nmo,Np[0],nmo,0.0,N2p[0],nmo);
    }
    return make_pair(N2,O);
}
std::pair<SharedMatrix, SharedVector> Prop::Nt_ao()
{
    std::pair<SharedMatrix, std::shared_ptr<Vector> > pair = Nt_so();
    SharedMatrix N = pair.first;
    std::shared_ptr<Vector> O = pair.second;

    SharedMatrix N2(new Matrix("Nt_ao", Cb_so_->nrow(), Cb_so_->ncol()));
    SharedMatrix N3(new Matrix("Nt_ao", Cb_so_->nrow(), Cb_so_->ncol()));
    std::shared_ptr<Vector> O2(new Vector("Total Occupation", Cb_so_->ncol()));

    int offset = 0;
    std::vector<std::pair<double,int> > index;
    for (int h = 0; h < Cb_so_->nirrep(); h++) {

        int ncol = Cb_so_->ncol();
        int nmo = Cb_so_->colspi()[h];
        int nso = AO2USO_->colspi()[h];
        int nao = AO2USO_->rowspi()[h];

        if (!nmo || !nso || !nao) continue;

        for (int i = 0; i < nmo; i++) {
            index.push_back(std::make_pair(O->get(h,i),i+offset));
        }

        double** Np = N->pointer(h);
        double** Up = AO2USO_->pointer(h);
        double** N2p = N2->pointer(h);

        C_DGEMM('N','N',nao,nmo,nso,1.0,Up[0],nso,Np[0],nmo,0.0,&N2p[0][offset],ncol);

        offset += nmo;
    }

    std::sort(index.begin(), index.end(), std::greater<std::pair<double,int> >());

    int nmo = N2->colspi()[0];
    int nao = N2->rowspi()[0];

    for (int i = 0; i < nmo; i++) {
        double occ = index[i].first;
        int ind    = index[i].second;
        O2->set(i,occ);

        C_DCOPY(nao, &(N2->pointer()[0][ind]), nmo, &(N3->pointer()[0][i]), nmo);
    }

    return make_pair(N3,O2);
}
SharedMatrix Prop::overlap_so()
{
    SharedMatrix S(new Matrix("S",Ca_so_->rowspi(),Ca_so_->rowspi()));
    std::shared_ptr<OneBodySOInt> Sint(integral_->so_overlap());
    Sint->compute(S);
    return S;
}

OEProp::OEProp(std::shared_ptr<Wavefunction> wfn) : Prop(wfn)
{
    common_init();
}
OEProp::~OEProp()
{
}

Vector3 OEProp::compute_center(const double *property) const
{
    std::shared_ptr<Molecule> mol = basisset_->molecule();
    int natoms = mol->natom();
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double sum = 0.0;
    for(int atom = 0; atom < natoms; ++atom){
        Vector3 xyz = mol->xyz(atom);
        double prop = property[atom];
        x += xyz[0] * prop;
        y += xyz[1] * prop;
        z += xyz[2] * prop;
    }
    x /= sum;
    y /= sum;
    z /= sum;
    return Vector3(x, y, z);
}

void OEProp::common_init()
{
    // See if the user specified the origin
    Options &options = Process::environment.options;

    print_ = options.get_int("PRINT");

    std::shared_ptr<Molecule> mol = basisset_->molecule();
    int natoms = mol->natom();
    if(options["PROPERTIES_ORIGIN"].has_changed()){
        int size = options["PROPERTIES_ORIGIN"].size();

        if(size == 1){
            double *property = new double[natoms];
            std::string str = options["PROPERTIES_ORIGIN"].to_string();
            if(str == "COM"){
                for(int atom = 0; atom < natoms; ++atom)
                    property[atom] = mol->mass(atom);
            }else if(str == "NUCLEAR_CHARGE"){
                for(int atom = 0; atom < natoms; ++atom)
                    property[atom] = mol->charge(atom);
            }else{
                throw PSIEXCEPTION("Invalid specification of PROPERTIES_ORIGIN.  Please consult the manual.");
            }
            origin_ = compute_center(property);
            delete [] property;
        }else if(size == 3){
            double x = options["PROPERTIES_ORIGIN"][0].to_double();
            double y = options["PROPERTIES_ORIGIN"][1].to_double();
            double z = options["PROPERTIES_ORIGIN"][2].to_double();
            bool convert = mol->units() == Molecule::Angstrom;
            if(convert){
                x /= pc_bohr2angstroms;
                y /= pc_bohr2angstroms;
                z /= pc_bohr2angstroms;
            }
            origin_ = Vector3(x, y, z);
        }else{
            throw PSIEXCEPTION("Invalid specification of PROPERTIES_ORIGIN.  Please consult the manual.");
        }
    }
    outfile->Printf( "\n\nProperties will be evaluated at %10.6f, %10.6f, %10.6f Bohr\n",
            origin_[0], origin_[1], origin_[2]);

    // Determine number of NOONs to print; default is 3
    if(options.get_str("PRINT_NOONS") == "ALL") max_noon_ = wfn_->nmo();
    else max_noon_ = to_integer(options.get_str("PRINT_NOONS"));

    /*
     * Now check the symmetry of the origin; if it's off-axis we can't use symmetry for multipoles anymore
     */
    CharacterTable ct = mol->point_group()->char_table();
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
            for (int xyz = 0; xyz < 3; ++xyz)
                t[xyz] += origin_[xyz]*rr(xyz, xyz) * gamma.character(G) / nirrep;
        }

        for (int xyz = 0; xyz < 3; ++xyz) {
            if(fabs(t[xyz]) > 1.0E-8){
                outfile->Printf( "The origin chosen breaks symmetry; multipoles will be computed without symmetry.\n");
                origin_preserves_symmetry_ = false;
            }
        }
    }

}

void OEProp::print_header()
{
    outfile->Printf( "\n OEPROP: One-electron properties/analyses.\n");
    outfile->Printf( "  by Rob Parrish and Justin Turney.\n");
    outfile->Printf( "  built on LIBMINTS.\n\n");
}

template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
}

void OEProp::compute()
{

    outfile->Printf( "\nProperties computed using the %s density matrix\n\n", title_.c_str());

    // Search for multipole strings, which are handled separately
    std::set<std::string>::const_iterator iter = tasks_.begin();
    std::regex mpoles("^MULTIPOLE(?:S)?\\s*\\((\\d+)\\)$");
    std::smatch matches;
    for(; iter !=tasks_.end(); ++iter){
        std::string str = *iter;
        if(std::regex_match(str, matches, mpoles)){
            int order;
            if(!from_string<int>(order, matches[1], std::dec))
                throw PSIEXCEPTION("Problem detemining multipole order!  Specify, e.g., MULTIPOLE(5)");
            compute_multipoles(order);
        }
    }
    // print_header();  // Not by default, happens too often -CDS
    if (tasks_.count("ESP_AT_NUCLEI"))
        compute_esp_at_nuclei();
    if (tasks_.count("DIPOLE"))
        compute_dipole(false);
    if (tasks_.count("QUADRUPOLE"))
        compute_quadrupole(false);
    if (tasks_.count("TRANSITION_DIPOLE"))
        compute_dipole(true);
    if (tasks_.count("TRANSITION_QUADRUPOLE"))
        compute_quadrupole(true);
    if (tasks_.count("MO_EXTENTS"))
        compute_mo_extents();
    if (tasks_.count("MULLIKEN_CHARGES"))
        compute_mulliken_charges();
    if (tasks_.count("LOWDIN_CHARGES"))
        compute_lowdin_charges();
    if (tasks_.count("MAYER_INDICES"))
        compute_mayer_indices();
    if (tasks_.count("WIBERG_LOWDIN_INDICES"))
        compute_wiberg_lowdin_indices();
    if (tasks_.count("NO_OCCUPATIONS"))
        compute_no_occupations();
    if (tasks_.count("GRID_FIELD"))
        compute_field_over_grid();
    if (tasks_.count("GRID_ESP"))
        compute_esp_over_grid();
}


void OEProp::compute_multipoles(int order, bool transition)
{
    std::shared_ptr<Molecule> mol = basisset_->molecule();

    SharedMatrix Da;
    SharedMatrix Db;

    std::vector<SharedMatrix> mp_ints;
    MultipoleSymmetry mpsymm (order, mol, integral_, factory_);

    if(origin_preserves_symmetry_){
        // We can use symmetry here
        std::shared_ptr<OneBodySOInt> sompOBI (integral_->so_multipoles(order));
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
    }else{
        // No symmetry
        mp_ints = mpsymm.create_matrices("", true);
        std::shared_ptr<OneBodyAOInt> aompOBI (integral_->ao_multipoles(order));
        aompOBI->set_origin(origin_);
        aompOBI->compute(mp_ints);

        if (same_dens_) {
            Da = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
            Db = Da;
        } else {
            Da = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
            Db = wfn_->D_subset_helper(Db_so_, Cb_so_, "AO");
        }
    }

    if(print_ > 4){
        std::vector<SharedMatrix>::iterator iter;
        for(iter = mp_ints.begin(); iter != mp_ints.end(); ++iter)
            iter->get()->print();
    }

    SharedVector nuclear_contributions = MultipoleInt::nuclear_contribution(mol, order, origin_);

    outfile->Printf("\n%s Multipole Moments:\n", transition ? "Transition" : "");
    outfile->Printf( "\n ------------------------------------------------------------------------------------\n");
    outfile->Printf( "     Multipole             Electric (a.u.)       Nuclear  (a.u.)        Total (a.u.)\n");
    outfile->Printf( " ------------------------------------------------------------------------------------\n\n");
    double convfac = pc_dipmom_au2debye;
    int address = 0;
    for(int l = 1; l <= order; ++l){
        int ncomponents = (l + 1) * (l + 2) / 2;
        std::stringstream ss;
        if(l > 1)
            ss << ".ang";
        if(l > 2)
            ss << "^" << l-1;
        std::string exp = ss.str();
        outfile->Printf( " L = %d.  Multiply by %.10f to convert to Debye%s\n", l, convfac, exp.c_str());
        for(int component = 0; component < ncomponents; ++component){
            SharedMatrix mpmat = mp_ints[address];
            std::string name = mpmat->name();
            double nuc = transition ? 0.0 : nuclear_contributions->get(address);
            double elec = Da->vector_dot(mpmat) + Db->vector_dot(mpmat);
            double tot = nuc + elec;
            outfile->Printf( " %-20s: %18.7f   %18.7f   %18.7f\n",
                    name.c_str(), elec, nuc, tot);
            std::string upper_name = to_upper_copy(name);
            /*- Process::environment.globals["DIPOLE X"] -*/
            /*- Process::environment.globals["DIPOLE Y"] -*/
            /*- Process::environment.globals["32-POLE XXXXX"] -*/
            /*- Process::environment.globals["32-POLE XXXXY"] -*/
            Process::environment.globals[upper_name] = tot;
            ++address;
        }
        outfile->Printf( "\n");
        convfac *= pc_bohr2angstroms;
    }
    outfile->Printf( " --------------------------------------------------------------------------------\n");


}


/**
 * @brief The GridIterator class:  A class to iterate over a user-provided grid for
 *                                 computing electrostatic properties.
 */
class GridIterator{
    std::ifstream gridfile_;
    Vector3 gridpoints_;
public:
    const Vector3& gridpoints() const {return gridpoints_;}
    GridIterator(const std::string &filename)
    {
        gridfile_.open(filename.c_str());
        if(!gridfile_)
            throw PSIEXCEPTION("Unable to open the grid.dat file.");
    }
    void first()
    {
        if(!gridfile_)
            throw PSIEXCEPTION("File not initialized in Griditer::first.");
        gridfile_.clear();
        gridfile_.seekg(0, std::ios::beg);
        next();
    }
    void next()
    {
        if(!gridfile_)
            throw PSIEXCEPTION("Griditer::next called before file stream was initialized.");
        if(!(gridfile_ >> gridpoints_[0])){
            if(gridfile_.eof()){
                // Hitting the end of file is OK in this situation.
                return;
            }else{
                throw PSIEXCEPTION("Problem reading x gridpoint from the grid file.");
            }
        }
        if(!(gridfile_ >> gridpoints_[1]))
            throw PSIEXCEPTION("Problem reading y gridpoint from the grid file.");
        if(!(gridfile_ >> gridpoints_[2]))
            throw PSIEXCEPTION("Problem reading z gridpoint from the grid file.");
    }
    bool last() const
    {
        return gridfile_.eof();
    }
    ~GridIterator()
    {
        gridfile_.close();
    }
};

void OEProp::compute_esp_over_grid()
{
    std::shared_ptr<Molecule> mol = basisset_->molecule();

    std::shared_ptr<ElectrostaticInt> epot(dynamic_cast<ElectrostaticInt*>(integral_->electrostatic()));

    outfile->Printf( "\n Electrostatic potential computed on the grid and written to grid_esp.dat\n");

    SharedMatrix Dtot = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
    if (same_dens_) {
        Dtot->scale(2.0);
    }else{
        Dtot->add(wfn_->D_subset_helper(Db_so_, Cb_so_, "AO"));
    }

    int nbf = basisset_->nbf();
    SharedMatrix ints(new Matrix("Ex integrals", nbf, nbf));

    Vvals_.clear();
    FILE *gridout = fopen("grid_esp.dat", "w");
    if(!gridout)
        throw PSIEXCEPTION("Unable to write to grid_esp.dat");
    GridIterator griditer("grid.dat");
    for(griditer.first(); !griditer.last(); griditer.next()){
        Vector3 origin(griditer.gridpoints());
        if(mol->units() == Molecule::Angstrom)
            origin /= pc_bohr2angstroms;
        ints->zero();
        epot->compute(ints, origin);
        double Velec = Dtot->vector_dot(ints);
        double Vnuc = 0.0;
        int natom = mol->natom();
        for(int i=0; i < natom; i++) {
            Vector3 dR = origin - mol->xyz(i);
            double r = dR.norm();
            if(r > 1.0E-8)
                Vnuc += mol->Z(i)/r;
        }
        Vvals_.push_back(Velec+Vnuc);
        fprintf(gridout, "%16.10f\n", Velec+Vnuc);
    }
    fclose(gridout);
}


void OEProp::compute_field_over_grid()
{
    std::shared_ptr<Molecule> mol = basisset_->molecule();

    std::shared_ptr<ElectrostaticInt> epot(dynamic_cast<ElectrostaticInt*>(integral_->electrostatic()));

    outfile->Printf( "\n Field computed on the grid and written to grid_field.dat\n");

    SharedMatrix Dtot = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
    if (same_dens_) {
        Dtot->scale(2.0);
    }else{
        Dtot->add(wfn_->D_subset_helper(Db_so_, Cb_so_, "AO"));
    }

    std::shared_ptr<ElectricFieldInt> field_ints(dynamic_cast<ElectricFieldInt*>(wfn_->integral()->electric_field()));

    int nbf = basisset_->nbf();
    std::vector<SharedMatrix> intmats;
    intmats.push_back(SharedMatrix(new Matrix("Ex integrals", nbf, nbf)));
    intmats.push_back(SharedMatrix(new Matrix("Ey integrals", nbf, nbf)));
    intmats.push_back(SharedMatrix(new Matrix("Ez integrals", nbf, nbf)));

    Exvals_.clear();
    Eyvals_.clear();
    Ezvals_.clear();

    FILE *gridout = fopen("grid_field.dat", "w");
    if(!gridout)
        throw PSIEXCEPTION("Unable to write to grid_field.dat");
    GridIterator griditer("grid.dat");
    for(griditer.first(); !griditer.last(); griditer.next()){
        Vector3 origin(griditer.gridpoints());
        if(mol->units() == Molecule::Angstrom)
            origin /= pc_bohr2angstroms;
        field_ints->set_origin(origin);
        for (int m=0; m<3; ++m)
            intmats[m]->zero();
        field_ints->compute(intmats);
        double Ex = Dtot->vector_dot(intmats[0]);
        double Ey = Dtot->vector_dot(intmats[1]);
        double Ez = Dtot->vector_dot(intmats[2]);
        Vector3 nuc = field_ints->nuclear_contribution(origin, mol);
        Exvals_.push_back(Ex+nuc[0]);
        Eyvals_.push_back(Ey+nuc[1]);
        Ezvals_.push_back(Ez+nuc[2]);
        fprintf(gridout, "%16.10f %16.10f %16.10f\n", Ex+nuc[0], Ey+nuc[1], Ez+nuc[2]);
    }
    fclose(gridout);
}


void OEProp::compute_esp_at_nuclei()
{
    std::shared_ptr<Molecule> mol = basisset_->molecule();

    std::shared_ptr<ElectrostaticInt> epot(dynamic_cast<ElectrostaticInt*>(integral_->electrostatic()));

    int nbf = basisset_->nbf();
    int natoms = mol->natom();

    SharedMatrix Dtot = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
    if (same_dens_) {
        Dtot->scale(2.0);
    }else{
        Dtot->add(wfn_->D_subset_helper(Db_so_, Cb_so_, "AO"));
    }

    Matrix dist = mol->distance_matrix();
    outfile->Printf( "\n Electrostatic potentials at the nuclear coordinates:\n");
    outfile->Printf( " ---------------------------------------------\n");
    outfile->Printf( "   Center     Electrostatic Potential (a.u.)\n");
    outfile->Printf( " ---------------------------------------------\n");
    for(int atom1 = 0; atom1 < natoms; ++atom1){
        std::stringstream s;
        s << "ESP AT CENTER " << atom1+1;
        SharedMatrix ints(new Matrix(s.str(), nbf, nbf));
        epot->compute(ints, mol->xyz(atom1));
        if(print_ > 2)
            ints->print();
        double elec = Dtot->vector_dot(ints);
        double nuc = 0.0;
        for(int atom2 = 0; atom2 < natoms; ++atom2){
            if(atom1 == atom2)
                continue;
            nuc += mol->Z(atom2) / dist[0][atom1][atom2];
        }
        outfile->Printf( "  %3d %2s           %16.12f\n",
                atom1+1, mol->label(atom1).c_str(), nuc+elec);
        /*- Process::environment.globals["ESP AT CENTER n"] -*/
        Process::environment.globals[s.str()] = nuc+elec;
    }
    outfile->Printf( " ---------------------------------------------\n");
}

void OEProp::compute_dipole(bool transition)
{
    std::shared_ptr<Molecule> mol = basisset_->molecule();

    Vector3 de;
    SharedMatrix Da;
    SharedMatrix Db;
    std::vector<SharedMatrix> dipole_ints;
    if (origin_preserves_symmetry_) {
        // Use symmetry
        OperatorSymmetry dipsymm(1, mol, integral_, factory_);
        dipole_ints = dipsymm.create_matrices("SO Dipole");
        std::shared_ptr<OneBodySOInt> sodOBI(integral_->so_dipole());
        sodOBI->ob()->set_origin(origin_);
        sodOBI->compute(dipole_ints);
        if (same_dens_) {
            Da = Da_so_;
            Db = Da;
        } else {
            Da = Da_so_;
            Db = Db_so_;
        }
    } else {
        // Don't use symmetry
        dipole_ints.push_back(SharedMatrix(
            new Matrix("AO Dipole X", basisset_->nbf(), basisset_->nbf())));
        dipole_ints.push_back(SharedMatrix(
            new Matrix("AO Dipole Y", basisset_->nbf(), basisset_->nbf())));
        dipole_ints.push_back(SharedMatrix(
            new Matrix("AO Dipole Z", basisset_->nbf(), basisset_->nbf())));
        std::shared_ptr<OneBodyAOInt> aodOBI(integral_->ao_dipole());
        aodOBI->set_origin(origin_);
        aodOBI->compute(dipole_ints);
        if (same_dens_) {
            Da = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
            Db = Da;
        } else {
            Da = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
            Db = wfn_->D_subset_helper(Db_so_, Cb_so_, "AO");
        }
    }

    if (print_ > 4){
        for (int n = 0; n < 3; ++n) dipole_ints[n]->print();
    }

    de[0] = Da->vector_dot(dipole_ints[0]) + Db->vector_dot(dipole_ints[0]);
    de[1] = Da->vector_dot(dipole_ints[1]) + Db->vector_dot(dipole_ints[1]);
    de[2] = Da->vector_dot(dipole_ints[2]) + Db->vector_dot(dipole_ints[2]);

    SharedVector ndip = DipoleInt::nuclear_contribution(mol, origin_);

    if (!transition) {

        outfile->Printf( "  Nuclear Dipole Moment: (a.u.)\n");
        outfile->Printf("     X: %10.4lf      Y: %10.4lf      Z: %10.4lf\n",
                ndip->get(0), ndip->get(1), ndip->get(2));
        outfile->Printf( "\n");
        outfile->Printf( "  Electronic Dipole Moment: (a.u.)\n");
        outfile->Printf("     X: %10.4lf      Y: %10.4lf      Z: %10.4lf\n",
                de[0], de[1], de[2]);
        outfile->Printf( "\n");

        de[0] += ndip->get(0, 0);
        de[1] += ndip->get(0, 1);
        de[2] += ndip->get(0, 2);
    }

    outfile->Printf("  %sDipole Moment: (a.u.)\n", (transition ? "Transition " : ""));
    outfile->Printf("     X: %10.4lf      Y: %10.4lf      Z: %10.4lf     Total: %10.4lf\n",
       de[0], de[1], de[2], de.norm());
    outfile->Printf( "\n");

    double dfac = pc_dipmom_au2debye;
    outfile->Printf("  %sDipole Moment: (Debye)\n", (transition ? "Transition " : ""));
    outfile->Printf("     X: %10.4lf      Y: %10.4lf      Z: %10.4lf     Total: %10.4lf\n",
       de[0]*dfac, de[1]*dfac, de[2]*dfac, de.norm()*dfac);
    outfile->Printf( "\n");

    // Dipole components in Debye
    std::stringstream s;
    s << title_ << " DIPOLE X";
    Process::environment.globals[s.str()] = de[0]*dfac;
    s.str(std::string());
    s << title_ << " DIPOLE Y";
    Process::environment.globals[s.str()] = de[1]*dfac;
    s.str(std::string());
    s << title_ << " DIPOLE Z";
    Process::environment.globals[s.str()] = de[2]*dfac;


}

void OEProp::compute_quadrupole(bool transition)
{
    std::shared_ptr<Molecule> mol = basisset_->molecule();
    SharedMatrix Da;
    SharedMatrix Db;

    std::vector<SharedMatrix> qpole_ints;
    if(origin_preserves_symmetry_){
        OperatorSymmetry quadsymm(2, mol, integral_, factory_);
        qpole_ints = quadsymm.create_matrices("SO Quadrupole");
        std::shared_ptr<OneBodySOInt> soqOBI(integral_->so_quadrupole());
        soqOBI->ob()->set_origin(origin_);
        soqOBI->compute(qpole_ints);
        if (same_dens_) {
            Da = Da_so_;
            Db = Da;
        } else {
            Da = Da_so_;
            Db = Db_so_;
        }
    }else{
        qpole_ints.push_back(SharedMatrix(new Matrix("AO Quadrupole XX", basisset_->nbf(), basisset_->nbf())));
        qpole_ints.push_back(SharedMatrix(new Matrix("AO Quadrupole XY", basisset_->nbf(), basisset_->nbf())));
        qpole_ints.push_back(SharedMatrix(new Matrix("AO Quadrupole XZ", basisset_->nbf(), basisset_->nbf())));
        qpole_ints.push_back(SharedMatrix(new Matrix("AO Quadrupole YY", basisset_->nbf(), basisset_->nbf())));
        qpole_ints.push_back(SharedMatrix(new Matrix("AO Quadrupole YZ", basisset_->nbf(), basisset_->nbf())));
        qpole_ints.push_back(SharedMatrix(new Matrix("AO Quadrupole ZZ", basisset_->nbf(), basisset_->nbf())));
        std::shared_ptr<OneBodyAOInt> aoqOBI(integral_->ao_quadrupole());
        aoqOBI->set_origin(origin_);
        aoqOBI->compute(qpole_ints);
        if (same_dens_) {
            Da = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
            Db = Da;
        } else {
            Da = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
            Db = wfn_->D_subset_helper(Db_so_, Cb_so_, "AO");
        }
    }

    if(print_ > 4)
        for(int n = 0; n < 6; ++n)
            qpole_ints[n]->print();

    // Each multipole integral needs to be dotted with the SO Density
    Vector qe(6);
    qe[0] = Da->vector_dot(qpole_ints[0]) + Db->vector_dot(qpole_ints[0]);
    qe[1] = Da->vector_dot(qpole_ints[1]) + Db->vector_dot(qpole_ints[1]);
    qe[2] = Da->vector_dot(qpole_ints[2]) + Db->vector_dot(qpole_ints[2]);
    qe[3] = Da->vector_dot(qpole_ints[3]) + Db->vector_dot(qpole_ints[3]);
    qe[4] = Da->vector_dot(qpole_ints[4]) + Db->vector_dot(qpole_ints[4]);
    qe[5] = Da->vector_dot(qpole_ints[5]) + Db->vector_dot(qpole_ints[5]);

    // Add in nuclear contribution
    if (!transition) {
        SharedVector nquad = QuadrupoleInt::nuclear_contribution(mol, origin_);
        qe[0] += nquad->get(0, 0);
        qe[1] += nquad->get(0, 1);
        qe[2] += nquad->get(0, 2);
        qe[3] += nquad->get(0, 3);
        qe[4] += nquad->get(0, 4);
        qe[5] += nquad->get(0, 5);
    }

    // Print multipole components
    double dfac = pc_dipmom_au2debye * pc_bohr2angstroms;
    outfile->Printf( "  %sQuadrupole Moment: (Debye Ang)\n", (transition ? "Transition " : ""));
    outfile->Printf( "    XX: %10.4lf     YY: %10.4lf     ZZ: %10.4lf\n", \
       qe[0]*dfac, qe[3]*dfac, qe[5]*dfac);
    outfile->Printf( "    XY: %10.4lf     XZ: %10.4lf     YZ: %10.4lf\n", \
       qe[1]*dfac, qe[2]*dfac, qe[4]*dfac);
    outfile->Printf( "\n");

    double dtrace = (1.0 / 3.0) * (qe[0] + qe[3] + qe[5]);
    outfile->Printf( "  Traceless %sQuadrupole Moment: (Debye Ang)\n", (transition ? "Transition " : ""));
    outfile->Printf( "    XX: %10.4lf     YY: %10.4lf     ZZ: %10.4lf\n", \
       (qe[0]-dtrace)*dfac, (qe[3]-dtrace)*dfac, (qe[5]-dtrace)*dfac);
    outfile->Printf( "    XY: %10.4lf     XZ: %10.4lf     YZ: %10.4lf\n", \
       qe[1]*dfac, qe[2]*dfac, qe[4]*dfac);
    outfile->Printf( "\n");

    // Quadrupole components in Debye Ang
    std::stringstream s;
    s << title_ << " QUADRUPOLE XX";
    Process::environment.globals[s.str()] = qe[0]*dfac;
    s.str(std::string());
    s << title_ << " QUADRUPOLE YY";
    Process::environment.globals[s.str()] = qe[3]*dfac;
    s.str(std::string());
    s << title_ << " QUADRUPOLE ZZ";
    Process::environment.globals[s.str()] = qe[5]*dfac;
    s.str(std::string());
    s << title_ << " QUADRUPOLE XY";
    Process::environment.globals[s.str()] = qe[1]*dfac;
    s.str(std::string());
    s << title_ << " QUADRUPOLE XZ";
    Process::environment.globals[s.str()] = qe[2]*dfac;
    s.str(std::string());
    s << title_ << " QUADRUPOLE YZ";
    Process::environment.globals[s.str()] = qe[4]*dfac;


}
void OEProp::compute_mo_extents()
{
    std::shared_ptr<Molecule> mol = basisset_->molecule();
    SharedMatrix Ca;
    SharedMatrix Cb;

    if (same_orbs_) {
        Ca = Ca_ao();
        Cb = Ca;
    } else {
        Ca = Ca_ao();
        Cb = Cb_ao();
    }

    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> ao_Dpole;
    std::vector<SharedMatrix> ao_Qpole;

    ao_Dpole.push_back(SharedMatrix(new Matrix("Dipole X", basisset_->nbf(), basisset_->nbf())));
    ao_Dpole.push_back(SharedMatrix(new Matrix("Dipole Y", basisset_->nbf(), basisset_->nbf())));
    ao_Dpole.push_back(SharedMatrix(new Matrix("Dipole Z", basisset_->nbf(), basisset_->nbf())));

    ao_Qpole.push_back(SharedMatrix(new Matrix("Quadrupole XX", basisset_->nbf(), basisset_->nbf())));
    ao_Qpole.push_back(SharedMatrix(new Matrix("Quadrupole XY", basisset_->nbf(), basisset_->nbf())));
    ao_Qpole.push_back(SharedMatrix(new Matrix("Quadrupole XZ", basisset_->nbf(), basisset_->nbf())));
    ao_Qpole.push_back(SharedMatrix(new Matrix("Quadrupole YY", basisset_->nbf(), basisset_->nbf())));
    ao_Qpole.push_back(SharedMatrix(new Matrix("Quadrupole YZ", basisset_->nbf(), basisset_->nbf())));
    ao_Qpole.push_back(SharedMatrix(new Matrix("Quadrupole ZZ", basisset_->nbf(), basisset_->nbf())));

    // Form the one-electron integral objects from the integral factory
    std::shared_ptr<OneBodyAOInt> aodOBI(integral_->ao_dipole());
    std::shared_ptr<OneBodyAOInt> aoqOBI(integral_->ao_quadrupole());

    aodOBI->set_origin(origin_);
    aoqOBI->set_origin(origin_);

    // Compute multipole moment integrals
    aodOBI->compute(ao_Dpole);
    aoqOBI->compute(ao_Qpole);

    aodOBI.reset();
    aoqOBI.reset();

    std::vector<SharedVector> dipole;
    dipole.push_back(SharedVector(new Vector("Orbital Dipole X", Ca->ncol())));
    dipole.push_back(SharedVector(new Vector("Orbital Dipole Y", Ca->ncol())));
    dipole.push_back(SharedVector(new Vector("Orbital Dipole Z", Ca->ncol())));

    std::vector<SharedVector> quadrupole;
    quadrupole.push_back(SharedVector(new Vector("Orbital Quadrupole XX", Ca->ncol())));
    quadrupole.push_back(SharedVector(new Vector("Orbital Quadrupole YY", Ca->ncol())));
    quadrupole.push_back(SharedVector(new Vector("Orbital Quadrupole ZZ", Ca->ncol())));

    SharedMatrix temp(new Matrix("Temp", Ca->nrow(), Ca->ncol()));

    int nao = Ca->nrow();
    int nmo = Ca->ncol();

    if (same_orbs_) {

        // Dipoles
        C_DGEMM('T','N',nmo,nao,nao,1.0,Ca->pointer()[0],nmo,ao_Dpole[0]->pointer()[0],nao,0.0,temp->pointer()[0],nao);
        for (int i = 0; i < nmo; i++) {
            dipole[0]->set(0,i,C_DDOT(nao,Ca->pointer()[i],nmo,temp->pointer()[i],1));
        }
        C_DGEMM('T','N',nmo,nao,nao,1.0,Ca->pointer()[0],nmo,ao_Dpole[1]->pointer()[0],nao,0.0,temp->pointer()[0],nao);
        for (int i = 0; i < nmo; i++) {
            dipole[1]->set(0,i,C_DDOT(nao,Ca->pointer()[i],nmo,temp->pointer()[i],1));
        }
        C_DGEMM('T','N',nmo,nao,nao,1.0,Ca->pointer()[0],nmo,ao_Dpole[2]->pointer()[0],nao,0.0,temp->pointer()[0],nao);
        for (int i = 0; i < nmo; i++) {
            dipole[2]->set(0,i,C_DDOT(nao,Ca->pointer()[i],nmo,temp->pointer()[i],1));
        }
        // Quadrupoles
#if 0
        C_DGEMM('T','N',nmo,nao,nao,1.0,Ca->pointer()[0],nmo,ao_Qpole[0]->pointer()[0],nao,0.0,temp->pointer()[0],nao);
        for (int i = 0; i < nmo; i++) {
            quadrupole[0]->set(0,i,C_DDOT(nao,Ca->pointer()[i],nmo,temp->pointer()[i],1));
        }
        C_DGEMM('T','N',nmo,nao,nao,1.0,Ca->pointer()[0],nmo,ao_Qpole[3]->pointer()[0],nao,0.0,temp->pointer()[0],nao);
        for (int i = 0; i < nmo; i++) {
            quadrupole[1]->set(0,i,C_DDOT(nao,Ca->pointer()[i],nmo,temp->pointer()[i],1));
        }
        C_DGEMM('T','N',nmo,nao,nao,1.0,Ca->pointer()[0],nmo,ao_Qpole[5]->pointer()[0],nao,0.0,temp->pointer()[0],nao);
        for (int i = 0; i < nmo; i++) {
            quadrupole[2]->set(0,i,C_DDOT(nao,Ca->pointer()[i],nmo,temp->pointer()[i],1));
        }
#else
        for (int i=0; i<Ca->ncol(); ++i) {
            double sumx=0.0, sumy=0.0, sumz=0.0;
            for (int k=0; k<Ca->nrow(); ++k) {
                for (int l=0; l<Ca->nrow(); ++l) {
                    double tmp = Ca->get(0, k, i) * Ca->get(0, l, i);
                    sumx += ao_Qpole[0]->get(0, k, l) * tmp;
                    sumy += ao_Qpole[3]->get(0, k, l) * tmp;
                    sumz += ao_Qpole[5]->get(0, k, l) * tmp;
                }
            }

            quadrupole[0]->set(0, i, fabs(sumx));
            quadrupole[1]->set(0, i, fabs(sumy));
            quadrupole[2]->set(0, i, fabs(sumz));
        }
#endif
        char** labels = basisset_->molecule()->irrep_labels();
        std::vector<std::tuple<double,int,int> > metric;
        for (int h = 0; h < epsilon_a_->nirrep(); h++) {
            for (int i = 0; i < epsilon_a_->dimpi()[h]; i++) {
                metric.push_back(std::tuple<double,int,int>(epsilon_a_->get(h,i),i,h));
            }
        }
        std::sort(metric.begin(),metric.end());

        outfile->Printf( "\n  Orbital extents (a.u.):\n");
        outfile->Printf( "    %10s%15s%15s%15s%15s\n", "MO", "<x^2>", "<y^2>", "<z^2>", "<r^2>");

        for (int i = 0; i < nmo; i++) {
            int n = std::get<1>(metric[i]);
            int h = std::get<2>(metric[i]);

            double xx = quadrupole[0]->get(0, i),
                   yy = quadrupole[1]->get(0, i),
                   zz = quadrupole[2]->get(0, i);
            outfile->Printf( "    %4d%3s%3d%15.10f%15.10f%15.10f%15.10f\n",
                    i,
                    labels[h],
                    n,
                    fabs(quadrupole[0]->get(0, i)),
                    fabs(quadrupole[1]->get(0, i)),
                    fabs(quadrupole[2]->get(0, i)),
                    fabs(xx + yy + zz));
        }

        outfile->Printf( "\n");
        for(int h = 0; h < epsilon_a_->nirrep(); h++) free(labels[h]); free(labels);


    } else {

        // TODO: Both alpha and beta orbitals are reported separately
        // This helps identify symmetry breaking
    }
}

void OEProp::compute_mulliken_charges()
{
    outfile->Printf( "  Mulliken Charges: (a.u.)\n");

    std::shared_ptr<Molecule> mol = basisset_->molecule();
    std::shared_ptr<std::vector<double>> apcs(new std::vector<double>(mol->natom()));
    double* Qa = new double[mol->natom()];
    double* PSa = new double[basisset_->nbf()];
    double suma = 0.0;

    double* Qb = new double[mol->natom()];
    double* PSb = new double[basisset_->nbf()];
    double sumb = 0.0;

    ::memset(Qa, '\0', mol->natom()*sizeof(double));
    ::memset(Qb, '\0', mol->natom()*sizeof(double));

    SharedMatrix Da;
    SharedMatrix Db;

//    Get the Density Matrices for alpha and beta spins
    if (same_dens_) {
        Da = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
        Db = Da;
    } else {
        Da = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
        Db = wfn_->D_subset_helper(Db_so_, Cb_so_, "AO");
    }

//    Compute the overlap matrix

    std::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    SharedMatrix S(new Matrix("S",basisset_->nbf(),basisset_->nbf()));
    overlap->compute(S);

//    Form the idempotent D*S matrix

    SharedMatrix PSam(new Matrix("PSa",basisset_->nbf(),basisset_->nbf()));
    PSam->gemm(false,false,1.0,Da,S,0.0);
    SharedMatrix PSbm(new Matrix("PSb",basisset_->nbf(),basisset_->nbf()));
    PSbm->gemm(false,false,1.0,Db,S,0.0);

//     Accumulate registers

    for (int mu = 0; mu < basisset_->nbf(); mu++) {
        PSa[mu] = PSam->get(0,mu,mu);
        PSb[mu] = PSbm->get(0,mu,mu);

        int shell = basisset_->function_to_shell(mu);
        int A = basisset_->shell_to_center(shell);

        Qa[A] += PSa[mu];
        Qb[A] += PSb[mu];

        suma += PSa[mu];
        sumb += PSb[mu];
    }

//    Print out the Mulliken populations and charges

    outfile->Printf( "   Center  Symbol    Alpha    Beta     Spin     Total\n");
    double nuc = 0.0;
    for (int A = 0; A < mol->natom(); A++) {
        double Qs = Qa[A] - Qb[A];
        double Qt = mol->Z(A) - (Qa[A] + Qb[A]);
        (*apcs)[A]=Qt;
        outfile->Printf("   %5d    %2s    %8.5f %8.5f %8.5f %8.5f\n", A+1,mol->label(A).c_str(), \
            Qa[A], Qb[A], Qs, Qt);
        nuc += (double) mol->Z(A);
   }

    outfile->Printf( "\n   Total alpha = %8.5f, Total beta = %8.5f, Total charge = %8.5f\n", \
        suma, sumb, nuc - suma - sumb);
    wfn_->set_atomic_point_charges(apcs);
//    Free memory
    delete[] Qa;
    delete[] Qb;
    delete[] PSa;
    delete[] PSb;

    outfile->Printf( "\n");

}
void OEProp::compute_lowdin_charges()
{
    outfile->Printf( "  Lowdin Charges: (a.u.)\n");

    std::shared_ptr<Molecule> mol = basisset_->molecule();
    std::shared_ptr<std::vector<double>> apcs(new std::vector<double>(mol->natom()));
    double* Qa = new double[mol->natom()];
    double suma = 0.0;

    double* Qb = new double[mol->natom()];
    double sumb = 0.0;

    ::memset(Qa, '\0', mol->natom()*sizeof(double));
    ::memset(Qb, '\0', mol->natom()*sizeof(double));

    SharedMatrix Da;
    SharedMatrix Db;
    SharedMatrix evecs(new Matrix("Eigenvectors of S matrix",basisset_->nbf(),basisset_->nbf()));
    SharedMatrix temp(new Matrix("Temporary matrix",basisset_->nbf(),basisset_->nbf()));
    SharedMatrix SDSa(new Matrix("S_12 * D * S_12 alpha matrix",basisset_->nbf(),basisset_->nbf()));
    SharedMatrix SDSb(new Matrix("S_12 * D * S_12 beta matrix",basisset_->nbf(),basisset_->nbf()));
    std::shared_ptr<Vector> evals(new Vector(basisset_->nbf()));

//    Get the Density Matrices for alpha and beta spins
    if (same_dens_) {
        Da = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
        Db = Da;
    } else {
        Da = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
        Db = wfn_->D_subset_helper(Db_so_, Cb_so_, "AO");
    }

//    Compute the overlap matrix

    std::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    SharedMatrix S(new Matrix("S",basisset_->nbf(),basisset_->nbf()));
    overlap->compute(S);

//    Form the S^(1/2) matrix

    S->power(1.0/2.0);

//    Compute the S^(1/2)*D*S^(1/2) matrix

    temp->gemm(false,false,1.0,Da,S,0.0);
    SDSa->gemm(false,false,1.0,S,temp,0.0);
    temp->gemm(false,false,1.0,Db,S,0.0);
    SDSb->gemm(false,false,1.0,S,temp,0.0);

//    Accumulate AO populations for each atom

    for (int mu = 0; mu < basisset_->nbf(); mu++) {
        int shell = basisset_->function_to_shell(mu);
        int A = basisset_->shell_to_center(shell);

        Qa[A] += SDSa->get(0,mu,mu);
        Qb[A] += SDSb->get(0,mu,mu);

        suma += SDSa->get(0,mu,mu);
        sumb += SDSb->get(0,mu,mu);
    }

//    Print out the populations and charges

    outfile->Printf( "   Center  Symbol    Alpha    Beta     Spin     Total\n");
    double nuc = 0.0;
    for (int A = 0; A < mol->natom(); A++) {
        double Qs = Qa[A] - Qb[A];
        double Qt = mol->Z(A) - (Qa[A] + Qb[A]);
        (*apcs)[A]=Qt;
        outfile->Printf("   %5d    %2s    %8.5f %8.5f %8.5f %8.5f\n", A+1,mol->label(A).c_str(), \
            Qa[A], Qb[A], Qs, Qt);
        nuc += (double) mol->Z(A);
    }

    outfile->Printf( "\n  Total alpha = %8.5f, Total beta = %8.5f, Total charge = %8.5f\n", \
        suma, sumb, nuc - suma - sumb);
    wfn_->set_atomic_point_charges(apcs);
    delete[] Qa;
    delete[] Qb;


}
void OEProp::compute_mayer_indices()
{
    outfile->Printf( "\n\n  Mayer Bond Indices:\n\n");

    std::shared_ptr<Molecule> mol = basisset_->molecule();

    int nbf = basisset_->nbf();

    SharedMatrix Da;      // Density matrix for alpha spin
    SharedMatrix Db;      // Density matrix for beta spin

    SharedMatrix DSa(new Matrix("D * S alpha matrix",nbf,nbf));
    SharedMatrix DSb(new Matrix("D * S beta matrix",nbf,nbf));

//    Get the Density Matrices for alpha and beta spins
    if (same_dens_) {
        Da = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
        Db = Da;
    } else {
        Da = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
        Db = wfn_->D_subset_helper(Db_so_, Cb_so_, "AO");
    }

//    Compute the overlap matrix

    std::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    SharedMatrix S(new Matrix("S matrix",nbf,nbf));
    overlap->compute(S);

//    Form the idempotent D*S matrix

    DSa->gemm(false,false,1.0,Da,S,0.0);
    DSb->gemm(false,false,1.0,Db,S,0.0);

//     Compute Mayer bond indices

    int natom = mol->natom();

    SharedMatrix MBI_total(new Matrix(natom,natom));
    SharedMatrix MBI_alpha;
    SharedMatrix MBI_beta;

    if (!same_dens_) {
        MBI_alpha = SharedMatrix (new Matrix(natom,natom));
        MBI_beta = SharedMatrix (new Matrix(natom,natom));
    }

    for (int mu = 0; mu < nbf; mu++) {
        for (int nu = 0; nu < mu; nu++) {
            int shell_mu = basisset_->function_to_shell(mu);
            int shell_nu = basisset_->function_to_shell(nu);
            int atom_mu = basisset_->shell_to_center(shell_mu);
            int atom_nu = basisset_->shell_to_center(shell_nu);
            if (atom_mu == atom_nu) continue;

            double alpha = DSa->get(0, mu, nu) * DSa->get(0, nu, mu) ;
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

    std::shared_ptr<Vector> MBI_valence(new Vector(natom));

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
//    the MBI value will be underestimated. For example, the MBI value for the H–H bond of H2+
//    is calculated to be 0.25 using the NBO program. The equation coded above gives the correct value of 0.5.
//    For reference, see IJQC 29 (1986) P. 73 and IJQC 29 (1986) P. 477.

//    A nicer output is needed ...

    if (same_dens_) {
        MBI_total->print();
        outfile->Printf( "  Atomic Valences: \n");
        MBI_valence->print();
    }
    else {
        outfile->Printf( "  Total Bond Index: \n");
        MBI_total->print();
        outfile->Printf( "  Alpha Contribution: \n");
        MBI_alpha->print();
        outfile->Printf( "  Beta Contribution: \n");
        MBI_beta->print();
        outfile->Printf( "  Atomic Valences: \n");
        MBI_valence->print();
    }


}
void OEProp::compute_wiberg_lowdin_indices()
{
    outfile->Printf( "\n\n  Wiberg Bond Indices using Orthogonal Lowdin Orbitals:\n\n");

//    We may wanna get rid of these if we have NAOs...

    std::shared_ptr<Molecule> mol = basisset_->molecule();

    int nbf = basisset_->nbf();

    SharedMatrix Da;
    SharedMatrix Db;
    SharedMatrix evecs(new Matrix("Eigenvectors of S matrix",nbf,nbf));
    SharedMatrix temp(new Matrix("Temporary matrix",nbf,nbf));
    SharedMatrix SDSa(new Matrix("S_12 * D * S_12 alpha matrix",nbf,nbf));
    SharedMatrix SDSb(new Matrix("S_12 * D * S_12 beta matrix",nbf,nbf));
    std::shared_ptr<Vector> evals(new Vector(nbf));

//    Get the Density Matrices for alpha and beta spins
    if (same_dens_) {
        Da = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
        Db = Da;
    } else {
        Da = wfn_->D_subset_helper(Da_so_, Ca_so_, "AO");
        Db = wfn_->D_subset_helper(Db_so_, Cb_so_, "AO");
    }

//    Compute the overlap matrix

    std::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    SharedMatrix S(new Matrix("S",basisset_->nbf(),basisset_->nbf()));
    overlap->compute(S);

//    Form the S^(1/2) matrix

    S->diagonalize(evecs,evals);
    S->zero();
    for (int p = 0; p < basisset_->nbf(); ++p) S->set(0, p, p, sqrt(evals->get(0,p)));
    S->back_transform(evecs);

//    Compute the S^(1/2)*D*S^(1/2) matrix

    temp->gemm(false,false,1.0,Da,S,0.0);
    SDSa->gemm(false,false,1.0,S,temp,0.0);
    temp->gemm(false,false,1.0,Db,S,0.0);
    SDSb->gemm(false,false,1.0,S,temp,0.0);

//    Compute Wiberg bond indices

    int natom = mol->natom();

    SharedMatrix WBI_total(new Matrix(natom,natom));
    SharedMatrix WBI_alpha;
    SharedMatrix WBI_beta;

    if (!same_dens_) {
        WBI_alpha = SharedMatrix (new Matrix(natom,natom));
        WBI_beta = SharedMatrix (new Matrix(natom,natom));
    }

    for (int mu = 0; mu < nbf; mu++) {
        for (int nu = 0; nu < mu; nu++) {
            int shell_mu = basisset_->function_to_shell(mu);
            int shell_nu = basisset_->function_to_shell(nu);
            int atom_mu = basisset_->shell_to_center(shell_mu);
            int atom_nu = basisset_->shell_to_center(shell_nu);
            if (atom_mu == atom_nu) continue;

            double alpha = SDSa->get(0, mu, nu) * SDSa->get(0, nu, mu) ;
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

        std::shared_ptr<Vector> WBI_valence(new Vector(natom));

        for (int iat = 0; iat < natom; iat++) {
            for (int jat = 0; jat < natom; jat++) {
                double valence = WBI_valence->get(0, iat);
                WBI_valence->set(0, iat, valence + WBI_total->get(0, iat, jat));
            }
        }

//    Print out the bond index matrix
//    A nicer output is needed ...

    if (same_dens_) {
        WBI_total->print();
        outfile->Printf( "  Atomic Valences: \n");
        WBI_valence->print();
    }
    else {
        outfile->Printf( "  Total Bond Index: \n");
        WBI_total->print();
        outfile->Printf( "  Alpha Contribution: \n");
        WBI_alpha->print();
        outfile->Printf( "  Beta Contribution: \n");
        WBI_beta->print();
        outfile->Printf( "  Atomic Valences: \n");
        WBI_valence->print();
    }
}

void OEProp::compute_no_occupations()
{
    char** labels = basisset_->molecule()->irrep_labels();

    outfile->Printf( "  Natural Orbital Occupations:\n\n");

    if (!same_dens_) {

        SharedVector Oa;
        SharedVector Ob;
        if (same_dens_) {
            std::pair<SharedMatrix,SharedVector> vals = Na_mo();
            Oa = vals.second;
            Ob = vals.second;
        } else {
            std::pair<SharedMatrix,SharedVector> vals = Na_mo();
            Oa = vals.second;
            std::pair<SharedMatrix,SharedVector> vals2 = Nb_mo();
            Ob = vals2.second;
        }
        std::vector<std::tuple<double, int, int> > metric_a;
        for (int h = 0; h < Oa->nirrep(); h++) {
            for (int i = 0; i < Oa->dimpi()[h]; i++) {
                metric_a.push_back(std::tuple<double,int,int>(Oa->get(h,i), i ,h));
            }
        }

        std::sort(metric_a.begin(), metric_a.end(), std::greater<std::tuple<double,int,int> >());
        int offset_a = wfn_->nalpha();
        int start_occ_a = offset_a - max_noon_;
        start_occ_a = (start_occ_a < 0 ? 0 : start_occ_a);
        int stop_vir_a = offset_a + max_noon_ + 1;
        stop_vir_a = (int)((size_t)stop_vir_a >= metric_a.size() ? metric_a.size()  : stop_vir_a);

        outfile->Printf( "  Alpha Occupations:\n");
        for (int index = start_occ_a; index < stop_vir_a; index++) {
            if (index < offset_a) {
                outfile->Printf( "  HONO-%-2d: %4d%3s %8.3f\n", offset_a - index - 1,
                std::get<1>(metric_a[index])+1,labels[std::get<2>(metric_a[index])],
                std::get<0>(metric_a[index]));
            } else {
                outfile->Printf( "  LUNO+%-2d: %4d%3s %8.3f\n", index - offset_a,
                std::get<1>(metric_a[index])+1,labels[std::get<2>(metric_a[index])],
                std::get<0>(metric_a[index]));
            }
        }
        outfile->Printf( "\n");

        std::vector<std::tuple<double, int, int> > metric_b;
        for (int h = 0; h < Ob->nirrep(); h++) {
            for (int i = 0; i < Ob->dimpi()[h]; i++) {
                metric_b.push_back(std::tuple<double,int,int>(Ob->get(h,i), i ,h));
            }
        }

        std::sort(metric_b.begin(), metric_b.end(), std::greater<std::tuple<double,int,int> >());

        int offset_b = wfn_->nbeta();
        int start_occ_b = offset_b - max_noon_;
        start_occ_b = (start_occ_b < 0 ? 0 : start_occ_b);
        int stop_vir_b = offset_b + max_noon_ + 1;
        stop_vir_b = (int)((size_t)stop_vir_b >= metric_b.size() ? metric_b.size()  : stop_vir_b);

        outfile->Printf( "  Beta Occupations:\n");
        for (int index = start_occ_b; index < stop_vir_b; index++) {
            if (index < offset_b) {
                outfile->Printf( "  HONO-%-2d: %4d%3s %8.3f\n", offset_b - index - 1,
                std::get<1>(metric_b[index])+1,labels[std::get<2>(metric_b[index])],
                std::get<0>(metric_b[index]));
            } else {
                outfile->Printf( "  LUNO+%-2d: %4d%3s %8.3f\n", index - offset_b,
                std::get<1>(metric_b[index])+1,labels[std::get<2>(metric_b[index])],
                std::get<0>(metric_b[index]));
            }
        }
        outfile->Printf( "\n");

    }

    std::pair<SharedMatrix,SharedVector> vals = Nt_mo();
    SharedVector Ot = vals.second;

    std::vector<std::tuple<double, int, int> > metric;
    for (int h = 0; h < Ot->nirrep(); h++) {
        for (int i = 0; i < Ot->dimpi()[h]; i++) {
            metric.push_back(std::tuple<double,int,int>(Ot->get(h,i), i ,h));
        }
    }

    std::sort(metric.begin(), metric.end(), std::greater<std::tuple<double,int,int> >());

    int offset = wfn_->nbeta();
    int start_occ = offset - max_noon_;
    start_occ = (start_occ < 0 ? 0 : start_occ);
    int stop_vir = offset + max_noon_ + 1;
    stop_vir = (int)((size_t)stop_vir >= metric.size() ? metric.size()  : stop_vir);

    outfile->Printf( "  Total Occupations:\n");
    for (int index = start_occ; index < stop_vir; index++) {
        if (index < offset) {
            outfile->Printf( "  HONO-%-2d: %4d%3s %8.3f\n", offset - index - 1,
            std::get<1>(metric[index])+1,labels[std::get<2>(metric[index])],
            std::get<0>(metric[index]));
        } else {
            outfile->Printf( "  LUNO+%-2d: %4d%3s %8.3f\n", index - offset,
            std::get<1>(metric[index])+1,labels[std::get<2>(metric[index])],
            std::get<0>(metric[index]));
        }
    }
    outfile->Printf( "\n");

    //for(int h = 0; h < epsilon_a_->nirrep(); h++) free(labels[h]); free(labels);

}

//GridProp::GridProp(std::shared_ptr<Wavefunction> wfn) : filename_("out.grid"), Prop(wfn)
//{
//    common_init();
//}
//GridProp::GridProp() : filename_("out.grid"), Prop(Process::environment.legacy_wavefunction())
//{
//    common_init();
//}
//GridProp::~GridProp()
//{
//    reset();
//    free_block(temp_tens_);
//}
//void GridProp::common_init()
//{
//    initialized_ = false;
//    format_ = "DF3";
//
//    n_[0] = 40;
//    n_[1] = 40;
//    n_[2] = 40;
//
//    l_[0] = 5.0;
//    l_[1] = 5.0;
//    l_[2] = 5.0;
//
//    o_[0] = 0.0;
//    o_[1] = 0.0;
//    o_[2] = 0.0;
//
//    block_size_= 5000;
//
//    irrep_offsets_[0] = 0;
//    for (int h = 0; h < Ca_so_->nirrep() - 1; h++)
//        irrep_offsets_[h + 1] = irrep_offsets_[h] + Ca_so_->colspi()[h];
//
//    temp_tens_ = block_matrix(block_size_, basisset_->nbf());
//}
//void GridProp::add_alpha_mo(int irrep, int index)
//{
//    alpha_mos_.push_back(make_pair(irrep,index));
//}
//void GridProp::add_beta_mo(int irrep, int index)
//{
//    beta_mos_.push_back(make_pair(irrep,index));
//}
//void GridProp::add_basis_fun(int irrep, int index)
//{
//    basis_funs_.push_back(make_pair(irrep,index));
//}
//void GridProp::print_header()
//{
//    outfile->Printf( "\n GRIDPROP: One-electron grid properties.\n");
//    outfile->Printf( "  by Rob Parrish and Justin Turney.\n");
//    outfile->Printf( "  built on LIBMINTS.\n\n");
//}
//double*** GridProp::block_grid(int nx, int ny, int nz)
//{
//    double*** grid = new double**[nx];
//
//    double** pointers = new double*[nx*(unsigned long int)ny];
//
//    double* memory = new double[nx*(unsigned long int)ny*nz];
//    memset(static_cast<void*>(memory), '\0', sizeof(double)*nx*ny*nz);
//
//    for (int i = 0; i < nx; i++)
//        for (int j = 0; j < ny; j++)
//            pointers[i*(unsigned long int)ny + j] = &memory[i*(unsigned long int)ny*nz + j*(unsigned long int)nz];
//
//    for (int i = 0; i < nx; i++)
//        grid[i] = &pointers[i*(unsigned long int)ny];
//
//    return grid;
//}
//void GridProp::free_grid(double*** grid)
//{
//    delete[] grid[0][0];
//    delete[] grid[0];
//    delete[] grid;
//}
//void GridProp::build_grid_overages(double over)
//{
//    std::shared_ptr<Molecule> mol = basisset_->molecule();
//
//    double min_x = mol->x(0);
//    double min_y = mol->y(0);
//    double min_z = mol->z(0);
//    double max_x = mol->x(0);
//    double max_y = mol->y(0);
//    double max_z = mol->z(0);
//
//    for (int A = 0; A < mol->natom(); A++) {
//        if (mol->x(A) <= min_x)
//            min_x = mol->x(A);
//        if (mol->x(A) >= max_x)
//            max_x = mol->x(A);
//        if (mol->y(A) <= min_y)
//            min_y = mol->y(A);
//        if (mol->y(A) >= max_y)
//            max_y = mol->y(A);
//        if (mol->z(A) <= min_z)
//            min_z = mol->z(A);
//        if (mol->z(A) >= max_z)
//            max_z = mol->z(A);
//    }
//
//    min_x -= over;
//    min_y -= over;
//    min_z -= over;
//    max_x += over;
//    max_y += over;
//    max_z += over;
//
//    o_[0] = 0.5*(min_x + max_x);
//    o_[1] = 0.5*(min_y + max_y);
//    o_[2] = 0.5*(min_z + max_z);
//    l_[0] = (-min_x + max_x);
//    l_[1] = (-min_y + max_y);
//    l_[2] = (-min_z + max_z);
//
//    caxis_[0] = 0.0;
//    caxis_[1] = 1.0;
//
//    build_grid();
//}
//void GridProp::build_grid()
//{
//    int nx = n_[0] + 1;
//    int ny = n_[1] + 1;
//    int nz = n_[2] + 1;
//
//    grid_["x"] = block_grid(nx,ny,nz);
//    grid_["y"] = block_grid(nx,ny,nz);
//    grid_["z"] = block_grid(nx,ny,nz);
//
//    double* x = new double[nx];
//    double* y = new double[nx];
//    double* z = new double[nx];
//
//    double*** xg = grid_["x"];
//    double*** yg = grid_["y"];
//    double*** zg = grid_["z"];
//
//    if (nx == 0)
//        x[0] = 0.0;
//    else
//        for (int i = 0; i < nx; i++)
//            x[i] = ((double) i) / (double (nx - 1));
//    if (ny == 0)
//        y[0] = 0.0;
//    else
//        for (int i = 0; i < ny; i++)
//           y[i] = ((double) i) / (double (ny - 1));
//    if (nz == 0)
//        z[0] = 0.0;
//    else
//        for (int i = 0; i < nz; i++)
//           z[i] = ((double) i) / (double (nz - 1));
//
//    for (int i = 0; i < nx; i++) {
//        x[i] = l_[0] * (x[i] - 0.5) + o_[0];
//    }
//    for (int i = 0; i < ny; i++) {
//        y[i] = l_[1] * (y[i] - 0.5) + o_[1];
//    }
//    for (int i = 0; i < nz; i++) {
//        z[i] = l_[2] * (z[i] - 0.5) + o_[2];
//    }
//
//    for (int i = 0; i < nx; i++)
//        for (int j = 0; j < ny; j++)
//            for (int k = 0; k < nz; k++) {
//                xg[i][j][k] = x[i];
//                yg[i][j][k] = y[j];
//                zg[i][j][k] = z[k];
//            }
//
//    delete[] x;
//    delete[] y;
//    delete[] z;
//}
//void GridProp::allocate_arrays()
//{
//    int nx = n_[0] + 1;
//    int ny = n_[1] + 1;
//    int nz = n_[2] + 1;
//
//    for (std::set<std::string>::iterator it = tasks_.begin(); it != tasks_.end(); it++) {
//        if ((*it) == "MOS") {
//            // TODO
//        } else if ((*it) == "BASIS_FUNS") {
//            // Also TODO
//        } else {
//            grid_[(*it)] = block_grid(nx,ny,nz);
//        }
//    }
//}
//void GridProp::compute()
//{
//#if 0
//    reset();
//
//    initialized_ = true;
//
//    Da_ao_ = Da_ao();
//    Ca_ao_ = Ca_ao();
//    if (restricted_) {
//        Db_ao_ = Da_ao_;
//        Cb_ao_ = Ca_ao_;
//    } else {
//        Db_ao_ = Db_ao();
//        Cb_ao_ = Cb_ao();
//    }
//
//    print_header();
//    build_grid();
//    allocate_arrays();
//
//    int nx = n_[0] + 1;
//    int ny = n_[1] + 1;
//    int nz = n_[2] + 1;
//    ULI ngrid = nx*(ULI)ny*nz;
//    int nblock = ngrid / block_size_;
//    if (ngrid % block_size_ != 0)
//        nblock++;
//
//    double*** xp = grid_["x"];
//    double*** yp = grid_["y"];
//    double*** zp = grid_["z"];
//
//    // Basis points object (heavy lifting)
//    points_ = std::shared_ptr<BasisPoints>(new BasisPoints(basisset_, block_size_));
//    if (tasks_.count("GAMMA_AA") || tasks_.count("GAMMA_BB") || tasks_.count("GAMMA_AB") \
//        || tasks_.count("TAU_A") || tasks_.count("TAU_B"))
//        points_->setToComputeGradients(true);
//
//    // Grid block traversal object
//    std::shared_ptr<GridBlock> gridblock(new GridBlock());
//    gridblock->setMaxPoints(block_size_);
//
//    for (int block = 0; block < nblock; block++) {
//        // Indexing
//        int size = block_size_;
//        if (block*(ULI)block_size_ >= ngrid)
//            size = ngrid - block*(ULI)block_size_;
//
//        ULI offset = block*(ULI)block_size_;
//
//        // Line up gridblock pointers
//        // Last xp is a dirty hack b/c w is not needed for points
//        gridblock->setGrid(&xp[0][0][offset],&yp[0][0][offset],&zp[0][0][offset],&xp[0][0][offset]);
//        gridblock->setTruePoints(size);
//
//        // Compute basis functions/gradients
//        points_->computePoints(gridblock);
//
//        // Call compute routines
//        if (tasks_.count("MOS"))
//            compute_mos(gridblock, offset);
//        if (tasks_.count("BASIS_FUNS"))
//            compute_basis_funs(gridblock, offset);
//        if (tasks_.count("RHO"))
//            compute_rho(gridblock, &grid_["RHO"][0][0][offset]);
//        if (tasks_.count("RHO_S"))
//            compute_rho_s(gridblock, &grid_["RHO_S"][0][0][offset]);
//        if (tasks_.count("RHO_A"))
//            compute_rho_a(gridblock, &grid_["RHO_A"][0][0][offset]);
//        if (tasks_.count("RHO_B"))
//            compute_rho_b(gridblock, &grid_["RHO_B"][0][0][offset]);
//        if (tasks_.count("GAMMA_AA"))
//            compute_gamma_aa(gridblock, &grid_["GAMMA_AA"][0][0][offset]);
//        if (tasks_.count("GAMMA_AB"))
//            compute_gamma_ab(gridblock, &grid_["GAMMA_AB"][0][0][offset]);
//        if (tasks_.count("GAMMA_BB"))
//            compute_gamma_bb(gridblock, &grid_["GAMMA_BB"][0][0][offset]);
//        if (tasks_.count("TAU_A"))
//            compute_rho_b(gridblock, &grid_["TAU_A"][0][0][offset]);
//        if (tasks_.count("TAU_B"))
//            compute_rho_b(gridblock, &grid_["TAU_B"][0][0][offset]);
//    }
//
//    // ESP is special we think
//    if (tasks_.count("ESP"))
//        compute_ESP();
//
//    if (format_ == "DF3")
//        write_df3_grid();
//    else
//        write_data_grid();
//
//#endif
//}
//void GridProp::write_data_grid()
//{
//    int nx = n_[0] + 1;
//    int ny = n_[1] + 1;
//    int nz = n_[2] + 1;
//
//    for (std::map<std::string, double***>::iterator it = grid_.begin(); it != grid_.end(); it++) {
//        std::string key = (*it).first;
//        double*** data = (*it).second;
//
//        /* Write it to a file */
//        int i,j,k;
//        std::string file = filename_ + "." + key + ".dat";
//        FILE* fptr = fopen(file.c_str(),"w");
//        outfile->Printf(fptr,"%d %d %d\n\n", nx,ny,nz);
//        for (k=0;k<nz;k++) {
//           for (j=0;j<ny;j++) {
//              for (i=0;i<nx;i++) {
//                    outfile->Printf(fptr,"%24.16f ", data[i][j][k]);
//                }
//            outfile->Printf(fptr,"\n");
//           }
//           outfile->Printf(fptr,"\n");
//        }
//        fclose(fptr);
//    }
//}
//void GridProp::write_df3_grid()
//{
//    int nx = n_[0] + 1;
//    int ny = n_[1] + 1;
//    int nz = n_[2] + 1;
//
//    for (std::map<std::string, double***>::iterator it = grid_.begin(); it != grid_.end(); it++) {
//        std::string key = (*it).first;
//        double*** data = (*it).second;
//
//        double v;
//        double themin = data[0][0][0];
//        double themax = data[0][0][0];
//
//        /* Write it to a file */
//        std::string file = filename_ + "." + key + ".df3";
//        FILE* fptr = fopen(file.c_str(),"w");
//        fputc(nx >> 8,fptr);
//        fputc(nx & 0xff,fptr);
//        fputc(ny >> 8,fptr);
//        fputc(ny & 0xff,fptr);
//        fputc(nz >> 8,fptr);
//        fputc(nz & 0xff,fptr);
//        int i,j,k;
//        for (k=0;k<nz;k++) {
//           for (j=0;j<ny;j++) {
//              for (i=0;i<nx;i++) {
//                 if (data[i][j][k] > caxis_[1] )
//                    v = 255;
//                 else if (data[i][j][k] < caxis_[0])
//                    v = 0;
//                 else
//                    v = 255 * (data[i][j][k]-caxis_[0])/(caxis_[1] - caxis_[0]);
//                 fputc((int)v,fptr);
//              }
//           }
//        }
//        fclose(fptr);
//    }
//}
//void GridProp::reset()
//{
//    if (!initialized_)
//        return;
//
//    // Free the points object
//    //points_.reset();
//
//    // Free the grids
//    for (std::map<std::string, double***>::iterator it = grid_.begin(); it != grid_.end(); ++it) {
//        free_grid((*it).second);
//    }
//}
//#if 0
//void GridProp::compute_mos(std::shared_ptr<GridBlock> g, ULI offset)
//{
//    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
//}
//void GridProp::compute_basis_funs(std::shared_ptr<GridBlock> g, ULI offset)
//{
//    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
//}
//void GridProp::compute_rho(std::shared_ptr<GridBlock> g, double* results)
//{
//    int npoints = g->getTruePoints();
//    int nbf = basisset_->nbf();
//    double** points = points_->getPoints();
//    double** Da = Da_ao_->pointer();
//    double** Db = Db_ao_->pointer();
//
//    // rho_a_
//    // rho_a^Q = phi_m^Q * Da_mn * phi_n^Q
//    C_DGEMM('N', 'N', npoints, nbf, nbf, 1.0, &points[0][0], nbf, &Da[0][0], nbf, \
//        0.0, &temp_tens_[0][0], nbf);
//
//    for (int Q = 0; Q < npoints; Q++) {
//        results[Q] = C_DDOT(nbf, &temp_tens_[Q][0], 1, &points[Q][0], 1);
//        //printf(" Q = %d, rho = %14.10E\n", Q, rho_a_[Q]);
//    }
//
//    if (!restricted_) {
//
//        // rho_b^Q = phi_m^Q * Db_mn * phi_n^Q
//        C_DGEMM('N', 'N', npoints, nbf, nbf, 1.0, &points[0][0], nbf, &Db[0][0], nbf, \
//            0.0, &temp_tens_[0][0], nbf);
//
//        for (int Q = 0; Q < npoints; Q++) {
//            results[Q] += C_DDOT(nbf, &temp_tens_[Q][0], 1, &points[Q][0], 1);
//            //printf(" Q = %d, rho = %14.10E\n", Q, rho_b_[Q]);
//        }
//
//    } else {
//        C_DSCAL(npoints,2.0,results,1);
//    }
//}
//void GridProp::compute_rho_s(std::shared_ptr<GridBlock> g, double* results)
//{
//    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
//}
//void GridProp::compute_rho_a(std::shared_ptr<GridBlock> g, double* results)
//{
//    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
//}
//void GridProp::compute_rho_b(std::shared_ptr<GridBlock> g, double* results)
//{
//    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
//}
//void GridProp::compute_gamma_aa(std::shared_ptr<GridBlock> g, double* results)
//{
//    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
//}
//void GridProp::compute_gamma_ab(std::shared_ptr<GridBlock> g, double* results)
//{
//    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
//}
//void GridProp::compute_gamma_bb(std::shared_ptr<GridBlock> g, double* results)
//{
//    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
//}
//void GridProp::compute_tau_a(std::shared_ptr<GridBlock> g, double* results)
//{
//    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
//}
//void GridProp::compute_tau_b(std::shared_ptr<GridBlock> g, double* results)
//{
//    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
//}
//#endif
//void GridProp::compute_ESP()
//{
//    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
//}

}
