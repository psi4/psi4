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

#include "fisapt.h"

#include <algorithm>
#include <ctime>
#include <functional>
#include <set>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/psi4-dec.h"
#include "psi4/physconst.h"

#include "psi4/lib3index/dfhelper.h"
#include "psi4/libcubeprop/csg.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libfock/jk.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/extern.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/potential.h"
#include "psi4/libmints/vector.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/mintshelper.h"

#include "local2.h"

namespace psi {

namespace fisapt {

FISAPT::FISAPT(SharedWavefunction scf) : options_(Process::environment.options), reference_(scf) { common_init(); }
FISAPT::FISAPT(SharedWavefunction scf, Options& options) : options_(options), reference_(scf) { common_init(); }
FISAPT::~FISAPT() {}

void FISAPT::common_init() {
    reference_->set_module("fisapt");

    primary_ = reference_->basisset();
    doubles_ = Process::environment.get_memory() / sizeof(double) * options_.get_double("FISAPT_MEM_SAFETY_FACTOR");

    matrices_["Cfocc"] = reference_->Ca_subset("AO", "FROZEN_OCC");

    vectors_["eps_all"] = reference_->epsilon_a_subset("AO", "ALL");

    matrices_["Call"] = reference_->Ca_subset("AO", "ALL");
    matrices_["Cocc"] = reference_->Ca_subset("AO", "OCC");
    matrices_["Cvir"] = reference_->Ca_subset("AO", "VIR");

    vectors_["eps_occ"] = reference_->epsilon_a_subset("AO", "OCC");
    vectors_["eps_vir"] = reference_->epsilon_a_subset("AO", "VIR");

    matrices_["Caocc"] = reference_->Ca_subset("AO", "ACTIVE_OCC");
    matrices_["Cavir"] = reference_->Ca_subset("AO", "ACTIVE_VIR");
    matrices_["Cfvir"] = reference_->Ca_subset("AO", "FROZEN_VIR");

    vectors_["eps_focc"] = reference_->epsilon_a_subset("AO", "FROZEN_OCC");
    vectors_["eps_aocc"] = reference_->epsilon_a_subset("AO", "ACTIVE_OCC");
    vectors_["eps_avir"] = reference_->epsilon_a_subset("AO", "ACTIVE_VIR");
    vectors_["eps_fvir"] = reference_->epsilon_a_subset("AO", "FROZEN_VIR");
}

void FISAPT::print_header() {
    outfile->Printf("\t --------------------------------------------\n");
    outfile->Printf("\t                    FISAPT0                  \n");
    outfile->Printf("\t                  Rob Parrish                \n");
    outfile->Printf("\t --------------------------------------------\n");
    outfile->Printf("\n");

    outfile->Printf("    Do F-SAPT = %11s\n", options_.get_bool("FISAPT_DO_FSAPT") ? "Yes" : "No");
    outfile->Printf("    Do Plot   = %11s\n", options_.get_bool("FISAPT_DO_PLOT") ? "Yes" : "No");
    outfile->Printf("    Memory    = %11.3f [GiB]\n", (doubles_ * 8) / (1024.0 * 1024.0 * 1024.0));
    outfile->Printf("\n");
}

void FISAPT::localize() {

    outfile->Printf("  ==> Localization (IBO) <==\n\n");

    std::shared_ptr<Matrix> Focc =
        std::make_shared<Matrix>("Focc", vectors_["eps_occ"]->dimpi()[0], vectors_["eps_occ"]->dimpi()[0]);
    Focc->set_diagonal(vectors_["eps_occ"]);

    std::vector<int> ranges;
    ranges.push_back(0);
    ranges.push_back(vectors_["eps_focc"]->dimpi()[0]);
    ranges.push_back(vectors_["eps_occ"]->dimpi()[0]);

    std::shared_ptr<fisapt::IBOLocalizer2> local =
        fisapt::IBOLocalizer2::build(primary_, reference_->get_basisset("MINAO"), matrices_["Cocc"], options_);
    local->print_header();
    std::map<std::string, std::shared_ptr<Matrix> > ret = local->localize(matrices_["Cocc"], Focc, ranges);

    matrices_["Locc"] = ret["L"];
    matrices_["Qocc"] = ret["Q"];
    matrices_["IAO"] = ret["A"];
}

void FISAPT::partition() {
    outfile->Printf("  ==> Partitioning <==\n\n");

    // => Sizing <= //

    std::shared_ptr<Molecule> mol = primary_->molecule();
    int nA = mol->natom();
    int na = matrices_["Locc"]->colspi()[0];

    // => Monomer Atoms <= //

    const std::vector<std::pair<int, int> >& fragment_list = mol->get_fragments();
    if (!(fragment_list.size() == 2 || fragment_list.size() == 3)) {
        throw PSIEXCEPTION("FISAPT: Molecular system must have 2 (A+B) or 3 (A+B+C) fragments");
    }

    std::vector<int> indA;
    std::vector<int> indB;
    std::vector<int> indC;

    for (int ind = fragment_list[0].first; ind < fragment_list[0].second; ind++) {
        indA.push_back(ind);
    }
    for (int ind = fragment_list[1].first; ind < fragment_list[1].second; ind++) {
        indB.push_back(ind);
    }
    if (fragment_list.size() == 3) {
        for (int ind = fragment_list[2].first; ind < fragment_list[2].second; ind++) {
            indC.push_back(ind);
        }
    }

    // => Lookup table for atoms (1-A, 2-B, 3-C) <= //
    vectors_["FRAG"] = std::make_shared<Vector>("FRAG", nA);
    double* FRAGp = vectors_["FRAG"]->pointer();
    for (int ind = 0; ind < indA.size(); ind++) {
        FRAGp[indA[ind]] = 1;
    }
    for (int ind = 0; ind < indB.size(); ind++) {
        FRAGp[indB[ind]] = 2;
    }
    for (int ind = 0; ind < indC.size(); ind++) {
        FRAGp[indC[ind]] = 3;
    }
    outfile->Printf(" Fragment lookup table:\n");
    vectors_["FRAG"]->print();

    outfile->Printf("   => Atomic Partitioning <= \n\n");
    outfile->Printf("    Monomer A: %3zu atoms\n", indA.size());
    outfile->Printf("    Monomer B: %3zu atoms\n", indB.size());
    outfile->Printf("    Monomer C: %3zu atoms\n", indC.size());
    outfile->Printf("\n");

    // => Fragment Orbital Charges <= //

    auto QF = std::make_shared<Matrix>("QF", 3, na);
    double** QFp = QF->pointer();
    double** Qp = matrices_["Qocc"]->pointer();

    for (int ind = 0; ind < indA.size(); ind++) {
        for (int a = 0; a < na; a++) {
            QFp[0][a] += Qp[indA[ind]][a];
        }
    }

    for (int ind = 0; ind < indB.size(); ind++) {
        for (int a = 0; a < na; a++) {
            QFp[1][a] += Qp[indB[ind]][a];
        }
    }

    for (int ind = 0; ind < indC.size(); ind++) {
        for (int a = 0; a < na; a++) {
            QFp[2][a] += Qp[indC[ind]][a];
        }
    }

    // => Link Identification <= //

    std::string link_selection = options_.get_str("FISAPT_LINK_SELECTION");
    outfile->Printf("   => Link Bond Identification <=\n\n");
    outfile->Printf("    Link Bond Selection = %s\n\n", link_selection.c_str());

    std::vector<int> link_orbs;
    std::vector<std::pair<int, int> > link_atoms;
    std::vector<std::string> link_types;

    if (link_selection == "AUTOMATIC") {
        double delta = options_.get_double("FISAPT_CHARGE_COMPLETENESS");
        outfile->Printf("    Charge Completeness = %5.3f\n\n", delta);
        for (int a = 0; a < na; a++) {
            if (QFp[0][a] > delta || QFp[1][a] > delta || QFp[2][a] > delta) {
                ;
            } else if (QFp[0][a] + QFp[2][a] > delta) {
                link_orbs.push_back(a);
                link_types.push_back("AC");
            } else if (QFp[1][a] + QFp[2][a] > delta) {
                link_orbs.push_back(a);
                link_types.push_back("BC");
            } else if (QFp[0][a] + QFp[1][a] > delta) {
                link_orbs.push_back(a);
                link_types.push_back("AB");
            } else if (QFp[0][a] + QFp[1][a] > delta) {
                ;
            } else {
                throw PSIEXCEPTION("FISAPT: A, B, and C are bonded?! 3c-2e bonds are not cool.");
            }
        }

        for (int ind = 0; ind < link_orbs.size(); ind++) {
            int a = link_orbs[ind];
            std::vector<std::pair<double, int> > Qvals;
            for (int A = 0; A < nA; A++) {
                Qvals.push_back(std::pair<double, int>(Qp[A][a], A));
            }
            std::sort(Qvals.begin(), Qvals.end(), std::greater<std::pair<double, int> >());
            int A1 = Qvals[0].second;
            int A2 = Qvals[1].second;
            if (A2 < A1) std::swap(A1, A2);
            link_atoms.push_back(std::pair<int, int>(A1, A2));
        }

    } else if (link_selection == "MANUAL") {
        for (int ind = 0; ind < options_["FISAPT_MANUAL_LINKS"].size(); ind++) {
            link_atoms.push_back(std::pair<int, int>(options_["FISAPT_MANUAL_LINKS"][ind][0].to_integer() - 1,
                                                     options_["FISAPT_MANUAL_LINKS"][ind][1].to_integer() - 1));
        }

        for (int ind = 0; ind < link_atoms.size(); ind++) {
            int A1 = link_atoms[ind].first;
            int A2 = link_atoms[ind].second;
            if (A2 < A1) std::swap(A1, A2);

            double Qmax = 0.0;
            int aind = -1;
            for (int a = 0; a < na; a++) {
                double Qval = Qp[A1][a] * Qp[A2][a];
                if (Qmax < Qval) {
                    Qmax = Qval;
                    aind = a;
                }
            }
            link_orbs.push_back(aind);

            if (std::find(indA.begin(), indA.end(), A1) != indA.end() &&
                std::find(indC.begin(), indC.end(), A2) != indC.end()) {
                link_types.push_back("AC");
            } else if (std::find(indB.begin(), indB.end(), A1) != indB.end() &&
                       std::find(indC.begin(), indC.end(), A2) != indC.end()) {
                link_types.push_back("BC");
            } else if (std::find(indA.begin(), indA.end(), A1) != indA.end() &&
                       std::find(indB.begin(), indB.end(), A2) != indB.end()) {
                link_types.push_back("AB");
            } else {
                throw PSIEXCEPTION("FISAPT: FISAPT_MANUAL_LINKS contains a bond which is not AB, AC, or BC");
            }
        }

    } else {
        throw PSIEXCEPTION("FISAPT: Unrecognized FISAPT_LINK_SELECTION option.");
    }

    outfile->Printf("    Total Link Bonds = %zu\n\n", link_orbs.size());

    if (link_orbs.size()) {
        outfile->Printf("    --------------------------\n");
        outfile->Printf("    %-4s %4s %4s %5s %5s\n", "N", "Orb", "Type", "Aind1", "Aind2");
        outfile->Printf("    --------------------------\n");
        for (int ind = 0; ind < link_orbs.size(); ind++) {
            outfile->Printf("    %-4d %4d %4s %5d %5d\n", ind + 1, link_orbs[ind], link_types[ind].c_str(),
                            link_atoms[ind].first + 1, link_atoms[ind].second + 1);
        }
        outfile->Printf("    --------------------------\n");
        outfile->Printf("\n");
    }

    // => Nuclear Charge Targets <= //

    vectors_["ZA"] = std::make_shared<Vector>("ZA", nA);
    vectors_["ZB"] = std::make_shared<Vector>("ZB", nA);
    vectors_["ZC"] = std::make_shared<Vector>("ZC", nA);
    vectors_["ZA_orig"] = std::make_shared<Vector>("ZA_orig", nA);
    vectors_["ZB_orig"] = std::make_shared<Vector>("ZB_orig", nA);
    vectors_["ZC_orig"] = std::make_shared<Vector>("ZC_orig", nA);

    double* ZAp = vectors_["ZA"]->pointer();
    double* ZBp = vectors_["ZB"]->pointer();
    double* ZCp = vectors_["ZC"]->pointer();
    double* ZA_origp = vectors_["ZA_orig"]->pointer();
    double* ZB_origp = vectors_["ZB_orig"]->pointer();
    double* ZC_origp = vectors_["ZC_orig"]->pointer();

    for (int ind = 0; ind < indA.size(); ind++) {
        ZAp[indA[ind]] = mol->Z(indA[ind]);
    }
    for (int ind = 0; ind < indB.size(); ind++) {
        ZBp[indB[ind]] = mol->Z(indB[ind]);
    }
    for (int ind = 0; ind < indC.size(); ind++) {
        ZCp[indC[ind]] = mol->Z(indC[ind]);
    }

    vectors_["ZA_orig"]->copy(*vectors_["ZA"]);
    vectors_["ZB_orig"]->copy(*vectors_["ZB"]);
    vectors_["ZC_orig"]->copy(*vectors_["ZC"]);

    // => Local Orbital Targets <= //

    std::vector<int> orbsA;
    std::vector<int> orbsB;
    std::vector<int> orbsC;
    std::vector<int> orbsL;

    // => Assign Links <= //

    std::string link_assignment = options_.get_str("FISAPT_LINK_ASSIGNMENT");
    if (!(link_assignment == "AB" || link_assignment == "C" || link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2" )) {
        throw PSIEXCEPTION("FISAPT: FISAPT_LINK_ASSIGNMENT not recognized");
    }

    outfile->Printf("   => Link Bond Assignment <=\n\n");
    outfile->Printf("    Link Bond Assignment      = %s\n", link_assignment.c_str());
    outfile->Printf("\n");

    if (link_assignment == "C" || link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") {
        for (int ind = 0; ind < link_orbs.size(); ind++) {
            int a = link_orbs[ind];
            int A1 = link_atoms[ind].first;
            int A2 = link_atoms[ind].second;
            std::string type = link_types[ind];

            typesL.push_back(type);
            if (type == "AC") {
                ZAp[A1] -= 1.0;
                ZCp[A1] += 1.0;
                orbsC.push_back(a);
                orbsL.push_back(a);
            } else if (type == "BC") {
                ZBp[A1] -= 1.0;
                ZCp[A1] += 1.0;
                orbsC.push_back(a);
                orbsL.push_back(a);
            } else if (type == "AB") {
                ZAp[A1] -= 1.0;
                ZCp[A1] += 1.0;
                ZBp[A2] -= 1.0;
                ZCp[A2] += 1.0;
                orbsC.push_back(a);
                orbsL.push_back(a);
            }
        }

    } else if (link_assignment == "AB") {
        for (int ind = 0; ind < link_orbs.size(); ind++) {
            int a = link_orbs[ind];
            int A1 = link_atoms[ind].first;
            int A2 = link_atoms[ind].second;
            std::string type = link_types[ind];

            if (type == "AC") {
                ZAp[A1] += 1.0;
                ZCp[A1] -= 1.0;
                orbsA.push_back(a);
            } else if (type == "BC") {
                ZBp[A1] += 1.0;
                ZCp[A1] -= 1.0;
                orbsB.push_back(a);
            } else if (type == "AB") {
                throw PSIEXCEPTION("FISAPT: AB link requires LINK_ASSIGNMENT C");
            }
        }
    }
//We altered the nuclear charges for the C algorithm in vectors_["ZA"]-vectors_["ZC"].
//We keep the unaltered charges for the SAOn/SIAOn algorithms in vectors_["ZA_orig"]-vectors_["ZC_orig"].

    // => Remaining Orbitals <= //

    const std::vector<int>& fragment_charges = mol->get_fragment_charges();
    int CA2 = fragment_charges[0];
    int CB2 = fragment_charges[1];
    int CC2 = (fragment_charges.size() == 3 ? fragment_charges[2] : 0);

    double ZA2 = 0.0;
    double ZB2 = 0.0;
    double ZC2 = 0.0;
    for (int ind = 0; ind < nA; ind++) {
        ZA2 += ZAp[ind];
        ZB2 += ZBp[ind];
        ZC2 += ZCp[ind];
    }

    int EA2 = round(ZA2) - CA2;  // Number of needed electrons
    int EB2 = round(ZB2) - CB2;  // Number of needed electrons
    int EC2 = round(ZC2) - CC2;  // Number of needed electrons

    if (EA2 % 2 != 0) throw PSIEXCEPTION("FISAPT: Charge on A is incompatible with singlet");
    if (EB2 % 2 != 0) throw PSIEXCEPTION("FISAPT: Charge on B is incompatible with singlet");
    if (EC2 % 2 != 0) throw PSIEXCEPTION("FISAPT: Charge on C is incompatible with singlet");

    int NA2 = EA2 / 2;
    int NB2 = EB2 / 2;
    int NC2 = EC2 / 2;

    if (NA2 + NB2 + NC2 != na)
        throw PSIEXCEPTION("FISAPT: Sum of charges is incompatible with total number of electrons.");

    int RA2 = NA2 - orbsA.size();
    int RB2 = NB2 - orbsB.size();
    int RC2 = NC2 - orbsC.size();

    std::set<int> taken_orbs;
    for (int x : orbsA) taken_orbs.insert(x);
    for (int x : orbsB) taken_orbs.insert(x);
    for (int x : orbsC) taken_orbs.insert(x);

    std::vector<std::pair<double, int> > QCvals;
    for (int a = 0; a < na; a++) {
        if (taken_orbs.count(a)) continue;
        QCvals.push_back(std::pair<double, int>(QFp[2][a], a));
    }
    std::sort(QCvals.begin(), QCvals.end(), std::greater<std::pair<double, int> >());
    for (int ind = 0; ind < RC2; ind++) {
        int a = QCvals[ind].second;
        orbsC.push_back(a);
        taken_orbs.insert(a);
    }

    std::vector<std::pair<double, int> > QAvals;
    for (int a = 0; a < na; a++) {
        if (taken_orbs.count(a)) continue;
        QAvals.push_back(std::pair<double, int>(QFp[0][a], a));
    }
    std::sort(QAvals.begin(), QAvals.end(), std::greater<std::pair<double, int> >());
    for (int ind = 0; ind < RA2; ind++) {
        int a = QAvals[ind].second;
        orbsA.push_back(a);
        taken_orbs.insert(a);
    }

    std::vector<std::pair<double, int> > QBvals;
    for (int a = 0; a < na; a++) {
        if (taken_orbs.count(a)) continue;
        QBvals.push_back(std::pair<double, int>(QFp[1][a], a));
    }
    std::sort(QBvals.begin(), QBvals.end(), std::greater<std::pair<double, int> >());
    for (int ind = 0; ind < RB2; ind++) {
        int a = QBvals[ind].second;
        orbsB.push_back(a);
        taken_orbs.insert(a);
    }

    // => Orbital Subsets <= //

    std::sort(orbsA.begin(), orbsA.end());
    std::sort(orbsB.begin(), orbsB.end());
    std::sort(orbsC.begin(), orbsC.end());
    if (link_orbs.size() > 1) {
        if (orbsL[0] > orbsL[1]) {
            int orbswap = orbsL[1];
            orbsL[1] = orbsL[0];
            orbsL[0] = orbswap;
            std::string typeswap = typesL[1];
            typesL[1] = typesL[0];
            typesL[0] = typeswap;
        }
    }

    // We save matrices with the localized orbital (IBO) coefficients for A, B, C, and the two link orbitals only.
    // We need the latter for later link reassignment in the SAOn/SIAOn methods.
    matrices_["LoccA"] = FISAPT::extract_columns(orbsA, matrices_["Locc"]);
    matrices_["LoccB"] = FISAPT::extract_columns(orbsB, matrices_["Locc"]);
    matrices_["LoccC"] = FISAPT::extract_columns(orbsC, matrices_["Locc"]);
    matrices_["LoccL"] = FISAPT::extract_columns(orbsL, matrices_["Locc"]);

    matrices_["LoccA"]->set_name("LoccA");
    matrices_["LoccB"]->set_name("LoccB");
    matrices_["LoccC"]->set_name("LoccC");
    matrices_["LoccL"]->set_name("LoccL");

    // matrices_["LoccA"]->print();
    // matrices_["LoccB"]->print();
    // matrices_["LoccC"]->print();

    // => Summary <= //

    double ZA = 0.0;
    double ZB = 0.0;
    double ZC = 0.0;
    for (int ind = 0; ind < nA; ind++) {
        ZA += ZAp[ind];
        ZB += ZBp[ind];
        ZC += ZCp[ind];
    }
    ZA = round(ZA);
    ZB = round(ZB);
    ZC = round(ZC);

    double YA = round(2.0 * orbsA.size());
    double YB = round(2.0 * orbsB.size());
    double YC = round(2.0 * orbsC.size());

    outfile->Printf("   => Partition Summary <=\n\n");
    outfile->Printf("    Monomer A: %2d charge, %3d protons, %3d electrons, %3zu docc\n", (int)(ZA - YA), (int)ZA,
                    (int)YA, orbsA.size());
    outfile->Printf("    Monomer B: %2d charge, %3d protons, %3d electrons, %3zu docc\n", (int)(ZB - YB), (int)ZB,
                    (int)YB, orbsB.size());
    outfile->Printf("    Monomer C: %2d charge, %3d protons, %3d electrons, %3zu docc\n", (int)(ZC - YC), (int)ZC,
                    (int)YC, orbsC.size());
    outfile->Printf("\n");
}

void FISAPT::overlap() {
    outfile->Printf("  ==> Overlap Integrals <==\n\n");

    int nm = primary_->nbf();
    auto Tfact = std::make_shared<IntegralFactory>(primary_);
    std::shared_ptr<OneBodyAOInt> Tint = std::shared_ptr<OneBodyAOInt>(Tfact->ao_overlap());
    matrices_["S"] = std::make_shared<Matrix>("S", nm, nm);
    Tint->compute(matrices_["S"]);
}

void FISAPT::kinetic() {
    outfile->Printf("  ==> Kinetic Integrals <==\n\n");

    int nm = primary_->nbf();
    auto Tfact = std::make_shared<IntegralFactory>(primary_);
    std::shared_ptr<OneBodyAOInt> Tint = std::shared_ptr<OneBodyAOInt>(Tfact->ao_kinetic());
    matrices_["T"] = std::make_shared<Matrix>("T", nm, nm);
    Tint->compute(matrices_["T"]);
}

// void sapt_test(){
//     printf("Hello from sapt_test\n");
// }
// void sapt_test(){
//     printf("Hello from sapt_test\n");
// }


void FISAPT::nuclear() {
    outfile->Printf("  ==> Nuclear Integrals <==\n\n");

    // => Sizing <= //

    std::shared_ptr<Molecule> mol = primary_->molecule();
    int nA = mol->natom();
    int nm = primary_->nbf();

    // => Nuclear Potentials <= //

    auto Vfact = std::make_shared<IntegralFactory>(primary_);
    std::shared_ptr<PotentialInt> Vint;
    Vint = std::shared_ptr<PotentialInt>(static_cast<PotentialInt*>(Vfact->ao_potential().release()));

    // => Dipole Integrals, Needed for Testing the Dipole Magnitude <= //

    std::shared_ptr<OneBodyAOInt> dipole(Vfact->ao_dipole());
    std::vector<SharedMatrix> dipole_ints;
    dipole_ints.push_back(std::make_shared<Matrix>("Dipole X", nm, nm));
    dipole_ints.push_back(std::make_shared<Matrix>("Dipole Y", nm, nm));
    dipole_ints.push_back(std::make_shared<Matrix>("Dipole Z", nm, nm));
    dipole->compute(dipole_ints);
    matrices_["X DIP"] = dipole_ints[0];
    matrices_["Y DIP"] = dipole_ints[1];
    matrices_["Z DIP"] = dipole_ints[2];
    auto dip_nuca = std::make_shared<Vector>("DIP NUCA", 3);
    auto dip_nucb = std::make_shared<Vector>("DIP NUCB", 3);
    auto dip_nucc = std::make_shared<Vector>("DIP NUCC", 3);
    double* dip_nucap = dip_nuca->pointer();
    double* dip_nucbp = dip_nucb->pointer();
    double* dip_nuccp = dip_nucc->pointer();

    // > Molecular Centers < //

    std::vector<std::pair<double, std::array<double, 3>>> Zxyz;

    // > A < //
    double* ZAp = vectors_["ZA"]->pointer();
    for (int A = 0; A < nA; A++) {
        Zxyz.push_back({ZAp[A], {{mol->x(A), mol->y(A), mol->z(A)}}});
    }
    Vint->set_charge_field(Zxyz);

    matrices_["VA"] = std::make_shared<Matrix>("VA", nm, nm);
    Vint->compute(matrices_["VA"]);

    // > B < //

    double* ZBp = vectors_["ZB"]->pointer();
    for (int A = 0; A < nA; A++) {
        Zxyz[A].first = ZBp[A];
    }
    Vint->set_charge_field(Zxyz);

    matrices_["VB"] = std::make_shared<Matrix>("VB", nm, nm);
    Vint->compute(matrices_["VB"]);

    // > C < //

    double* ZCp = vectors_["ZC"]->pointer();
    for (int A = 0; A < nA; A++) {
        Zxyz[A].first = ZCp[A];
    }
    Vint->set_charge_field(Zxyz);

    matrices_["VC"] = std::make_shared<Matrix>("VC", nm, nm);
    Vint->compute(matrices_["VC"]);

    // => Nuclear Repulsions <= //

    auto Zs = std::make_shared<Matrix>("Zs", nA, 3);
    double** Zsp = Zs->pointer();

    auto Rinv = std::make_shared<Matrix>("Rinv", nA, nA);
    double** Rinvp = Rinv->pointer();

    for (int A = 0; A < nA; A++) {
        Zsp[A][0] = ZAp[A];
        Zsp[A][1] = ZBp[A];
        Zsp[A][2] = ZCp[A];
        dip_nucap[0] += ZAp[A]*(mol->x(A));
        dip_nucap[1] += ZAp[A]*(mol->y(A));
        dip_nucap[2] += ZAp[A]*(mol->z(A));
        dip_nucbp[0] += ZBp[A]*(mol->x(A));
        dip_nucbp[1] += ZBp[A]*(mol->y(A));
        dip_nucbp[2] += ZBp[A]*(mol->z(A));
        dip_nuccp[0] += ZCp[A]*(mol->x(A));
        dip_nuccp[1] += ZCp[A]*(mol->y(A));
        dip_nuccp[2] += ZCp[A]*(mol->z(A));
    }

    vectors_["DIP NUCA"] = dip_nuca;
    vectors_["DIP NUCB"] = dip_nucb;
    vectors_["DIP NUCC"] = dip_nucc;

    for (int A = 0; A < nA; A++) {
        for (int B = 0; B < nA; B++) {
            if (A == B) continue;
            Rinvp[A][B] = 1.0 / mol->xyz(A).distance(mol->xyz(B));
        }
    }

    /// Nuclear repulsion for A, B, C,
    Zs->print_out();
    Rinv->print_out();
    std::shared_ptr<Matrix> Enucs = linalg::triplet(Zs, Rinv, Zs, true, false, false);
    Enucs->scale(0.5);
    Enucs->set_name("E Nuc");
    Enucs->print_out();
    matrices_["E NUC"] = Enucs;

    double** Enucsp = Enucs->pointer();
    double Etot = 0.0;
    for (int A = 0; A < 3; A++) {
        for (int B = 0; B < 3; B++) {
            Etot += Enucsp[A][B];
        }
    }
    outfile->Printf("    Nuclear Repulsion Tot: %24.16E [Eh]\n", Etot);


    // are there any external potentials? If so, we need all QM atoms to feel their effects in the nuclear 
    // repulsion energy.  
    // Apparently, we were using C full strength, but the others I think get scaled by 0.5 because Rob counts
    // A->B and B->A separately and adds them (see a few lines up this fn...maybe due to how FSAPT files are written) 
    matrices_["Enucs"] = Enucs;
    Enucs->print();
    Etot += psi::sapt_nuclear_external_potential_matrix(
        reference_,
        matrices_,
        options_
    );
    // Enucs->print();
    Enucs = matrices_["Enucsp"];

    // => Print <= //

    // Zs->print();
    //Enucs->print();

    outfile->Printf("    Nuclear Repulsion Tot: %24.16E [Eh]\n", Etot);
    outfile->Printf("\n");
}


void FISAPT::coulomb() {
    outfile->Printf("  ==> Coulomb Integrals <==\n\n");

    // => Global JK Object <= //

    jk_ = JK::build_JK(primary_, reference_->get_basisset("DF_BASIS_SCF"), options_, false, doubles_);
    jk_->set_memory(doubles_);

    // => Build J and K for embedding <= //

    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    Cl.clear();
    Cr.clear();

    // => Prevent Failure if C, D, or E are empty <= //

    if (matrices_["LoccC"]->colspi()[0] > 0) {
        Cl.push_back(matrices_["LoccC"]);
        Cr.push_back(matrices_["LoccC"]);
    }

    jk_->set_do_J(true);
    jk_->set_do_K(true);
    jk_->initialize();
    jk_->print_header();

    jk_->compute();

    int nn = primary_->nbf();
    matrices_["JC"] = std::make_shared<Matrix>("JC", nn, nn);
    matrices_["KC"] = std::make_shared<Matrix>("KC", nn, nn);
    if (matrices_["LoccC"]->colspi()[0] > 0) {
        matrices_["JC"]->copy(J[0]);
        matrices_["KC"]->copy(K[0]);
    }
}

void FISAPT::scf() {
    outfile->Printf("  ==> Relaxed SCF Equations <==\n\n");

    // => Restricted Basis Sets with C Projected <= //

    std::vector<std::shared_ptr<Matrix> > Xs;
    Xs.push_back(matrices_["LoccA"]);
    Xs.push_back(matrices_["LoccB"]);
    Xs.push_back(matrices_["Cvir"]);
    matrices_["XC"] = linalg::horzcat(Xs);
    matrices_["XC"]->set_name("XC");

    // => Embedding Potential for C <= //

    std::shared_ptr<Matrix> WC(matrices_["VC"]->clone());
    WC->copy(matrices_["VC"]);
    WC->add(matrices_["JC"]);
    WC->add(matrices_["JC"]);
    WC->subtract(matrices_["KC"]);
    matrices_["WC"] = WC;

    // => A <= //

    outfile->Printf("  ==> SCF A: <==\n\n");
    std::shared_ptr<Matrix> VA_SCF(matrices_["VA"]->clone());
    VA_SCF->copy(matrices_["VA"]);
    if (reference_->has_potential_variable("C")) VA_SCF->add(matrices_["VE"]);
    std::shared_ptr<FISAPTSCF> scfA =
        std::make_shared<FISAPTSCF>(jk_, matrices_["E NUC"]->get(0, 0), matrices_["S"], matrices_["XC"], matrices_["T"],
                                    VA_SCF, matrices_["WC"], matrices_["LoccA"], options_);
    scfA->compute_energy();

    scalars_["E0 A"] = scfA->scalars()["E SCF"];
    matrices_["Cocc0A"] = scfA->matrices()["Cocc"];
    matrices_["Cvir0A"] = scfA->matrices()["Cvir"];
    matrices_["J0A"] = scfA->matrices()["J"];
    matrices_["K0A"] = scfA->matrices()["K"];
    vectors_["eps_occ0A"] = scfA->vectors()["eps_occ"];
    vectors_["eps_vir0A"] = scfA->vectors()["eps_vir"];

    // => B <= //

    outfile->Printf("  ==> SCF B: <==\n\n");
    std::shared_ptr<Matrix> VB_SCF(matrices_["VB"]->clone());
    VB_SCF->copy(matrices_["VB"]);
    if (reference_->has_potential_variable("C")) VB_SCF->add(matrices_["VE"]);
    std::shared_ptr<FISAPTSCF> scfB =
        std::make_shared<FISAPTSCF>(jk_, matrices_["E NUC"]->get(1, 1), matrices_["S"], matrices_["XC"], matrices_["T"],
                                    VB_SCF, matrices_["WC"], matrices_["LoccB"], options_);
    scfB->compute_energy();

    scalars_["E0 B"] = scfB->scalars()["E SCF"];
    matrices_["Cocc0B"] = scfB->matrices()["Cocc"];
    matrices_["Cvir0B"] = scfB->matrices()["Cvir"];
    matrices_["J0B"] = scfB->matrices()["J"];
    matrices_["K0B"] = scfB->matrices()["K"];
    vectors_["eps_occ0B"] = scfB->vectors()["eps_occ"];
    vectors_["eps_vir0B"] = scfB->vectors()["eps_vir"];
}

// Prep the computation to handle frozen core orbitals
void FISAPT::freeze_core() {
    outfile->Printf("  ==> Frozen Core <==\n\n");

    // => Frozen Core (for Disp) <= //

    std::shared_ptr<Molecule> mol = primary_->molecule();

    std::vector<int> none;
    std::vector<int> zero;
    zero.push_back(0);
    std::vector<int> one;
    one.push_back(1);
    std::shared_ptr<Molecule> molA = mol->extract_subsets(zero, none);
    std::shared_ptr<Molecule> molB = mol->extract_subsets(one, none);

    // std::shared_ptr<Molecule> molA = mol->extract_subsets({0},{});
    // std::shared_ptr<Molecule> molB = mol->extract_subsets({1},{});

    int nfocc0A = reference_->basisset()->n_frozen_core(options_.get_str("FREEZE_CORE"), molA);
    int nfocc0B = reference_->basisset()->n_frozen_core(options_.get_str("FREEZE_CORE"), molB);

    int nbf = matrices_["Cocc0A"]->rowspi()[0];
    int nocc0A = matrices_["Cocc0A"]->colspi()[0];
    int nocc0B = matrices_["Cocc0B"]->colspi()[0];
    int naocc0A = nocc0A - nfocc0A;
    int naocc0B = nocc0B - nfocc0B;
    int nvir0A = matrices_["Cvir0A"]->colspi()[0];
    int nvir0B = matrices_["Cvir0B"]->colspi()[0];
    int nmoA = nocc0A + nvir0A;
    int nmoB = nocc0B + nvir0B;

    outfile->Printf("\n");
    outfile->Printf("    ------------------\n");
    outfile->Printf("    %-6s %5s %5s\n", "Range", "A", "B");
    outfile->Printf("    ------------------\n");
    outfile->Printf("    %-6s %5d %5d\n", "nbf", nbf, nbf);
    outfile->Printf("    %-6s %5d %5d\n", "nmo", nmoA, nmoB);
    outfile->Printf("    %-6s %5d %5d\n", "nocc", nocc0A, nocc0B);
    outfile->Printf("    %-6s %5d %5d\n", "nvir", nvir0A, nvir0B);
    outfile->Printf("    %-6s %5d %5d\n", "nfocc", nfocc0A, nfocc0B);
    outfile->Printf("    %-6s %5d %5d\n", "naocc", naocc0A, naocc0B);
    outfile->Printf("    %-6s %5d %5d\n", "navir", nvir0A, nvir0B);
    outfile->Printf("    %-6s %5d %5d\n", "nfvir", 0, 0);
    outfile->Printf("    ------------------\n");
    outfile->Printf("\n");

    matrices_["Cfocc0A"] = std::make_shared<Matrix>("Cfocc0A", nbf, nfocc0A);
    matrices_["Caocc0A"] = std::make_shared<Matrix>("Caocc0A", nbf, naocc0A);
    matrices_["Cfocc0B"] = std::make_shared<Matrix>("Cfocc0B", nbf, nfocc0B);
    matrices_["Caocc0B"] = std::make_shared<Matrix>("Caocc0B", nbf, naocc0B);

    vectors_["eps_focc0A"] = std::make_shared<Vector>("eps_focc0A", nfocc0A);
    vectors_["eps_aocc0A"] = std::make_shared<Vector>("eps_aocc0A", naocc0A);
    vectors_["eps_focc0B"] = std::make_shared<Vector>("eps_focc0B", nfocc0B);
    vectors_["eps_aocc0B"] = std::make_shared<Vector>("eps_aocc0B", naocc0B);

    double** Cocc0Ap = matrices_["Cocc0A"]->pointer();
    double** Cocc0Bp = matrices_["Cocc0B"]->pointer();
    double** Cfocc0Ap = matrices_["Cfocc0A"]->pointer();
    double** Caocc0Ap = matrices_["Caocc0A"]->pointer();
    double** Cfocc0Bp = matrices_["Cfocc0B"]->pointer();
    double** Caocc0Bp = matrices_["Caocc0B"]->pointer();

    double* eps_occ0Ap = vectors_["eps_occ0A"]->pointer();
    double* eps_occ0Bp = vectors_["eps_occ0B"]->pointer();
    double* eps_focc0Ap = vectors_["eps_focc0A"]->pointer();
    double* eps_aocc0Ap = vectors_["eps_aocc0A"]->pointer();
    double* eps_focc0Bp = vectors_["eps_focc0B"]->pointer();
    double* eps_aocc0Bp = vectors_["eps_aocc0B"]->pointer();

    for (int m = 0; m < nbf; m++) {
        for (int a = 0; a < nfocc0A; a++) {
            Cfocc0Ap[m][a] = Cocc0Ap[m][a];
        }
        for (int a = 0; a < naocc0A; a++) {
            Caocc0Ap[m][a] = Cocc0Ap[m][a + nfocc0A];
        }
        for (int a = 0; a < nfocc0B; a++) {
            Cfocc0Bp[m][a] = Cocc0Bp[m][a];
        }
        for (int a = 0; a < naocc0B; a++) {
            Caocc0Bp[m][a] = Cocc0Bp[m][a + nfocc0B];
        }
    }

    for (int a = 0; a < nfocc0A; a++) {
        eps_focc0Ap[a] = eps_occ0Ap[a];
    }
    for (int a = 0; a < naocc0A; a++) {
        eps_aocc0Ap[a] = eps_occ0Ap[a + nfocc0A];
    }
    for (int a = 0; a < nfocc0B; a++) {
        eps_focc0Bp[a] = eps_occ0Bp[a];
    }
    for (int a = 0; a < naocc0B; a++) {
        eps_aocc0Bp[a] = eps_occ0Bp[a + nfocc0B];
    }

    // vectors_["eps_occ0A"]->print();
    // vectors_["eps_occ0B"]->print();
    // vectors_["eps_focc0A"]->print();
    // vectors_["eps_focc0B"]->print();
    // vectors_["eps_aocc0A"]->print();
    // vectors_["eps_aocc0B"]->print();
}


void FISAPT::unify() {
    outfile->Printf("  ==> Unification <==\n\n");

    // Matrices Cocc0A, Cocc0B come from the SCF calculation for noninteracting fragments A,B embedded in C
    // Matrices LoccA,LoccB,LoccC,LoccL come from the SCF calculation for entire molecule, and then localization.
    std::shared_ptr<Matrix> Cocc_A = matrices_["Cocc0A"];
    std::shared_ptr<Matrix> Cocc_B = matrices_["Cocc0B"];
    std::shared_ptr<Matrix> Cocc_C = matrices_["LoccC"];
    std::shared_ptr<Matrix> Cocc_L = matrices_["LoccL"];
    std::shared_ptr<Matrix> Cvir = matrices_["Cvir"];
    std::shared_ptr<Matrix> Sao = matrices_["S"];
    std::shared_ptr<Matrix> LinkA(Cocc_L->clone());
    std::shared_ptr<Matrix> LinkB(Cocc_L->clone());
    std::shared_ptr<Matrix> LinkC(Cocc_L->clone());
    std::shared_ptr<Matrix> DC_A(Sao->clone());  
    std::shared_ptr<Matrix> DC_B(Sao->clone());  
    std::shared_ptr<Matrix> DC_C(Sao->clone());  
    LinkA->copy(Cocc_L);
    LinkB->copy(Cocc_L);
    LinkC->copy(Cocc_L);
    double* FRAGp = vectors_["FRAG"]->pointer();
    double** LinkAp = LinkA->pointer();
    double** LinkBp = LinkB->pointer();
    double** LinkCp = LinkC->pointer();
    std::shared_ptr<Molecule> mol = primary_->molecule();
    int nA = mol->natom();
    std::shared_ptr<Vector> ZA = std::make_shared<Vector>("ZA", nA);
    std::shared_ptr<Vector> ZB = std::make_shared<Vector>("ZB", nA);
    std::shared_ptr<Vector> ZC = std::make_shared<Vector>("ZC", nA);

    std::string link_assignment = options_.get_str("FISAPT_LINK_ASSIGNMENT");
    std::string link_ortho = options_.get_str("FISAPT_LINK_ORTHO");

    // For the SAOn/SIAOn link assignments, the nuclear charge assignment to fragments is altered (one proton from the connecting atoms on A,B is no longer counted as part of C). Therefore, we need to revert to the unmodified nuclear charges saved earlier.
    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2" ) {
        ZA->copy(*vectors_["ZA_orig"]);
        ZB->copy(*vectors_["ZB_orig"]);
        ZC->copy(*vectors_["ZC_orig"]);
    }
    else {
        ZA->copy(*vectors_["ZA"]);
        ZB->copy(*vectors_["ZB"]);
        ZC->copy(*vectors_["ZC"]);
    }
    double* ZAp = ZA->pointer();
    double* ZBp = ZB->pointer();
    double* ZCp = ZC->pointer();
    // Now we will identify the link intrinsic hybrid orbitals that will be reassigned (with one electron each) from C to A,B.
    // These link IHOs will be kept in matrices thislinkA/thislinkB.

    std::vector<int> link_orbs;
    int nm = primary_->nbf();
    Dimension zero(1);
    Dimension dimnm(1);
    Dimension di(1);
    Dimension dj(1);
    dimnm[0] = nm;
    auto thislinkA = std::make_shared<Matrix>("thislinkA", nm, 1);
    auto thislinkB = std::make_shared<Matrix>("thislinkB", nm, 1);
    std::shared_ptr<Matrix> thislinkA0(thislinkA->clone());
    std::shared_ptr<Matrix> thislinkB0(thislinkB->clone());
    thislinkA->zero();
    thislinkB->zero();
    Slice rowslice(zero,dimnm);
    int na = matrices_["Locc"]->colspi()[0];
    int norbA = matrices_["Cocc0A"]->colspi()[0];
    int norbB = matrices_["Cocc0B"]->colspi()[0];
    int norbC = matrices_["LoccC"]->colspi()[0];
    int nlink = matrices_["LoccL"]->colspi()[0];
    int nvir = matrices_["Call"]->colspi()[0] - na;

    if (link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") {
//In the SIAO algorithms, we project the link IBOs onto IAOs of the respective fragment.
//IAOs have tails, but every IAO has a primary, dominant  center so we can easily assign them to fragments A/B/C.
        outfile->Printf("  ==> Link bond redistribution (SIAOn) <==\n\n");
        int natoms = primary_->molecule()->natom();  
        std::shared_ptr<Matrix> Iao(matrices_["IAO"]->clone());
        double** Iaop = Iao->pointer();
        std::shared_ptr<BasisSet> minao = reference_->get_basisset("MINAO");
        int nminao = minao->nbf();
        outfile->Printf("\n Number of MINAO functions: %d Number of IBOs: %d",nminao,na);
        auto Sminao = reference_->mintshelper()->ao_overlap(minao,minao);
//IAOs themselves are localized on atoms - find out which atoms and how well localized. The matrix atiao will have the weight of atom k in the IAO number l.
        auto atiao = std::make_shared<Matrix>("ATIAO", nminao, natoms);
        double** atiaop = atiao->pointer();
        atiao->zero();
        for (int k = 0; k < nminao; k++) {
            for (int l = 0; l < nm; l++) {
                int atind = primary_->function_to_center(l);
                atiaop[k][atind] += Iaop[l][k] * Iaop[l][k];
            }
        }
//      atiao->print();
//it seems that IAOs are ordered and localized in the same way as the MINAO basis functions, so 
//minao->function_to_center(k) returns the center for IAO number k.
//We can now express the two link IBOs as linear combinations of IAOs.
        std::shared_ptr<Matrix> CLiao = linalg::triplet(matrices_["IAO"], Sao, matrices_["LoccL"], true, false, false);
        double** CLiaop = CLiao->pointer();
        std::shared_ptr<Matrix> CLiaoA(CLiao->clone());
        std::shared_ptr<Matrix> CLiaoB(CLiao->clone());
        std::shared_ptr<Matrix> CLiaoC(CLiao->clone());
        double** CLiaoAp = CLiaoA->pointer();
        double** CLiaoBp = CLiaoB->pointer();
        double** CLiaoCp = CLiaoC->pointer();
//The FRAG matrix has 1/2/3 when the atom belongs to A/B/C, respectively.
        for (int k = 0; k < nminao; k++) {
          if (FRAGp[minao->function_to_center(k)] == 1.0) {
            for (int l = 0; l < 2; l++) {
              CLiaoBp[k][l] = 0.0 ;
              CLiaoCp[k][l] = 0.0 ;
            }
          } else
          if (FRAGp[minao->function_to_center(k)] == 2.0) {
            for (int l = 0; l < 2; l++) {
              CLiaoAp[k][l] = 0.0 ;
              CLiaoCp[k][l] = 0.0 ;
            }
          } else
          if (FRAGp[minao->function_to_center(k)] == 3.0) {
            for (int l = 0; l < 2; l++) {
              CLiaoAp[k][l] = 0.0 ;
              CLiaoBp[k][l] = 0.0 ;
            }
          } else
          {
          outfile->Printf(" Fragmentation problem! %3d %5.1f\n",k,FRAGp[minao->function_to_center(k)]);
          }
        }
// We no longer need ISAPT(C) nuclear charges and can overwrite them.
        vectors_["ZA"] = ZA;
        vectors_["ZB"] = ZB;
        vectors_["ZC"] = ZC;
        std::shared_ptr<Matrix> ortho_targetA = matrices_["LoccA"];
        std::shared_ptr<Matrix> ortho_targetB = matrices_["LoccB"];
        if (link_ortho == "FRAGMENT") {
            ortho_targetA = matrices_["Cocc0A"];
            ortho_targetB = matrices_["Cocc0B"];
            }
// At this point, the matrices CLiaoA,CLiaoB,CLiaoC have the 2 link orbitals expanded in IAOs of one fragment only, the rest are zeroed out.
// We now recast the projected vectors back to the AO basis, identify which one is A and which is B, 
// and orthogonalize to the occupied orbitals of a given fragment.
// Orthogonalization options:
//     FRAGMENT (default): orthogonalize to occupied orbitals of the noninteracting fragment A/B embedded in C
//     WHOLE: orthogonalize to occupied IBOs of the whole molecule localized on A/B
//     NONE: do not orthogonalize
        std::shared_ptr<Matrix> CLiaoA_AO = linalg::doublet(matrices_["IAO"], CLiaoA, false, false);
        std::shared_ptr<Matrix> CLiaoB_AO = linalg::doublet(matrices_["IAO"], CLiaoB, false, false);
        for (int l = 0; l < 2; l++) {
            if (typesL[l] == "AC") {
                outfile->Printf("\n Link orbital %3d belongs to A.\n",l);
                di[0] = l;
                dj[0] = l + 1;
                Slice colslice(di,dj);
                thislinkA = CLiaoA_AO->get_block(rowslice,colslice);
                thislinkA0->copy(thislinkA);
                if (link_ortho == "FRAGMENT" || link_ortho == "WHOLE") {
                    std::shared_ptr<Matrix> CoccLA_ortho = linalg::triplet(ortho_targetA, Sao, thislinkA, true, false, false);
                    // CoccLA_ortho->print();
                    double** CoccLA_orthop = CoccLA_ortho->pointer();
                    for (int j = 0; j < norbA; j++) {
                        di[0] = j;
                        dj[0] = j + 1;
                        Slice colslice(di,dj);
                        thislinkA->axpy(-CoccLA_orthop[j][0],ortho_targetA->get_block(rowslice,colslice));
                    }
                }
                std::shared_ptr<Matrix> thisnormA = linalg::triplet(thislinkA, Sao, thislinkA, true, false, false);
                double** thisnormAp = thisnormA->pointer();
                thislinkA->scale(1.0/std::sqrt(thisnormAp[0][0]));
                matrices_["thislinkAbak"] = thislinkA->clone();
            }
            else if (typesL[l] == "BC") {
                outfile->Printf("\n Link orbital %3d belongs to B.\n",l);
                di[0] = l;
                dj[0] = l + 1;
                Slice colslice(di,dj);
                thislinkB = CLiaoB_AO->get_block(rowslice,colslice);
                thislinkB0->copy(thislinkB);
                if (link_ortho == "FRAGMENT" || link_ortho == "WHOLE") {
                    std::shared_ptr<Matrix> CoccLB_ortho = linalg::triplet(ortho_targetB, Sao, thislinkB, true, false, false);
                    // CoccLB_ortho->print();
                    double** CoccLB_orthop = CoccLB_ortho->pointer();
                    for (int j = 0; j < norbB; j++) {
                        di[0] = j;
                        dj[0] = j + 1;
                        Slice colslice(di,dj);
                        thislinkB->axpy(-CoccLB_orthop[j][0],ortho_targetB->get_block(rowslice,colslice));
                    }
                }
                std::shared_ptr<Matrix> thisnormB = linalg::triplet(thislinkB, Sao, thislinkB, true, false, false);
                double** thisnormBp = thisnormB->pointer();
                thislinkB->scale(1.0/std::sqrt(thisnormBp[0][0]));
                matrices_["thislinkBbak"] = thislinkB->clone();
            }
            else {
                outfile->Printf("PROBLEM assigning link orbital %3d.",l);
            }
        } 
    }
    else if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" ) {
//In the SAO algorithms, we project the link IBOs onto AO basis functions of the respective fragment.
        outfile->Printf("  ==> Link bond redistribution (SAOn) <==\n\n");
// We no longer need ISAPT(C) nuclear charges and can overwrite them.
        vectors_["ZA"] = ZA;
        vectors_["ZB"] = ZB;
        vectors_["ZC"] = ZC;
        vectors_["ZA"]->print();
        vectors_["ZB"]->print();
        vectors_["ZC"]->print();
        for (int i = 0; i < nlink; i++) {
          link_orbs.push_back(0);
          for (int k = 0; k < nm; k++) {
            if (FRAGp[primary_->function_to_center(k)] == 1.0) {
              LinkBp[k][i] = 0.0 ;
              LinkCp[k][i] = 0.0 ;
            } else
            if (FRAGp[primary_->function_to_center(k)] == 2.0) {
              LinkAp[k][i] = 0.0 ;
              LinkCp[k][i] = 0.0 ;
            } else
            if (FRAGp[primary_->function_to_center(k)] == 3.0) {
              LinkAp[k][i] = 0.0 ;
              LinkBp[k][i] = 0.0 ;
            } else
            {
            outfile->Printf(" Fragmentation problem! %3d %5.1f\n",k,FRAGp[primary_->function_to_center(k)]);
            }
          }
        }
        std::shared_ptr<Matrix> LinkAnorm = linalg::triplet(LinkA, Sao, LinkA, true, false, false);
        outfile->Printf("\n Norm of link vectors truncated to A: \n");
        LinkAnorm->print();
        std::shared_ptr<Matrix> LinkBnorm = linalg::triplet(LinkB, Sao, LinkB, true, false, false);
        outfile->Printf("\n Norm of link vectors truncated to B: \n");
        LinkBnorm->print();
        std::shared_ptr<Matrix> LinkCnorm = linalg::triplet(LinkC, Sao, LinkC, true, false, false);
        outfile->Printf("\n Norm of link vectors truncated to C: \n");
        LinkCnorm->print();
        std::shared_ptr<Matrix> LinkAllnorm = linalg::triplet(Cocc_L, Sao, Cocc_L, true, false, false);
        outfile->Printf("\n Norm of link vectors untruncated: \n");
        LinkAllnorm->print();

// At this point, the matrices LinkA,LinkB,LinkC have the 2 link orbitals expanded in AOs of one fragment only, the rest are zeroed out.
// We now identify which one is A and which is B, and orthogonalize to the occupied orbitals of a given fragment.
// Orthogonalization options:
//     FRAGMENT (default): orthogonalize to occupied orbitals of the noninteracting fragment A/B embedded in C
//     WHOLE: orthogonalize to occupied IBOs of the whole molecule localized on A/B
//     NONE: do not orthogonalize
        std::shared_ptr<Matrix> ortho_targetA = matrices_["LoccA"];
        std::shared_ptr<Matrix> ortho_targetB = matrices_["LoccB"];
        if (link_ortho == "FRAGMENT") {
            ortho_targetA = matrices_["Cocc0A"];
            ortho_targetB = matrices_["Cocc0B"];
            }
        for (int i = 0; i < nlink; i++) {
            if (typesL[i] == "AC") {
                outfile->Printf("\n Link orbital %3d belongs to A.\n",i);
                di[0] = i;
                dj[0] = i + 1;
                Slice colslice(di,dj);
                thislinkA = LinkA->get_block(rowslice,colslice);
                thislinkA0->copy(thislinkA);
                if (link_ortho == "FRAGMENT" || link_ortho == "WHOLE") {
                    std::shared_ptr<Matrix> CoccLA_ortho = linalg::triplet(ortho_targetA, Sao, thislinkA, true, false, false);
                    // CoccLA_ortho->print();
                    double** CoccLA_orthop = CoccLA_ortho->pointer();
                    for (int j = 0; j < norbA; j++) {
                        di[0] = j;
                        dj[0] = j + 1;
                        Slice colslice(di,dj);
                        thislinkA->axpy(-CoccLA_orthop[j][0],ortho_targetA->get_block(rowslice,colslice));
                    }
                }
                std::shared_ptr<Matrix> thisnormA = linalg::triplet(thislinkA, Sao, thislinkA, true, false, false);
                double** thisnormAp = thisnormA->pointer();
                thislinkA->scale(1.0/std::sqrt(thisnormAp[0][0]));
                matrices_["thislinkAbak"] = thislinkA->clone();
            }
            else if (typesL[i] == "BC") {
                outfile->Printf("\n Link orbital %3d belongs to B.\n",i);
                di[0] = i;
                dj[0] = i + 1;
                Slice colslice(di,dj);
                Slice rowslice(zero,dimnm);
                thislinkB = LinkB->get_block(rowslice,colslice);
                thislinkB0->copy(thislinkB);
                if (link_ortho == "FRAGMENT" || link_ortho == "WHOLE") {
                    std::shared_ptr<Matrix> CoccLB_ortho = linalg::triplet(ortho_targetB, Sao, thislinkB, true, false, false);
                    // CoccLB_ortho->print();
                    double** CoccLB_orthop = CoccLB_ortho->pointer();
                    for (int j = 0; j < norbB; j++) {
                        di[0] = j;
                        dj[0] = j + 1;
                        Slice colslice(di,dj);
                        thislinkB->axpy(-CoccLB_orthop[j][0],ortho_targetB->get_block(rowslice,colslice));
                    }
                }
                std::shared_ptr<Matrix> thisnormB = linalg::triplet(thislinkB, Sao, thislinkB, true, false, false);
                double** thisnormBp = thisnormB->pointer();
                thislinkB->scale(1.0/std::sqrt(thisnormBp[0][0]));
                matrices_["thislinkBbak"] = thislinkB->clone();
            }
            else {
                outfile->Printf("PROBLEM assigning link orbital %3d.",i);
            }
        }
    }

    std::shared_ptr<Matrix> D_A = linalg::doublet(Cocc_A, Cocc_A, false, true);
    std::shared_ptr<Matrix> D_B = linalg::doublet(Cocc_B, Cocc_B, false, true);
    std::shared_ptr<Matrix> D_C(D_A->clone());
    D_C->zero();
    if (Cocc_C->colspi()[0] > 0) {
        D_C = linalg::doublet(Cocc_C, Cocc_C, false, true);
    }
    matrices_["D_A0"] = D_A->clone();
    matrices_["D_B0"] = D_B->clone();
    matrices_["D_C0"] = D_C->clone();

// At this point, for all SAOn/SIAOn variants, matrices thislinkA and thislinkB contain the link hybrid orbitals that are being reassigned from C to A/B.
// We can now introduce this reassignment to the fragment density matrices D_A, D_B, D_C.
// Note the multiplication by 1/sqrt{2} which correspond to multiplying the density matrix contribution by 1/2 (the orbital is singly occupied).
    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2" ) {
        outfile->Printf("\n Trace of D_A before link reassignment: %10.6f \n",D_A->vector_dot(Sao));
        outfile->Printf("\n Trace of D_B before link reassignment: %10.6f \n",D_B->vector_dot(Sao));
        outfile->Printf("\n Trace of D_C before link reassignment: %10.6f \n",D_C->vector_dot(Sao));
        thislinkA->scale(0.7071067811865475244);
        thislinkB->scale(0.7071067811865475244);
        D_A->axpy(1.0,linalg::doublet(thislinkA, thislinkA, false, true)); 
        D_B->axpy(1.0,linalg::doublet(thislinkB, thislinkB, false, true)); 
        if (Cocc_C->colspi()[0] > 0) {
            D_C->axpy(-1.0,linalg::doublet(thislinkA, thislinkA, false, true)); 
            D_C->axpy(-1.0,linalg::doublet(thislinkB, thislinkB, false, true)); 
        }
        outfile->Printf("\n Trace of D_A after link reassignment: %10.6f \n",D_A->vector_dot(Sao));
        outfile->Printf("\n Trace of D_B after link reassignment: %10.6f \n",D_B->vector_dot(Sao));
        outfile->Printf("\n Trace of D_C after link reassignment: %10.6f \n",D_C->vector_dot(Sao));
    }

    matrices_["D_A"] = D_A;
    matrices_["D_B"] = D_B;
    matrices_["D_C"] = D_C;

    // Incorrect for this application: C is not frozen in these orbitals
    // std::shared_ptr<Matrix> P_A = linalg::doublet(matrices_["Cvir0A"], matrices_["Cvir0A"], false, true);
    // std::shared_ptr<Matrix> P_B = linalg::doublet(matrices_["Cvir0B"], matrices_["Cvir0B"], false, true);

    // PA and PB are used only to define the complement of DA and DB in the DCBS
    std::shared_ptr<Matrix> P_A =
        linalg::doublet(reference_->Ca_subset("AO", "ALL"), reference_->Ca_subset("AO", "ALL"), false, true);
    P_A->subtract(D_A);

    std::shared_ptr<Matrix> P_B =
        linalg::doublet(reference_->Ca_subset("AO", "ALL"), reference_->Ca_subset("AO", "ALL"), false, true);
    P_B->subtract(D_B);

    matrices_["P_A"] = P_A;
    matrices_["P_B"] = P_B;

    matrices_["Cocc_A"] = matrices_["Cocc0A"];
    matrices_["Cocc_B"] = matrices_["Cocc0B"];
    matrices_["Cocc_C"] = matrices_["Cocc0C"];

    matrices_["V_A"] = matrices_["VA"];
    matrices_["V_B"] = matrices_["VB"];
    matrices_["V_C"] = matrices_["VC"];

    std::shared_ptr<Matrix> J_A(matrices_["J0A"]->clone());
    std::shared_ptr<Matrix> K_A(matrices_["K0A"]->clone());
    std::shared_ptr<Matrix> J_B(matrices_["J0B"]->clone());
    std::shared_ptr<Matrix> K_B(matrices_["K0B"]->clone());
    std::shared_ptr<Matrix> J_C(matrices_["JC"]->clone());
    std::shared_ptr<Matrix> K_C(matrices_["KC"]->clone());

// If density has been reassigned between fragments (the SAOn/SIAOn variants), we need to redo the Coulomb and exchange operators now.
    SharedMatrix AlloccA(Cocc_A->clone());
    SharedMatrix AlloccB(Cocc_B->clone());
    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2" ) {
        AlloccA = linalg::horzcat({Cocc_A,thislinkA});
        AlloccB = linalg::horzcat({Cocc_B,thislinkB});

        std::vector<SharedMatrix>& Cl = jk_->C_left();
        std::vector<SharedMatrix>& Cr = jk_->C_right();
        const std::vector<SharedMatrix>& J = jk_->J();
        const std::vector<SharedMatrix>& K = jk_->K();

        Cl.clear();
        Cr.clear();
        Cl.push_back(thislinkA);
        Cr.push_back(thislinkA);
        Cl.push_back(thislinkB);
        Cr.push_back(thislinkB);

        jk_->compute();

        J_A->add(J[0]);
        K_A->add(K[0]);
        J_B->add(J[1]);
        K_B->add(K[1]);
        J_C->subtract(J[0]);
        K_C->subtract(K[0]);
        J_C->subtract(J[1]);
        K_C->subtract(K[1]);
        matrices_["JLA"] = J[0]->clone(); 
        matrices_["KLA"] = K[0]->clone(); 
        matrices_["JLB"] = J[1]->clone(); 
        matrices_["KLB"] = K[1]->clone(); 
    }
    else {
// dummy zero matrices for other link assignments
        matrices_["JLA"] = matrices_["VC"]->clone();
        matrices_["KLA"] = matrices_["VC"]->clone();
        matrices_["JLB"] = matrices_["VC"]->clone();
        matrices_["KLB"] = matrices_["VC"]->clone();
        matrices_["JLA"]->zero();
        matrices_["JLB"]->zero();
        matrices_["KLA"]->zero();
        matrices_["KLB"]->zero();
    }

    matrices_["AlloccA"] = AlloccA;
    matrices_["AlloccB"] = AlloccB;
    matrices_["thislinkA"] = thislinkA;
    matrices_["thislinkB"] = thislinkB;
    matrices_["thislinkA0"] = thislinkA0;
    matrices_["thislinkB0"] = thislinkB0;
    matrices_["E NUC0"] = matrices_["E NUC"];
    matrices_["J_A"] = J_A;
    matrices_["J_B"] = J_B;
    matrices_["K_A"] = K_A;
    matrices_["K_B"] = K_B;
    matrices_["J_C"] = J_C;
    matrices_["K_C"] = K_C;
    vectors_["DIP NUCA0"] = vectors_["DIP NUCA"];
    vectors_["DIP NUCB0"] = vectors_["DIP NUCB"];
    vectors_["DIP NUCC0"] = vectors_["DIP NUCC"];
}

void FISAPT::unify_part2() {
// What must happen between unify() and unify_part2() is a recalculation of one-electron integrals with updated nuclear charges,
// that is, another call to nuclear().
    outfile->Printf("  ==> Unification (part 2) <==\n\n");

    std::string link_assignment = options_.get_str("FISAPT_LINK_ASSIGNMENT");
    std::string link_ortho = options_.get_str("FISAPT_LINK_ORTHO");
    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2" ) {

//Now redo monomer SCF calculations with the hybrid orbital on B removed from the embedding potential for A
//and vice versa.
//This is the FIRST iteration towards self-consistency.
//If zero iterations requested (SAO0/SIAO0), we still come in here but do not save the SCF vectors. MAybe that is unnecessary? 
//Anyway, SAO0/SIAO0 is bad, you should always do at least one iteration.
//old nuclear potential, with +1 charges on the connecting atoms of A/B belonging to C
    matrices_["V_A0"] = matrices_["V_A"];
    matrices_["V_B0"] = matrices_["V_B"];
    matrices_["V_C0"] = matrices_["V_C"];
//new, updated nuclear potential, with full charges on A/B belonging to A/B
    matrices_["V_A"] = matrices_["VA"];
    matrices_["V_B"] = matrices_["VB"];
    matrices_["V_C"] = matrices_["VC"];

    outfile->Printf("  ==> Redo Relaxed SCF Equations <==\n\n");

    // => Restricted Basis Sets with C Projected <= //

    std::vector<std::shared_ptr<Matrix> > Xs;
    Xs.push_back(matrices_["LoccA"]);
    Xs.push_back(matrices_["LoccB"]);
    Xs.push_back(matrices_["Cvir"]);
    matrices_["XC"] = linalg::horzcat(Xs);
    matrices_["XC"]->set_name("XC");

    // => Embedding Potential for C <= //

    std::shared_ptr<Matrix> WC(matrices_["V_C"]->clone());
    WC->copy(matrices_["V_C"]);
    WC->add(matrices_["J_C"]);
    WC->add(matrices_["J_C"]);
    WC->subtract(matrices_["K_C"]);
    std::shared_ptr<Matrix> WCA(WC->clone());
    WCA->copy(WC);
    WCA->add(matrices_["JLA"]) ;
    WCA->add(matrices_["JLA"]) ;
    WCA->subtract(matrices_["KLA"]) ;
    std::shared_ptr<Matrix> WCB(WC->clone());
    WCB->copy(WC);
    WCB->add(matrices_["JLB"]) ;
    WCB->add(matrices_["JLB"]) ;
    WCB->subtract(matrices_["KLB"]) ;
    matrices_["WCA"] = WCA;
    matrices_["WCB"] = WCB;

    // => A <= //

    outfile->Printf("  ==> SCF A: <==\n\n");
    outfile->Printf("  Nuclear repulsion: %10.6f\n",matrices_["E NUC"]->get(0, 0));
    std::shared_ptr<Matrix> VA_SCF(matrices_["V_A"]->clone());
    VA_SCF->copy(matrices_["V_A"]);
    if (reference_->has_potential_variable("C")) VA_SCF->add(matrices_["VE"]);
    double nucrep = matrices_["E NUC"]->get(0, 0);
    nucrep += 2.0*linalg::triplet(matrices_["thislinkA"], matrices_["T"], matrices_["thislinkA"], true, false, false)->get(0,0);
    nucrep += 2.0*linalg::triplet(matrices_["thislinkA"], VA_SCF, matrices_["thislinkA"], true, false, false)->get(0,0);
    nucrep += 2.0*linalg::triplet(matrices_["thislinkA"], matrices_["V_C"], matrices_["thislinkA"], true, false, false)->get(0,0);
    nucrep += 4.0*matrices_["D_C"]->vector_dot(matrices_["JLA"]) - 2.0*matrices_["D_C"]->vector_dot(matrices_["KLA"] );
    // Seems ultra-weird but apparently these "self-interaction" terms need to be included.
    nucrep += 2.0*linalg::triplet(matrices_["thislinkA"], matrices_["JLA"], matrices_["thislinkA"], true, false, false)->get(0,0);
    nucrep -= linalg::triplet(matrices_["thislinkA"], matrices_["KLA"], matrices_["thislinkA"], true, false, false)->get(0,0);
    std::shared_ptr<FISAPTSCF> scfA =
        std::make_shared<FISAPTSCF>(jk_, nucrep, matrices_["S"], matrices_["XC"], matrices_["T"],
                                    VA_SCF, matrices_["WCA"], matrices_["LoccA"], options_);
    scfA->compute_energy();

    scalars_["E0 A"] = scfA->scalars()["E SCF"];
    if (link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO1" || link_assignment == "SIAO2" ) {
        matrices_["Cocc0A"] = scfA->matrices()["Cocc"];
        matrices_["Cvir0A"] = scfA->matrices()["Cvir"];
        matrices_["J0A"] = scfA->matrices()["J"];
        matrices_["K0A"] = scfA->matrices()["K"];
        vectors_["eps_occ0A"] = scfA->vectors()["eps_occ"];
        vectors_["eps_vir0A"] = scfA->vectors()["eps_vir"];
    }

    // => B <= //

    outfile->Printf("  ==> SCF B: <==\n\n");
    outfile->Printf("  Nuclear repulsion: %10.6f\n",matrices_["E NUC"]->get(1, 1));
    std::shared_ptr<Matrix> VB_SCF(matrices_["V_B"]->clone());
    VB_SCF->copy(matrices_["V_B"]);
    if (reference_->has_potential_variable("C")) VB_SCF->add(matrices_["VE"]);
    nucrep = matrices_["E NUC"]->get(1, 1);
    nucrep += 2.0*linalg::triplet(matrices_["thislinkB"], matrices_["T"], matrices_["thislinkB"], true, false, false)->get(0,0);
    nucrep += 2.0*linalg::triplet(matrices_["thislinkB"], VB_SCF, matrices_["thislinkB"], true, false, false)->get(0,0);
    nucrep += 2.0*linalg::triplet(matrices_["thislinkB"], matrices_["V_C"], matrices_["thislinkB"], true, false, false)->get(0,0);
    nucrep += 4.0*matrices_["D_C"]->vector_dot(matrices_["JLB"]) - 2.0*matrices_["D_C"]->vector_dot(matrices_["KLB"] );
    // Seems ultra-weird but apparently these "self-interaction" terms need to be included.
    nucrep += 2.0*linalg::triplet(matrices_["thislinkB"], matrices_["JLB"], matrices_["thislinkB"], true, false, false)->get(0,0);
    nucrep -= linalg::triplet(matrices_["thislinkB"], matrices_["KLB"], matrices_["thislinkB"], true, false, false)->get(0,0);
    std::shared_ptr<FISAPTSCF> scfB =
        std::make_shared<FISAPTSCF>(jk_, nucrep, matrices_["S"], matrices_["XC"], matrices_["T"],
                                    VB_SCF, matrices_["WCB"], matrices_["LoccB"], options_);
    scfB->compute_energy();

    scalars_["E0 B"] = scfB->scalars()["E SCF"];
    if (link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO1" || link_assignment == "SIAO2" ) {
        matrices_["Cocc0B"] = scfB->matrices()["Cocc"];
        matrices_["Cvir0B"] = scfB->matrices()["Cvir"];
        matrices_["J0B"] = scfB->matrices()["J"];
        matrices_["K0B"] = scfB->matrices()["K"];
        vectors_["eps_occ0B"] = scfB->vectors()["eps_occ"];
        vectors_["eps_vir0B"] = scfB->vectors()["eps_vir"];
    }

//Now reorthogonalize the truncated link orbital to the new A/B SCF space, completing the FIRST iteration to self-consistency.
//thislinkA0/thislinkB0 are the original unorthogonalized link IHOs computed before.
    if (link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO1" || link_assignment == "SIAO2") {
        std::shared_ptr<Matrix> thislinkA(matrices_["thislinkA0"]->clone());
        std::shared_ptr<Matrix> thislinkB(matrices_["thislinkB0"]->clone());
        int norbA = matrices_["Cocc0A"]->colspi()[0];
        int norbB = matrices_["Cocc0B"]->colspi()[0];
        int nm = primary_->nbf();
        Dimension zero(1);
        Dimension dimnm(1);
        Dimension di(1);
        Dimension dj(1);
        dimnm[0] = nm;
        Slice rowslice(zero,dimnm);
        std::shared_ptr<Matrix> ortho_targetA = matrices_["LoccA"];
        std::shared_ptr<Matrix> ortho_targetB = matrices_["LoccB"];
        if (link_ortho == "FRAGMENT") {
            ortho_targetA = matrices_["Cocc0A"];
            ortho_targetB = matrices_["Cocc0B"];
            }
        if (link_ortho == "FRAGMENT" || link_ortho == "WHOLE") {
            std::shared_ptr<Matrix> CoccLA_ortho = linalg::triplet(ortho_targetA, matrices_["S"], thislinkA, true, false, false);
            // CoccLA_ortho->print();
            double** CoccLA_orthop = CoccLA_ortho->pointer();
            for (int j = 0; j < norbA; j++) {
                di[0] = j;
                dj[0] = j + 1;
                Slice colslice(di,dj);
                thislinkA->axpy(-CoccLA_orthop[j][0],ortho_targetA->get_block(rowslice,colslice));
            }
        }
        std::shared_ptr<Matrix> thisnormA = linalg::triplet(thislinkA, matrices_["S"], thislinkA, true, false, false);
        double** thisnormAp = thisnormA->pointer();
        thislinkA->scale(1.0/std::sqrt(thisnormAp[0][0]));
        if (link_ortho == "FRAGMENT" || link_ortho == "WHOLE") {
            std::shared_ptr<Matrix> CoccLB_ortho = linalg::triplet(ortho_targetB, matrices_["S"], thislinkB, true, false, false);
            // CoccLB_ortho->print();
            double** CoccLB_orthop = CoccLB_ortho->pointer();
            for (int j = 0; j < norbB; j++) {
                di[0] = j;
                dj[0] = j + 1;
                Slice colslice(di,dj);
                thislinkB->axpy(-CoccLB_orthop[j][0],ortho_targetB->get_block(rowslice,colslice));
            }
        }
        std::shared_ptr<Matrix> thisnormB = linalg::triplet(thislinkB, matrices_["S"], thislinkB, true, false, false);
        double** thisnormBp = thisnormB->pointer();
        thislinkB->scale(1.0/std::sqrt(thisnormBp[0][0]));

// We can now use the link IHOs from the first iteration to update the fragment density matrices D_A, D_B, D_C.
// Note the multiplication by 1/sqrt{2} which correspond to multiplying the density matrix contribution by 1/2 (the orbital is singly occupied).
        std::shared_ptr<Matrix> Cocc_C = matrices_["LoccC"];
        std::shared_ptr<Matrix> D_A = linalg::doublet(matrices_["Cocc0A"], matrices_["Cocc0A"], false, true);
        std::shared_ptr<Matrix> D_B = linalg::doublet(matrices_["Cocc0B"], matrices_["Cocc0B"], false, true);
        std::shared_ptr<Matrix> D_C(D_A->clone());
        D_C->zero();
        if (Cocc_C->colspi()[0] > 0) {
            D_C = linalg::doublet(Cocc_C, Cocc_C, false, true);
        }
 
        thislinkA->scale(0.7071067811865475244);
        thislinkB->scale(0.7071067811865475244);
        D_A->axpy(1.0,linalg::doublet(thislinkA, thislinkA, false, true)); 
        D_B->axpy(1.0,linalg::doublet(thislinkB, thislinkB, false, true)); 
        if (Cocc_C->colspi()[0] > 0) {
            D_C->axpy(-1.0,linalg::doublet(thislinkA, thislinkA, false, true)); 
            D_C->axpy(-1.0,linalg::doublet(thislinkB, thislinkB, false, true)); 
        }
 
        matrices_["D_A"] = D_A;
        matrices_["D_B"] = D_B;
        matrices_["D_C"] = D_C;
 
        // PA and PB are used only to define the complement of DA and DB in the DCBS
        std::shared_ptr<Matrix> P_A =
            linalg::doublet(reference_->Ca_subset("AO", "ALL"), reference_->Ca_subset("AO", "ALL"), false, true);
        P_A->subtract(D_A);
 
        std::shared_ptr<Matrix> P_B =
            linalg::doublet(reference_->Ca_subset("AO", "ALL"), reference_->Ca_subset("AO", "ALL"), false, true);
        P_B->subtract(D_B);
 
        matrices_["P_A"] = P_A;
        matrices_["P_B"] = P_B;
 
        matrices_["Cocc_A"] = matrices_["Cocc0A"];
        matrices_["Cocc_B"] = matrices_["Cocc0B"];

        std::shared_ptr<Matrix> J_A(matrices_["J0A"]->clone());
        std::shared_ptr<Matrix> K_A(matrices_["K0A"]->clone());
        std::shared_ptr<Matrix> J_B(matrices_["J0B"]->clone());
        std::shared_ptr<Matrix> K_B(matrices_["K0B"]->clone());
        std::shared_ptr<Matrix> J_C(matrices_["JC"]->clone());
        std::shared_ptr<Matrix> K_C(matrices_["KC"]->clone());

// Now redo the Coulomb and exchange operators with the new densities. We do just the link part and save it separately as well.
        SharedMatrix AlloccA(matrices_["Cocc0A"]->clone());
        SharedMatrix AlloccB(matrices_["Cocc0B"]->clone());
        AlloccA = linalg::horzcat({matrices_["Cocc0A"],thislinkA});
        AlloccB = linalg::horzcat({matrices_["Cocc0B"],thislinkB});
 
        std::vector<SharedMatrix>& Cl = jk_->C_left();
        std::vector<SharedMatrix>& Cr = jk_->C_right();
 
        const std::vector<SharedMatrix>& J = jk_->J();
        const std::vector<SharedMatrix>& K = jk_->K();
 
        Cl.clear();
        Cr.clear();
 
        Cl.push_back(thislinkA);
        Cr.push_back(thislinkA);
        Cl.push_back(thislinkB);
        Cr.push_back(thislinkB);
 
        jk_->compute();
 
        J_A->add(J[0]);
        K_A->add(K[0]);
        J_B->add(J[1]);
        K_B->add(K[1]);
        J_C->subtract(J[0]);
        K_C->subtract(K[0]);
        J_C->subtract(J[1]);
        K_C->subtract(K[1]);
        matrices_["JLA"] = J[0]->clone(); 
        matrices_["KLA"] = K[0]->clone(); 
        matrices_["JLB"] = J[1]->clone(); 
        matrices_["KLB"] = K[1]->clone();
 
        matrices_["AlloccA"] = AlloccA;
        matrices_["AlloccB"] = AlloccB;
        matrices_["thislinkA"] = thislinkA;
        matrices_["thislinkB"] = thislinkB;
        matrices_["J_A"] = J_A;
        matrices_["J_B"] = J_B;
        matrices_["K_A"] = K_A;
        matrices_["K_B"] = K_B;
        matrices_["J_C"] = J_C;
        matrices_["K_C"] = K_C;
     
    if (link_assignment == "SAO2" || link_assignment == "SIAO2" ) {

//Now redo monomer SCF calculations with the hybrid orbital on B removed from the embedding potential for A
//and vice versa.
//This is the SECOND iteration towards self-consistency, performed only when SAO2/SIAO2 is requested.
    outfile->Printf("  ==> Redo Relaxed SCF Equations (SAO2/SIAO2) <==\n\n");

    // => Embedding Potential for C <= //

    WC->copy(matrices_["VC"]);
    WC->add(matrices_["J_C"]);
    WC->add(matrices_["J_C"]);
    WC->subtract(matrices_["K_C"]);
    WCA->copy(WC);
    WCA->add(matrices_["JLA"]) ;
    WCA->add(matrices_["JLA"]) ;
    WCA->subtract(matrices_["KLA"]) ;
    WCB->copy(WC);
    WCB->add(matrices_["JLB"]) ;
    WCB->add(matrices_["JLB"]) ;
    WCB->subtract(matrices_["KLB"]) ;
    matrices_["WCA"] = WCA;
    matrices_["WCB"] = WCB;

    // => A <= //

    outfile->Printf("  ==> SCF A: <==\n\n");
    outfile->Printf("  Nuclear repulsion: %10.6f\n",matrices_["E NUC"]->get(0, 0));
    nucrep = matrices_["E NUC"]->get(0, 0);
    nucrep += 2.0*linalg::triplet(matrices_["thislinkA"], matrices_["T"], matrices_["thislinkA"], true, false, false)->get(0,0);
    nucrep += 2.0*linalg::triplet(matrices_["thislinkA"], VA_SCF, matrices_["thislinkA"], true, false, false)->get(0,0);
    nucrep += 2.0*linalg::triplet(matrices_["thislinkA"], matrices_["VC"], matrices_["thislinkA"], true, false, false)->get(0,0);
    nucrep += 4.0*matrices_["D_C"]->vector_dot(matrices_["JLA"]) - 2.0*matrices_["D_C"]->vector_dot(matrices_["KLA"] );
    nucrep += 2.0*linalg::triplet(matrices_["thislinkA"], matrices_["JLA"], matrices_["thislinkA"], true, false, false)->get(0,0);
    nucrep -= linalg::triplet(matrices_["thislinkA"], matrices_["KLA"], matrices_["thislinkA"], true, false, false)->get(0,0);
    std::shared_ptr<FISAPTSCF> scfA2 =
        std::make_shared<FISAPTSCF>(jk_, nucrep, matrices_["S"], matrices_["XC"], matrices_["T"],
                                    VA_SCF, matrices_["WCA"], matrices_["LoccA"], options_);
    scfA2->compute_energy();

    scalars_["E0 A"] = scfA2->scalars()["E SCF"];
    matrices_["Cocc0A"] = scfA2->matrices()["Cocc"];
    matrices_["Cvir0A"] = scfA2->matrices()["Cvir"];
    matrices_["J0A"] = scfA2->matrices()["J"];
    matrices_["K0A"] = scfA2->matrices()["K"];
    vectors_["eps_occ0A"] = scfA2->vectors()["eps_occ"];
    vectors_["eps_vir0A"] = scfA2->vectors()["eps_vir"];

    // => B <= //

    outfile->Printf("  ==> SCF B: <==\n\n");
    outfile->Printf("  Nuclear repulsion: %10.6f\n",matrices_["E NUC"]->get(1, 1));
    nucrep = matrices_["E NUC"]->get(1, 1);
    nucrep += 2.0*linalg::triplet(matrices_["thislinkB"], matrices_["T"], matrices_["thislinkB"], true, false, false)->get(0,0);
    nucrep += 2.0*linalg::triplet(matrices_["thislinkB"], VB_SCF, matrices_["thislinkB"], true, false, false)->get(0,0);
    nucrep += 2.0*linalg::triplet(matrices_["thislinkB"], matrices_["VC"], matrices_["thislinkB"], true, false, false)->get(0,0);
    nucrep += 4.0*matrices_["D_C"]->vector_dot(matrices_["JLB"]) - 2.0*matrices_["D_C"]->vector_dot(matrices_["KLB"] );
    nucrep += 2.0*linalg::triplet(matrices_["thislinkB"], matrices_["JLB"], matrices_["thislinkB"], true, false, false)->get(0,0);
    nucrep -= linalg::triplet(matrices_["thislinkB"], matrices_["KLB"], matrices_["thislinkB"], true, false, false)->get(0,0);
    std::shared_ptr<FISAPTSCF> scfB2 =
        std::make_shared<FISAPTSCF>(jk_, nucrep, matrices_["S"], matrices_["XC"], matrices_["T"],
                                    VB_SCF, matrices_["WCB"], matrices_["LoccB"], options_);
    scfB2->compute_energy();

    scalars_["E0 B"] = scfB2->scalars()["E SCF"];
    matrices_["Cocc0B"] = scfB2->matrices()["Cocc"];
    matrices_["Cvir0B"] = scfB2->matrices()["Cvir"];
    matrices_["J0B"] = scfB2->matrices()["J"];
    matrices_["K0B"] = scfB2->matrices()["K"];
    vectors_["eps_occ0B"] = scfB2->vectors()["eps_occ"];
    vectors_["eps_vir0B"] = scfB2->vectors()["eps_vir"];

//Now reorthogonalize the truncated link orbital to the new A/B SCF space, completing the SECOND iteration to self-consistency.
//thislinkA0/thislinkB0 are the original unorthogonalized link IHOs computed before.
    thislinkA->copy(matrices_["thislinkA0"]);
    thislinkB->copy(matrices_["thislinkB0"]);
    if (link_ortho == "FRAGMENT") {
        ortho_targetA = matrices_["Cocc0A"];
        ortho_targetB = matrices_["Cocc0B"];
        }
    else {
        ortho_targetA = matrices_["LoccA"];
        ortho_targetB = matrices_["LoccB"];
        } 
    if (link_ortho == "FRAGMENT" || link_ortho == "WHOLE") {
        std::shared_ptr<Matrix> CoccLA_ortho = linalg::triplet(ortho_targetA, matrices_["S"], thislinkA, true, false, false);
        // CoccLA_ortho->print();
        double** CoccLA_orthop = CoccLA_ortho->pointer();
        for (int j = 0; j < norbA; j++) {
            di[0] = j;
            dj[0] = j + 1;
            Slice colslice(di,dj);
            thislinkA->axpy(-CoccLA_orthop[j][0],ortho_targetA->get_block(rowslice,colslice));
        }
    }
    thisnormA = linalg::triplet(thislinkA, matrices_["S"], thislinkA, true, false, false);
    double** thisnormAp = thisnormA->pointer();
    thislinkA->scale(1.0/std::sqrt(thisnormAp[0][0]));

    if (link_ortho == "FRAGMENT" || link_ortho == "WHOLE") {
        std::shared_ptr<Matrix> CoccLB_ortho = linalg::triplet(ortho_targetB, matrices_["S"], thislinkB, true, false, false);
        // CoccLB_ortho->print();
        double** CoccLB_orthop = CoccLB_ortho->pointer();
        for (int j = 0; j < norbB; j++) {
            di[0] = j;
            dj[0] = j + 1;
            Slice colslice(di,dj);
            thislinkB->axpy(-CoccLB_orthop[j][0],ortho_targetB->get_block(rowslice,colslice));
        }
    }
    thisnormB = linalg::triplet(thislinkB, matrices_["S"], thislinkB, true, false, false);
    double** thisnormBp = thisnormB->pointer();
    thislinkB->scale(1.0/std::sqrt(thisnormBp[0][0]));

// We can now use the link IHOs from the second iteration to update the fragment density matrices D_A, D_B, D_C.
// Note the multiplication by 1/sqrt{2} which correspond to multiplying the density matrix contribution by 1/2 (the orbital is singly occupied).
    D_A = linalg::doublet(matrices_["Cocc0A"], matrices_["Cocc0A"], false, true);
    D_B = linalg::doublet(matrices_["Cocc0B"], matrices_["Cocc0B"], false, true);
    D_C->zero();
    if (Cocc_C->colspi()[0] > 0) {
        D_C = linalg::doublet(Cocc_C, Cocc_C, false, true);
    }
 
    thislinkA->scale(0.7071067811865475244);
    thislinkB->scale(0.7071067811865475244);
    D_A->axpy(1.0,linalg::doublet(thislinkA, thislinkA, false, true)); 
    D_B->axpy(1.0,linalg::doublet(thislinkB, thislinkB, false, true)); 
    if (Cocc_C->colspi()[0] > 0) {
        D_C->axpy(-1.0,linalg::doublet(thislinkA, thislinkA, false, true)); 
        D_C->axpy(-1.0,linalg::doublet(thislinkB, thislinkB, false, true)); 
    }
 
    matrices_["D_A"] = D_A;
    matrices_["D_B"] = D_B;
    matrices_["D_C"] = D_C;
 
    // PA and PB are used only to define the complement of DA and DB in the DCBS
    P_A =
        linalg::doublet(reference_->Ca_subset("AO", "ALL"), reference_->Ca_subset("AO", "ALL"), false, true);
    P_A->subtract(D_A);
 
    P_B =
        linalg::doublet(reference_->Ca_subset("AO", "ALL"), reference_->Ca_subset("AO", "ALL"), false, true);
    P_B->subtract(D_B);
 
    matrices_["P_A"] = P_A;
    matrices_["P_B"] = P_B;
 
    matrices_["Cocc_A"] = matrices_["Cocc0A"];
    matrices_["Cocc_B"] = matrices_["Cocc0B"];

    J_A->copy(matrices_["J0A"]);
    K_A->copy(matrices_["K0A"]);
    J_B->copy(matrices_["J0B"]);
    K_B->copy(matrices_["K0B"]);
    J_C->copy(matrices_["JC"]);
    K_C->copy(matrices_["KC"]);

// Finally, redo the Coulomb and exchange operators one last time.
    AlloccA = linalg::horzcat({matrices_["Cocc0A"],thislinkA});
    AlloccB = linalg::horzcat({matrices_["Cocc0B"],thislinkB});
 
    Cl = jk_->C_left();
    Cr = jk_->C_right();
 
    Cl.clear();
    Cr.clear();
 
    Cl.push_back(thislinkA);
    Cr.push_back(thislinkA);
    Cl.push_back(thislinkB);
    Cr.push_back(thislinkB);
 
    jk_->compute();
 
    J_A->add(J[0]);
    K_A->add(K[0]);
    J_B->add(J[1]);
    K_B->add(K[1]);
    J_C->subtract(J[0]);
    K_C->subtract(K[0]);
    J_C->subtract(J[1]);
    K_C->subtract(K[1]);
    matrices_["JLA"] = J[0]->clone(); 
    matrices_["KLA"] = K[0]->clone(); 
    matrices_["JLB"] = J[1]->clone(); 
    matrices_["KLB"] = K[1]->clone(); 
 
    matrices_["AlloccA"] = AlloccA;
    matrices_["AlloccB"] = AlloccB;
    matrices_["thislinkA"] = thislinkA;
    matrices_["thislinkB"] = thislinkB;
    matrices_["J_A"] = J_A;
    matrices_["J_B"] = J_B;
    matrices_["K_A"] = K_A;
    matrices_["K_B"] = K_B;
    matrices_["J_C"] = J_C;
    matrices_["K_C"] = K_C;
    }
    }
    }
// Now we are done even for SAO2/SIAO2.
}

// Create cube files for displaying link orbitals and densities.
// Three kinds of cubes are supported: 
// FISAPT_CUBE_LINKIBOS==true - compute the two link IBOs
// FISAPT_CUBE_LINKIHOS==true - compute the two link IHOs obtained by projecting link IBOs onto a suitable fragment space
// FISAPT_CUBE_DENSMAT==true - compute the three densities for fragments A/B/C
void FISAPT::do_cubes() {

    if ((!(options_.get_bool("FISAPT_CUBE_LINKIHOS"))) && (!(options_.get_bool("FISAPT_CUBE_LINKIBOS"))) && (!(options_.get_bool("FISAPT_CUBE_DENSMAT")))) {
    return;
    }
    std::shared_ptr<BasisSet> auxiliary = reference_->get_basisset("DF_BASIS_SAPT");
    std::shared_ptr<CubicScalarGrid> grid_ = std::make_shared<CubicScalarGrid>(primary_, options_);
    grid_->set_filepath(options_.get_str("CUBEPROP_FILEPATH"));
    grid_->set_auxiliary_basis(auxiliary);

    if (options_.get_bool("FISAPT_CUBE_LINKIHOS")) {
        std::vector<int> indices;
        indices.push_back(0);
        std::vector<std::string> labels1;
        labels1.push_back("LinkA");
        std::vector<std::string> labels2;
        labels2.push_back("LinkB");
        grid_->compute_orbitals(matrices_["thislinkA"], indices, labels1, "LinkA");
        grid_->compute_orbitals(matrices_["thislinkB"], indices, labels2, "LinkB");
    }

    if (options_.get_bool("FISAPT_CUBE_LINKIBOS")) {
        std::vector<std::string> labels3;
        std::vector<int> indices2;
        for (int k = 0; k < matrices_["LoccL"]->colspi()[0]; k++) {
            labels3.push_back("LoccL");
            indices2.push_back(k);
        }
        grid_->compute_orbitals(matrices_["LoccL"], indices2, labels3, "LoccL");
    }

    if (options_.get_bool("FISAPT_CUBE_DENSMAT")) {    
        std::shared_ptr<Matrix> D_A = matrices_["D_A"];
        std::shared_ptr<Matrix> D_B = matrices_["D_B"];
        std::shared_ptr<Matrix> D_C = matrices_["D_C"];
        grid_->compute_density(D_A, "DA");
        grid_->compute_density(D_B, "DB");
        grid_->compute_density(D_C, "DC");
    }

}

// Compute deltaHF contribution (counted as induction)
void FISAPT::dHF() {
    outfile->Printf("  ==> dHF <==\n\n");

    // => Pointers <= //

    std::shared_ptr<Matrix> T = matrices_["T"];

    std::shared_ptr<Matrix> D_A = matrices_["D_A"];
    std::shared_ptr<Matrix> D_B = matrices_["D_B"];
    std::shared_ptr<Matrix> D_C = matrices_["D_C"];
    std::shared_ptr<Matrix> V_A = matrices_["V_A"];
    std::shared_ptr<Matrix> V_B = matrices_["V_B"];
    std::shared_ptr<Matrix> V_C = matrices_["V_C"];
    std::shared_ptr<Matrix> J_A = matrices_["J_A"];
    std::shared_ptr<Matrix> J_B = matrices_["J_B"];
    std::shared_ptr<Matrix> J_C = matrices_["J_C"];
    std::shared_ptr<Matrix> K_A = matrices_["K_A"];
    std::shared_ptr<Matrix> K_B = matrices_["K_B"];
    std::shared_ptr<Matrix> K_C = matrices_["K_C"];

    std::string link_assignment = options_.get_str("FISAPT_LINK_ASSIGNMENT");
    double** Enuc2p = matrices_["E NUC"]->pointer();

    // => Dimer HF (Already done) <= //

    double EABC = reference_->energy();

    // => Monomer AC Energy <= //

    double EAC = 0.0;
    EAC += Enuc2p[0][0];
    EAC += Enuc2p[2][2];
    EAC += Enuc2p[0][2];
    EAC += Enuc2p[2][0];

    std::shared_ptr<Matrix> H_AC(D_A->clone());
    H_AC->copy(T);
    H_AC->add(V_A);
    H_AC->add(V_C);
    if (reference_->has_potential_variable("C")) H_AC->add(matrices_["VE"]);

    std::shared_ptr<Matrix> F_AC(D_A->clone());
    F_AC->copy(H_AC);
    F_AC->add(J_A);
    F_AC->add(J_A);
    F_AC->add(J_C);
    F_AC->add(J_C);
    F_AC->subtract(K_A);
    F_AC->subtract(K_C);

    std::shared_ptr<Matrix> D_AC(D_A->clone());
    D_AC->copy(D_A);
    D_AC->add(D_C);
    outfile->Printf("\n Trace of D_AC: %10.6f \n",D_AC->vector_dot(matrices_["S"]));

    EAC += D_AC->vector_dot(H_AC) + D_AC->vector_dot(F_AC);

    H_AC.reset();
    F_AC.reset();
    D_AC.reset();

    // => Monomer BC Energy <= //

    double EBC = 0.0;
    EBC += Enuc2p[1][1];
    EBC += Enuc2p[2][2];
    EBC += Enuc2p[1][2];
    EBC += Enuc2p[2][1];

    std::shared_ptr<Matrix> H_BC(D_B->clone());
    H_BC->copy(T);
    H_BC->add(V_B);
    H_BC->add(V_C);
    if (reference_->has_potential_variable("C")) H_BC->add(matrices_["VE"]);

    std::shared_ptr<Matrix> F_BC(D_B->clone());
    F_BC->copy(H_BC);
    F_BC->add(J_B);
    F_BC->add(J_B);
    F_BC->add(J_C);
    F_BC->add(J_C);
    F_BC->subtract(K_B);
    F_BC->subtract(K_C);

    std::shared_ptr<Matrix> D_BC(D_B->clone());
    D_BC->copy(D_B);
    D_BC->add(D_C);
    outfile->Printf("\n Trace of D_BC: %10.6f \n",D_BC->vector_dot(matrices_["S"]));

    EBC += D_BC->vector_dot(H_BC) + D_BC->vector_dot(F_BC);

    H_BC.reset();
    F_BC.reset();
    D_BC.reset();

    // We also compute the energies of the isolated A and B
    // fragments, so that the orbital deformation energy is available

    // => Monomer A Energy <= //

    double EA = 0.0;
    EA += Enuc2p[0][0];

    std::shared_ptr<Matrix> H_A(D_A->clone());
    H_A->copy(T);
    H_A->add(V_A);
    if (reference_->has_potential_variable("C")) H_A->add(matrices_["VE"]);

    std::shared_ptr<Matrix> F_A(D_A->clone());
    F_A->copy(H_A);
    F_A->add(J_A);
    F_A->add(J_A);
    F_A->subtract(K_A);

    outfile->Printf("\n Trace of D_A: %10.6f \n",D_A->vector_dot(matrices_["S"]));
    EA += D_A->vector_dot(H_A) + D_A->vector_dot(F_A);

    // => Monomer B Energy <= //

    double EB = 0.0;
    EB += Enuc2p[1][1];

    std::shared_ptr<Matrix> H_B(D_B->clone());
    H_B->copy(T);
    H_B->add(V_B);
    if (reference_->has_potential_variable("C")) H_B->add(matrices_["VE"]);

    std::shared_ptr<Matrix> F_B(D_B->clone());
    F_B->copy(H_B);
    F_B->add(J_B);
    F_B->add(J_B);
    F_B->subtract(K_B);

    outfile->Printf("\n Trace of D_B: %10.6f \n",D_B->vector_dot(matrices_["S"]));
    EB += D_B->vector_dot(H_B) + D_B->vector_dot(F_B);

    // => Monomer C Energy <= //

    double EC = 0.0;
    EC += Enuc2p[2][2];

    std::shared_ptr<Matrix> H_C(D_C->clone());
    H_C->copy(T);
    H_C->add(V_C);
    if (reference_->has_potential_variable("C")) H_C->add(matrices_["VE"]);

    std::shared_ptr<Matrix> F_C(D_C->clone());
    F_C->copy(H_C);
    F_C->add(J_C);
    F_C->add(J_C);
    F_C->subtract(K_C);

    outfile->Printf("\n Trace of D_C: %10.6f \n",D_C->vector_dot(matrices_["S"]));
    EC += D_C->vector_dot(H_C) + D_C->vector_dot(F_C);

    // => Delta HF <= //

    double EHF = EABC - EAC - EBC + EC;
 
    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") {

        // => Now do analysis of dipole moment contributions <= //
        std::shared_ptr<Matrix> D_ABC(D_A->clone());
        std::shared_ptr<Matrix> D_ABCA(D_A->clone());
        std::shared_ptr<Matrix> D_ABCB(D_A->clone());
        std::shared_ptr<Matrix> D_ABCLA(D_A->clone());
        std::shared_ptr<Matrix> D_ABCLB(D_A->clone());
        D_ABC = linalg::doublet(matrices_["Locc"], matrices_["Locc"], false, true);
        D_ABCA = linalg::doublet(matrices_["LoccA"], matrices_["LoccA"], false, true);
        D_ABCB = linalg::doublet(matrices_["LoccB"], matrices_["LoccB"], false, true);
        D_ABCLA = linalg::doublet(matrices_["thislinkA"], matrices_["thislinkA"], false, true);
        D_ABCLB = linalg::doublet(matrices_["thislinkB"], matrices_["thislinkB"], false, true);
        std::shared_ptr<Matrix> D0A = linalg::doublet(matrices_["Cocc0A"], matrices_["Cocc0A"], false, true);
        std::shared_ptr<Matrix> D0B = linalg::doublet(matrices_["Cocc0B"], matrices_["Cocc0B"], false, true);

        double dipabc_x = 2.0*D_ABC->vector_dot(matrices_["X DIP"]) + vectors_["DIP NUCA"]->get(0) + vectors_["DIP NUCB"]->get(0) + vectors_["DIP NUCC"]->get(0);
        double dipabc_y = 2.0*D_ABC->vector_dot(matrices_["Y DIP"]) + vectors_["DIP NUCA"]->get(1) + vectors_["DIP NUCB"]->get(1) + vectors_["DIP NUCC"]->get(1);
        double dipabc_z = 2.0*D_ABC->vector_dot(matrices_["Z DIP"]) + vectors_["DIP NUCA"]->get(2) + vectors_["DIP NUCB"]->get(2) + vectors_["DIP NUCC"]->get(2);
        outfile->Printf("\n HF dipole moment for ABC: (%12.8f, %12.8f, %12.8f) length: %12.8f", dipabc_x, dipabc_y, dipabc_z, std::sqrt(dipabc_x*dipabc_x+dipabc_y*dipabc_y+dipabc_z*dipabc_z));
        dipabc_x = 2.0*D_ABCA->vector_dot(matrices_["X DIP"]) + vectors_["DIP NUCA0"]->get(0);
        dipabc_y = 2.0*D_ABCA->vector_dot(matrices_["Y DIP"]) + vectors_["DIP NUCA0"]->get(1);
        dipabc_z = 2.0*D_ABCA->vector_dot(matrices_["Z DIP"]) + vectors_["DIP NUCA0"]->get(2);
        outfile->Printf("\n HF dipole moment for A before link reassignment: (%12.8f, %12.8f, %12.8f) length: %12.8f", dipabc_x, dipabc_y, dipabc_z, std::sqrt(dipabc_x*dipabc_x+dipabc_y*dipabc_y+dipabc_z*dipabc_z));
        dipabc_x = 2.0*D_ABCB->vector_dot(matrices_["X DIP"]) + vectors_["DIP NUCB0"]->get(0);
        dipabc_y = 2.0*D_ABCB->vector_dot(matrices_["Y DIP"]) + vectors_["DIP NUCB0"]->get(1);
        dipabc_z = 2.0*D_ABCB->vector_dot(matrices_["Z DIP"]) + vectors_["DIP NUCB0"]->get(2);
        outfile->Printf("\n HF dipole moment for B before link reassignment: (%12.8f, %12.8f, %12.8f) length: %12.8f", dipabc_x, dipabc_y, dipabc_z, std::sqrt(dipabc_x*dipabc_x+dipabc_y*dipabc_y+dipabc_z*dipabc_z));
        dipabc_x = 2.0*D0A->vector_dot(matrices_["X DIP"]) + 2.0*D_ABCLA->vector_dot(matrices_["X DIP"]) + vectors_["DIP NUCA"]->get(0);
        dipabc_y = 2.0*D0A->vector_dot(matrices_["Y DIP"]) + 2.0*D_ABCLA->vector_dot(matrices_["Y DIP"]) + vectors_["DIP NUCA"]->get(1);
        dipabc_z = 2.0*D0A->vector_dot(matrices_["Z DIP"]) + 2.0*D_ABCLA->vector_dot(matrices_["Z DIP"]) + vectors_["DIP NUCA"]->get(2);
        outfile->Printf("\n HF dipole moment for A after link reassignment:  (%12.8f, %12.8f, %12.8f) length: %12.8f", dipabc_x, dipabc_y, dipabc_z, std::sqrt(dipabc_x*dipabc_x+dipabc_y*dipabc_y+dipabc_z*dipabc_z));
        dipabc_x = 2.0*D0B->vector_dot(matrices_["X DIP"]) + 2.0*D_ABCLB->vector_dot(matrices_["X DIP"]) + vectors_["DIP NUCB"]->get(0);
        dipabc_y = 2.0*D0B->vector_dot(matrices_["Y DIP"]) + 2.0*D_ABCLB->vector_dot(matrices_["Y DIP"]) + vectors_["DIP NUCB"]->get(1);
        dipabc_z = 2.0*D0B->vector_dot(matrices_["Z DIP"]) + 2.0*D_ABCLB->vector_dot(matrices_["Z DIP"]) + vectors_["DIP NUCB"]->get(2);
        outfile->Printf("\n HF dipole moment for B after link reassignment:  (%12.8f, %12.8f, %12.8f) length: %12.8f\n\n", dipabc_x, dipabc_y, dipabc_z, std::sqrt(dipabc_x*dipabc_x+dipabc_y*dipabc_y+dipabc_z*dipabc_z));
    }

    // => Monomer and dimer energies in the original ABC full system <= //

    // Compute density from A HF localized orbitals
    std::shared_ptr<Matrix> LoccA = matrices_["LoccA"];
    std::shared_ptr<Matrix> LoccB = matrices_["LoccB"];
    std::shared_ptr<Matrix> LD_A = linalg::doublet(LoccA, LoccA, false, true);
    std::shared_ptr<Matrix> LD_B = linalg::doublet(LoccB, LoccB, false, true);

    // Get J and K from A and B HF localized orbitals while we are at it
    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    Cl.clear();
    Cr.clear();

    Cl.push_back(LoccA);
    Cr.push_back(LoccA);
    Cl.push_back(LoccB);
    Cr.push_back(LoccB);

    jk_->compute();

    std::shared_ptr<Matrix> LJ_A(J[0]->clone());
    std::shared_ptr<Matrix> LK_A(K[0]->clone());
    std::shared_ptr<Matrix> LJ_B(J[1]->clone());
    std::shared_ptr<Matrix> LK_B(K[1]->clone());

    // We have all the ingredients, now we build everything
    // Monomer A localised energy
    double LE_A = 0.0;
    LE_A += Enuc2p[0][0];
    std::shared_ptr<Matrix> LH_A(T->clone());
    LH_A->copy(T);
    LH_A->add(V_A);
    if (reference_->has_potential_variable("C")) LH_A->add(matrices_["VE"]);
    std::shared_ptr<Matrix> LF_A(LH_A->clone());
    LF_A->copy(LH_A);
    LF_A->add(LJ_A);
    LF_A->add(LJ_A);
    LF_A->subtract(LK_A);
    LE_A += LD_A->vector_dot(LH_A) + LD_A->vector_dot(LF_A);
    LH_A.reset();
    LF_A.reset();

    // Monomer B localised energy
    double LE_B = 0.0;
    LE_B += Enuc2p[1][1];
    std::shared_ptr<Matrix> LH_B(T->clone());
    LH_B->copy(T);
    LH_B->add(V_B);
    if (reference_->has_potential_variable("C")) LH_B->add(matrices_["VE"]);
    std::shared_ptr<Matrix> LF_B(LH_B->clone());
    LF_B->copy(LH_B);
    LF_B->add(LJ_B);
    LF_B->add(LJ_B);
    LF_B->subtract(LK_B);
    LE_B += LD_B->vector_dot(LH_B) + LD_B->vector_dot(LF_B);
    LH_B.reset();
    LF_B.reset();

    // Dimer AC localised energy
    double LE_AC = 0.0;
    LE_AC += Enuc2p[0][0];
    LE_AC += Enuc2p[0][2];
    LE_AC += Enuc2p[2][2];
    LE_AC += Enuc2p[2][0];

    std::shared_ptr<Matrix> LH_AC(T->clone());
    LH_AC->copy(T);
    LH_AC->add(V_A);
    LH_AC->add(V_C);
    if (reference_->has_potential_variable("C")) LH_AC->add(matrices_["VE"]);
    std::shared_ptr<Matrix> LF_AC(LH_AC->clone());
    LF_AC->copy(LH_AC);
    LF_AC->add(J_C);
    LF_AC->add(J_C);
    LF_AC->add(LJ_A);
    LF_AC->add(LJ_A);
    LF_AC->subtract(LK_A);
    LF_AC->subtract(K_C);
    std::shared_ptr<Matrix> LD_AC(LD_A->clone());
    LD_AC->copy(LD_A);
    LD_AC->add(D_C);

    LE_AC += LD_AC->vector_dot(LH_AC) + LD_AC->vector_dot(LF_AC);
    LD_AC.reset();
    LH_AC.reset();
    LF_AC.reset();

    // Dimer BC localised energy
    double LE_BC = 0.0;
    LE_BC += Enuc2p[1][1];
    LE_BC += Enuc2p[1][2];
    LE_BC += Enuc2p[2][2];
    LE_BC += Enuc2p[2][1];

    std::shared_ptr<Matrix> LH_BC(T->clone());
    LH_BC->copy(T);
    LH_BC->add(V_B);
    LH_BC->add(V_C);
    if (reference_->has_potential_variable("C")) LH_BC->add(matrices_["VE"]);
    std::shared_ptr<Matrix> LF_BC(LH_BC->clone());
    LF_BC->copy(LH_BC);
    LF_BC->add(J_C);
    LF_BC->add(J_C);
    LF_BC->add(LJ_B);
    LF_BC->add(LJ_B);
    LF_BC->subtract(LK_B);
    LF_BC->subtract(K_C);
    std::shared_ptr<Matrix> LD_BC(LD_B->clone());
    LD_BC->copy(LD_B);
    LD_BC->add(D_C);

    LE_BC += LD_BC->vector_dot(LH_BC) + LD_BC->vector_dot(LF_BC);
    LD_BC.reset();
    LH_BC.reset();
    LF_BC.reset();

    // Dimer AB localised energy
    double LE_BA = 0.0;
    LE_BA += Enuc2p[1][1];
    LE_BA += Enuc2p[1][0];
    LE_BA += Enuc2p[0][0];
    LE_BA += Enuc2p[0][1];

    std::shared_ptr<Matrix> LH_BA(T->clone());
    LH_BA->copy(T);
    LH_BA->add(V_B);
    LH_BA->add(V_A);
    if (reference_->has_potential_variable("C")) LH_BA->add(matrices_["VE"]);
    std::shared_ptr<Matrix> LF_BA(LH_BA->clone());
    LF_BA->copy(LH_BA);
    LF_BA->add(LJ_A);
    LF_BA->add(LJ_A);
    LF_BA->add(LJ_B);
    LF_BA->add(LJ_B);
    LF_BA->subtract(LK_B);
    LF_BA->subtract(LK_A);
    std::shared_ptr<Matrix> LD_BA(LD_B->clone());
    LD_BA->copy(LD_B);
    LD_BA->add(LD_A);

    LE_BA += LD_BA->vector_dot(LH_BA) + LD_BA->vector_dot(LF_BA);
    LD_BA.reset();
    LH_BA.reset();
    LF_BA.reset();

    // => Print <= //

    outfile->Printf("    E ABC(HF) = %24.16E [Eh]\n", EABC);
    outfile->Printf("    E AC(0)   = %24.16E [Eh]\n", EAC);
    outfile->Printf("    E BC(0)   = %24.16E [Eh]\n", EBC);
    outfile->Printf("    E A(0)    = %24.16E [Eh]\n", EA);
    outfile->Printf("    E B(0)    = %24.16E [Eh]\n", EB);
    outfile->Printf("    E AC(HF)  = %24.16E [Eh]\n", LE_AC);
    outfile->Printf("    E BC(HF)  = %24.16E [Eh]\n", LE_BC);
    outfile->Printf("    E AB(HF)  = %24.16E [Eh]\n", LE_BA);
    outfile->Printf("    E A(HF)   = %24.16E [Eh]\n", LE_A);
    outfile->Printf("    E B(HF)   = %24.16E [Eh]\n", LE_B);
    outfile->Printf("    E C       = %24.16E [Eh]\n", EC);
    outfile->Printf("    E HF      = %24.16E [Eh]\n", EHF);
    outfile->Printf("\n");

    scalars_["HF"] = EHF;
    scalars_["E_A"] = EA;
    scalars_["E_B"] = EB;

    // Export all components of dHF as Psi4 variables
    scalars_["E_C"] = EC;
    scalars_["E_AC"] = EAC;
    scalars_["E_BC"] = EBC;
    scalars_["E_ABC_HF"] = EABC;
    scalars_["E_AC_HF"] = LE_AC;
    scalars_["E_BC_HF"] = LE_BC;
    scalars_["E_AB_HF"] = LE_BA;
    scalars_["E_A_HF"] = LE_A;
    scalars_["E_B_HF"] = LE_B;
}

// Compute total electrostatics contribution
// This definition works just as well for the SAOn/SIAOn methods because the densities and Coulomb operators have been adjusted.
void FISAPT::elst() {
    outfile->Printf("  ==> Electrostatics <==\n\n");

    std::shared_ptr<Matrix> D_A = matrices_["D_A"];
    std::shared_ptr<Matrix> D_B = matrices_["D_B"];
    std::shared_ptr<Matrix> V_A = matrices_["V_A"];
    std::shared_ptr<Matrix> V_B = matrices_["V_B"];
    std::shared_ptr<Matrix> J_A = matrices_["J_A"];
    std::shared_ptr<Matrix> J_B = matrices_["J_B"];

    double Enuc = 0.0;
    double** Enuc2p = matrices_["E NUC"]->pointer();
    Enuc += 2.0 * Enuc2p[0][1];  // A - B

    double Elst10 = 0.0;
    double Elst10_sum = 0.0;
    std::vector<double> Elst10_terms;
    Elst10_terms.resize(4);
    Elst10_terms[0] += 2.0 * D_A->vector_dot(V_B);
    Elst10_sum += 2.0 * D_A->vector_dot(V_B);
    outfile->Printf("Elst10: %18.12lf\n", Elst10_sum);
    Elst10_terms[1] += 2.0 * D_B->vector_dot(V_A);
    Elst10_sum += 2.0 * D_B->vector_dot(V_A);
    outfile->Printf("Elst10: %18.12lf\n", Elst10_sum);
    Elst10_terms[2] += 4.0 * D_A->vector_dot(J_B);
    Elst10_sum += 4.0 * D_A->vector_dot(J_B);
    outfile->Printf("Elst10: %18.12lf\n", Elst10_sum);
    Elst10_terms[3] += 1.0 * Enuc;
    Elst10_sum += 1.0 * Enuc;
    outfile->Printf("Elst10: %18.12lf\n", Elst10_sum);
    // correct
    outfile->Printf("2.0 * D_A->vector_dot(V_B)  (HF): %18.12lf\n",2.0 * D_A->vector_dot(V_B));
    outfile->Printf("2.0 * D_B->vector_dot(V_A)  (HF): %18.12lf\n",2.0 * D_B->vector_dot(V_A));
    // maybe?
    outfile->Printf("4.0 * D_A->vector_dot(J_B): %18.12lf\n",4.0 * D_A->vector_dot(J_B));
    outfile->Printf("Enuc_t: %18.12lf\n", Enuc);
    for (int k = 0; k < Elst10_terms.size(); k++) {
        Elst10 += Elst10_terms[k];
    }
    // for (int k = 0; k < Elst10_terms.size(); k++) {
    //    outfile->Printf("    Elst10,r (%1d)        = %18.12lf [Eh]\n",k+1,Elst10_terms[k]);
    //}

    scalars_["Elst10,r"] = Elst10;
    outfile->Printf("    Elst10,r            = %18.12lf [Eh]\n", Elst10);

    // External potential interactions
    if (reference_->has_potential_variable("A") && reference_->has_potential_variable("B")) {
        // Add the interaction between the external potentials  in A and B
        // Multiply by 2 to get the full A-B + B-A interaction energy
        scalars_["Extern-Extern"] = matrices_["extern_extern_IE"]->get(0, 1)*2.0; 
        outfile->Printf("    Extern-Extern       = %18.12lf [Eh]\n", scalars_["Extern-Extern"]);
    }

   outfile->Printf("\n");
    // fflush(outfile);
}

// Compute total exchange contribution
// Special code (https://doi.org/10.1021/acs.jpca.2c06465) is used for the SAOn/SIAOn variants - see below.
void FISAPT::exch() {
    outfile->Printf("  ==> Exchange <==\n\n");

    // => Density and Potential Matrices <= //

    std::shared_ptr<Matrix> S = matrices_["S"];

    std::shared_ptr<Matrix> D_A = matrices_["D_A"];
    std::shared_ptr<Matrix> D_B = matrices_["D_B"];
    std::shared_ptr<Matrix> P_A = matrices_["P_A"];
    std::shared_ptr<Matrix> P_B = matrices_["P_B"];
    std::shared_ptr<Matrix> V_A = matrices_["V_A"];
    std::shared_ptr<Matrix> V_B = matrices_["V_B"];
    std::shared_ptr<Matrix> J_A = matrices_["J_A"];
    std::shared_ptr<Matrix> J_B = matrices_["J_B"];
    std::shared_ptr<Matrix> K_A = matrices_["K_A"];
    std::shared_ptr<Matrix> K_B = matrices_["K_B"];

    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    // ==> Exchange Terms (S^2, MCBS or DCBS) <== //

    std::string link_assignment = options_.get_str("FISAPT_LINK_ASSIGNMENT");
    if (link_assignment == "C" || link_assignment == "AB" ) {

    std::shared_ptr<Matrix> Cocc_A = matrices_["Cocc_A"];   
    std::shared_ptr<Matrix> Cocc_B = matrices_["Cocc_B"];

    std::shared_ptr<Matrix> C_O = linalg::triplet(D_B, S, Cocc_A);
    Cl.clear();
    Cr.clear();
    Cl.push_back(Cocc_A);
    Cr.push_back(C_O);
    jk_->compute();
    std::shared_ptr<Matrix> K_O = K[0];

    double Exch10_2M = 0.0;
    std::vector<double> Exch10_2M_terms;
    Exch10_2M_terms.resize(6);
    Exch10_2M_terms[0] -= 2.0 * D_A->vector_dot(K_B);
    Exch10_2M_terms[1] -= 2.0 * linalg::triplet(D_A, S, D_B)->vector_dot(V_A);
    Exch10_2M_terms[1] -= 4.0 * linalg::triplet(D_A, S, D_B)->vector_dot(J_A);
    Exch10_2M_terms[1] += 2.0 * linalg::triplet(D_A, S, D_B)->vector_dot(K_A);
    Exch10_2M_terms[2] -= 2.0 * linalg::triplet(D_B, S, D_A)->vector_dot(V_B);
    Exch10_2M_terms[2] -= 4.0 * linalg::triplet(D_B, S, D_A)->vector_dot(J_B);
    Exch10_2M_terms[2] += 2.0 * linalg::triplet(D_B, S, D_A)->vector_dot(K_B);
    Exch10_2M_terms[3] += 2.0 * linalg::triplet(linalg::triplet(D_B, S, D_A), S, D_B)->vector_dot(V_A);
    Exch10_2M_terms[3] += 4.0 * linalg::triplet(linalg::triplet(D_B, S, D_A), S, D_B)->vector_dot(J_A);
    Exch10_2M_terms[4] += 2.0 * linalg::triplet(linalg::triplet(D_A, S, D_B), S, D_A)->vector_dot(V_B);
    Exch10_2M_terms[4] += 4.0 * linalg::triplet(linalg::triplet(D_A, S, D_B), S, D_A)->vector_dot(J_B);
    Exch10_2M_terms[5] -= 2.0 * linalg::triplet(D_A, S, D_B)->vector_dot(K_O);
    for (int k = 0; k < Exch10_2M_terms.size(); k++) {
        Exch10_2M += Exch10_2M_terms[k];
    }
    // for (int k = 0; k < Exch10_2M_terms.size(); k++) {
    //     outfile->Printf("    Exch10(S^2) (%1d)     = %18.12lf [Eh]\n",k+1,Exch10_2M_terms[k]);
    // }
    // scalars_["Exch10(S^2)"] = Exch10_2;
    outfile->Printf("    Exch10(S^2) [MCBS]  = %18.12lf [Eh]\n",Exch10_2M);
    outfile->Printf("    Exch10(S^2)         = %18.12lf [Eh]\n",Exch10_2M);
    // fflush(outfile);

    // ==> Exchange Terms (S^2, DCBS only) <== //
    // DCBS is redundant because we have already computed MCBS. Commenting out code that computes DCBS below.
    // => K_AS <= //

    //std::shared_ptr<Matrix> C_AS = linalg::triplet(P_B, S, Cocc_A);
    //Cl.clear();
    //Cr.clear();
    //Cl.push_back(Cocc_A);
    //Cr.push_back(C_AS);
    //jk_->compute();
    //std::shared_ptr<Matrix> K_AS = K[0];

    // => Accumulation <= //

    //double Exch10_2 = 0.0;
    //std::vector<double> Exch10_2_terms;
    //Exch10_2_terms.resize(3);
    //Exch10_2_terms[0] -= 2.0 * linalg::triplet(linalg::triplet(D_A, S, D_B), S, P_A)->vector_dot(V_B);
    //Exch10_2_terms[0] -= 4.0 * linalg::triplet(linalg::triplet(D_A, S, D_B), S, P_A)->vector_dot(J_B);
    //Exch10_2_terms[1] -= 2.0 * linalg::triplet(linalg::triplet(D_B, S, D_A), S, P_B)->vector_dot(V_A);
    //Exch10_2_terms[1] -= 4.0 * linalg::triplet(linalg::triplet(D_B, S, D_A), S, P_B)->vector_dot(J_A);
    //Exch10_2_terms[2] -= 2.0 * linalg::triplet(P_A, S, D_B)->vector_dot(K_AS);
    //for (int k = 0; k < Exch10_2_terms.size(); k++) {
    //    Exch10_2 += Exch10_2_terms[k];
    //}
    //for (int k = 0; k < Exch10_2_terms.size(); k++) {
    //   outfile->Printf("    Exch10(S^2) (%1d)     = %18.12lf [Eh]\n",k+1,Exch10_2_terms[k]);
    //}
    //scalars_["Exch10(S^2)"] = Exch10_2;
    //outfile->Printf("    Exch10(S^2) [DCBS]  = %18.12lf [Eh]\n",Exch10_2);
    //outfile->Printf("    Exch10(S^2)         = %18.12lf [Eh]\n", Exch10_2);
    // fflush(outfile);

    // ==> Exchange Terms (S^\infty, MCBS or DCBS) <== //

    // => T Matrix <= //

    int na = matrices_["Cocc0A"]->colspi()[0];
    int nb = matrices_["Cocc0B"]->colspi()[0];
    int nbf = matrices_["Cocc0A"]->rowspi()[0];

    std::shared_ptr<Matrix> Sab = linalg::triplet(matrices_["Cocc0A"], S, matrices_["Cocc0B"], true, false, false);
    double** Sabp = Sab->pointer();
    auto T = std::make_shared<Matrix>("T", na + nb, na + nb);
    T->identity();
    double** Tp = T->pointer();
    for (int a = 0; a < na; a++) {
        for (int b = 0; b < nb; b++) {
            Tp[a][b + na] = Tp[b + na][a] = Sabp[a][b];
        }
    }
    // T->print();
    T->power(-1.0, 1.0E-12);
    Tp = T->pointer();
    for (int a = 0; a < na + nb; a++) {
        Tp[a][a] -= 1.0;
    }
    // T->print();

    auto C_T_A_n = std::make_shared<Matrix>("C_T_A_n", nbf, na);
    auto C_T_B_n = std::make_shared<Matrix>("C_T_A_n", nbf, nb);
    auto C_T_BA_n = std::make_shared<Matrix>("C_T_BA_n", nbf, nb);
    auto C_T_AB_n = std::make_shared<Matrix>("C_T_AB_n", nbf, na);

    C_DGEMM('N', 'N', nbf, na, na, 1.0, matrices_["Cocc0A"]->pointer()[0], na, &Tp[0][0], na + nb, 0.0,
            C_T_A_n->pointer()[0], na);
    C_DGEMM('N', 'N', nbf, nb, nb, 1.0, matrices_["Cocc0B"]->pointer()[0], nb, &Tp[na][na], na + nb, 0.0,
            C_T_B_n->pointer()[0], nb);
    C_DGEMM('N', 'N', nbf, nb, na, 1.0, matrices_["Cocc0A"]->pointer()[0], na, &Tp[0][na], na + nb, 0.0,
            C_T_BA_n->pointer()[0], nb);
    C_DGEMM('N', 'N', nbf, na, nb, 1.0, matrices_["Cocc0B"]->pointer()[0], nb, &Tp[na][0], na + nb, 0.0,
            C_T_AB_n->pointer()[0], na);

    // => K Terms <= //

    Cl.clear();
    Cr.clear();
    // J/K[T^A, S^\infty]
    Cl.push_back(matrices_["Cocc0A"]);
    Cr.push_back(C_T_A_n);
    // J/K[T^AB, S^\infty]
    Cl.push_back(matrices_["Cocc0A"]);
    Cr.push_back(C_T_AB_n);

    jk_->compute();

    std::shared_ptr<Matrix> J_T_A_n = J[0];
    std::shared_ptr<Matrix> K_T_A_n = K[0];
    std::shared_ptr<Matrix> J_T_AB_n = J[1];
    std::shared_ptr<Matrix> K_T_AB_n = K[1];

    std::shared_ptr<Matrix> T_A_n = linalg::doublet(matrices_["Cocc0A"], C_T_A_n, false, true);
    std::shared_ptr<Matrix> T_B_n = linalg::doublet(matrices_["Cocc0B"], C_T_B_n, false, true);
    std::shared_ptr<Matrix> T_BA_n = linalg::doublet(matrices_["Cocc0B"], C_T_BA_n, false, true);
    std::shared_ptr<Matrix> T_AB_n = linalg::doublet(matrices_["Cocc0A"], C_T_AB_n, false, true);

    double Exch10_n = 0.0;
    std::vector<double> Exch10_n_terms;
    Exch10_n_terms.resize(9);
    Exch10_n_terms[0] -= 2.0 * D_A->vector_dot(K_B);  // This needs to be the full D_A
    Exch10_n_terms[1] += 2.0 * T_A_n->vector_dot(V_B);
    Exch10_n_terms[1] += 4.0 * T_A_n->vector_dot(J_B);
    Exch10_n_terms[1] -= 2.0 * T_A_n->vector_dot(K_B);
    Exch10_n_terms[2] += 2.0 * T_B_n->vector_dot(V_A);
    Exch10_n_terms[2] += 4.0 * T_B_n->vector_dot(J_A);
    Exch10_n_terms[2] -= 2.0 * T_B_n->vector_dot(K_A);
    Exch10_n_terms[3] += 2.0 * T_AB_n->vector_dot(V_A);
    Exch10_n_terms[3] += 4.0 * T_AB_n->vector_dot(J_A);
    Exch10_n_terms[3] -= 2.0 * T_AB_n->vector_dot(K_A);
    Exch10_n_terms[4] += 2.0 * T_AB_n->vector_dot(V_B);
    Exch10_n_terms[4] += 4.0 * T_AB_n->vector_dot(J_B);
    Exch10_n_terms[4] -= 2.0 * T_AB_n->vector_dot(K_B);
    Exch10_n_terms[5] += 4.0 * T_B_n->vector_dot(J_T_AB_n);
    Exch10_n_terms[5] -= 2.0 * T_B_n->vector_dot(K_T_AB_n);
    Exch10_n_terms[6] += 4.0 * T_A_n->vector_dot(J_T_AB_n);
    Exch10_n_terms[6] -= 2.0 * T_A_n->vector_dot(K_T_AB_n);
    Exch10_n_terms[7] += 4.0 * T_B_n->vector_dot(J_T_A_n);
    Exch10_n_terms[7] -= 2.0 * T_B_n->vector_dot(K_T_A_n);
    Exch10_n_terms[8] += 4.0 * T_AB_n->vector_dot(J_T_AB_n);
    Exch10_n_terms[8] -= 2.0 * T_AB_n->vector_dot(K_T_AB_n);
    for (int k = 0; k < Exch10_n_terms.size(); k++) {
        Exch10_n += Exch10_n_terms[k];
    }
    // for (int k = 0; k < Exch10_n_terms.size(); k++) {
    //    outfile->Printf("    Exch10 (%1d)          = %18.12lf [Eh]\n",k+1,Exch10_n_terms[k]);
    //}
    scalars_["Exch10"] = Exch10_n;
    // outfile->Printf("    Exch10      [MCBS]  = %18.12lf [Eh]\n",Exch10_n);
    outfile->Printf("    Exch10              = %18.12lf [Eh]\n", Exch10_n);
    outfile->Printf("\n");
    // fflush(outfile);
    } // (link_assignment == "C" || link_assignment == "AB" )

// Now come the first-order exchange energy expressions appropriate for the SAOn/SIAOn link assignments.
// The expressions depend on the (parallel, perpendicular, or averaged) spin coupling between the two singly occupied link orbitals.
// If FISAPT_EXCH_PARPERP == true, all 3 spin couplings are calculated.
// If FISAPT_EXCH_PARPERP == false, only the averaged spin coupling is computed, avoiding the computation of a bunch of terms that cancel out in that case.
// See https://doi.org/10.1021/acs.jpca.2c06465 for details - equation numbers below refer to this paper and its SI.
    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2" ) {

// first let's do E(10)exch(S^2) for the parallel and perpendicular spin coupling of link orbitals - Eq. (5)
    std::shared_ptr<Matrix> Cocc_A = matrices_["AlloccA"];
    std::shared_ptr<Matrix> Cocc_B = matrices_["AlloccB"];
    std::shared_ptr<Matrix> Sab = linalg::triplet(matrices_["AlloccA"], S, matrices_["AlloccB"], true, false, false);
    double** Sabp = Sab->pointer();
    std::shared_ptr<Matrix> D_X(D_A->clone());
    std::shared_ptr<Matrix> D_Y(D_A->clone());
// here we need the link-only parts of density matrices D_X,D_Y and their corresponding Coulomb and exchange matrices (computed earlier)
    D_X = linalg::doublet(matrices_["thislinkA"], matrices_["thislinkA"], false, true);
    D_Y = linalg::doublet(matrices_["thislinkB"], matrices_["thislinkB"], false, true);
    std::shared_ptr<Matrix> J_X = matrices_["JLA"];
    std::shared_ptr<Matrix> K_X = matrices_["KLA"];
    std::shared_ptr<Matrix> J_Y = matrices_["JLB"];
    std::shared_ptr<Matrix> K_Y = matrices_["KLB"];
// a few extra J/K matrices still need to be computed here
    std::shared_ptr<Matrix> C_AOB = linalg::triplet(D_B, S, Cocc_A);
    std::shared_ptr<Matrix> C_XOB = linalg::triplet(D_B, S, matrices_["thislinkA"]);
    std::shared_ptr<Matrix> C_AOY = linalg::triplet(D_Y, S, Cocc_A);
    std::shared_ptr<Matrix> C_XOY = linalg::triplet(D_Y, S, matrices_["thislinkA"]);
    Cl.clear();
    Cr.clear();
    Cl.push_back(Cocc_A);
    Cl.push_back(matrices_["thislinkA"]);
    Cl.push_back(Cocc_A);
    Cl.push_back(matrices_["thislinkA"]);
    Cr.push_back(C_AOB);
    Cr.push_back(C_XOB);
    Cr.push_back(C_AOY);
    Cr.push_back(C_XOY);
    jk_->compute();
    std::shared_ptr<Matrix> K_AOB = K[0];
    std::shared_ptr<Matrix> K_XOB = K[1];
    std::shared_ptr<Matrix> K_AOY = K[2];
    std::shared_ptr<Matrix> K_XOY = K[3];

    if (options_.get_bool("FISAPT_EXCH_PARPERP")) {
// we calculate the parallel and perpendicular results separately - Eq. (5) with the upper and lower signs
        double Exch10_S2PAR = 0.0;
        std::vector<double> Exch10_S2PAR_terms;
        Exch10_S2PAR_terms.resize(6);
        double Exch10_S2PERP = 0.0;
        std::vector<double> Exch10_S2PERP_terms;
        Exch10_S2PERP_terms.resize(6);
 
        Exch10_S2PAR_terms[0] -= 2.0 * D_A->vector_dot(K_B);
        Exch10_S2PAR_terms[0] -= 2.0 * D_X->vector_dot(K_Y);
 
        Exch10_S2PAR_terms[1] -= 2.0 * linalg::triplet(D_A, S, D_B)->vector_dot(V_A);
        Exch10_S2PAR_terms[1] -= 2.0 * linalg::triplet(D_X, S, D_Y)->vector_dot(V_A);
        Exch10_S2PAR_terms[1] -= 4.0 * linalg::triplet(D_A, S, D_B)->vector_dot(J_A);
        Exch10_S2PAR_terms[1] -= 4.0 * linalg::triplet(D_X, S, D_Y)->vector_dot(J_A);
        Exch10_S2PAR_terms[1] += 2.0 * linalg::triplet(D_A, S, D_B)->vector_dot(K_A);
        Exch10_S2PAR_terms[1] += 2.0 * linalg::triplet(D_X, S, D_B)->vector_dot(K_X);
        Exch10_S2PAR_terms[1] += 2.0 * linalg::triplet(D_A, S, D_Y)->vector_dot(K_X);
        Exch10_S2PAR_terms[1] += 2.0 * linalg::triplet(D_X, S, D_Y)->vector_dot(K_A);
 
        Exch10_S2PAR_terms[2] -= 2.0 * linalg::triplet(D_B, S, D_A)->vector_dot(V_B);
        Exch10_S2PAR_terms[2] -= 2.0 * linalg::triplet(D_Y, S, D_X)->vector_dot(V_B);
        Exch10_S2PAR_terms[2] -= 4.0 * linalg::triplet(D_B, S, D_A)->vector_dot(J_B);
        Exch10_S2PAR_terms[2] -= 4.0 * linalg::triplet(D_Y, S, D_X)->vector_dot(J_B);
        Exch10_S2PAR_terms[2] += 2.0 * linalg::triplet(D_B, S, D_A)->vector_dot(K_B);
        Exch10_S2PAR_terms[2] += 2.0 * linalg::triplet(D_B, S, D_X)->vector_dot(K_Y);
        Exch10_S2PAR_terms[2] += 2.0 * linalg::triplet(D_Y, S, D_A)->vector_dot(K_Y);
        Exch10_S2PAR_terms[2] += 2.0 * linalg::triplet(D_Y, S, D_X)->vector_dot(K_B);

        Exch10_S2PAR_terms[3] += 2.0 * linalg::triplet(linalg::triplet(D_B, S, D_A), S, D_B)->vector_dot(V_A);
        Exch10_S2PAR_terms[3] += 2.0 * linalg::triplet(linalg::triplet(D_B, S, D_X), S, D_Y)->vector_dot(V_A);
        Exch10_S2PAR_terms[3] += 2.0 * linalg::triplet(linalg::triplet(D_Y, S, D_X), S, D_B)->vector_dot(V_A);
        Exch10_S2PAR_terms[3] += 2.0 * linalg::triplet(linalg::triplet(D_Y, S, D_A), S, D_Y)->vector_dot(V_A);
        Exch10_S2PAR_terms[3] += 4.0 * linalg::triplet(linalg::triplet(D_B, S, D_A), S, D_B)->vector_dot(J_A);
        Exch10_S2PAR_terms[3] += 4.0 * linalg::triplet(linalg::triplet(D_Y, S, D_X), S, D_B)->vector_dot(J_A);
        Exch10_S2PAR_terms[3] += 4.0 * linalg::triplet(linalg::triplet(D_B, S, D_X), S, D_Y)->vector_dot(J_A);
        Exch10_S2PAR_terms[3] += 4.0 * linalg::triplet(linalg::triplet(D_Y, S, D_A), S, D_Y)->vector_dot(J_A);
      
        Exch10_S2PAR_terms[4] += 2.0 * linalg::triplet(linalg::triplet(D_A, S, D_B), S, D_A)->vector_dot(V_B);
        Exch10_S2PAR_terms[4] += 2.0 * linalg::triplet(linalg::triplet(D_X, S, D_Y), S, D_A)->vector_dot(V_B);
        Exch10_S2PAR_terms[4] += 2.0 * linalg::triplet(linalg::triplet(D_A, S, D_Y), S, D_X)->vector_dot(V_B);
        Exch10_S2PAR_terms[4] += 2.0 * linalg::triplet(linalg::triplet(D_X, S, D_B), S, D_X)->vector_dot(V_B);
        Exch10_S2PAR_terms[4] += 4.0 * linalg::triplet(linalg::triplet(D_A, S, D_B), S, D_A)->vector_dot(J_B);
        Exch10_S2PAR_terms[4] += 4.0 * linalg::triplet(linalg::triplet(D_X, S, D_Y), S, D_A)->vector_dot(J_B);
        Exch10_S2PAR_terms[4] += 4.0 * linalg::triplet(linalg::triplet(D_A, S, D_Y), S, D_X)->vector_dot(J_B);
        Exch10_S2PAR_terms[4] += 4.0 * linalg::triplet(linalg::triplet(D_X, S, D_B), S, D_X)->vector_dot(J_B);
      
        Exch10_S2PAR_terms[5] -= 2.0 * linalg::triplet(D_A, S, D_B)->vector_dot(K_AOB);
        Exch10_S2PAR_terms[5] -= 2.0 * linalg::triplet(D_X, S, D_B)->vector_dot(K_XOB);
        Exch10_S2PAR_terms[5] -= 2.0 * linalg::triplet(D_A, S, D_Y)->vector_dot(K_XOB);
        Exch10_S2PAR_terms[5] -= 2.0 * linalg::triplet(D_A, S, D_B)->vector_dot(K_XOY);
        Exch10_S2PAR_terms[5] -= 2.0 * linalg::triplet(D_X, S, D_Y)->vector_dot(K_AOB);
        Exch10_S2PAR_terms[5] -= 2.0 * linalg::triplet(D_X, S, D_B)->vector_dot(K_AOY);
        Exch10_S2PAR_terms[5] -= 2.0 * linalg::triplet(D_A, S, D_Y)->vector_dot(K_AOY);
        Exch10_S2PAR_terms[5] -= 2.0 * linalg::triplet(D_X, S, D_Y)->vector_dot(K_XOY);
       
        Exch10_S2PERP_terms[0] -= 2.0 * D_A->vector_dot(K_B);
        Exch10_S2PERP_terms[0] += 2.0 * D_X->vector_dot(K_Y);
       
        Exch10_S2PERP_terms[1] -= 2.0 * linalg::triplet(D_A, S, D_B)->vector_dot(V_A);
        Exch10_S2PERP_terms[1] += 2.0 * linalg::triplet(D_X, S, D_Y)->vector_dot(V_A);
        Exch10_S2PERP_terms[1] -= 4.0 * linalg::triplet(D_A, S, D_B)->vector_dot(J_A);
        Exch10_S2PERP_terms[1] += 4.0 * linalg::triplet(D_X, S, D_Y)->vector_dot(J_A);
        Exch10_S2PERP_terms[1] += 2.0 * linalg::triplet(D_A, S, D_B)->vector_dot(K_A);
        Exch10_S2PERP_terms[1] += 2.0 * linalg::triplet(D_X, S, D_B)->vector_dot(K_X);
        Exch10_S2PERP_terms[1] -= 2.0 * linalg::triplet(D_A, S, D_Y)->vector_dot(K_X);
        Exch10_S2PERP_terms[1] -= 2.0 * linalg::triplet(D_X, S, D_Y)->vector_dot(K_A);
       
        Exch10_S2PERP_terms[2] -= 2.0 * linalg::triplet(D_B, S, D_A)->vector_dot(V_B);
        Exch10_S2PERP_terms[2] += 2.0 * linalg::triplet(D_Y, S, D_X)->vector_dot(V_B);
        Exch10_S2PERP_terms[2] -= 4.0 * linalg::triplet(D_B, S, D_A)->vector_dot(J_B);
        Exch10_S2PERP_terms[2] += 4.0 * linalg::triplet(D_Y, S, D_X)->vector_dot(J_B);
        Exch10_S2PERP_terms[2] += 2.0 * linalg::triplet(D_B, S, D_A)->vector_dot(K_B);
        Exch10_S2PERP_terms[2] -= 2.0 * linalg::triplet(D_B, S, D_X)->vector_dot(K_Y);
        Exch10_S2PERP_terms[2] += 2.0 * linalg::triplet(D_Y, S, D_A)->vector_dot(K_Y);
        Exch10_S2PERP_terms[2] -= 2.0 * linalg::triplet(D_Y, S, D_X)->vector_dot(K_B);
       
        Exch10_S2PERP_terms[3] += 2.0 * linalg::triplet(linalg::triplet(D_B, S, D_A), S, D_B)->vector_dot(V_A);
        Exch10_S2PERP_terms[3] -= 2.0 * linalg::triplet(linalg::triplet(D_B, S, D_X), S, D_Y)->vector_dot(V_A);
        Exch10_S2PERP_terms[3] -= 2.0 * linalg::triplet(linalg::triplet(D_Y, S, D_X), S, D_B)->vector_dot(V_A);
        Exch10_S2PERP_terms[3] += 2.0 * linalg::triplet(linalg::triplet(D_Y, S, D_A), S, D_Y)->vector_dot(V_A);
        Exch10_S2PERP_terms[3] += 4.0 * linalg::triplet(linalg::triplet(D_B, S, D_A), S, D_B)->vector_dot(J_A);
        Exch10_S2PERP_terms[3] -= 4.0 * linalg::triplet(linalg::triplet(D_Y, S, D_X), S, D_B)->vector_dot(J_A);
        Exch10_S2PERP_terms[3] -= 4.0 * linalg::triplet(linalg::triplet(D_B, S, D_X), S, D_Y)->vector_dot(J_A);
        Exch10_S2PERP_terms[3] += 4.0 * linalg::triplet(linalg::triplet(D_Y, S, D_A), S, D_Y)->vector_dot(J_A);
       
        Exch10_S2PERP_terms[4] += 2.0 * linalg::triplet(linalg::triplet(D_A, S, D_B), S, D_A)->vector_dot(V_B);
        Exch10_S2PERP_terms[4] -= 2.0 * linalg::triplet(linalg::triplet(D_X, S, D_Y), S, D_A)->vector_dot(V_B);
        Exch10_S2PERP_terms[4] -= 2.0 * linalg::triplet(linalg::triplet(D_A, S, D_Y), S, D_X)->vector_dot(V_B);
        Exch10_S2PERP_terms[4] += 2.0 * linalg::triplet(linalg::triplet(D_X, S, D_B), S, D_X)->vector_dot(V_B);
        Exch10_S2PERP_terms[4] += 4.0 * linalg::triplet(linalg::triplet(D_A, S, D_B), S, D_A)->vector_dot(J_B);
        Exch10_S2PERP_terms[4] -= 4.0 * linalg::triplet(linalg::triplet(D_X, S, D_Y), S, D_A)->vector_dot(J_B);
        Exch10_S2PERP_terms[4] -= 4.0 * linalg::triplet(linalg::triplet(D_A, S, D_Y), S, D_X)->vector_dot(J_B);
        Exch10_S2PERP_terms[4] += 4.0 * linalg::triplet(linalg::triplet(D_X, S, D_B), S, D_X)->vector_dot(J_B);
       
        Exch10_S2PERP_terms[5] -= 2.0 * linalg::triplet(D_A, S, D_B)->vector_dot(K_AOB);
        Exch10_S2PERP_terms[5] -= 2.0 * linalg::triplet(D_X, S, D_B)->vector_dot(K_XOB);
        Exch10_S2PERP_terms[5] += 2.0 * linalg::triplet(D_A, S, D_Y)->vector_dot(K_XOB);
        Exch10_S2PERP_terms[5] += 2.0 * linalg::triplet(D_A, S, D_B)->vector_dot(K_XOY);
        Exch10_S2PERP_terms[5] += 2.0 * linalg::triplet(D_X, S, D_Y)->vector_dot(K_AOB);
        Exch10_S2PERP_terms[5] += 2.0 * linalg::triplet(D_X, S, D_B)->vector_dot(K_AOY);
        Exch10_S2PERP_terms[5] -= 2.0 * linalg::triplet(D_A, S, D_Y)->vector_dot(K_AOY);
        Exch10_S2PERP_terms[5] -= 2.0 * linalg::triplet(D_X, S, D_Y)->vector_dot(K_XOY);
   
        for (int k = 0; k < Exch10_S2PAR_terms.size(); k++) {
            Exch10_S2PAR += Exch10_S2PAR_terms[k];
        }
        for (int k = 0; k < Exch10_S2PAR_terms.size(); k++) {
            outfile->Printf("    Exch10(S2PAR %1d)     = %18.12lf [Eh]\n",k+1,Exch10_S2PAR_terms[k]);
        }
        outfile->Printf("    Exch10(S^2) [PAR]   = %18.12lf [Eh]\n",Exch10_S2PAR);
       
        for (int k = 0; k < Exch10_S2PERP_terms.size(); k++) {
            Exch10_S2PERP += Exch10_S2PERP_terms[k];
        }
        for (int k = 0; k < Exch10_S2PERP_terms.size(); k++) {
            outfile->Printf("    Exch10(S2PERP %1d)    = %18.12lf [Eh]\n",k+1,Exch10_S2PERP_terms[k]);
        }
        outfile->Printf("    Exch10(S^2) [PERP]  = %18.12lf [Eh]\n",Exch10_S2PERP);
        outfile->Printf("    Exch10(S^2) [avg]   = %18.12lf [Eh]\n",0.5*(Exch10_S2PAR+Exch10_S2PERP));
        scalars_["Exch10(S^2)"] = 0.5*(Exch10_S2PAR+Exch10_S2PERP);
    }
    else {
// we only calculate the average of parallel and perpendicular results, so some terms are not needed at all
// this means we compute Eq. (5) but skip all the terms with a +- or -+ in front of them - these terms cancel out
        double Exch10_S2AVG = 0.0;
        std::vector<double> Exch10_S2AVG_terms;
        Exch10_S2AVG_terms.resize(6);
 
        Exch10_S2AVG_terms[0] -= 2.0 * D_A->vector_dot(K_B);
 
        Exch10_S2AVG_terms[1] -= 2.0 * linalg::triplet(D_A, S, D_B)->vector_dot(V_A);
        Exch10_S2AVG_terms[1] -= 4.0 * linalg::triplet(D_A, S, D_B)->vector_dot(J_A);
        Exch10_S2AVG_terms[1] += 2.0 * linalg::triplet(D_A, S, D_B)->vector_dot(K_A);
        Exch10_S2AVG_terms[1] += 2.0 * linalg::triplet(D_X, S, D_B)->vector_dot(K_X);

        Exch10_S2AVG_terms[2] -= 2.0 * linalg::triplet(D_B, S, D_A)->vector_dot(V_B);
        Exch10_S2AVG_terms[2] -= 4.0 * linalg::triplet(D_B, S, D_A)->vector_dot(J_B);
        Exch10_S2AVG_terms[2] += 2.0 * linalg::triplet(D_B, S, D_A)->vector_dot(K_B);
        Exch10_S2AVG_terms[2] += 2.0 * linalg::triplet(D_Y, S, D_A)->vector_dot(K_Y);

        Exch10_S2AVG_terms[3] += 2.0 * linalg::triplet(linalg::triplet(D_B, S, D_A), S, D_B)->vector_dot(V_A);
        Exch10_S2AVG_terms[3] += 2.0 * linalg::triplet(linalg::triplet(D_Y, S, D_A), S, D_Y)->vector_dot(V_A);
        Exch10_S2AVG_terms[3] += 4.0 * linalg::triplet(linalg::triplet(D_B, S, D_A), S, D_B)->vector_dot(J_A);
        Exch10_S2AVG_terms[3] += 4.0 * linalg::triplet(linalg::triplet(D_Y, S, D_A), S, D_Y)->vector_dot(J_A);
     
        Exch10_S2AVG_terms[4] += 2.0 * linalg::triplet(linalg::triplet(D_A, S, D_B), S, D_A)->vector_dot(V_B);
        Exch10_S2AVG_terms[4] += 2.0 * linalg::triplet(linalg::triplet(D_X, S, D_B), S, D_X)->vector_dot(V_B);
        Exch10_S2AVG_terms[4] += 4.0 * linalg::triplet(linalg::triplet(D_A, S, D_B), S, D_A)->vector_dot(J_B);
        Exch10_S2AVG_terms[4] += 4.0 * linalg::triplet(linalg::triplet(D_X, S, D_B), S, D_X)->vector_dot(J_B);
      
        Exch10_S2AVG_terms[5] -= 2.0 * linalg::triplet(D_A, S, D_B)->vector_dot(K_AOB);
        Exch10_S2AVG_terms[5] -= 2.0 * linalg::triplet(D_X, S, D_B)->vector_dot(K_XOB);
        Exch10_S2AVG_terms[5] -= 2.0 * linalg::triplet(D_A, S, D_Y)->vector_dot(K_AOY);
        Exch10_S2AVG_terms[5] -= 2.0 * linalg::triplet(D_X, S, D_Y)->vector_dot(K_XOY);
   
        for (int k = 0; k < Exch10_S2AVG_terms.size(); k++) {
            Exch10_S2AVG += Exch10_S2AVG_terms[k];
        }
        for (int k = 0; k < Exch10_S2AVG_terms.size(); k++) {
            outfile->Printf("    Exch10(S2AVG %1d)     = %18.12lf [Eh]\n",k+1,Exch10_S2AVG_terms[k]);
        }
        outfile->Printf("    Exch10(S^2) [AVG]   = %18.12lf [Eh]\n",Exch10_S2AVG);
        scalars_["Exch10(S^2)"] = Exch10_S2AVG;
    }
// end of parallel- and perpendicular-link-coupling E(10)exch(S^2)

// Now, for the nonapproximated E(10)exch with the SAOn/SIAOn link assignments, there is no shortcut.
// We have to calculate both spin couplings and only average the final number regardless of the 
// FISAPT_EXCH_PARPERP value (we tried averaging out the matrices and the results were very bad).
// First comes the parallel spin coupling between singly occupied link orbitals
    int na = matrices_["AlloccA"]->colspi()[0] - 1;
    int nb = matrices_["AlloccB"]->colspi()[0] - 1;
    std::shared_ptr<Matrix> Saa = linalg::triplet(matrices_["AlloccA"], S, matrices_["AlloccA"], true, false, false);
    std::shared_ptr<Matrix> Sbb = linalg::triplet(matrices_["AlloccB"], S, matrices_["AlloccB"], true, false, false);
// We now construct and invert the extended overlap matrix - Eq. (12)
    auto Tbig = std::make_shared<Matrix>("T", 2*na + 2*nb + 2, 2*na + 2*nb + 2);
    Tbig->identity();
    double** Tbigp = Tbig->pointer();
    for (int a = 0; a < na; a++) {
        for (int b = 0; b < nb; b++) {
            Tbigp[a][b + 2*na + 1] = Tbigp[b + 2*na + 1][a] = Sabp[a][b];
            Tbigp[a + na][b + 2*na + nb + 1] = Tbigp[b + 2*na + nb + 1][a + na] = Sabp[a][b];
        }
// The 0.70710678118655 factor already sits in Sabp
        Tbigp[a][2*na + 2*nb + 1] = Tbigp[2*na + 2*nb + 1][a] = Sabp[a][nb];
        Tbigp[a + na][2*na + 2*nb + 1] = Tbigp[2*na + 2*nb + 1][a + na] = Sabp[a][nb];
    }
    for (int b = 0; b < nb; b++) {
        Tbigp[2*na][b + 2*na + 1] = Tbigp[b + 2*na + 1][2*na] = Sabp[na][b];
        Tbigp[2*na][b + 2*na + nb + 1] = Tbigp[b + 2*na + nb + 1][2*na] = Sabp[na][b];
    }
// Now we need to offset the 1/2 factor in Sabp
    Tbigp[2*na][2*na + 2*nb + 1] = Tbigp[2*na + 2*nb + 1][2*na] = 2 * Sabp[na][nb];
    outfile->Printf("\n Sxy: %12.8f\n",2 * Sabp[na][nb]);
//now let's slice out the relevant pieces of this monster
    Dimension zero(1);
    Dimension nadim(1);
    Dimension naadim(1);
    Dimension naaxdim(1);
    Dimension naaxbdim(1);
    Dimension naaxbbdim(1);
    Dimension naaxbbydim(1);
    zero[0] = 0;
    nadim[0] = na;
    naadim[0] = 2 * na;
    naaxdim[0] = 2 * na + 1;
    naaxbdim[0] = 2 * na + 1 + nb;
    naaxbbdim[0] = 2 * na + 1 + 2*nb;
    naaxbbydim[0] = 2 * na + 2 + 2*nb;
    Tbig->power(-1.0, 1.0E-14);
    Tbigp = Tbig->pointer();
//note that the slices should contain x and y as well
    auto Dii_ss = Tbig->get_block(Slice(nadim, naaxdim), Slice(nadim, naaxdim));
    auto Dij_ss = Tbig->get_block(Slice(nadim, naaxdim), Slice(naaxbdim, naaxbbydim));
    auto Dji_ss = Tbig->get_block(Slice(naaxbdim, naaxbbydim), Slice(nadim, naaxdim));
    auto Djj_ss = Tbig->get_block(Slice(naaxbdim, naaxbbydim), Slice(naaxbdim, naaxbbydim));
//for the opposite-spin case, we have to pick a spin-up slice that is noncontiguous with x/y
    std::vector<SharedMatrix> hold_temp;
    
    hold_temp.push_back(Tbig->get_block(Slice(zero, nadim), Slice(nadim, naaxdim)));
    hold_temp.push_back(Tbig->get_block(Slice(naadim, naaxdim), Slice(nadim, naaxdim)));
    auto Dii_os = linalg::vertcat(hold_temp);
    hold_temp.clear();
    
    hold_temp.push_back(Tbig->get_block(Slice(zero, nadim), Slice(naaxbdim, naaxbbydim)));
    hold_temp.push_back(Tbig->get_block(Slice(naadim, naaxdim), Slice(naaxbdim, naaxbbydim)));
    auto Dij_os = linalg::vertcat(hold_temp);
    hold_temp.clear();
    
    hold_temp.push_back(Tbig->get_block(Slice(naaxdim, naaxbdim), Slice(nadim, naaxdim)));
    hold_temp.push_back(Tbig->get_block(Slice(naaxbbdim, naaxbbydim), Slice(nadim, naaxdim)));
    auto Dji_os = linalg::vertcat(hold_temp);
    hold_temp.clear();

    hold_temp.push_back(Tbig->get_block(Slice(naaxdim, naaxbdim), Slice(naaxbdim, naaxbbydim)));
    hold_temp.push_back(Tbig->get_block(Slice(naaxbbdim, naaxbbydim), Slice(naaxbdim, naaxbbydim)));
    auto Djj_os = linalg::vertcat(hold_temp);
    hold_temp.clear();

//half-transformed matrices for the JK operators
    auto Dsj_ssh = linalg::doublet(matrices_["AlloccA"], Dij_ss, false, false);
    Dsj_ssh->add(linalg::doublet(matrices_["AlloccB"], Djj_ss, false, false));
    auto Dsj_osh = linalg::doublet(matrices_["AlloccA"], Dij_os, false, false);
    Dsj_osh->add(linalg::doublet(matrices_["AlloccB"], Djj_os, false, false));

    Cl.clear();
    Cr.clear();
    Cl.push_back(Dsj_ssh);
    Cr.push_back(matrices_["AlloccB"]);
    Cl.push_back(Dsj_osh);
    Cr.push_back(matrices_["AlloccB"]);
    
    jk_->compute();

    std::shared_ptr<Matrix> J_ss = J[0];
    std::shared_ptr<Matrix> K_ss = K[0];
    std::shared_ptr<Matrix> J_os = J[1];
    std::shared_ptr<Matrix> K_os = K[1];

    K_ss->transpose_this();
    K_os->transpose_this();

//finish the transformation of the D matrices
    auto Dsj_ss = linalg::doublet(Dsj_ssh, matrices_["AlloccB"], false, true);
    auto Dsj_os = linalg::doublet(Dsj_osh, matrices_["AlloccB"], false, true);
    auto Dri_ss = linalg::triplet(matrices_["AlloccA"], Dii_ss, matrices_["AlloccA"], false, false, true);
    Dri_ss->add(linalg::triplet(matrices_["AlloccB"], Dji_ss, matrices_["AlloccA"], false, false, true));
    auto Dri_os = linalg::triplet(matrices_["AlloccA"], Dii_os, matrices_["AlloccA"], false, false, true);
    Dri_os->add(linalg::triplet(matrices_["AlloccB"], Dji_os, matrices_["AlloccA"], false, false, true));
    
    double Enuc = 0.0;
    double** Enuc2p = matrices_["E NUC"]->pointer();
    Enuc += 2.0 * Enuc2p[0][1];  // A - B

//The full first-order energy (elst+exch) is given by Eq. (21)
    double E10_full_ss = Enuc;
    E10_full_ss += 2.0 * Dri_ss->vector_dot(V_B);
    E10_full_ss += 2.0 * Dsj_ss->vector_dot(V_A);
    E10_full_ss += 4.0 * Dri_ss->vector_dot(J_ss);
    E10_full_ss -= 2.0 * Dri_ss->vector_dot(K_ss);
    double E10_full_os = -2.0 * Dri_os->vector_dot(K_os);
    double E10_full_par = E10_full_ss + E10_full_os;
    
    double E10_elst = Enuc;
    E10_elst += 2.0 * D_A->vector_dot(V_B);
    E10_elst += 2.0 * D_B->vector_dot(V_A);
    E10_elst += 4.0 * D_A->vector_dot(J_B);

    outfile->Printf("\n E(10)(full) value (PAR)  : %16.12f ",E10_full_par);
    outfile->Printf("\n SS/OS parts:           %16.12f %16.12f",E10_full_ss,E10_full_os);
    outfile->Printf("\n E(10)(elst) part (PAR)  : %16.12f ",E10_elst);
    outfile->Printf("\n E(10)(exch) part (PAR)  : %16.12f \n\n",E10_full_par-E10_elst);

//Now we do the same but with the perpendicular spin coupling between singly occupied link orbitals
//We now construct and invert the extended overlap matrix - Eq. (13)
    Tbig->identity();
    for (int a = 0; a < na; a++) {
        for (int b = 0; b < nb; b++) {
            Tbigp[a][b + 2*na + 1] = Tbigp[b + 2*na + 1][a] = Sabp[a][b];
            Tbigp[a + na][b + 2*na + nb + 1] = Tbigp[b + 2*na + nb + 1][a + na] = Sabp[a][b];
        }
// The 0.70710678118655 factor already sits in Sabp
        Tbigp[a][2*na + 2*nb + 1] = Tbigp[2*na + 2*nb + 1][a] = -Sabp[a][nb];
        Tbigp[a + na][2*na + 2*nb + 1] = Tbigp[2*na + 2*nb + 1][a + na] = Sabp[a][nb];
// We reversed the signs above (reversed the beta/alpha order) to have contiguous slices below
    }
    for (int b = 0; b < nb; b++) {
        Tbigp[2*na][b + 2*na + 1] = Tbigp[b + 2*na + 1][2*na] = Sabp[na][b];
        Tbigp[2*na][b + 2*na + nb + 1] = Tbigp[b + 2*na + nb + 1][2*na] = Sabp[na][b];
    }
    Tbigp[2*na][2*na + 2*nb + 1] = Tbigp[2*na + 2*nb + 1][2*na] = 0.0;
//now let's slice out the relevant pieces of this monster
    Tbig->power(-1.0, 1.0E-14);
    Tbigp = Tbig->pointer();
//note that the slices should contain x and y as well
    Dii_ss = Tbig->get_block(Slice(nadim, naaxdim), Slice(nadim, naaxdim));
    Dij_ss = Tbig->get_block(Slice(nadim, naaxdim), Slice(naaxbdim, naaxbbydim));
    Dji_ss = Tbig->get_block(Slice(naaxbdim, naaxbbydim), Slice(nadim, naaxdim));
    Djj_ss = Tbig->get_block(Slice(naaxbdim, naaxbbydim), Slice(naaxbdim, naaxbbydim));
//for the opposite-spin case, we have to pick a spin-up slice that is noncontiguous with x/y
    hold_temp.push_back(Tbig->get_block(Slice(zero, nadim), Slice(nadim, naaxdim)));
    hold_temp.push_back(Tbig->get_block(Slice(naadim, naaxdim), Slice(nadim, naaxdim)));
    Dii_os = linalg::vertcat(hold_temp);
    hold_temp.clear();

    hold_temp.push_back(Tbig->get_block(Slice(zero, nadim), Slice(naaxbdim, naaxbbydim)));
    hold_temp.push_back(Tbig->get_block(Slice(naadim, naaxdim), Slice(naaxbdim, naaxbbydim)));
    Dij_os = linalg::vertcat(hold_temp);
    hold_temp.clear();

    hold_temp.push_back(Tbig->get_block(Slice(naaxdim, naaxbdim), Slice(nadim, naaxdim)));
    hold_temp.push_back(Tbig->get_block(Slice(naaxbbdim, naaxbbydim), Slice(nadim, naaxdim)));
    Dji_os = linalg::vertcat(hold_temp);
    hold_temp.clear();

    hold_temp.push_back(Tbig->get_block(Slice(naaxdim, naaxbdim), Slice(naaxbdim, naaxbbydim)));
    hold_temp.push_back(Tbig->get_block(Slice(naaxbbdim, naaxbbydim), Slice(naaxbdim, naaxbbydim)));
    Djj_os = linalg::vertcat(hold_temp);
    hold_temp.clear();

// now we have some signs to adjust in D_os terms containing y
    double** Dij_osp = Dij_os->pointer();
    double** Dji_osp = Dji_os->pointer();
    double** Djj_osp = Djj_os->pointer();
    for (int a = 0; a < na; a++) {
        Dji_osp [nb][a] *= -1;
    }
    for (int b = 0; b < nb; b++) {
        Djj_osp [nb][b] *= -1;  
    }
    Djj_osp [nb][nb] *= -1;

// half-transformed matrices for the JK operators
    Dsj_ssh = linalg::doublet(matrices_["AlloccA"], Dij_ss, false, false);
    Dsj_ssh->add(linalg::doublet(matrices_["AlloccB"], Djj_ss, false, false));
    Dsj_osh = linalg::doublet(matrices_["AlloccA"], Dij_os, false, false);
    Dsj_osh->add(linalg::doublet(matrices_["AlloccB"], Djj_os, false, false));

    Cl.clear();
    Cr.clear();
    Cl.push_back(Dsj_ssh);
    Cr.push_back(matrices_["AlloccB"]);
    Cl.push_back(Dsj_osh);
    Cr.push_back(matrices_["AlloccB"]);
    
    jk_->compute();

    std::shared_ptr<Matrix> J2_ss = J[0];
    std::shared_ptr<Matrix> K2_ss = K[0];
    std::shared_ptr<Matrix> J2_os = J[1];
    std::shared_ptr<Matrix> K2_os = K[1];

    K2_ss->transpose_this();
    K2_os->transpose_this();

// finish the transformation of the D matrices
    Dsj_ss = linalg::doublet(Dsj_ssh, matrices_["AlloccB"], false, true);
    Dsj_os = linalg::doublet(Dsj_osh, matrices_["AlloccB"], false, true);
    Dri_ss = linalg::triplet(matrices_["AlloccA"], Dii_ss, matrices_["AlloccA"], false, false, true);
    Dri_ss->add(linalg::triplet(matrices_["AlloccB"], Dji_ss, matrices_["AlloccA"], false, false, true));
    Dri_os = linalg::triplet(matrices_["AlloccA"], Dii_os, matrices_["AlloccA"], false, false, true);
    Dri_os->add(linalg::triplet(matrices_["AlloccB"], Dji_os, matrices_["AlloccA"], false, false, true));

// The full first-order energy (elst+exch) is given by Eq. (21)
    E10_full_ss = Enuc;
    E10_full_ss += 2.0 * Dri_ss->vector_dot(V_B);
    E10_full_ss += 2.0 * Dsj_ss->vector_dot(V_A);
    E10_full_ss += 4.0 * Dri_ss->vector_dot(J2_ss);
    E10_full_ss -= 2.0 * Dri_ss->vector_dot(K2_ss);
    E10_full_os = -2.0 * Dri_os->vector_dot(K2_os);
    double E10_full_perp = E10_full_ss + E10_full_os;
    
    outfile->Printf("\n E(10)(full) value (PERP) : %16.12f ",E10_full_perp);
    outfile->Printf("\n SS/OS parts:            %16.12f %16.12f",E10_full_ss,E10_full_os);
    outfile->Printf("\n E(10)(exch) part (PERP): : %16.12f \n\n",E10_full_perp-E10_elst);
    scalars_["Exch10"] = 0.5*(E10_full_par+E10_full_perp)-E10_elst;
//  The algorithm which only calculates the average of parallel and perpendicular matrices proved to be
//  inaccurate (and not invariant with respect to the A<->B interchange) so the code has been removed.

//end of the new E(10)exch(full) code
    } // new link assignments

    if (options_.get_bool("SSAPT0_SCALE")) {
        sSAPT0_scale_ = scalars_["Exch10"] / scalars_["Exch10(S^2)"];
        sSAPT0_scale_ = pow(sSAPT0_scale_, 3.0);
        outfile->Printf("    Scaling F-SAPT Exch-Ind and Exch-Disp by %11.3E \n\n", sSAPT0_scale_);
    }

}

// Compute the total induction contribution
void FISAPT::ind() {
    outfile->Printf("  ==> Induction <==\n\n");

    // => Pointers <= //

    std::shared_ptr<Matrix> S = matrices_["S"];

    std::shared_ptr<Matrix> D_A = matrices_["D_A"];
    std::shared_ptr<Matrix> D_B = matrices_["D_B"];
    std::shared_ptr<Matrix> P_A = matrices_["P_A"];
    std::shared_ptr<Matrix> P_B = matrices_["P_B"];
    std::shared_ptr<Matrix> V_A = matrices_["V_A"];
    std::shared_ptr<Matrix> V_B = matrices_["V_B"];
    std::shared_ptr<Matrix> J_A = matrices_["J_A"];
    std::shared_ptr<Matrix> J_B = matrices_["J_B"];
    std::shared_ptr<Matrix> K_A = matrices_["K_A"];
    std::shared_ptr<Matrix> K_B = matrices_["K_B"];

    std::shared_ptr<Matrix> Cocc0A = matrices_["Cocc0A"];
    std::shared_ptr<Matrix> Cocc0B = matrices_["Cocc0B"];
    std::shared_ptr<Matrix> Cvir0A = matrices_["Cvir0A"];
    std::shared_ptr<Matrix> Cvir0B = matrices_["Cvir0B"];

    std::shared_ptr<Vector> eps_occ0A = vectors_["eps_occ0A"];
    std::shared_ptr<Vector> eps_occ0B = vectors_["eps_occ0B"];
    std::shared_ptr<Vector> eps_vir0A = vectors_["eps_vir0A"];
    std::shared_ptr<Vector> eps_vir0B = vectors_["eps_vir0B"];

    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    int na = eps_occ0A->dimpi()[0];
    int nb = eps_occ0B->dimpi()[0];
    int nr = eps_vir0A->dimpi()[0];
    int ns = eps_vir0B->dimpi()[0];

    auto uA = std::make_shared<Matrix>("uA", nb, ns);
    auto wA = std::make_shared<Matrix>("wA", nb, ns);
    auto uApar = std::make_shared<Matrix>("uApar", nb, ns);
    auto uAperp = std::make_shared<Matrix>("uAperp", nb, ns);
    auto uB = std::make_shared<Matrix>("uB", na, nr);
    auto wB = std::make_shared<Matrix>("wB", na, nr);
    auto uBpar = std::make_shared<Matrix>("uBpar", na, nr);
    auto uBperp = std::make_shared<Matrix>("uBperp", na, nr);

// Now come the second-order exchange-induction energy expressions appropriate for the SAOn/SIAOn link assignments.
// The expressions depend on the (parallel, perpendicular, or averaged) spin coupling between the two singly occupied link orbitals.
// If FISAPT_EXCH_PARPERP == true, all 3 spin couplings are calculated.
// If FISAPT_EXCH_PARPERP == false, only the averaged spin coupling is computed, avoiding the computation of a bunch of terms that cancel out in that case.
// See https://doi.org/10.1021/acs.jpca.2c06465 for details - equation numbers below refer to this paper and its SI.
    std::string link_assignment = options_.get_str("FISAPT_LINK_ASSIGNMENT");
    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2" ) {

//here we need the link-only parts of density matrices D_X,D_Y and their corresponding Coulomb and exchange matrices (computed earlier)
        std::shared_ptr<Matrix> D_X(D_A->clone());
        std::shared_ptr<Matrix> D_Y(D_A->clone());
        D_X = linalg::doublet(matrices_["thislinkA"], matrices_["thislinkA"], false, true);
        D_Y = linalg::doublet(matrices_["thislinkB"], matrices_["thislinkB"], false, true);
        auto J_X(matrices_["JLA"]->clone());
        auto K_X(matrices_["KLA"]->clone());
        auto J_Y(matrices_["JLB"]->clone());
        auto K_Y(matrices_["KLB"]->clone());
//we need a couple extra J/K matrices as well
        std::shared_ptr<Matrix> C_AOB = linalg::triplet(D_B, S, matrices_["AlloccA"]);
        std::shared_ptr<Matrix> C_XOB = linalg::triplet(D_B, S, matrices_["thislinkA"]);
        std::shared_ptr<Matrix> C_AOY = linalg::triplet(D_Y, S, matrices_["AlloccA"]);
        std::shared_ptr<Matrix> C_XOY = linalg::triplet(D_Y, S, matrices_["thislinkA"]);
        Cl.clear();
        Cr.clear();
        Cl.push_back(matrices_["AlloccA"]);
        Cl.push_back(matrices_["thislinkA"]);
        Cl.push_back(matrices_["AlloccA"]);
        Cl.push_back(matrices_["thislinkA"]);
        Cr.push_back(C_AOB);
        Cr.push_back(C_XOB);
        Cr.push_back(C_AOY);
        Cr.push_back(C_XOY);
        jk_->compute();
        std::shared_ptr<Matrix> K_AOB = K[0];
        std::shared_ptr<Matrix> K_XOB = K[1];
        std::shared_ptr<Matrix> K_AOY = K[2];
        std::shared_ptr<Matrix> K_XOY = K[3];
        std::shared_ptr<Matrix> J_XOY = J[3];

    // => ExchInd perturbations <= //

        std::shared_ptr<Matrix> C_O_A = linalg::triplet(D_B, S, matrices_["AlloccA"]);
        std::shared_ptr<Matrix> C_P_A = linalg::triplet(linalg::triplet(D_B, S, D_A), S, matrices_["AlloccB"]);
        std::shared_ptr<Matrix> C_P_BXY = linalg::triplet(linalg::triplet(D_B, S, D_X), S, matrices_["thislinkB"]);
        std::shared_ptr<Matrix> C_P_YXB = linalg::triplet(linalg::triplet(D_Y, S, D_X), S, matrices_["AlloccB"]);
        std::shared_ptr<Matrix> C_P_YAY = linalg::triplet(linalg::triplet(D_Y, S, D_A), S, matrices_["thislinkB"]);
        std::shared_ptr<Matrix> C_P_B = linalg::triplet(linalg::triplet(D_A, S, D_B), S, matrices_["AlloccA"]);
        std::shared_ptr<Matrix> C_P_AYX = linalg::triplet(linalg::triplet(D_A, S, D_Y), S, matrices_["thislinkA"]);
        std::shared_ptr<Matrix> C_P_XYA = linalg::triplet(linalg::triplet(D_X, S, D_Y), S, matrices_["AlloccA"]);
        std::shared_ptr<Matrix> C_P_XBX = linalg::triplet(linalg::triplet(D_X, S, D_B), S, matrices_["thislinkA"]);
       
        Cl.clear();
        Cr.clear();

        // J/K[O]
        Cl.push_back(matrices_["AlloccA"]);
        Cr.push_back(C_O_A);
        // J/K[P_B]
        Cl.push_back(matrices_["AlloccA"]);
        Cr.push_back(C_P_B);
        // J/K[P_B]
        Cl.push_back(matrices_["AlloccB"]);
        Cr.push_back(C_P_A);
        Cl.push_back(matrices_["thislinkB"]);
        Cr.push_back(C_P_BXY);
        Cl.push_back(matrices_["AlloccB"]);
        Cr.push_back(C_P_YXB);
        Cl.push_back(matrices_["thislinkB"]);
        Cr.push_back(C_P_YAY);
        Cl.push_back(matrices_["thislinkA"]);
        Cr.push_back(C_P_AYX);
        Cl.push_back(matrices_["AlloccA"]);
        Cr.push_back(C_P_XYA);
        Cl.push_back(matrices_["thislinkA"]);
        Cr.push_back(C_P_XBX);
      
        // => Compute the JK matrices <= //
      
        jk_->compute();
      
        // => Unload the JK Object <= //
      
        std::shared_ptr<Matrix> J_O = J[0];
        std::shared_ptr<Matrix> J_P_B = J[1];
        std::shared_ptr<Matrix> J_P_A = J[2];
        std::shared_ptr<Matrix> J_P_BXY = J[3];
        std::shared_ptr<Matrix> J_P_YXB = J[4];
        std::shared_ptr<Matrix> J_P_YAY = J[5];
        std::shared_ptr<Matrix> J_P_AYX = J[6];
        std::shared_ptr<Matrix> J_P_XYA = J[7];
        std::shared_ptr<Matrix> J_P_XBX = J[8];
      
        std::shared_ptr<Matrix> K_O = K[0];
        std::shared_ptr<Matrix> K_P_B = K[1];
        std::shared_ptr<Matrix> K_P_A = K[2];
       
        // ==> Generalized ESP (Flat and Exchange) <== //
       
        std::map<std::string, std::shared_ptr<Matrix> > mapA;
        mapA["Cocc_A"] = Cocc0A;
        mapA["Cvir_A"] = Cvir0A;
        mapA["S"] = S;
        mapA["D_A"] = D_A;
        mapA["V_A"] = V_A;
        mapA["J_A"] = J_A;
        mapA["K_A"] = K_A;
        mapA["D_B"] = D_B;
        mapA["V_B"] = V_B;
        mapA["J_B"] = J_B;
        mapA["K_B"] = K_B;
        mapA["D_X"] = D_X;
        mapA["J_X"] = J_X;
        mapA["K_X"] = K_X;
        mapA["D_Y"] = D_Y;
        mapA["J_Y"] = J_Y;
        mapA["K_Y"] = K_Y;
        mapA["J_O"] = J_O;
        mapA["K_O"] = K_O;
        mapA["J_XOY"] = J_XOY;
        mapA["K_XOB"] = K_XOB;
        mapA["K_AOY"] = K_AOY;
        mapA["K_XOY"] = K_XOY;
        mapA["J_P"] = J_P_A;
        mapA["J_PBXY"] = J_P_BXY;
        mapA["J_PYXB"] = J_P_YXB;
        mapA["J_PYAY"] = J_P_YAY;

        wB = build_ind_pot(mapA);
        if (options_.get_bool("FISAPT_EXCH_PARPERP")) {
            uBpar = build_exch_ind_pot_par(mapA);
            uBperp = build_exch_ind_pot_perp(mapA);
        }
        else {
            uB = build_exch_ind_pot_avg(mapA);
        }

        K_O->transpose_this();
        K_AOY->transpose_this();
        K_XOB->transpose_this();
        K_XOY->transpose_this();

        std::map<std::string, std::shared_ptr<Matrix> > mapB;
        mapB["Cocc_A"] = Cocc0B;
        mapB["Cvir_A"] = Cvir0B;
        mapB["S"] = S;
        mapB["D_A"] = D_B;
        mapB["V_A"] = V_B;
        mapB["J_A"] = J_B;
        mapB["K_A"] = K_B;
        mapB["D_B"] = D_A;
        mapB["V_B"] = V_A;
        mapB["J_B"] = J_A;
        mapB["K_B"] = K_A;
        mapB["D_X"] = D_Y;
        mapB["J_X"] = J_Y;
        mapB["K_X"] = K_Y;
        mapB["D_Y"] = D_X;
        mapB["J_Y"] = J_X;
        mapB["K_Y"] = K_X;
        mapB["J_O"] = J_O;
        mapB["K_O"] = K_O;
        mapB["J_XOY"] = J_XOY;
        mapB["K_XOB"] = K_AOY;
        mapB["K_AOY"] = K_XOB;
        mapB["K_XOY"] = K_XOY;
        mapB["J_P"] = J_P_B;
        mapB["J_PBXY"] = J_P_AYX;
        mapB["J_PYXB"] = J_P_XYA;
        mapB["J_PYAY"] = J_P_XBX;

        wA = build_ind_pot(mapB);
        if (options_.get_bool("FISAPT_EXCH_PARPERP")) {
            uApar = build_exch_ind_pot_par(mapB);
            uAperp = build_exch_ind_pot_perp(mapB);
        }
        else {
            uA = build_exch_ind_pot_avg(mapB);
        }

        K_O->transpose_this();
        K_AOY->transpose_this();
        K_XOB->transpose_this();
        K_XOY->transpose_this();

        // => Stash for ExchDisp <= //

        matrices_["J_O"] = J_O;
        matrices_["K_O"] = K_O;
        matrices_["J_P_A"] = J_P_A;
        matrices_["J_P_B"] = J_P_B;
        matrices_["J_P_YAY"] = J_P_YAY;
        matrices_["J_P_XBX"] = J_P_XBX;
        matrices_["K_XOB"] = K_XOB;
        matrices_["K_AOY"] = K_AOY;
        matrices_["K_XOY"] = K_XOY;
        matrices_["D_X"] = D_X;
        matrices_["D_Y"] = D_Y;
    }
    else {

    // => ExchInd perturbations <= //

        std::shared_ptr<Matrix> C_O_A = linalg::triplet(D_B, S, matrices_["Cocc_A"]);
        std::shared_ptr<Matrix> C_P_A = linalg::triplet(linalg::triplet(D_B, S, D_A), S, matrices_["Cocc_B"]);
        std::shared_ptr<Matrix> C_P_B = linalg::triplet(linalg::triplet(D_A, S, D_B), S, matrices_["Cocc_A"]);
       
        Cl.clear();
        Cr.clear();

        // J/K[O]
        Cl.push_back(matrices_["Cocc_A"]);
        Cr.push_back(C_O_A);
        // J/K[P_B]
        Cl.push_back(matrices_["Cocc_A"]);
        Cr.push_back(C_P_B);
        // J/K[P_B]
        Cl.push_back(matrices_["Cocc_B"]);
        Cr.push_back(C_P_A);
      
        // => Compute the JK matrices <= //
      
        jk_->compute();
      
        // => Unload the JK Object <= //
      
        std::shared_ptr<Matrix> J_O = J[0];
        std::shared_ptr<Matrix> J_P_B = J[1];
        std::shared_ptr<Matrix> J_P_A = J[2];
      
        std::shared_ptr<Matrix> K_O = K[0];
        std::shared_ptr<Matrix> K_P_B = K[1];
        std::shared_ptr<Matrix> K_P_A = K[2];
       
        // ==> Generalized ESP (Flat and Exchange) <== //
       
        std::map<std::string, std::shared_ptr<Matrix> > mapA;
        mapA["Cocc_A"] = Cocc0A;
        mapA["Cvir_A"] = Cvir0A;
        mapA["S"] = S;
        mapA["D_A"] = D_A;
        mapA["V_A"] = V_A;
        mapA["J_A"] = J_A;
        mapA["K_A"] = K_A;
        mapA["D_B"] = D_B;
        mapA["V_B"] = V_B;
        mapA["J_B"] = J_B;
        mapA["K_B"] = K_B;
        mapA["J_O"] = J_O;
        mapA["K_O"] = K_O;
        mapA["J_P"] = J_P_A;

        wB = build_ind_pot(mapA);
        uB = build_exch_ind_pot(mapA);

        K_O->transpose_this();

        std::map<std::string, std::shared_ptr<Matrix> > mapB;
        mapB["Cocc_A"] = Cocc0B;
        mapB["Cvir_A"] = Cvir0B;
        mapB["S"] = S;
        mapB["D_A"] = D_B;
        mapB["V_A"] = V_B;
        mapB["J_A"] = J_B;
        mapB["K_A"] = K_B;
        mapB["D_B"] = D_A;
        mapB["V_B"] = V_A;
        mapB["J_B"] = J_A;
        mapB["K_B"] = K_A;
        mapB["J_O"] = J_O;
        mapB["K_O"] = K_O;
        mapB["J_P"] = J_P_B;

        wA = build_ind_pot(mapB);
        uA = build_exch_ind_pot(mapB);

        K_O->transpose_this();
    
        // => Stash for ExchDisp <= //

        matrices_["J_O"] = J_O;
        matrices_["K_O"] = K_O;
        matrices_["J_P_A"] = J_P_A;
        matrices_["J_P_B"] = J_P_B;
    }

    // ==> Uncoupled Induction <== //

    std::shared_ptr<Matrix> xuA(wB->clone());
    std::shared_ptr<Matrix> xuB(wA->clone());

    {
        double** xuAp = xuA->pointer();
        double** xuBp = xuB->pointer();
        double** wAp = wA->pointer();
        double** wBp = wB->pointer();
        double* eap = eps_occ0A->pointer();
        double* erp = eps_vir0A->pointer();
        double* ebp = eps_occ0B->pointer();
        double* esp = eps_vir0B->pointer();

        for (int a = 0; a < na; a++) {
            for (int r = 0; r < nr; r++) {
                xuAp[a][r] = wBp[a][r] / (eap[a] - erp[r]);
            }
        }

        for (int b = 0; b < nb; b++) {
            for (int s = 0; s < ns; s++) {
                xuBp[b][s] = wAp[b][s] / (ebp[b] - esp[s]);
            }
        }
    }

    // ==> Induction <== //

    double Ind20u_AB = 2.0 * xuA->vector_dot(wB);
    double Ind20u_BA = 2.0 * xuB->vector_dot(wA);
    double Ind20u = Ind20u_AB + Ind20u_BA;
    scalars_["Ind20,u (A<-B)"] = Ind20u_AB;
    scalars_["Ind20,u (B<-A)"] = Ind20u_BA;
    scalars_["Ind20,u"] = Ind20u;
    outfile->Printf("    Ind20,u (A<-B)      = %18.12lf [Eh]\n", Ind20u_AB);
    outfile->Printf("    Ind20,u (B<-A)      = %18.12lf [Eh]\n", Ind20u_BA);
    outfile->Printf("    Ind20,u             = %18.12lf [Eh]\n", Ind20u);
    // fflush(outfile);

    // => Exchange-Induction <= //
    double ExchInd20u_ABpar;
    double ExchInd20u_BApar;
    double ExchInd20upar;
    double ExchInd20u_ABperp;
    double ExchInd20u_BAperp;
    double ExchInd20uperp;
    double ExchInd20u_AB;
    double ExchInd20u_BA;
    double ExchInd20u;
    double ExchInd20r_ABpar;
    double ExchInd20r_BApar;
    double ExchInd20rpar;
    double ExchInd20r_ABperp;
    double ExchInd20r_BAperp;
    double ExchInd20rperp;
    double ExchInd20r_AB;
    double ExchInd20r_BA;
    double ExchInd20r;

    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2" ) {
        if (options_.get_bool("FISAPT_EXCH_PARPERP")) {
            ExchInd20u_ABpar = 2.0 * xuA->vector_dot(uBpar);
            ExchInd20u_BApar = 2.0 * xuB->vector_dot(uApar);
            ExchInd20upar = ExchInd20u_ABpar + ExchInd20u_BApar;
            outfile->Printf("    Exch-Ind20,u (A<-B) PAR  = %18.12lf [Eh]\n", ExchInd20u_ABpar);
            outfile->Printf("    Exch-Ind20,u (B<-A) PAR  = %18.12lf [Eh]\n", ExchInd20u_BApar);
            outfile->Printf("    Exch-Ind20,u        PAR  = %18.12lf [Eh]\n", ExchInd20upar);
            outfile->Printf("\n");
            ExchInd20u_ABperp = 2.0 * xuA->vector_dot(uBperp);
            ExchInd20u_BAperp = 2.0 * xuB->vector_dot(uAperp);
            ExchInd20uperp = ExchInd20u_ABperp + ExchInd20u_BAperp;
            outfile->Printf("    Exch-Ind20,u (A<-B) PERP = %18.12lf [Eh]\n", ExchInd20u_ABperp);
            outfile->Printf("    Exch-Ind20,u (B<-A) PERP = %18.12lf [Eh]\n", ExchInd20u_BAperp);
            outfile->Printf("    Exch-Ind20,u        PERP = %18.12lf [Eh]\n", ExchInd20uperp);
            outfile->Printf("\n");
            ExchInd20u_AB = 0.5*(ExchInd20u_ABpar+ExchInd20u_ABperp);
            ExchInd20u_BA = 0.5*(ExchInd20u_BApar+ExchInd20u_BAperp);
            ExchInd20u = 0.5*(ExchInd20upar+ExchInd20uperp);
            outfile->Printf("    Exch-Ind20,u        AVG  = %18.12lf [Eh]\n", ExchInd20u);
            outfile->Printf("\n");
        }
        else { 
            ExchInd20u_AB = 2.0 * xuA->vector_dot(uB);
            ExchInd20u_BA = 2.0 * xuB->vector_dot(uA);
            ExchInd20u = ExchInd20u_AB + ExchInd20u_BA;
            outfile->Printf("    Exch-Ind20,u (A<-B) AVG  = %18.12lf [Eh]\n", ExchInd20u_AB);
            outfile->Printf("    Exch-Ind20,u (B<-A) AVG  = %18.12lf [Eh]\n", ExchInd20u_BA);
            outfile->Printf("    Exch-Ind20,u        AVG  = %18.12lf [Eh]\n", ExchInd20u);
            outfile->Printf("\n");
        }
    }
    else { 
        ExchInd20u_AB = 2.0 * xuA->vector_dot(uB);
        ExchInd20u_BA = 2.0 * xuB->vector_dot(uA);
        ExchInd20u = ExchInd20u_AB + ExchInd20u_BA;
        outfile->Printf("    Exch-Ind20,u (A<-B) = %18.12lf [Eh]\n", ExchInd20u_AB);
        outfile->Printf("    Exch-Ind20,u (B<-A) = %18.12lf [Eh]\n", ExchInd20u_BA);
        outfile->Printf("    Exch-Ind20,u        = %18.12lf [Eh]\n", ExchInd20u);
        outfile->Printf("\n");
    }
    if (options_.get_bool("SSAPT0_SCALE")) {
        double scale = sSAPT0_scale_;
        double sExchInd20u_AB = scale * ExchInd20u_AB;
        double sExchInd20u_BA = scale * ExchInd20u_BA;
        double sExchInd20u = sExchInd20u_AB + sExchInd20u_BA;
        outfile->Printf("    sExch-Ind20,u (A<-B) = %18.12lf [Eh]\n", sExchInd20u_AB);
        outfile->Printf("    sExch-Ind20,u (B<-A) = %18.12lf [Eh]\n", sExchInd20u_BA);
        outfile->Printf("    sExch-Ind20,u        = %18.12lf [Eh]\n", sExchInd20u);
        outfile->Printf("\n");
        scalars_["sExch-Ind20,u (A<-B)"] = sExchInd20u_AB;
        scalars_["sExch-Ind20,u (B<-A)"] = sExchInd20u_BA;
        scalars_["sExch-Ind20,u"] = sExchInd20u_AB + sExchInd20u_BA;
    }

    scalars_["Exch-Ind20,u (A<-B)"] = ExchInd20u_AB;
    scalars_["Exch-Ind20,u (B<-A)"] = ExchInd20u_BA;
    scalars_["Exch-Ind20,u"] = ExchInd20u;

    // => Coupled Induction <= //

    auto cphf = std::make_shared<CPHF_FISAPT>();

    // Effective constructor
    cphf->delta_ = options_.get_double("D_CONVERGENCE");
    cphf->maxiter_ = options_.get_int("MAXITER");
    cphf->jk_ = jk_;

    cphf->w_A_ = wB;  // Reversal of convention
    cphf->Cocc_A_ = Cocc0A;
    cphf->Cvir_A_ = Cvir0A;
    cphf->eps_occ_A_ = eps_occ0A;
    cphf->eps_vir_A_ = eps_vir0A;

    cphf->w_B_ = wA;  // Reversal of convention
    cphf->Cocc_B_ = Cocc0B;
    cphf->Cvir_B_ = Cvir0B;
    cphf->eps_occ_B_ = eps_occ0B;
    cphf->eps_vir_B_ = eps_vir0B;

    // Gogo CPKS
    cphf->compute_cphf();

    std::shared_ptr<Matrix> xA = cphf->x_A_;
    std::shared_ptr<Matrix> xB = cphf->x_B_;

    // Backward in Ed's convention
    xA->scale(-1.0);
    xB->scale(-1.0);

    // => Induction <= //

    double Ind20r_AB = 2.0 * xA->vector_dot(wB);
    double Ind20r_BA = 2.0 * xB->vector_dot(wA);
    double Ind20r = Ind20r_AB + Ind20r_BA;
    scalars_["Ind20,r (A<-B)"] = Ind20r_AB;
    scalars_["Ind20,r (B<-A)"] = Ind20r_BA;
    scalars_["Ind20,r"] = Ind20r;
    outfile->Printf("    Ind20,r (A<-B)      = %18.12lf [Eh]\n", Ind20r_AB);
    outfile->Printf("    Ind20,r (B<-A)      = %18.12lf [Eh]\n", Ind20r_BA);
    outfile->Printf("    Ind20,r             = %18.12lf [Eh]\n", Ind20r);
    // fflush(outfile);

    // => Exchange-Induction <= //

    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2" ) {
        if (options_.get_bool("FISAPT_EXCH_PARPERP")) {
            ExchInd20r_ABpar = 2.0 * xA->vector_dot(uBpar);
            ExchInd20r_BApar = 2.0 * xB->vector_dot(uApar);
            ExchInd20rpar = ExchInd20r_ABpar + ExchInd20r_BApar;
            outfile->Printf("    Exch-Ind20,r (A<-B) PAR  = %18.12lf [Eh]\n", ExchInd20r_ABpar);
            outfile->Printf("    Exch-Ind20,r (B<-A) PAR  = %18.12lf [Eh]\n", ExchInd20r_BApar);
            outfile->Printf("    Exch-Ind20,r        PAR  = %18.12lf [Eh]\n", ExchInd20rpar);
            outfile->Printf("\n");
            ExchInd20r_ABperp = 2.0 * xA->vector_dot(uBperp);
            ExchInd20r_BAperp = 2.0 * xB->vector_dot(uAperp);
            ExchInd20rperp = ExchInd20r_ABperp + ExchInd20r_BAperp;
            outfile->Printf("    Exch-Ind20,r (A<-B) PERP = %18.12lf [Eh]\n", ExchInd20r_ABperp);
            outfile->Printf("    Exch-Ind20,r (B<-A) PERP = %18.12lf [Eh]\n", ExchInd20r_BAperp);
            outfile->Printf("    Exch-Ind20,r        PERP = %18.12lf [Eh]\n", ExchInd20rperp);
            outfile->Printf("\n");
            ExchInd20r_AB = 0.5*(ExchInd20r_ABpar+ExchInd20r_ABperp);
            ExchInd20r_BA = 0.5*(ExchInd20r_BApar+ExchInd20r_BAperp);
            ExchInd20r = 0.5*(ExchInd20rpar+ExchInd20rperp);
            outfile->Printf("    Exch-Ind20,r        AVG  = %18.12lf [Eh]\n", ExchInd20r);
            outfile->Printf("\n");
        }
        else { 
            ExchInd20r_AB = 2.0 * xA->vector_dot(uB);
            ExchInd20r_BA = 2.0 * xB->vector_dot(uA);
            ExchInd20r = ExchInd20r_AB + ExchInd20r_BA;
            outfile->Printf("    Exch-Ind20,r (A<-B) AVG  = %18.12lf [Eh]\n", ExchInd20r_AB);
            outfile->Printf("    Exch-Ind20,r (B<-A) AVG  = %18.12lf [Eh]\n", ExchInd20r_BA);
            outfile->Printf("    Exch-Ind20,r        AVG  = %18.12lf [Eh]\n", ExchInd20r);
            outfile->Printf("\n");
        }
    }
    else { 
        ExchInd20r_AB = 2.0 * xA->vector_dot(uB);
        ExchInd20r_BA = 2.0 * xB->vector_dot(uA);
        ExchInd20r = ExchInd20r_AB + ExchInd20r_BA;
        outfile->Printf("    Exch-Ind20,r (A<-B) = %18.12lf [Eh]\n", ExchInd20r_AB);
        outfile->Printf("    Exch-Ind20,r (B<-A) = %18.12lf [Eh]\n", ExchInd20r_BA);
        outfile->Printf("    Exch-Ind20,r        = %18.12lf [Eh]\n", ExchInd20r);
        outfile->Printf("\n");
    }
    if (options_.get_bool("SSAPT0_SCALE")) {
        double scale = sSAPT0_scale_;
        double sExchInd20r_AB = scale * ExchInd20r_AB;
        double sExchInd20r_BA = scale * ExchInd20r_BA;
        double sExchInd20r = sExchInd20r_AB + sExchInd20r_BA;
        outfile->Printf("    sExch-Ind20,r (A<-B) = %18.12lf [Eh]\n", sExchInd20r_AB);
        outfile->Printf("    sExch-Ind20,r (B<-A) = %18.12lf [Eh]\n", sExchInd20r_BA);
        outfile->Printf("    sExch-Ind20,r        = %18.12lf [Eh]\n", sExchInd20r);
        outfile->Printf("\n");
        scalars_["sExch-Ind20,r (A<-B)"] = sExchInd20r_AB;
        scalars_["sExch-Ind20,r (B<-A)"] = sExchInd20r_BA;
        scalars_["sExch-Ind20,r"] = sExchInd20r_AB + sExchInd20r_BA;
    }

    scalars_["Exch-Ind20,r (A<-B)"] = ExchInd20r_AB;
    scalars_["Exch-Ind20,r (B<-A)"] = ExchInd20r_BA;
    scalars_["Exch-Ind20,r"] = ExchInd20r;

    scalars_["delta HF,r (2)"] = 0.0;
    if (scalars_["HF"] != 0.0) {
        scalars_["delta HF,r (2)"] =
            scalars_["HF"] - scalars_["Elst10,r"] - scalars_["Exch10"] - scalars_["Ind20,r"] - scalars_["Exch-Ind20,r"];
    }

    // => Kill the JK Object <= //

    jk_.reset();
}

// build potential for induction contribution
std::shared_ptr<Matrix> FISAPT::build_ind_pot(std::map<std::string, std::shared_ptr<Matrix> >& vars) {
    std::shared_ptr<Matrix> Ca = vars["Cocc_A"];
    std::shared_ptr<Matrix> Cr = vars["Cvir_A"];
    std::shared_ptr<Matrix> V_B = vars["V_B"];
    std::shared_ptr<Matrix> J_B = vars["J_B"];

    std::shared_ptr<Matrix> W(V_B->clone());
    W->copy(J_B);
    W->scale(2.0);
    W->add(V_B);

    return linalg::triplet(Ca, W, Cr, true, false, false);
}

// build exchange-induction potential
std::shared_ptr<Matrix> FISAPT::build_exch_ind_pot(std::map<std::string, std::shared_ptr<Matrix> >& vars) {
    std::shared_ptr<Matrix> Ca = vars["Cocc_A"];
    std::shared_ptr<Matrix> Cr = vars["Cvir_A"];

    std::shared_ptr<Matrix> S = vars["S"];

    std::shared_ptr<Matrix> D_A = vars["D_A"];
    std::shared_ptr<Matrix> J_A = vars["J_A"];
    std::shared_ptr<Matrix> K_A = vars["K_A"];
    std::shared_ptr<Matrix> V_A = vars["V_A"];
    std::shared_ptr<Matrix> D_B = vars["D_B"];
    std::shared_ptr<Matrix> J_B = vars["J_B"];
    std::shared_ptr<Matrix> K_B = vars["K_B"];
    std::shared_ptr<Matrix> V_B = vars["V_B"];

    std::shared_ptr<Matrix> J_O = vars["J_O"];  // J[D^A S D^B]
    std::shared_ptr<Matrix> K_O = vars["K_O"];  // K[D^A S D^B]
    std::shared_ptr<Matrix> J_P = vars["J_P"];  // J[D^B S D^A S D^B]

    std::shared_ptr<Matrix> W(K_B->clone());
    std::shared_ptr<Matrix> T;

    // 1
    W->copy(K_B);
    W->scale(-1.0);

    // 2
    T = linalg::triplet(S, D_B, J_A);
    T->scale(-2.0);
    W->add(T);

    // 3
    T->copy(K_O);
    T->scale(1.0);
    W->add(T);

    // 4
    T->copy(J_O);
    T->scale(-2.0);
    W->add(T);

    // 5
    T = linalg::triplet(S, D_B, K_A);
    T->scale(1.0);
    W->add(T);

    // 6
    T = linalg::triplet(J_B, D_B, S);
    T->scale(-2.0);
    W->add(T);

    // 7
    T = linalg::triplet(K_B, D_B, S);
    T->scale(1.0);
    W->add(T);

    // 8
    T = linalg::triplet(linalg::triplet(S, D_B, J_A), D_B, S);
    T->scale(2.0);
    W->add(T);

    // 9
    T = linalg::triplet(linalg::triplet(J_B, D_A, S), D_B, S);
    T->scale(2.0);
    W->add(T);

    // 10
    T = linalg::triplet(K_O, D_B, S);
    T->scale(-1.0);
    W->add(T);

    // 11
    T->copy(J_P);
    T->scale(2.0);
    W->add(T);

    // 12
    T = linalg::triplet(linalg::triplet(S, D_B, S), D_A, J_B);
    T->scale(2.0);
    W->add(T);

    // 13
    T = linalg::triplet(S, D_B, K_O, false, false, true);
    T->scale(-1.0);
    W->add(T);

    // 14
    T = linalg::triplet(S, D_B, V_A);
    T->scale(-1.0);
    W->add(T);

    // 15
    T = linalg::triplet(V_B, D_B, S);
    T->scale(-1.0);
    W->add(T);

    // 16
    T = linalg::triplet(linalg::triplet(S, D_B, V_A), D_B, S);
    T->scale(1.0);
    W->add(T);

    // 17
    T = linalg::triplet(linalg::triplet(V_B, D_A, S), D_B, S);
    T->scale(1.0);
    W->add(T);

    // 18
    T = linalg::triplet(linalg::triplet(S, D_B, S), D_A, V_B);
    T->scale(1.0);
    W->add(T);

    return linalg::triplet(Ca, W, Cr, true, false, false);
}

// Build the exchange-induction potential - SAOn/SIAOn variants, parallel spin coupling
// This is Eq. (9) in https://doi.org/10.1021/acs.jpca.2c06465 with upper signs picked for all +- and -+
std::shared_ptr<Matrix> FISAPT::build_exch_ind_pot_par(std::map<std::string, std::shared_ptr<Matrix> >& vars) {
    std::shared_ptr<Matrix> Ca = vars["Cocc_A"];
    std::shared_ptr<Matrix> Cr = vars["Cvir_A"];

    std::shared_ptr<Matrix> S = vars["S"];

    std::shared_ptr<Matrix> D_A = vars["D_A"];
    std::shared_ptr<Matrix> J_A = vars["J_A"];
    std::shared_ptr<Matrix> K_A = vars["K_A"];
    std::shared_ptr<Matrix> V_A = vars["V_A"];
    std::shared_ptr<Matrix> D_B = vars["D_B"];
    std::shared_ptr<Matrix> J_B = vars["J_B"];
    std::shared_ptr<Matrix> K_B = vars["K_B"];
    std::shared_ptr<Matrix> V_B = vars["V_B"];
    std::shared_ptr<Matrix> D_X = vars["D_X"];
    std::shared_ptr<Matrix> J_X = vars["J_X"];
    std::shared_ptr<Matrix> K_X = vars["K_X"];
    std::shared_ptr<Matrix> D_Y = vars["D_Y"];
    std::shared_ptr<Matrix> J_Y = vars["J_Y"];
    std::shared_ptr<Matrix> K_Y = vars["K_Y"];

    std::shared_ptr<Matrix> J_O = vars["J_O"];  // J[D^A S D^B]
    std::shared_ptr<Matrix> K_O = vars["K_O"];  // K[D^A S D^B]
    std::shared_ptr<Matrix> J_XOB = vars["J_XOB"];  // J[D^X S D^B]
    std::shared_ptr<Matrix> J_AOY = vars["J_AOY"];  // J[D^A S D^Y]
    std::shared_ptr<Matrix> J_XOY = vars["J_XOY"];  // J[D^X S D^Y]
    std::shared_ptr<Matrix> K_XOB = vars["K_XOB"];  // K[D^X S D^B]
    std::shared_ptr<Matrix> K_AOY = vars["K_AOY"];  // K[D^A S D^Y]
    std::shared_ptr<Matrix> K_XOY = vars["K_XOY"];  // K[D^X S D^Y]

    std::shared_ptr<Matrix> J_P = vars["J_P"];  // J[D^B S D^A S D^B]
    std::shared_ptr<Matrix> J_PYXB = vars["J_PYXB"];  // J[D^Y S D^X S D^B]
    std::shared_ptr<Matrix> J_PBXY = vars["J_PBXY"];  // J[D^B S D^X S D^Y]
    std::shared_ptr<Matrix> J_PYAY = vars["J_PYAY"];  // J[D^Y S D^A S D^Y]

    std::shared_ptr<Matrix> W(K_B->clone());
    std::shared_ptr<Matrix> T;

    // 1
    W->copy(K_B);
    W->scale(-1.0);

    // 2
    T = linalg::triplet(S, D_B, J_A);
    T->scale(-2.0);
    W->add(T);

    // 3
    T->copy(K_O);
    T->scale(1.0);
    W->add(T);
    T->copy(K_XOY);
    T->scale(1.0);
    W->add(T);

    // 4
    T->copy(J_O);
    T->scale(-2.0);
    W->add(T);
    T->copy(J_XOY);
    T->scale(-2.0);
    W->add(T);

    // 5
    T = linalg::triplet(S, D_B, K_A);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(S, D_Y, K_X);
    T->scale(1.0);
    W->add(T);

    // 6
    T = linalg::triplet(J_B, D_B, S);
    T->scale(-2.0);
    W->add(T);

    // 7
    T = linalg::triplet(K_B, D_B, S);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(K_Y, D_Y, S);
    T->scale(1.0);
    W->add(T);

    // 8
    T = linalg::triplet(linalg::triplet(S, D_B, J_A), D_B, S);
    T->scale(2.0);
    W->add(T);
    T = linalg::triplet(linalg::triplet(S, D_Y, J_A), D_Y, S);
    T->scale(2.0);
    W->add(T);

    // 9
    T = linalg::triplet(linalg::triplet(J_B, D_A, S), D_B, S);
    T->scale(2.0);
    W->add(T);
    T = linalg::triplet(linalg::triplet(J_B, D_X, S), D_Y, S);
    T->scale(2.0);
    W->add(T);

    // 10
    T = linalg::triplet(K_O, D_B, S);
    T->scale(-1.0);
    W->add(T);
    T = linalg::triplet(K_XOB, D_Y, S);
    T->scale(-1.0);
    W->add(T);
    T = linalg::triplet(K_XOY, D_B, S);
    T->scale(-1.0);
    W->add(T);
    T = linalg::triplet(K_AOY, D_Y, S);
    T->scale(-1.0);
    W->add(T);

    // 11
    T->copy(J_P);
    T->scale(2.0);
    W->add(T);
    T->copy(J_PYXB);
    T->scale(2.0);
    W->add(T);
    T->copy(J_PBXY);
    T->scale(2.0);
    W->add(T);
    T->copy(J_PYAY);
    T->scale(2.0);
    W->add(T);

    // 12
    T = linalg::triplet(linalg::triplet(S, D_B, S), D_A, J_B);
    T->scale(2.0);
    W->add(T);
    T = linalg::triplet(linalg::triplet(S, D_Y, S), D_X, J_B);
    T->scale(2.0);
    W->add(T);

    // 13
    T = linalg::triplet(S, D_B, K_O, false, false, true);
    T->scale(-1.0);
    W->add(T);
    T = linalg::triplet(S, D_B, K_XOY, false, false, true);
    T->scale(-1.0);
    W->add(T);
    T = linalg::triplet(S, D_Y, K_XOB, false, false, true);
    T->scale(-1.0);
    W->add(T);
    T = linalg::triplet(S, D_Y, K_AOY, false, false, true);
    T->scale(-1.0);
    W->add(T);

    // 14
    T = linalg::triplet(S, D_B, V_A);
    T->scale(-1.0);
    W->add(T);

    // 15
    T = linalg::triplet(V_B, D_B, S);
    T->scale(-1.0);
    W->add(T);

    // 16
    T = linalg::triplet(linalg::triplet(S, D_B, V_A), D_B, S);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(linalg::triplet(S, D_Y, V_A), D_Y, S);
    T->scale(1.0);
    W->add(T);

    // 17
    T = linalg::triplet(linalg::triplet(V_B, D_A, S), D_B, S);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(linalg::triplet(V_B, D_X, S), D_Y, S);
    T->scale(1.0);
    W->add(T);

    // 18
    T = linalg::triplet(linalg::triplet(S, D_B, S), D_A, V_B);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(linalg::triplet(S, D_Y, S), D_X, V_B);
    T->scale(1.0);
    W->add(T);

    return linalg::triplet(Ca, W, Cr, true, false, false);
}

// Build the exchange-induction potential - SAOn/SIAOn variants, perpendicular spin coupling
// This is Eq. (9) in https://doi.org/10.1021/acs.jpca.2c06465 with lower signs picked for all +- and -+
std::shared_ptr<Matrix> FISAPT::build_exch_ind_pot_perp(std::map<std::string, std::shared_ptr<Matrix> >& vars) {
    std::shared_ptr<Matrix> Ca = vars["Cocc_A"];
    std::shared_ptr<Matrix> Cr = vars["Cvir_A"];

    std::shared_ptr<Matrix> S = vars["S"];

    std::shared_ptr<Matrix> D_A = vars["D_A"];
    std::shared_ptr<Matrix> J_A = vars["J_A"];
    std::shared_ptr<Matrix> K_A = vars["K_A"];
    std::shared_ptr<Matrix> V_A = vars["V_A"];
    std::shared_ptr<Matrix> D_B = vars["D_B"];
    std::shared_ptr<Matrix> J_B = vars["J_B"];
    std::shared_ptr<Matrix> K_B = vars["K_B"];
    std::shared_ptr<Matrix> V_B = vars["V_B"];
    std::shared_ptr<Matrix> D_X = vars["D_X"];
    std::shared_ptr<Matrix> J_X = vars["J_X"];
    std::shared_ptr<Matrix> K_X = vars["K_X"];
    std::shared_ptr<Matrix> D_Y = vars["D_Y"];
    std::shared_ptr<Matrix> J_Y = vars["J_Y"];
    std::shared_ptr<Matrix> K_Y = vars["K_Y"];

    std::shared_ptr<Matrix> J_O = vars["J_O"];  // J[D^A S D^B]
    std::shared_ptr<Matrix> K_O = vars["K_O"];  // K[D^A S D^B]
    std::shared_ptr<Matrix> J_XOB = vars["J_XOB"];  // J[D^X S D^B]
    std::shared_ptr<Matrix> J_AOY = vars["J_AOY"];  // J[D^A S D^Y]
    std::shared_ptr<Matrix> J_XOY = vars["J_XOY"];  // J[D^X S D^Y]
    std::shared_ptr<Matrix> K_XOB = vars["K_XOB"];  // K[D^X S D^B]
    std::shared_ptr<Matrix> K_AOY = vars["K_AOY"];  // K[D^A S D^Y]
    std::shared_ptr<Matrix> K_XOY = vars["K_XOY"];  // K[D^X S D^Y]

    std::shared_ptr<Matrix> J_P = vars["J_P"];  // J[D^B S D^A S D^B]
    std::shared_ptr<Matrix> J_PYXB = vars["J_PYXB"];  // J[D^Y S D^X S D^B]
    std::shared_ptr<Matrix> J_PBXY = vars["J_PBXY"];  // J[D^B S D^X S D^Y]
    std::shared_ptr<Matrix> J_PYAY = vars["J_PYAY"];  // J[D^Y S D^A S D^Y]

    std::shared_ptr<Matrix> W(K_B->clone());
    std::shared_ptr<Matrix> T;

    // 1
    W->copy(K_B);
    W->scale(-1.0);

    // 2
    T = linalg::triplet(S, D_B, J_A);
    T->scale(-2.0);
    W->add(T);

    // 3
    T->copy(K_O);
    T->scale(1.0);
    W->add(T);
    T->copy(K_XOY);
    T->scale(-1.0);
    W->add(T);

    // 4
    T->copy(J_O);
    T->scale(-2.0);
    W->add(T);
    T->copy(J_XOY);
    T->scale(2.0);
    W->add(T);

    // 5
    T = linalg::triplet(S, D_B, K_A);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(S, D_Y, K_X);
    T->scale(-1.0);
    W->add(T);

    // 6
    T = linalg::triplet(J_B, D_B, S);
    T->scale(-2.0);
    W->add(T);

    // 7
    T = linalg::triplet(K_B, D_B, S);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(K_Y, D_Y, S);
    T->scale(1.0);
    W->add(T);

    // 8
    T = linalg::triplet(linalg::triplet(S, D_B, J_A), D_B, S);
    T->scale(2.0);
    W->add(T);
    T = linalg::triplet(linalg::triplet(S, D_Y, J_A), D_Y, S);
    T->scale(2.0);
    W->add(T);

    // 9
    T = linalg::triplet(linalg::triplet(J_B, D_A, S), D_B, S);
    T->scale(2.0);
    W->add(T);
    T = linalg::triplet(linalg::triplet(J_B, D_X, S), D_Y, S);
    T->scale(-2.0);
    W->add(T);

    // 10
    T = linalg::triplet(K_O, D_B, S);
    T->scale(-1.0);
    W->add(T);
    T = linalg::triplet(K_XOB, D_Y, S);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(K_XOY, D_B, S);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(K_AOY, D_Y, S);
    T->scale(-1.0);
    W->add(T);

    // 11
    T->copy(J_P);
    T->scale(2.0);
    W->add(T);
    T->copy(J_PYXB);
    T->scale(-2.0);
    W->add(T);
    T->copy(J_PBXY);
    T->scale(-2.0);
    W->add(T);
    T->copy(J_PYAY);
    T->scale(2.0);
    W->add(T);

    // 12
    T = linalg::triplet(linalg::triplet(S, D_B, S), D_A, J_B);
    T->scale(2.0);
    W->add(T);
    T = linalg::triplet(linalg::triplet(S, D_Y, S), D_X, J_B);
    T->scale(-2.0);
    W->add(T);

    // 13
    T = linalg::triplet(S, D_B, K_O, false, false, true);
    T->scale(-1.0);
    W->add(T);
    T = linalg::triplet(S, D_B, K_XOY, false, false, true);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(S, D_Y, K_XOB, false, false, true);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(S, D_Y, K_AOY, false, false, true);
    T->scale(-1.0);
    W->add(T);

    // 14
    T = linalg::triplet(S, D_B, V_A);
    T->scale(-1.0);
    W->add(T);

    // 15
    T = linalg::triplet(V_B, D_B, S);
    T->scale(-1.0);
    W->add(T);

    // 16
    T = linalg::triplet(linalg::triplet(S, D_B, V_A), D_B, S);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(linalg::triplet(S, D_Y, V_A), D_Y, S);
    T->scale(1.0);
    W->add(T);

    // 17
    T = linalg::triplet(linalg::triplet(V_B, D_A, S), D_B, S);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(linalg::triplet(V_B, D_X, S), D_Y, S);
    T->scale(-1.0);
    W->add(T);

    // 18
    T = linalg::triplet(linalg::triplet(S, D_B, S), D_A, V_B);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(linalg::triplet(S, D_Y, S), D_X, V_B);
    T->scale(-1.0);
    W->add(T);

    return linalg::triplet(Ca, W, Cr, true, false, false);
}

// Build the exchange-induction potential - SAOn/SIAOn variants, averaged spin coupling
// This is Eq. (9) in https://doi.org/10.1021/acs.jpca.2c06465 with all the +- and -+ omitted as they cancel out when averaging
std::shared_ptr<Matrix> FISAPT::build_exch_ind_pot_avg(std::map<std::string, std::shared_ptr<Matrix> >& vars) {
    std::shared_ptr<Matrix> Ca = vars["Cocc_A"];
    std::shared_ptr<Matrix> Cr = vars["Cvir_A"];

    std::shared_ptr<Matrix> S = vars["S"];

    std::shared_ptr<Matrix> D_A = vars["D_A"];
    std::shared_ptr<Matrix> J_A = vars["J_A"];
    std::shared_ptr<Matrix> K_A = vars["K_A"];
    std::shared_ptr<Matrix> V_A = vars["V_A"];
    std::shared_ptr<Matrix> D_B = vars["D_B"];
    std::shared_ptr<Matrix> J_B = vars["J_B"];
    std::shared_ptr<Matrix> K_B = vars["K_B"];
    std::shared_ptr<Matrix> V_B = vars["V_B"];
    std::shared_ptr<Matrix> D_X = vars["D_X"];
    std::shared_ptr<Matrix> J_X = vars["J_X"];
    std::shared_ptr<Matrix> K_X = vars["K_X"];
    std::shared_ptr<Matrix> D_Y = vars["D_Y"];
    std::shared_ptr<Matrix> J_Y = vars["J_Y"];
    std::shared_ptr<Matrix> K_Y = vars["K_Y"];

    std::shared_ptr<Matrix> J_O = vars["J_O"];  // J[D^A S D^B]
    std::shared_ptr<Matrix> K_O = vars["K_O"];  // K[D^A S D^B]
    std::shared_ptr<Matrix> K_AOY = vars["K_AOY"];  // K[D^A S D^Y]

    std::shared_ptr<Matrix> J_P = vars["J_P"];  // J[D^B S D^A S D^B]
    std::shared_ptr<Matrix> J_PYAY = vars["J_PYAY"];  // J[D^Y S D^A S D^Y]

    std::shared_ptr<Matrix> W(K_B->clone());
    std::shared_ptr<Matrix> T;

    // 1
    W->copy(K_B);
    W->scale(-1.0);

    // 2
    T = linalg::triplet(S, D_B, J_A);
    T->scale(-2.0);
    W->add(T);

    // 3
    T->copy(K_O);
    T->scale(1.0);
    W->add(T);

    // 4
    T->copy(J_O);
    T->scale(-2.0);
    W->add(T);

    // 5
    T = linalg::triplet(S, D_B, K_A);
    T->scale(1.0);
    W->add(T);

    // 6
    T = linalg::triplet(J_B, D_B, S);
    T->scale(-2.0);
    W->add(T);

    // 7
    T = linalg::triplet(K_B, D_B, S);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(K_Y, D_Y, S);
    T->scale(1.0);
    W->add(T);

    // 8
    T = linalg::triplet(linalg::triplet(S, D_B, J_A), D_B, S);
    T->scale(2.0);
    W->add(T);
    T = linalg::triplet(linalg::triplet(S, D_Y, J_A), D_Y, S);
    T->scale(2.0);
    W->add(T);

    // 9
    T = linalg::triplet(linalg::triplet(J_B, D_A, S), D_B, S);
    T->scale(2.0);
    W->add(T);

    // 10
    T = linalg::triplet(K_O, D_B, S);
    T->scale(-1.0);
    W->add(T);
    T = linalg::triplet(K_AOY, D_Y, S);
    T->scale(-1.0);
    W->add(T);

    // 11
    T->copy(J_P);
    T->scale(2.0);
    W->add(T);
    T->copy(J_PYAY);
    T->scale(2.0);
    W->add(T);

    // 12
    T = linalg::triplet(linalg::triplet(S, D_B, S), D_A, J_B);
    T->scale(2.0);
    W->add(T);

    // 13
    T = linalg::triplet(S, D_B, K_O, false, false, true);
    T->scale(-1.0);
    W->add(T);
    T = linalg::triplet(S, D_Y, K_AOY, false, false, true);
    T->scale(-1.0);
    W->add(T);

    // 14
    T = linalg::triplet(S, D_B, V_A);
    T->scale(-1.0);
    W->add(T);

    // 15
    T = linalg::triplet(V_B, D_B, S);
    T->scale(-1.0);
    W->add(T);

    // 16
    T = linalg::triplet(linalg::triplet(S, D_B, V_A), D_B, S);
    T->scale(1.0);
    W->add(T);
    T = linalg::triplet(linalg::triplet(S, D_Y, V_A), D_Y, S);
    T->scale(1.0);
    W->add(T);

    // 17
    T = linalg::triplet(linalg::triplet(V_B, D_A, S), D_B, S);
    T->scale(1.0);
    W->add(T);

    // 18
    T = linalg::triplet(linalg::triplet(S, D_B, S), D_A, V_B);
    T->scale(1.0);
    W->add(T);

    return linalg::triplet(Ca, W, Cr, true, false, false);
}

// Compute total dispersion contribution
// This definition is NOT appropriate for the SAOn/SIAOn ISAPT variants, but it is not called in an ISAPT workflow (FISAPT::fdisp is).
void FISAPT::disp(std::map<std::string, SharedMatrix> matrix_cache, std::map<std::string, SharedVector> vector_cache,
                  bool do_print) {
    if (do_print) {
        outfile->Printf("  ==> Dispersion <==\n\n");
    }

    // => Auxiliary Basis Set <= //
    std::shared_ptr<BasisSet> auxiliary = reference_->get_basisset("DF_BASIS_SAPT");

    // => Pointers <= //

    std::shared_ptr<Matrix> Cocc0A = matrix_cache["Caocc0A"];
    std::shared_ptr<Matrix> Cocc0B = matrix_cache["Caocc0B"];
    std::shared_ptr<Matrix> Cvir0A = matrix_cache["Cvir0A"];
    std::shared_ptr<Matrix> Cvir0B = matrix_cache["Cvir0B"];

    std::shared_ptr<Vector> eps_occ0A = vector_cache["eps_aocc0A"];
    std::shared_ptr<Vector> eps_occ0B = vector_cache["eps_aocc0B"];
    std::shared_ptr<Vector> eps_vir0A = vector_cache["eps_vir0A"];
    std::shared_ptr<Vector> eps_vir0B = vector_cache["eps_vir0B"];

    //For the SAOn/SIAOn variants, we are NOT calling freeze_core() to regenerate the matrices Caocc0A, Caocc0B after
    //redoing the orbitals, so we should not be loading old orbital coefficients.
    std::string link_assignment = options_.get_str("FISAPT_LINK_ASSIGNMENT");
    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2" ) {
        Cocc0A = matrix_cache["Cocc0A"];
        Cocc0B = matrix_cache["Cocc0B"];
        eps_occ0A = vector_cache["eps_occ0A"];
        eps_occ0B = vector_cache["eps_occ0B"];
    }

    // => Sizing <= //

    int nn = primary_->nbf();

    int na = Cocc0A->colspi()[0];
    int nb = Cocc0B->colspi()[0];
    int nr = Cvir0A->colspi()[0];
    int ns = Cvir0B->colspi()[0];
    int nQ = auxiliary->nbf();
    size_t nrQ = nr * (size_t)nQ;
    size_t nsQ = ns * (size_t)nQ;

    int nT = 1;
#ifdef _OPENMP
    nT = Process::environment.get_n_threads();
#endif

    // => Stashed Variables <= //

    std::shared_ptr<Matrix> S = matrix_cache["S"];
    std::shared_ptr<Matrix> D_A = matrix_cache["D_A"];
    std::shared_ptr<Matrix> P_A = matrix_cache["P_A"];
    std::shared_ptr<Matrix> V_A = matrix_cache["V_A"];
    std::shared_ptr<Matrix> J_A = matrix_cache["J_A"];
    std::shared_ptr<Matrix> K_A = matrix_cache["K_A"];
    std::shared_ptr<Matrix> D_B = matrix_cache["D_B"];
    std::shared_ptr<Matrix> P_B = matrix_cache["P_B"];
    std::shared_ptr<Matrix> V_B = matrix_cache["V_B"];
    std::shared_ptr<Matrix> J_B = matrix_cache["J_B"];
    std::shared_ptr<Matrix> K_B = matrix_cache["K_B"];
    std::shared_ptr<Matrix> K_O = matrix_cache["K_O"];

    // => Auxiliary C matrices <= //

    std::shared_ptr<Matrix> Cr1 = linalg::triplet(D_B, S, Cvir0A);
    Cr1->scale(-1.0);
    Cr1->add(Cvir0A);
    std::shared_ptr<Matrix> Cs1 = linalg::triplet(D_A, S, Cvir0B);
    Cs1->scale(-1.0);
    Cs1->add(Cvir0B);
    std::shared_ptr<Matrix> Ca2 = linalg::triplet(D_B, S, Cocc0A);
    std::shared_ptr<Matrix> Cb2 = linalg::triplet(D_A, S, Cocc0B);
    std::shared_ptr<Matrix> Cr3 = linalg::triplet(D_B, S, Cvir0A);
    std::shared_ptr<Matrix> CrX = linalg::triplet(linalg::triplet(D_A, S, D_B), S, Cvir0A);
    Cr3->subtract(CrX);
    Cr3->scale(2.0);
    std::shared_ptr<Matrix> Cs3 = linalg::triplet(D_A, S, Cvir0B);
    std::shared_ptr<Matrix> CsX = linalg::triplet(linalg::triplet(D_B, S, D_A), S, Cvir0B);
    Cs3->subtract(CsX);
    Cs3->scale(2.0);
    std::shared_ptr<Matrix> Ca4 = linalg::triplet(linalg::triplet(D_A, S, D_B), S, Cocc0A);
    Ca4->scale(-2.0);
    std::shared_ptr<Matrix> Cb4 = linalg::triplet(linalg::triplet(D_B, S, D_A), S, Cocc0B);
    Cb4->scale(-2.0);

    // => Auxiliary V matrices <= //

    std::shared_ptr<Matrix> Jbr = linalg::triplet(Cocc0B, J_A, Cvir0A, true, false, false);
    Jbr->scale(2.0);
    std::shared_ptr<Matrix> Kbr = linalg::triplet(Cocc0B, K_A, Cvir0A, true, false, false);
    Kbr->scale(-1.0);

    std::shared_ptr<Matrix> Jas = linalg::triplet(Cocc0A, J_B, Cvir0B, true, false, false);
    Jas->scale(2.0);
    std::shared_ptr<Matrix> Kas = linalg::triplet(Cocc0A, K_B, Cvir0B, true, false, false);
    Kas->scale(-1.0);

    std::shared_ptr<Matrix> KOas = linalg::triplet(Cocc0A, K_O, Cvir0B, true, false, false);
    KOas->scale(1.0);
    std::shared_ptr<Matrix> KObr = linalg::triplet(Cocc0B, K_O, Cvir0A, true, true, false);
    KObr->scale(1.0);

    std::shared_ptr<Matrix> JBas = linalg::triplet(linalg::triplet(Cocc0A, S, D_B, true, false, false), J_A, Cvir0B);
    JBas->scale(-2.0);
    std::shared_ptr<Matrix> JAbr = linalg::triplet(linalg::triplet(Cocc0B, S, D_A, true, false, false), J_B, Cvir0A);
    JAbr->scale(-2.0);

    std::shared_ptr<Matrix> Jbs = linalg::triplet(Cocc0B, J_A, Cvir0B, true, false, false);
    Jbs->scale(4.0);
    std::shared_ptr<Matrix> Jar = linalg::triplet(Cocc0A, J_B, Cvir0A, true, false, false);
    Jar->scale(4.0);

    std::shared_ptr<Matrix> JAas = linalg::triplet(linalg::triplet(Cocc0A, J_B, D_A, true, false, false), S, Cvir0B);
    JAas->scale(-2.0);
    std::shared_ptr<Matrix> JBbr = linalg::triplet(linalg::triplet(Cocc0B, J_A, D_B, true, false, false), S, Cvir0A);
    JBbr->scale(-2.0);

    // Get your signs right Hesselmann!
    std::shared_ptr<Matrix> Vbs = linalg::triplet(Cocc0B, V_A, Cvir0B, true, false, false);
    Vbs->scale(2.0);
    std::shared_ptr<Matrix> Var = linalg::triplet(Cocc0A, V_B, Cvir0A, true, false, false);
    Var->scale(2.0);
    std::shared_ptr<Matrix> VBas = linalg::triplet(linalg::triplet(Cocc0A, S, D_B, true, false, false), V_A, Cvir0B);
    VBas->scale(-1.0);
    std::shared_ptr<Matrix> VAbr = linalg::triplet(linalg::triplet(Cocc0B, S, D_A, true, false, false), V_B, Cvir0A);
    VAbr->scale(-1.0);
    std::shared_ptr<Matrix> VRas = linalg::triplet(linalg::triplet(Cocc0A, V_B, P_A, true, false, false), S, Cvir0B);
    VRas->scale(1.0);
    std::shared_ptr<Matrix> VSbr = linalg::triplet(linalg::triplet(Cocc0B, V_A, P_B, true, false, false), S, Cvir0A);
    VSbr->scale(1.0);

    std::shared_ptr<Matrix> Sas = linalg::triplet(Cocc0A, S, Cvir0B, true, false, false);
    std::shared_ptr<Matrix> Sbr = linalg::triplet(Cocc0B, S, Cvir0A, true, false, false);

    std::shared_ptr<Matrix> Qbr(Jbr->clone());
    Qbr->zero();
    Qbr->add(Jbr);
    Qbr->add(Kbr);
    Qbr->add(KObr);
    Qbr->add(JAbr);
    Qbr->add(JBbr);
    Qbr->add(VAbr);
    Qbr->add(VSbr);

    std::shared_ptr<Matrix> Qas(Jas->clone());
    Qas->zero();
    Qas->add(Jas);
    Qas->add(Kas);
    Qas->add(KOas);
    Qas->add(JAas);
    Qas->add(JBas);
    Qas->add(VBas);
    Qas->add(VRas);

    std::shared_ptr<Matrix> SBar = linalg::triplet(linalg::triplet(Cocc0A, S, D_B, true, false, false), S, Cvir0A);
    std::shared_ptr<Matrix> SAbs = linalg::triplet(linalg::triplet(Cocc0B, S, D_A, true, false, false), S, Cvir0B);

    std::shared_ptr<Matrix> Qar(Jar->clone());
    Qar->zero();
    Qar->add(Jar);
    Qar->add(Var);

    std::shared_ptr<Matrix> Qbs(Jbs->clone());
    Qbs->zero();
    Qbs->add(Jbs);
    Qbs->add(Vbs);

    Jbr.reset();
    Kbr.reset();
    Jas.reset();
    Kas.reset();
    KOas.reset();
    KObr.reset();
    JBas.reset();
    JAbr.reset();
    Jbs.reset();
    Jar.reset();
    JAas.reset();
    JBbr.reset();
    Vbs.reset();
    Var.reset();
    VBas.reset();
    VAbr.reset();
    VRas.reset();
    VSbr.reset();

    // => Memory <= //

    std::vector<std::shared_ptr<Matrix> > Cs;
    Cs.push_back(Cocc0A);
    Cs.push_back(Cvir0A);
    Cs.push_back(Cocc0B);
    Cs.push_back(Cvir0B);
    Cs.push_back(Cr1);
    Cs.push_back(Cs1);
    Cs.push_back(Ca2);
    Cs.push_back(Cb2);
    Cs.push_back(Cr3);
    Cs.push_back(Cs3);
    Cs.push_back(Ca4);
    Cs.push_back(Cb4);

    size_t max_MO = 0, ncol = 0;
    for (auto& mat : Cs) {
        max_MO = std::max(max_MO, (size_t)mat->ncol());
        ncol += (size_t)mat->ncol();
    }

    // => Get integrals from DFHelper <= //
    auto dfh(std::make_shared<DFHelper>(primary_, auxiliary));
    dfh->set_memory(doubles_ - Cs[0]->nrow() * ncol);
    dfh->set_method("DIRECT_iaQ");
    dfh->set_nthreads(nT);
    dfh->initialize();
    dfh->print_header();

    dfh->add_space("a", Cs[0]);
    dfh->add_space("r", Cs[1]);
    dfh->add_space("b", Cs[2]);
    dfh->add_space("s", Cs[3]);
    dfh->add_space("r1", Cs[4]);
    dfh->add_space("s1", Cs[5]);
    dfh->add_space("a2", Cs[6]);
    dfh->add_space("b2", Cs[7]);
    dfh->add_space("r3", Cs[8]);
    dfh->add_space("s3", Cs[9]);
    dfh->add_space("a4", Cs[10]);
    dfh->add_space("b4", Cs[11]);

    dfh->add_transformation("Aar", "a", "r");
    dfh->add_transformation("Abs", "b", "s");
    dfh->add_transformation("Bas", "a", "s1");
    dfh->add_transformation("Bbr", "b", "r1");
    dfh->add_transformation("Cas", "a2", "s");
    dfh->add_transformation("Cbr", "b2", "r");
    dfh->add_transformation("Dar", "a", "r3");
    dfh->add_transformation("Dbs", "b", "s3");
    dfh->add_transformation("Ear", "a4", "r");
    dfh->add_transformation("Ebs", "b4", "s");

    dfh->transform();

    Cr1.reset();
    Cs1.reset();
    Ca2.reset();
    Cb2.reset();
    Cr3.reset();
    Cs3.reset();
    Ca4.reset();
    Cb4.reset();
    Cs.clear();
    dfh->clear_spaces();

    // => Blocking ... figure out how big a tensor slice to handle at a time <= //

    long int overhead = 0L;
    overhead += 2L * nT * nr * ns; // Thread work arrays Trs and Vrs below
    overhead += 2L * na * ns + 2L * nb * nr + 2L * na * nr + 2L * nb * ns; // Sas, Sbr, sBar, sAbs, Qas, Qbr, Qar, Qbs
    // account for a few of the smaller matrices already defined, but not exhaustively
    overhead += 12L * nn * nn; // D, V, J, K, P, and C matrices for A and B (neglecting C)

    long int rem = doubles_ - overhead;

    outfile->Printf("    %ld doubles - %ld overhead leaves %ld for dispersion\n", doubles_, overhead, rem);

    if (rem < 0L) {
        throw PSIEXCEPTION("Too little static memory for DFTSAPT::mp2_terms");
    }

    long int cost_a = 2L * nr * nQ + 2L * ns * nQ; // how much mem for Aar, Bas, Cas, Dar for a single a
    // cost_b would be the same value, and would be how much mem for Abs, Bbr, Cbr, Dbs for a single b
    long int max_a_l = rem / (2L * cost_a);
    long int max_b_l = max_a_l;
    int max_a = (max_a_l > na ? na : (int) max_a_l);
    int max_b = (max_b_l > nb ? nb : (int) max_b_l);
    if (max_a < 1 || max_b < 1) {
        throw PSIEXCEPTION("Too little dynamic memory for DFTSAPT::mp2_terms");
    }
    int nablocks = (na / max_a);
    if (na % max_a) nablocks++;
    int nbblocks = (nb / max_b);
    if (nb % max_b) nbblocks++;
    outfile->Printf("    Processing a single (a,b) pair requires %ld doubles\n", cost_a * 2L);
    outfile->Printf("    %d values of a processed in %d blocks of %d\n", na, nablocks, max_a);
    outfile->Printf("    %d values of b processed in %d blocks of %d\n\n", nb, nbblocks, max_b);

    // => Tensor Slices <= //

    auto Aar = std::make_shared<Matrix>("Aar", max_a * nr, nQ);
    auto Abs = std::make_shared<Matrix>("Abs", max_b * ns, nQ);
    auto Bas = std::make_shared<Matrix>("Bas", max_a * ns, nQ);
    auto Bbr = std::make_shared<Matrix>("Bbr", max_b * nr, nQ);
    auto Cas = std::make_shared<Matrix>("Cas", max_a * ns, nQ);
    auto Cbr = std::make_shared<Matrix>("Cbr", max_b * nr, nQ);
    auto Dar = std::make_shared<Matrix>("Dar", max_a * nr, nQ);
    auto Dbs = std::make_shared<Matrix>("Dbs", max_b * ns, nQ);

    // => Thread Work Arrays <= //

    std::vector<std::shared_ptr<Matrix> > Trs;
    std::vector<std::shared_ptr<Matrix> > Vrs;
    for (int t = 0; t < nT; t++) {
        Trs.push_back(std::make_shared<Matrix>("Trs", nr, ns));
        Vrs.push_back(std::make_shared<Matrix>("Vrs", nr, ns));
    }

    // => Pointers <= //

    double** Aarp = Aar->pointer();
    double** Absp = Abs->pointer();
    double** Basp = Bas->pointer();
    double** Bbrp = Bbr->pointer();
    double** Casp = Cas->pointer();
    double** Cbrp = Cbr->pointer();
    double** Darp = Dar->pointer();
    double** Dbsp = Dbs->pointer();

    double** Sasp = Sas->pointer();
    double** Sbrp = Sbr->pointer();
    double** SBarp = SBar->pointer();
    double** SAbsp = SAbs->pointer();

    double** Qasp = Qas->pointer();
    double** Qbrp = Qbr->pointer();
    double** Qarp = Qar->pointer();
    double** Qbsp = Qbs->pointer();

    double* eap = eps_occ0A->pointer();
    double* ebp = eps_occ0B->pointer();
    double* erp = eps_vir0A->pointer();
    double* esp = eps_vir0B->pointer();

    // => Slice D + E -> D <= //

    dfh->add_disk_tensor("Far", std::make_tuple(na, nr, nQ));

    for (size_t astart = 0; astart < na; astart += max_a) {
        size_t nablock = (astart + max_a >= na ? na - astart : max_a);

        dfh->fill_tensor("Dar", Dar, {astart, astart + nablock});
        dfh->fill_tensor("Ear", Aar, {astart, astart + nablock});

        double* D2p = Darp[0];
        double* A2p = Aarp[0];
        for (long int arQ = 0L; arQ < nablock * nrQ; arQ++) {
            (*D2p++) += (*A2p++);
        }
        dfh->write_disk_tensor("Far", Dar, {astart, astart + nablock});
    }

    dfh->add_disk_tensor("Fbs", std::make_tuple(nb, ns, nQ));

    for (size_t bstart = 0; bstart < nb; bstart += max_b) {
        size_t nbblock = (bstart + max_b >= nb ? nb - bstart : max_b);

        dfh->fill_tensor("Dbs", Dbs, {bstart, bstart + nbblock});
        dfh->fill_tensor("Ebs", Abs, {bstart, bstart + nbblock});

        double* D2p = Dbsp[0];
        double* A2p = Absp[0];
        for (long int bsQ = 0L; bsQ < nbblock * nsQ; bsQ++) {
            (*D2p++) += (*A2p++);
        }
        dfh->write_disk_tensor("Fbs", Dbs, {bstart, bstart + nbblock});
    }

    // => Targets <= //

    double Disp20 = 0.0;
    double ExchDisp20 = 0.0;

    // ==> Master Loop <== //

    for (size_t astart = 0; astart < na; astart += max_a) {
        size_t nablock = (astart + max_a >= na ? na - astart : max_a);

        dfh->fill_tensor("Aar", Aar, {astart, astart + nablock});
        dfh->fill_tensor("Bas", Bas, {astart, astart + nablock});
        dfh->fill_tensor("Cas", Cas, {astart, astart + nablock});
        dfh->fill_tensor("Far", Dar, {astart, astart + nablock});

        for (size_t bstart = 0; bstart < nb; bstart += max_b) {
            size_t nbblock = (bstart + max_b >= nb ? nb - bstart : max_b);

            dfh->fill_tensor("Abs", Abs, {bstart, bstart + nbblock});
            dfh->fill_tensor("Bbr", Bbr, {bstart, bstart + nbblock});
            dfh->fill_tensor("Cbr", Cbr, {bstart, bstart + nbblock});
            dfh->fill_tensor("Fbs", Dbs, {bstart, bstart + nbblock});

            long int nab = nablock * nbblock;

#pragma omp parallel for schedule(dynamic) reduction(+ : Disp20, ExchDisp20)
            for (long int ab = 0L; ab < nab; ab++) {
                int a = ab / nbblock;
                int b = ab % nbblock;

                int thread = 0;
#ifdef _OPENMP
                thread = omp_get_thread_num();
#endif

                double** Trsp = Trs[thread]->pointer();
                double** Vrsp = Vrs[thread]->pointer();

                // => Amplitudes, Disp20 <= //

                C_DGEMM('N', 'T', nr, ns, nQ, 1.0, Aarp[(a)*nr], nQ, Absp[(b)*ns], nQ, 0.0, Vrsp[0], ns);

                for (int r = 0; r < nr; r++) {
                    for (int s = 0; s < ns; s++) {
                        Trsp[r][s] = Vrsp[r][s] / (eap[a + astart] + ebp[b + bstart] - erp[r] - esp[s]);
                        Disp20 += 4.0 * Trsp[r][s] * Vrsp[r][s];
                    }
                }

                // => Exch-Disp20 <= //

                // > Q1-Q3 < //

                C_DGEMM('N', 'T', nr, ns, nQ, 1.0, Bbrp[(b)*nr], nQ, Basp[(a)*ns], nQ, 0.0, Vrsp[0], ns);
                C_DGEMM('N', 'T', nr, ns, nQ, 1.0, Cbrp[(b)*nr], nQ, Casp[(a)*ns], nQ, 1.0, Vrsp[0], ns);
                C_DGEMM('N', 'T', nr, ns, nQ, 1.0, Aarp[(a)*nr], nQ, Dbsp[(b)*ns], nQ, 1.0, Vrsp[0], ns);
                C_DGEMM('N', 'T', nr, ns, nQ, 1.0, Darp[(a)*nr], nQ, Absp[(b)*ns], nQ, 1.0, Vrsp[0], ns);

                // > V,J,K < //

                C_DGER(nr, ns, 1.0, Qbrp[b + bstart], 1, Sasp[a + astart], 1, Vrsp[0], ns);
                C_DGER(nr, ns, 1.0, Sbrp[b + bstart], 1, Qasp[a + astart], 1, Vrsp[0], ns);
                C_DGER(nr, ns, 1.0, Qarp[a + astart], 1, SAbsp[b + bstart], 1, Vrsp[0], ns);
                C_DGER(nr, ns, 1.0, SBarp[a + astart], 1, Qbsp[b + bstart], 1, Vrsp[0], ns);

                for (int r = 0; r < nr; r++) {
                    for (int s = 0; s < ns; s++) {
                        ExchDisp20 -= 2.0 * Trsp[r][s] * Vrsp[r][s];
                    }
                }
            }
        }
    }

    scalars_["Disp20"] = Disp20;
    scalars_["Exch-Disp20"] = ExchDisp20;
    if (do_print) {
        outfile->Printf("    Disp20              = %18.12lf [Eh]\n", Disp20);
        outfile->Printf("    Exch-Disp20         = %18.12lf [Eh]\n", ExchDisp20);
        outfile->Printf("\n");
    }
}

// Compute total dispersion energy and S-infinity version of total exchange-dispersion
// This code applies to regular SAPT, not FISAPT, and is called from the SAPT workflow when DO_DISP_EXCH_SINF == true.
void FISAPT::sinf_disp(std::map<std::string, SharedMatrix> matrix_cache, std::map<std::string, SharedVector> vector_cache,
                        bool do_print) {
    if (do_print) {
        outfile->Printf("  ==> Dispersion <==\n\n");
    }

    // => Pointers <= //

    std::shared_ptr<Matrix> Cocc0A = matrix_cache["Caocc0A"];
    std::shared_ptr<Matrix> Cocc0B = matrix_cache["Caocc0B"];
    std::shared_ptr<Matrix> Cvir0A = matrix_cache["Cvir0A"];
    std::shared_ptr<Matrix> Cvir0B = matrix_cache["Cvir0B"];

    std::shared_ptr<Vector> eps_occ0A = vector_cache["eps_aocc0A"];
    std::shared_ptr<Vector> eps_occ0B = vector_cache["eps_aocc0B"];
    std::shared_ptr<Vector> eps_vir0A = vector_cache["eps_vir0A"];
    std::shared_ptr<Vector> eps_vir0B = vector_cache["eps_vir0B"];

    // => Auxiliary Basis Set <= //
    std::shared_ptr<BasisSet> auxiliary = reference_->get_basisset("DF_BASIS_SAPT");

    // => Sizing <= //

    int nn = primary_->nbf();
    int na = Cocc0A->colspi()[0];
    int nb = Cocc0B->colspi()[0];
    int nr = Cvir0A->colspi()[0];
    int ns = Cvir0B->colspi()[0];
    int nQ = auxiliary->nbf();
    size_t nrQ = nr * (size_t)nQ;
    size_t nsQ = ns * (size_t)nQ;

    int nT = 1;
#ifdef _OPENMP
    nT = Process::environment.get_n_threads();
#endif

    // => Stashed Variables <= //

    std::shared_ptr<Matrix> S = matrix_cache["S"];
    std::shared_ptr<Matrix> V_A = matrix_cache["V_A"];
    std::shared_ptr<Matrix> V_B = matrix_cache["V_B"];

    // => Intermolecular overlap matrix and inverse <= //
    std::shared_ptr<Matrix> Sab = linalg::triplet(Cocc0A, S, Cocc0B, true, false, false);
    double** Sabp = Sab->pointer();
    auto D = std::make_shared<Matrix>("D", na + nb, na + nb);
    D->identity();
    double** Dp = D->pointer();
    for (int a = 0; a < na; a++) {
        for (int b = 0; b < nb; b++) {
            Dp[a][b + na] = Dp[b + na][a] = Sabp[a][b];
        }
    }
    D->power(-1.0, 1.0E-12);
    Dp = D->pointer();

    // => New Stuff <= //
    // Start with T's
    std::shared_ptr<Matrix> Sbr = linalg::triplet(Cocc0B, S, Cvir0A, true, false, false);
    std::shared_ptr<Matrix> Sas = linalg::triplet(Cocc0A, S, Cvir0B, true, false, false);
    auto Tar = std::make_shared<Matrix>("Tar", na, nr);
    auto Tbr = std::make_shared<Matrix>("Tbr", nb, nr);
    auto Tas = std::make_shared<Matrix>("Tas", na, ns);
    auto Tbs = std::make_shared<Matrix>("Tbs", nb, ns);

    C_DGEMM('N', 'N', na, nr, nb, 1.0, &Dp[0][na], na + nb, Sbr->pointer()[0], nr, 0.0,
            Tar->pointer()[0], nr);
    C_DGEMM('N', 'N', nb, nr, nb, 1.0, &Dp[na][na], na + nb, Sbr->pointer()[0], nr, 0.0,
            Tbr->pointer()[0], nr);
    C_DGEMM('N', 'N', na, ns, na, 1.0, &Dp[0][0], na + nb, Sas->pointer()[0], ns, 0.0,
            Tas->pointer()[0], ns);
    C_DGEMM('N', 'N', nb, ns, na, 1.0, &Dp[na][0], na + nb, Sas->pointer()[0], ns, 0.0,
            Tbs->pointer()[0], ns);

    // C1's and C2's from D's and C's.
    // C1's are C times D diagonal blocks.
    // C2's are times off-diagonal blocks.
    auto C1a = std::make_shared<Matrix>("C1a", nn, na);
    auto C1b = std::make_shared<Matrix>("C1b", nn, nb);
    auto C2a = std::make_shared<Matrix>("C2a", nn, na);
    auto C2b = std::make_shared<Matrix>("C2b", nn, nb);

    C_DGEMM('N', 'N', nn, na, na, 1.0, Cocc0A->pointer()[0], na, &Dp[0][0], na + nb, 0.0,
            C1a->pointer()[0], na);
    C_DGEMM('N', 'N', nn, nb, nb, 1.0, Cocc0B->pointer()[0], nb, &Dp[na][na], na + nb, 0.0,
            C1b->pointer()[0], nb);
    C_DGEMM('N', 'N', nn, na, nb, 1.0, Cocc0B->pointer()[0], nb, &Dp[na][0], na + nb, 0.0,
            C2a->pointer()[0], na);
    C_DGEMM('N', 'N', nn, nb, na, 1.0, Cocc0A->pointer()[0], na, &Dp[0][na], na + nb, 0.0,
            C2b->pointer()[0], nb);

    // Coeffs for all occupied
    std::vector<std::shared_ptr<Matrix> > hold_these;
    hold_these.push_back(Cocc0A);
    hold_these.push_back(Cocc0B);

    auto Cocc0AB = linalg::horzcat(hold_these);
    hold_these.clear();

    // Half transform D_ia and D_ib for JK
    auto D_Ni_a = std::make_shared<Matrix>("D_Ni_a", nn, na + nb);
    auto D_Ni_b = std::make_shared<Matrix>("D_Ni_b", nn, na + nb);

    C_DGEMM('N', 'N', nn, na + nb, na, 1.0, Cocc0A->pointer()[0], na, &Dp[0][0], na + nb, 0.0,
            D_Ni_a->pointer()[0], na + nb);
    C_DGEMM('N', 'N', nn, na + nb, nb, 1.0, Cocc0B->pointer()[0], nb, &Dp[na][0], na + nb, 0.0,
            D_Ni_b->pointer()[0], na + nb);

    // Make JK's
    jk_ = JK::build_JK(primary_, reference_->get_basisset("DF_BASIS_SCF"), options_, false, doubles_);
    jk_->set_memory(doubles_);
    jk_->set_do_J(true);
    jk_->set_do_K(true);
    jk_->initialize();
    jk_->print_header();

    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    Cl.clear();
    Cr.clear();
    Cl.push_back(Cocc0AB);
    Cr.push_back(D_Ni_a);
    Cl.push_back(Cocc0AB);
    Cr.push_back(D_Ni_b);
    jk_->compute();

    std::shared_ptr<Matrix> J_D_ia = J[0];
    std::shared_ptr<Matrix> K_D_ia = K[0];

    std::shared_ptr<Matrix> J_D_ib = J[1];
    std::shared_ptr<Matrix> K_D_ib = K[1];

    // Finish D_ia and D_ib transformation to make tilded C's
    auto D_ia = linalg::doublet(Cocc0AB, D_Ni_a, false, true);
    auto D_ib = linalg::doublet(Cocc0AB, D_Ni_b, false, true);

    auto Ct_Kr = linalg::triplet(D_ib, S, Cvir0A, false, false, false);
    Ct_Kr->scale(-1);
    Ct_Kr->add(Cvir0A);
    auto Ct_Ks = linalg::triplet(D_ia, S, Cvir0B, false, false, false);
    Ct_Ks->scale(-1);
    Ct_Ks->add(Cvir0B);

    // Make omega cores
    std::shared_ptr<Matrix> AJK(J_D_ia->clone());
    AJK->zero();
    AJK->add(J_D_ia);
    AJK->scale(2);
    AJK->add(V_A);
    AJK->subtract(K_D_ia);

    auto AJK_ar = linalg::triplet(C2a, AJK, Ct_Kr, true, false, false);
    auto AJK_as = linalg::triplet(C2a, AJK, Ct_Ks, true, false, false);
    auto AJK_br = linalg::triplet(C1b, AJK, Ct_Kr, true, false, false);
    auto AJK_bs = linalg::triplet(C1b, AJK, Ct_Ks, true, false, false);

    std::shared_ptr<Matrix> BJK(J_D_ib->clone());
    BJK->zero();
    BJK->add(J_D_ib);
    BJK->scale(2);
    BJK->add(V_B);
    BJK->subtract(K_D_ib);

    auto BJK_ar = linalg::triplet(C1a, BJK, Ct_Kr, true, false, false);
    auto BJK_as = linalg::triplet(C1a, BJK, Ct_Ks, true, false, false);
    auto BJK_br = linalg::triplet(C2b, BJK, Ct_Kr, true, false, false);
    auto BJK_bs = linalg::triplet(C2b, BJK, Ct_Ks, true, false, false);

    // Finish omega terms
    std::shared_ptr<Matrix> omega_ar(AJK_ar->clone());
    omega_ar->zero();
    omega_ar->add(AJK_ar);
    omega_ar->add(BJK_ar);
    omega_ar->scale(4);

    std::shared_ptr<Matrix> omega_as(AJK_as->clone());
    omega_as->zero();
    omega_as->add(AJK_as);
    omega_as->add(BJK_as);
    omega_as->scale(2);

    std::shared_ptr<Matrix> omega_br(AJK_br->clone());
    omega_br->zero();
    omega_br->add(AJK_br);
    omega_br->add(BJK_br);
    omega_br->scale(2);

    std::shared_ptr<Matrix> omega_bs(AJK_bs->clone());
    omega_bs->zero();
    omega_bs->add(AJK_bs);
    omega_bs->add(BJK_bs);
    omega_bs->scale(4);

    D.reset();
    Sbr.reset();
    Sas.reset();
    Cocc0AB.reset();
    D_Ni_a.reset();
    D_Ni_b.reset();
    J_D_ia.reset();
    K_D_ia.reset();
    J_D_ib.reset();
    K_D_ib.reset();
    D_ia.reset();
    D_ib.reset();
    AJK.reset();
    BJK.reset();
    AJK_ar.reset();
    AJK_as.reset();
    AJK_br.reset();
    AJK_bs.reset();
    BJK_ar.reset();
    BJK_as.reset();
    BJK_br.reset();
    BJK_bs.reset();

    // => Memory <= //

    std::vector<std::shared_ptr<Matrix> > Cs;
    Cs.push_back(Cocc0A);
    Cs.push_back(Cvir0A);
    Cs.push_back(Cocc0B);
    Cs.push_back(Cvir0B);
    Cs.push_back(C1a);
    Cs.push_back(C1b);
    Cs.push_back(C2a);
    Cs.push_back(C2b);
    Cs.push_back(Ct_Kr);
    Cs.push_back(Ct_Ks);

    size_t max_MO = 0, ncol = 0;
    for (auto& mat : Cs) {
        max_MO = std::max(max_MO, (size_t)mat->ncol());
        ncol += (size_t)mat->ncol();
    }

    // => Get integrals from DFHelper <= //
    auto dfh(std::make_shared<DFHelper>(primary_, auxiliary));
    dfh->set_memory(doubles_ - Cs[0]->nrow() * ncol);
    dfh->set_method("DIRECT_iaQ");
    dfh->set_nthreads(nT);
    dfh->initialize();
    dfh->print_header();

    dfh->add_space("a", Cs[0]);
    dfh->add_space("r", Cs[1]);
    dfh->add_space("b", Cs[2]);
    dfh->add_space("s", Cs[3]);
    dfh->add_space("a1", Cs[4]);
    dfh->add_space("b1", Cs[5]);
    dfh->add_space("a2", Cs[6]);
    dfh->add_space("b2", Cs[7]);
    dfh->add_space("r1", Cs[8]);
    dfh->add_space("s1", Cs[9]);

    dfh->add_transformation("Aar", "a", "r");
    dfh->add_transformation("Abs", "b", "s");
    dfh->add_transformation("Bar", "a1", "r1");
    dfh->add_transformation("Bbs", "b1", "s1");
    dfh->add_transformation("Cbr", "b2", "r1");
    dfh->add_transformation("Cas", "a2", "s1");
    dfh->add_transformation("Das", "a1", "s1");
    dfh->add_transformation("Dbr", "b1", "r1");
    dfh->add_transformation("Ebs", "b2", "s1");
    dfh->add_transformation("Ear", "a2", "r1");

    dfh->transform();

    C1a.reset();
    C1b.reset();
    C2a.reset();
    C2b.reset();
    Ct_Kr.reset();
    Ct_Ks.reset();
    Cs.clear();
    dfh->clear_spaces();

    // => Blocking <= //

    long int overhead = 0L;
    overhead += 2L * nT * nr * ns;
    overhead += 2L * na * ns + 2L * nb * nr + 2L * na * nr + 2L * nb * ns;
    long int rem = doubles_ - overhead;

    if (rem < 0L) {
        throw PSIEXCEPTION("Too little static memory for DFTSAPT::mp2_terms");
    }

    long int cost_a = 2L * nr * nQ + 2L * ns * nQ;
    long int max_a = rem / (2L * cost_a);
    long int max_b = max_a;
    max_a = (max_a > na ? na : max_a);
    max_b = (max_b > nb ? nb : max_b);
    if (max_a < 1L || max_b < 1L) {
        throw PSIEXCEPTION("Too little dynamic memory for DFTSAPT::mp2_terms");
    }

    // => Tensor Slices <= //

    auto Aar = std::make_shared<Matrix>("Aar", max_a * nr, nQ);
    auto Abs = std::make_shared<Matrix>("Abs", max_b * ns, nQ);
    auto Bar = std::make_shared<Matrix>("Bar", max_a * nr, nQ);
    auto Bbs = std::make_shared<Matrix>("Bbs", max_b * ns, nQ);
    auto Cbr = std::make_shared<Matrix>("Cbr", max_b * nr, nQ);
    auto Cas = std::make_shared<Matrix>("Cas", max_a * ns, nQ);
    auto Das = std::make_shared<Matrix>("Das", max_a * ns, nQ);
    auto Dbr = std::make_shared<Matrix>("Dbr", max_b * nr, nQ);
    auto Ebs = std::make_shared<Matrix>("Ebs", max_b * ns, nQ);
    auto Ear = std::make_shared<Matrix>("Ear", max_a * nr, nQ);

    // => Thread Work Arrays <= //

    std::vector<std::shared_ptr<Matrix> > Trs;
    std::vector<std::shared_ptr<Matrix> > Vrs;
    for (int t = 0; t < nT; t++) {
        Trs.push_back(std::make_shared<Matrix>("Trs", nr, ns));
        Vrs.push_back(std::make_shared<Matrix>("Vrs", nr, ns));
    }

    // => Pointers <= //
    double** Aarp = Aar->pointer();
    double** Absp = Abs->pointer();
    double** Barp = Bar->pointer();
    double** Bbsp = Bbs->pointer();
    double** Cbrp = Cbr->pointer();
    double** Casp = Cas->pointer();
    double** Dasp = Das->pointer();
    double** Dbrp = Dbr->pointer();
    double** Ebsp = Ebs->pointer();
    double** Earp = Ear->pointer();

    double** Tarp = Tar->pointer();
    double** Tbrp = Tbr->pointer();
    double** Tasp = Tas->pointer();
    double** Tbsp = Tbs->pointer();

    double** omega_arp = omega_ar->pointer();
    double** omega_asp = omega_as->pointer();
    double** omega_brp = omega_br->pointer();
    double** omega_bsp = omega_bs->pointer();

    double* eap = eps_occ0A->pointer();
    double* ebp = eps_occ0B->pointer();
    double* erp = eps_vir0A->pointer();
    double* esp = eps_vir0B->pointer();

    // => Targets <= //

    double Disp20 = 0.0;
    double CompleteDisp20 = 0.0;

    // ==> Master Loop <== //

    for (size_t astart = 0; astart < na; astart += max_a) {
        size_t nablock = (astart + max_a >= na ? na - astart : max_a);

        dfh->fill_tensor("Aar", Aar, {astart, astart + nablock});
        dfh->fill_tensor("Bar", Bar, {astart, astart + nablock});
        dfh->fill_tensor("Cas", Cas, {astart, astart + nablock});
        dfh->fill_tensor("Das", Das, {astart, astart + nablock});
        dfh->fill_tensor("Ear", Ear, {astart, astart + nablock});

        for (size_t bstart = 0; bstart < nb; bstart += max_b) {
            size_t nbblock = (bstart + max_b >= nb ? nb - bstart : max_b);

            dfh->fill_tensor("Abs", Abs, {bstart, bstart + nbblock});
            dfh->fill_tensor("Bbs", Bbs, {bstart, bstart + nbblock});
            dfh->fill_tensor("Cbr", Cbr, {bstart, bstart + nbblock});
            dfh->fill_tensor("Dbr", Dbr, {bstart, bstart + nbblock});
            dfh->fill_tensor("Ebs", Ebs, {bstart, bstart + nbblock});

            long int nab = nablock * nbblock;

#pragma omp parallel for schedule(dynamic) reduction(+ : Disp20, CompleteDisp20)
            for (long int ab = 0L; ab < nab; ab++) {
                int a = ab / nbblock;
                int b = ab % nbblock;

                int thread = 0;
#ifdef _OPENMP
                thread = omp_get_thread_num();
#endif

                double** Trsp = Trs[thread]->pointer();
                double** Vrsp = Vrs[thread]->pointer();

                // => Amplitudes, Disp20 <= //

                C_DGEMM('N', 'T', nr, ns, nQ, 1.0, Aarp[(a)*nr], nQ, Absp[(b)*ns], nQ, 0.0, Vrsp[0], ns);

                for (int r = 0; r < nr; r++) {
                    for (int s = 0; s < ns; s++) {
                        Trsp[r][s] = Vrsp[r][s] / (eap[a + astart] + ebp[b + bstart] - erp[r] - esp[s]);
                        Disp20 += 4.0 * Trsp[r][s] * Vrsp[r][s];
                    }
                }

                // => Exch-Disp20 <= //

                // > DF-Part < //

                C_DGEMM('N', 'T', nr, ns, nQ, 4.0, Barp[(a)*nr], nQ, Bbsp[(b)*ns], nQ, 0.0, Vrsp[0], ns);
                C_DGEMM('N', 'T', nr, ns, nQ, -2.0, Cbrp[(b)*nr], nQ, Casp[(a)*ns], nQ, 1.0, Vrsp[0], ns);
                C_DGEMM('N', 'T', nr, ns, nQ, -2.0, Dbrp[(b)*nr], nQ, Dasp[(a)*ns], nQ, 1.0, Vrsp[0], ns);
                C_DGEMM('N', 'T', nr, ns, nQ, 4.0, Earp[(a)*nr], nQ, Ebsp[(b)*ns], nQ, 1.0, Vrsp[0], ns);

                // > AO-Part < //

                C_DGER(nr, ns, 1.0, Tarp[a + astart], 1, omega_bsp[b + bstart], 1, Vrsp[0], ns);
                C_DGER(nr, ns, -1.0, omega_brp[b + bstart], 1, Tasp[a + astart], 1, Vrsp[0], ns);
                C_DGER(nr, ns, 1.0, omega_arp[a + astart], 1, Tbsp[b + bstart], 1, Vrsp[0], ns);
                C_DGER(nr, ns, -1.0, Tbrp[b + bstart], 1, omega_asp[a + astart], 1, Vrsp[0], ns);
                for (int r = 0; r < nr; r++) {
                    for (int s = 0; s < ns; s++) {
                        CompleteDisp20 += Trsp[r][s] * Vrsp[r][s];
                    }
                }
            }
        }
    }

    double ExchDisp20 = CompleteDisp20 - Disp20;

    scalars_["Disp20"] = Disp20;
    scalars_["Exch-Disp20 (S^inf)"] = ExchDisp20;
    Process::environment.globals["SAPT EXCH-DISP20(S^INF) ENERGY"] = scalars_["Exch-Disp20 (S^inf)"];
    if (do_print) {
        outfile->Printf("    Disp20              = %18.12lf [Eh]\n", Disp20);
        outfile->Printf("    Exch-Disp20 (S^inf) = %18.12lf [Eh]\n", ExchDisp20);
        outfile->Printf("\n");
    }
}

void FISAPT::print_trailer() {
    scalars_["Electrostatics"] = scalars_["Elst10,r"];
    if (reference_->has_potential_variable("A") && reference_->has_potential_variable("B")) {
        scalars_["Electrostatics"] += scalars_["Extern-Extern"];
    }
    scalars_["Exchange"] = scalars_["Exch10"];
    scalars_["Induction"] = scalars_["Ind20,r"] + scalars_["Exch-Ind20,r"] + scalars_["delta HF,r (2)"];
    scalars_["sInduction"] = scalars_["Ind20,r"] + scalars_["sExch-Ind20,r"] + scalars_["delta HF,r (2)"];
    scalars_["Dispersion"] = scalars_["Disp20"] + scalars_["Exch-Disp20"];
    scalars_["sDispersion"] = scalars_["Disp20"] + scalars_["sExch-Disp20"];
    scalars_["SAPT"] =
        scalars_["Electrostatics"] + scalars_["Exchange"] + scalars_["Induction"] + scalars_["Dispersion"];
    scalars_["sSAPT"] =
        scalars_["Electrostatics"] + scalars_["Exchange"] + scalars_["sInduction"] + scalars_["sDispersion"];

    double Sdelta = scalars_["Induction"] / (scalars_["Ind20,r"] + scalars_["Exch-Ind20,r"]);
    scalars_["Induction (A<-B)"] = Sdelta * (scalars_["Ind20,r (A<-B)"] + scalars_["Exch-Ind20,r (A<-B)"]);
    scalars_["Induction (B<-A)"] = Sdelta * (scalars_["Ind20,r (B<-A)"] + scalars_["Exch-Ind20,r (B<-A)"]);

    outfile->Printf("  ==> Results <==\n\n");

    outfile->Printf("\n    SAPT Results  \n");
    std::string scaled = "   ";
    outfile->Printf(
        "  --------------------------------------------------------------------------------------------------------\n");
    outfile->Printf("    Electrostatics            %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scalars_["Electrostatics"] * 1000.0, scalars_["Electrostatics"] * pc_hartree2kcalmol,
                    scalars_["Electrostatics"] * pc_hartree2kJmol);
    outfile->Printf("      Elst10,r                %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scalars_["Elst10,r"] * 1000.0, scalars_["Elst10,r"] * pc_hartree2kcalmol,
                    scalars_["Elst10,r"] * pc_hartree2kJmol);

    if (reference_->has_potential_variable("A") && reference_->has_potential_variable("B")) {
        outfile->Printf("      Extern-Extern           %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                        scalars_["Extern-Extern"] * 1000.0, scalars_["Extern-Extern"] * pc_hartree2kcalmol,
                        scalars_["Extern-Extern"] * pc_hartree2kJmol);
    }
    outfile->Printf("\n");

    outfile->Printf("    Exchange %3s              %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n", scaled.c_str(),
                    scalars_["Exchange"] * 1000.0, scalars_["Exchange"] * pc_hartree2kcalmol,
                    scalars_["Exchange"] * pc_hartree2kJmol);
    outfile->Printf("      Exch10                  %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scalars_["Exch10"] * 1000.0, scalars_["Exch10"] * pc_hartree2kcalmol,
                    scalars_["Exch10"] * pc_hartree2kJmol);
    outfile->Printf("      Exch10(S^2)             %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n\n",
                    scalars_["Exch10(S^2)"] * 1000.0, scalars_["Exch10(S^2)"] * pc_hartree2kcalmol,
                    scalars_["Exch10(S^2)"] * pc_hartree2kJmol);

    outfile->Printf("    Induction %3s             %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n", scaled.c_str(),
                    scalars_["Induction"] * 1000.0, scalars_["Induction"] * pc_hartree2kcalmol,
                    scalars_["Induction"] * pc_hartree2kJmol);
    outfile->Printf("      Ind20,r                 %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scalars_["Ind20,r"] * 1000.0, scalars_["Ind20,r"] * pc_hartree2kcalmol,
                    scalars_["Ind20,r"] * pc_hartree2kJmol);
    outfile->Printf("      Exch-Ind20,r %3s        %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n", scaled.c_str(),
                    scalars_["Exch-Ind20,r"] * 1000.0, scalars_["Exch-Ind20,r"] * pc_hartree2kcalmol,
                    scalars_["Exch-Ind20,r"] * pc_hartree2kJmol);
    outfile->Printf("      delta HF,r (2) %3s      %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n", scaled.c_str(),
                    scalars_["delta HF,r (2)"] * 1000.0, scalars_["delta HF,r (2)"] * pc_hartree2kcalmol,
                    scalars_["delta HF,r (2)"] * pc_hartree2kJmol);
    outfile->Printf("      Induction (A<-B) %3s    %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n", scaled.c_str(),
                    scalars_["Induction (A<-B)"] * 1000.0, scalars_["Induction (A<-B)"] * pc_hartree2kcalmol,
                    scalars_["Induction (A<-B)"] * pc_hartree2kJmol);
    outfile->Printf("      Induction (B<-A) %3s    %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n\n",
                    scaled.c_str(), scalars_["Induction (B<-A)"] * 1000.0,
                    scalars_["Induction (B<-A)"] * pc_hartree2kcalmol, scalars_["Induction (B<-A)"] * pc_hartree2kJmol);

    outfile->Printf("    Dispersion %3s            %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n", scaled.c_str(),
                    scalars_["Dispersion"] * 1000.0, scalars_["Dispersion"] * pc_hartree2kcalmol,
                    scalars_["Dispersion"] * pc_hartree2kJmol);
    outfile->Printf("      Disp20                  %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scalars_["Disp20"] * 1000.0, scalars_["Disp20"] * pc_hartree2kcalmol,
                    scalars_["Disp20"] * pc_hartree2kJmol);
    outfile->Printf("      Exch-Disp20 %3s         %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n\n",
                    scaled.c_str(), scalars_["Exch-Disp20"] * 1000.0, scalars_["Exch-Disp20"] * pc_hartree2kcalmol,
                    scalars_["Exch-Disp20"] * pc_hartree2kJmol);

    outfile->Printf("  Total HF                    %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                    scalars_["HF"] * 1000.0, scalars_["HF"] * pc_hartree2kcalmol, scalars_["HF"] * pc_hartree2kJmol);
    outfile->Printf("  Total SAPT0 %3s             %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n", scaled.c_str(),
                    scalars_["SAPT"] * 1000.0, scalars_["SAPT"] * pc_hartree2kcalmol,
                    scalars_["SAPT"] * pc_hartree2kJmol);
    if (options_.get_bool("SSAPT0_SCALE")) {
        outfile->Printf("  Total sSAPT0 %3s            %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
                        scaled.c_str(), scalars_["sSAPT"] * 1000.0, scalars_["sSAPT"] * pc_hartree2kcalmol,
                        scalars_["sSAPT"] * pc_hartree2kJmol);
    }
    outfile->Printf("\n");
    outfile->Printf(
        "  --------------------------------------------------------------------------------------------------------\n");

    outfile->Printf("    Han Solo: This is *not* gonna work.\n");
    outfile->Printf("    Luke Skywalker: Why didn't you say so before?\n");
    outfile->Printf("    Han Solo: I *did* say so before.\n");

    Process::environment.globals["SAPT ELST ENERGY"] = scalars_["Electrostatics"];
    Process::environment.globals["SAPT ELST10,R ENERGY"] = scalars_["Elst10,r"];
    if (reference_->has_potential_variable("A") && reference_->has_potential_variable("B")) {
        Process::environment.globals["SAPT ELST EXTERN-EXTERN ENERGY"] = scalars_["Extern-Extern"];
    }

    Process::environment.globals["SAPT EXCH ENERGY"] = scalars_["Exchange"];
    Process::environment.globals["SAPT EXCH10 ENERGY"] = scalars_["Exch10"];
    Process::environment.globals["SAPT EXCH10(S^2) ENERGY"] = scalars_["Exch10(S^2)"];

    Process::environment.globals["SAPT IND ENERGY"] = scalars_["Induction"];
    Process::environment.globals["SAPT IND20,R ENERGY"] = scalars_["Ind20,r"];
    Process::environment.globals["SAPT EXCH-IND20,R ENERGY"] = scalars_["Exch-Ind20,r"];
    Process::environment.globals["SAPT IND20,U ENERGY"] = scalars_["Ind20,u"];
    Process::environment.globals["SAPT EXCH-IND20,U ENERGY"] = scalars_["Exch-Ind20,u"];

    Process::environment.globals["SAPT DISP ENERGY"] = scalars_["Dispersion"];
    Process::environment.globals["SAPT DISP20 ENERGY"] = scalars_["Disp20"];
    Process::environment.globals["SAPT EXCH-DISP20 ENERGY"] = scalars_["Exch-Disp20"];

    Process::environment.globals["SAPT0 TOTAL ENERGY"] = scalars_["SAPT"];
    Process::environment.globals["SAPT TOTAL ENERGY"] = scalars_["SAPT"];
    Process::environment.globals["CURRENT ENERGY"] = Process::environment.globals["SAPT TOTAL ENERGY"];

    // Export the components of dHF to Psi4 variables
    Process::environment.globals["SAPT HF(2) ENERGY ABC(HF)"] = scalars_["E_ABC_HF"];
    Process::environment.globals["SAPT HF(2) ENERGY AC(0)"] = scalars_["E_AC"];
    Process::environment.globals["SAPT HF(2) ENERGY BC(0)"] = scalars_["E_BC"];
    Process::environment.globals["SAPT HF(2) ENERGY A(0)"] = scalars_["E_A"];
    Process::environment.globals["SAPT HF(2) ENERGY B(0)"] = scalars_["E_B"];
    Process::environment.globals["SAPT HF(2) ENERGY AC(HF)"] = scalars_["E_AC_HF"];
    Process::environment.globals["SAPT HF(2) ENERGY BC(HF)"] = scalars_["E_BC_HF"];
    Process::environment.globals["SAPT HF(2) ENERGY AB(HF)"] = scalars_["E_AB_HF"];
    Process::environment.globals["SAPT HF(2) ENERGY A(HF)"] = scalars_["E_A_HF"];
    Process::environment.globals["SAPT HF(2) ENERGY B(HF)"] = scalars_["E_B_HF"];
    Process::environment.globals["SAPT HF(2) ENERGY C"] = scalars_["E_C"];
    Process::environment.globals["SAPT HF(2) ENERGY HF"] = scalars_["HF"];
}

void FISAPT::raw_plot(const std::string& filepath) {
// Konrad: there is at the moment some redundancy between do_cubes and raw_plot that should be eliminated. My bad.
    outfile->Printf("  ==> Scalar Field Plots <==\n\n");

    outfile->Printf("    F-SAPT Plot Filepath = %s\n\n", filepath.c_str());

    auto csg = std::make_shared<CubicScalarGrid>(primary_, options_);
    csg->set_filepath(filepath);
    csg->print_header();
    csg->set_auxiliary_basis(reference_->get_basisset("DF_BASIS_SCF"));

    /// Zeroth-order wavefunctions
    std::shared_ptr<Matrix> D_A = matrices_["D_A"];
    std::shared_ptr<Matrix> D_B = matrices_["D_B"];
    std::shared_ptr<Matrix> D_C = matrices_["D_C"];

    /// Fully interacting wavefunctions
    std::shared_ptr<Matrix> DFA = linalg::doublet(matrices_["LoccA"], matrices_["LoccA"], false, true);
    std::shared_ptr<Matrix> DFB = linalg::doublet(matrices_["LoccB"], matrices_["LoccB"], false, true);

    // => Density Fields <= //

    csg->compute_density(D_A, "DA");
    csg->compute_density(D_B, "DB");
    csg->compute_density(D_C, "DC");
    csg->compute_density(DFA, "DFA");
    csg->compute_density(DFB, "DFB");

    // => Difference Density Fields <= //

    DFA->subtract(D_A);
    DFB->subtract(D_B);

    csg->compute_density(DFA, "dDA");
    csg->compute_density(DFB, "dDB");

    // => ESP Fields <= //

    double* ZAp = vectors_["ZA"]->pointer();
    double* ZBp = vectors_["ZB"]->pointer();
    double* ZCp = vectors_["ZC"]->pointer();

    std::shared_ptr<Molecule> mol = primary_->molecule();

    std::vector<double> w_A(mol->natom());
    std::vector<double> w_B(mol->natom());
    std::vector<double> w_C(mol->natom());

    for (int A = 0; A < mol->natom(); A++) {
        w_A[A] = (ZAp[A]) / mol->Z(A);
        w_B[A] = (ZBp[A]) / mol->Z(A);
        w_C[A] = (ZCp[A]) / mol->Z(A);
    }

    D_A->scale(2.0);
    D_B->scale(2.0);
    D_C->scale(2.0);

    csg->compute_esp(D_A, w_A, "VA");
    csg->compute_esp(D_B, w_B, "VB");
    csg->compute_esp(D_C, w_C, "VC");

    D_A->scale(0.5);
    D_B->scale(0.5);
    D_C->scale(0.5);
}

void FISAPT::flocalize() {
    outfile->Printf("  ==> F-SAPT Localization (IBO) <==\n\n");

    // Currently always separating core and valence
    std::string link_assignment = options_.get_str("FISAPT_LINK_ASSIGNMENT");
    {
        outfile->Printf("  Local Orbitals for Monomer A:\n\n");

        int nn = matrices_["Caocc0A"]->rowspi()[0];
        int nf = matrices_["Cfocc0A"]->colspi()[0];
        int na = matrices_["Caocc0A"]->colspi()[0];
        int nm = nf + na;

        std::vector<int> ranges;
        ranges.push_back(0);
        ranges.push_back(nf);
        ranges.push_back(nm);

        std::shared_ptr<Matrix> Focc(
            new Matrix("Focc", vectors_["eps_occ0A"]->dimpi()[0], vectors_["eps_occ0A"]->dimpi()[0]));
        Focc->set_diagonal(vectors_["eps_occ0A"]);

        std::shared_ptr<fisapt::IBOLocalizer2> local =
            fisapt::IBOLocalizer2::build(primary_, reference_->get_basisset("MINAO"), matrices_["Cocc0A"], options_);
        local->print_header();
        std::map<std::string, std::shared_ptr<Matrix> > ret = local->localize(matrices_["Cocc0A"], Focc, ranges);

        matrices_["Locc0A"] = ret["L"];
        matrices_["Uocc0A"] = ret["U"];
        matrices_["Qocc0A"] = ret["Q"];

        matrices_["Lfocc0A"] = std::make_shared<Matrix>("Lfocc0A", nn, nf);
        matrices_["Laocc0A"] = std::make_shared<Matrix>("Laocc0A", nn, na);
        matrices_["Ufocc0A"] = std::make_shared<Matrix>("Ufocc0A", nf, nf);
        matrices_["Uaocc0A"] = std::make_shared<Matrix>("Uaocc0A", na, na);

        double** Lp = matrices_["Locc0A"]->pointer();
        double** Lfp = matrices_["Lfocc0A"]->pointer();
        double** Lap = matrices_["Laocc0A"]->pointer();
        double** Up = matrices_["Uocc0A"]->pointer();
        double** Ufp = matrices_["Ufocc0A"]->pointer();
        double** Uap = matrices_["Uaocc0A"]->pointer();

        for (int n = 0; n < nn; n++) {
            for (int i = 0; i < nf; i++) {
                Lfp[n][i] = Lp[n][i];
            }
            for (int i = 0; i < na; i++) {
                Lap[n][i] = Lp[n][i + nf];
            }
        }

        for (int i = 0; i < nf; i++) {
            for (int j = 0; j < nf; j++) {
                Ufp[i][j] = Up[i][j];
            }
        }

        for (int i = 0; i < na; i++) {
            for (int j = 0; j < na; j++) {
                Uap[i][j] = Up[i + nf][j + nf];
            }
        }

        matrices_["Locc0A"]->set_name("Locc0A");
        matrices_["Lfocc0A"]->set_name("Lfocc0A");
        matrices_["Laocc0A"]->set_name("Laocc0A");
        matrices_["Uocc0A"]->set_name("Uocc0A");
        matrices_["Ufocc0A"]->set_name("Ufocc0A");
        matrices_["Uaocc0A"]->set_name("Uaocc0A");
        matrices_["Qocc0A"]->set_name("Qocc0A");

        if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") {
// Extra matrices of localized orbitals augmented by the link orbital (normalized to 1/2, not 1)
            matrices_["LLocc0A"] = std::make_shared<Matrix>("LLocc0A", nn, nm+1);
            double** LLp = matrices_["LLocc0A"]->pointer();
            for (int n = 0; n < nn; n++) {
                for (int i = 0; i < nm; i++) {
                    LLp[n][i] = Lp[n][i];
                }
            }
            double** thislinkAp = matrices_["thislinkA"]->pointer();
            for (int n = 0; n < nn; n++) {
                LLp[n][nm] = thislinkAp[n][0];
            }
        }
    }

    {
        outfile->Printf("  Local Orbitals for Monomer B:\n\n");

        int nn = matrices_["Caocc0B"]->rowspi()[0];
        int nf = matrices_["Cfocc0B"]->colspi()[0];
        int na = matrices_["Caocc0B"]->colspi()[0];
        int nm = nf + na;

        std::vector<int> ranges;
        ranges.push_back(0);
        ranges.push_back(nf);
        ranges.push_back(nm);

        std::shared_ptr<Matrix> Focc(
            new Matrix("Focc", vectors_["eps_occ0B"]->dimpi()[0], vectors_["eps_occ0B"]->dimpi()[0]));
        Focc->set_diagonal(vectors_["eps_occ0B"]);

        std::shared_ptr<fisapt::IBOLocalizer2> local =
            fisapt::IBOLocalizer2::build(primary_, reference_->get_basisset("MINAO"), matrices_["Cocc0B"], options_);
        local->print_header();
        std::map<std::string, std::shared_ptr<Matrix> > ret = local->localize(matrices_["Cocc0B"], Focc, ranges);

        matrices_["Locc0B"] = ret["L"];
        matrices_["Uocc0B"] = ret["U"];
        matrices_["Qocc0B"] = ret["Q"];

        matrices_["Lfocc0B"] = std::make_shared<Matrix>("Lfocc0B", nn, nf);
        matrices_["Laocc0B"] = std::make_shared<Matrix>("Laocc0B", nn, na);
        matrices_["Ufocc0B"] = std::make_shared<Matrix>("Ufocc0B", nf, nf);
        matrices_["Uaocc0B"] = std::make_shared<Matrix>("Uaocc0B", na, na);

        double** Lp = matrices_["Locc0B"]->pointer();
        double** Lfp = matrices_["Lfocc0B"]->pointer();
        double** Lap = matrices_["Laocc0B"]->pointer();
        double** Up = matrices_["Uocc0B"]->pointer();
        double** Ufp = matrices_["Ufocc0B"]->pointer();
        double** Uap = matrices_["Uaocc0B"]->pointer();

        for (int n = 0; n < nn; n++) {
            for (int i = 0; i < nf; i++) {
                Lfp[n][i] = Lp[n][i];
            }
            for (int i = 0; i < na; i++) {
                Lap[n][i] = Lp[n][i + nf];
            }
        }

        for (int i = 0; i < nf; i++) {
            for (int j = 0; j < nf; j++) {
                Ufp[i][j] = Up[i][j];
            }
        }

        for (int i = 0; i < na; i++) {
            for (int j = 0; j < na; j++) {
                Uap[i][j] = Up[i + nf][j + nf];
            }
        }

        matrices_["Locc0B"]->set_name("Locc0B");
        matrices_["Lfocc0B"]->set_name("Lfocc0B");
        matrices_["Laocc0B"]->set_name("Laocc0B");
        matrices_["Uocc0B"]->set_name("Uocc0B");
        matrices_["Ufocc0B"]->set_name("Ufocc0B");
        matrices_["Uaocc0B"]->set_name("Uaocc0B");
        matrices_["Qocc0B"]->set_name("Qocc0B");

        if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") {
// Extra matrices of localized orbitals augmented by the link orbital (normalized to 1/2, not 1)
            matrices_["LLocc0B"] = std::make_shared<Matrix>("LLocc0B", nn, nm+1);
            double** LLp = matrices_["LLocc0B"]->pointer();
            for (int n = 0; n < nn; n++) {
                for (int i = 0; i < nm; i++) {
                    LLp[n][i] = Lp[n][i];
                }
            }
            double** thislinkBp = matrices_["thislinkB"]->pointer();
            for (int n = 0; n < nn; n++) {
                LLp[n][nm] = thislinkBp[n][0];
            }
        }
    }
}

// Compute fragment-fragment partitioning of electrostatic contribution
void FISAPT::felst() {
    outfile->Printf("  ==> F-SAPT Electrostatics <==\n\n");

    // => Sizing <= //

    std::string link_assignment = options_.get_str("FISAPT_LINK_ASSIGNMENT");
    std::shared_ptr<Molecule> mol = primary_->molecule();
    int nn = primary_->nbf();
    int nA = mol->natom();
    int nB = mol->natom();

// Some special considerations for the SAOn/SIAOn variants of ISAPT below
    int na = matrices_["Locc0A"]->colspi()[0];
    int nb = matrices_["Locc0B"]->colspi()[0];
    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") {
        na = matrices_["LLocc0A"]->colspi()[0];
        nb = matrices_["LLocc0B"]->colspi()[0];
    }
    std::shared_ptr<Matrix> L0A = std::make_shared<Matrix>("L0A", nn, na);
    std::shared_ptr<Matrix> L0B = std::make_shared<Matrix>("L0B", nn, nb);
    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") {
        L0A->copy(matrices_["LLocc0A"]);
        L0B->copy(matrices_["LLocc0B"]);
    }
    else {
        L0A->copy(matrices_["Locc0A"]);
        L0B->copy(matrices_["Locc0B"]);
    }
// for the SAOn/SIAOn variants, we need to do the same things but with the Locc0A, Locc0B matrices augmented by the link orbital. 
// The local matrices L0A/L0B hold the correct form, augmented or not, as needed.

    // => Targets <= //

    double Elst10 = 0.0;
    std::vector<double> Elst10_terms;
    Elst10_terms.resize(4);

    // This matrix contains the decomposition of the F-SAPT electrostatics energy to nuclear, orbital, and external
    // point charges contributions. nA and nB are the total number of atoms in subsystems A and B. na(nb) is the number
    // of occupied oritals for subsystem A(B). One additional entry is added to include the effect of external 
    // point charges in A and B (if present; otherwise it will be zero columns/rows). This matrix is saved
    // in the Elst.dat file generated after the F-SAPT analysis. The fsapt.py script reads this file and the
    // user-defined functional group partitioning and analyzes the F-SAPT interaction energy in terms of functional
    // group contributions.
    //
    // For the nuclei-nuclei interactions (entries [0:nA-1, 0:nB-1]), the total number of A-B atoms is used for the rows
    // and columns but only the entries where nuclei A interact with nuclei B are actually nonzero. Entries
    // [nA:nA+na-1, nB:nB+nb-1] represent the interactions between the local occupied orbitals in A and the local occupied
    // orbitals in B. Lastly, entry [nA, nB] represent the interaction between the external point charges in A and B. 
    // Cross terms represent the interaction between nuclei and electrons, nuclei and external point charges, 
    // and electrons and external point charges.
    //
    // Similar matrices are used for the other SAPT components.
    matrices_["Elst_AB"] = std::make_shared<Matrix>("Elst_AB", nA + na + 1, nB + nb + 1); // Add one entry for external
                                                                                          // potentials in A and B
    double** Ep = matrices_["Elst_AB"]->pointer();

    // => A <-> B <= //

    double* ZAp = vectors_["ZA"]->pointer();
    double* ZBp = vectors_["ZB"]->pointer();
    for (int A = 0; A < nA; A++) {
        for (int B = 0; B < nB; B++) {
            if (A == B) continue;
            double E = ZAp[A] * ZBp[B] / mol->xyz(A).distance(mol->xyz(B));
            Ep[A][B] += E;
            Elst10_terms[3] += E;
        }
    }

    // Extern A-atom B interaction
    // Compute the interaction between the external potenrial and each atom in the fragment
    if (reference_->has_potential_variable("A")) {
        double conv = 1;
        if (mol->units() == Molecule::Angstrom) {
            conv *= pc_bohr2angstroms;
        }
        double E = 0.0;
        for (int B = 0; B < nB; B++) {
            auto atom = std::make_shared<Molecule>();
            atom->set_units(mol->units());
            atom->add_atom(ZBp[B], mol->x(B)*conv, mol->y(B)*conv, mol->z(B)*conv);
            double interaction = reference_->potential_variable("A")->computeNuclearEnergy(atom);
            Ep[nA + na][B] = interaction;
            E += interaction;
        }
        Elst10_terms[3] += E;
    }

    // Extern B-atom A interaction
    // Compute the interaction between the external potenrial and each atom in the fragment
    if (reference_->has_potential_variable("B")) {
        double conv = 1;
        if (mol->units() == Molecule::Angstrom) {
            conv *= pc_bohr2angstroms;
        }
        double E = 0.0;
        for (int A = 0; A < nA; A++) {
            auto atom = std::make_shared<Molecule>();
            atom->set_units(mol->units());
            atom->add_atom(ZAp[A], mol->x(A)*conv, mol->y(A)*conv, mol->z(A)*conv);
            double interaction = reference_->potential_variable("B")->computeNuclearEnergy(atom);
            Ep[A][nB + nb] = interaction;
            E += interaction;
        }
        Elst10_terms[3] += E;
    }

    // => a <-> b <= //

    int nT = 1;
#ifdef _OPENMP
    nT = Process::environment.get_n_threads();
#endif

    // => Get integrals from DFHelper <= //
    dfh_ = std::make_shared<DFHelper>(primary_, reference_->get_basisset("DF_BASIS_SCF"));
    dfh_->set_memory(doubles_);
    dfh_->set_method("DIRECT_iaQ");
    dfh_->set_nthreads(nT);
    dfh_->initialize();
    dfh_->print_header();

    dfh_->add_space("a", L0A);
    dfh_->add_space("b", L0B);

    dfh_->add_transformation("Aaa", "a", "a");
    dfh_->add_transformation("Abb", "b", "b");

    dfh_->transform();

    size_t nQ = dfh_->get_naux();
    auto QaC = std::make_shared<Matrix>("QaC", na, nQ);
    double** QaCp = QaC->pointer();
    for (size_t a = 0; a < na; a++) {
        dfh_->fill_tensor("Aaa", QaCp[a], {a, a + 1}, {a, a + 1});
    }

    auto QbC = std::make_shared<Matrix>("QbC", nb, nQ);
    double** QbCp = QbC->pointer();
    for (size_t b = 0; b < nb; b++) {
        dfh_->fill_tensor("Abb", QbCp[b], {b, b + 1}, {b, b + 1});
    }

    std::shared_ptr<Matrix> Elst10_3 = linalg::doublet(QaC, QbC, false, true);
    double** Elst10_3p = Elst10_3->pointer();
    for (int a = 0; a < na; a++) {
        for (int b = 0; b < nb; b++) {
            double E = 4.0 * Elst10_3p[a][b];
            Elst10_terms[2] += E;
            Ep[a + nA][b + nB] += E;
        }
    }

    matrices_["Vlocc0A"] = QaC;
    matrices_["Vlocc0B"] = QbC;

    // => Nuclear Part (PITA) <= //
    auto Vfact2 = std::make_shared<IntegralFactory>(primary_);
    std::shared_ptr<PotentialInt> Vint2(static_cast<PotentialInt*>(Vfact2->ao_potential().release()));
    auto Vtemp2 = std::make_shared<Matrix>("Vtemp2", nn, nn);


    // => A <-> b <= //

    for (int A = 0; A < nA; A++) {
        if (ZAp[A] == 0.0) continue;
        Vtemp2->zero();
        Vint2->set_charge_field({{ZAp[A], {mol->x(A), mol->y(A), mol->z(A)}}});
        Vint2->compute(Vtemp2);
        std::shared_ptr<Matrix> Vbb =
            linalg::triplet(L0B, Vtemp2, L0B, true, false, false);
        double** Vbbp = Vbb->pointer();
        for (int b = 0; b < nb; b++) {
            double E = 2.0 * Vbbp[b][b];
            Elst10_terms[1] += E;
            Ep[A][b + nB] += E;
        }
    }

    // Add Extern-A - Orbital b interaction
    if (reference_->has_potential_variable("A")) {
        std::shared_ptr<Matrix> Vbb =
            linalg::triplet(L0B, matrices_["VA_extern"], L0B, true, false, false);
        double** Vbbp = Vbb->pointer();
        for (int b = 0; b < nb; b++) {
            double E = 2.0 * Vbbp[b][b];
            Elst10_terms[1] += E;
            Ep[nA + na][b + nB] += E;
        }
    }

    // => a <-> B <= //

    for (int B = 0; B < nB; B++) {
        if (ZBp[B] == 0.0) continue;
        Vtemp2->zero();
        Vint2->set_charge_field({{ZBp[B], {mol->x(B), mol->y(B), mol->z(B)}}});
        Vint2->compute(Vtemp2);
        std::shared_ptr<Matrix> Vaa =
            linalg::triplet(L0A, Vtemp2, L0A, true, false, false);
        double** Vaap = Vaa->pointer();
        for (int a = 0; a < na; a++) {
            double E = 2.0 * Vaap[a][a];
            Elst10_terms[0] += E;
            Ep[a + nA][B] += E;
        }
    }

    // Add Extern-B - Orbital a interaction
    if (reference_->has_potential_variable("B")) {
        std::shared_ptr<Matrix> Vaa =
            linalg::triplet(L0A, matrices_["VB_extern"], L0A, true, false, false);
        double** Vaap = Vaa->pointer();
        for (int a = 0; a < na; a++) {
            double E = 2.0 * Vaap[a][a];
            Elst10_terms[0] += E;
            Ep[a + nA][nB + nb] += E;
        }
    }


    // Prepare DFHelper object for the next module
    dfh_->clear_spaces();

    // => Summation <= //

    for (int k = 0; k < Elst10_terms.size(); k++) {
        Elst10 += Elst10_terms[k];
    }
    // for (int k = 0; k < Elst10_terms.size(); k++) {
    //    outfile->Printf("    Elst10,r (%1d)        = %18.12lf [Eh]\n",k+1,Elst10_terms[k]);
    //}
    // scalars_["Elst10,r"] = Elst10;
    outfile->Printf("    Elst10,r            = %18.12lf [Eh]\n", Elst10);

    // Add extern-extern contribution
    if (reference_->has_potential_variable("A") && reference_->has_potential_variable("B")) {
        // Add the interaction between the external potenitals in A and B
        // Multiply by 2 to get the full A-B + B-A interaction energy
        double ext = matrices_["extern_extern_IE"]->get(0, 1)*2.0; 
        Ep[nA + na][nB + nb] += ext;
        outfile->Printf("    Extern-Extern       = %18.12lf [Eh]\n", ext);
    }

    outfile->Printf("\n");

    // fflush(outfile);
}

// Compute fragment-fragment partitioning of exchange contribution
void FISAPT::fexch() {
    outfile->Printf("  ==> F-SAPT Exchange <==\n\n");

    // => Sizing <= //

    std::shared_ptr<Molecule> mol = primary_->molecule();
    int nn = primary_->nbf();
    int nA = mol->natom();
    int nB = mol->natom();
    int na = matrices_["Locc0A"]->colspi()[0];
    int nb = matrices_["Locc0B"]->colspi()[0];
    int nr = matrices_["Cvir0A"]->colspi()[0];
    int ns = matrices_["Cvir0B"]->colspi()[0];

//for the SAOn/SIAOn variants, we sometimes need na1 = na+1 (with link orbital) and sometimes na (without) - be careful with this!
    std::string link_assignment = options_.get_str("FISAPT_LINK_ASSIGNMENT");
    int na1 = na;
    int nb1 = nb;
    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") {
        na1 = na + 1;
        nb1 = nb + 1;
    }

    // => Targets <= //

    double Exch10_2 = 0.0;
    std::vector<double> Exch10_2_terms;
    Exch10_2_terms.resize(3);

    matrices_["Exch_AB"] = std::make_shared<Matrix>("Exch_AB", nA + na1 + 1, nB + nb1 + 1); // Add one entry for external
                                                                                          // potentials in A and B
    double** Ep = matrices_["Exch_AB"]->pointer();

    // ==> Stack Variables <== //

    std::shared_ptr<Matrix> S = matrices_["S"];
    std::shared_ptr<Matrix> V_A = matrices_["V_A"];
    std::shared_ptr<Matrix> J_A = matrices_["J_A"];
    std::shared_ptr<Matrix> V_B = matrices_["V_B"];
    std::shared_ptr<Matrix> J_B = matrices_["J_B"];

    std::shared_ptr<Matrix> LoccA = matrices_["Locc0A"];
    std::shared_ptr<Matrix> LoccB = matrices_["Locc0B"];
    std::shared_ptr<Matrix> CvirA = matrices_["Cvir0A"];
    std::shared_ptr<Matrix> CvirB = matrices_["Cvir0B"];

    // ==> DF ERI Setup (JKFIT Type, in Full Basis) <== //

    int nT = 1;
#ifdef _OPENMP
    nT = Process::environment.get_n_threads();
#endif

    std::vector<std::shared_ptr<Matrix> > Cs;
    Cs.push_back(LoccA);
    Cs.push_back(CvirA);
    Cs.push_back(LoccB);
    Cs.push_back(CvirB);

    size_t max_MO = 0;
    for (auto& mat : Cs) max_MO = std::max(max_MO, (size_t)mat->ncol());

    dfh_->add_space("a", Cs[0]);
    dfh_->add_space("r", Cs[1]);
    dfh_->add_space("b", Cs[2]);
    dfh_->add_space("s", Cs[3]);

    dfh_->add_transformation("Aar", "a", "r");
    dfh_->add_transformation("Abs", "b", "s");

    dfh_->transform();

    // ==> Electrostatic Potentials <== //

    std::shared_ptr<Matrix> W_A(J_A->clone());
    W_A->copy(J_A);
    W_A->scale(2.0);
    W_A->add(V_A);

    std::shared_ptr<Matrix> W_B(J_B->clone());
    W_B->copy(J_B);
    W_B->scale(2.0);
    W_B->add(V_B);

    std::shared_ptr<Matrix> WAbs = linalg::triplet(LoccB, W_A, CvirB, true, false, false);
    std::shared_ptr<Matrix> WBar = linalg::triplet(LoccA, W_B, CvirA, true, false, false);
    double** WBarp = WBar->pointer();
    double** WAbsp = WAbs->pointer();

    W_A.reset();
    W_B.reset();

    // ==> Exchange S^2 Computation <== //

    std::shared_ptr<Matrix> Sab = linalg::triplet(LoccA, S, LoccB, true, false, false);
    std::shared_ptr<Matrix> Sba = linalg::triplet(LoccB, S, LoccA, true, false, false);
    std::shared_ptr<Matrix> Sas = linalg::triplet(LoccA, S, CvirB, true, false, false);
    std::shared_ptr<Matrix> Sbr = linalg::triplet(LoccB, S, CvirA, true, false, false);
    double** Sabp = Sab->pointer();
    double** Sbap = Sba->pointer();
    double** Sasp = Sas->pointer();
    double** Sbrp = Sbr->pointer();

    auto WBab = std::make_shared<Matrix>("WBab", na, nb);
    double** WBabp = WBab->pointer();
    auto WAba = std::make_shared<Matrix>("WAba", nb, na);
    double** WAbap = WAba->pointer();

    C_DGEMM('N', 'T', na, nb, nr, 1.0, WBarp[0], nr, Sbrp[0], nr, 0.0, WBabp[0], nb);
    C_DGEMM('N', 'T', nb, na, ns, 1.0, WAbsp[0], ns, Sasp[0], ns, 0.0, WAbap[0], na);

    auto E_exch1 = std::make_shared<Matrix>("E_exch [a <x- b]", na, nb);
    double** E_exch1p = E_exch1->pointer();
    auto E_exch2 = std::make_shared<Matrix>("E_exch [a -x> b]", na, nb);
    double** E_exch2p = E_exch2->pointer();

    for (int a = 0; a < na; a++) {
        for (int b = 0; b < nb; b++) {
            E_exch1p[a][b] -= 2.0 * Sabp[a][b] * WBabp[a][b];
            E_exch2p[a][b] -= 2.0 * Sbap[b][a] * WAbap[b][a];
        }
    }

    // E_exch1->print();
    // E_exch2->print();

    size_t nQ = dfh_->get_naux();
    auto TrQ = std::make_shared<Matrix>("TrQ", nr, nQ);
    double** TrQp = TrQ->pointer();
    auto TsQ = std::make_shared<Matrix>("TsQ", ns, nQ);
    double** TsQp = TsQ->pointer();
    auto TbQ = std::make_shared<Matrix>("TbQ", nb, nQ);
    double** TbQp = TbQ->pointer();
    auto TaQ = std::make_shared<Matrix>("TaQ", na, nQ);
    double** TaQp = TaQ->pointer();

    dfh_->add_disk_tensor("Bab", std::make_tuple(na, nb, nQ));

    for (size_t a = 0; a < na; a++) {
        dfh_->fill_tensor("Aar", TrQ, {a, a + 1});
        C_DGEMM('N', 'N', nb, nQ, nr, 1.0, Sbrp[0], nr, TrQp[0], nQ, 0.0, TbQp[0], nQ);
        dfh_->write_disk_tensor("Bab", TbQ, {a, a + 1});
    }

    dfh_->add_disk_tensor("Bba", std::make_tuple(nb, na, nQ));

    for (size_t b = 0; b < nb; b++) {
        dfh_->fill_tensor("Abs", TsQ, {b, b + 1});
        C_DGEMM('N', 'N', na, nQ, ns, 1.0, Sasp[0], ns, TsQp[0], nQ, 0.0, TaQp[0], nQ);
        dfh_->write_disk_tensor("Bba", TaQ, {b, b + 1});
    }

    auto E_exch3 = std::make_shared<Matrix>("E_exch [a <x-x> b]", na, nb);
    double** E_exch3p = E_exch3->pointer();

    for (size_t a = 0; a < na; a++) {
        dfh_->fill_tensor("Bab", TbQ, {a, a + 1});
        for (size_t b = 0; b < nb; b++) {
            dfh_->fill_tensor("Bba", TaQ, {b, b + 1}, {a, a + 1});
            E_exch3p[a][b] -= 2.0 * C_DDOT(nQ, TbQp[b], 1, TaQp[0], 1);
        }
    }

    // E_exch3->print();

    // => Totals <= //

    for (int a = 0; a < na; a++) {
        for (int b = 0; b < nb; b++) {
            Ep[a + nA][b + nB] = E_exch1p[a][b] + E_exch2p[a][b] + E_exch3p[a][b];
            Exch10_2_terms[0] += E_exch1p[a][b];
            Exch10_2_terms[1] += E_exch2p[a][b];
            Exch10_2_terms[2] += E_exch3p[a][b];
        }
    }

    for (int k = 0; k < Exch10_2_terms.size(); k++) {
        Exch10_2 += Exch10_2_terms[k];
    }
    // for (int k = 0; k < Exch10_2_terms.size(); k++) {
    //    outfile->Printf("    Exch10(S^2) (%1d)     = %18.12lf [Eh]\n",k+1,Exch10_2_terms[k]);
    //}
    // scalars_["Exch10(S^2)"] = Exch10_2;
    outfile->Printf("    Exch10(S^2)         = %18.12lf [Eh]\n", Exch10_2);
    outfile->Printf("\n");
    // fflush(outfile);

    // => Exchange scaling <= //

    if (options_.get_bool("FISAPT_FSAPT_EXCH_SCALE")) {
    // KP: the following change in the scaling (see the commented out line) is essential for ISAPT(SAOn) and ISAPT(SIAOn),
    // and should not mess up the other variants. The problem is that exch() calculates E(10)exch(S^2) in the DM formalism
    // (good for SAO/SIAO because no virtual orbitals in formulas, bad for FSAPT), while fexch() calculates this term
    // in the second quantization formalism (good for FSAPT, bad for SAO/SIAO). In effect, the scaling for SAO/SIAO accounts
    // for two effects at once: terms beyond S^2 and linking terms computed by exch() but not here.
    //
    // Not changing the scaling for sSAPT0 - if you try using sSAPT0 with ISAPT(SAO/SIAO), you've only got yourself to blame.
    //  double scale = scalars_["Exch10"] / scalars_["Exch10(S^2)"];
        double scale = scalars_["Exch10"] / Exch10_2;
        matrices_["Exch_AB"]->scale(scale);
        outfile->Printf("    Scaling F-SAPT Exch10(S^2) by %11.3E to match Exch10\n\n", scale);
    }
    if (options_.get_bool("SSAPT0_SCALE")) {
        sSAPT0_scale_ = scalars_["Exch10"] / scalars_["Exch10(S^2)"];
        sSAPT0_scale_ = pow(sSAPT0_scale_, 3.0);
        outfile->Printf("    Scaling F-SAPT Exch-Ind and Exch-Disp by %11.3E \n\n", sSAPT0_scale_);
    }

    // Prepare DFHelper object for the next module
    dfh_->clear_spaces();
}

// Compute fragment-fragment partitioning of induction contribution
void FISAPT::find() {
    outfile->Printf("  ==> F-SAPT Induction <==\n\n");

    // => Options <= //

    bool ind_resp = options_.get_bool("FISAPT_FSAPT_IND_RESPONSE");
    bool ind_scale = options_.get_bool("FISAPT_FSAPT_IND_SCALE");
    std::string link_assignment = options_.get_str("FISAPT_LINK_ASSIGNMENT");

    // => Sizing <= //

    std::shared_ptr<Molecule> mol = primary_->molecule();
    int nn = primary_->nbf();
    int nA = mol->natom();
    int nB = mol->natom();
    int na = matrices_["Locc0A"]->colspi()[0];
    int nb = matrices_["Locc0B"]->colspi()[0];
    int nr = matrices_["Cvir0A"]->colspi()[0];
    int ns = matrices_["Cvir0B"]->colspi()[0];

//for the SAOn/SIAOn variants, we sometimes need na1 = na+1 (with link orbital) and sometimes na (without) - be careful with this!
    int na1 = na;
    int nb1 = nb;
    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") {
        na1 = na + 1;
        nb1 = nb + 1;
    }

    // => Pointers <= //

    std::shared_ptr<Matrix> Locc_A = matrices_["Locc0A"];
    std::shared_ptr<Matrix> Locc_B = matrices_["Locc0B"];

    std::shared_ptr<Matrix> Uocc_A = matrices_["Uocc0A"];
    std::shared_ptr<Matrix> Uocc_B = matrices_["Uocc0B"];

    std::shared_ptr<Matrix> Cocc_A = matrices_["Cocc0A"];
    std::shared_ptr<Matrix> Cocc_B = matrices_["Cocc0B"];
    std::shared_ptr<Matrix> Cvir_A = matrices_["Cvir0A"];
    std::shared_ptr<Matrix> Cvir_B = matrices_["Cvir0B"];

    std::shared_ptr<Vector> eps_occ_A = vectors_["eps_occ0A"];
    std::shared_ptr<Vector> eps_occ_B = vectors_["eps_occ0B"];
    std::shared_ptr<Vector> eps_vir_A = vectors_["eps_vir0A"];
    std::shared_ptr<Vector> eps_vir_B = vectors_["eps_vir0B"];

    // => DFHelper = DF + disk tensors <= //

    int nT = 1;
#ifdef _OPENMP
    nT = Process::environment.get_n_threads();
#endif

    size_t nQ = dfh_->get_naux();

    // => ESPs <= //

    dfh_->add_disk_tensor("WBar", std::make_tuple(nB + nb1 + 1, na, nr)); // add entry for external potentials
    dfh_->add_disk_tensor("WAbs", std::make_tuple(nA + na1 + 1, nb, ns)); // add entry for external potentials

    // => Nuclear Part (PITA) <= //

    auto Vfact2 = std::make_shared<IntegralFactory>(primary_);
    std::shared_ptr<PotentialInt> Vint2(static_cast<PotentialInt*>(Vfact2->ao_potential().release()));
    auto Vtemp2 = std::make_shared<Matrix>("Vtemp2", nn, nn);

    double* ZAp = vectors_["ZA"]->pointer();
    for (size_t A = 0; A < nA; A++) {
        Vtemp2->zero();
        Vint2->set_charge_field({{ZAp[A], {mol->x(A), mol->y(A), mol->z(A)}}});
        Vint2->compute(Vtemp2);
        std::shared_ptr<Matrix> Vbs = linalg::triplet(Cocc_B, Vtemp2, Cvir_B, true, false, false);
        dfh_->write_disk_tensor("WAbs", Vbs, {A, A + 1});
    }

    double* ZBp = vectors_["ZB"]->pointer();
    for (size_t B = 0; B < nB; B++) {
        Vtemp2->zero();
        Vint2->set_charge_field({{ZBp[B], {mol->x(B), mol->y(B), mol->z(B)}}});
        Vint2->compute(Vtemp2);
        std::shared_ptr<Matrix> Var = linalg::triplet(Cocc_A, Vtemp2, Cvir_A, true, false, false);
        dfh_->write_disk_tensor("WBar", Var, {B, B + 1});
    }

    // ==> DFHelper Setup (JKFIT Type, in Full Basis) <== //

    std::vector<std::shared_ptr<Matrix> > Cs;
    Cs.push_back(Cocc_A);
    Cs.push_back(Cvir_A);
    Cs.push_back(Cocc_B);
    Cs.push_back(Cvir_B);

    size_t max_MO = 0;
    for (auto& mat : Cs) max_MO = std::max(max_MO, (size_t)mat->ncol());

    dfh_->add_space("a", Cs[0]);
    dfh_->add_space("r", Cs[1]);
    dfh_->add_space("b", Cs[2]);
    dfh_->add_space("s", Cs[3]);

    dfh_->add_transformation("Aar", "a", "r");
    dfh_->add_transformation("Abs", "b", "s");

    dfh_->transform();

    // => Electronic Part (Massive PITA) <= //

    double** RaCp = matrices_["Vlocc0A"]->pointer();
    double** RbDp = matrices_["Vlocc0B"]->pointer();
    // for SAOn/SIAOn, the above two matrices contain the link orbital contribution

    auto TsQ = std::make_shared<Matrix>("TsQ", ns, nQ);
    auto T1As = std::make_shared<Matrix>("T1As", na1, ns);
    double** TsQp = TsQ->pointer();
    double** T1Asp = T1As->pointer();
    for (size_t b = 0; b < nb; b++) {
        dfh_->fill_tensor("Abs", TsQ, {b, b + 1});
        C_DGEMM('N', 'T', na1, ns, nQ, 2.0, RaCp[0], nQ, TsQp[0], nQ, 0.0, T1Asp[0], ns);
        for (size_t a = 0; a < na1; a++) {
            dfh_->write_disk_tensor("WAbs", T1Asp[a], {nA + a, nA + a + 1}, {b, b + 1});
        }
    }

    auto TrQ = std::make_shared<Matrix>("TrQ", nr, nQ);
    auto T1Br = std::make_shared<Matrix>("T1Br", nb1, nr);
    double** TrQp = TrQ->pointer();
    double** T1Brp = T1Br->pointer();
    for (size_t a = 0; a < na; a++) {
        dfh_->fill_tensor("Aar", TrQ, {a, a + 1});
        C_DGEMM('N', 'T', nb1, nr, nQ, 2.0, RbDp[0], nQ, TrQp[0], nQ, 0.0, T1Brp[0], nr);
        for (size_t b = 0; b < nb1; b++) {
            dfh_->write_disk_tensor("WBar", T1Brp[b], {nB + b, nB + b + 1}, {a, a + 1});
        }
    }

    // ==> Stack Variables <== //

    double* eap = eps_occ_A->pointer();
    double* ebp = eps_occ_B->pointer();
    double* erp = eps_vir_A->pointer();
    double* esp = eps_vir_B->pointer();

    std::shared_ptr<Matrix> S = matrices_["S"];
    std::shared_ptr<Matrix> D_A = matrices_["D_A"];
    std::shared_ptr<Matrix> V_A = matrices_["V_A"];
    std::shared_ptr<Matrix> J_A = matrices_["J_A"];
    std::shared_ptr<Matrix> K_A = matrices_["K_A"];
    std::shared_ptr<Matrix> D_B = matrices_["D_B"];
    std::shared_ptr<Matrix> V_B = matrices_["V_B"];
    std::shared_ptr<Matrix> J_B = matrices_["J_B"];
    std::shared_ptr<Matrix> K_B = matrices_["K_B"];
    std::shared_ptr<Matrix> J_O = matrices_["J_O"];
    std::shared_ptr<Matrix> K_O = matrices_["K_O"];
    std::shared_ptr<Matrix> J_P_A = matrices_["J_P_A"];
    std::shared_ptr<Matrix> J_P_B = matrices_["J_P_B"];

    // ==> MO Amplitudes/Sources (by source atom) <== //

    auto xA = std::make_shared<Matrix>("xA", na, nr);
    auto xB = std::make_shared<Matrix>("xB", nb, ns);
    double** xAp = xA->pointer();
    double** xBp = xB->pointer();

    auto wB = std::make_shared<Matrix>("wB", na, nr);
    auto wA = std::make_shared<Matrix>("wA", nb, ns);
    double** wBp = wB->pointer();
    double** wAp = wA->pointer();

    auto uAT = std::make_shared<Matrix>("uAT", nb, ns);
    auto wAT = std::make_shared<Matrix>("wAT", nb, ns);
    auto uBT = std::make_shared<Matrix>("uBT", na, nr);
    auto wBT = std::make_shared<Matrix>("wBT", na, nr);

// Now come the second-order exchange-induction energy expressions appropriate for the SAOn/SIAOn link assignments.
// Here, we compute only the averaged spin coupling, like when FISAPT_EXCH_PARPERP == false.
// See https://doi.org/10.1021/acs.jpca.2c06465 for details - equation numbers below refer to this paper and its SI.
    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2" ) {

//here we need the link-only parts of density matrices D_X,D_Y and their corresponding Coulomb and exchange matrices (computed earlier)
        std::shared_ptr<Matrix> D_X(D_A->clone());
        std::shared_ptr<Matrix> D_Y(D_A->clone());
        D_X = linalg::doublet(matrices_["thislinkA"], matrices_["thislinkA"], false, true);
        D_Y = linalg::doublet(matrices_["thislinkB"], matrices_["thislinkB"], false, true);
        auto J_X(matrices_["JLA"]->clone());
        auto K_X(matrices_["KLA"]->clone());
        auto J_Y(matrices_["JLB"]->clone());
        auto K_Y(matrices_["KLB"]->clone());
       
        // ==> Generalized ESP (Flat and Exchange) <== //
        std::shared_ptr<Matrix> K_AOY = matrices_["K_AOY"];
        std::shared_ptr<Matrix> K_XOB = matrices_["K_XOB"];
        K_XOB->transpose_this();
        std::shared_ptr<Matrix> J_P_YAY = matrices_["J_P_YAY"];
        std::shared_ptr<Matrix> J_P_XBX = matrices_["J_P_XBX"];
       
        std::map<std::string, std::shared_ptr<Matrix> > mapA;
        mapA["Cocc_A"] = Locc_A;
        mapA["Cvir_A"] = Cvir_A;
        mapA["S"] = S;
        mapA["D_A"] = D_A;
        mapA["V_A"] = V_A;
        mapA["J_A"] = J_A;
        mapA["K_A"] = K_A;
        mapA["D_B"] = D_B;
        mapA["V_B"] = V_B;
        mapA["J_B"] = J_B;
        mapA["K_B"] = K_B;
        mapA["D_X"] = D_X;
        mapA["J_X"] = J_X;
        mapA["K_X"] = K_X;
        mapA["D_Y"] = D_Y;
        mapA["J_Y"] = J_Y;
        mapA["K_Y"] = K_Y;
        mapA["J_O"] = J_O;
        mapA["K_O"] = K_O;
        mapA["K_AOY"] = K_AOY;
        mapA["J_P"] = J_P_A;
        mapA["J_PYAY"] = J_P_YAY;

        wBT = build_ind_pot(mapA);
        uBT = build_exch_ind_pot_avg(mapA);

        K_O->transpose_this();

        std::map<std::string, std::shared_ptr<Matrix> > mapB;
        mapB["Cocc_A"] = Locc_B;
        mapB["Cvir_A"] = Cvir_B;
        mapB["S"] = S;
        mapB["D_A"] = D_B;
        mapB["V_A"] = V_B;
        mapB["J_A"] = J_B;
        mapB["K_A"] = K_B;
        mapB["D_B"] = D_A;
        mapB["V_B"] = V_A;
        mapB["J_B"] = J_A;
        mapB["K_B"] = K_A;
        mapB["D_X"] = D_Y;
        mapB["J_X"] = J_Y;
        mapB["K_X"] = K_Y;
        mapB["D_Y"] = D_X;
        mapB["J_Y"] = J_X;
        mapB["K_Y"] = K_X;
        mapB["J_O"] = J_O;
        mapB["K_O"] = K_O;
        mapB["K_AOY"] = K_XOB;
        mapB["J_P"] = J_P_B;
        mapB["J_PYAY"] = J_P_XBX;

        wAT = build_ind_pot(mapB);
        uAT = build_exch_ind_pot_avg(mapB);

        K_O->transpose_this();

    }
    else {

    // ==> Generalized ESP (Flat and Exchange) <== //

        std::map<std::string, std::shared_ptr<Matrix> > mapA;
        mapA["Cocc_A"] = Locc_A;
        mapA["Cvir_A"] = Cvir_A;
        mapA["Cocc_B"] = Locc_B;
        mapA["Cvir_B"] = Cvir_B;
        mapA["S"] = S;
        mapA["D_A"] = D_A;
        mapA["V_A"] = V_A;
        mapA["J_A"] = J_A;
        mapA["K_A"] = K_A;
        mapA["D_B"] = D_B;
        mapA["V_B"] = V_B;
        mapA["J_B"] = J_B;
        mapA["K_B"] = K_B;
        mapA["J_O"] = J_O;
        mapA["K_O"] = K_O;
        mapA["J_P"] = J_P_A;

        wBT = build_ind_pot(mapA);
        uBT = build_exch_ind_pot(mapA);

        K_O->transpose_this();

        std::map<std::string, std::shared_ptr<Matrix> > mapB;
        mapB["Cocc_A"] = Locc_B;
        mapB["Cvir_A"] = Cvir_B;
        mapB["Cocc_B"] = Locc_A;
        mapB["Cvir_B"] = Cvir_A;
        mapB["S"] = S;
        mapB["D_A"] = D_B;
        mapB["V_A"] = V_B;
        mapB["J_A"] = J_B;
        mapB["K_A"] = K_B;
        mapB["D_B"] = D_A;
        mapB["V_B"] = V_A;
        mapB["J_B"] = J_A;
        mapB["K_B"] = K_A;
        mapB["J_O"] = J_O;
        mapB["K_O"] = K_O;
        mapB["J_P"] = J_P_B;

        wAT = build_ind_pot(mapB);
        uAT = build_exch_ind_pot(mapB);

        K_O->transpose_this();
    }
    double** wATp = wAT->pointer();
    double** uATp = uAT->pointer();
    double** wBTp = wBT->pointer();
    double** uBTp = uBT->pointer();

    // ==> Uncoupled Targets <== //

    auto Ind20u_AB_terms = std::make_shared<Matrix>("Ind20 [A<-B] (a x B)", na, nB + nb1 + 1); // add one entry for external potential 
    auto Ind20u_BA_terms = std::make_shared<Matrix>("Ind20 [B<-A] (A x b)", nA + na1 + 1, nb); // add one entry for external potential
    double** Ind20u_AB_termsp = Ind20u_AB_terms->pointer();
    double** Ind20u_BA_termsp = Ind20u_BA_terms->pointer();

    double Ind20u_AB = 0.0;
    double Ind20u_BA = 0.0;

    auto ExchInd20u_AB_terms = std::make_shared<Matrix>("ExchInd20 [A<-B] (a x B)", na, nB + nb1 + 1); // add one for external potential 
    auto ExchInd20u_BA_terms = std::make_shared<Matrix>("ExchInd20 [B<-A] (A x b)", nA + na1 + 1, nb); // add one for external potential
    double** ExchInd20u_AB_termsp = ExchInd20u_AB_terms->pointer();
    double** ExchInd20u_BA_termsp = ExchInd20u_BA_terms->pointer();

    double ExchInd20u_AB = 0.0;
    double ExchInd20u_BA = 0.0;

    int sna = 0;
    int snB = 0;
    int snb = 0;
    int snA = 0;

    if (options_.get_bool("SSAPT0_SCALE")) {
    // This will NOT work with ISAPT(SAOn/SIAOn) and I (Konrad) think that's OK.
        sna = na;
        snB = nB;
        snb = nb;
        snA = nA;
    }

    std::shared_ptr<Matrix> sExchInd20u_AB_terms =
        std::make_shared<Matrix>("sExchInd20 [A<-B] (a x B)", sna, snB + snb + 1); // add one entry for external potential
    std::shared_ptr<Matrix> sExchInd20u_BA_terms =
        std::make_shared<Matrix>("sExchInd20 [B<-A] (A x b)", snA + sna + 1, snb); // add one entry for external potential
    double** sExchInd20u_AB_termsp = sExchInd20u_AB_terms->pointer();
    double** sExchInd20u_BA_termsp = sExchInd20u_BA_terms->pointer();

    double sExchInd20u_AB = 0.0;
    double sExchInd20u_BA = 0.0;

    auto Indu_AB_terms = std::make_shared<Matrix>("Ind [A<-B] (a x B)", na, nB + nb1 + 1); // add one entry for external potential
    auto Indu_BA_terms = std::make_shared<Matrix>("Ind [B<-A] (A x b)", nA + na1 + 1, nb); // add one entry for external potential
    double** Indu_AB_termsp = Indu_AB_terms->pointer();
    double** Indu_BA_termsp = Indu_BA_terms->pointer();

    double Indu_AB = 0.0;
    double Indu_BA = 0.0;

    auto sIndu_AB_terms = std::make_shared<Matrix>("sInd [A<-B] (a x B)", sna, snB + snb + 1); // add one entry for external potential 
    auto sIndu_BA_terms = std::make_shared<Matrix>("sInd [B<-A] (A x b)", snA + sna + 1, snb); // add one entry for external potential
    double** sIndu_AB_termsp = sIndu_AB_terms->pointer();
    double** sIndu_BA_termsp = sIndu_BA_terms->pointer();

    double sIndu_AB = 0.0;
    double sIndu_BA = 0.0;

    // ==> A <- B Uncoupled <== //

    // Add the external potential
    if (reference_->has_potential_variable("B")) {
        std::shared_ptr<Matrix> Var = linalg::triplet(Cocc_A, matrices_["VB_extern"], Cvir_A, true, false, false);
        dfh_->write_disk_tensor("WBar", Var, {(size_t) nB + nb1, (size_t) nB + nb1 + 1});
    }

    else { // Add empty matrix
        std::shared_ptr<Matrix> Var = std::make_shared<Matrix>("zero", na, nr);
        Var->zero();
        dfh_->write_disk_tensor("WBar", Var, {(size_t) nB + nb1, (size_t) nB + nb1 + 1});
    }

    for (size_t B = 0; B < nB + nb1 + 1; B++) { // add one for external potential
        // ESP
        dfh_->fill_tensor("WBar", wB, {B, B + 1});

        // Uncoupled amplitude
        for (int a = 0; a < na; a++) {
            for (int r = 0; r < nr; r++) {
                xAp[a][r] = wBp[a][r] / (eap[a] - erp[r]);
            }
        }

        // Backtransform the amplitude to LO
        std::shared_ptr<Matrix> x2A = linalg::doublet(Uocc_A, xA, true, false);
        double** x2Ap = x2A->pointer();

        // Zip up the Ind20 contributions
        for (int a = 0; a < na; a++) {
            double Jval = 2.0 * C_DDOT(nr, x2Ap[a], 1, wBTp[a], 1);
            double Kval = 2.0 * C_DDOT(nr, x2Ap[a], 1, uBTp[a], 1);
            Ind20u_AB_termsp[a][B] = Jval;
            Ind20u_AB += Jval;
            ExchInd20u_AB_termsp[a][B] = Kval;
            ExchInd20u_AB += Kval;
            if (options_.get_bool("SSAPT0_SCALE")) {
                sExchInd20u_AB_termsp[a][B] = Kval;
                sExchInd20u_AB += Kval;
                sIndu_AB_termsp[a][B] = Jval + Kval;
                sIndu_AB += Jval + Kval;
            }

            Indu_AB_termsp[a][B] = Jval + Kval;
            Indu_AB += Jval + Kval;
        }
    }

    // ==> B <- A Uncoupled <== //

    // Add the external potential
    if (reference_->has_potential_variable("A")) {
        std::shared_ptr<Matrix> Vbs = linalg::triplet(Cocc_B, matrices_["VA_extern"], Cvir_B, true, false, false);
        dfh_->write_disk_tensor("WAbs", Vbs, {(size_t) nA + na1, (size_t) nA + na1 + 1});
    }

    else { // add empty matrix
        std::shared_ptr<Matrix> Vbs = std::make_shared<Matrix>("zero", nb, ns);
        Vbs->zero();
        dfh_->write_disk_tensor("WAbs", Vbs, {(size_t) nA + na1, (size_t) nA + na1 + 1});
    }


    for (size_t A = 0; A < nA + na1 + 1; A++) { // add one for extenral potential
        // ESP
        dfh_->fill_tensor("WAbs", wA, {A, A + 1});

        // Uncoupled amplitude
        for (int b = 0; b < nb; b++) {
            for (int s = 0; s < ns; s++) {
                xBp[b][s] = wAp[b][s] / (ebp[b] - esp[s]);
            }
        }

        // Backtransform the amplitude to LO
        std::shared_ptr<Matrix> x2B = linalg::doublet(Uocc_B, xB, true, false);
        double** x2Bp = x2B->pointer();

        // Zip up the Ind20 contributions
        for (int b = 0; b < nb; b++) {
            double Jval = 2.0 * C_DDOT(ns, x2Bp[b], 1, wATp[b], 1);
            double Kval = 2.0 * C_DDOT(ns, x2Bp[b], 1, uATp[b], 1);
            Ind20u_BA_termsp[A][b] = Jval;
            Ind20u_BA += Jval;
            ExchInd20u_BA_termsp[A][b] = Kval;
            ExchInd20u_BA += Kval;
            if (options_.get_bool("SSAPT0_SCALE")) {
                sExchInd20u_BA_termsp[A][b] = Kval;
                sExchInd20u_BA += Kval;
                sIndu_BA_termsp[A][b] = Jval + Kval;
                sIndu_BA += Jval + Kval;
            }
            Indu_BA_termsp[A][b] = Jval + Kval;
            Indu_BA += Jval + Kval;
        }
    }

    double Ind20u = Ind20u_AB + Ind20u_BA;
    outfile->Printf("    Ind20,u (A<-B)      = %18.12lf [Eh]\n", Ind20u_AB);
    outfile->Printf("    Ind20,u (B<-A)      = %18.12lf [Eh]\n", Ind20u_BA);
    outfile->Printf("    Ind20,u             = %18.12lf [Eh]\n", Ind20u);
    // fflush(outfile);

    double ExchInd20u = ExchInd20u_AB + ExchInd20u_BA;
    outfile->Printf("    Exch-Ind20,u (A<-B) = %18.12lf [Eh]\n", ExchInd20u_AB);
    outfile->Printf("    Exch-Ind20,u (B<-A) = %18.12lf [Eh]\n", ExchInd20u_BA);
    outfile->Printf("    Exch-Ind20,u        = %18.12lf [Eh]\n", ExchInd20u);
    outfile->Printf("\n");
    // fflush(outfile);
    if (options_.get_bool("SSAPT0_SCALE")) {
        double sExchInd20u = sExchInd20u_AB + sExchInd20u_BA;
        outfile->Printf("    sExch-Ind20,u (A<-B) = %18.12lf [Eh]\n", sExchInd20u_AB);
        outfile->Printf("    sExch-Ind20,u (B<-A) = %18.12lf [Eh]\n", sExchInd20u_BA);
        outfile->Printf("    sExch-Ind20,u        = %18.12lf [Eh]\n", sExchInd20u);
        outfile->Printf("\n");
    }

    double Ind = Ind20u + ExchInd20u;
    std::shared_ptr<Matrix> Ind_AB_terms = Indu_AB_terms;
    std::shared_ptr<Matrix> Ind_BA_terms = Indu_BA_terms;
    std::shared_ptr<Matrix> sInd_AB_terms = sIndu_AB_terms;
    std::shared_ptr<Matrix> sInd_BA_terms = sIndu_BA_terms;

    if (ind_resp) {
        outfile->Printf("  COUPLED INDUCTION (You asked for it!):\n\n");

        // ==> Coupled Targets <== //

        auto Ind20r_AB_terms = std::make_shared<Matrix>("Ind20 [A<-B] (a x B)", na, nB + nb1 + 1); // add one for external potential 
        auto Ind20r_BA_terms = std::make_shared<Matrix>("Ind20 [B<-A] (A x b)", nA + na1 + 1, nb); // add one for external potential
        double** Ind20r_AB_termsp = Ind20r_AB_terms->pointer();
        double** Ind20r_BA_termsp = Ind20r_BA_terms->pointer();

        double Ind20r_AB = 0.0;
        double Ind20r_BA = 0.0;

        auto ExchInd20r_AB_terms = std::make_shared<Matrix>("ExchInd20 [A<-B] (a x B)", na, nB + nb1 + 1); // add one for external
                                                                                                          // potential
        auto ExchInd20r_BA_terms = std::make_shared<Matrix>("ExchInd20 [B<-A] (A x b)", nA + na1 + 1, nb); // add one for external
                                                                                                          // potential
        double** ExchInd20r_AB_termsp = ExchInd20r_AB_terms->pointer();
        double** ExchInd20r_BA_termsp = ExchInd20r_BA_terms->pointer();

        double ExchInd20r_AB = 0.0;
        double ExchInd20r_BA = 0.0;

        auto Indr_AB_terms = std::make_shared<Matrix>("Ind [A<-B] (a x B)", na, nB + nb1 + 1); // add one entry for external potential 
        auto Indr_BA_terms = std::make_shared<Matrix>("Ind [B<-A] (A x b)", nA + na1 + 1, nb); // add one entry for external potential
        double** Indr_AB_termsp = Indr_AB_terms->pointer();
        double** Indr_BA_termsp = Indr_BA_terms->pointer();

        double Indr_AB = 0.0;
        double Indr_BA = 0.0;

        // => JK Object <= //

        // TODO: Account for 2-index overhead in memory
        auto nso = primary_->nbf();
        auto jk_memory = (long int)doubles_;
        jk_memory -= 24 * nso * nso;
        jk_memory -= 4 * na * nso;
        jk_memory -= 4 * nb * nso;
        if (jk_memory < 0L) {
            throw PSIEXCEPTION("Too little static memory for FISAPT::induction");
        }

        std::shared_ptr<JK> jk =
            JK::build_JK(primary_, reference_->get_basisset("DF_BASIS_SCF"), options_, false, (size_t)jk_memory);

        jk->set_memory((size_t)jk_memory);
        jk->set_do_J(true);
        jk->set_do_K(true);
        jk->initialize();
        jk->print_header();

        // ==> Master Loop over perturbing atoms <== //

        int nC = std::max(nA + na1 + 1, nB + nb1 + 1); // add one for external potential

        for (size_t C = 0; C < nC; C++) {
            if (C < nB + nb1) dfh_->fill_tensor("WBar", wB, {C, C + 1});
            if (C < nA + na1) dfh_->fill_tensor("WAbs", wB, {C, C + 1});

            outfile->Printf("    Responses for (A <- Source B = %3zu) and (B <- Source A = %3zu)\n\n",
                            (C < nB + nb1 ? C : nB + nb1), (C < nA + na1 ? C : nA + na1));

            auto cphf = std::make_shared<CPHF_FISAPT>();

            // Effective constructor
            cphf->delta_ = options_.get_double("CPHF_R_CONVERGENCE");
            cphf->maxiter_ = options_.get_int("MAXITER");
            cphf->jk_ = jk;

            cphf->w_A_ = wB;  // Reversal of convention
            cphf->Cocc_A_ = Cocc_A;
            cphf->Cvir_A_ = Cvir_A;
            cphf->eps_occ_A_ = eps_occ_A;
            cphf->eps_vir_A_ = eps_vir_A;

            cphf->w_B_ = wA;  // Reversal of convention
            cphf->Cocc_B_ = Cocc_B;
            cphf->Cvir_B_ = Cvir_B;
            cphf->eps_occ_B_ = eps_occ_B;
            cphf->eps_vir_B_ = eps_vir_B;

            // Gogo CPKS
            cphf->compute_cphf();

            xA = cphf->x_A_;
            xB = cphf->x_B_;

            xA->scale(-1.0);
            xB->scale(-1.0);

            if (C < nB + nb1 + 1) {
                // Backtransform the amplitude to LO
                std::shared_ptr<Matrix> x2A = linalg::doublet(Uocc_A, xA, true, false);
                double** x2Ap = x2A->pointer();

                // Zip up the Ind20 contributions
                for (int a = 0; a < na; a++) {
                    double Jval = 2.0 * C_DDOT(nr, x2Ap[a], 1, wBTp[a], 1);
                    double Kval = 2.0 * C_DDOT(nr, x2Ap[a], 1, uBTp[a], 1);
                    Ind20r_AB_termsp[a][C] = Jval;
                    Ind20r_AB += Jval;
                    ExchInd20r_AB_termsp[a][C] = Kval;
                    ExchInd20r_AB += Kval;
                    Indr_AB_termsp[a][C] = Jval + Kval;
                    Indr_AB += Jval + Kval;
                }
            }

            if (C < nA + na1 + 1) {
                // Backtransform the amplitude to LO
                std::shared_ptr<Matrix> x2B = linalg::doublet(Uocc_B, xB, true, false);
                double** x2Bp = x2B->pointer();

                // Zip up the Ind20 contributions
                for (int b = 0; b < nb; b++) {
                    double Jval = 2.0 * C_DDOT(ns, x2Bp[b], 1, wATp[b], 1);
                    double Kval = 2.0 * C_DDOT(ns, x2Bp[b], 1, uATp[b], 1);
                    Ind20r_BA_termsp[C][b] = Jval;
                    Ind20r_BA += Jval;
                    ExchInd20r_BA_termsp[C][b] = Kval;
                    ExchInd20r_BA += Kval;
                    Indr_BA_termsp[C][b] = Jval + Kval;
                    Indr_BA += Jval + Kval;
                }
            }
        }

        double Ind20r = Ind20r_AB + Ind20r_BA;
        outfile->Printf("    Ind20,r (A<-B)      = %18.12lf [Eh]\n", Ind20r_AB);
        outfile->Printf("    Ind20,r (B<-A)      = %18.12lf [Eh]\n", Ind20r_BA);
        outfile->Printf("    Ind20,r             = %18.12lf [Eh]\n", Ind20r);
        // fflush(outfile);

        double ExchInd20r = ExchInd20r_AB + ExchInd20r_BA;
        outfile->Printf("    Exch-Ind20,r (A<-B) = %18.12lf [Eh]\n", ExchInd20r_AB);
        outfile->Printf("    Exch-Ind20,r (B<-A) = %18.12lf [Eh]\n", ExchInd20r_BA);
        outfile->Printf("    Exch-Ind20,r        = %18.12lf [Eh]\n", ExchInd20r);
        outfile->Printf("\n");
        // fflush(outfile);

        Ind = Ind20r + ExchInd20r;
        Ind_AB_terms = Indr_AB_terms;
        Ind_BA_terms = Indr_BA_terms;
    }

    // => Induction scaling <= //

    if (ind_scale) {
        double dHF = 0.0;
        if (scalars_["HF"] != 0.0) {
            dHF = scalars_["HF"] - scalars_["Elst10,r"] - scalars_["Exch10"] - scalars_["Ind20,r"] -
                  scalars_["Exch-Ind20,r"];
        }
        double IndHF = scalars_["Ind20,r"] + scalars_["Exch-Ind20,r"] + dHF;
        double IndSAPT0 = scalars_["Ind20,r"] + scalars_["Exch-Ind20,r"];

        double Sdelta = IndHF / IndSAPT0;
        double SrAB = (ind_resp ? 1.0
                                : (scalars_["Ind20,r (A<-B)"] + scalars_["Exch-Ind20,r (A<-B)"]) /
                                      (scalars_["Ind20,u (A<-B)"] + scalars_["Exch-Ind20,u (A<-B)"]));
        double SrBA = (ind_resp ? 1.0
                                : (scalars_["Ind20,r (B<-A)"] + scalars_["Exch-Ind20,r (B<-A)"]) /
                                      (scalars_["Ind20,u (B<-A)"] + scalars_["Exch-Ind20,u (B<-A)"]));

        double sIndHF = scalars_["Ind20,r"] + scalars_["sExch-Ind20,r"] + dHF;
        double sIndSAPT0 = scalars_["Ind20,r"] + scalars_["sExch-Ind20,r"];

        double sSdelta = sIndHF / IndSAPT0;

        double sSrAB = (ind_resp ? 1.0
                                 : (scalars_["Ind20,r (A<-B)"] + scalars_["sExch-Ind20,r (A<-B)"]) /
                                       (scalars_["Ind20,u (A<-B)"] + scalars_["sExch-Ind20,u (A<-B)"]));
        double sSrBA = (ind_resp ? 1.0
                                 : (scalars_["Ind20,r (B<-A)"] + scalars_["sExch-Ind20,r (B<-A)"]) /
                                       (scalars_["Ind20,u (B<-A)"] + scalars_["sExch-Ind20,u (B<-A)"]));

        outfile->Printf("    Scaling for delta HF        = %11.3E\n", Sdelta);
        outfile->Printf("    Scaling for response (A<-B) = %11.3E\n", SrAB);
        outfile->Printf("    Scaling for response (B<-A) = %11.3E\n", SrBA);
        outfile->Printf("    Scaling for total (A<-B)    = %11.3E\n", Sdelta * SrAB);
        outfile->Printf("    Scaling for total (B<-A)    = %11.3E\n", Sdelta * SrBA);
        outfile->Printf("\n");

        Ind_AB_terms->scale(Sdelta * SrAB);
        Ind_BA_terms->scale(Sdelta * SrBA);
        Ind20u_AB_terms->scale(Sdelta * SrAB);
        ExchInd20u_AB_terms->scale(Sdelta * SrAB);
        Ind20u_BA_terms->scale(Sdelta * SrBA);
        ExchInd20u_BA_terms->scale(Sdelta * SrBA);
        sInd_AB_terms->scale(sSdelta * SrAB);
        sInd_BA_terms->scale(sSdelta * SrBA);
    }

    matrices_["IndAB_AB"] = std::make_shared<Matrix>("IndAB_AB", nA + na1 + 1, nB + nb1 + 1); // add one entry for external
                                                                                            // potentials in A and B
    matrices_["IndBA_AB"] = std::make_shared<Matrix>("IndBA_AB", nA + na1 + 1, nB + nb1 + 1); // add one entry for external
                                                                                            // potentials in A and B
    matrices_["Ind20u_AB_terms"] = std::make_shared<Matrix>("Ind20uAB_AB", nA + na1 + 1, nB + nb1 + 1); // add one for external pot
    matrices_["ExchInd20u_AB_terms"] = std::make_shared<Matrix>("ExchInd20uAB_AB", nA + na1 + 1, nB + nb1 + 1); // add one for external
    matrices_["Ind20u_BA_terms"] = std::make_shared<Matrix>("Ind20uBA_AB", nA + na1 + 1, nB + nb1 + 1); // add one for external
    matrices_["ExchInd20u_BA_terms"] = std::make_shared<Matrix>("ExchInd20uBA_AB", nA + na1 + 1, nB + nb1 + 1); // add one for external
    double** EABp = matrices_["IndAB_AB"]->pointer();
    double** EBAp = matrices_["IndBA_AB"]->pointer();
    double** Ind20ABp = matrices_["Ind20u_AB_terms"]->pointer();
    double** ExchInd20ABp = matrices_["ExchInd20u_AB_terms"]->pointer();
    double** Ind20BAp = matrices_["Ind20u_BA_terms"]->pointer();
    double** ExchInd20BAp = matrices_["ExchInd20u_BA_terms"]->pointer();
    double** EAB2p = Ind_AB_terms->pointer();
    double** EBA2p = Ind_BA_terms->pointer();
    double** Ind20AB2p = Ind20u_AB_terms->pointer();
    double** ExchInd20AB2p = ExchInd20u_AB_terms->pointer();
    double** Ind20BA2p = Ind20u_BA_terms->pointer();
    double** ExchInd20BA2p = ExchInd20u_BA_terms->pointer();

    for (int a = 0; a < na; a++) {
        for (int B = 0; B < nB + nb1 + 1; B++) { // add one for external potential
            EABp[a + nA][B] = EAB2p[a][B];
            Ind20ABp[a + nA][B] = Ind20AB2p[a][B];
            ExchInd20ABp[a + nA][B] = ExchInd20AB2p[a][B];
        }
    }

    for (int A = 0; A < nA + na1 + 1; A++) { // add one for external potential
        for (int b = 0; b < nb; b++) {
            EBAp[A][b + nB] = EBA2p[A][b];
            Ind20BAp[A][b + nB] = Ind20BA2p[A][b];
            ExchInd20BAp[A][b + nB] = ExchInd20BA2p[A][b];
        }
    }

    matrices_["sIndAB_AB"] = std::make_shared<Matrix>("sIndAB_AB", snA + sna + 1, snB + snb + 1); // add one entry for external
                                                                                                  // potentials in A and B
    matrices_["sIndBA_AB"] = std::make_shared<Matrix>("sIndBA_AB", snA + sna + 1, snB + snb + 1); // add one entry for external
                                                                                                  // potentials in A and B
    double** sEABp = matrices_["sIndAB_AB"]->pointer();
    double** sEBAp = matrices_["sIndBA_AB"]->pointer();
    double** sEAB2p = sInd_AB_terms->pointer();
    double** sEBA2p = sInd_BA_terms->pointer();

    for (int a = 0; a < sna; a++) {
        for (int B = 0; B < snB + snb + 1; B++) { // add one for external potential
            sEABp[a + snA][B] = sEAB2p[a][B];
        }
    }

    for (int A = 0; A < snA + sna + 1; A++) { // add one for external potential
        for (int b = 0; b < snb; b++) {
            sEBAp[A][b + snB] = sEBA2p[A][b];
        }
    }
    // We're done with dfh_'s integrals
    dfh_->clear_all();
}

// Compute fragment-fragment partitioning of dispersion contribution
// This is also where the modified exchange-dispersion expressions for the SAOn/SIAOn ISAPT algorithms are coded.
void FISAPT::fdisp() {
    outfile->Printf("  ==> F-SAPT Dispersion <==\n\n");

    // => Auxiliary Basis Set <= //

    std::shared_ptr<BasisSet> auxiliary = reference_->get_basisset("DF_BASIS_SAPT");

    // => Sizing <= //

    std::shared_ptr<Molecule> mol = primary_->molecule();
    int nn = primary_->nbf();
    int nA = mol->natom();
    int nB = mol->natom();
    int na = matrices_["Laocc0A"]->colspi()[0];
    int nb = matrices_["Laocc0B"]->colspi()[0];
    int nr = matrices_["Cvir0A"]->colspi()[0];
    int ns = matrices_["Cvir0B"]->colspi()[0];
    int nQ = auxiliary->nbf();
    size_t naQ = na * (size_t)nQ;
    size_t nbQ = nb * (size_t)nQ;

//for the SAOn/SIAOn variants, we sometimes need na1 = na+1 (with link orbital) and sometimes na (without) - be careful with this!
    std::string link_assignment = options_.get_str("FISAPT_LINK_ASSIGNMENT");
    int na1 = na;
    int nb1 = nb;
    if (link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") {
        na1 = na + 1;
        nb1 = nb + 1;
    }

    int nfa = matrices_["Lfocc0A"]->colspi()[0];
    int nfb = matrices_["Lfocc0B"]->colspi()[0];

    int nT = 1;
#ifdef _OPENMP
    nT = Process::environment.get_n_threads();
#endif

    // => Targets <= //

    matrices_["Disp_AB"] = std::make_shared<Matrix>("Disp_AB", nA + nfa + na1 + 1, nB + nfb + nb1 + 1); // add one entry for external
                                                                                                      // potentials in A and B
    double** Ep = matrices_["Disp_AB"]->pointer();

    int snA = 0;
    int snfa = 0;
    int sna = 0;
    int snB = 0;
    int snfb = 0;
    int snb = 0;

    if (options_.get_bool("SSAPT0_SCALE")) {
        snA = nA;
        snfa = nfa;
        sna = na;
        snB = nB;
        snfb = nfb;
        snb = nb;
    }

    matrices_["sDisp_AB"] = std::make_shared<Matrix>("Disp_AB", snA + snfa + sna + 1, snB + snfb + snb + 1); // add one entry for
                                                                                                             // external potenrials
                                                                                                             // in A and B
    double** sEp = matrices_["sDisp_AB"]->pointer();

    // => Stashed Variables <= //

    std::shared_ptr<Matrix> S = matrices_["S"];
    std::shared_ptr<Matrix> D_A = matrices_["D_A"];
    std::shared_ptr<Matrix> P_A = matrices_["P_A"];
    std::shared_ptr<Matrix> V_A = matrices_["V_A"];
    std::shared_ptr<Matrix> J_A = matrices_["J_A"];
    std::shared_ptr<Matrix> K_A = matrices_["K_A"];
    std::shared_ptr<Matrix> D_B = matrices_["D_B"];
    std::shared_ptr<Matrix> P_B = matrices_["P_B"];
    std::shared_ptr<Matrix> V_B = matrices_["V_B"];
    std::shared_ptr<Matrix> J_B = matrices_["J_B"];
    std::shared_ptr<Matrix> K_B = matrices_["K_B"];
    std::shared_ptr<Matrix> K_O = matrices_["K_O"];

    bool parperp = options_.get_bool("FISAPT_EXCH_PARPERP");

    std::shared_ptr<Matrix> Caocc_A = matrices_["Caocc0A"];
    std::shared_ptr<Matrix> Caocc_B = matrices_["Caocc0B"];
    std::shared_ptr<Matrix> Cavir_A = matrices_["Cvir0A"];
    std::shared_ptr<Matrix> Cavir_B = matrices_["Cvir0B"];

    std::shared_ptr<Vector> eps_aocc_A = vectors_["eps_aocc0A"];
    std::shared_ptr<Vector> eps_aocc_B = vectors_["eps_aocc0B"];
    std::shared_ptr<Vector> eps_avir_A = vectors_["eps_vir0A"];
    std::shared_ptr<Vector> eps_avir_B = vectors_["eps_vir0B"];

    std::shared_ptr<Matrix> Uaocc_A = matrices_["Uaocc0A"];
    std::shared_ptr<Matrix> Uaocc_B = matrices_["Uaocc0B"];

    // => Auxiliary C matrices <= //

    std::shared_ptr<Matrix> Cr1 = linalg::triplet(D_B, S, Cavir_A);
    Cr1->scale(-1.0);
    Cr1->add(Cavir_A);
    std::shared_ptr<Matrix> Cs1 = linalg::triplet(D_A, S, Cavir_B);
    Cs1->scale(-1.0);
    Cs1->add(Cavir_B);
    std::shared_ptr<Matrix> Ca2 = linalg::triplet(D_B, S, Caocc_A);
    std::shared_ptr<Matrix> Cb2 = linalg::triplet(D_A, S, Caocc_B);

    std::shared_ptr<Matrix> Cr3 = linalg::triplet(D_B, S, Cavir_A);
    std::shared_ptr<Matrix> CrX = linalg::triplet(linalg::triplet(D_A, S, D_B), S, Cavir_A);
    Cr3->subtract(CrX);
    Cr3->scale(2.0);
    std::shared_ptr<Matrix> Cs3 = linalg::triplet(D_A, S, Cavir_B);
    std::shared_ptr<Matrix> CsX = linalg::triplet(linalg::triplet(D_B, S, D_A), S, Cavir_B);
    Cs3->subtract(CsX);
    Cs3->scale(2.0);

    std::shared_ptr<Matrix> Ca4 = linalg::triplet(linalg::triplet(D_A, S, D_B), S, Caocc_A);
    Ca4->scale(-2.0);
    std::shared_ptr<Matrix> Cb4 = linalg::triplet(linalg::triplet(D_B, S, D_A), S, Caocc_B);
    Cb4->scale(-2.0);

    // => Auxiliary V matrices <= //

    std::shared_ptr<Matrix> Jbr = linalg::triplet(Caocc_B, J_A, Cavir_A, true, false, false);
    Jbr->scale(2.0);
    std::shared_ptr<Matrix> Kbr = linalg::triplet(Caocc_B, K_A, Cavir_A, true, false, false);
    Kbr->scale(-1.0);

    std::shared_ptr<Matrix> Jas = linalg::triplet(Caocc_A, J_B, Cavir_B, true, false, false);
    Jas->scale(2.0);
    std::shared_ptr<Matrix> Kas = linalg::triplet(Caocc_A, K_B, Cavir_B, true, false, false);
    Kas->scale(-1.0);

    std::shared_ptr<Matrix> KOas = linalg::triplet(Caocc_A, K_O, Cavir_B, true, false, false);
    KOas->scale(1.0);
    std::shared_ptr<Matrix> KObr = linalg::triplet(Caocc_B, K_O, Cavir_A, true, true, false);
    KObr->scale(1.0);

    std::shared_ptr<Matrix> JBas = linalg::triplet(linalg::triplet(Caocc_A, S, D_B, true, false, false), J_A, Cavir_B);
    JBas->scale(-2.0);
    std::shared_ptr<Matrix> JAbr = linalg::triplet(linalg::triplet(Caocc_B, S, D_A, true, false, false), J_B, Cavir_A);
    JAbr->scale(-2.0);

    std::shared_ptr<Matrix> Jbs = linalg::triplet(Caocc_B, J_A, Cavir_B, true, false, false);
    Jbs->scale(4.0);
    std::shared_ptr<Matrix> Jar = linalg::triplet(Caocc_A, J_B, Cavir_A, true, false, false);
    Jar->scale(4.0);

    std::shared_ptr<Matrix> JAas = linalg::triplet(linalg::triplet(Caocc_A, J_B, D_A, true, false, false), S, Cavir_B);
    JAas->scale(-2.0);
    std::shared_ptr<Matrix> JBbr = linalg::triplet(linalg::triplet(Caocc_B, J_A, D_B, true, false, false), S, Cavir_A);
    JBbr->scale(-2.0);

    // Get your signs right Hesselmann!
    std::shared_ptr<Matrix> Vbs = linalg::triplet(Caocc_B, V_A, Cavir_B, true, false, false);
    Vbs->scale(2.0);
    std::shared_ptr<Matrix> Var = linalg::triplet(Caocc_A, V_B, Cavir_A, true, false, false);
    Var->scale(2.0);
    std::shared_ptr<Matrix> VBas = linalg::triplet(linalg::triplet(Caocc_A, S, D_B, true, false, false), V_A, Cavir_B);
    VBas->scale(-1.0);
    std::shared_ptr<Matrix> VAbr = linalg::triplet(linalg::triplet(Caocc_B, S, D_A, true, false, false), V_B, Cavir_A);
    VAbr->scale(-1.0);
    std::shared_ptr<Matrix> VRas = linalg::triplet(linalg::triplet(Caocc_A, V_B, P_A, true, false, false), S, Cavir_B);
    VRas->scale(1.0);
    std::shared_ptr<Matrix> VSbr = linalg::triplet(linalg::triplet(Caocc_B, V_A, P_B, true, false, false), S, Cavir_A);
    VSbr->scale(1.0);

    std::shared_ptr<Matrix> Sas = linalg::triplet(Caocc_A, S, Cavir_B, true, false, false);
    std::shared_ptr<Matrix> Sbr = linalg::triplet(Caocc_B, S, Cavir_A, true, false, false);

    std::shared_ptr<Matrix> Qbr(Jbr->clone());
    Qbr->zero();
    Qbr->add(Jbr);
    Qbr->add(Kbr);
    Qbr->add(KObr);
    Qbr->add(JAbr);
    Qbr->add(JBbr);
    Qbr->add(VAbr);
    Qbr->add(VSbr);

    std::shared_ptr<Matrix> Qas(Jas->clone());
    Qas->zero();
    Qas->add(Jas);
    Qas->add(Kas);
    Qas->add(KOas);
    Qas->add(JAas);
    Qas->add(JBas);
    Qas->add(VBas);
    Qas->add(VRas);

    std::shared_ptr<Matrix> SBar = linalg::triplet(linalg::triplet(Caocc_A, S, D_B, true, false, false), S, Cavir_A);
    std::shared_ptr<Matrix> SAbs = linalg::triplet(linalg::triplet(Caocc_B, S, D_A, true, false, false), S, Cavir_B);

    std::shared_ptr<Matrix> Qar(Jar->clone());
    Qar->zero();
    Qar->add(Jar);
    Qar->add(Var);

    std::shared_ptr<Matrix> Qbs(Jbs->clone());
    Qbs->zero();
    Qbs->add(Jbs);
    Qbs->add(Vbs);

    std::shared_ptr<Matrix> KXOYas(Jas->clone());
    std::shared_ptr<Matrix> KXOYbr(Jbr->clone());
    double** KXOYasp = KXOYas->pointer();
    double** KXOYbrp = KXOYbr->pointer();

    Jbr.reset();
    Kbr.reset();
    Jas.reset();
    Kas.reset();
    KOas.reset();
    KObr.reset();
    JBas.reset();
    JAbr.reset();
    Jbs.reset();
    Jar.reset();
    JAas.reset();
    JBbr.reset();
    Vbs.reset();
    Var.reset();
    VBas.reset();
    VAbr.reset();
    VRas.reset();
    VSbr.reset();

    // => Integrals from DFHelper <= //

    std::vector<std::shared_ptr<Matrix> > Cs;
    Cs.push_back(Caocc_A);
    Cs.push_back(Cavir_A);
    Cs.push_back(Caocc_B);
    Cs.push_back(Cavir_B);
    Cs.push_back(Cr1);
    Cs.push_back(Cs1);
    Cs.push_back(Ca2);
    Cs.push_back(Cb2);
    Cs.push_back(Cr3);
    Cs.push_back(Cs3);
    Cs.push_back(Ca4);
    Cs.push_back(Cb4);

    size_t max_MO = 0, ncol = 0;
    for (auto& mat : Cs) {
        max_MO = std::max(max_MO, (size_t)mat->ncol());
        ncol += (size_t)mat->ncol();
    }

    auto dfh(std::make_shared<DFHelper>(primary_, auxiliary));
    dfh->set_memory(doubles_ - Cs[0]->nrow() * ncol);
    dfh->set_method("DIRECT_iaQ");
    dfh->set_nthreads(nT);
    dfh->initialize();
    dfh->print_header();

    dfh->add_space("a", Cs[0]);
    dfh->add_space("r", Cs[1]);
    dfh->add_space("b", Cs[2]);
    dfh->add_space("s", Cs[3]);
    dfh->add_space("r1", Cs[4]);
    dfh->add_space("s1", Cs[5]);
    dfh->add_space("a2", Cs[6]);
    dfh->add_space("b2", Cs[7]);
    dfh->add_space("r3", Cs[8]);
    dfh->add_space("s3", Cs[9]);
    dfh->add_space("a4", Cs[10]);
    dfh->add_space("b4", Cs[11]);

    dfh->add_transformation("Aar", "r", "a");
    dfh->add_transformation("Abs", "s", "b");
    dfh->add_transformation("Bas", "s1", "a");
    dfh->add_transformation("Bbr", "r1", "b");
    dfh->add_transformation("Cas", "s", "a2");
    dfh->add_transformation("Cbr", "r", "b2");
    dfh->add_transformation("Dar", "r3", "a");
    dfh->add_transformation("Dbs", "s3", "b");
    dfh->add_transformation("Ear", "r", "a4");
    dfh->add_transformation("Ebs", "s", "b4");

    // => Additional quantities needed for parallel/perpendicular link orbital spin coupling (but not for their average)  <= //

    if ((link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") && parperp) {
        std::shared_ptr<Matrix> K_XOY = matrices_["K_XOY"];
        std::shared_ptr<Matrix> D_X = matrices_["D_X"];
        std::shared_ptr<Matrix> J_X = matrices_["J_X"];
        std::shared_ptr<Matrix> K_X = matrices_["K_X"];
        std::shared_ptr<Matrix> D_Y = matrices_["D_Y"];
        std::shared_ptr<Matrix> J_Y = matrices_["J_Y"];
        std::shared_ptr<Matrix> K_Y = matrices_["K_Y"];
        std::shared_ptr<Matrix> Cx = matrices_["thislinkA"];
        std::shared_ptr<Matrix> Cy = matrices_["thislinkB"];
        std::shared_ptr<Matrix> Cx1 = linalg::triplet(D_Y, S, Cavir_A);
        Cx1->scale(-1.0);
        std::shared_ptr<Matrix> Cy1 = linalg::triplet(D_X, S, Cavir_B);
        Cy1->scale(-1.0);
        std::shared_ptr<Matrix> Cx2 = linalg::triplet(D_Y, S, Caocc_A);
        std::shared_ptr<Matrix> Cy2 = linalg::triplet(D_X, S, Caocc_B);
        std::shared_ptr<Matrix> Cx3 = linalg::triplet(linalg::triplet(D_X, S, D_Y), S, Cavir_A);
        Cx3->scale(-2.0);
        std::shared_ptr<Matrix> Cy3 = linalg::triplet(linalg::triplet(D_Y, S, D_X), S, Cavir_B);
        Cy3->scale(-2.0);
        std::shared_ptr<Matrix> Cx4 = linalg::triplet(linalg::triplet(D_X, S, D_Y), S, Caocc_A);
        Cx4->scale(-2.0);
        std::shared_ptr<Matrix> Cy4 = linalg::triplet(linalg::triplet(D_Y, S, D_X), S, Caocc_B);
        Cy4->scale(-2.0);
        KXOYas = linalg::triplet(Caocc_A, K_XOY, Cavir_B, true, false, false);
        KXOYas->scale(1.0);
        KXOYbr = linalg::triplet(Caocc_B, K_XOY, Cavir_A, true, true, false);
        KXOYbr->scale(1.0);

        Cs.push_back(Cx1);
        Cs.push_back(Cy1);
        Cs.push_back(Cx2);
        Cs.push_back(Cy2);
        Cs.push_back(Cx3);
        Cs.push_back(Cy3);
        Cs.push_back(Cx4);
        Cs.push_back(Cy4);

        dfh->add_space("x1", Cs[12]);
        dfh->add_space("y1", Cs[13]);
        dfh->add_space("x2", Cs[14]);
        dfh->add_space("y2", Cs[15]);
        dfh->add_space("x3", Cs[16]);
        dfh->add_space("y3", Cs[17]);
        dfh->add_space("x4", Cs[18]);
        dfh->add_space("y4", Cs[19]);

        dfh->add_transformation("BYas", "y1", "a");
        dfh->add_transformation("BXbr", "x1", "b");
        dfh->add_transformation("CXas", "s", "x2");
        dfh->add_transformation("CYbr", "r", "y2");
        dfh->add_transformation("DXar", "x3", "a");
        dfh->add_transformation("DYbs", "y3", "b");
        dfh->add_transformation("EXar", "r", "x4");
        dfh->add_transformation("EYbs", "s", "y4");

// now ready for DF transformation in both cases
        dfh->transform();

        Cx1.reset();
        Cy1.reset();
        Cx2.reset();
        Cy2.reset();
        Cx3.reset();
        Cy3.reset();
        Cx4.reset();
        Cy4.reset();

    }
    else {
        dfh->transform();
 
    }

    Cr1.reset();
    Cs1.reset();
    Ca2.reset();
    Cb2.reset();
    Cr3.reset();
    Cs3.reset();
    Ca4.reset();
    Cb4.reset();
    Cs.clear();
    dfh->clear_spaces();

    // => Blocking ... figure out how big a tensor slice to handle at a time <= //

    long int overhead = 0L;
    overhead += 5L * nT * na * nb; // Tab, Vab, T2ab, V2ab, and Iab work arrays below
    if ((link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") && parperp) {
        overhead += 4L * na * ns + 4L * nb * nr + 4L * na * nr + 4L * nb * ns; // Sas, Sbr, sBar, sAbs, Qas, Qbr, Qar, Qbs
    }
    else {
        overhead += 2L * na * ns + 2L * nb * nr + 2L * na * nr + 2L * nb * ns; // Sas, Sbr, sBar, sAbs, Qas, Qbr, Qar, Qbs
    }
    // the next few matrices allocated here don't take too much room (but might if large numbers of threads)
    overhead += 2L * na * nb * (nT + 1); // E_disp20 and E_exch_disp20 thread work and final matrices
    overhead += 1L * sna * snb * (nT + 1); // sE_exch_disp20 thread work and final matrices
    overhead += 1L * (nA + nfa + na) * (nB + nfb + nb); // Disp_AB
    overhead += 1L * (snA + snfa + sna) * (snB + snfb + snb); // sDisp_AB
    // account for a few of the smaller matrices already defined, but not exhaustively
    overhead += 12L * nn * nn; // D, V, J, K, P, and C matrices for A and B (neglecting C)
    long int rem = doubles_ - overhead;

    outfile->Printf("    %ld doubles - %ld overhead leaves %ld for dispersion\n", doubles_, overhead, rem);
    
    if (rem < 0L) {
        throw PSIEXCEPTION("Too little static memory for DFTSAPT::mp2_terms");
    }

    // cost_r is how much mem for Aar, Bbr, Cbr, Dar for a single r
    // cost_s would be the same value, and is the mem requirement for Abs, Bas, Cas, and Dbs for single s
    long int cost_r = 0L;
    if ((link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") && parperp) {
        cost_r = 4L * na * nQ + 4L * nb * nQ; 
    }
    else {
        cost_r = 2L * na * nQ + 2L * nb * nQ;
    } 
    long int max_r_l = rem / (2L * cost_r); // 2 b/c need to hold both an r and an s
    long int max_s_l = max_r_l;
    int max_r = (max_r_l > nr ? nr : (int) max_r_l);
    int max_s = (max_s_l > ns ? ns : (int) max_s_l);
    if (max_r < 1 || max_s < 1) {
        throw PSIEXCEPTION("Too little dynamic memory for DFTSAPT::mp2_terms");
    }
    int nrblocks = (nr / max_r);
    if (nr % max_r) nrblocks++;
    int nsblocks = (ns / max_s);
    if (ns % max_s) nsblocks++;
    outfile->Printf("    Processing a single (r,s) pair requires %ld doubles\n", cost_r * 2L);
    outfile->Printf("    %d values of r processed in %d blocks of %d\n", nr, nrblocks, max_r);
    outfile->Printf("    %d values of s processed in %d blocks of %d\n\n", ns, nsblocks, max_s);

    // => Tensor Slices <= //

    auto Aar = std::make_shared<Matrix>("Aar", max_r * na, nQ);
    auto Abs = std::make_shared<Matrix>("Abs", max_s * nb, nQ);
    auto Bas = std::make_shared<Matrix>("Bas", max_s * na, nQ);
    auto Bbr = std::make_shared<Matrix>("Bbr", max_r * nb, nQ);
    auto Cas = std::make_shared<Matrix>("Cas", max_s * na, nQ);
    auto Cbr = std::make_shared<Matrix>("Cbr", max_r * nb, nQ);
    auto Dar = std::make_shared<Matrix>("Dar", max_r * na, nQ);
    auto Dbs = std::make_shared<Matrix>("Dbs", max_s * nb, nQ);

    // => Thread Work Arrays <= //

    std::vector<std::shared_ptr<Matrix> > Tab;
    std::vector<std::shared_ptr<Matrix> > Vab;
    std::vector<std::shared_ptr<Matrix> > T2ab;
    std::vector<std::shared_ptr<Matrix> > V2ab;
    std::vector<std::shared_ptr<Matrix> > Iab;
    for (int t = 0; t < nT; t++) {
        Tab.push_back(std::make_shared<Matrix>("Tab", na, nb));
        Vab.push_back(std::make_shared<Matrix>("Vab", na, nb));
        T2ab.push_back(std::make_shared<Matrix>("T2ab", na, nb));
        V2ab.push_back(std::make_shared<Matrix>("V2ab", na, nb));
        Iab.push_back(std::make_shared<Matrix>("Iab", na, nb));
    }

    // => Pointers <= //

    double** Aarp = Aar->pointer();
    double** Absp = Abs->pointer();
    double** Basp = Bas->pointer();
    double** Bbrp = Bbr->pointer();
    double** Casp = Cas->pointer();
    double** Cbrp = Cbr->pointer();
    double** Darp = Dar->pointer();
    double** Dbsp = Dbs->pointer();

    double** Sasp = Sas->pointer();
    double** Sbrp = Sbr->pointer();
    double** SBarp = SBar->pointer();
    double** SAbsp = SAbs->pointer();

    double** Qasp = Qas->pointer();
    double** Qbrp = Qbr->pointer();
    double** Qarp = Qar->pointer();
    double** Qbsp = Qbs->pointer();

    double* eap = eps_aocc_A->pointer();
    double* ebp = eps_aocc_B->pointer();
    double* erp = eps_avir_A->pointer();
    double* esp = eps_avir_B->pointer();

    // => Slice D + E -> D <= //

    dfh->add_disk_tensor("Far", std::make_tuple(nr, na, nQ));

    for (size_t rstart = 0; rstart < nr; rstart += max_r) {
        size_t nrblock = (rstart + max_r >= nr ? nr - rstart : max_r);

        dfh->fill_tensor("Dar", Dar, {rstart, rstart + nrblock});
        dfh->fill_tensor("Ear", Aar, {rstart, rstart + nrblock});

        double* D2p = Darp[0];
        double* A2p = Aarp[0];
        for (long int arQ = 0L; arQ < nrblock * naQ; arQ++) {
            (*D2p++) += (*A2p++);
        }
        dfh->write_disk_tensor("Far", Dar, {rstart, rstart + nrblock});
    }

    dfh->add_disk_tensor("Fbs", std::make_tuple(ns, nb, nQ));

    for (size_t sstart = 0; sstart < ns; sstart += max_s) {
        size_t nsblock = (sstart + max_s >= ns ? ns - sstart : max_s);

        dfh->fill_tensor("Dbs", Dbs, {sstart, sstart + nsblock});
        dfh->fill_tensor("Ebs", Abs, {sstart, sstart + nsblock});

        double* D2p = Dbsp[0];
        double* A2p = Absp[0];
        for (long int bsQ = 0L; bsQ < nsblock * nbQ; bsQ++) {
            (*D2p++) += (*A2p++);
        }
        dfh->write_disk_tensor("Fbs", Dbs, {sstart, sstart + nsblock});
    }

    // => Targets <= //

    double Disp20 = 0.0;
    double ExchDisp20 = 0.0;
    double sExchDisp20 = 0.0;
    double par_ExchDisp20 = 0.0;

    // => Local Targets <= //

    std::vector<std::shared_ptr<Matrix> > E_disp20_threads;
    std::vector<std::shared_ptr<Matrix> > E_exch_disp20_threads;
    std::vector<std::shared_ptr<Matrix> > sE_exch_disp20_threads;
    std::vector<std::shared_ptr<Matrix> > par_E_exch_disp20_threads;
    for (int t = 0; t < nT; t++) {
        E_disp20_threads.push_back(std::make_shared<Matrix>("E_disp20", na, nb));
        E_exch_disp20_threads.push_back(std::make_shared<Matrix>("E_exch_disp20", na, nb));
        sE_exch_disp20_threads.push_back(std::make_shared<Matrix>("sE_exch_disp20", sna, snb));
        par_E_exch_disp20_threads.push_back(std::make_shared<Matrix>("par_E_exch_disp20", na, nb));
    }

    // => MO => LO Transform <= //

    double** UAp = Uaocc_A->pointer();
    double** UBp = Uaocc_B->pointer();

    double scale = 1.0;
    if (options_.get_bool("SSAPT0_SCALE")) {
        scale = sSAPT0_scale_;
    }

    // ==> Master Loop <== //

// For the SAOn/SIAOn algorithms if parallel and perpendicular spin couplings are requested for E(20)exch-disp, 
// we need to compute the expression from https://doi.org/10.1021/acs.jpca.2c06465, Eq. (8) in the SI.
// If only averaged spin coupling is requested, the standard SAPT0 formula works and we don't do anything extra.
    if ((link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") && parperp) {

        auto BYas = std::make_shared<Matrix>("BYas", max_s * na, nQ);
        auto BXbr = std::make_shared<Matrix>("BXbr", max_r * nb, nQ);
        auto CXas = std::make_shared<Matrix>("CXas", max_s * na, nQ);
        auto CYbr = std::make_shared<Matrix>("CYbr", max_r * nb, nQ);
        auto DXar = std::make_shared<Matrix>("DXar", max_r * na, nQ);
        auto DYbs = std::make_shared<Matrix>("DYbs", max_s * nb, nQ);
        double** BYasp = BYas->pointer();
        double** BXbrp = BXbr->pointer();
        double** CXasp = CXas->pointer();
        double** CYbrp = CYbr->pointer();
        double** DXarp = DXar->pointer();
        double** DYbsp = DYbs->pointer();
        double** KXOYasp = KXOYas->pointer();
        double** KXOYbrp = KXOYbr->pointer();

        dfh->add_disk_tensor("FXar", std::make_tuple(nr, na, nQ));

        for (size_t rstart = 0; rstart < nr; rstart += max_r) {
            size_t nrblock = (rstart + max_r >= nr ? nr - rstart : max_r);

            dfh->fill_tensor("DXar", DXar, {rstart, rstart + nrblock});
            dfh->fill_tensor("EXar", Aar, {rstart, rstart + nrblock});

            double* DX2p = DXarp[0];
            double* A2p = Aarp[0];
            for (long int arQ = 0L; arQ < nrblock * naQ; arQ++) {
                (*DX2p++) += (*A2p++);
            }
            dfh->write_disk_tensor("FXar", DXar, {rstart, rstart + nrblock});
        }

        dfh->add_disk_tensor("FYbs", std::make_tuple(ns, nb, nQ));

        for (size_t sstart = 0; sstart < ns; sstart += max_s) {
            size_t nsblock = (sstart + max_s >= ns ? ns - sstart : max_s);

            dfh->fill_tensor("DYbs", DYbs, {sstart, sstart + nsblock});
            dfh->fill_tensor("EYbs", Abs, {sstart, sstart + nsblock});

            double* DY2p = DYbsp[0];
            double* A2p = Absp[0];
            for (long int bsQ = 0L; bsQ < nsblock * nbQ; bsQ++) {
                (*DY2p++) += (*A2p++);
            }
            dfh->write_disk_tensor("FYbs", DYbs, {sstart, sstart + nsblock});
        }
    
        for (size_t rstart = 0; rstart < nr; rstart += max_r) {
            size_t nrblock = (rstart + max_r >= nr ? nr - rstart : max_r);
    
            dfh->fill_tensor("Aar", Aar, {rstart, rstart + nrblock});
            dfh->fill_tensor("Far", Dar, {rstart, rstart + nrblock});
            dfh->fill_tensor("Bbr", Bbr, {rstart, rstart + nrblock});
            dfh->fill_tensor("Cbr", Cbr, {rstart, rstart + nrblock});
            dfh->fill_tensor("FXar", DXar, {rstart, rstart + nrblock});
            dfh->fill_tensor("BXbr", BXbr, {rstart, rstart + nrblock});
            dfh->fill_tensor("CYbr", CYbr, {rstart, rstart + nrblock});
    
            for (size_t sstart = 0; sstart < ns; sstart += max_s) {
                size_t nsblock = (sstart + max_s >= ns ? ns - sstart : max_s);
    
                dfh->fill_tensor("Abs", Abs, {sstart, sstart + nsblock});
                dfh->fill_tensor("Fbs", Dbs, {sstart, sstart + nsblock});
                dfh->fill_tensor("Bas", Bas, {sstart, sstart + nsblock});
                dfh->fill_tensor("Cas", Cas, {sstart, sstart + nsblock});
                dfh->fill_tensor("FYbs", DYbs, {sstart, sstart + nsblock});
                dfh->fill_tensor("BYas", BYas, {sstart, sstart + nsblock});
                dfh->fill_tensor("CXas", CXas, {sstart, sstart + nsblock});

                long int nrs = nrblock * nsblock;

#pragma omp parallel for schedule(dynamic) reduction(+ : Disp20, ExchDisp20, sExchDisp20, par_ExchDisp20)
                for (long int rs = 0L; rs < nrs; rs++) {
                    int r = rs / nsblock;
                    int s = rs % nsblock;

                    int thread = 0;
#ifdef _OPENMP
                    thread = omp_get_thread_num();
#endif

                    double** E_disp20Tp = E_disp20_threads[thread]->pointer();
                    double** E_exch_disp20Tp = E_exch_disp20_threads[thread]->pointer();
                    double** sE_exch_disp20Tp = sE_exch_disp20_threads[thread]->pointer();
                    double** par_E_exch_disp20Tp = par_E_exch_disp20_threads[thread]->pointer();
 
                    double** Tabp = Tab[thread]->pointer();
                    double** Vabp = Vab[thread]->pointer();
                    double** T2abp = T2ab[thread]->pointer();
                    double** V2abp = V2ab[thread]->pointer();
                    double** Iabp = Iab[thread]->pointer();
 
                    // => Amplitudes, Disp20 <= //
 
                    C_DGEMM('N', 'T', na, nb, nQ, 1.0, Aarp[(r)*na], nQ, Absp[(s)*nb], nQ, 0.0, Vabp[0], nb);
                    for (int a = 0; a < na; a++) {
                        for (int b = 0; b < nb; b++) {
                            Tabp[a][b] = Vabp[a][b] / (eap[a] + ebp[b] - erp[r + rstart] - esp[s + sstart]);
                        }
                    }
 
                    C_DGEMM('N', 'N', na, nb, nb, 1.0, Tabp[0], nb, UBp[0], nb, 0.0, Iabp[0], nb);
                    C_DGEMM('T', 'N', na, nb, na, 1.0, UAp[0], na, Iabp[0], nb, 0.0, T2abp[0], nb);
                    C_DGEMM('N', 'N', na, nb, nb, 1.0, Vabp[0], nb, UBp[0], nb, 0.0, Iabp[0], nb);
                    C_DGEMM('T', 'N', na, nb, na, 1.0, UAp[0], na, Iabp[0], nb, 0.0, V2abp[0], nb);
 
                    for (int a = 0; a < na; a++) {
                        for (int b = 0; b < nb; b++) {
                            E_disp20Tp[a][b] += 4.0 * T2abp[a][b] * V2abp[a][b];
                            Disp20 += 4.0 * T2abp[a][b] * V2abp[a][b];
                        }
                    }

                    // => Exch-Disp20 <= //

                    // > Q1-Q3 < //

                    C_DGEMM('N', 'T', na, nb, nQ, 1.0, Basp[(s)*na], nQ, Bbrp[(r)*nb], nQ, 0.0, Vabp[0], nb);
                    C_DGEMM('N', 'T', na, nb, nQ, 1.0, Casp[(s)*na], nQ, Cbrp[(r)*nb], nQ, 1.0, Vabp[0], nb);
                    C_DGEMM('N', 'T', na, nb, nQ, 1.0, Aarp[(r)*na], nQ, Dbsp[(s)*nb], nQ, 1.0, Vabp[0], nb);
                    C_DGEMM('N', 'T', na, nb, nQ, 1.0, Darp[(r)*na], nQ, Absp[(s)*nb], nQ, 1.0, Vabp[0], nb);

                    // > V,J,K < //

                    C_DGER(na, nb, 1.0, &Sasp[0][s + sstart], ns, &Qbrp[0][r + rstart], nr, Vabp[0], nb);
                    C_DGER(na, nb, 1.0, &Qasp[0][s + sstart], ns, &Sbrp[0][r + rstart], nr, Vabp[0], nb);
                    C_DGER(na, nb, 1.0, &Qarp[0][r + rstart], nr, &SAbsp[0][s + sstart], ns, Vabp[0], nb);
                    C_DGER(na, nb, 1.0, &SBarp[0][r + rstart], nr, &Qbsp[0][s + sstart], ns, Vabp[0], nb);

                    C_DGEMM('N', 'N', na, nb, nb, 1.0, Vabp[0], nb, UBp[0], nb, 0.0, Iabp[0], nb);
                    C_DGEMM('T', 'N', na, nb, na, 1.0, UAp[0], na, Iabp[0], nb, 0.0, V2abp[0], nb);

                    for (int a = 0; a < na; a++) {
                        for (int b = 0; b < nb; b++) {
                            E_exch_disp20Tp[a][b] -= 2.0 * T2abp[a][b] * V2abp[a][b];
                            if (options_.get_bool("SSAPT0_SCALE"))
                                sE_exch_disp20Tp[a][b] -= scale * 2.0 * T2abp[a][b] * V2abp[a][b];
                            ExchDisp20 -= 2.0 * T2abp[a][b] * V2abp[a][b];
                            sExchDisp20 -= scale * 2.0 * T2abp[a][b] * V2abp[a][b];
                        }
                    }

                    // now, additional term for parallel/perpendicular spin coupling
                    // > Q1-Q3 < //

                    C_DGEMM('N', 'T', na, nb, nQ, 1.0, BYasp[(s)*na], nQ, BXbrp[(r)*nb], nQ, 0.0, Vabp[0], nb);
                    C_DGEMM('N', 'T', na, nb, nQ, 1.0, CXasp[(s)*na], nQ, CYbrp[(r)*nb], nQ, 1.0, Vabp[0], nb);
                    C_DGEMM('N', 'T', na, nb, nQ, 1.0, Aarp[(r)*na], nQ, DYbsp[(s)*nb], nQ, 1.0, Vabp[0], nb);
                    C_DGEMM('N', 'T', na, nb, nQ, 1.0, DXarp[(r)*na], nQ, Absp[(s)*nb], nQ, 1.0, Vabp[0], nb);

                    // > V,J,K < //

                    C_DGER(na, nb, 1.0, &Sasp[0][s + sstart], ns, &KXOYbrp[0][r + rstart], nr, Vabp[0], nb);
                    C_DGER(na, nb, 1.0, &KXOYasp[0][s + sstart], ns, &Sbrp[0][r + rstart], nr, Vabp[0], nb);

                    C_DGEMM('N', 'N', na, nb, nb, 1.0, Vabp[0], nb, UBp[0], nb, 0.0, Iabp[0], nb);
                    C_DGEMM('T', 'N', na, nb, na, 1.0, UAp[0], na, Iabp[0], nb, 0.0, V2abp[0], nb);

                    for (int a = 0; a < na; a++) {
                        for (int b = 0; b < nb; b++) {
                            par_E_exch_disp20Tp[a][b] -= 2.0 * T2abp[a][b] * V2abp[a][b];
                            par_ExchDisp20 -= 2.0 * T2abp[a][b] * V2abp[a][b];
                        }
                    }
                }
            }
        }
    }
    else {
        for (size_t rstart = 0; rstart < nr; rstart += max_r) {
            size_t nrblock = (rstart + max_r >= nr ? nr - rstart : max_r);

            dfh->fill_tensor("Aar", Aar, {rstart, rstart + nrblock});
            dfh->fill_tensor("Far", Dar, {rstart, rstart + nrblock});
            dfh->fill_tensor("Bbr", Bbr, {rstart, rstart + nrblock});
            dfh->fill_tensor("Cbr", Cbr, {rstart, rstart + nrblock});

            for (size_t sstart = 0; sstart < ns; sstart += max_s) {
                size_t nsblock = (sstart + max_s >= ns ? ns - sstart : max_s);

                dfh->fill_tensor("Abs", Abs, {sstart, sstart + nsblock});
                dfh->fill_tensor("Fbs", Dbs, {sstart, sstart + nsblock});
                dfh->fill_tensor("Bas", Bas, {sstart, sstart + nsblock});
                dfh->fill_tensor("Cas", Cas, {sstart, sstart + nsblock});

                long int nrs = nrblock * nsblock;

#pragma omp parallel for schedule(dynamic) reduction(+ : Disp20, ExchDisp20, sExchDisp20, par_ExchDisp20)
                for (long int rs = 0L; rs < nrs; rs++) {
                    int r = rs / nsblock;
                    int s = rs % nsblock;

                    int thread = 0;
#ifdef _OPENMP
                    thread = omp_get_thread_num();
#endif

                    double** E_disp20Tp = E_disp20_threads[thread]->pointer();
                    double** E_exch_disp20Tp = E_exch_disp20_threads[thread]->pointer();
                    double** sE_exch_disp20Tp = sE_exch_disp20_threads[thread]->pointer();
                    double** par_E_exch_disp20Tp = par_E_exch_disp20_threads[thread]->pointer();

                    double** Tabp = Tab[thread]->pointer();
                    double** Vabp = Vab[thread]->pointer();
                    double** T2abp = T2ab[thread]->pointer();
                    double** V2abp = V2ab[thread]->pointer();
                    double** Iabp = Iab[thread]->pointer();

                    // => Amplitudes, Disp20 <= //

                    C_DGEMM('N', 'T', na, nb, nQ, 1.0, Aarp[(r)*na], nQ, Absp[(s)*nb], nQ, 0.0, Vabp[0], nb);
                    for (int a = 0; a < na; a++) {
                        for (int b = 0; b < nb; b++) {
                            Tabp[a][b] = Vabp[a][b] / (eap[a] + ebp[b] - erp[r + rstart] - esp[s + sstart]);
                        }
                    }

                    C_DGEMM('N', 'N', na, nb, nb, 1.0, Tabp[0], nb, UBp[0], nb, 0.0, Iabp[0], nb);
                    C_DGEMM('T', 'N', na, nb, na, 1.0, UAp[0], na, Iabp[0], nb, 0.0, T2abp[0], nb);
                    C_DGEMM('N', 'N', na, nb, nb, 1.0, Vabp[0], nb, UBp[0], nb, 0.0, Iabp[0], nb);
                    C_DGEMM('T', 'N', na, nb, na, 1.0, UAp[0], na, Iabp[0], nb, 0.0, V2abp[0], nb);

                    for (int a = 0; a < na; a++) {
                        for (int b = 0; b < nb; b++) {
                            E_disp20Tp[a][b] += 4.0 * T2abp[a][b] * V2abp[a][b];
                            Disp20 += 4.0 * T2abp[a][b] * V2abp[a][b];
                        }
                    }

                    // => Exch-Disp20 <= //

                    // > Q1-Q3 < //

                    C_DGEMM('N', 'T', na, nb, nQ, 1.0, Basp[(s)*na], nQ, Bbrp[(r)*nb], nQ, 0.0, Vabp[0], nb);
                    C_DGEMM('N', 'T', na, nb, nQ, 1.0, Casp[(s)*na], nQ, Cbrp[(r)*nb], nQ, 1.0, Vabp[0], nb);
                    C_DGEMM('N', 'T', na, nb, nQ, 1.0, Aarp[(r)*na], nQ, Dbsp[(s)*nb], nQ, 1.0, Vabp[0], nb);
                    C_DGEMM('N', 'T', na, nb, nQ, 1.0, Darp[(r)*na], nQ, Absp[(s)*nb], nQ, 1.0, Vabp[0], nb);

                    // > V,J,K < //

                    C_DGER(na, nb, 1.0, &Sasp[0][s + sstart], ns, &Qbrp[0][r + rstart], nr, Vabp[0], nb);
                    C_DGER(na, nb, 1.0, &Qasp[0][s + sstart], ns, &Sbrp[0][r + rstart], nr, Vabp[0], nb);
                    C_DGER(na, nb, 1.0, &Qarp[0][r + rstart], nr, &SAbsp[0][s + sstart], ns, Vabp[0], nb);
                    C_DGER(na, nb, 1.0, &SBarp[0][r + rstart], nr, &Qbsp[0][s + sstart], ns, Vabp[0], nb);

                    C_DGEMM('N', 'N', na, nb, nb, 1.0, Vabp[0], nb, UBp[0], nb, 0.0, Iabp[0], nb);
                    C_DGEMM('T', 'N', na, nb, na, 1.0, UAp[0], na, Iabp[0], nb, 0.0, V2abp[0], nb);

                    for (int a = 0; a < na; a++) {
                        for (int b = 0; b < nb; b++) {
                            E_exch_disp20Tp[a][b] -= 2.0 * T2abp[a][b] * V2abp[a][b];
                            if (options_.get_bool("SSAPT0_SCALE"))
                                sE_exch_disp20Tp[a][b] -= scale * 2.0 * T2abp[a][b] * V2abp[a][b];
                            ExchDisp20 -= 2.0 * T2abp[a][b] * V2abp[a][b];
                            sExchDisp20 -= scale * 2.0 * T2abp[a][b] * V2abp[a][b];
                        }
                    }
                }
            }
        }
    }

    auto E_disp20 = std::make_shared<Matrix>("E_disp20", na, nb);
    auto E_exch_disp20 = std::make_shared<Matrix>("E_exch_disp20", na, nb);
    double** E_disp20p = E_disp20->pointer();
    double** E_exch_disp20p = E_exch_disp20->pointer();

    for (int t = 0; t < nT; t++) {
        E_disp20->add(E_disp20_threads[t]);
        E_exch_disp20->add(E_exch_disp20_threads[t]);
    }

    for (int a = 0; a < na; a++) {
        for (int b = 0; b < nb; b++) {
            Ep[a + nfa + nA][b + nfb + nB] = E_disp20p[a][b] + E_exch_disp20p[a][b];
        }
    }

    if (options_.get_bool("SSAPT0_SCALE")) {
        auto sE_exch_disp20 = std::make_shared<Matrix>("sE_exch_disp20", na, nb);
        sE_exch_disp20->copy(E_exch_disp20);
        double** sE_exch_disp20p = sE_exch_disp20->pointer();
        sE_exch_disp20->scale(sSAPT0_scale_);

        for (int a = 0; a < na; a++) {
            for (int b = 0; b < nb; b++) {
                sEp[a + nfa + nA][b + nfb + nB] = E_disp20p[a][b] + sE_exch_disp20p[a][b];
            }
        }
    }

    // E_disp20->print();
    // E_exch_disp20->print();

    scalars_["Disp20"] = Disp20;
    scalars_["Exch-Disp20"] = ExchDisp20;
    if (options_.get_bool("SSAPT0_SCALE")) scalars_["sExch-Disp20"] = sExchDisp20;
    outfile->Printf("    Disp20              = %18.12lf [Eh]\n", Disp20);
    outfile->Printf("    Exch-Disp20         = %18.12lf [Eh]\n", ExchDisp20);
    if ((link_assignment == "SAO0" || link_assignment == "SAO1" || link_assignment == "SAO2" || link_assignment == "SIAO0" || link_assignment == "SIAO1" || link_assignment == "SIAO2") && parperp) {
        outfile->Printf("    Exch-Disp20 (PAR)   = %18.12lf [Eh]\n", ExchDisp20 + par_ExchDisp20);
        outfile->Printf("    Exch-Disp20 (PERP)  = %18.12lf [Eh]\n", ExchDisp20 - par_ExchDisp20);
    }
    if (options_.get_bool("SSAPT0_SCALE")) outfile->Printf("    sExch-Disp20         = %18.12lf [Eh]\n", sExchDisp20);
    outfile->Printf("\n");
    // fflush(outfile);
}

std::shared_ptr<Matrix> FISAPT::extract_columns(const std::vector<int>& cols, std::shared_ptr<Matrix> A) {
    int nm = A->rowspi()[0];
    int na = A->colspi()[0];
    int ni = cols.size();

    auto A2 = std::make_shared<Matrix>("A2", nm, ni);
    double** Ap = A->pointer();
    double** A2p = A2->pointer();

    for (int m = 0; m < nm; m++) {
        for (int i = 0; i < ni; i++) {
            A2p[m][i] = Ap[m][cols[i]];
        }
    }

    return A2;
}

FISAPTSCF::FISAPTSCF(std::shared_ptr<JK> jk, double enuc, std::shared_ptr<Matrix> S, std::shared_ptr<Matrix> X,
                     std::shared_ptr<Matrix> T, std::shared_ptr<Matrix> V, std::shared_ptr<Matrix> W,
                     std::shared_ptr<Matrix> C, Options& options)
    : options_(options), jk_(jk) {
    scalars_["E NUC"] = enuc;
    matrices_["S"] = S;
    matrices_["X"] = X;
    matrices_["T"] = T;
    matrices_["V"] = V;
    matrices_["W"] = W;
    matrices_["C0"] = C;
}

FISAPTSCF::~FISAPTSCF() {}
void FISAPTSCF::compute_energy() {
    // => Sizing <= //
    int nbf = matrices_["X"]->rowspi()[0];
    int nmo = matrices_["X"]->colspi()[0];
    int nocc = matrices_["C0"]->colspi()[0];
    int nvir = nmo - nocc;

    // => One-electron potential <= //

    matrices_["H"] = std::shared_ptr<Matrix>(matrices_["T"]->clone());
    matrices_["H"]->set_name("H");
    matrices_["H"]->copy(matrices_["T"]);
    matrices_["H"]->add(matrices_["V"]);
    // matrices_["H"]->add(matrices_["W"]);

    // => Fock Matrix <= //

    matrices_["F"] = std::shared_ptr<Matrix>(matrices_["T"]->clone());
    matrices_["F"]->set_name("F");

    // => For Convenience <= //

    auto H = matrices_["H"];
    auto F = matrices_["F"];
    auto S = matrices_["S"];
    auto X = matrices_["X"];
    auto W = matrices_["W"];

    // => Guess <= //

    std::shared_ptr<Matrix> Cocc2(matrices_["C0"]->clone());
    Cocc2->copy(matrices_["C0"]);

    // => Convergence Criteria <= //

    int maxiter = options_.get_int("MAXITER");
    double Etol = options_.get_double("E_CONVERGENCE");
    double Gtol = options_.get_double("D_CONVERGENCE");
    bool converged = false;
    double Eold = 0.0;

    outfile->Printf("    Maxiter = %11d\n", maxiter);
    outfile->Printf("    E Tol   = %11.3E\n", Etol);
    outfile->Printf("    D Tol   = %11.3E\n", Gtol);
    outfile->Printf("\n");

    // => DIIS Setup <= //

    int max_diis_vectors = options_.get_int("DIIS_MAX_VECS");
    outfile->Printf("    Max DIIS Vectors = %d\n", max_diis_vectors);
    outfile->Printf("\n");

    bool diised = false;
    auto Gsize = std::make_shared<Matrix>("Gsize", nmo, nmo);
    DIISManager diis(max_diis_vectors, "FISAPT DIIS");
    diis.set_error_vector_size(Gsize.get());
    diis.set_vector_size(F.get());
    Gsize.reset();

    // ==> Master Loop <== //

    outfile->Printf("    Iter %3s: %24s %11s %11s\n", "N", "E", "dE", "|D|");
    for (int iter = 1; iter <= maxiter; iter++) {
        // => Compute Density Matrix <= //

        auto D = linalg::doublet(Cocc2, Cocc2, false, true);

        // => Compute Fock Matrix <= //

        auto& Cl = jk_->C_left();
        auto& Cr = jk_->C_right();

        const auto& Js = jk_->J();
        const auto& Ks = jk_->K();

        Cl.clear();
        Cr.clear();

        Cl.push_back(Cocc2);
        Cr.push_back(Cocc2);

        jk_->compute();

        auto J = Js[0];
        auto K = Ks[0];

        F->copy(H);
        F->add(W);
        F->add(J);
        F->add(J);
        F->subtract(K);

        // => Compute Energy <= //

        double E = scalars_["E NUC"] + D->vector_dot(H) + D->vector_dot(F) + D->vector_dot(W);
        double Ediff = E - Eold;
        scalars_["E SCF"] = E;

        // => Compute Orbital Gradient <= //

        auto G1 = linalg::triplet(F, D, S);
        auto G2 = G1->transpose();
        G1->subtract(G2);
        G1->transform(X);
        double Gnorm = G1->rms();

        // => Print and Check Convergence <= //

        outfile->Printf("    Iter %3d: %24.16E %11.3E %11.3E %s\n", iter, E, Ediff, Gnorm, (diised ? "DIIS" : ""));

        if (std::fabs(Ediff) < Etol && std::fabs(Gnorm) < Gtol) {
            converged = true;
            break;
        }

        Eold = E;

        // => DIIS <= //

        diis.add_entry(G1.get(), F.get());
        diised = diis.extrapolate(F.get());

        // => Diagonalize Fock Matrix <= //

        auto F2 = linalg::triplet(X, F, X, true, false, false);
        auto U2 = std::make_shared<Matrix>("C", nmo, nmo);
        auto e2 = std::make_shared<Vector>("eps", nmo);
        F2->diagonalize(U2, e2, ascending);
        auto C = linalg::doublet(X, U2, false, false);

        // => Assign New Orbitals <= //

        double** Coccp = Cocc2->pointer();
        double** Cp = C->pointer();
        for (int m = 0; m < nbf; m++) {
            for (int i = 0; i < nocc; i++) {
                Coccp[m][i] = Cp[m][i];
            }
        }

        matrices_["C"] = C;
        vectors_["eps"] = e2;
    }
    outfile->Printf("\n");

    if (converged) {
        outfile->Printf("    FISAPTSCF Converged.\n\n");
    } else {
        outfile->Printf("    FISAPTSCF Failed.\n\n");
    }

    // => Post Results <= //

    auto eps = vectors_["eps"];
    auto eps_occ = std::make_shared<Vector>("eps_occ", nocc);
    auto eps_vir = std::make_shared<Vector>("eps_vir", nvir);

    auto ep = eps->pointer();
    auto eop = eps_occ->pointer();
    auto evp = eps_vir->pointer();

    for (int i = 0; i < nocc; i++) {
        eop[i] = ep[i];
    }

    for (int a = 0; a < nvir; a++) {
        evp[a] = ep[a + nocc];
    }

    vectors_["eps_occ"] = eps_occ;
    vectors_["eps_vir"] = eps_vir;

    auto C = matrices_["C"];
    auto Cocc = std::make_shared<Matrix>("Cocc", nbf, nocc);
    auto Cvir = std::make_shared<Matrix>("Cvir", nbf, nvir);

    auto Cp = C->pointer();
    auto Cop = Cocc->pointer();
    auto Cvp = Cvir->pointer();

    for (int m = 0; m < nbf; m++) {
        for (int i = 0; i < nocc; i++) {
            Cop[m][i] = Cp[m][i];
        }
    }

    for (int m = 0; m < nbf; m++) {
        for (int a = 0; a < nvir; a++) {
            Cvp[m][a] = Cp[m][a + nocc];
        }
    }

    matrices_["Cocc"] = Cocc;
    matrices_["Cvir"] = Cvir;

    const auto& Js = jk_->J();
    const auto& Ks = jk_->K();

    matrices_["J"] = std::shared_ptr<Matrix>(Js[0]->clone());
    matrices_["K"] = std::shared_ptr<Matrix>(Ks[0]->clone());
    matrices_["J"]->copy(Js[0]);
    matrices_["K"]->copy(Ks[0]);
    matrices_["J"]->set_name("J");
    matrices_["K"]->set_name("K");

    // => Print Final Info <= //

    outfile->Printf("    Final SCF Energy: %24.16E [Eh]\n\n", scalars_["E SCF"]);

    print_orbitals("Occupied Orbital Energies", 1, eps_occ);
    print_orbitals("Virtual Orbital Energies", nocc + 1, eps_vir);
}
void FISAPTSCF::print_orbitals(const std::string& header, int start, std::shared_ptr<Vector> eps) {
    outfile->Printf("   => %s <=\n\n", header.c_str());
    outfile->Printf("    ");
    int n = eps->dimpi()[0];
    double* ep = eps->pointer();
    int count = 0;
    for (int i = 0; i < n; i++) {
        outfile->Printf("%4d %11.6f  ", i + start, ep[i]);
        if (count++ % 3 == 2 && count != n) outfile->Printf("\n    ");
    }
    outfile->Printf("\n\n");
}

CPHF_FISAPT::CPHF_FISAPT() {}
CPHF_FISAPT::~CPHF_FISAPT() {}
void CPHF_FISAPT::compute_cphf() {
    // Allocate
    x_A_ = std::shared_ptr<Matrix>(w_A_->clone());
    x_B_ = std::shared_ptr<Matrix>(w_B_->clone());
    x_A_->zero();
    x_B_->zero();

    std::shared_ptr<Matrix> r_A(w_A_->clone());
    std::shared_ptr<Matrix> z_A(w_A_->clone());
    std::shared_ptr<Matrix> p_A(w_A_->clone());
    std::shared_ptr<Matrix> r_B(w_B_->clone());
    std::shared_ptr<Matrix> z_B(w_B_->clone());
    std::shared_ptr<Matrix> p_B(w_B_->clone());

    // Initialize (x_0 = 0)
    r_A->copy(w_A_);
    r_B->copy(w_B_);

    preconditioner(r_A, z_A, eps_occ_A_, eps_vir_A_);
    preconditioner(r_B, z_B, eps_occ_B_, eps_vir_B_);

    // Uncoupled value
    // outfile->Printf("(A<-B): %24.16E\n", -2.0 * z_A->vector_dot(w_A_));
    // outfile->Printf("(B<-A): %24.16E\n", -2.0 * z_B->vector_dot(w_B_));

    p_A->copy(z_A);
    p_B->copy(z_B);

    double zr_old_A = z_A->vector_dot(r_A);
    double zr_old_B = z_B->vector_dot(r_B);

    double r2A = 1.0;
    double r2B = 1.0;

    double b2A = sqrt(w_A_->vector_dot(w_A_));
    double b2B = sqrt(w_B_->vector_dot(w_B_));

    outfile->Printf("  ==> CPHF Iterations <==\n\n");

    outfile->Printf("    Maxiter     = %11d\n", maxiter_);
    outfile->Printf("    Convergence = %11.3E\n", delta_);
    outfile->Printf("\n");

    std::time_t start;
    std::time_t stop;

    start = std::time(nullptr);

    outfile->Printf("    -----------------------------------------\n");
    outfile->Printf("    %-4s %11s  %11s  %10s\n", "Iter", "Monomer A", "Monomer B", "Time [s]");
    outfile->Printf("    -----------------------------------------\n");
    // fflush(outfile);

    int iter;
    for (iter = 0; iter < maxiter_; iter++) {
        std::map<std::string, std::shared_ptr<Matrix> > b;
        if (r2A > delta_) {
            b["A"] = p_A;
        }
        if (r2B > delta_) {
            b["B"] = p_B;
        }

        std::map<std::string, std::shared_ptr<Matrix> > s = product(b);

        if (r2A > delta_) {
            std::shared_ptr<Matrix> s_A = s["A"];
            double alpha = r_A->vector_dot(z_A) / p_A->vector_dot(s_A);
            if (alpha < 0.0) {
                throw PSIEXCEPTION("Monomer A: A Matrix is not SPD");
            }
            size_t no = x_A_->nrow();
            size_t nv = x_A_->ncol();
            double** xp = x_A_->pointer();
            double** rp = r_A->pointer();
            double** pp = p_A->pointer();
            double** sp = s_A->pointer();
            C_DAXPY(no * nv, alpha, pp[0], 1, xp[0], 1);
            C_DAXPY(no * nv, -alpha, sp[0], 1, rp[0], 1);
            r2A = sqrt(C_DDOT(no * nv, rp[0], 1, rp[0], 1)) / b2A;
        }

        if (r2B > delta_) {
            std::shared_ptr<Matrix> s_B = s["B"];
            double alpha = r_B->vector_dot(z_B) / p_B->vector_dot(s_B);
            if (alpha < 0.0) {
                throw PSIEXCEPTION("Monomer B: A Matrix is not SPD");
            }
            int no = x_B_->nrow();
            int nv = x_B_->ncol();
            double** xp = x_B_->pointer();
            double** rp = r_B->pointer();
            double** pp = p_B->pointer();
            double** sp = s_B->pointer();
            C_DAXPY(no * nv, alpha, pp[0], 1, xp[0], 1);
            C_DAXPY(no * nv, -alpha, sp[0], 1, rp[0], 1);
            r2B = sqrt(C_DDOT(no * nv, rp[0], 1, rp[0], 1)) / b2B;
        }

        stop = std::time(nullptr);
        outfile->Printf("    %-4d %11.3E%1s %11.3E%1s %10ld\n", iter + 1, r2A, (r2A < delta_ ? "*" : " "), r2B,
                        (r2B < delta_ ? "*" : " "), stop - start);
        // fflush(outfile);

        if (r2A <= delta_ && r2B <= delta_) {
            break;
        }

        if (r2A > delta_) {
            preconditioner(r_A, z_A, eps_occ_A_, eps_vir_A_);
            double zr_new = z_A->vector_dot(r_A);
            double beta = zr_new / zr_old_A;
            zr_old_A = zr_new;
            int no = x_A_->nrow();
            int nv = x_A_->ncol();
            double** pp = p_A->pointer();
            double** zp = z_A->pointer();
            C_DSCAL(no * nv, beta, pp[0], 1);
            C_DAXPY(no * nv, 1.0, zp[0], 1, pp[0], 1);
        }

        if (r2B > delta_) {
            preconditioner(r_B, z_B, eps_occ_B_, eps_vir_B_);
            double zr_new = z_B->vector_dot(r_B);
            double beta = zr_new / zr_old_B;
            zr_old_B = zr_new;
            int no = x_B_->nrow();
            int nv = x_B_->ncol();
            double** pp = p_B->pointer();
            double** zp = z_B->pointer();
            C_DSCAL(no * nv, beta, pp[0], 1);
            C_DAXPY(no * nv, 1.0, zp[0], 1, pp[0], 1);
        }
    }

    outfile->Printf("    -----------------------------------------\n");
    outfile->Printf("\n");
    // fflush(outfile);

    if (iter == maxiter_) throw PSIEXCEPTION("CPHF did not converge.");
}

void CPHF_FISAPT::preconditioner(std::shared_ptr<Matrix> r, std::shared_ptr<Matrix> z, std::shared_ptr<Vector> o,
                                 std::shared_ptr<Vector> v) {
    int no = o->dim();
    int nv = v->dim();

    double** rp = r->pointer();
    double** zp = z->pointer();

    double* op = o->pointer();
    double* vp = v->pointer();

    for (int i = 0; i < no; i++) {
        for (int a = 0; a < nv; a++) {
            zp[i][a] = rp[i][a] / (vp[a] - op[i]);
        }
    }
}

std::map<std::string, std::shared_ptr<Matrix> > CPHF_FISAPT::product(
    std::map<std::string, std::shared_ptr<Matrix> > b) {
    std::map<std::string, std::shared_ptr<Matrix> > s;

    bool do_A = b.count("A");
    bool do_B = b.count("B");

    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    if (do_A) {
        Cl.push_back(Cocc_A_);
        auto T = linalg::doublet(Cvir_A_, b["A"], false, true);
        Cr.push_back(T);
    }

    if (do_B) {
        Cl.push_back(Cocc_B_);
        auto T = linalg::doublet(Cvir_B_, b["B"], false, true);
        Cr.push_back(T);
    }

    jk_->compute();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    int indA = 0;
    int indB = (do_A ? 1 : 0);

    if (do_A) {
        std::shared_ptr<Matrix> Jv = J[indA];
        std::shared_ptr<Matrix> Kv = K[indA];
        Jv->scale(4.0);
        Jv->subtract(Kv);
        Jv->subtract(Kv->transpose());

        int no = b["A"]->nrow();
        int nv = b["A"]->ncol();
        int nso = Cvir_A_->nrow();
        s["A"] = linalg::triplet(Cocc_A_, Jv, Cvir_A_, true, false, false);
        auto Sp = s["A"]->pointer();

        auto bp = b["A"]->pointer();
        auto op = eps_occ_A_->pointer();
        auto vp = eps_vir_A_->pointer();
        for (int i = 0; i < no; i++) {
            for (int a = 0; a < nv; a++) {
                Sp[i][a] += bp[i][a] * (vp[a] - op[i]);
            }
        }
    }

    if (do_B) {
        std::shared_ptr<Matrix> Jv = J[indB];
        std::shared_ptr<Matrix> Kv = K[indB];
        Jv->scale(4.0);
        Jv->subtract(Kv);
        Jv->subtract(Kv->transpose());

        int no = b["B"]->nrow();
        int nv = b["B"]->ncol();
        s["B"] = linalg::triplet(Cocc_B_, Jv, Cvir_B_, true, false, false);
        auto Sp = s["B"]->pointer();

        auto bp = b["B"]->pointer();
        auto op = eps_occ_B_->pointer();
        auto vp = eps_vir_B_->pointer();
        for (int i = 0; i < no; i++) {
            for (int a = 0; a < nv; a++) {
                Sp[i][a] += bp[i][a] * (vp[a] - op[i]);
            }
        }
    }

    return s;
}

}  // Namespace fisapt

double sapt_nuclear_external_potential_matrix(
    std::shared_ptr<Wavefunction> reference_,
    std::map<std::string, std::shared_ptr<Matrix>>& matrices_,
    Options& options_
) {
    // => External potential <= //
    // for every key in matrices_ print matrix name
    double** Enucsp = matrices_["Enucs"]->pointer();

    std::shared_ptr<BasisSet> primary_ = reference_->basisset();
    std::shared_ptr<Molecule> mol = reference_->molecule();
    std::vector<std::shared_ptr<ExternalPotential>> pot_list;
    std::vector<int> pot_ids;
    double Etot = 0.0;


    // this is where any "external potential C" goes... any common potential felt by all subsystems

    if (reference_->has_potential_variable("C")) {
        pot_list.push_back(reference_->potential_variable("C"));
        pot_ids.push_back(2); // code 2 == C

        if (options_.get_bool("EXTERNAL_POTENTIAL_SYMMETRY") == false && reference_->nirrep() != 1)
            throw PSIEXCEPTION("SCF: External Fields are not consistent with symmetry. Set symmetry c1.");

        std::shared_ptr<Matrix> V_extern = reference_->potential_variable("C")->computePotentialMatrix(primary_);

        if (options_.get_bool("EXTERNAL_POTENTIAL_SYMMETRY")) {
            // Attempt to apply symmetry. No error checking is performed.
            std::shared_ptr<Matrix> V_extern_sym = reference_->matrix_factory()->create_shared_matrix("External Potential");
            V_extern_sym->apply_symmetry(V_extern, reference_->aotoso());
            V_extern = V_extern_sym;
        }

        if (reference_->get_print()) {
            reference_->potential_variable("C")->set_print(reference_->get_print());
            outfile->Printf("  External Potential C:\n");
            reference_->potential_variable("C")->print();
        }

        if (reference_->get_print() > 3) V_extern->print();


        // Save external potential to add to one-electron SCF potential
        matrices_["VE"] = V_extern;
    }

    std::vector<std::string> subsystem_labels = {"A", "B"};
    std::vector<int> potential_ids = {0, 1};

    for (int i=0; i<2; i++) {
        if (reference_->has_potential_variable(subsystem_labels[i])) {
            outfile->Printf("  External Potential %s:\n", subsystem_labels[i].c_str());
            
            pot_list.push_back(reference_->potential_variable(subsystem_labels[i]));
            pot_ids.push_back(potential_ids[i]);

            if (options_.get_bool("EXTERNAL_POTENTIAL_SYMMETRY") == false && reference_->nirrep() != 1)
                throw PSIEXCEPTION("SCF: External Fields are not consistent with symmetry. Set symmetry c1.");

            std::shared_ptr<Matrix> V_extern = reference_->potential_variable(subsystem_labels[i])->computePotentialMatrix(primary_);
            // if (matrices_["V" + subsystem_labels[i]]) {
            //     V_extern = matrices_["V" + subsystem_labels[i]];
            // }
            // V_extern->ncol()
            printf("V_extern: ncol %d, nrow: %d", V_extern->ncol(), V_extern->nrow());

            if (options_.get_bool("EXTERNAL_POTENTIAL_SYMMETRY")) {
                // Attempt to apply symmetry. No error checking is performed.
                std::shared_ptr<Matrix> V_extern_sym = reference_->matrix_factory()->create_shared_matrix("External Potential");
                V_extern_sym->apply_symmetry(V_extern, reference_->aotoso());
                V_extern = V_extern_sym;
            }

            if (reference_->get_print()) {
                reference_->potential_variable(subsystem_labels[i])->set_print(reference_->get_print());
                outfile->Printf("  External Potential " + subsystem_labels[i] + ":\n");
                reference_->potential_variable(subsystem_labels[i])->print();
            }
            if (reference_->get_print() > 3) V_extern->print();

            matrices_["V" + subsystem_labels[i]]->add(V_extern);
            matrices_["V" + subsystem_labels[i] + "_extern"] = V_extern;
        } else{
            outfile->Printf("  No External Potential %s\n", subsystem_labels[i].c_str());
        }
    }
    // matrices_["VA"]->print();

    // are there any external potentials? If so, we need all QM atoms to feel their effects in the nuclear 
    // repulsion energy.  
    // Apparently, we were using C full strength, but the others I think get scaled by 0.5 because Rob counts
    // A->B and B->A separately and adds them (see a few lines up this fn...maybe due to how FSAPT files are written) 

    if (pot_list.size() > 0) { 

        for (int A = 0; A < mol->nfragments(); A++) {
            outfile->Printf("           Old Nuclear Repulsion %c: %24.16E [Eh]\n", 'A'+A, Enucsp[A][A]);
        }

        // First compute mol-extern interaction
        // a maximum of 3 external potentials and 3 fragments
        auto extern_mol_IE_mat = std::make_shared<Matrix>("extern_mol_IE", 3, 3);
        auto extern_mol_IE_matp = extern_mol_IE_mat->pointer();

        std::vector<int> none;
        std::vector<int> frag_list(1);
        double Enuc_extern;
        char frag = '@'; // Next characters are 'A', 'B', 'C'
        for (int p = 0; p < pot_list.size(); p++) {
            for (int A = 0; A < mol->nfragments(); A++) {
                frag++;
                frag_list[0] = A;
                auto mol_frag = mol->extract_subsets(frag_list, none);
                Enuc_extern = pot_list[p]->computeNuclearEnergy(mol_frag);
                Etot += Enuc_extern;

                Enucsp[A][pot_ids[p]] += Enuc_extern * 0.5;
                Enucsp[pot_ids[p]][A] += Enuc_extern * 0.5;

                extern_mol_IE_matp[pot_ids[p]][A] = Enuc_extern;

            } // end frag loop
        } // end pot loop

        for (int A = 0; A < mol->nfragments(); A++) {
            outfile->Printf("       Updated Nuclear Repulsion %c: %24.16E [Eh]\n", 'A'+A, Enucsp[A][A]);
        }

        outfile->Printf("\n");

        matrices_["extern_mol_IE"] = extern_mol_IE_mat; // save interaction energy between external potential and molecule

        // Now compute extern-extern interaction
        // We store the interaction energy between the external potentials in A, B, and C.
        auto extern_extern_IE_mat = std::make_shared<Matrix>("extern_extern_IE", 3, 3);
        double** extern_extern_IE_matp = extern_extern_IE_mat->pointer();
        double extern_extern_IE = 0.0;
        for (int p1 = 0; p1 < pot_list.size(); p1++) {
            for (int p2 = p1+1; p2 < pot_list.size(); p2++) {

                double IE = pot_list[p1]->computeExternExternInteraction(pot_list[p2]);
                // store half the interaction so that Eij + Eji = Etotal
                extern_extern_IE_matp[pot_ids[p1]][pot_ids[p2]] = IE * 0.5;
                extern_extern_IE_matp[pot_ids[p2]][pot_ids[p1]] = IE * 0.5;
                extern_extern_IE += IE;

                outfile->Printf("    Interaction Energy between External Potentials %c and %c: %24.16E [Eh]\n",
                                'A'+pot_ids[p1], 'A'+pot_ids[p2], IE);
            }
        }

        outfile->Printf("\n");
        outfile->Printf("    Total Interaction Energy between External Potentials: %24.16E [Eh]\n", extern_extern_IE);
        matrices_["extern_extern_IE"] = extern_extern_IE_mat;
    }
    // matrices_["VE"]->print();
    // for (auto const& x : matrices_) {
    //     std::string key = x.first;
    //     std::shared_ptr<Matrix> matrix = x.second;
    //     outfile->Printf("  Matrix %s\n", key.c_str());
    // }
    outfile->Printf("       Total Nuclear Repulsion: %24.16E [Eh]\n", Etot);
    return Etot;
}
}  // Namespace psi
