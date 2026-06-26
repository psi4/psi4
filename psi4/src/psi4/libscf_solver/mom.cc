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

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <utility>

#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"

#include "hf.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {
namespace scf {

void HF::MOM_start() {
    // Perhaps no MOM (or at least no MOM_start())?
    if (iteration_ != options_.get_int("MOM_START") || options_.get_int("MOM_START") == 0 || MOM_performed_) return;

    // If we're here, we're doing MOM of some kind
    MOM_performed_ = true;  // Gets printed next iteration

    // Build Ca_old_ matrices
    Ca_old_ = std::make_shared<Matrix>("C Alpha Old (SO Basis)", nsopi_, nmopi_);
    if (!same_a_b_orbs()) {
        Cb_old_ = std::make_shared<Matrix>("C Beta Old (SO Basis)", nsopi_, nmopi_);
    } else {
        Cb_old_ = Ca_old_;
    }

    Ca_old_->copy(Ca_);
    Cb_old_->copy(Cb_);

    // If no excitation requested, it's a stabilizer MOM, nothing fancy, don't print
    if (!options_["MOM_OCC"].size()) return;

    // If we're here, its an exciting MOM
    outfile->Printf("\n");
    print_orbitals();
    outfile->Printf("\n  ==> MOM Excited-State Iterations <==\n\n");

    // Reset DIIS (will automagically restart)
    if (initialized_diis_manager_) {
        diis_manager_.attr("delete_diis_file")();
        diis_manager_ = py::none();
        initialized_diis_manager_ = false;
    }

    // Sort MOs by (energy, (irrep, # within irrep))
    std::vector<std::pair<double, std::pair<int, int> > > orbs_a;
    for (int h = 0; h < nirrep_; h++) {
        int nmo = nmopi_[h];
        if (nmo == 0) continue;
        double* eps = epsilon_a_->pointer(h);
        for (int a = 0; a < nmo; a++) orbs_a.push_back(std::make_pair(eps[a], std::make_pair(h, a)));
    }
    std::vector<std::pair<double, std::pair<int, int> > > orbs_b;
    for (int h = 0; h < nirrep_; h++) {
        int nmo = nmopi_[h];
        if (nmo == 0) continue;
        double* eps = epsilon_b_->pointer(h);
        for (int a = 0; a < nmo; a++) orbs_b.push_back(std::make_pair(eps[a], std::make_pair(h, a)));
    }
    sort(orbs_a.begin(), orbs_a.end());
    sort(orbs_b.begin(), orbs_b.end());

    // Validate the orbitals to excite to/from
    if (options_["MOM_OCC"].size() != options_["MOM_VIR"].size())
        throw PSIEXCEPTION("SCF: MOM_OCC and MOM_VIR are not the same size");

    std::set<int> tot_check;

    for (int n = 0; n < options_["MOM_OCC"].size(); n++) {
        tot_check.insert(options_["MOM_OCC"][n].to_integer());
        tot_check.insert(options_["MOM_VIR"][n].to_integer());
    }

    if (tot_check.size() != 2 * options_["MOM_OCC"].size())
        throw PSIEXCEPTION("SCF::MOM_start: Occupied/Virtual index array elements are not unique");

    CharacterTable ct = molecule_->point_group()->char_table();

    outfile->Printf("  Excitations:\n");

    for (int n = 0; n < options_["MOM_OCC"].size(); n++) {
        int i = options_["MOM_OCC"][n].to_integer();
        int a = options_["MOM_VIR"][n].to_integer();

        bool si = (i > 0);
        bool sa = (a > 0);

        i = std::abs(i) - 1;
        a = std::abs(a) - 1;

        if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {
            if (!si || !sa)
                throw PSIEXCEPTION(
                    "SCF::MOM_start: RHF/RKS requires + -> + in input, as only double excitations are performed");

            int hi = orbs_a[i].second.first;
            int ha = orbs_a[a].second.first;

            if (hi == ha) {
                // Same irrep
                int pi = orbs_a[i].second.second;
                int pa = orbs_a[a].second.second;

                int nso = nsopi_[hi];
                int nmo = nmopi_[hi];

                double** Ca = Ca_->pointer(hi);
                double* Ct = new double[nso];
                double* eps = epsilon_a_->pointer(hi);
                double epst;

                // Swap eigvals
                epst = eps[pi];
                eps[pi] = eps[pa];
                eps[pa] = epst;

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pa], nmo, &Ca[0][pi], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pa], nmo);

                delete[] Ct;

                outfile->Printf("   %8s: %4d%-4s -> %4d%-4s \n", "AB -> AB", pi + 1, ct.gamma(hi).symbol(), pa + 1,
                                ct.gamma(ha).symbol());
            } else {
                // Different irrep
                // Occ -> Vir
                int pi = orbs_a[i].second.second;
                int pi2 = nalphapi_[hi] - 1;

                int nso = nsopi_[hi];
                int nmo = nmopi_[hi];

                double** Ca = Ca_->pointer(hi);
                double* Ct = new double[nso];
                double* eps = epsilon_a_->pointer(hi);
                double epst;

                // Swap eigvals
                epst = eps[pi];
                eps[pi] = eps[pi2];
                eps[pi2] = epst;

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);

                delete[] Ct;

                // Redo indexing
                nalphapi_[hi]--;
                nbetapi_[hi]--;

                // Vir -> Occ
                int pa = orbs_a[a].second.second;
                int pa2 = nalphapi_[ha];

                nso = nsopi_[ha];
                nmo = nmopi_[ha];

                Ca = Ca_->pointer(ha);
                Ct = new double[nso];
                eps = epsilon_a_->pointer(ha);

                // Swap eigvals
                epst = eps[pa];
                eps[pa] = eps[pa2];
                eps[pa2] = epst;

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);

                delete[] Ct;

                // Redo indexing
                nalphapi_[ha]++;
                nbetapi_[ha]++;

                outfile->Printf("   %8s: %4d%-4s -> %4d%-4s \n", "AB -> AB", pi + 1, ct.gamma(hi).symbol(), pa + 1,
                                ct.gamma(ha).symbol());
            }

        } else if (options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS") {
            if (si && sa) {
                // Alpha - alpha
                int hi = orbs_a[i].second.first;
                int ha = orbs_a[a].second.first;

                if (hi == ha) {
                    // Same irrep
                    int pi = orbs_a[i].second.second;
                    int pa = orbs_a[a].second.second;

                    int nso = nsopi_[hi];
                    int nmo = nmopi_[hi];

                    double** Ca = Ca_->pointer(hi);
                    double* Ct = new double[nso];
                    double* eps = epsilon_a_->pointer(hi);
                    double epst;

                    // Swap eigvals
                    epst = eps[pi];
                    eps[pi] = eps[pa];
                    eps[pa] = epst;

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pa], nmo, &Ca[0][pi], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pa], nmo);

                    delete[] Ct;

                    outfile->Printf("   %8s: %4d%-4s -> %4d%-4s \n", "A  -> A ", pi + 1, ct.gamma(hi).symbol(), pa + 1,
                                    ct.gamma(ha).symbol());
                } else {
                    // Different irrep
                    // Occ -> Vir
                    int pi = orbs_a[i].second.second;
                    int pi2 = nalphapi_[hi] - 1;

                    int nso = nsopi_[hi];
                    int nmo = nmopi_[hi];

                    double** Ca = Ca_->pointer(hi);
                    double* Ct = new double[nso];
                    double* eps = epsilon_a_->pointer(hi);
                    double epst;

                    // Swap eigvals
                    epst = eps[pi];
                    eps[pi] = eps[pi2];
                    eps[pi2] = epst;

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);

                    delete[] Ct;

                    // Redo indexing
                    nalphapi_[hi]--;

                    // Vir -> Occ
                    int pa = orbs_a[a].second.second;
                    int pa2 = nalphapi_[ha];

                    nso = nsopi_[ha];
                    nmo = nmopi_[ha];

                    Ca = Ca_->pointer(ha);
                    Ct = new double[nso];
                    eps = epsilon_a_->pointer(ha);

                    // Swap eigvals
                    epst = eps[pa];
                    eps[pa] = eps[pa2];
                    eps[pa2] = epst;

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);

                    delete[] Ct;

                    // Redo indexing
                    nalphapi_[ha]++;

                    outfile->Printf("   %8s: %4d%-4s -> %4d%-4s \n", "A  -> A ", pi + 1, ct.gamma(hi).symbol(), pa + 1,
                                    ct.gamma(ha).symbol());
                }
            } else if (!si && !sa) {
                // Beta->Beta
                int hi = orbs_b[i].second.first;
                int ha = orbs_b[a].second.first;

                if (hi == ha) {
                    // Same irrep
                    int pi = orbs_b[i].second.second;
                    int pa = orbs_b[a].second.second;

                    int nso = nsopi_[hi];
                    int nmo = nmopi_[hi];

                    double** Ca = Cb_->pointer(hi);
                    double* Ct = new double[nso];
                    double* eps = epsilon_b_->pointer(hi);
                    double epst;

                    // Swap eigvals
                    epst = eps[pi];
                    eps[pi] = eps[pa];
                    eps[pa] = epst;

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pa], nmo, &Ca[0][pi], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pa], nmo);

                    delete[] Ct;

                    outfile->Printf("   %8s: %4d%-4s -> %4d%-4s \n", "B  -> B ", pi + 1, ct.gamma(hi).symbol(), pa + 1,
                                    ct.gamma(ha).symbol());
                } else {
                    // Different irrep
                    // Occ -> Vir
                    int pi = orbs_b[i].second.second;
                    int pi2 = nbetapi_[hi] - 1;

                    int nso = nsopi_[hi];
                    int nmo = nmopi_[hi];

                    double** Ca = Cb_->pointer(hi);
                    double* Ct = new double[nso];
                    double* eps = epsilon_b_->pointer(hi);
                    double epst;

                    // Swap eigvals
                    epst = eps[pi];
                    eps[pi] = eps[pi2];
                    eps[pi2] = epst;

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);

                    delete[] Ct;

                    // Redo indexing
                    nbetapi_[hi]--;

                    // Vir -> Occ
                    int pa = orbs_b[a].second.second;
                    int pa2 = nbetapi_[ha];

                    nso = nsopi_[ha];
                    nmo = nmopi_[ha];

                    Ca = Cb_->pointer(ha);
                    Ct = new double[nso];
                    eps = epsilon_b_->pointer(ha);

                    // Swap eigvals
                    epst = eps[pa];
                    eps[pa] = eps[pa2];
                    eps[pa2] = epst;

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);

                    delete[] Ct;

                    // Redo indexing
                    nbetapi_[ha]++;

                    outfile->Printf("   %8s: %4d%-4s -> %4d%-4s \n", "B  -> B ", pi + 1, ct.gamma(hi).symbol(), pa + 1,
                                    ct.gamma(ha).symbol());
                }
            } else if (!si && sa) {
                // Beta->Alpha
                int hi = orbs_b[i].second.first;
                int ha = orbs_a[a].second.first;

                // Different irrep
                // Occ -> Vir
                int pi = orbs_b[i].second.second;
                int pi2 = nbetapi_[hi] - 1;

                int nso = nsopi_[hi];
                int nmo = nmopi_[hi];

                double** Ca = Cb_->pointer(hi);
                double* Ct = new double[nso];
                double* eps = epsilon_b_->pointer(hi);
                double epst;

                // Swap eigvals
                epst = eps[pi];
                eps[pi] = eps[pi2];
                eps[pi2] = epst;

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);

                delete[] Ct;

                // Redo indexing
                nbetapi_[hi]--;
                nbeta_--;

                // Vir -> Occ
                int pa = orbs_a[a].second.second;
                int pa2 = nalphapi_[ha];

                nso = nsopi_[ha];
                nmo = nmopi_[ha];

                Ca = Ca_->pointer(ha);
                Ct = new double[nso];
                eps = epsilon_a_->pointer(ha);

                // Swap eigvals
                epst = eps[pa];
                eps[pa] = eps[pa2];
                eps[pa2] = epst;

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);

                delete[] Ct;

                // Redo indexing
                nalphapi_[ha]++;
                nalpha_++;

                outfile->Printf("   %8s: %4d%-4s -> %4d%-4s \n", "B  -> A ", pi + 1, ct.gamma(hi).symbol(), pa + 1,
                                ct.gamma(ha).symbol());
            } else if (si && !sa) {
                // Alpha->Beta
                int hi = orbs_a[i].second.first;
                int ha = orbs_b[a].second.first;

                // Different irrep
                // Occ -> Vir
                int pi = orbs_a[i].second.second;
                int pi2 = nalphapi_[hi] - 1;

                int nso = nsopi_[hi];
                int nmo = nmopi_[hi];

                double** Ca = Ca_->pointer(hi);
                double* Ct = new double[nso];
                double* eps = epsilon_a_->pointer(hi);
                double epst;

                // Swap eigvals
                epst = eps[pi];
                eps[pi] = eps[pi2];
                eps[pi2] = epst;

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);

                delete[] Ct;

                // Redo indexing
                nalphapi_[hi]--;
                nalpha_--;

                // Vir -> Occ
                int pa = orbs_b[a].second.second;
                int pa2 = nbetapi_[ha];

                nso = nsopi_[ha];
                nmo = nmopi_[ha];

                Ca = Cb_->pointer(ha);
                Ct = new double[nso];
                eps = epsilon_b_->pointer(ha);

                // Swap eigvals
                epst = eps[pa];
                eps[pa] = eps[pa2];
                eps[pa2] = epst;

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);

                delete[] Ct;

                // Redo indexing
                nbetapi_[ha]++;
                nbeta_++;

                outfile->Printf("   %8s: %4d%-4s -> %4d%-4s \n", "A  -> B ", pi + 1, ct.gamma(hi).symbol(), pa + 1,
                                ct.gamma(ha).symbol());
            }
            if (nalpha_ < nbeta_)
                throw PSIEXCEPTION("PSI::MOM_start: Nbeta ends up being less than Nalpha, this is not supported");

        } else if (options_.get_str("REFERENCE") == "ROHF" || options_.get_str("REFERENCE") == "ROKS") {
            // For ROHF/ROKS, orbital position determines spin channel (sign is ignored):
            // DOCC->SOCC: beta electron, SOCC->VIR: alpha electron, DOCC->VIR: both.

            int hi = orbs_a[i].second.first;
            int ha = orbs_a[a].second.first;
            int pi = orbs_a[i].second.second;
            int pa = orbs_a[a].second.second;

            auto in_docc = [&](int h, int p) { return p < nbetapi_[h]; };
            auto in_socc = [&](int h, int p) { return p >= nbetapi_[h] && p < nalphapi_[h]; };
            auto in_vir  = [&](int h, int p) { return p >= nalphapi_[h]; };

            auto swap_mo = [&](SharedMatrix C, SharedVector eps, int h, int p, int q) {
                if (p == q) return;
                int nso = nsopi_[h];
                int nmo = nmopi_[h];
                if (nso == 0 || nmo == 0) return;
                double** Cp = C->pointer(h);
                double* ep = eps->pointer(h);
                auto* Ct = new double[nso];
                double epst = ep[p];
                ep[p] = ep[q];
                ep[q] = epst;
                C_DCOPY(nso, &Cp[0][p], nmo, Ct, 1);
                C_DCOPY(nso, &Cp[0][q], nmo, &Cp[0][p], nmo);
                C_DCOPY(nso, Ct, 1, &Cp[0][q], nmo);
                delete[] Ct;
            };

            enum Case { DOCC_TO_SOCC, SOCC_TO_VIR, DOCC_TO_VIR };
            Case excitation;
            const char* tag = "";
            if (in_docc(hi, pi) && in_socc(ha, pa)) {
                excitation = DOCC_TO_SOCC;
                tag = "D  -> S ";
            } else if (in_socc(hi, pi) && in_vir(ha, pa)) {
                excitation = SOCC_TO_VIR;
                tag = "S  -> V ";
            } else if (in_docc(hi, pi) && in_vir(ha, pa)) {
                excitation = DOCC_TO_VIR;
                tag = "D  -> V ";
            } else {
                throw PSIEXCEPTION(
                    "SCF::MOM_start: ROHF excitation must be DOCC->SOCC, SOCC->VIR, or DOCC->VIR");
            }

            if (excitation == DOCC_TO_SOCC) {
                if (hi == ha) {
                    swap_mo(Ca_, epsilon_a_, hi, pi, pa);
                } else {
                    swap_mo(Ca_, epsilon_a_, hi, pi, nbetapi_[hi] - 1);
                    nbetapi_[hi]--;
                    nbeta_--;
                    swap_mo(Ca_, epsilon_a_, ha, pa, nbetapi_[ha]);
                    nbetapi_[ha]++;
                    nbeta_++;
                }
            } else if (excitation == SOCC_TO_VIR) {
                if (hi == ha) {
                    swap_mo(Ca_, epsilon_a_, hi, pi, pa);
                } else {
                    swap_mo(Ca_, epsilon_a_, hi, pi, nalphapi_[hi] - 1);
                    nalphapi_[hi]--;
                    nalpha_--;
                    swap_mo(Ca_, epsilon_a_, ha, pa, nalphapi_[ha]);
                    nalphapi_[ha]++;
                    nalpha_++;
                }
            } else {  // DOCC_TO_VIR
                if (hi == ha) {
                    swap_mo(Ca_, epsilon_a_, hi, pi, pa);
                } else {
                    // DOCC -> Vir via SOCC slot
                    swap_mo(Ca_, epsilon_a_, hi, pi, nbetapi_[hi] - 1);
                    nbetapi_[hi]--;
                    nbeta_--;
                    swap_mo(Ca_, epsilon_a_, hi, nbetapi_[hi], nalphapi_[hi] - 1);
                    nalphapi_[hi]--;
                    nalpha_--;
                    // Vir -> DOCC via SOCC slot
                    swap_mo(Ca_, epsilon_a_, ha, pa, nalphapi_[ha]);
                    nalphapi_[ha]++;
                    nalpha_++;
                    swap_mo(Ca_, epsilon_a_, ha, nalphapi_[ha] - 1, nbetapi_[ha]);
                    nbetapi_[ha]++;
                    nbeta_++;
                }
            }

            outfile->Printf("   %8s: %4d%-4s -> %4d%-4s \n", tag, pi + 1, ct.gamma(hi).symbol(), pa + 1,
                            ct.gamma(ha).symbol());

            if (nalpha_ < nbeta_)
                throw PSIEXCEPTION("PSI::MOM_start: Nbeta ends up being less than Nalpha, this is not supported");
        }
    }
    Ca_old_->copy(Ca_);
    Cb_old_->copy(Cb_);

    outfile->Printf("\n                        Total Energy        Delta E      Density RMS\n\n");
}
void HF::MOM() {
    // ROHF/ROKS: track DOCC and SOCC via separate overlap pools so the DOCC/SOCC
    // boundary is preserved when orbitals are near-degenerate.
    const bool is_rohf_type =
        (options_.get_str("REFERENCE") == "ROHF" || options_.get_str("REFERENCE") == "ROKS");

    // Alpha
    for (int h = 0; h < nirrep_; h++) {
        // Indexing
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int nalpha = nalphapi_[h];
        int nbeta = nbetapi_[h];

        if (nso == 0 || nmo == 0 || nalpha == 0) continue;

        double** Cold = Ca_old_->pointer(h);
        double** Cnew = Ca_->pointer(h);
        double* eps = epsilon_a_->pointer(h);
        double** S = S_->pointer(h);

        if (is_rohf_type && nbeta > 0 && nalpha > nbeta) {
            int nsocc = nalpha - nbeta;
            int nvir = nmo - nalpha;

            // P[a,i] = (C_new^T S C_old)[a,i]; p_docc[a] = sum_{i<nbeta} P[a,i]^2,
            // p_socc[a] = sum_{nbeta<=i<nalpha} P[a,i]^2.
            double** SC = block_matrix(nso, nmo);
            C_DGEMM('N', 'N', nso, nmo, nso, 1.0, S[0], nso, Cold[0], nmo, 0.0, SC[0], nmo);
            double** P = block_matrix(nmo, nmo);
            C_DGEMM('T', 'N', nmo, nmo, nso, 1.0, Cnew[0], nmo, SC[0], nmo, 0.0, P[0], nmo);

            auto* p_docc = new double[nmo];
            auto* p_socc = new double[nmo];
            for (int a = 0; a < nmo; a++) {
                double sd = 0.0, ss = 0.0;
                for (int j = 0; j < nbeta; j++) sd += P[a][j] * P[a][j];
                for (int j = nbeta; j < nalpha; j++) ss += P[a][j] * P[a][j];
                p_docc[a] = sd;
                p_socc[a] = ss;
            }
            free_block(SC);
            free_block(P);

            // Top nbeta by |p_docc| -> DOCC
            std::vector<std::pair<double, int> > pvec_docc(nmo);
            for (int a = 0; a < nmo; a++) pvec_docc[a] = std::make_pair(std::fabs(p_docc[a]), a);
            sort(pvec_docc.begin(), pvec_docc.end(), std::greater<std::pair<double, int> >());

            std::vector<bool> used(nmo, false);
            std::vector<int> docc_idx(nbeta);
            for (int a = 0; a < nbeta; a++) {
                docc_idx[a] = pvec_docc[a].second;
                used[pvec_docc[a].second] = true;
            }

            // Top nsocc by |p_socc| among remaining -> SOCC
            std::vector<std::pair<double, int> > pvec_socc;
            pvec_socc.reserve(nmo - nbeta);
            for (int a = 0; a < nmo; a++) {
                if (!used[a]) pvec_socc.push_back(std::make_pair(std::fabs(p_socc[a]), a));
            }
            sort(pvec_socc.begin(), pvec_socc.end(), std::greater<std::pair<double, int> >());

            std::vector<int> socc_idx(nsocc);
            for (int a = 0; a < nsocc; a++) {
                socc_idx[a] = pvec_socc[a].second;
                used[pvec_socc[a].second] = true;
            }

            std::vector<int> vir_idx;
            vir_idx.reserve(nvir);
            for (int a = 0; a < nmo; a++) {
                if (!used[a]) vir_idx.push_back(a);
            }

            // Order each group by energy
            std::vector<std::pair<double, int> > docc_vec(nbeta);
            for (int a = 0; a < nbeta; a++) docc_vec[a] = std::make_pair(eps[docc_idx[a]], docc_idx[a]);
            sort(docc_vec.begin(), docc_vec.end());

            std::vector<std::pair<double, int> > socc_vec(nsocc);
            for (int a = 0; a < nsocc; a++) socc_vec[a] = std::make_pair(eps[socc_idx[a]], socc_idx[a]);
            sort(socc_vec.begin(), socc_vec.end());

            std::vector<std::pair<double, int> > vir_vec(nvir);
            for (int a = 0; a < nvir; a++) vir_vec[a] = std::make_pair(eps[vir_idx[a]], vir_idx[a]);
            sort(vir_vec.begin(), vir_vec.end());

            double** Ct = block_matrix(nso, nmo);
            memcpy(static_cast<void*>(Ct[0]), static_cast<void*>(Cnew[0]), sizeof(double) * nso * nmo);

            for (int a = 0; a < nbeta; a++) {
                eps[a] = docc_vec[a].first;
                C_DCOPY(nso, &Ct[0][docc_vec[a].second], nmo, &Cnew[0][a], nmo);
            }
            for (int a = 0; a < nsocc; a++) {
                eps[a + nbeta] = socc_vec[a].first;
                C_DCOPY(nso, &Ct[0][socc_vec[a].second], nmo, &Cnew[0][a + nbeta], nmo);
            }
            for (int a = 0; a < nvir; a++) {
                eps[a + nalpha] = vir_vec[a].first;
                C_DCOPY(nso, &Ct[0][vir_vec[a].second], nmo, &Cnew[0][a + nalpha], nmo);
            }

            free_block(Ct);
            delete[] p_docc;
            delete[] p_socc;
            continue;
        }

        auto* c = new double[nso];
        auto* d = new double[nso];
        auto* p = new double[nmo];

        memset(static_cast<void*>(c), '\0', sizeof(double) * nso);

        for (int i = 0; i < nalpha; i++) C_DAXPY(nso, 1.0, &Cold[0][i], nmo, c, 1);

        C_DGEMV('N', nso, nso, 1.0, S[0], nso, c, 1, 0.0, d, 1);
        C_DGEMV('T', nso, nmo, 1.0, Cnew[0], nmo, d, 1, 0.0, p, 1);

        // outfile->Printf("  P_a:\n");
        // for (int a = 0; a < nmo; a++)
        //    outfile->Printf("   a = %3d: %14.10f\n", a + 1, p[a]);

        // Find the largest contributions
        std::vector<std::pair<double, int> > pvec;
        pvec.resize(nmo);
        for (int a = 0; a < nmo; a++) pvec[a] = std::make_pair(std::fabs(p[a]), a);
        sort(pvec.begin(), pvec.end(), std::greater<std::pair<double, int> >());

        // outfile->Printf("  P_a sorted:\n");
        // for (int a = 0; a < nmo; a++)
        //    outfile->Printf("   a = %3d: Index = %3d, %14.10f\n", a + 1, pvec[a].second, pvec[a].first);

        // Now order the mos in each group
        std::vector<std::pair<double, int> > occvec;
        occvec.resize(nalpha);
        for (int a = 0; a < nalpha; a++) occvec[a] = std::make_pair(eps[pvec[a].second], pvec[a].second);
        sort(occvec.begin(), occvec.end());

        // outfile->Printf("  P_a_occ sorted:\n");
        // for (int a = 0; a < nalpha; a++)
        //    outfile->Printf("   a = %3d: Index = %3d, %14.10f\n", a + 1, occvec[a].second, occvec[a].first);

        std::vector<std::pair<double, int> > virvec;
        virvec.resize(nmo - nalpha);
        for (int a = 0; a < nmo - nalpha; a++)
            virvec[a] = std::make_pair(eps[pvec[a + nalpha].second], pvec[a + nalpha].second);
        sort(virvec.begin(), virvec.end());

        // outfile->Printf("  P_a_vir sorted:\n");
        // for (int a = 0; a < nmo - nalpha; a++)
        //    outfile->Printf("   a = %3d: Index = %3d, %14.10f\n", a + 1, virvec[a].second, virvec[a].first);

        double** Ct = block_matrix(nso, nmo);

        // Use Cold and p as a buffer
        memcpy(static_cast<void*>(Ct[0]), static_cast<void*>(Cnew[0]), sizeof(double) * nso * nmo);
        memcpy(static_cast<void*>(p), static_cast<void*>(eps), sizeof(double) * nmo);

        for (int a = 0; a < nalpha; a++) {
            eps[a] = occvec[a].first;
            C_DCOPY(nso, &Ct[0][occvec[a].second], nmo, &Cnew[0][a], nmo);
        }

        for (int a = 0; a < nmo - nalpha; a++) {
            eps[a + nalpha] = virvec[a].first;
            C_DCOPY(nso, &Ct[0][virvec[a].second], nmo, &Cnew[0][a + nalpha], nmo);
        }

        free_block(Ct);

        delete[] c;
        delete[] d;
        delete[] p;
    }

    // Beta
    if (same_a_b_orbs()) return;

    for (int h = 0; h < nirrep_; h++) {
        // Indexing
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int nbeta = nbetapi_[h];

        if (nso == 0 || nmo == 0 || nbeta == 0) continue;

        double** Cold = Cb_old_->pointer(h);
        double** Cnew = Cb_->pointer(h);
        double* eps = epsilon_b_->pointer(h);
        double** S = S_->pointer(h);

        auto* c = new double[nso];
        auto* d = new double[nso];
        auto* p = new double[nmo];

        memset(static_cast<void*>(c), '\0', sizeof(double) * nso);

        for (int i = 0; i < nbeta; i++) C_DAXPY(nso, 1.0, &Cold[0][i], nmo, c, 1);

        C_DGEMV('N', nso, nso, 1.0, S[0], nso, c, 1, 0.0, d, 1);
        C_DGEMV('T', nso, nmo, 1.0, Cnew[0], nmo, d, 1, 0.0, p, 1);

        // outfile->Printf("  P_a:\n");
        // for (int a = 0; a < nmo; a++)
        //    outfile->Printf("   a = %3d: %14.10f\n", a + 1, p[a]);

        // Find the largest contributions
        std::vector<std::pair<double, int> > pvec;
        pvec.resize(nmo);
        for (int a = 0; a < nmo; a++) pvec[a] = std::make_pair(std::fabs(p[a]), a);
        sort(pvec.begin(), pvec.end(), std::greater<std::pair<double, int> >());

        // outfile->Printf("  P_a sorted:\n");
        // for (int a = 0; a < nmo; a++)
        //    outfile->Printf("   a = %3d: Index = %3d, %14.10f\n", a + 1, pvec[a].second, pvec[a].first);

        // Now order the mos in each group
        std::vector<std::pair<double, int> > occvec;
        occvec.resize(nbeta);
        for (int a = 0; a < nbeta; a++) occvec[a] = std::make_pair(eps[pvec[a].second], pvec[a].second);
        sort(occvec.begin(), occvec.end());

        // outfile->Printf("  P_a_occ sorted:\n");
        // for (int a = 0; a < nbeta; a++)
        //    outfile->Printf("   a = %3d: Index = %3d, %14.10f\n", a + 1, occvec[a].second, occvec[a].first);

        std::vector<std::pair<double, int> > virvec;
        virvec.resize(nmo - nbeta);
        for (int a = 0; a < nmo - nbeta; a++)
            virvec[a] = std::make_pair(eps[pvec[a + nbeta].second], pvec[a + nbeta].second);
        sort(virvec.begin(), virvec.end());

        // outfile->Printf("  P_a_vir sorted:\n");
        // for (int a = 0; a < nmo - nalpha; a++)
        //    outfile->Printf("   a = %3d: Index = %3d, %14.10f\n", a + 1, virvec[a].second, virvec[a].first);

        double** Ct = block_matrix(nso, nmo);

        // Use Cold and p as a buffer
        memcpy(static_cast<void*>(Ct[0]), static_cast<void*>(Cnew[0]), sizeof(double) * nso * nmo);
        memcpy(static_cast<void*>(p), static_cast<void*>(eps), sizeof(double) * nmo);

        for (int a = 0; a < nbeta; a++) {
            eps[a] = occvec[a].first;
            C_DCOPY(nso, &Ct[0][occvec[a].second], nmo, &Cnew[0][a], nmo);
        }

        for (int a = 0; a < nmo - nbeta; a++) {
            eps[a + nbeta] = virvec[a].first;
            C_DCOPY(nso, &Ct[0][virvec[a].second], nmo, &Cnew[0][a + nbeta], nmo);
        }

        free_block(Ct);

        delete[] c;
        delete[] d;
        delete[] p;
    }
    Cb_old_->copy(Cb_);
}
}  // namespace scf
}  // namespace psi
