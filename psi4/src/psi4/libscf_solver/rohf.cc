/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <utility>
#include <vector>

#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/matrix.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"

#include "rohf.h"

namespace psi {
namespace scf {

ROHF::ROHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func)
    : HF(ref_wfn, func, Process::environment.options, PSIO::shared_object()) {
    common_init();
}

ROHF::ROHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func, Options& options,
           std::shared_ptr<PSIO> psio)
    : HF(ref_wfn, func, options, psio) {
    common_init();
}

ROHF::~ROHF() {}

void ROHF::common_init() {
    name_ = "ROHF";
    Fa_ = SharedMatrix(factory_->create_matrix("Alpha Fock Matrix"));
    Fb_ = SharedMatrix(factory_->create_matrix("Beta Fock Matrix"));
    moFeff_ = SharedMatrix(factory_->create_matrix("F effective (MO basis)"));
    soFeff_ = SharedMatrix(factory_->create_matrix("F effective (orthogonalized SO basis)"));
    Ct_ = SharedMatrix(factory_->create_matrix("Orthogonalized Molecular orbitals"));
    Ca_ = SharedMatrix(factory_->create_matrix("MO coefficients (C)"));
    Cb_ = Ca_;
    Da_ = SharedMatrix(factory_->create_matrix("SCF alpha density"));
    Db_ = SharedMatrix(factory_->create_matrix("SCF beta density"));
    Lagrangian_ = SharedMatrix(factory_->create_matrix("Lagrangian matrix"));
    Ka_ = SharedMatrix(factory_->create_matrix("K alpha"));
    Kb_ = SharedMatrix(factory_->create_matrix("K beta"));
    Ga_ = SharedMatrix(factory_->create_matrix("G alpha"));
    Gb_ = SharedMatrix(factory_->create_matrix("G beta"));
    Dt_ = SharedMatrix(factory_->create_matrix("Total SCF density"));
    Dt_old_ = SharedMatrix(factory_->create_matrix("Old total SCF density"));
    Da_old_ = SharedMatrix(factory_->create_matrix("Old alpha SCF density"));
    Db_old_ = SharedMatrix(factory_->create_matrix("Old beta SCF density"));
    moFa_ = SharedMatrix(factory_->create_matrix("MO alpha Fock Matrix (MO basis)"));
    moFb_ = SharedMatrix(factory_->create_matrix("MO beta Fock Matrix (MO basis)"));

    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_a_->set_name("orbital energies");
    epsilon_b_ = epsilon_a_;
    same_a_b_dens_ = false;
    same_a_b_orbs_ = true;

    subclass_init();
}

void ROHF::format_guess() {
    // Need to build Ct_
    Ct_ = linalg::triplet(X_, S_, Ca_, true, false, false);
}

void ROHF::semicanonicalize() {
    // Quick sanity checks
    if (Fa_ == 0 || Fb_ == 0)
        throw PSIEXCEPTION("Wavefunction: Semicanonicalize called, but Fock matrices are not initialized.");
    if (Ca_ == 0 || Cb_ == 0)
        throw PSIEXCEPTION("Wavefunction: Semicanonicalize called, but orbitals are not initialized.");
    if (Ca_ != Cb_) throw PSIEXCEPTION("Wavefunction: Semicanonicalize called, but orbitals are not the same.");
    if (Fa_ == Fb_) throw PSIEXCEPTION("Wavefunction: Semicanonicalize called, but Fock matrices are the same.");
    if (epsilon_a_ != epsilon_b_)
        throw PSIEXCEPTION("Wavefunction: Semicanonicalize called, but orbital energies are not the same.");
    // Now, make space for the new orbitals
    Cb_ = SharedMatrix(Ca_->clone());
    epsilon_b_ = std::make_shared<Vector>(std::move(epsilon_b_->clone()));
    Ca_->set_name("Alpha semicanonical orbitals");
    Cb_->set_name("Beta semicanonical orbitals");
    epsilon_a_->set_name("Alpha semicanonical orbital energies");
    epsilon_b_->set_name("Beta semicanonical orbital energies");
    same_a_b_orbs_ = false;

    SharedMatrix Crohf(Ca_->clone());
    SharedMatrix evecs;
    SharedVector evals;

    // Transform the Fock matrix to the MO basis
    SharedMatrix moFa(Fa_->clone());
    SharedMatrix moFb(Fb_->clone());
    moFa->transform(Ca_);
    moFb->transform(Ca_);

    // Pick out occ-occ, and vir-vir subsets of the Fock matrices
    Dimension dim_zero(nirrep_);
    auto aoccpi = nalphapi_;
    auto boccpi = nbetapi_;
    auto avirpi = nmopi_ - aoccpi;
    auto bvirpi = nmopi_ - boccpi;
    Slice aocc_slice(dim_zero, aoccpi);
    Slice bocc_slice(dim_zero, boccpi);
    Slice avir_slice(aoccpi, nmopi_);
    Slice bvir_slice(boccpi, nmopi_);
    auto aFOO = moFa->get_block(aocc_slice, aocc_slice);
    auto aFVV = moFa->get_block(avir_slice, avir_slice);
    auto bFOO = moFb->get_block(bocc_slice, bocc_slice);
    auto bFVV = moFb->get_block(bvir_slice, bvir_slice);

    // Canonicalize the Alpha occ-occ block
    evecs = std::make_shared<Matrix>(aoccpi, aoccpi);
    evals = std::make_shared<Vector>(aoccpi);
    aFOO->diagonalize(evecs, evals);
    for (int h = 0; h < nirrep_; ++h) {
        auto pC = Crohf->pointer(h);
        auto pCa = Ca_->pointer(h);
        auto pR = evecs->pointer(h);
        if (aoccpi[h]) {
            C_DGEMM('n', 'n', nsopi_[h], aoccpi[h], aoccpi[h], 1.0, pC[0], nmopi_[h], pR[0], aoccpi[h], 0.0, pCa[0],
                    nmopi_[h]);
            for (int p = 0; p < aoccpi[h]; ++p) {
                double epsilon = evals->get(h, p);
                epsilon_a_->set(h, p, epsilon);
                for (int q = 0; q < aoccpi[h]; ++q) {
                    moFa->set(h, p, q, p == q ? epsilon : 0.0);
                }
            }
        }
    }
    // Canonicalize the Alpha vir-vir block
    evecs = std::make_shared<Matrix>(avirpi, avirpi);
    evals = std::make_shared<Vector>(avirpi);
    aFVV->diagonalize(evecs, evals);
    for (int h = 0; h < nirrep_; ++h) {
        auto pC = Crohf->pointer(h);
        auto pCa = Ca_->pointer(h);
        auto pR = evecs->pointer(h);
        if (avirpi[h]) {
            C_DGEMM('n', 'n', nsopi_[h], avirpi[h], avirpi[h], 1.0, &(pC[0][aoccpi[h]]), nmopi_[h], pR[0], avirpi[h],
                    0.0, &(pCa[0][aoccpi[h]]), nmopi_[h]);
            for (int p = 0; p < avirpi[h]; ++p) {
                double epsilon = evals->get(h, p);
                epsilon_a_->set(h, p + aoccpi[h], epsilon);
                for (int q = 0; q < avirpi[h]; ++q) {
                    moFa->set(h, p + aoccpi[h], q + aoccpi[h], p == q ? epsilon : 0.0);
                }
            }
        }
    }
    // Canonicalize the Beta occ-occ block
    evecs = std::make_shared<Matrix>(boccpi, boccpi);
    evals = std::make_shared<Vector>(boccpi);
    bFOO->diagonalize(evecs, evals);
    for (int h = 0; h < nirrep_; ++h) {
        auto pC = Crohf->pointer(h);
        auto pCb = Cb_->pointer(h);
        auto pR = evecs->pointer(h);
        if (boccpi[h]) {
            C_DGEMM('n', 'n', nsopi_[h], boccpi[h], boccpi[h], 1.0, pC[0], nmopi_[h], pR[0], boccpi[h], 0.0, pCb[0],
                    nmopi_[h]);
            for (int p = 0; p < boccpi[h]; ++p) {
                double epsilon = evals->get(h, p);
                epsilon_b_->set(h, p, epsilon);
                for (int q = 0; q < boccpi[h]; ++q) {
                    moFb->set(h, p, q, p == q ? epsilon : 0.0);
                }
            }
        }
    }
    // Canonicalize the Beta vir-vir block
    evecs = std::make_shared<Matrix>(bvirpi, bvirpi);
    evals = std::make_shared<Vector>(bvirpi);
    bFVV->diagonalize(evecs, evals);
    for (int h = 0; h < nirrep_; ++h) {
        auto pC = Crohf->pointer(h);
        auto pCb = Cb_->pointer(h);
        auto pR = evecs->pointer(h);
        if (bvirpi[h]) {
            C_DGEMM('n', 'n', nsopi_[h], bvirpi[h], bvirpi[h], 1.0, &(pC[0][boccpi[h]]), nmopi_[h], pR[0], bvirpi[h],
                    0.0, &(pCb[0][boccpi[h]]), nmopi_[h]);
            for (int p = 0; p < bvirpi[h]; ++p) {
                double epsilon = evals->get(h, p);
                epsilon_b_->set(h, p + boccpi[h], epsilon);
                for (int q = 0; q < bvirpi[h]; ++q) {
                    moFb->set(h, p + boccpi[h], q + boccpi[h], p == q ? epsilon : 0.0);
                }
            }
        }
    }
    // Given the invariance w.r.t. occ-occ rotations, the existing Fa and Fb matrices
    // are still valid, because the densities used to construct them in the old basis
    // are equivalent to those in the new basis.  If the user forward transforms them
    // using the semicanonical orbitals, the correct semicanonical basis Fa and Fb
    // will be obtained.
}

void ROHF::finalize() {
    // Form Lagrangian
    //
    // In HF, the Lagrangian is Xpi = Fpi. For RHF and UHF, this reduces
    // to Xii = ei. In the AO basis (where we want it), Xmn = Cmi ei Cni.
    // For ROHF, the effective Fock matrix is diagonal, not Fa and Fb (as
    // in UHF). So we need to form the Lagrangian as: Xmn = Cmp Fpi Cni.
    //
    // --EGH
    //
    // Let's build the MO Lagrangian in moFeff_
    // ...I'm assuming these are square
    moFeff_->zero();
    moFa_->transform(Fa_, Ca_);
    moFb_->transform(Fb_, Ca_);
    for (int h = 0; h < nirrep_; ++h) {
        for (int m = 0; m < moFeff_->rowdim(h); ++m) {
            double tval;
            for (int i = 0; i < nbetapi_[h]; ++i) {
                tval = moFa_->get(h, m, i);
                tval += moFb_->get(h, m, i);
                moFeff_->set(h, m, i, tval);
            }
            for (int i = nbetapi_[h]; i < nalphapi_[h]; ++i) {
                tval = moFa_->get(h, m, i);
                moFeff_->set(h, m, i, tval);
            }
        }
    }
    Lagrangian_->back_transform(moFeff_, Ca_);

    moFeff_.reset();
    Ka_.reset();
    Kb_.reset();
    Ga_.reset();
    Gb_.reset();
    Da_old_.reset();
    Db_old_.reset();
    Dt_old_.reset();
    Dt_.reset();
    moFa_.reset();
    moFb_.reset();

    HF::finalize();
}

void ROHF::save_density_and_energy() {
    Da_old_->copy(Da_);
    Db_old_->copy(Db_);
    Dt_old_->copy(Dt_);
}

void ROHF::form_initial_F() {
    // Form the initial Fock matrix, closed and open variants
    Fa_->copy(H_);
    Fa_->add(Ga_);
    for (const auto& Vext : external_potentials_) {
        Fa_->add(Vext);
    }
    Fb_->copy(Fa_);

    if (debug_) {
        outfile->Printf("Initial alpha Fock matrix:\n");
        Fa_->print("outfile");
        outfile->Printf("Initial beta Fock matrix:\n");
        Fb_->print("outfile");
    }
}

void ROHF::form_F() {
    // Start by constructing the standard Fa and Fb matrices encountered in UHF
    Fa_->copy(H_);
    Fa_->add(Ga_);
    for (const auto& Vext : external_potentials_) {
        Fa_->add(Vext);
    }

    Fb_->copy(H_);
    Fb_->add(Gb_);
    for (const auto& Vext : external_potentials_) {
        Fb_->add(Vext);
    }

    moFa_->transform(Fa_, Ca_);
    moFb_->transform(Fb_, Ca_);

    /*
     * Fo = open-shell fock matrix = 0.5 Fa
     * Fc = closed-shell fock matrix = 0.5 (Fa + Fb)
     *
     * Therefore
     *
     * 2(Fc-Fo) = Fb
     * 2Fo = Fa
     *
     * Form the effective Fock matrix, too
     * The effective Fock matrix has the following structure
     *          |  closed     open    virtual
     *  ----------------------------------------
     *  closed  |    Fc     2(Fc-Fo)    Fc
     *  open    | 2(Fc-Fo)     Fc      2Fo
     *  virtual |    Fc       2Fo       Fc
     */
    moFeff_->copy(moFa_);
    moFeff_->add(moFb_);
    moFeff_->scale(0.5);
    for (int h = 0; h < nirrep_; ++h) {
        for (int i = nbetapi_[h]; i < nalphapi_[h]; ++i) {
            // Set the open/closed portion
            for (int j = 0; j < nbetapi_[h]; ++j) {
                double val = moFb_->get(h, i, j);
                moFeff_->set(h, i, j, val);
                moFeff_->set(h, j, i, val);
            }
            // Set the open/virtual portion
            for (int j = nalphapi_[h]; j < nmopi_[h]; ++j) {
                double val = moFa_->get(h, i, j);
                moFeff_->set(h, i, j, val);
                moFeff_->set(h, j, i, val);
            }
        }
    }

    // Form the orthogonalized SO basis moFeff matrix, for use in DIIS
    soFeff_->back_transform(moFeff_, Ct_);

    if (debug_) {
        Fa_->print();
        Fb_->print();
        moFa_->print();
        moFb_->print();
        moFeff_->print();
        soFeff_->print();
    }
}

void ROHF::form_C(double shift) {
    if (shift == 0.0) {
        soFeff_->diagonalize(Ct_, epsilon_a_);
    } else {
        // Implementation of the shifting scheme from M.F. Guest &
        // V. R. Saunders (1974), "On methods for converging open-shell
        // Hartree-Fock wave-functions", Molecular Physics, 28:3,
        // 819-828, DOI: 10.1080/00268977400102171. Virtuals are shifted
        // up by shift, open-shell orbitals are shifted up by shift/2.

        // Shifted Fock operator
        auto shifted_F = SharedMatrix(factory_->create_matrix("F"));

        // Zero dimensions
        Dimension dim_zero(nirrep_);

        // Shift open-shell orbitals by 0.5*shift
        auto ospi = nalphapi_ - nbetapi_;
        auto Cos = Ct_->get_block({dim_zero, nmopi_}, {nbetapi_, nbetapi_ + ospi});
        Cos->set_name("Cos");
        shifted_F->gemm(false, true, 0.5 * shift, Cos, Cos, 0.0);

        // Shift virtuals by shift
        auto virpi = nmopi_ - nalphapi_;
        auto Cvir = Ct_->get_block({dim_zero, nmopi_}, {nalphapi_, nalphapi_ + virpi});
        Cvir->set_name("Cvir");
        shifted_F->gemm(false, true, shift, Cvir, Cvir, 1.0);

        // Add in the Fock itself and diagonalize
        shifted_F->add(soFeff_);
        shifted_F->diagonalize(Ct_, epsilon_a_);
    }
    // Form C = XC'
    Ca_->gemm(false, false, 1.0, X_, Ct_, 0.0);

    find_occupation();

    if (debug_) {
        Ca_->print("outfile");
        outfile->Printf("In ROHF::form_C:\n");
        Ct_->eivprint(epsilon_a_);
    }
}

void ROHF::prepare_canonical_orthogonalization() {
    // Some matrix size changes if we canonical orthogonalization
    Ct_->init(nmopi_, nmopi_);
    moFa_->init(nmopi_, nmopi_);
    moFb_->init(nmopi_, nmopi_);
    moFeff_->init(nmopi_, nmopi_);
    soFeff_->init(nmopi_, nmopi_);  // This is in the "Orthogonalized SO" basis
}
void ROHF::form_initial_C() {
    // In ROHF the creation of the C matrix depends on the previous iteration's C
    // matrix. Here we use Fa to generate the first C, where Fa was set by guess()
    // to either H or the GWH Hamiltonian.

    // Form F' = X'FX for canonical orthogonalization
    auto diag_F_temp = linalg::triplet(X_, Fa_, X_, true, false, false);

    // Form C' = eig(F')
    diag_F_temp->diagonalize(Ct_, epsilon_a_);

    // Form C = XC'
    Ca_->gemm(false, false, 1.0, X_, Ct_, 0.0);

    find_occupation();

    if (debug_) {
        Ca_->print("outfile");
        outfile->Printf("In ROHF::form_initial_C:\n");
        Ct_->eivprint(epsilon_a_);
    }
}

void ROHF::form_D() {
    Da_->zero();
    Db_->zero();

    for (int h = 0; h < nirrep_; ++h) {
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int na = nalphapi_[h];
        int nb = nbetapi_[h];

        if (nso == 0 || nmo == 0) continue;

        auto Ca = Ca_->pointer(h);
        auto Da = Da_->pointer(h);
        auto Db = Db_->pointer(h);

        C_DGEMM('N', 'T', nso, nso, na, 1.0, Ca[0], nmo, Ca[0], nmo, 0.0, Da[0], nso);
        C_DGEMM('N', 'T', nso, nso, nb, 1.0, Ca[0], nmo, Ca[0], nmo, 0.0, Db[0], nso);
    }

    Dt_->copy(Da_);
    Dt_->add(Db_);

    if (debug_) {
        outfile->Printf("in ROHF::form_D:\n");
        Da_->print();
        Db_->print();
    }
}

double ROHF::compute_initial_E() { return nuclearrep_ + Dt_->vector_dot(H_); }

double ROHF::compute_E() {
    double one_electron_E = Dt_->vector_dot(H_);
    double kinetic_E = Dt_->vector_dot(T_);
    double two_electron_E = 0.5 * (Da_->vector_dot(Ga_) + Db_->vector_dot(Gb_));

    energies_["Nuclear"] = nuclearrep_;
    energies_["Kinetic"] = kinetic_E;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = two_electron_E;
    energies_["XC"] = 0.0;
    energies_["VV10_E"] = 0.0;
    energies_["-D"] = 0.0;

    double Eelec = one_electron_E + two_electron_E;
    double Etotal = nuclearrep_ + Eelec;
    return Etotal;
}

void ROHF::Hx(SharedMatrix x, SharedMatrix ret) {
    if (functional_->needs_xc()) {
        throw PSIEXCEPTION("SCF: Cannot yet compute DFT Hessian-vector prodcuts.\n");
    }
    // Index reference
    // left = IAJB + IAjb, right = iajb + iaJB
    // i = docc, a = socc, p = pure virtual
    // o = docc + socc, v = socc + vir

    // Spaces
    auto occpi = nalphapi_;
    auto virpi = nmopi_ - nbetapi_;
    auto pvir = nmopi_ - nalphapi_;
    auto docc = doccpi();
    auto socc = soccpi();

    if (ret->rowspi() != occpi) {
        throw PSIEXCEPTION("ROHF:Hx First dimension of rotation matrix is not correct.");
    }
    if (ret->colspi() != virpi) {
        throw PSIEXCEPTION("ROHF:Hx Second dimension of rotation matrix is not correct.");
    }

    // => Effective one electron part <= //
    auto Hx_left = std::make_shared<Matrix>("Partial Hx tensor left", ret->rowspi(), ret->colspi());
    auto Hx_right = std::make_shared<Matrix>("Partial Hx tensor right", ret->rowspi(), ret->colspi());

    // Passing these guys is annoying, pretty cheap to rebuild
    auto dim_zero = Dimension(nirrep_, "Zero Dim");

    auto Cocc = Ca_->get_block({dim_zero, nsopi_}, {dim_zero, ret->rowspi()});
    Cocc->set_name("Cocc");
    auto Cvir = Ca_->get_block({dim_zero, nsopi_}, {docc, docc + ret->colspi()});
    Cvir->set_name("Cvir");

    for (size_t h = 0; h < nirrep_; h++) {
        if (!occpi[h] || !virpi[h]) continue;

        auto leftp = Hx_left->pointer(h);
        auto rightp = Hx_right->pointer(h);
        auto xp = x->pointer(h);
        auto Fap = moFa_->pointer(h);
        auto Fbp = moFb_->pointer(h);

        // left_ov += 0.5 * x_op Fa_pv
        C_DGEMM('N', 'N', occpi[h], virpi[h], pvir[h], 0.5, (xp[0] + socc[h]), virpi[h],
                (Fap[occpi[h]] + docc[h]), nmopi_[h], 0.0, leftp[0], virpi[h]);

        // left_ov -= Fa_oo x_ov
        C_DGEMM('N', 'N', occpi[h], virpi[h], occpi[h], -0.5, Fap[0], nmopi_[h], xp[0], virpi[h], 1.0, leftp[0],
                virpi[h]);

        // right_ov += 0.5 * x_ov Fb_vv
        C_DGEMM('N', 'N', occpi[h], virpi[h], virpi[h], 0.5, xp[0], virpi[h], (Fbp[docc[h]] + docc[h]), nmopi_[h],
                0.0, rightp[0], virpi[h]);

        // right_ov -= Fb_oi x_iv
        C_DGEMM('N', 'N', occpi[h], virpi[h], docc[h], -0.5, Fbp[0], nmopi_[h], xp[0], virpi[h], 1.0, rightp[0],
                virpi[h]);
        if (socc[h]) {
            // Socc terms
            // left_av += 0.5 * x_oa.T Fb_ov
            C_DGEMM('T', 'N', socc[h], virpi[h], occpi[h], 0.5, xp[0], virpi[h], (Fbp[0] + docc[h]), nmopi_[h],
                    1.0, leftp[docc[h]], virpi[h]);

            // right_oa += 0.5 * Fb_op x_ap.T
            C_DGEMM('N', 'T', occpi[h], socc[h], pvir[h], 0.5, (Fbp[0] + occpi[h]), nmopi_[h],
                    (xp[docc[h]] + socc[h]), virpi[h], 1.0, rightp[0], virpi[h]);
        }
    }

    // => Two electron part <= //
    auto& Cl = jk_->C_left();
    auto& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    // If scf_type is DF we can do some extra JK voodo
    if ((options_.get_str("SCF_TYPE").find("DF") != std::string::npos) || (options_.get_str("SCF_TYPE") == "CD")) {
        auto Cdocc = Ca_->get_block({dim_zero, nsopi_}, {dim_zero, docc});
        Cdocc->set_name("Cdocc");

        auto Csocc = Ca_->get_block({dim_zero, nsopi_}, {docc, nalphapi_});
        Csocc->set_name("Csocc");

        auto Cr_i = std::make_shared<Matrix>("Cright for docc", nsopi_, docc);
        auto Cr_a = std::make_shared<Matrix>("Cright for socc", nsopi_, socc);
        auto Cl_a = std::make_shared<Matrix>("Cleft for socc", nsopi_, socc);

        for (size_t h = 0; h < nirrep_; h++) {
            if (!nsopi_[h]) continue;

            auto xp = x->pointer(h);
            auto Cp = Ca_->pointer(h);
            auto Cr_ip = Cr_i->pointer(h);
            auto Cr_ap = Cr_a->pointer(h);
            auto Cl_ap = Cl_a->pointer(h);

            if (docc[h] && pvir[h]) {
                // C_lambda,i = C_lambda,p x_ip.T
                C_DGEMM('N', 'T', nsopi_[h], docc[h], pvir[h], 1.0, Cp[0] + occpi[h], nmopi_[h],
                        (xp[0] + socc[h]), virpi[h], 0.0, Cr_ip[0], docc[h]);
            }

            if (socc[h] && pvir[h]) {
                // C_lambda,a = C_lambda,p x_ap.T
                C_DGEMM('N', 'T', nsopi_[h], socc[h], pvir[h], 1.0, Cp[0] + occpi[h], nmopi_[h],
                        (xp[docc[h]] + socc[h]), virpi[h], 0.0, Cr_ap[0], socc[h]);
            }

            if (socc[h] && docc[h]) {
                // C_lambda,a = C_lambda,i x_ia
                C_DGEMM('N', 'N', nsopi_[h], socc[h], docc[h], 1.0, Cp[0], nmopi_[h], xp[0], virpi[h], 0.0,
                        Cl_ap[0], socc[h]);
            }
        }

        Cl.push_back(Cdocc);
        Cl.push_back(Csocc);
        Cl.push_back(Cl_a);

        Cr.push_back(Cr_i);
        Cr.push_back(Cr_a);
        Cr.push_back(Csocc);

        jk_->compute();

        // Just in case someone only clears out Cleft and gets very strange errors
        Cl.clear();
        Cr.clear();
        Cdocc.reset();
        Csocc.reset();
        Cr_i.reset();
        Cr_a.reset();
        Cl_a.reset();

        const auto& J = jk_->J();
        const auto& K = jk_->K();

        // Collect left terms
        J[0]->scale(2.0);
        J[0]->add(J[1]);
        J[0]->add(J[2]);
        J[1]->copy(J[0]);

        K[2]->add(K[0]);
        K[0]->add(K[1]);

        K[0]->scale(0.5);
        J[0]->subtract(K[0]);
        J[0]->subtract(K[0]->transpose());

        // Collect right terms
        K[2]->scale(0.5);
        J[1]->subtract(K[2]);
        J[1]->subtract(K[2]->transpose());

        // Transform to MO basis and add to exsisting
        auto half_trans = std::make_shared<Matrix>("half_trans temp space", occpi, nsopi_);

        half_trans->gemm(true, false, 1.0, Cocc, J[0], 0.0);
        Hx_left->gemm(false, false, 1.0, half_trans, Cvir, 1.0);

        half_trans->gemm(true, false, 1.0, Cocc, J[1], 0.0);
        Hx_right->gemm(false, false, 1.0, half_trans, Cvir, 1.0);

        half_trans.reset();

    } else {
        auto Cdocc = Ca_->get_block({dim_zero, nsopi_}, {dim_zero, docc});
        Cdocc->set_name("Cdocc");

        Cl.push_back(Cocc);
        Cl.push_back(Cdocc);

        auto Cr_a = std::make_shared<Matrix>("Cright for alpha", nsopi_, occpi);
        auto Cr_b = std::make_shared<Matrix>("Cright for beta", nsopi_, docc);

        for (size_t h = 0; h < nirrep_; h++) {
            if (!nsopi_[h]) continue;

            auto xp = x->pointer(h);
            auto Cp = Ca_->pointer(h);
            auto Cr_ap = Cr_a->pointer(h);
            auto Cr_bp = Cr_b->pointer(h);

            if (occpi[h] && pvir[h]) {
                // C_lambda,o = C_lambda,p x_op.T
                C_DGEMM('N', 'T', nsopi_[h], occpi[h], pvir[h], 1.0, Cp[0] + occpi[h], nmopi_[h], (xp[0] + socc[h]),
                        virpi[h], 0.0, Cr_ap[0], occpi[h]);
            }

            if (docc[h] && virpi[h]) {
                // C_lambda,i = C_lambda,p x_ip.T
                C_DGEMM('N', 'T', nsopi_[h], docc[h], virpi[h], 1.0, Cp[0] + docc[h], nmopi_[h], xp[0], virpi[h],
                        0.0, Cr_bp[0], docc[h]);
            }
        }

        Cr.push_back(Cr_a);
        Cr.push_back(Cr_b);

        jk_->compute();

        // Just in case someone only clears out Cleft and gets very strange errors
        Cl.clear();
        Cr.clear();
        Cr_a.reset();
        Cr_b.reset();
        Cdocc.reset();

        const auto& J = jk_->J();
        const auto& K = jk_->K();

        // Collect left terms
        J[0]->add(J[1]);
        J[1]->copy(J[0]);

        K[0]->scale(0.5);
        J[0]->subtract(K[0]);
        J[0]->subtract(K[0]->transpose());

        // Collect right terms
        K[1]->scale(0.5);
        J[1]->subtract(K[1]);
        J[1]->subtract(K[1]->transpose());

        // Transform to MO basis and add to exsisting
        auto half_trans = std::make_shared<Matrix>("half_trans temp space", occpi, nsopi_);

        half_trans->gemm(true, false, 1.0, Cocc, J[0], 0.0);
        Hx_left->gemm(false, false, 1.0, half_trans, Cvir, 1.0);

        half_trans->gemm(true, false, 1.0, Cocc, J[1], 0.0);
        Hx_right->gemm(false, false, 1.0, half_trans, Cvir, 1.0);

        half_trans.reset();
    }

    // Zero out socc-socc terms
    for (size_t h = 0; h < nirrep_; h++) {
        if (!occpi[h] || !virpi[h]) continue;

        auto leftp = Hx_left->pointer(h);
        auto rightp = Hx_right->pointer(h);

        for (size_t i = 0; i < socc[h]; i++) {
            for (size_t j = 0; j < occpi[h]; j++) {
                leftp[j][i] = 0.0;
            }
            for (size_t j = 0; j < virpi[h]; j++) {
                rightp[docc[h] + i][j] = 0.0;
            }
        }
    }

    ret->zero();
    ret->add(Hx_left);
    ret->add(Hx_right);
    ret->scale(4.0);

    Hx_left.reset();
    Hx_right.reset();
    Cocc.reset();
    Cvir.reset();
}

void ROHF::damping_update(double damping_percentage) {
    Da_->scale(1.0 - damping_percentage);
    Da_->axpy(damping_percentage, Da_old_);
    Db_->scale(1.0 - damping_percentage);
    Db_->axpy(damping_percentage, Db_old_);
    Dt_->copy(Da_);
    Dt_->add(Db_);
}

int ROHF::soscf_update(double soscf_conv, int soscf_min_iter, int soscf_max_iter, bool soscf_print) {
    std::time_t start, stop;
    start = std::time(nullptr);

    // => Build gradient and preconditioner <= //

    // Only the inact-act, inact-vir, and act-vir rotations are non-redundant
    auto dim_zero = Dimension(nirrep_, "Zero Dim");
    auto occpi = nalphapi_;
    auto virpi = nmopi_ - nbetapi_;
    auto docc = doccpi();
    auto socc = soccpi();

    auto Gradient = moFeff_->get_block({dim_zero, occpi}, {docc, nmopi_});
    Gradient->scale(-4.0);
    auto Precon = std::make_shared<Matrix>("Precon", occpi, virpi);

    for (size_t h = 0; h < nirrep_; h++) {
        if (!occpi[h] || !virpi[h]) continue;

        auto preconp = Precon->pointer(h);
        auto gradientp = Gradient->pointer(h);
        auto feffp = moFeff_->pointer(h);

        // Precon
        for (size_t i = 0; i < occpi[h]; i++) {
            for (size_t j = 0; j < virpi[h]; j++) {
                // Should be scaled but 3.5 works better (4 * 0.875 or GWH coef)
                preconp[i][j] = -3.5 * (feffp[i][i] - feffp[docc[h] + j][docc[h] + j]);
            }
        }

        // Divide docc-socc and socc-vir by two
        // Fix socc-socc term to avoid divide by zero errors
        for (size_t i = 0; i < socc[h]; i++) {
            for (size_t j = 0; j < occpi[h]; j++) {
                gradientp[j][i] *= 0.5;
                preconp[j][i] *= 0.5;
            }
            for (size_t j = 0; j < virpi[h]; j++) {
                gradientp[docc[h] + i][j] *= 0.5;
                preconp[docc[h] + i][j] *= 0.5;
            }
            for (size_t j = 0; j < socc[h]; j++) {
                preconp[docc[h] + i][j] = 1.0;
                gradientp[docc[h] + i][j] = 0.0;
            }
        }
    }

    // Make sure the MO gradient is reasonably small
    if (Gradient->absmax() > 0.3) {
        if (print_ > 1) {
            outfile->Printf("    Gradient element too large for SOSCF, using DIIS.\n");
        }
        return 0;
    }

    if (soscf_print) {
        outfile->Printf("\n");
        outfile->Printf("    ==> SOROHF Iterations <==\n");
        outfile->Printf("    Maxiter     = %11d\n", soscf_max_iter);
        outfile->Printf("    Miniter     = %11d\n", soscf_min_iter);
        outfile->Printf("    Convergence = %11.3E\n", soscf_conv);
        outfile->Printf("    ---------------------------------------\n");
        outfile->Printf("    %-4s   %11s     %10s\n", "Iter", "Residual RMS", "Time [s]");
        outfile->Printf("    ---------------------------------------\n");
    }

    // => Initial CG guess <= //
    auto x = Gradient->clone();
    x->set_name("Current ROHF CG Guess");
    x->apply_denominator(Precon);

    // Calc hessian vector product, find residual and conditioned residual
    auto r = Gradient->clone();
    auto Ap = std::make_shared<Matrix>("Ap", occpi, virpi);
    Hx(x, Ap);
    r->subtract(Ap);

    // Print iteration 0 timings and rms
    double rconv = r->vector_dot(r);
    double grad_rms = Gradient->vector_dot(Gradient);
    if (grad_rms < 1.e-14) {
        grad_rms = 1.e-14;  // Prevent rel denom from being too small
    }
    double rms = std::sqrt(rconv / grad_rms);
    stop = std::time(nullptr);
    if (soscf_print) {
        outfile->Printf("    %-5s %11.3E %10ld\n", "Guess", rms, stop - start);
    }

    // Build new p and z vectors
    SharedMatrix z = r->clone();
    z->apply_denominator(Precon);
    SharedMatrix p = z->clone();

    // => CG iterations <= //
    int fock_builds = 1;
    for (int cg_iter = 1; cg_iter < soscf_max_iter; cg_iter++) {
        // Calc hessian vector product
        Hx(p, Ap);
        fock_builds += 1;

        // Find factors and scale
        double rzpre = r->vector_dot(z);
        double alpha = rzpre / p->vector_dot(Ap);
        if (std::isnan(alpha)) {
            outfile->Printf("ROHF::SOSCF Warning CG alpha is zero/nan. Stopping CG.\n");
            alpha = 0.0;
        }

        x->axpy(alpha, p);
        r->axpy(-alpha, Ap);

        // Get residual
        double rconv = r->sum_of_squares();
        double rms = std::sqrt(rconv / grad_rms);
        stop = std::time(nullptr);
        if (soscf_print) {
            outfile->Printf("    %-5d %11.3E %10ld\n", cg_iter, rms, stop - start);
        }

        // Check convergence
        if (((rms < soscf_conv) && (cg_iter >= soscf_min_iter)) || (alpha == 0.0)) {
            break;
        }

        // Update p and z
        z->copy(r);
        z->apply_denominator(Precon);

        double beta = r->vector_dot(z) / rzpre;

        p->scale(beta);
        p->add(z);
    }

    if (soscf_print) {
        outfile->Printf("    ---------------------------------------\n");
        outfile->Printf("\n");
    }

    // => Rotate orbitals <= //
    rotate_orbitals(Ca_, x);
    rotate_orbitals(Ct_, x);

    // // => Cleanup <= //
    Precon.reset();
    Gradient.reset();
    Ap.reset();
    z.reset();
    r.reset();
    p.reset();

    return fock_builds;
}

void ROHF::form_G() {
    Dimension dim_zero(nirrep_, "Zero Dim");

    auto& C = jk_->C_left();
    C.clear();

    // Push back docc orbitals
    SharedMatrix Cdocc = Ca_->get_block({dim_zero, nsopi_}, {dim_zero, nbetapi_});
    C.push_back(Cdocc);

    // Push back socc orbitals
    SharedMatrix Csocc = Ca_->get_block({dim_zero, nsopi_}, {nbetapi_, nalphapi_});
    C.push_back(Csocc);

    // Run the JK object
    jk_->compute();

    // Pull the J and K matrices off
    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();
    Ga_->copy(J[0]);
    Ga_->scale(2.0);
    Ga_->add(J[1]);

    Ka_->copy(K[0]);
    Ka_->add(K[1]);
    Kb_ = K[0];

    Gb_->copy(Ga_);
    Ga_->subtract(Ka_);
    Gb_->subtract(Kb_);
}

bool ROHF::stability_analysis() {
    if (functional_->needs_xc()) {
        throw PSIEXCEPTION("Stability analysis not yet supported for XC functionals.");
    }
    if (scf_type_ == "DF" || scf_type_ == "CD") {
        throw PSIEXCEPTION("Stability analysis has not been implemented for density fitted wavefunctions yet.");
    } else {
        auto nvirpi = nmopi_ - nalphapi_;
        auto zero = Dimension(nirrep_);
        auto docc = doccpi();
        auto socc = soccpi();
        auto navir = nmopi_ - nalphapi_;
        auto nbvir = nmopi_ - nbetapi_;

        auto FIJ = std::make_shared<Matrix>("Alpha occupied MO basis Fock matrix", nalphapi_, nalphapi_);
        auto Fij = std::make_shared<Matrix>("Beta occupied MO basis Fock matrix", nalphapi_, nalphapi_);
        auto FAB = std::make_shared<Matrix>("Alpha virtual MO basis Fock matrix", nbvir, nbvir);
        auto Fab = std::make_shared<Matrix>("Beta virtual MO basis Fock matrix", nbvir, nbvir);
        auto FIA = std::make_shared<Matrix>("Alpha occ-vir MO basis Fock matrix", nalphapi_, nbvir);
        auto Fia = std::make_shared<Matrix>("Beta occ-vir MO basis Fock matrix", nalphapi_, nbvir);

        SharedMatrix Cocc = Ca_->get_block({zero, nsopi_}, {zero, nalphapi_});
        std::vector<SharedMatrix> virandsoc;
        virandsoc.push_back(Ca_->get_block({zero, nsopi_}, {nalphapi_, nmopi_}));
        virandsoc.push_back(Ca_->get_block({zero, nsopi_}, {docc, nalphapi_}));
        SharedMatrix Cvir = linalg::horzcat(virandsoc);
        FIJ->transform(Fa_, Cocc);
        Fij->transform(Fb_, Cocc);

        FIA->transform(Cocc, Fa_, Cvir);
        Fia->transform(Cocc, Fb_, Cvir);

        FAB->transform(Fa_, Cvir);
        Fab->transform(Fb_, Cvir);

        // Zero out alpha socc terms corresponding to vir indices
        for (int h = 0; h < nirrep_; ++h) {
            for (int mo1 = navir[h]; mo1 < nbvir[h]; ++mo1) {
                for (int mo2 = navir[h]; mo2 < nbvir[h]; ++mo2) {
                    FAB->set(h, mo1, mo2, 0.0);
                }
            }
            for (int mo1 = 0; mo1 < nalphapi_[h]; ++mo1) {
                for (int mo2 = navir[h]; mo2 < nbvir[h]; ++mo2) {
                    FIA->set(h, mo1, mo2, 0.0);
                }
            }
        }
        // Zero out beta socc terms corresponding to occ indices
        for (int h = 0; h < nirrep_; ++h) {
            for (int mo1 = docc[h]; mo1 < nalphapi_[h]; ++mo1) {
                for (int mo2 = docc[h]; mo2 < nalphapi_[h]; ++mo2) {
                    Fij->set(h, mo1, mo2, 0.0);
                }
            }
            for (int mo1 = docc[h]; mo1 < nalphapi_[h]; ++mo1) {
                for (int mo2 = 0; mo2 < nbvir[h]; ++mo2) {
                    Fia->set(h, mo1, mo2, 0.0);
                }
            }
        }

        // Some indexing arrays to allow us to compare occ and vir spatial orbitals, below.
        int occ_offsets[8];
        occ_offsets[0] = 0;
        for (int h = 1; h < nirrep_; ++h) occ_offsets[h] = occ_offsets[h - 1] + docc[h - 1];

        int vir_offsets[8];
        vir_offsets[0] = docc[0] - navir[0];
        for (int h = 1; h < nirrep_; ++h) vir_offsets[h] = vir_offsets[h - 1] + docc[h] - navir[h];

        std::vector<std::shared_ptr<MOSpace> > spaces;
        spaces.push_back(MOSpace::occ);
        spaces.push_back(MOSpace::vir);
#define ID(x) ints.DPD_ID(x)
        IntegralTransform ints(shared_from_this(), spaces, IntegralTransform::TransformationType::Restricted,
                               IntegralTransform::OutputType::DPDOnly, IntegralTransform::MOOrdering::QTOrder,
                               IntegralTransform::FrozenOrbitals::None);
        ints.set_keep_dpd_so_ints(true);
        ints.transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
        ints.transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir);
        dpd_set_default(ints.get_dpd_id());
        dpdbuf4 Aaa, Aab, Aba, Abb, I, A;
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "MO Ints (OV|OV)");
        // A_IA_JB = (IA|JB)
        global_dpd_->buf4_scmcopy(&I, PSIF_LIBTRANS_DPD, "ROHF Hessian (IA|JB)", 1.0);
        // A_IA_jb = (IA|jb)
        global_dpd_->buf4_scmcopy(&I, PSIF_LIBTRANS_DPD, "ROHF Hessian (IA|jb)", 1.0);

        // A_IA_JB -= 0.5 (IB|JA)
        global_dpd_->buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, psrq, ID("[O,V]"), ID("[O,V]"), "ROHF Hessian (IA|JB)",
                                    -0.5);
        global_dpd_->buf4_close(&I);

        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>=O]+"), ID("[V>=V]+"), 0,
                               "MO Ints (OO|VV)");
        // A_IA_JB -= 0.5 (IJ|AB)
        global_dpd_->buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "ROHF Hessian (IA|JB)",
                                    -0.5);
        global_dpd_->buf4_close(&I);

        global_dpd_->buf4_init(&Aaa, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "ROHF Hessian (IA|JB)");

        // A_ia_jb = A_IA_JB
        global_dpd_->buf4_copy(&Aaa, PSIF_LIBTRANS_DPD, "ROHF Hessian (ia|jb)");
        for (int h = 0; h < Aaa.params->nirreps; ++h) {
            global_dpd_->buf4_mat_irrep_init(&Aaa, h);
            global_dpd_->buf4_mat_irrep_rd(&Aaa, h);
            for (int ia = 0; ia < Aaa.params->rowtot[h]; ++ia) {
                int iabs = Aaa.params->roworb[h][ia][0];
                int aabs = Aaa.params->roworb[h][ia][1];
                int isym = Aaa.params->psym[iabs];
                int asym = Aaa.params->qsym[aabs];
                int irel = iabs - Aaa.params->poff[isym];
                int arel = aabs - Aaa.params->qoff[asym];
                for (int jb = 0; jb < Aaa.params->coltot[h]; ++jb) {
                    int jabs = Aaa.params->colorb[h][jb][0];
                    int babs = Aaa.params->colorb[h][jb][1];
                    int jsym = Aaa.params->rsym[jabs];
                    int bsym = Aaa.params->ssym[babs];
                    int jrel = jabs - Aaa.params->roff[jsym];
                    int brel = babs - Aaa.params->soff[bsym];
                    double val = Aaa.matrix[h][ia][jb];
                    // A_IA_JB += 0.5 delta_IJ F_AB - 0.5 delta_AB F_IJ
                    if ((iabs == jabs) && (asym == bsym)) val += 0.5 * FAB->get(asym, arel, brel);
                    if ((aabs == babs) && (isym == jsym)) val -= 0.5 * FIJ->get(isym, irel, jrel);
                    // Zero out any socc-socc terms
                    if (arel >= nvirpi[asym] || brel >= nvirpi[bsym]) val = 0.0;
                    Aaa.matrix[h][ia][jb] = val;
                }
            }
            global_dpd_->buf4_mat_irrep_wrt(&Aaa, h);
        }
        // global_dpd_->buf4_sort(&Aaa, PSIF_LIBTRANS_DPD, qpsr, ID("[V,O]"), ID("[V,O]"), "tmp");
        global_dpd_->buf4_close(&Aaa);

        // global_dpd_->buf4_init(&Aaa, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[V,O]"),
        //              ID("[V,O]"), ID("[V,O]"), 0, "tmp");
        // global_dpd_->buf4_print(&Aaa, "outfile", 1);
        // global_dpd_->buf4_close(&Aaa);
        // exit(1);

        global_dpd_->buf4_init(&Abb, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "ROHF Hessian (ia|jb)");
        for (int h = 0; h < Abb.params->nirreps; ++h) {
            global_dpd_->buf4_mat_irrep_init(&Abb, h);
            global_dpd_->buf4_mat_irrep_rd(&Abb, h);
            for (int ia = 0; ia < Abb.params->rowtot[h]; ++ia) {
                int iabs = Abb.params->roworb[h][ia][0];
                int aabs = Abb.params->roworb[h][ia][1];
                int isym = Abb.params->psym[iabs];
                int asym = Abb.params->qsym[aabs];
                int irel = iabs - Abb.params->poff[isym];
                int arel = aabs - Abb.params->qoff[asym];
                for (int jb = 0; jb < Abb.params->coltot[h]; ++jb) {
                    int jabs = Abb.params->colorb[h][jb][0];
                    int babs = Abb.params->colorb[h][jb][1];
                    int jsym = Abb.params->rsym[jabs];
                    int bsym = Abb.params->ssym[babs];
                    int jrel = jabs - Abb.params->roff[jsym];
                    int brel = babs - Abb.params->soff[bsym];
                    double val = Abb.matrix[h][ia][jb];
                    // A_ia_jb += 0.5 delta_ij F_ab - 0.5 delta_ab F_ij
                    if ((iabs == jabs) && (asym == bsym)) {
                        val += 0.5 * Fab->get(asym, arel, brel);
                    }
                    if ((aabs == babs) && (isym == jsym)) val -= 0.5 * Fij->get(isym, irel, jrel);
                    // Zero out any socc-socc terms
                    if (irel >= docc[isym] || jrel >= docc[jsym]) val = 0.0;
                    Abb.matrix[h][ia][jb] = val;
                }
            }
            global_dpd_->buf4_mat_irrep_wrt(&Abb, h);
        }
        global_dpd_->buf4_close(&Abb);

        global_dpd_->buf4_init(&Aab, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "ROHF Hessian (IA|jb)");
        for (int h = 0; h < Aab.params->nirreps; ++h) {
            global_dpd_->buf4_mat_irrep_init(&Aab, h);
            global_dpd_->buf4_mat_irrep_rd(&Aab, h);
            for (int ia = 0; ia < Aab.params->rowtot[h]; ++ia) {
                int iabs = Aab.params->roworb[h][ia][0];
                int aabs = Aab.params->roworb[h][ia][1];
                int isym = Aab.params->psym[iabs];
                int asym = Aab.params->qsym[aabs];
                int irel = iabs - Aab.params->poff[isym];
                int arel = aabs - Aab.params->qoff[asym];
                for (int jb = 0; jb < Aab.params->coltot[h]; ++jb) {
                    int jabs = Aab.params->colorb[h][jb][0];
                    int babs = Aab.params->colorb[h][jb][1];
                    int jsym = Aab.params->rsym[jabs];
                    int bsym = Aab.params->ssym[babs];
                    int jrel = jabs - Aab.params->roff[jsym];
                    int brel = babs - Aab.params->soff[bsym];
                    // A_IA_jb += 0.5 delta_Ib F(beta)_jA
                    // The delta_Ib is actually comparing spatial orbitals, so we have to use
                    // offsets to make the comparison properly.
                    if (((irel + occ_offsets[isym]) == (brel + vir_offsets[bsym])) && (jsym == asym))
                        Aab.matrix[h][ia][jb] += 0.5 * Fia->get(jsym, jrel, arel);
                    // Zero out any socc-socc terms
                    if (jrel >= docc[jsym] || arel >= navir[asym]) Aab.matrix[h][ia][jb] = 0.0;
                }
            }
            global_dpd_->buf4_mat_irrep_wrt(&Aab, h);
        }

        // A_ia_JB = A_IA_jb
        global_dpd_->buf4_sort(&Aab, PSIF_LIBTRANS_DPD, rspq, ID("[O,V]"), ID("[O,V]"), "ROHF Hessian (ia|JB)");
        global_dpd_->buf4_close(&Aab);

        // Alpha-Alpha
        global_dpd_->buf4_init(&Aaa, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "ROHF Hessian (IA|JB)");
        global_dpd_->buf4_copy(&Aaa, PSIF_LIBTRANS_DPD, "ROHF Hessian");
        global_dpd_->buf4_close(&Aaa);
        global_dpd_->buf4_init(&A, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "ROHF Hessian");
        // Alpha-Beta
        global_dpd_->buf4_init(&Aab, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "ROHF Hessian (IA|jb)");
        global_dpd_->buf4_axpy(&Aab, &A, 1.0);
        global_dpd_->buf4_close(&Aab);
        // Beta-Alpha
        global_dpd_->buf4_init(&Aba, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "ROHF Hessian (ia|JB)");
        global_dpd_->buf4_axpy(&Aba, &A, 1.0);
        global_dpd_->buf4_close(&Aba);
        // Beta-Beta
        global_dpd_->buf4_init(&Abb, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "ROHF Hessian (ia|jb)");
        global_dpd_->buf4_axpy(&Abb, &A, 1.0);
        global_dpd_->buf4_close(&Abb);
        global_dpd_->buf4_close(&A);

        /*
         *  Perform the stability analysis
         */
        std::vector<std::pair<double, int> > eval_sym;
        global_dpd_->buf4_init(&A, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                               "ROHF Hessian");

        // Allocate some space for the eigenvalues
        int nsave = options_.get_int("SOLVER_N_ROOT");
        std::vector<int> onevec(nirrep_, 1);
        std::vector<int> dimvec(nirrep_, nsave);
        Dimension ones(onevec);
        Dimension evalsdim(dimvec);
        auto stabvals = std::make_shared<Matrix>("Eigenvalues from ROHF stability calculation", evalsdim, ones);

        for (int h = 0; h < A.params->nirreps; ++h) {
            double** pEvals = stabvals->pointer(h);
            int npairs = A.params->rowtot[h];
            if (npairs == 0) continue;

            // Store the row indices, for convenience
            size_t rank = 0;
            double** U = block_matrix(npairs, npairs);
            for (int ia = 0; ia < npairs; ++ia) {
                int iabs = Aab.params->roworb[h][ia][0];
                int aabs = Aab.params->roworb[h][ia][1];
                int isym = Aab.params->psym[iabs];
                int asym = Aab.params->qsym[aabs];
                int irel = iabs - Aab.params->poff[isym];
                int arel = aabs - Aab.params->qoff[asym];
                if ((arel >= nvirpi[asym]) && (irel >= docc[isym])) {
                    U[ia][ia] = 0.0;
                } else {
                    U[ia][ia] = 1.0;
                    rank++;
                }
            }
            if (rank == 0) continue;
            int lastcol = npairs - 1;
            for (int ia = 0; ia < npairs; ++ia) {
                if (U[ia][ia] == 0.0) {
                    while (U[lastcol][lastcol] == 0.0) lastcol--;
                    if (lastcol > ia) {
                        U[lastcol][ia] = U[ia][lastcol] = 1.0;
                        U[lastcol][lastcol] = 0.0;
                    }
                }
            }

            global_dpd_->buf4_mat_irrep_init(&A, h);
            global_dpd_->buf4_mat_irrep_rd(&A, h);

            // Use the transformation matrix to rearrange the columns
            double** temp = block_matrix(npairs, npairs);
            C_DGEMM('n', 'n', npairs, npairs, npairs, 1.0, A.matrix[h][0], npairs, U[0], npairs, 0.0, temp[0], npairs);
            C_DGEMM('n', 'n', npairs, npairs, npairs, 1.0, U[0], npairs, temp[0], npairs, 0.0, A.matrix[h][0], npairs);
            free_block(temp);

            auto* evals = new double[rank];
            double** evecs = block_matrix(rank, rank);

            if (DSYEV_ascending(rank, A.matrix[h], evals, evecs) != 0){
                throw PSIEXCEPTION("DSYEV diagonalizer failed in ROHF stability check!");
            }
            global_dpd_->buf4_mat_irrep_close(&A, h);
            int mindim = rank < 15 ? rank : 15;
            for (int i = 0; i < mindim; i++) eval_sym.push_back(std::make_pair(evals[i], h));
            for (int i = 0; i < nsave; ++i) pEvals[i][0] = evals[i];

            free_block(evecs);
            delete[] evals;
        }
        global_dpd_->buf4_close(&A);
        Process::environment.arrays["SCF STABILITY EIGENVALUES"] = stabvals;

        outfile->Printf("    Lowest ROHF->ROHF stability eigenvalues:\n");
        print_stability_analysis(eval_sym);
        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }
    // FOLLOW not implemented yet for ROHF
    return false;
}

std::shared_ptr<ROHF> ROHF::c1_deep_copy(std::shared_ptr<BasisSet> basis) {
    auto wfn = Wavefunction::c1_deep_copy(basis);
    auto hf_wfn = std::make_shared<ROHF>(wfn, functional_, wfn->options(), wfn->psio());

    // now just have to copy the matrices that RHF initializes
    // include only those that are not temporary (some deleted in finalize())
    if (Ca_) {
        hf_wfn->Ca_ = Ca_subset("AO", "ALL");
        hf_wfn->Cb_ = hf_wfn->Ca_;
    }
    if (Da_) hf_wfn->Da_ = Da_subset("AO");
    if (Db_) hf_wfn->Db_ = Db_subset("AO");
    if (Fa_) hf_wfn->Fa_ = Fa_subset("AO");
    if (Fb_) hf_wfn->Fb_ = Fb_subset("AO");
    if (epsilon_a_) {
        hf_wfn->epsilon_a_ = epsilon_subset_helper(epsilon_a_, nalphapi_, "AO", "ALL");
        hf_wfn->epsilon_b_ = hf_wfn->epsilon_a_;
    }
    // H_ ans X_ reset in the HF constructor, copy them over here
    auto SO2AO = aotoso()->transpose();
    if (H_) hf_wfn->H_->remove_symmetry(H_, SO2AO);
    if (X_) hf_wfn->X_->remove_symmetry(X_, SO2AO);

    return hf_wfn;
}

void ROHF::compute_SAD_guess(bool natorb) {
    // Form the SAD guess
    HF::compute_SAD_guess(natorb);
    if (natorb) {
        // Need to build Ct
        Ct_ = linalg::triplet(X_, S_, Ca_);
    } else {
        // Form the total density matrix that's used in energy evaluation
        Dt_->copy(Da_);
        Dt_->add(Db_);
    }
}

void ROHF::setup_potential() {
    if (functional_->needs_xc()) {
        throw PSIEXCEPTION("ROHF: Cannot compute XC components!");
    }
}

}  // namespace scf
}  // namespace psi
