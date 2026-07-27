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

#include "cghf.h"

#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/v.h"

#include <Einsums/TensorBase/TensorBase.hpp>
#include <Einsums/TensorBase/Common.hpp>
#include <Einsums/Tensor/Tensor.hpp>

namespace psi {
namespace scf {

std::shared_ptr<BlockTensorView<double, 2>> share_matrix_to_einsums(SharedMatrix M) {
    // Heap-allocated holder keeps both the source Matrix and the view tensor alive; the
    // returned shared_ptr aliases the holder's refcount while pointing at the tensor.
    struct Holder {
        SharedMatrix matrix;
        BlockTensorView<double, 2> tensor;
    };
    auto holder = std::make_shared<Holder>();
    holder->matrix = M;
    for (size_t h = 0; h < M->nirrep(); h++) {
        auto dim = einsums::Dim<2>{M->rowspi(h), M->colspi(h)};
        einsums::TensorView<double, 2> t{M->get_pointer(h), dim};
        holder->tensor.push_block(t);
    }
    return std::shared_ptr<BlockTensorView<double, 2>>(holder, &holder->tensor);
}

CGHF::CGHF(std::shared_ptr<ComplexWavefunction> ref_wfn, std::shared_ptr<SuperFunctional> func)
    : CGHF(ref_wfn, func, Process::environment.options, PSIO::shared_object()) {}

CGHF::CGHF(std::shared_ptr<ComplexWavefunction> ref_wfn, std::shared_ptr<SuperFunctional> func, Options& options,
           std::shared_ptr<PSIO> psio)
    : ComplexWavefunction(options), BaseHF(func) {
    // Analogous to HF::HF's shallow_copy(ref_wfn); full ComplexWavefunction::shallow_copy comes later.
    molecule_ = ref_wfn->molecule();
    basisset_ = ref_wfn->basisset();
    psio_ = psio;
    // Options-only base ctor skipped init; run it now that molecule/basis are set.
    ComplexWavefunction::common_init();
    common_init();
}

CGHF::~CGHF() {}

void CGHF::common_init() {
    BaseHF::common_init(options_, module_, molecule_, dipole_field_strength_);
    name_ = "CGHF";
    // Prefer BaseHF::nelectron_ over BaseWavefunction's; copy from ComplexWavefunction.
    nelectron_ = nelec_;

    if (molecule_->schoenflies_symbol() != "c1") {
        throw PSIEXCEPTION("CGHF currently supports only C1 symmetry. Set symmetry c1 in the molecule block.");
    }

    // DFT stuff (would typically go in subclass_init)
    setup_potential();
}

void CGHF::setup_potential() {
    if (functional_->needs_xc()) {
        throw PSIEXCEPTION("DFT functionals are not supported with CGHF.");
    } else {
        potential_ = nullptr;
    }
}

void CGHF::set_jk(std::shared_ptr<BaseJK> jk) {
    // Cheap basis check
    int jk_nbf = jk->basisset()->nbf();
    int hf_nbf = basisset_->nbf();
    if (hf_nbf != jk_nbf) {
        throw PSIEXCEPTION("Tried setting a JK object whos number of basis functions does not match HF's!");
    }

    jk_ = jk;
}

std::shared_ptr<BaseJK> CGHF::build_jk(size_t memory) const {
    // auto jk = ComplexJK::build_JK(get_basisset("ORBITAL"), get_basisset("DF_BASIS_SCF"), Process::environment.options, false, memory);
    // return jk;
    return jk_;
}

}  // namespace scf
}  // namespace psi
