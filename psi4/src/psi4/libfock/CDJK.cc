/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/integral.h"
#include "psi4/lib3index/cholesky.h"
#include "psi4/libpsi4util/process.h"

#include "jk.h"

#include <sstream>
#include "psi4/libpsi4util/PsiOutStream.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {

CDJK::CDJK(std::shared_ptr<BasisSet> primary, double cholesky_tolerance)
    : DiskDFJK(primary, primary), cholesky_tolerance_(cholesky_tolerance) {}
CDJK::~CDJK() {}
void CDJK::initialize_JK_disk() { throw PsiException("Disk algorithm for CD JK not implemented.", __FILE__, __LINE__); }
size_t CDJK::memory_estimate() {
    // Size is unknown until actual evaluation
    size_t nbf = primary_->nbf();

    // Assume cholesky index is ~4x of nbf.
    return nbf * nbf * nbf * 4;
}

void CDJK::initialize_JK_core() {
    timer_on("CD: cholesky decomposition");
    auto integral = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);
    cderi_ = std::shared_ptr<TwoBodyAOInt>(integral->eri());
    int ntri = cderi_->function_pairs().size();
    /// If user asks to read integrals from disk, just read them from disk.
    /// Qmn is only storing upper triangle.
    /// Ugur needs ncholesky_ in NAUX (SCF), but it can also be read from disk
    if (df_ints_io_ == "LOAD") {
        psio_->open(unit_, PSIO_OPEN_OLD);
        psio_->read_entry(unit_, "length", (char*)&ncholesky_, sizeof(long int));
        Qmn_ = std::make_shared<Matrix>("Qmn (CD Integrals)", ncholesky_, ntri);
        double** Qmnp = Qmn_->pointer();
        psio_->read_entry(unit_, "(Q|mn) Integrals", (char*)Qmnp[0], sizeof(double) * ntri * ncholesky_);
        psio_->close(unit_, 1);
        Process::environment.globals["NAUX (SCF)"] = ncholesky_;
        timer_off("CD: cholesky decomposition");
        return;
    }

    /// If user does not want to read from disk, recompute the cholesky integrals
    auto Ch = std::make_shared<CholeskyERI>(cderi_, 0.0, cholesky_tolerance_,
                                            memory_);
    Ch->choleskify();
    ncholesky_ = Ch->Q();
    size_t three_memory = ncholesky_ * ntri;
    size_t nbf = primary_->nbf();

    /// Kinda silly to check for memory after you perform CD.
    /// Most likely redundant as cholesky also checks for memory.
    if (memory_ < ((size_t)sizeof(double) * three_memory + (size_t)sizeof(double) * ncholesky_ * nbf * nbf))
        throw PsiException("Not enough memory for CD.", __FILE__, __LINE__);

    std::shared_ptr<Matrix> L = Ch->L();
    double** Lp = L->pointer();
    timer_off("CD: cholesky decomposition");

    Qmn_ = std::make_shared<Matrix>("Qmn (CD Integrals)", ncholesky_, ntri);

    double** Qmnp = Qmn_->pointer();

    const std::vector<long int>& schwarz_fun_pairs = cderi_->function_pairs_to_dense();

    timer_on("CD: schwarz");
    for (size_t mu = 0; mu < nbf; mu++) {
        for (size_t nu = mu; nu < nbf; nu++) {
            if (schwarz_fun_pairs[nu * (nu + 1) / 2 + mu] < 0) continue;
            for (long int P = 0; P < ncholesky_; P++) {
                Qmnp[P][schwarz_fun_pairs[nu * (nu + 1) / 2 + mu]] = Lp[P][mu * nbf + nu];
            }
        }
    }
    timer_off("CD: schwarz");
    if (df_ints_io_ == "SAVE") {
        psio_->open(unit_, PSIO_OPEN_NEW);
        psio_->write_entry(unit_, "length", (char*)&ncholesky_, sizeof(long int));
        psio_->write_entry(unit_, "(Q|mn) Integrals", (char*)Qmnp[0], sizeof(double) * ntri * ncholesky_);
        psio_->close(unit_, 1);
        // stick ncholesky in process environment for other codes that may use the integrals
        // Not sure if this should really be here.  It is here because Ugur uses this option to get the number of
        // cholesky vectors.
        Process::environment.globals["NAUX (SCF)"] = ncholesky_;
    }
}
void CDJK::manage_JK_core() {
    for (int Q = 0; Q < ncholesky_; Q += max_rows_) {
        int naux = (ncholesky_ - Q <= max_rows_ ? ncholesky_ - Q : max_rows_);
        if (do_J_) {
            timer_on("JK: J");
            block_J(&Qmn_->pointer()[Q], naux);
            timer_off("JK: J");
        }
        if (do_K_) {
            timer_on("JK: K");
            block_K(&Qmn_->pointer()[Q], naux);
            timer_off("JK: K");
        }
    }
}
void CDJK::print_header() const {
    if (print_) {
        outfile->Printf("  ==> CDJK: Cholesky-decomposed J/K Matrices <==\n\n");

        outfile->Printf("    J tasked:             %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:             %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:            %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) {
            throw PsiException("no wk for scf_type cd.", __FILE__, __LINE__);
            // outfile->Printf( "    Omega:                %11.3E\n", omega_);
        }
        outfile->Printf("    OpenMP threads:       %11d\n", omp_nthread_);
        outfile->Printf("    Integrals threads:    %11d\n", df_ints_num_threads_);
        outfile->Printf("    Memory [MiB]:         %11ld\n", (memory_ * 8L) / (1024L * 1024L));
        outfile->Printf("    Algorithm:            %11s\n", (is_core_ ? "Core" : "Disk"));
        outfile->Printf("    Integral Cache:       %11s\n", df_ints_io_.c_str());
        outfile->Printf("    Schwarz Cutoff:       %11.0E\n", cutoff_);
        outfile->Printf("    Cholesky tolerance:   %11.2E\n", cholesky_tolerance_);
        outfile->Printf("    No. Cholesky vectors: %11li\n\n", ncholesky_);
    }
}
}  // namespace psi
