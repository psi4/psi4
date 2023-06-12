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

#include "jk.h"

#include "psi4/lib3index/3index.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/lib3index/cholesky.h"
#include "psi4/lib3index/dfhelper.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/integral.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#include <sstream>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {

template <class T>
void _set_dfjk_options(std::shared_ptr<T> jk, Options& options) {
    if (options["INTS_TOLERANCE"].has_changed()) jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
    if (options["PRINT"].has_changed()) jk->set_print(options.get_int("PRINT"));
    if (options["DEBUG"].has_changed()) jk->set_debug(options.get_int("DEBUG"));
    if (options["BENCH"].has_changed()) jk->set_bench(options.get_int("BENCH"));
    jk->set_condition(options.get_double("DF_FITTING_CONDITION"));
    if (options["DF_INTS_NUM_THREADS"].has_changed())
        jk->set_df_ints_num_threads(options.get_int("DF_INTS_NUM_THREADS"));
}

JK::JK(std::shared_ptr<BasisSet> primary) : primary_(primary) { common_init(); }
JK::~JK() {}
std::shared_ptr<JK> JK::build_JK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                                 Options& options, std::string jk_type) {

    // check if algorithm is composite
    std::array<std::string, 3> composite_algos = { "DFDIRJ", "COSX", "LINK" };
    bool is_composite = std::any_of(
      composite_algos.cbegin(),
      composite_algos.cend(),
      [&](std::string composite_algo) { return jk_type.find(composite_algo) != std::string::npos; }
    );

    // exit calculation if density screening is selected for incompatible JK algo
    bool do_density_screen = options.get_str("SCREENING") == "DENSITY";

    std::array<std::string, 3> can_do_density_screen = { "DIRECT", "DFDIRJ+LINK", "DFDIRJ" };
    bool is_compatible_density_screen = std::any_of(
        can_do_density_screen.cbegin(),
        can_do_density_screen.cend(),
        [&](std::string jk_algo) { return jk_type == jk_algo; }
    ); 

    std::array<std::string, 3> can_be_in_df_scf_guess = { "DIRECT", "MEM_DF", "DISK_DF" };
    bool is_in_df_scf_guess = std::any_of(
        can_be_in_df_scf_guess.cbegin(),
        can_be_in_df_scf_guess.cend(),
        [&](std::string jk_algo) { return jk_type == jk_algo; }
    ); 
    bool do_df_scf_guess = options.get_bool("DF_SCF_GUESS") && is_in_df_scf_guess;

    bool is_incompatible_density_screen = !(is_compatible_density_screen || do_df_scf_guess);
    
    if (do_density_screen && is_incompatible_density_screen) {
        std::string error_message = "SCREENING=DENSITY has not been implemented for ";
        error_message += (do_df_scf_guess) ? "DF_SCF_GUESS" : jk_type;
        error_message += ".";
   
        throw PSIEXCEPTION(error_message);
    }
    
    // exit calculation if no screening is selected for incompatible JK algo
    bool do_no_screen = options.get_str("SCREENING") == "NONE";
    
    std::array<std::string, 3> cant_do_no_screen = { "PK", "DISK_DF", "DIRECT" };
    bool is_incompatible_no_screen = std::any_of(
        cant_do_no_screen.cbegin(),
        cant_do_no_screen.cend(),
        [&](std::string jk_algo) { return jk_type == jk_algo; }
    ); 
    is_incompatible_no_screen |= is_composite; 
    
    if (do_no_screen && is_incompatible_no_screen) {
        std::string error_message = "SCREENING=NONE has not been implemented for ";
        error_message += jk_type;
        error_message += ".";
 
        throw PSIEXCEPTION(error_message);
    }

    // exit calculation for other incompatible JK + SCREENING combos 
    std::string screening_type = options.get_str("SCREENING");
    if (jk_type == "DFDIRJ+LINK" && ((screening_type == "SCHWARZ") || screening_type == "CSAM" )) {
        std::string error_message = "SCREENING=";
        error_message += screening_type;
        error_message += " has not been implemented for ";
        error_message += jk_type;
        error_message += ".";
 
        throw PSIEXCEPTION(error_message);
    }
   
    // Throw small DF warning
    if (jk_type == "DF") {
        outfile->Printf("\n  Warning: JK type 'DF' found in simple constructor, defaulting to DiskDFJK.\n");
        outfile->Printf("           Please use the build_JK(primary, auxiliary, options, do_wK, memory)\n");
        outfile->Printf("           constructor as DiskDFJK has non-optimal performance for many workloads.\n\n");
        jk_type = "DISK_DF";
    }

    if (jk_type == "CD") {
        auto jk = std::make_shared<CDJK>(primary, options, options.get_double("CHOLESKY_TOLERANCE"));

        if (options["INTS_TOLERANCE"].has_changed()) jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["SCREENING"].has_changed()) jk->set_csam(options.get_str("SCREENING") == "CSAM");
        if (options["PRINT"].has_changed()) jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed()) jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed()) jk->set_bench(options.get_int("BENCH"));
        if (options["DF_INTS_IO"].has_changed()) jk->set_df_ints_io(options.get_str("DF_INTS_IO"));
        jk->set_condition(options.get_double("DF_FITTING_CONDITION"));
        if (options["DF_INTS_NUM_THREADS"].has_changed())
            jk->set_df_ints_num_threads(options.get_int("DF_INTS_NUM_THREADS"));

        return jk;

    } else if (jk_type == "DISK_DF") {
        auto jk = std::make_shared<DiskDFJK>(primary, auxiliary, options);
        _set_dfjk_options<DiskDFJK>(jk, options);
        if (options["DF_INTS_IO"].has_changed()) jk->set_df_ints_io(options.get_str("DF_INTS_IO"));

        return jk;

    } else if (jk_type == "MEM_DF") {
        auto jk = std::make_shared<MemDFJK>(primary, auxiliary, options);
        // TODO: re-enable after fixing all bugs
        jk->set_wcombine(false);
        _set_dfjk_options<MemDFJK>(jk, options);
        if (options["WCOMBINE"].has_changed()) { jk->set_wcombine(options.get_bool("WCOMBINE")); }

        return jk;
    } else if (jk_type == "PK") {
        auto jk = std::make_shared<PKJK>(primary, options);

        if (options["INTS_TOLERANCE"].has_changed()) jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["SCREENING"].has_changed()) jk->set_csam(options.get_str("SCREENING") == "CSAM");
        if (options["PRINT"].has_changed()) jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed()) jk->set_debug(options.get_int("DEBUG"));

        return jk;

    } else if (jk_type == "OUT_OF_CORE") {
        auto jk = std::make_shared<DiskJK>(primary, options);

        if (options["INTS_TOLERANCE"].has_changed()) jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["SCREENING"].has_changed()) jk->set_csam(options.get_str("SCREENING") == "CSAM");
        if (options["PRINT"].has_changed()) jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed()) jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed()) jk->set_bench(options.get_int("BENCH"));

        return jk;

    } else if (jk_type == "DIRECT") {
        auto jk = std::make_shared<DirectJK>(primary, options);

        if (options["INTS_TOLERANCE"].has_changed()) jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["SCREENING"].has_changed()) jk->set_csam(options.get_str("SCREENING") == "CSAM");
        if (options["PRINT"].has_changed()) jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed()) jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed()) jk->set_bench(options.get_int("BENCH"));
        if (options["DF_INTS_NUM_THREADS"].has_changed())
            jk->set_df_ints_num_threads(options.get_int("DF_INTS_NUM_THREADS"));

        return jk;

    /// handle composite methods
    } else if (is_composite) {
        auto jk = std::make_shared<CompositeJK>(primary, auxiliary, options);

        //if (options["INTS_TOLERANCE"].has_changed()) jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["SCREENING"].has_changed()) jk->set_csam(options.get_str("SCREENING") == "CSAM");
        if (options["PRINT"].has_changed()) jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed()) jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed()) jk->set_bench(options.get_int("BENCH"));

        return jk;

    } else {
        std::stringstream message;
        message << "JK::build_JK: Unkown SCF Type '" << jk_type << "'" << std::endl;
        throw PSIEXCEPTION(message.str());
    }
}
std::shared_ptr<JK> JK::build_JK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                                 Options& options) {
    // if SCF_TYPE == DF, you are using the wrong constructor and get an error next constructor in
    return build_JK(primary, auxiliary, options, options.get_str("SCF_TYPE"));
}
std::shared_ptr<JK> JK::build_JK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                                 Options& options, bool do_wK, size_t doubles) {
    std::string jk_type = options.get_str("SCF_TYPE");
    if (jk_type == "DF") {
        // Build exact estimate via Schwarz metrics
        auto jk = build_JK(primary, auxiliary, options, "MEM_DF");
        jk->set_do_wK(do_wK);
        if (jk->memory_estimate() < doubles) {
            return jk;
        }
        jk.reset();

        // Use Disk DFJK
        return build_JK(primary, auxiliary, options, "DISK_DF");

    } else {  // otherwise it has already been set
        return build_JK(primary, auxiliary, options, options.get_str("SCF_TYPE"));
    }

    // I am not passing wK and doubles to the next constructor FIXME??
    // instead, I will let the already existing sets do their job
    // this requires do_wK and doubles to be passed here and set
}
SharedVector JK::iaia(SharedMatrix /*Ci*/, SharedMatrix /*Ca*/) {
    throw PSIEXCEPTION("JK: (ia|ia) integrals not implemented");
}

const std::unordered_map<std::string, std::vector<size_t> >& JK::computed_shells_per_iter() {
    return computed_shells_per_iter_;
}

const std::vector<size_t>& JK::computed_shells_per_iter(const std::string& n_let) {
    return computed_shells_per_iter_[n_let];
}

void JK::common_init() {
    print_ = 1;
    debug_ = 0;
    bench_ = 0;

    // 256 MB default
    memory_ = 32000000L;
    omp_nthread_ = 1;
#ifdef _OPENMP
    omp_nthread_ = Process::environment.get_n_threads();
#endif
    cutoff_ = 1.0E-12;
    do_csam_ = false;

    do_J_ = true;
    do_K_ = true;
    do_wK_ = false;
    wcombine_ = false;
    lr_symmetric_ = false;
    omega_ = 0.0;
    omega_alpha_ = 1.0;
    omega_beta_ = 0.0;

    num_computed_shells_ = 0L;
    computed_shells_per_iter_ = {};

    std::shared_ptr<IntegralFactory> integral =
        std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);
    auto pet = std::make_shared<PetiteList>(primary_, integral);
    AO2USO_ = SharedMatrix(pet->aotoso());
}
size_t JK::memory_overhead() const {
    size_t mem = 0L;
    int JKwKD_factor = 1;
    if (do_J_) JKwKD_factor++;
    if (do_K_) JKwKD_factor++;
    if (do_wK_) JKwKD_factor++;

    int C_factor = 1;
    if (!lr_symmetric_) C_factor++;

    // USO Quantities
    for (size_t N = 0; N < C_left_.size(); N++) {
        int symml = C_left_[N]->symmetry();
        for (int h = 0; h < C_left_[N]->nirrep(); h++) {
            int nbfl = C_left_[N]->rowspi()[h];
            int nbfr = C_right_[N]->rowspi()[h];
            int nocc = C_left_[N]->colspi()[symml ^ h];

            mem += C_factor * (size_t)nocc * (nbfl + nbfr) / 2L + JKwKD_factor * (size_t)nbfl * nbfr;
        }
    }

    // AO Copies
    if (C1() && C_left_.size() && C_left_[0]->nirrep() != 1) {
        int nbf = primary_->nbf();
        for (size_t N = 0; N < C_left_.size(); N++) {
            int nocc = 0;
            for (int h = 0; h < C_left_[N]->nirrep(); h++) {
                nocc += C_left_[N]->colspi()[h];
            }
            mem += C_factor * (size_t)nocc * nbf + JKwKD_factor * (size_t)nbf * nbf;
        }
    }

    return mem;
}
void JK::compute_D() {
    /// Make sure the memory is there
    bool same = true;
    if (C_left_.size() != D_.size()) {
        same = false;
    } else {
        for (size_t N = 0; N < D_.size(); N++) {
            if (D_[N]->symmetry() != (C_left_[N]->symmetry() ^ C_right_[N]->symmetry())) same = false;
        }
    }

    if (!same) {
        D_.clear();
        for (size_t N = 0; N < C_left_.size(); ++N) {
            std::stringstream s;
            s << "D " << N << " (SO)";
            D_.push_back(std::make_shared<Matrix>(s.str(), C_left_[N]->nirrep(), C_left_[N]->rowspi(),
                                                  C_right_[N]->rowspi(),
                                                  C_left_[N]->symmetry() ^ C_right_[N]->symmetry()));
        }
    }

    // Form the density, differs from dou
    for (size_t N = 0; N < D_.size(); ++N) {
        int symm = D_[N]->symmetry();
        D_[N]->zero();
        for (int h = 0; h < D_[N]->nirrep(); ++h) {
            int nsol = C_left_[N]->rowspi()[h ^ C_left_[N]->symmetry()];
            int nocc = C_left_[N]->colspi()[h];
            int nsor = C_right_[N]->rowspi()[h ^ symm];

            if (!nsol || !nsor || !nocc) continue;

            double** Dp = D_[N]->pointer(h ^ symm);
            double** Clp = C_left_[N]->pointer(h);
            double** Crp = C_right_[N]->pointer(h ^ symm);

            C_DGEMM('N', 'T', nsol, nsor, nocc, 1.0, Clp[0], nocc, Crp[0], nocc, 0.0, Dp[0], nsor);
        }
    }
}
void JK::allocate_JK() {
    // Allocate J/K in the case that the algorithm uses USOs, so AO2USO will not allocate.
    bool same = true;
    if (J_.size() != D_.size()) {
        same = false;
    } else {
        for (size_t N = 0; N < D_.size(); N++) {
            if (D_[N]->symmetry() != J_[N]->symmetry()) same = false;
        }
    }

    if (!same) {
        J_.clear();
        K_.clear();
        wK_.clear();
        for (size_t N = 0; N < D_.size() && do_J_; ++N) {
            std::stringstream s;
            s << "J " << N << " (SO)";
            J_.push_back(std::make_shared<Matrix>(s.str(), D_[N]->nirrep(), D_[N]->rowspi(), D_[N]->rowspi(),
                                                  D_[N]->symmetry()));
        }
        for (size_t N = 0; N < D_.size() && do_K_; ++N) {
            std::stringstream s;
            s << "K " << N << " (SO)";
            K_.push_back(std::make_shared<Matrix>(s.str(), D_[N]->nirrep(), D_[N]->rowspi(), D_[N]->rowspi(),
                                                  D_[N]->symmetry()));
        }
        for (size_t N = 0; N < D_.size() && do_wK_; ++N) {
            std::stringstream s;
            s << "wK " << N << " (SO)";
            wK_.push_back(std::make_shared<Matrix>(s.str(), D_[N]->nirrep(), D_[N]->rowspi(), D_[N]->rowspi(),
                                                   D_[N]->symmetry()));
        }
    }
}
void JK::USO2AO() {
    allocate_JK();

    // If C1, C_ao and D_ao are equal to C and D
    if (AO2USO_->nirrep() == 1) {
        C_left_ao_ = C_left_;
        C_right_ao_ = C_right_;
        D_ao_ = D_;
        J_ao_ = J_;
        K_ao_ = K_;
        wK_ao_ = wK_;
        return;
    }

    if (J_ao_.size() != D_.size()) {
        J_ao_.clear();
        K_ao_.clear();
        wK_ao_.clear();
        D_ao_.clear();
        int nao = AO2USO_->rowspi()[0];
        for (size_t N = 0; N < D_.size() && do_J_; ++N) {
            std::stringstream s;
            s << "J " << N << " (AO)";
            J_ao_.push_back(std::make_shared<Matrix>(s.str(), nao, nao));
        }
        for (size_t N = 0; N < D_.size() && do_K_; ++N) {
            std::stringstream s;
            s << "K " << N << " (AO)";
            K_ao_.push_back(std::make_shared<Matrix>(s.str(), nao, nao));
        }
        for (size_t N = 0; N < D_.size() && do_wK_; ++N) {
            std::stringstream s;
            s << "wK " << N << " (AO)";
            wK_ao_.push_back(std::make_shared<Matrix>(s.str(), nao, nao));
        }
        for (size_t N = 0; N < D_.size(); ++N) {
            std::stringstream s;
            s << "D " << N << " (AO)";
            D_ao_.push_back(std::make_shared<Matrix>(s.str(), nao, nao));
        }
    }

    // Always reallocate C matrices, the occupations are tricky
    C_left_ao_.clear();
    C_right_ao_.clear();
    for (size_t N = 0; N < D_.size(); ++N) {
        std::stringstream s;
        s << "C Left " << N << " (AO)";
        int ncol = C_left_[N]->colspi().sum();
        C_left_ao_.push_back(std::make_shared<Matrix>(s.str(), AO2USO_->rowspi()[0], ncol));
    }
    for (size_t N = 0; (N < D_.size()) && (!lr_symmetric_); ++N) {
        std::stringstream s;
        s << "C Right " << N << " (AO)";
        int ncol = C_right_[N]->colspi().sum();
        C_right_ao_.push_back(std::make_shared<Matrix>(s.str(), AO2USO_->rowspi()[0], ncol));
    }

    // Alias pointers if lr_symmetric_
    if (lr_symmetric_) {
        C_right_ao_ = C_left_ao_;
    }

    // Transform D
    double* temp = new double[AO2USO_->max_ncol() * AO2USO_->max_nrow()];
    for (size_t N = 0; N < D_.size(); ++N) {
        // Input is already C1
        if (!input_symmetry_cast_map_[N]) {
            D_ao_[N]->copy(D_[N]);
            continue;
        }

        if (D_[N]->nirrep() != AO2USO_->nirrep()) {
            throw PSIEXCEPTION("JK::AO2USO: Dimensions of C and D do not match AO2USO!\n");
        }
        D_ao_[N]->zero();
        int symm = D_[N]->symmetry();
        for (int h = 0; h < AO2USO_->nirrep(); ++h) {
            int nao = AO2USO_->rowspi()[0];
            int nsol = AO2USO_->colspi()[h];
            int nsor = AO2USO_->colspi()[h ^ symm];
            if (!nsol || !nsor) continue;
            double** Ulp = AO2USO_->pointer(h);
            double** Urp = AO2USO_->pointer(h ^ symm);
            double** DSOp = D_[N]->pointer(h ^ symm);
            double** DAOp = D_ao_[N]->pointer();
            C_DGEMM('N', 'T', nsol, nao, nsor, 1.0, DSOp[0], nsor, Urp[0], nsor, 0.0, temp, nao);
            C_DGEMM('N', 'N', nao, nao, nsol, 1.0, Ulp[0], nsol, temp, nao, 1.0, DAOp[0], nao);
        }
    }
    delete[] temp;

    // Transform the left-index of all C matrices from SO basis to AO basis.

    // Transform C_left. Assumed totally symmetric.
    for (size_t N = 0; N < D_.size(); ++N) {
        // Input is already C1
        if (!input_symmetry_cast_map_[N]) {
            C_left_ao_[N]->copy(C_left_[N]);
            continue;
        }

        int offset = 0;
        for (int h = 0; h < AO2USO_->nirrep(); ++h) {
            int nao = AO2USO_->rowspi()[0];
            int nso = AO2USO_->colspi()[h];
            int ncol = C_left_ao_[N]->colspi()[0];
            int ncolspi = C_left_[N]->colspi()[h];
            if (nso == 0 || ncolspi == 0) continue;
            double** Up = AO2USO_->pointer(h);
            double** CAOp = C_left_ao_[N]->pointer();
            double** CSOp = C_left_[N]->pointer(h);
            C_DGEMM('N', 'N', nao, ncolspi, nso, 1.0, Up[0], nso, CSOp[0], ncolspi, 0.0, &CAOp[0][offset], ncol);
            offset += ncolspi;
        }
    }

    // Transform C_right. Not assumed totally symmetric.
    for (size_t N = 0; (N < D_.size()) && (!lr_symmetric_); ++N) {
        // Input is already C1
        if (!input_symmetry_cast_map_[N]) {
            C_right_ao_[N]->copy(C_right_[N]);
            continue;
        }

        int offset = 0;
        int symm = D_[N]->symmetry();
        for (int h = 0; h < AO2USO_->nirrep(); ++h) {
            // We MUST pack columns in the order in which they appear for totally symmetric C_left.
            // This means we transform in order of h ^ symm, not in order of h.
            int nao = AO2USO_->rowspi()[0];
            int nso = AO2USO_->colspi()[h ^ symm];
            int ncol = C_right_ao_[N]->colspi()[0];
            // Remember: colspi_[h] describes not the orbitals of block h, but the orbitals that transform as h.
            int ncolspi = C_right_[N]->colspi()[h];
            if (nso == 0 || ncolspi == 0) continue;
            double** Up = AO2USO_->pointer(h ^ symm);
            double** CAOp = C_right_ao_[N]->pointer();
            double** CSOp = C_right_[N]->pointer(h ^ symm);
            C_DGEMM('N', 'N', nao, ncolspi, nso, 1.0, Up[0], nso, CSOp[0], ncolspi, 0.0, &CAOp[0][offset], ncol);
            offset += ncolspi;
        }
    }

}
void JK::AO2USO() {
    // If already C1, J/K are J_ao/K_ao, pointers are already aliased
    if (AO2USO_->nirrep() == 1) {
        return;
    }

    // If not C1, J/K/wK are already allocated

    // Transform
    double* temp = new double[AO2USO_->max_ncol() * AO2USO_->max_nrow()];
    for (size_t N = 0; N < D_.size(); ++N) {
        // Input was desymmetrized, return as same
        if (!input_symmetry_cast_map_[N]) {
            if (do_J_) {
                J_[N]->copy(J_ao_[N]);
            }
            if (do_K_) {
                K_[N]->copy(K_ao_[N]);
            }
            if (do_wK_) {
                wK_[N]->copy(wK_ao_[N]);
            }
            continue;
        }

        int symm = D_[N]->symmetry();
        for (int h = 0; h < AO2USO_->nirrep(); ++h) {
            int nao = AO2USO_->rowspi()[0];
            int nsol = AO2USO_->colspi()[h];
            int nsor = AO2USO_->colspi()[h ^ symm];

            if (!nsol || !nsor) continue;

            double** Ulp = AO2USO_->pointer(h);
            double** Urp = AO2USO_->pointer(h ^ symm);

            if (do_J_) {
                double** JAOp = J_ao_[N]->pointer();
                double** JSOp = J_[N]->pointer(h);
                C_DGEMM('N', 'N', nao, nsor, nao, 1.0, JAOp[0], nao, Urp[0], nsor, 0.0, temp, nsor);
                C_DGEMM('T', 'N', nsol, nsor, nao, 1.0, Ulp[0], nsol, temp, nsor, 0.0, JSOp[0], nsor);
            }
            if (do_K_) {
                double** KAOp = K_ao_[N]->pointer();
                double** KSOp = K_[N]->pointer(h);
                C_DGEMM('N', 'N', nao, nsor, nao, 1.0, KAOp[0], nao, Urp[0], nsor, 0.0, temp, nsor);
                C_DGEMM('T', 'N', nsol, nsor, nao, 1.0, Ulp[0], nsol, temp, nsor, 0.0, KSOp[0], nsor);
            }
            if (do_wK_) {
                double** wKAOp = wK_ao_[N]->pointer();
                double** wKSOp = wK_[N]->pointer(h);
                C_DGEMM('N', 'N', nao, nsor, nao, 1.0, wKAOp[0], nao, Urp[0], nsor, 0.0, temp, nsor);
                C_DGEMM('T', 'N', nsol, nsor, nao, 1.0, Ulp[0], nsol, temp, nsor, 0.0, wKSOp[0], nsor);
            }
        }
    }
    delete[] temp;
}
void JK::initialize() { preiterations(); }

void JK::compute() {
    // Is this density symmetric?
    if (C_left_.size() && !C_right_.size()) {
        lr_symmetric_ = true;
        C_right_ = C_left_;
    } else {
        lr_symmetric_ = false;
    }

    // Figure out the symmetry and which codes will stay in C1 symmetry
    input_symmetry_cast_map_.clear();
    for (size_t i = 0; i < C_left_.size(); i++) {
        // Make sure they have the same symmetry
        if (C_left_[i]->nirrep() != C_right_[i]->nirrep()) {
            throw PSIEXCEPTION("JK: C_left/C_right irrep mismatch!");
        }

        // Make sure they have the same zip index
        if (C_left_[i]->colspi() != C_right_[i]->colspi()) {
            throw PSIEXCEPTION("JK: C_left/C_right MO zip index size mismatch!");
        }

        // Figure out if we need to convert or not
        if ((AO2USO_->nirrep() == 1) && (C_left_[i]->nirrep() == 1)) {
            // Everything in C1, nothing to do
            input_symmetry_cast_map_.push_back(false);
        } else if (C_left_[i]->nirrep() == AO2USO_->nirrep()) {
            // We match symmetry, does this code uses C1?
            if (C1()) {
                input_symmetry_cast_map_.push_back(true);
            } else {
                input_symmetry_cast_map_.push_back(false);
            }
        } else if ((C_left_[i]->nirrep() == 1) && C1()) {
            // Code uses C1, nothing to do
            input_symmetry_cast_map_.push_back(false);
        } else {
            // No other cases, throw
            throw PSIEXCEPTION("JK: Input orbital irrep mismatch!");
        }
    }

    // Construct the densities
    timer_on("JK: D");
    compute_D();
    timer_off("JK: D");

    if (C1()) {
        timer_on("JK: USO2AO");
        USO2AO();
        timer_off("JK: USO2AO");
    } else {
        allocate_JK();
    }

    timer_on("JK: JK");
    compute_JK();
    timer_off("JK: JK");

    if (C1()) {
        timer_on("JK: AO2USO");
        AO2USO();
        timer_off("JK: AO2USO");
    }

    if (debug_ > 6) {
        outfile->Printf("   > JK <\n\n");
        for (size_t N = 0; N < C_left_.size(); N++) {
            if (C1() && AO2USO_->nirrep() != 1) {
                C_left_ao_[N]->print("outfile");
                C_right_ao_[N]->print("outfile");
                D_ao_[N]->print("outfile");
                J_ao_[N]->print("outfile");
                K_ao_[N]->print("outfile");
            }
            C_left_[N]->print("outfile");
            C_right_[N]->print("outfile");
            D_[N]->print("outfile");
            J_[N]->print("outfile");
            K_[N]->print("outfile");
        }
    }

    if (lr_symmetric_) {
        C_right_.clear();
    }
}

void JK::set_wcombine(bool wcombine) {
    wcombine_ = wcombine;
    if (wcombine) {
        throw PSIEXCEPTION("To combine exchange terms, use MemDFJK\n");
    }
}

void JK::zero() {
    if (do_J_) {
        for(auto J : J_) J->zero();
        for(auto J : J_ao_) J->zero();
    }
    if (do_K_) {
        for(auto K : K_) K->zero();
        for(auto K : K_ao_) K->zero();
    }
    if (do_wK_) {
        for(auto wK : wK_) wK->zero();
        for(auto wK : wK_ao_) wK->zero();
    }
}

size_t JK::num_computed_shells() {
    outfile->Printf("WARNING: JK::num_computed_shells() was called, but benchmarking is disabled for the chosen JK algorithm.");
    outfile->Printf(" Returning 0 as computed shells count.\n");

    return 0;
}

void JK::finalize() { postiterations(); }
}  // namespace psi
