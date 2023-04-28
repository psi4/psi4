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
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/lib3index/dftensor.h"

#include "jk.h"

#include <sstream>
#include "psi4/libpsi4util/PsiOutStream.h"
#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

using namespace psi;

namespace psi {

DiskDFJK::DiskDFJK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options)
    : JK(primary), auxiliary_(auxiliary), options_(options) {
    common_init();
}
DiskDFJK::~DiskDFJK() {}
void DiskDFJK::common_init() {
    df_ints_num_threads_ = 1;
#ifdef _OPENMP
    df_ints_num_threads_ = Process::environment.get_n_threads();
#endif
    df_ints_io_ = "NONE";
    condition_ = 1.0E-12;
    unit_ = PSIF_DFSCF_BJ;
    if (options_["SCF_SUBTYPE"].has_changed()) set_subalgo(options_.get_str("SCF_SUBTYPE"));
    is_core_ = true;
    psio_ = PSIO::shared_object();
    // We need to make an integral object so we can figure out memory requirements from the sieve
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory =
        std::make_shared<IntegralFactory>(auxiliary_, zero, primary_, primary_);
    auto tmperi = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    n_function_pairs_ = tmperi->function_pairs().size();
}
size_t DiskDFJK::memory_estimate() {
    size_t three_memory = ((size_t)auxiliary_->nbf()) * n_function_pairs_;
    size_t two_memory = 2 * ((size_t)auxiliary_->nbf()) * auxiliary_->nbf();

    if (do_wK_) { three_memory *= 3; }

    size_t memory = three_memory + two_memory;
    memory += memory_overhead();
    memory += memory_temp();

    return memory;
}
SharedVector DiskDFJK::iaia(SharedMatrix Ci, SharedMatrix Ca) {
    // Target quantity
    Dimension dim(Ci->nirrep());
    for (int symm = 0; symm < Ci->nirrep(); symm++) {
        int rank = 0;
        for (int h = 0; h < Ci->nirrep(); h++) {
            rank += Ci->colspi()[h] * Ca->colspi()[h ^ symm];
        }
        dim[symm] = rank;
    }

    auto Iia = std::make_shared<Vector>("(ia|ia)", dim);

    // AO-basis quantities
    int nirrep = Ci->nirrep();
    int nocc = Ci->ncol();
    int nvir = Ca->ncol();
    int nso = AO2USO_->rowspi()[0];

    auto Ci_ao = std::make_shared<Matrix>("Ci AO", nso, nocc);
    auto Ca_ao = std::make_shared<Matrix>("Ca AO", nso, nvir);
    auto Iia_ao = std::make_shared<Vector>("(ia|ia) AO", nocc * (size_t)nvir);

    int offset = 0;
    for (int h = 0; h < nirrep; h++) {
        int ni = Ci->colspi()[h];
        int nm = Ci->rowspi()[h];
        if (!ni || !nm) continue;
        double** Cip = Ci->pointer(h);
        double** Cp = Ci_ao->pointer();
        double** Up = AO2USO_->pointer(h);
        C_DGEMM('N', 'N', nso, ni, nm, 1.0, Up[0], nm, Cip[0], ni, 0.0, &Cp[0][offset], nocc);
        offset += ni;
    }

    offset = 0;
    for (int h = 0; h < nirrep; h++) {
        int ni = Ca->colspi()[h];
        int nm = Ca->rowspi()[h];
        if (!ni || !nm) continue;
        double** Cip = Ca->pointer(h);
        double** Cp = Ca_ao->pointer();
        double** Up = AO2USO_->pointer(h);
        C_DGEMM('N', 'N', nso, ni, nm, 1.0, Up[0], nm, Cip[0], ni, 0.0, &Cp[0][offset], nvir);
        offset += ni;
    }

    // Memory size
    int naux = auxiliary_->nbf();
    int maxrows = max_rows();

    const std::vector<std::pair<int, int> >& function_pairs = eri_.front()->function_pairs();
    const std::vector<long int>& function_pairs_to_dense = eri_.front()->function_pairs_to_dense();
    size_t num_nm = function_pairs.size();

// Temps
#ifdef _OPENMP
    int temp_nthread = Process::environment.get_n_threads();
    omp_set_num_threads(omp_nthread_);
    C_temp_.resize(omp_nthread_);
    Q_temp_.resize(omp_nthread_);
#pragma omp parallel
    {
        C_temp_[omp_get_thread_num()] = std::make_shared<Matrix>("Ctemp", nocc, nso);
        Q_temp_[omp_get_thread_num()] = std::make_shared<Matrix>("Qtemp", maxrows, nso);
    }
    omp_set_num_threads(temp_nthread);
#else
    for (int thread = 0; thread < omp_nthread_; thread++) {
        C_temp_.push_back(std::make_shared<Matrix>("Ctemp", nocc, nso));
        Q_temp_.push_back(std::make_shared<Matrix>("Qtemp", maxrows, nso));
    }
#endif

    E_left_ = std::make_shared<Matrix>("E_left", nso, maxrows * nocc);
    E_right_ = std::make_shared<Matrix>("E_right", nvir, maxrows * nocc);

    // Disk overhead
    psio_address addr = PSIO_ZERO;
    if (!is_core()) {
        Qmn_ = std::make_shared<Matrix>("(Q|mn) Block", maxrows, num_nm);
        psio_->open(unit_, PSIO_OPEN_OLD);
    }

    // Blocks of Q
    double** Qmnp;
    double** Clp = Ci_ao->pointer();
    double** Crp = Ca_ao->pointer();
    double** Elp = E_left_->pointer();
    double** Erp = E_right_->pointer();
    double* Iiap = Iia_ao->pointer();
    for (int Q = 0; Q < naux; Q += maxrows) {
        // Read block of (Q|mn) in
        int rows = (naux - Q <= maxrows ? naux - Q : maxrows);
        if (is_core()) {
            Qmnp = &Qmn_->pointer()[Q];
        } else {
            Qmnp = Qmn_->pointer();
            psio_->read(unit_, "(Q|mn) Integrals", (char*)(Qmn_->pointer()[0]), sizeof(double) * naux * num_nm, addr,
                        &addr);
        }

// (mi|Q)
#pragma omp parallel for schedule(dynamic)
        for (int m = 0; m < nso; m++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            double** Ctp = C_temp_[thread]->pointer();
            double** QSp = Q_temp_[thread]->pointer();

            const std::vector<int>& pairs = eri_.front()->significant_partners_per_function()[m];
            int mrows = pairs.size();

            for (int i = 0; i < mrows; i++) {
                int n = pairs[i];
                long int ij = function_pairs_to_dense[(m >= n ? (m * (m + 1L) >> 1) + n : (n * (n + 1L) >> 1) + m)];
                C_DCOPY(rows, &Qmnp[0][ij], num_nm, &QSp[0][i], nso);
                C_DCOPY(nocc, Clp[n], 1, &Ctp[0][i], nso);
            }

            C_DGEMM('N', 'T', nocc, rows, mrows, 1.0, Ctp[0], nso, QSp[0], nso, 0.0, &Elp[0][m * (size_t)nocc * rows],
                    rows);
        }

        // (ai|Q)
        C_DGEMM('T', 'N', nvir, nocc * (size_t)rows, nso, 1.0, Crp[0], nvir, Elp[0], nocc * (size_t)rows, 0.0, Erp[0],
                nocc * (size_t)rows);

        // (ia|Q)(Q|ia)
        for (int i = 0; i < nocc; i++) {
            for (int a = 0; a < nvir; a++) {
                double* Ep = &Erp[0][a * static_cast<size_t>(nocc) * rows + i * rows];
                Iiap[i * nvir + a] += C_DDOT(rows, Ep, 1, Ep, 1);
            }
        }
    }

    // Free disk overhead
    if (!is_core()) {
        Qmn_.reset();
        psio_->close(unit_, 1);
    }

    // Free Temps
    E_left_.reset();
    E_right_.reset();
    C_temp_.clear();
    Q_temp_.clear();

    // SO-basis (ia|ia)
    Dimension i_offsets(Ci->nirrep());
    Dimension a_offsets(Ci->nirrep());
    for (int h = 1; h < Ci->nirrep(); h++) {
        i_offsets[h] = i_offsets[h - 1] + Ci->colspi()[h - 1];
        a_offsets[h] = a_offsets[h - 1] + Ca->colspi()[h - 1];
    }

    for (int symm = 0; symm < Ci->nirrep(); symm++) {
        offset = 0;
        for (int h = 0; h < Ci->nirrep(); h++) {
            int ni = Ci->colspi()[h];
            int na = Ca->colspi()[h ^ symm];
            int ioff = i_offsets[h];
            int aoff = a_offsets[h ^ symm];
            for (int i = 0; i < ni; i++) {
                for (int a = 0; a < na; a++) {
                    Iia->set(symm, i * na + a + offset, Iiap[(ioff + i) * nvir + (aoff + a)]);
                }
            }
            offset += ni * na;
        }
    }

    return Iia;
}
void DiskDFJK::print_header() const {
    if (print_) {
        outfile->Printf("  ==> DiskDFJK: Density-Fitted J/K Matrices <==\n\n");

        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("    Omega:             %11.3E\n", omega_);
        outfile->Printf("    OpenMP threads:    %11d\n", omp_nthread_);
        outfile->Printf("    Integrals threads: %11d\n", df_ints_num_threads_);
        outfile->Printf("    Memory [MiB]:      %11ld\n", (memory_ * 8L) / (1024L * 1024L));
        outfile->Printf("    Algorithm:         %11s\n", (is_core_ ? "Core" : "Disk"));
        outfile->Printf("    Integral Cache:    %11s\n", df_ints_io_.c_str());
        outfile->Printf("    Schwarz Cutoff:    %11.0E\n", cutoff_);
        outfile->Printf("    Fitting Condition: %11.0E\n\n", condition_);

        outfile->Printf("   => Auxiliary Basis Set <=\n\n");
        auxiliary_->print_by_level("outfile", print_);
    }
}
bool DiskDFJK::is_core() {
    auto do_core = is_core_;

    // determine do_core either automatically...
    if (subalgo_ == "AUTO") {
        do_core = memory_estimate() < memory_;

    // .. or forcibly disable do_core if user specifies ...
    } else if (subalgo_ == "OUT_OF_CORE") {
        if (print_ > 0) {
            outfile->Printf("  SCF_SUBTYPE = OUT_OF_CORE selected. Out-of-core DISK_DF algorithm will be used.\n\n");
        }

        do_core = false;

   // .. or force do_core if user specifies
    } else if (subalgo_ == "INCORE") {
        if (memory_estimate() > memory_) {
            throw PSIEXCEPTION("SCF_SUBTYPE=INCORE was specified, but there is not enough memory to do in-core! Increase the amount of memory allocated to Psi4 or allow for out-of-core to be used.\n");
        } else {
            if (print_ > 0) {
                outfile->Printf("  SCF_SUBTYPE=INCORE selected. In-core DISK_DF algorithm will be used.\n\n");
            }
            do_core = true;
        }
    } else {
        throw PSIEXCEPTION("Invalid SCF_SUBTYPE option! The choices for SCF_SUBTYPE are AUTO, INCORE, and OUT_OF_CORE.");
    }

    return do_core;
}

size_t DiskDFJK::memory_temp() const {
    size_t mem = 0L;

    // J Overhead (Jtri, Dtri, d)
    mem += 2L * n_function_pairs_ + auxiliary_->nbf();
    // K Overhead (C_temp, Q_temp)
    mem += omp_nthread_ * (size_t)primary_->nbf() * (auxiliary_->nbf() + max_nocc());

    return mem;
}
int DiskDFJK::max_rows() const {
    // Start with all memory
    size_t mem = memory_;
    // Subtract J/K/wK/C/D overhead
    mem -= memory_overhead();
    // Subtract threading temp overhead
    mem -= memory_temp();

    // How much will each row cost?
    size_t row_cost = 0L;
    // Copies of E tensor
    row_cost += (lr_symmetric_ ? 1L : 2L) * max_nocc() * primary_->nbf();
    // Slices of Qmn tensor, including AIO buffer (NOTE: AIO not implemented yet)
    row_cost += (is_core_ ? 1L : 1L) * n_function_pairs_;

    size_t max_rows = mem / row_cost;

    if (max_rows > (size_t)auxiliary_->nbf()) max_rows = (size_t)auxiliary_->nbf();
    if (max_rows < 1L) max_rows = 1L;

    return (int)max_rows;
}
int DiskDFJK::max_nocc() const {
    int max_nocc = 0;
    for (size_t N = 0; N < C_left_ao_.size(); N++) {
        max_nocc = (C_left_ao_[N]->colspi()[0] > max_nocc ? C_left_ao_[N]->colspi()[0] : max_nocc);
    }
    return max_nocc;
}
void DiskDFJK::initialize_temps() {
    J_temp_ = std::make_shared<Vector>("Jtemp", n_function_pairs_);
    D_temp_ = std::make_shared<Vector>("Dtemp", n_function_pairs_);
    d_temp_ = std::make_shared<Vector>("dtemp", max_rows_);

#ifdef _OPENMP
    int temp_nthread = Process::environment.get_n_threads();
    omp_set_num_threads(omp_nthread_);
    C_temp_.resize(omp_nthread_);
    Q_temp_.resize(omp_nthread_);
#pragma omp parallel
    {
        C_temp_[omp_get_thread_num()] = std::make_shared<Matrix>("Ctemp", max_nocc_, primary_->nbf());
        Q_temp_[omp_get_thread_num()] = std::make_shared<Matrix>("Qtemp", max_rows_, primary_->nbf());
    }
    omp_set_num_threads(temp_nthread);
#else
    for (int thread = 0; thread < omp_nthread_; thread++) {
        C_temp_.push_back(std::make_shared<Matrix>("Ctemp", max_nocc_, primary_->nbf()));
        Q_temp_.push_back(std::make_shared<Matrix>("Qtemp", max_rows_, primary_->nbf()));
    }
#endif

    E_left_ = std::make_shared<Matrix>("E_left", primary_->nbf(), max_rows_ * max_nocc_);
    if (lr_symmetric_)
        E_right_ = E_left_;
    else
        E_right_ = std::make_shared<Matrix>("E_right", primary_->nbf(), max_rows_ * max_nocc_);
}
void DiskDFJK::initialize_w_temps() {
    int max_rows_w = max_rows_ / 2;
    max_rows_w = (max_rows_w < 1 ? 1 : max_rows_w);

#ifdef _OPENMP
    int temp_nthread = Process::environment.get_n_threads();
    omp_set_num_threads(omp_nthread_);
    C_temp_.resize(omp_nthread_);
    Q_temp_.resize(omp_nthread_);
#pragma omp parallel
    {
        C_temp_[omp_get_thread_num()] = std::make_shared<Matrix>("Ctemp", max_nocc_, primary_->nbf());
        Q_temp_[omp_get_thread_num()] = std::make_shared<Matrix>("Qtemp", max_rows_w, primary_->nbf());
    }
    omp_set_num_threads(temp_nthread);
#else
    for (int thread = 0; thread < omp_nthread_; thread++) {
        C_temp_.push_back(std::make_shared<Matrix>("Ctemp", max_nocc_, primary_->nbf()));
        Q_temp_.push_back(std::make_shared<Matrix>("Qtemp", max_rows_w, primary_->nbf()));
    }
#endif

    E_left_ = std::make_shared<Matrix>("E_left", primary_->nbf(), max_rows_w * max_nocc_);
    E_right_ = std::make_shared<Matrix>("E_right", primary_->nbf(), max_rows_w * max_nocc_);
}
void DiskDFJK::free_temps() {
    J_temp_.reset();
    D_temp_.reset();
    d_temp_.reset();
    E_left_.reset();
    E_right_.reset();
    C_temp_.clear();
    Q_temp_.clear();
}
void DiskDFJK::free_w_temps() {
    E_left_.reset();
    E_right_.reset();
    C_temp_.clear();
    Q_temp_.clear();
}
void DiskDFJK::preiterations() {
    // Setup integral objects
    eri_.clear();
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory =
        std::make_shared<IntegralFactory>(auxiliary_, zero, primary_, primary_);
    eri_.emplace_back(rifactory->eri());
    for (int Q = 1; Q < df_ints_num_threads_; Q++) {
        eri_.emplace_back(eri_.front()->clone());
    }
    n_function_pairs_ = eri_.front()->function_pairs().size();
    // Setup erf integrals, if needed
    if (do_wK_) {
        erf_eri_.clear();
        std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
        std::shared_ptr<IntegralFactory> rifactory =
            std::make_shared<IntegralFactory>(auxiliary_, zero, primary_, primary_);
        erf_eri_.emplace_back(rifactory->erf_eri(omega_));
        for (int Q = 1; Q < df_ints_num_threads_; Q++) {
            erf_eri_.emplace_back(erf_eri_.front()->clone());
        }
    }

    // Core or disk?
    is_core_ = is_core();

    if (is_core_)
        initialize_JK_core();
    else
        initialize_JK_disk();

    if (do_wK_) {
        if (is_core_)
            initialize_wK_core();
        else
            initialize_wK_disk();
    }
}

void DiskDFJK::compute_JK() {

    // zero out J, K, and wK matrices
    zero();

    max_nocc_ = max_nocc();
    max_rows_ = max_rows();

    if (do_J_ || do_K_) {
        initialize_temps();
        if (is_core_)
            manage_JK_core();
        else
            manage_JK_disk();
        free_temps();
    }

    if (do_wK_) {
        initialize_w_temps();
        if (is_core_)
            manage_wK_core();
        else
            manage_wK_disk();
        free_w_temps();
        // Bring the wK matrices back to Hermitian
        if (lr_symmetric_) {
            for (size_t N = 0; N < wK_ao_.size(); N++) {
                wK_ao_[N]->hermitivitize();
            }
        }
    }
}
void DiskDFJK::postiterations() {
    Qmn_.reset();
    Qlmn_.reset();
    Qrmn_.reset();
}
void DiskDFJK::initialize_JK_core() {
    size_t three_memory = ((size_t)auxiliary_->nbf()) * n_function_pairs_;
    size_t two_memory = ((size_t)auxiliary_->nbf()) * auxiliary_->nbf();

    int nthread = 1;
#ifdef _OPENMP
    nthread = df_ints_num_threads_;
#endif

    Qmn_ = std::make_shared<Matrix>("Qmn (Fitted Integrals)", auxiliary_->nbf(), n_function_pairs_);
    double** Qmnp = Qmn_->pointer();

    // Try to load
    if (df_ints_io_ == "LOAD") {
        psio_->open(unit_, PSIO_OPEN_OLD);
        psio_->read_entry(unit_, "(Q|mn) Integrals", (char*)Qmnp[0], sizeof(double) * n_function_pairs_ * auxiliary_->nbf());
        psio_->close(unit_, 1);
        return;
    }

    const int primary_maxam = primary_->max_am();
    const int aux_maxam = auxiliary_->max_am();
    const int primary_nshell = primary_->nshell();
    const int aux_nshell = auxiliary_->nshell();

    const std::vector<long int>& schwarz_shell_pairs = eri_.front()->shell_pairs_to_dense();
    const std::vector<long int>& schwarz_fun_pairs = eri_.front()->function_pairs_to_dense();

    // shell pair blocks
    std::vector<ShellPairBlock> p_blocks = eri_[0]->get_blocks12();
    std::vector<ShellPairBlock> mn_blocks = eri_[0]->get_blocks34();

    timer_on("JK: (A|mn)");

// We are calculting integrals of (p|mn). Note that mn is on the ket side.
// Integrals may be vectorized over the ket shells, so we want more of the
// work to be on that side.
//
// Variable description:
// mn = mu nu (in old terminology)
// p_block_idx = index of what block of p we are doing
// mn_block_idx = index of what block of mn we are doing
// p_block = The block of p we are doing (a vector of std::pair<int, int>)
// mn_block = The block of mn we are doing
// mn_pair = a single pair of indicies, corresponding to the indicies of the
//           m and n shells in the primary basis set
// p_pair = a single pair of indicies, corresponding to the indicies of the
//          p shell in the auxiliary basis set
// num_m, num_n, num_p = number of basis functions in shell m, n, and p
// m_start, n_start, p_start = starting index of the basis functions of this shell
//                             relative to the beginning of their corresponding basis
//                             sets.
// im, in, ip = indices for looping over all the basis functions of the shells
// im_idx, in_idx, ip_idx = indices of a particular basis function in a shell
//                          relative to the start of the basis set
//                          (ie, m_start + im)
// sfp = screen function pair. The actual output index, taking into account
//                             whether or not functions have been screened.

#pragma omp parallel for schedule(dynamic) num_threads(nthread)
    for (size_t mn_block_idx = 0; mn_block_idx < mn_blocks.size(); mn_block_idx++) {
#ifdef _OPENMP
        const int rank = omp_get_thread_num();
#else
        const int rank = 0;
#endif
        const auto &buffers = eri_[rank]->buffers();
        const auto& mn_block = mn_blocks[mn_block_idx];

        // loop over all the blocks of P
        for (int p_block_idx = 0; p_block_idx < p_blocks.size(); ++p_block_idx) {
            // compute the
            eri_[rank]->compute_shell_blocks(p_block_idx, mn_block_idx);
            const auto* buffer = buffers[0];

            const auto& p_block = p_blocks[p_block_idx];

            for (const auto& mn_pair : mn_block) {
                const int m = mn_pair.first;
                const int n = mn_pair.second;
                const int num_m = primary_->shell(m).nfunction();
                const int num_n = primary_->shell(n).nfunction();
                const int m_start = primary_->shell(m).function_index();
                const int n_start = primary_->shell(n).function_index();
                const int num_mn = num_m * num_n;

                for (const auto& p_pair : p_block) {
                    // remember that this vector will only contain one shell pair
                    const int p = p_pair.first;
                    const int num_p = auxiliary_->shell(p).nfunction();
                    const int p_start = auxiliary_->shell(p).function_index();

                    for (int im = 0; im < num_m; ++im) {
                        const int im_idx = m_start + im;

                        for (int in = 0; in < num_n; ++in) {
                            const int in_idx = n_start + in;
                            const int imn_idx = im * num_n + in;
                            const int sfp_idx = im_idx > in_idx ? (im_idx * (im_idx + 1)) / 2 + in_idx :
                                               (in_idx * (in_idx + 1))/2 + im_idx;

                            int sfp;

                            // note the assignment in the following conditional
                            if ((sfp = schwarz_fun_pairs[sfp_idx]) > -1) {
                                for (int ip = 0; ip < num_p; ++ip) {
                                    const int ip_idx = p_start + ip;
                                    Qmnp[ip_idx][sfp] = buffer[ip * num_mn + imn_idx];
                                }
                            }

                        }
                    }

                    buffer += num_mn * num_p;
                }
            }
        }
    }
    timer_off("JK: (A|mn)");

    timer_on("JK: (A|Q)^-1/2");

    auto Jinv = std::make_shared<FittingMetric>(auxiliary_, true);
    Jinv->form_eig_inverse(condition_);
    double** Jinvp = Jinv->get_metric()->pointer();

    timer_off("JK: (A|Q)^-1/2");

    size_t max_cols = (memory_ - three_memory - two_memory) / auxiliary_->nbf();
    if (max_cols < 1) max_cols = 1;
    if (max_cols > n_function_pairs_) max_cols = n_function_pairs_;
    auto temp = std::make_shared<Matrix>("Qmn buffer", auxiliary_->nbf(), max_cols);
    double** tempp = temp->pointer();

    size_t nblocks = n_function_pairs_ / max_cols;
    if ((size_t)nblocks * max_cols != n_function_pairs_) nblocks++;

    size_t ncol = 0;
    size_t col = 0;

    timer_on("JK: (Q|mn)");

    for (size_t block = 0; block < nblocks; block++) {
        ncol = max_cols;
        if (col + ncol > n_function_pairs_) ncol = n_function_pairs_ - col;

        C_DGEMM('N', 'N', auxiliary_->nbf(), ncol, auxiliary_->nbf(), 1.0, Jinvp[0], auxiliary_->nbf(), &Qmnp[0][col],
                n_function_pairs_, 0.0, tempp[0], max_cols);

        for (int Q = 0; Q < auxiliary_->nbf(); Q++) {
            C_DCOPY(ncol, tempp[Q], 1, &Qmnp[Q][col], 1);
        }

        col += ncol;
    }

    timer_off("JK: (Q|mn)");
    // Qmn_->print();

    if (df_ints_io_ == "SAVE") {
        psio_->open(unit_, PSIO_OPEN_NEW);
        psio_->write_entry(unit_, "(Q|mn) Integrals", (char*)Qmnp[0], sizeof(double) * n_function_pairs_ * auxiliary_->nbf());
        psio_->close(unit_, 1);
    }
}
void DiskDFJK::initialize_JK_disk() {
    // Try to load
    if (df_ints_io_ == "LOAD") {
        return;
    }

    int nshell = primary_->nshell();
    int naux = auxiliary_->nbf();

    // ==> Schwarz Indexing <== //
    const std::vector<std::pair<int, int> >& schwarz_shell_pairs = eri_.front()->shell_pairs();
    const std::vector<std::pair<int, int> >& schwarz_fun_pairs = eri_.front()->function_pairs();
    int nshellpairs = schwarz_shell_pairs.size();
    const std::vector<long int>& schwarz_shell_pairs_r = eri_.front()->shell_pairs_to_dense();
    const std::vector<long int>& schwarz_fun_pairs_r = eri_.front()->function_pairs_to_dense();

    // ==> Memory Sizing <== //
    size_t two_memory = ((size_t)auxiliary_->nbf()) * auxiliary_->nbf();
    size_t buffer_memory =
        (memory_ > 2 * two_memory) ? memory_ - 2 * two_memory : 0;  // Two is for buffer space in fitting

    // outfile->Printf( "Buffer memory = %ld words\n", buffer_memory);

    // outfile->Printf("Schwarz Shell Pairs:\n");
    // for (int MN = 0; MN < nshellpairs; MN++) {
    //    outfile->Printf("  %3d: (%3d,%3d)\n", MN, schwarz_shell_pairs[2*MN], schwarz_shell_pairs[2*MN + 1]);
    //}

    // outfile->Printf("Schwarz Function Pairs:\n");
    // for (int MN = 0; MN < ntri; MN++) {
    //    outfile->Printf("  %3d: (%3d,%3d)\n", MN, schwarz_fun_pairs[2*MN], schwarz_fun_pairs[2*MN + 1]);
    //}

    // outfile->Printf("Schwarz Reverse Shell Pairs:\n");
    // for (int MN = 0; MN < primary_->nshell() * (primary_->nshell() + 1) / 2; MN++) {
    //    outfile->Printf("  %3d: %4ld\n", MN, schwarz_shell_pairs_r[MN]);
    //}

    // outfile->Printf("Schwarz Reverse Function Pairs:\n");
    // for (int MN = 0; MN < primary_->nbf() * (primary_->nbf() + 1) / 2; MN++) {
    //    outfile->Printf("  %3d: %4ld\n", MN, schwarz_fun_pairs_r[MN]);
    //}

    // Find out exactly how much memory per MN shell
    auto MN_mem = std::make_shared<IntVector>("Memory per MN pair", nshell * (nshell + 1) / 2);
    int* MN_memp = MN_mem->pointer();

    for (int mn = 0; mn < n_function_pairs_; mn++) {
        int m = schwarz_fun_pairs[mn].first;
        int n = schwarz_fun_pairs[mn].second;

        int M = primary_->function_to_shell(m);
        int N = primary_->function_to_shell(n);
        size_t addr = M > N ? M * (M + 1) / 2 + N : N * (N + 1) / 2 + M;
        MN_memp[addr] += naux;
    }

    // MN_mem->print(outfile);

    // Figure out exactly how much memory per M row
    size_t* M_memp = new size_t[nshell];
    memset(static_cast<void*>(M_memp), '\0', nshell * sizeof(size_t));

    for (int M = 0; M < nshell; M++) {
        for (int N = 0; N <= M; N++) {
            M_memp[M] += MN_memp[M * (M + 1) / 2 + N];
        }
    }

    // outfile->Printf("  # Memory per M row #\n\n");
    // for (int M = 0; M < nshell; M++)
    //    outfile->Printf("   %3d: %10ld\n", M+1,M_memp[M]);
    // outfile->Printf("\n");

    // Find and check the minimum required memory for this problem
    size_t min_mem = naux * (size_t)n_function_pairs_;
    for (int M = 0; M < nshell; M++) {
        if (min_mem > M_memp[M]) min_mem = M_memp[M];
    }

    if (min_mem > buffer_memory) {
        std::stringstream message;
        message << "SCF::DF: Disk based algorithm requires 2 (A|B) fitting metrics and an (A|mn) chunk on core."
                << std::endl;
        message << "         This is 2Q^2 + QNP doubles, where Q is the auxiliary basis size, N is the" << std::endl;
        message << "         primary basis size, and P is the maximum number of functions in a primary shell."
                << std::endl;
        message << "         For this problem, that is " << ((8L * (min_mem + 2 * two_memory)))
                << " bytes before taxes,";
        message << ((80L * (min_mem + 2 * two_memory) / 7L)) << " bytes after taxes. " << std::endl;

        throw PSIEXCEPTION(message.str());
    }

    // ==> Reduced indexing by M <== //

    // Figure out the MN start index per M row
    auto MN_start = std::make_shared<IntVector>("MUNU start per M row", nshell);
    int* MN_startp = MN_start->pointer();

    MN_startp[0] = schwarz_shell_pairs_r[0];
    int M_index = 1;
    for (int MN = 0; MN < nshellpairs; MN++) {
        if (std::max(schwarz_shell_pairs[MN].first, schwarz_shell_pairs[MN].second) == M_index) {
            MN_startp[M_index] = MN;
            M_index++;
        }
    }

    // Figure out the mn start index per M row
    auto mn_start = std::make_shared<IntVector>("munu start per M row", nshell);
    int* mn_startp = mn_start->pointer();

    mn_startp[0] = schwarz_fun_pairs[0].first;
    int m_index = 1;
    for (int mn = 0; mn < n_function_pairs_; mn++) {
        if (primary_->function_to_shell(schwarz_fun_pairs[mn].first) == m_index) {
            mn_startp[m_index] = mn;
            m_index++;
        }
    }

    // Figure out the MN columns per M row
    auto MN_col = std::make_shared<IntVector>("MUNU cols per M row", nshell);
    int* MN_colp = MN_col->pointer();

    for (int M = 1; M < nshell; M++) {
        MN_colp[M - 1] = MN_startp[M] - MN_startp[M - 1];
    }
    MN_colp[nshell - 1] = nshellpairs - MN_startp[nshell - 1];

    // Figure out the mn columns per M row
    auto mn_col = std::make_shared<IntVector>("munu cols per M row", nshell);
    int* mn_colp = mn_col->pointer();

    for (int M = 1; M < nshell; M++) {
        mn_colp[M - 1] = mn_startp[M] - mn_startp[M - 1];
    }
    mn_colp[nshell - 1] = n_function_pairs_ - mn_startp[nshell - 1];

    // MN_start->print(outfile);
    // MN_col->print(outfile);
    // mn_start->print(outfile);
    // mn_col->print(outfile);

    // ==> Block indexing <== //
    // Sizing by block
    std::vector<int> MN_start_b;
    std::vector<int> MN_col_b;
    std::vector<int> mn_start_b;
    std::vector<int> mn_col_b;

    // Determine MN and mn block starts
    // also MN and mn block cols
    int nblock = 1;
    size_t current_mem = 0L;
    MN_start_b.push_back(0);
    mn_start_b.push_back(0);
    MN_col_b.push_back(0);
    mn_col_b.push_back(0);
    for (int M = 0; M < nshell; M++) {
        if (current_mem + M_memp[M] > buffer_memory) {
            MN_start_b.push_back(MN_startp[M]);
            mn_start_b.push_back(mn_startp[M]);
            MN_col_b.push_back(0);
            mn_col_b.push_back(0);
            nblock++;
            current_mem = 0L;
        }
        MN_col_b[nblock - 1] += MN_colp[M];
        mn_col_b[nblock - 1] += mn_colp[M];
        current_mem += M_memp[M];
    }

    // outfile->Printf("Block, MN start, MN cols, mn start, mn cols\n");
    // for (int block = 0; block < nblock; block++) {
    //    outfile->Printf("  %3d: %12d %12d %12d %12d\n", block, MN_start_b[block], MN_col_b[block], mn_start_b[block],
    //    mn_col_b[block]);
    //}
    //

    // Full sizing not required any longer
    MN_mem.reset();
    MN_start.reset();
    MN_col.reset();
    mn_start.reset();
    mn_col.reset();
    delete[] M_memp;

    // ==> Buffer allocation <== //
    int max_cols = 0;
    for (int block = 0; block < nblock; block++) {
        if (max_cols < mn_col_b[block]) max_cols = mn_col_b[block];
    }

    // Primary buffer
    Qmn_ = std::make_shared<Matrix>("(Q|mn) (Disk Chunk)", naux, max_cols);
    // Fitting buffer
    auto Amn = std::make_shared<Matrix>("(Q|mn) (Buffer)", naux, naux);
    double** Qmnp = Qmn_->pointer();
    double** Amnp = Amn->pointer();

    // ==> Prestripe/Jinv <== //

    timer_on("JK: (A|Q)^-1");

    psio_->open(unit_, PSIO_OPEN_NEW);
    auto aio = std::make_shared<AIOHandler>(psio_);

    // Dispatch the prestripe
    aio->zero_disk(unit_, "(Q|mn) Integrals", naux, n_function_pairs_);

    // Form the J symmetric inverse
    auto Jinv = std::make_shared<FittingMetric>(auxiliary_, true);
    Jinv->form_eig_inverse(condition_);
    double** Jinvp = Jinv->get_metric()->pointer();

    // Synch up
    aio->synchronize();

    timer_off("JK: (A|Q)^-1");

    // ==> Thread setup <== //
    int nthread = 1;
#ifdef _OPENMP
    nthread = df_ints_num_threads_;
#endif

    // ==> Main loop <== //
    for (int block = 0; block < nblock; block++) {
        int MN_start_val = MN_start_b[block];
        int mn_start_val = mn_start_b[block];
        int MN_col_val = MN_col_b[block];
        int mn_col_val = mn_col_b[block];

        // ==> (A|mn) integrals <== //

        timer_on("JK: (A|mn)");

#pragma omp parallel for schedule(guided) num_threads(nthread)
        for (int MUNU = MN_start_val; MUNU < MN_start_val + MN_col_val; MUNU++) {
            int rank = 0;
#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif
            const auto &buffers = eri_[rank]->buffers();
            int MU = schwarz_shell_pairs[MUNU + 0].first;
            int NU = schwarz_shell_pairs[MUNU + 0].second;
            int nummu = primary_->shell(MU).nfunction();
            int numnu = primary_->shell(NU).nfunction();
            int mu = primary_->shell(MU).function_index();
            int nu = primary_->shell(NU).function_index();
            for (int P = 0; P < auxiliary_->nshell(); P++) {
                int nump = auxiliary_->shell(P).nfunction();
                int p = auxiliary_->shell(P).function_index();
                eri_[rank]->compute_shell(P, 0, MU, NU);
                const auto *buffer = buffers[0];
                for (int dm = 0; dm < nummu; dm++) {
                    int omu = mu + dm;
                    for (int dn = 0; dn < numnu; dn++) {
                        int onu = nu + dn;
                        size_t addr = omu > onu ? omu * (omu + 1) / 2 + onu :
                                                  onu * (onu + 1) / 2 + omu;
                        if (schwarz_fun_pairs_r[addr] >= 0) {
                            int delta = schwarz_fun_pairs_r[addr] - mn_start_val;
                            for (int dp = 0; dp < nump; dp++) {
                                int op = p + dp;
                                Qmnp[op][delta] = buffer[dp * nummu * numnu + dm * numnu + dn];
                            }
                        }
                    }
                }
            }
        }

        timer_off("JK: (A|mn)");

        // ==> (Q|mn) fitting <== //

        timer_on("JK: (Q|mn)");

        for (int mn = 0; mn < mn_col_val; mn += naux) {
            int cols = naux;
            if (mn + naux >= mn_col_val) cols = mn_col_val - mn;

            for (int Q = 0; Q < naux; Q++) C_DCOPY(cols, &Qmnp[Q][mn], 1, Amnp[Q], 1);

            C_DGEMM('N', 'N', naux, cols, naux, 1.0, Jinvp[0], naux, Amnp[0], naux, 0.0, &Qmnp[0][mn], max_cols);
        }

        timer_off("JK: (Q|mn)");

        // ==> Disk striping <== //

        timer_on("JK: (Q|mn) Write");

        psio_address addr;
        for (int Q = 0; Q < naux; Q++) {
            addr = psio_get_address(PSIO_ZERO, (Q * (size_t)n_function_pairs_ + mn_start_val) * sizeof(double));
            psio_->write(unit_, "(Q|mn) Integrals", (char*)Qmnp[Q], mn_col_val * sizeof(double), addr, &addr);
        }

        timer_off("JK: (Q|mn) Write");
    }

    // ==> Close out <== //
    Qmn_.reset();

    psio_->close(unit_, 1);
}
void DiskDFJK::initialize_wK_core() {
    int naux = auxiliary_->nbf();

    int nthread = 1;
#ifdef _OPENMP
    nthread = df_ints_num_threads_;
#endif
    int rank = 0;

    // Check that the right integrals are using the correct omega
    if (df_ints_io_ == "LOAD") {
        psio_->open(unit_, PSIO_OPEN_OLD);
        double check_omega;
        psio_->read_entry(unit_, "Omega", (char*)&check_omega, sizeof(double));
        if (check_omega != omega_) {
            rebuild_wK_disk();
        }
        psio_->close(unit_, 1);
    }

    Qlmn_ = std::make_shared<Matrix>("Qlmn (Fitted Integrals)", auxiliary_->nbf(), n_function_pairs_);
    double** Qmnp = Qlmn_->pointer();

    Qrmn_ = std::make_shared<Matrix>("Qrmn (Fitted Integrals)", auxiliary_->nbf(), n_function_pairs_);
    double** Qmn2p = Qrmn_->pointer();

    // Try to load
    if (df_ints_io_ == "LOAD") {
        psio_->open(unit_, PSIO_OPEN_OLD);
        psio_->read_entry(unit_, "Left (Q|w|mn) Integrals", (char*)Qmnp[0], sizeof(double) * n_function_pairs_ * auxiliary_->nbf());
        psio_->read_entry(unit_, "Right (Q|w|mn) Integrals", (char*)Qmn2p[0],
                          sizeof(double) * n_function_pairs_ * auxiliary_->nbf());
        psio_->close(unit_, 1);
        return;
    }

    // => Left Integrals <= //

    const std::vector<long int>& schwarz_shell_pairs = eri_.front()->shell_pairs_to_dense();
    const std::vector<long int>& schwarz_fun_pairs = eri_.front()->function_pairs_to_dense();

    int numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu;
    // The integrals (A|mn)

    timer_on("JK: (A|mn)^L");

#pragma omp parallel for private(numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, \
                                 rank) schedule(dynamic) num_threads(nthread)
    for (MU = 0; MU < primary_->nshell(); ++MU) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        const auto& buffers = eri_[rank]->buffers();
        nummu = primary_->shell(MU).nfunction();
        for (NU = 0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU).nfunction();
            size_t addr = MU > NU ? MU * (MU + 1) / 2 + NU : NU * (NU + 1) / 2 + MU;
            if (schwarz_shell_pairs[addr] > -1) {
                for (Pshell = 0; Pshell < auxiliary_->nshell(); ++Pshell) {
                    numP = auxiliary_->shell(Pshell).nfunction();
                    eri_[rank]->compute_shell(Pshell, 0, MU, NU);
                    const auto *buffer = buffers[0];
                    for (mu = 0; mu < nummu; ++mu) {
                        omu = primary_->shell(MU).function_index() + mu;
                        for (nu = 0; nu < numnu; ++nu) {
                            onu = primary_->shell(NU).function_index() + nu;
                            size_t addr = omu > onu ? omu * (omu + 1) / 2 + onu :
                                                      onu * (onu + 1) / 2 + omu;
                            if (schwarz_fun_pairs[addr] > -1) {
                                for (P = 0; P < numP; ++P) {
                                    PHI = auxiliary_->shell(Pshell).function_index() + P;
                                    Qmn2p[PHI][schwarz_fun_pairs[addr]] =
                                        buffer[P * nummu * numnu + mu * numnu + nu];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    timer_off("JK: (A|mn)^L");

    // => Fitting <= //

    timer_on("JK: (A|Q)^-1");

    // Fitting metric
    auto Jinv = std::make_shared<FittingMetric>(auxiliary_, true);
    Jinv->form_full_eig_inverse(condition_);
    double** Jinvp = Jinv->get_metric()->pointer();

    timer_off("JK: (A|Q)^-1");

    timer_on("JK: (Q|mn)^L");

    // Fitting in one GEMM (being a clever bastard)
    C_DGEMM('N', 'N', naux, n_function_pairs_, naux, 1.0, Jinvp[0], naux, Qmn2p[0], n_function_pairs_, 0.0, Qmnp[0], n_function_pairs_);

    timer_off("JK: (Q|mn)^L");

    // The integrals (A|w|mn)

    timer_on("JK: (A|mn)^R");

#pragma omp parallel for private(numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, \
                                 rank) schedule(dynamic) num_threads(nthread)
    for (MU = 0; MU < primary_->nshell(); ++MU) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        const auto& buffers = erf_eri_[rank]->buffers();
        nummu = primary_->shell(MU).nfunction();
        for (NU = 0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU).nfunction();
            size_t addr = MU > NU ? MU * (MU + 1) / 2 + NU : NU * (NU + 1) / 2 + MU;
            if (schwarz_shell_pairs[addr] > -1) {
                for (Pshell = 0; Pshell < auxiliary_->nshell(); ++Pshell) {
                    numP = auxiliary_->shell(Pshell).nfunction();
                    erf_eri_[rank]->compute_shell(Pshell, 0, MU, NU);
                    const auto *buffer = buffers[0];
                    for (mu = 0; mu < nummu; ++mu) {
                        omu = primary_->shell(MU).function_index() + mu;
                        for (nu = 0; nu < numnu; ++nu) {
                            onu = primary_->shell(NU).function_index() + nu;
                            size_t addr = omu > onu ? omu * (omu + 1) / 2 + onu :
                                                      onu * (onu + 1) / 2 + omu;
                            if (schwarz_fun_pairs[addr] > -1) {
                                for (P = 0; P < numP; ++P) {
                                    PHI = auxiliary_->shell(Pshell).function_index() + P;
                                    Qmn2p[PHI][schwarz_fun_pairs[addr]] =
                                        buffer[P * nummu * numnu + mu * numnu + nu];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    timer_off("JK: (A|mn)^R");

    // Try to save
    if (df_ints_io_ == "SAVE") {
        psio_->open(unit_, PSIO_OPEN_OLD);
        psio_->write_entry(unit_, "Left (Q|w|mn) Integrals", (char*)Qmnp[0], sizeof(double) * n_function_pairs_ * auxiliary_->nbf());
        psio_->write_entry(unit_, "Right (Q|w|mn) Integrals", (char*)Qmn2p[0],
                           sizeof(double) * n_function_pairs_ * auxiliary_->nbf());
        psio_->write_entry(unit_, "Omega", (char*)&omega_, sizeof(double));
        psio_->close(unit_, 1);
    }
}
void DiskDFJK::initialize_wK_disk() {
    // Try to load
    if (df_ints_io_ == "LOAD") {
        psio_->open(unit_, PSIO_OPEN_OLD);
        double check_omega;
        psio_->read_entry(unit_, "Omega", (char*)&check_omega, sizeof(double));
        if (check_omega != omega_) {
            rebuild_wK_disk();
        }
        psio_->close(unit_, 1);
    }

    size_t nshell = primary_->nshell();
    size_t naux = auxiliary_->nbf();

    // ==> Schwarz Indexing <== //
    const std::vector<std::pair<int, int> >& schwarz_shell_pairs = eri_.front()->shell_pairs();
    const std::vector<std::pair<int, int> >& schwarz_fun_pairs = eri_.front()->function_pairs();
    int nshellpairs = schwarz_shell_pairs.size();
    const std::vector<long int>& schwarz_shell_pairs_r = eri_.front()->shell_pairs_to_dense();
    const std::vector<long int>& schwarz_fun_pairs_r = eri_.front()->function_pairs_to_dense();

    // ==> Memory Sizing <== //
    size_t two_memory = ((size_t)auxiliary_->nbf()) * auxiliary_->nbf();
    size_t buffer_memory =
        (memory_ > 2 * two_memory) ? memory_ - 2 * two_memory : 0;  // Two is for buffer space in fitting

    // outfile->Printf( "Buffer memory = %ld words\n", buffer_memory);

    // outfile->Printf("Schwarz Shell Pairs:\n");
    // for (int MN = 0; MN < nshellpairs; MN++) {
    //    outfile->Printf("  %3d: (%3d,%3d)\n", MN, schwarz_shell_pairs[2*MN], schwarz_shell_pairs[2*MN + 1]);
    //}

    // outfile->Printf("Schwarz Function Pairs:\n");
    // for (int MN = 0; MN < ntri; MN++) {
    //    outfile->Printf("  %3d: (%3d,%3d)\n", MN, schwarz_fun_pairs[2*MN], schwarz_fun_pairs[2*MN + 1]);
    //}

    // outfile->Printf("Schwarz Reverse Shell Pairs:\n");
    // for (int MN = 0; MN < primary_->nshell() * (primary_->nshell() + 1) / 2; MN++) {
    //    outfile->Printf("  %3d: %4ld\n", MN, schwarz_shell_pairs_r[MN]);
    //}

    // outfile->Printf("Schwarz Reverse Function Pairs:\n");
    // for (int MN = 0; MN < primary_->nbf() * (primary_->nbf() + 1) / 2; MN++) {
    //    outfile->Printf("  %3d: %4ld\n", MN, schwarz_fun_pairs_r[MN]);
    //}

    // Find out exactly how much memory per MN shell
    auto MN_mem = std::make_shared<IntVector>("Memory per MN pair", nshell * (nshell + 1) / 2);
    int* MN_memp = MN_mem->pointer();

    for (int mn = 0; mn < n_function_pairs_; mn++) {
        int m = schwarz_fun_pairs[mn].first;
        int n = schwarz_fun_pairs[mn].second;

        int M = primary_->function_to_shell(m);
        int N = primary_->function_to_shell(n);
        size_t addr = M > N ? M * (M + 1) / 2 + N : N * (N + 1) / 2 + M;
        MN_memp[addr] += naux;
    }

    // MN_mem->print(outfile);

    // Figure out exactly how much memory per M row
    size_t* M_memp = new size_t[nshell];
    memset(static_cast<void*>(M_memp), '\0', nshell * sizeof(size_t));

    for (size_t M = 0; M < nshell; M++) {
        for (size_t N = 0; N <= M; N++) {
            M_memp[M] += MN_memp[M * (M + 1) / 2 + N];
        }
    }

    // outfile->Printf("  # Memory per M row #\n\n");
    // for (int M = 0; M < nshell; M++)
    //    outfile->Printf("   %3d: %10ld\n", M+1,M_memp[M]);
    // outfile->Printf("\n");

    // Find and check the minimum required memory for this problem
    size_t min_mem = naux * (size_t)n_function_pairs_;
    for (size_t M = 0; M < nshell; M++) {
        if (min_mem > M_memp[M]) min_mem = M_memp[M];
    }

    if (min_mem > buffer_memory) {
        std::stringstream message;
        message << "SCF::DF: Disk based algorithm requires 2 (A|B) fitting metrics and an (A|mn) chunk on core."
                << std::endl;
        message << "         This is 2Q^2 + QNP doubles, where Q is the auxiliary basis size, N is the" << std::endl;
        message << "         primary basis size, and P is the maximum number of functions in a primary shell."
                << std::endl;
        message << "         For this problem, that is " << ((8L * (min_mem + 2 * two_memory)))
                << " bytes before taxes,";
        message << ((80L * (min_mem + 2 * two_memory) / 7L)) << " bytes after taxes. " << std::endl;

        throw PSIEXCEPTION(message.str());
    }

    // ==> Reduced indexing by M <== //

    // Figure out the MN start index per M row
    auto MN_start = std::make_shared<IntVector>("MUNU start per M row", nshell);
    int* MN_startp = MN_start->pointer();

    MN_startp[0] = schwarz_shell_pairs_r[0];
    int M_index = 1;
    for (int MN = 0; MN < nshellpairs; MN++) {
        if (std::max(schwarz_shell_pairs[MN].first, schwarz_shell_pairs[MN].second) == M_index) {
            MN_startp[M_index] = MN;
            M_index++;
        }
    }

    // Figure out the mn start index per M row
    auto mn_start = std::make_shared<IntVector>("munu start per M row", nshell);
    int* mn_startp = mn_start->pointer();

    mn_startp[0] = schwarz_fun_pairs[0].first;
    int m_index = 1;
    for (int mn = 0; mn < n_function_pairs_; mn++) {
        int fn = std::max(schwarz_fun_pairs[mn].first, schwarz_fun_pairs[mn].second);
        if (primary_->function_to_shell(fn) == m_index) {
            mn_startp[m_index] = mn;
            m_index++;
        }
    }

    // Figure out the MN columns per M row
    auto MN_col = std::make_shared<IntVector>("MUNU cols per M row", nshell);
    int* MN_colp = MN_col->pointer();

    for (size_t M = 1; M < nshell; M++) {
        MN_colp[M - 1] = MN_startp[M] - MN_startp[M - 1];
    }
    MN_colp[nshell - 1] = nshellpairs - MN_startp[nshell - 1];

    // Figure out the mn columns per M row
    auto mn_col = std::make_shared<IntVector>("munu cols per M row", nshell);
    int* mn_colp = mn_col->pointer();

    for (size_t M = 1; M < nshell; M++) {
        mn_colp[M - 1] = mn_startp[M] - mn_startp[M - 1];
    }
    mn_colp[nshell - 1] = n_function_pairs_ - mn_startp[nshell - 1];

    // MN_start->print(outfile);
    // MN_col->print(outfile);
    // mn_start->print(outfile);
    // mn_col->print(outfile);

    // ==> Block indexing <== //
    // Sizing by block
    std::vector<int> MN_start_b;
    std::vector<int> MN_col_b;
    std::vector<int> mn_start_b;
    std::vector<int> mn_col_b;

    // Determine MN and mn block starts
    // also MN and mn block cols
    int nblock = 1;
    size_t current_mem = 0L;
    MN_start_b.push_back(0);
    mn_start_b.push_back(0);
    MN_col_b.push_back(0);
    mn_col_b.push_back(0);
    for (size_t M = 0; M < nshell; M++) {
        if (current_mem + M_memp[M] > buffer_memory) {
            MN_start_b.push_back(MN_startp[M]);
            mn_start_b.push_back(mn_startp[M]);
            MN_col_b.push_back(0);
            mn_col_b.push_back(0);
            nblock++;
            current_mem = 0L;
        }
        MN_col_b[nblock - 1] += MN_colp[M];
        mn_col_b[nblock - 1] += mn_colp[M];
        current_mem += M_memp[M];
    }

    // outfile->Printf("Block, MN start, MN cols, mn start, mn cols\n");
    // for (int block = 0; block < nblock; block++) {
    //    outfile->Printf("  %3d: %12d %12d %12d %12d\n", block, MN_start_b[block], MN_col_b[block], mn_start_b[block],
    //    mn_col_b[block]);
    //}
    //

    // Full sizing not required any longer
    MN_mem.reset();
    MN_start.reset();
    MN_col.reset();
    mn_start.reset();
    mn_col.reset();
    delete[] M_memp;

    // ==> Buffer allocation <== //
    int max_cols = 0;
    for (int block = 0; block < nblock; block++) {
        if (max_cols < mn_col_b[block]) max_cols = mn_col_b[block];
    }

    // Primary buffer
    Qmn_ = std::make_shared<Matrix>("(Q|mn) (Disk Chunk)", naux, max_cols);
    // Fitting buffer
    auto Amn = std::make_shared<Matrix>("(Q|mn) (Buffer)", naux, naux);
    double** Qmnp = Qmn_->pointer();
    double** Amnp = Amn->pointer();

    // ==> Prestripe/Jinv <== //
    psio_->open(unit_, PSIO_OPEN_OLD);
    auto aio = std::make_shared<AIOHandler>(psio_);

    // Dispatch the prestripe
    aio->zero_disk(unit_, "Left (Q|w|mn) Integrals", naux, n_function_pairs_);

    // Form the J full inverse
    auto Jinv = std::make_shared<FittingMetric>(auxiliary_, true);
    Jinv->form_full_eig_inverse(condition_);
    double** Jinvp = Jinv->get_metric()->pointer();

    // Synch up
    aio->synchronize();

    // ==> Thread setup <== //
    int nthread = 1;
#ifdef _OPENMP
    nthread = df_ints_num_threads_;
#endif

    // ==> Main loop <== //
    for (int block = 0; block < nblock; block++) {
        int MN_start_val = MN_start_b[block];
        int mn_start_val = mn_start_b[block];
        int MN_col_val = MN_col_b[block];
        int mn_col_val = mn_col_b[block];

        // ==> (A|mn) integrals <== //

        timer_on("JK: (A|mn)^L");

#pragma omp parallel for schedule(guided) num_threads(nthread)
        for (int MUNU = MN_start_val; MUNU < MN_start_val + MN_col_val; MUNU++) {
            int rank = 0;
#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif
            const auto& buffers = eri_[rank]->buffers();
            int MU = schwarz_shell_pairs[MUNU + 0].first;
            int NU = schwarz_shell_pairs[MUNU + 0].second;
            int nummu = primary_->shell(MU).nfunction();
            int numnu = primary_->shell(NU).nfunction();
            int mu = primary_->shell(MU).function_index();
            int nu = primary_->shell(NU).function_index();
            for (int P = 0; P < auxiliary_->nshell(); P++) {
                int nump = auxiliary_->shell(P).nfunction();
                int p = auxiliary_->shell(P).function_index();
                eri_[rank]->compute_shell(P, 0, MU, NU);
                const auto *buffer = buffers[0];
                for (int dm = 0; dm < nummu; dm++) {
                    int omu = mu + dm;
                    for (int dn = 0; dn < numnu; dn++) {
                        int onu = nu + dn;
                        size_t addr = omu > onu ? omu * (omu + 1) / 2 + onu :
                                                  onu * (onu + 1) / 2 + omu;
                        if (schwarz_fun_pairs_r[addr] >= 0) {
                            int delta = schwarz_fun_pairs_r[addr] - mn_start_val;
                            for (int dp = 0; dp < nump; dp++) {
                                int op = p + dp;
                                Qmnp[op][delta] = buffer[dp * nummu * numnu + dm * numnu + dn];
                            }
                        }
                    }
                }
            }
        }

        timer_off("JK: (A|mn)^L");

        // ==> (Q|mn) fitting <== //

        timer_on("JK: (Q|mn)^L");

        for (int mn = 0; mn < mn_col_val; mn += naux) {
            int cols = naux;
            if (mn + naux >= (size_t)mn_col_val) cols = mn_col_val - mn;

            for (size_t Q = 0; Q < naux; Q++) C_DCOPY(cols, &Qmnp[Q][mn], 1, Amnp[Q], 1);

            C_DGEMM('N', 'N', naux, cols, naux, 1.0, Jinvp[0], naux, Amnp[0], naux, 0.0, &Qmnp[0][mn], max_cols);
        }

        timer_off("JK: (Q|mn)^L");

        // ==> Disk striping <== //

        timer_on("JK: (Q|mn)^L Write");

        psio_address addr;
        for (size_t Q = 0; Q < naux; Q++) {
            addr = psio_get_address(PSIO_ZERO, (Q * (size_t)n_function_pairs_ + mn_start_val) * sizeof(double));
            psio_->write(unit_, "Left (Q|w|mn) Integrals", (char*)Qmnp[Q], mn_col_val * sizeof(double), addr, &addr);
        }

        timer_off("JK: (Q|mn)^L Write");
    }

    Qmn_.reset();

    // => Right Integrals <= //

    size_t maxP = auxiliary_->max_function_per_shell();
    size_t max_rows = memory_ / n_function_pairs_;
    max_rows = (max_rows > naux ? naux : max_rows);
    max_rows = (max_rows < maxP ? maxP : max_rows);

    // Block extents
    std::vector<int> block_Q_starts;
    size_t counter = 0;
    block_Q_starts.push_back(0);
    for (int Q = 0; Q < auxiliary_->nshell(); Q++) {
        int nQ = auxiliary_->shell(Q).nfunction();
        if (counter + nQ > max_rows) {
            counter = 0;
            block_Q_starts.push_back(Q);
        }
        counter += nQ;
    }
    block_Q_starts.push_back(auxiliary_->nshell());

    auto Amn2 = std::make_shared<Matrix>("(A|mn) Block", max_rows, n_function_pairs_);
    double** Amn2p = Amn2->pointer();
    psio_address next_AIA = PSIO_ZERO;

    const std::vector<std::pair<int, int> >& shell_pairs = eri_.front()->shell_pairs();
    const size_t npairs = shell_pairs.size();

    // Loop over blocks of Qshell
    for (size_t block = 0; block < block_Q_starts.size() - 1; block++) {
        // Block sizing/offsets
        int Qstart = block_Q_starts[block];
        int Qstop = block_Q_starts[block + 1];
        int qoff = auxiliary_->shell(Qstart).function_index();
        int nrows = (Qstop == auxiliary_->nshell()
                         ? auxiliary_->nbf() - auxiliary_->shell(Qstart).function_index()
                         : auxiliary_->shell(Qstop).function_index() - auxiliary_->shell(Qstart).function_index());

        // Compute TEI tensor block (A|mn)

        timer_on("JK: (Q|mn)^R");

#pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (size_t QMN = 0L; QMN < (Qstop - Qstart) * (size_t)npairs; QMN++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif
            const auto &buffers = erf_eri_[thread]->buffers();
            int Q = QMN / npairs + Qstart;
            int MN = QMN % npairs;

            std::pair<int, int> pair = shell_pairs[MN];
            int M = pair.first;
            int N = pair.second;

            int nq = auxiliary_->shell(Q).nfunction();
            int nm = primary_->shell(M).nfunction();
            int nn = primary_->shell(N).nfunction();

            int sq = auxiliary_->shell(Q).function_index();
            int sm = primary_->shell(M).function_index();
            int sn = primary_->shell(N).function_index();

            erf_eri_[thread]->compute_shell(Q, 0, M, N);
            const auto *buffer = buffers[0];

            for (int om = 0; om < nm; om++) {
                for (int on = 0; on < nn; on++) {
                    long int m = sm + om;
                    long int n = sn + on;
                    size_t addr = m > n ? m * (m + 1) / 2 + n : n * (n + 1) / 2 + m;
                    if (schwarz_fun_pairs_r[addr] >= 0) {
                        long int delta = schwarz_fun_pairs_r[addr];
                        for (int oq = 0; oq < nq; oq++) {
                            Amn2p[sq + oq - qoff][delta] = buffer[oq * nm * nn + om * nn + on];
                        }
                    }
                }
            }
        }

        timer_off("JK: (Q|mn)^R");

        // Dump block to disk
        timer_on("JK: (Q|mn)^R Write");

        psio_->write(unit_, "Right (Q|w|mn) Integrals", (char*)Amn2p[0], sizeof(double) * nrows * n_function_pairs_, next_AIA,
                     &next_AIA);

        timer_off("JK: (Q|mn)^R Write");
    }
    Amn2.reset();

    psio_->write_entry(unit_, "Omega", (char*)&omega_, sizeof(double));
    psio_->close(unit_, 1);
}
void DiskDFJK::rebuild_wK_disk() {
    // Already open
    outfile->Printf("    Rebuilding (Q|w|mn) Integrals (new omega)\n\n");

    size_t naux = auxiliary_->nbf();

    // ==> Schwarz Indexing <== //
    const std::vector<std::pair<int, int> >& schwarz_fun_pairs = eri_.front()->function_pairs();
    const std::vector<long int>& schwarz_fun_pairs_r = eri_.front()->function_pairs_to_dense();

    // ==> Thread setup <== //
    int nthread = 1;
#ifdef _OPENMP
    nthread = df_ints_num_threads_;
#endif

    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory =
        std::make_shared<IntegralFactory>(auxiliary_, zero, primary_, primary_);

    size_t maxP = auxiliary_->max_function_per_shell();
    size_t max_rows = memory_ / n_function_pairs_;
    max_rows = (max_rows > naux ? naux : max_rows);
    max_rows = (max_rows < maxP ? maxP : max_rows);

    // Block extents
    std::vector<int> block_Q_starts;
    int counter = 0;
    block_Q_starts.push_back(0);
    for (int Q = 0; Q < auxiliary_->nshell(); Q++) {
        size_t nQ = auxiliary_->shell(Q).nfunction();
        if (counter + nQ > max_rows) {
            counter = 0;
            block_Q_starts.push_back(Q);
        }
        counter += nQ;
    }
    block_Q_starts.push_back(auxiliary_->nshell());

    auto Amn2 = std::make_shared<Matrix>("(A|mn) Block", max_rows, n_function_pairs_);
    double** Amn2p = Amn2->pointer();
    psio_address next_AIA = PSIO_ZERO;

    const std::vector<std::pair<int, int> >& shell_pairs = eri_.front()->shell_pairs();
    const size_t npairs = shell_pairs.size();

    // Loop over blocks of Qshell
    for (size_t block = 0; block < block_Q_starts.size() - 1; block++) {
        // Block sizing/offsets
        int Qstart = block_Q_starts[block];
        int Qstop = block_Q_starts[block + 1];
        int qoff = auxiliary_->shell(Qstart).function_index();
        int nrows = (Qstop == auxiliary_->nshell()
                         ? auxiliary_->nbf() - auxiliary_->shell(Qstart).function_index()
                         : auxiliary_->shell(Qstop).function_index() - auxiliary_->shell(Qstart).function_index());

        // Compute TEI tensor block (A|mn)

        timer_on("JK: (Q|mn)^R");

#pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (size_t QMN = 0L; QMN < (Qstop - Qstart) * (size_t)npairs; QMN++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif
            const auto &buffers = erf_eri_[thread]->buffers();
            int Q = QMN / npairs + Qstart;
            int MN = QMN % npairs;

            std::pair<int, int> pair = shell_pairs[MN];
            int M = pair.first;
            int N = pair.second;

            int nq = auxiliary_->shell(Q).nfunction();
            int nm = primary_->shell(M).nfunction();
            int nn = primary_->shell(N).nfunction();

            int sq = auxiliary_->shell(Q).function_index();
            int sm = primary_->shell(M).function_index();
            int sn = primary_->shell(N).function_index();

            erf_eri_[thread]->compute_shell(Q, 0, M, N);
            const auto *buffer = buffers[0];

            for (int om = 0; om < nm; om++) {
                for (int on = 0; on < nn; on++) {
                    long int m = sm + om;
                    long int n = sn + on;
                    size_t addr = m > n ? m * (m + 1) / 2 + n : n * (n + 1) / 2 + m;
                    if (schwarz_fun_pairs_r[addr] >= 0) {
                        long int delta = schwarz_fun_pairs_r[addr];
                        for (int oq = 0; oq < nq; oq++) {
                            Amn2p[sq + oq - qoff][delta] = buffer[oq * nm * nn + om * nn + on];
                        }
                    }
                }
            }
        }

        timer_off("JK: (Q|mn)^R");

        // Dump block to disk
        timer_on("JK: (Q|mn)^R Write");

        psio_->write(unit_, "Right (Q|w|mn) Integrals", (char*)Amn2p[0], sizeof(double) * nrows * n_function_pairs_, next_AIA,
                     &next_AIA);

        timer_off("JK: (Q|mn)^R Write");
    }
    Amn2.reset();

    psio_->write_entry(unit_, "Omega", (char*)&omega_, sizeof(double));
    // No need to close
}
void DiskDFJK::manage_JK_core() {
    for (int Q = 0; Q < auxiliary_->nbf(); Q += max_rows_) {
        int naux = (auxiliary_->nbf() - Q <= max_rows_ ? auxiliary_->nbf() - Q : max_rows_);
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
void DiskDFJK::manage_JK_disk() {
    Qmn_ = std::make_shared<Matrix>("(Q|mn) Block", max_rows_, n_function_pairs_);
    psio_->open(unit_, PSIO_OPEN_OLD);
    for (int Q = 0; Q < auxiliary_->nbf(); Q += max_rows_) {
        int naux = (auxiliary_->nbf() - Q <= max_rows_ ? auxiliary_->nbf() - Q : max_rows_);
        psio_address addr = psio_get_address(PSIO_ZERO, (Q * (size_t)n_function_pairs_) * sizeof(double));

        timer_on("JK: (Q|mn) Read");
        psio_->read(unit_, "(Q|mn) Integrals", (char*)(Qmn_->pointer()[0]), sizeof(double) * naux * n_function_pairs_, addr, &addr);
        timer_off("JK: (Q|mn) Read");

        if (do_J_) {
            timer_on("JK: J");
            block_J(&Qmn_->pointer()[0], naux);
            timer_off("JK: J");
        }
        if (do_K_) {
            timer_on("JK: K");
            block_K(&Qmn_->pointer()[0], naux);
            timer_off("JK: K");
        }
    }
    psio_->close(unit_, 1);
    Qmn_.reset();
}
void DiskDFJK::manage_wK_core() {
    int max_rows_w = max_rows_ / 2;
    max_rows_w = (max_rows_w < 1 ? 1 : max_rows_w);
    for (int Q = 0; Q < auxiliary_->nbf(); Q += max_rows_w) {
        int naux = (auxiliary_->nbf() - Q <= max_rows_w ? auxiliary_->nbf() - Q : max_rows_w);

        timer_on("JK: wK");
        block_wK(&Qlmn_->pointer()[Q], &Qrmn_->pointer()[Q], naux);
        timer_off("JK: wK");
    }
}
void DiskDFJK::manage_wK_disk() {
    int max_rows_w = max_rows_ / 2;
    max_rows_w = (max_rows_w < 1 ? 1 : max_rows_w);
    Qlmn_ = std::make_shared<Matrix>("(Q|mn) Block", max_rows_w, n_function_pairs_);
    Qrmn_ = std::make_shared<Matrix>("(Q|mn) Block", max_rows_w, n_function_pairs_);
    psio_->open(unit_, PSIO_OPEN_OLD);
    for (int Q = 0; Q < auxiliary_->nbf(); Q += max_rows_w) {
        int naux = (auxiliary_->nbf() - Q <= max_rows_w ? auxiliary_->nbf() - Q : max_rows_w);
        psio_address addr = psio_get_address(PSIO_ZERO, (Q * (size_t)n_function_pairs_) * sizeof(double));

        timer_on("JK: (Q|mn)^L Read");
        psio_->read(unit_, "Left (Q|w|mn) Integrals", (char*)(Qlmn_->pointer()[0]), sizeof(double) * naux * n_function_pairs_, addr,
                    &addr);
        timer_off("JK: (Q|mn)^L Read");

        addr = psio_get_address(PSIO_ZERO, (Q * (size_t)n_function_pairs_) * sizeof(double));

        timer_on("JK: (Q|mn)^R Read");
        psio_->read(unit_, "Right (Q|w|mn) Integrals", (char*)(Qrmn_->pointer()[0]), sizeof(double) * naux * n_function_pairs_, addr,
                    &addr);
        timer_off("JK: (Q|mn)^R Read");

        timer_on("JK: wK");
        block_wK(&Qlmn_->pointer()[0], &Qrmn_->pointer()[0], naux);
        timer_off("JK: wK");
    }
    psio_->close(unit_, 1);
    Qlmn_.reset();
    Qrmn_.reset();
}
void DiskDFJK::block_J(double** Qmnp, int naux) {
    const std::vector<std::pair<int, int> >& function_pairs = eri_.front()->function_pairs();
    size_t num_nm = function_pairs.size();

    for (size_t N = 0; N < J_ao_.size(); N++) {
        double** Dp = D_ao_[N]->pointer();
        double** Jp = J_ao_[N]->pointer();
        double* J2p = J_temp_->pointer();
        double* D2p = D_temp_->pointer();
        double* dp = d_temp_->pointer();
        for (size_t mn = 0; mn < num_nm; ++mn) {
            int m = function_pairs[mn].first;
            int n = function_pairs[mn].second;
            D2p[mn] = (m == n ? Dp[m][n] : Dp[m][n] + Dp[n][m]);
        }

        timer_on("JK: J1");
        C_DGEMV('N', naux, num_nm, 1.0, Qmnp[0], num_nm, D2p, 1, 0.0, dp, 1);
        timer_off("JK: J1");

        timer_on("JK: J2");
        C_DGEMV('T', naux, num_nm, 1.0, Qmnp[0], num_nm, dp, 1, 0.0, J2p, 1);
        timer_off("JK: J2");
        for (size_t mn = 0; mn < num_nm; ++mn) {
            int m = function_pairs[mn].first;
            int n = function_pairs[mn].second;
            Jp[m][n] += J2p[mn];
            Jp[n][m] += (m == n ? 0.0 : J2p[mn]);
        }
    }
}
void DiskDFJK::block_K(double** Qmnp, int naux) {
    const std::vector<std::pair<int, int> >& function_pairs = eri_.front()->function_pairs();
    const std::vector<long int>& function_pairs_to_dense = eri_.front()->function_pairs_to_dense();
    size_t num_nm = function_pairs.size();

    for (size_t N = 0; N < K_ao_.size(); N++) {
        int nbf = C_left_ao_[N]->rowspi()[0];
        int nocc = C_left_ao_[N]->colspi()[0];

        if (!nocc) continue;

        double** Clp = C_left_ao_[N]->pointer();
        double** Crp = C_right_ao_[N]->pointer();
        double** Elp = E_left_->pointer();
        double** Erp = E_right_->pointer();
        double** Kp = K_ao_[N]->pointer();

        if (N == 0 || C_left_[N].get() != C_left_[N - 1].get()) {
            timer_on("JK: K1");

#pragma omp parallel for schedule(dynamic)
            for (int m = 0; m < nbf; m++) {
                int thread = 0;
#ifdef _OPENMP
                thread = omp_get_thread_num();
#endif

                double** Ctp = C_temp_[thread]->pointer();
                double** QSp = Q_temp_[thread]->pointer();

                const std::vector<int>& pairs = eri_.front()->significant_partners_per_function()[m];
                int rows = pairs.size();

                for (int i = 0; i < rows; i++) {
                    int n = pairs[i];
                    long int ij = function_pairs_to_dense[(m >= n ? (m * (m + 1L) >> 1) + n : (n * (n + 1L) >> 1) + m)];
                    C_DCOPY(naux, &Qmnp[0][ij], num_nm, &QSp[0][i], nbf);
                    C_DCOPY(nocc, Clp[n], 1, &Ctp[0][i], nbf);
                }

                C_DGEMM('N', 'T', nocc, naux, rows, 1.0, Ctp[0], nbf, QSp[0], nbf, 0.0,
                        &Elp[0][m * (size_t)nocc * naux], naux);
            }

            timer_off("JK: K1");
        }

        if (!lr_symmetric_ && (N == 0 || C_right_[N].get() != C_right_[N - 1].get())) {
            if (C_right_[N].get() == C_left_[N].get()) {
                ::memcpy((void*)Erp[0], (void*)Elp[0], sizeof(double) * naux * nocc * nbf);
            } else {
                timer_on("JK: K1");

#pragma omp parallel for schedule(dynamic)
                for (int m = 0; m < nbf; m++) {
                    int thread = 0;
#ifdef _OPENMP
                    thread = omp_get_thread_num();
#endif

                    double** Ctp = C_temp_[thread]->pointer();
                    double** QSp = Q_temp_[thread]->pointer();

                    const std::vector<int>& pairs = eri_.front()->significant_partners_per_function()[m];
                    int rows = pairs.size();

                    for (int i = 0; i < rows; i++) {
                        int n = pairs[i];
                        long int ij =
                            function_pairs_to_dense[(m >= n ? (m * (m + 1L) >> 1) + n : (n * (n + 1L) >> 1) + m)];
                        C_DCOPY(naux, &Qmnp[0][ij], num_nm, &QSp[0][i], nbf);
                        C_DCOPY(nocc, Crp[n], 1, &Ctp[0][i], nbf);
                    }

                    C_DGEMM('N', 'T', nocc, naux, rows, 1.0, Ctp[0], nbf, QSp[0], nbf, 0.0,
                            &Erp[0][m * (size_t)nocc * naux], naux);
                }

                timer_off("JK: K1");
            }
        }

        timer_on("JK: K2");
        C_DGEMM('N', 'T', nbf, nbf, naux * nocc, 1.0, Elp[0], naux * nocc, Erp[0], naux * nocc, 1.0, Kp[0], nbf);
        timer_off("JK: K2");
    }
}
void DiskDFJK::block_wK(double** Qlmnp, double** Qrmnp, int naux) {
    const std::vector<std::pair<int, int> >& function_pairs = eri_.front()->function_pairs();
    const std::vector<long int>& function_pairs_to_dense = eri_.front()->function_pairs_to_dense();
    size_t num_nm = function_pairs.size();

    for (size_t N = 0; N < wK_ao_.size(); N++) {
        int nbf = C_left_ao_[N]->rowspi()[0];
        int nocc = C_left_ao_[N]->colspi()[0];

        if (!nocc) continue;

        double** Clp = C_left_ao_[N]->pointer();
        double** Crp = C_right_ao_[N]->pointer();
        double** Elp = E_left_->pointer();
        double** Erp = E_right_->pointer();
        double** wKp = wK_ao_[N]->pointer();

        if (N == 0 || C_left_[N].get() != C_left_[N - 1].get()) {
            timer_on("JK: wK1");

#pragma omp parallel for schedule(dynamic)
            for (int m = 0; m < nbf; m++) {
                int thread = 0;
#ifdef _OPENMP
                thread = omp_get_thread_num();
#endif

                double** Ctp = C_temp_[thread]->pointer();
                double** QSp = Q_temp_[thread]->pointer();

                const std::vector<int>& pairs = eri_.front()->significant_partners_per_function()[m];
                int rows = pairs.size();

                for (int i = 0; i < rows; i++) {
                    int n = pairs[i];
                    long int ij = function_pairs_to_dense[(m >= n ? (m * (m + 1L) >> 1) + n : (n * (n + 1L) >> 1) + m)];
                    C_DCOPY(naux, &Qlmnp[0][ij], num_nm, &QSp[0][i], nbf);
                    C_DCOPY(nocc, Clp[n], 1, &Ctp[0][i], nbf);
                }

                C_DGEMM('N', 'T', nocc, naux, rows, 1.0, Ctp[0], nbf, QSp[0], nbf, 0.0,
                        &Elp[0][m * (size_t)nocc * naux], naux);
            }

            timer_off("JK: wK1");
        }

        timer_on("JK: wK1");

#pragma omp parallel for schedule(dynamic)
        for (int m = 0; m < nbf; m++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            double** Ctp = C_temp_[thread]->pointer();
            double** QSp = Q_temp_[thread]->pointer();

            const std::vector<int>& pairs = eri_.front()->significant_partners_per_function()[m];
            int rows = pairs.size();

            for (int i = 0; i < rows; i++) {
                int n = pairs[i];
                long int ij = function_pairs_to_dense[(m >= n ? (m * (m + 1L) >> 1) + n : (n * (n + 1L) >> 1) + m)];
                C_DCOPY(naux, &Qrmnp[0][ij], num_nm, &QSp[0][i], nbf);
                C_DCOPY(nocc, Crp[n], 1, &Ctp[0][i], nbf);
            }

            C_DGEMM('N', 'T', nocc, naux, rows, 1.0, Ctp[0], nbf, QSp[0], nbf, 0.0, &Erp[0][m * (size_t)nocc * naux],
                    naux);
        }

        timer_off("JK: wK1");

        timer_on("JK: wK2");
        C_DGEMM('N', 'T', nbf, nbf, naux * nocc, 1.0, Elp[0], naux * nocc, Erp[0], naux * nocc, 1.0, wKp[0], nbf);
        timer_off("JK: wK2");
    }
}
}  // namespace psi
