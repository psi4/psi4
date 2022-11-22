/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "dfhelper.h"

#include <algorithm>
#include <cstdlib>
#ifdef _MSC_VER
#include <process.h>
#define SYSTEM_GETPID ::_getpid
#else
#include <unistd.h>
#define SYSTEM_GETPID ::getpid
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include <memory>
#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libfock/jk.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"

#include "dftensor.h"

namespace psi {

DFHelper::DFHelper(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> aux) {
    primary_ = primary;
    aux_ = aux;

    nbf_ = primary_->nbf();
    naux_ = aux_->nbf();
    prepare_blocking();
}

DFHelper::~DFHelper() { clear_all(); }

void DFHelper::prepare_blocking() {
    Qshells_ = aux_->nshell();
    pshells_ = primary_->nshell();

    Qshell_aggs_.resize(Qshells_ + 1);
    pshell_aggs_.resize(pshells_ + 1);

    // Aux shell blocking
    Qshell_max_ = aux_->max_function_per_shell();
    Qshell_aggs_[0] = 0;
    for (size_t i = 0; i < Qshells_; i++) {
        Qshell_aggs_[i + 1] = Qshell_aggs_[i] + aux_->shell(i).nfunction();
    }

    // AO shell blocking
    pshell_aggs_[0] = 0;
    for (size_t i = 0; i < pshells_; i++) {
        pshell_aggs_[i + 1] = pshell_aggs_[i] + primary_->shell(i).nfunction();
    }
}

void DFHelper::AO_filename_maker(size_t i) {
    auto name = start_filename("dfh.AO" + std::to_string(i));
    AO_names_.push_back(name);
}

std::string DFHelper::start_filename(std::string start) {
    std::string name = PSIOManager::shared_object()->get_default_path();
    name += start + "." + std::to_string(SYSTEM_GETPID());
    name += "." + primary_->molecule()->name() + ".";
    name += std::to_string(rand()) + "." + ".dat";
    return name;
}

void DFHelper::filename_maker(std::string name, size_t Q, size_t p, size_t q, size_t op) {
    auto pfilename = start_filename("dfh.p" + name);
    auto filename = start_filename("dfh" + name);

    std::tuple<std::string, std::string> files(pfilename, filename);
    files_[name] = files;

    bool is_transf = transf_.count(name);

    // direct_iaQ is special, because it has two different sizes
    if (direct_iaQ_ && is_transf) {
        sizes_[pfilename] = std::make_tuple(Q, p, q);
        sizes_[filename] = std::make_tuple(p, q, Q);

    } else {
        std::tuple<size_t, size_t, size_t> sizes;
        if (op == 0) {
            sizes = std::make_tuple(Q, p, q);
        } else if (op == 1) {
            sizes = std::make_tuple(p, Q, q);
        } else {
            sizes = std::make_tuple(p, q, Q);
        }

        sizes_[pfilename] = sizes;
        sizes_[filename] = sizes;
    }
}

void DFHelper::initialize() {
    if (debug_) {
        outfile->Printf("Entering DFHelper::initialize\n");
    }
    timer_on("DFH: initialize()");

    // have the algorithm specified before init
    if (method_.compare("DIRECT") && method_.compare("STORE") && method_.compare("DIRECT_iaQ")) {
        std::stringstream error;
        error << "DFHelper:initialize: specified method (" << method_ << ") is incorrect";
        throw PSIEXCEPTION(error.str().c_str());
    }

    // workflow holders
    direct_iaQ_ = (!method_.compare("DIRECT_iaQ") ? true : false);
    direct_ = (!method_.compare("DIRECT") ? true : false);

    // did we get enough memory for at least the metric?
    if (naux_ * naux_ > memory_) {
        std::stringstream error;
        error << "DFHelper: The Coulomb metric requires at least " << naux_ * naux_ * 8 / (1024 * 1024 * 1024.0)
              << "[GiB].  We need that plus some more, but we only got " << memory_ * 8 / (1024 * 1024 * 1024.0)
              << "[GiB].";
        throw PSIEXCEPTION(error.str().c_str());
    }

    // if metric power is not zero, prepare it
    if (!(std::fabs(mpower_ - 0.0) < 1e-13)) (hold_met_ ? prepare_metric_core() : prepare_metric());

    // if metric power for omega integrals is not zero, prepare its metric
    if (do_wK_) {
        if (!(std::fabs(wmpower_ - 0.0) < 1e-13) && (std::fabs(mpower_ - 0.0) < 1e-13))
            (hold_met_ ? prepare_metric_core() : prepare_metric());
    }

    // prepare sparsity masks
    prepare_sparsity();

    // figure out AO_core
    AO_core();
    if (print_lvl_ > 0) {
        outfile->Printf("  DFHelper Memory: AOs need %.3f GiB; user supplied %.3f GiB. ",
                        (required_core_size_ * 8 / (1024 * 1024 * 1024.0)), (memory_ * 8 / (1024 * 1024 * 1024.0)));
        outfile->Printf("%s in-core AOs.\n\n", AO_core_ ? "Using" : "Turning off");
    }

    // prepare AOs for STORE method
    if (AO_core_) {
        if (do_wK_) {
            prepare_AO_wK_core();
        } else {  // It is possible to reformulate the expression for the
            //   coulomb matrix to save memory in case do_wK_ is
            //   is true, but do_K_ is false. This code isn't written
            prepare_AO_core();
        }
    } else if (!direct_ && !direct_iaQ_) {
        prepare_AO();
        if (do_wK_) {
            std::stringstream error;
            error << "DFHelper: not equipped to do wK out of core. \nPlease supply more memory or remove scf_type Mem_DF from the imput file";
            throw PSIEXCEPTION(error.str().c_str());
            //prepare_AO_wK();
        }
    }

    built_ = true;
    timer_off("DFH: initialize()");

    if (debug_) {
        outfile->Printf("Exiting DFHelper::initialize\n");
    }
}
void DFHelper::AO_core() {
    prepare_sparsity();

    if (direct_iaQ_) {
        // the direct_iaQ method does not use sparse storage
        // if do_wK added to code, the following will need to be changed to match
        required_core_size_ = naux_ * nbf_ * nbf_;
    } else {
        // total size of sparse AOs.
        // yes, a nested ternary operator
        required_core_size_ = (do_wK_ ? ( wcombine_ ? 2 * big_skips_[nbf_] : 3 * big_skips_[nbf_] ) : big_skips_[nbf_]);
    }

    // Auxiliary metric
    required_core_size_ += naux_ * naux_;

    // C_buffers (conservative estimate since I do not have max_nocc TODO)
    required_core_size_ += nthreads_ * nbf_ * nbf_;

    // Tmp buffers
    required_core_size_ += 3 * nbf_ * nbf_ * Qshell_max_;

    // a fraction of memory to use, do we want it as an option?
    AO_core_ = true;
    if (memory_ < required_core_size_) AO_core_ = false;
}
void DFHelper::print_header() {
    // Preps any required metadata, safe to call multiple times
    get_core_size();

    outfile->Printf("  ==> DFHelper <==\n");
    outfile->Printf("    NBF:                     %11ld\n", nbf_);
    outfile->Printf("    NAux:                    %11ld\n", naux_);
    outfile->Printf("    Schwarz Cutoff:          %11.0E\n", cutoff_);
    outfile->Printf("    Mask sparsity (%%):       %11.0f\n", 100. * ao_sparsity());
    outfile->Printf("    DFH Avail. Memory [GiB]: %11.3f\n", (memory_ * 8L) / ((double)(1024L * 1024L * 1024L)));
    outfile->Printf("    OpenMP threads:          %11zu\n", nthreads_);
    outfile->Printf("    Algorithm:               %11s\n", method_.c_str());
    outfile->Printf("    AO Core:                 %11s\n", (AO_core_ ? "True" : "False"));
    outfile->Printf("    MO Core:                 %11s\n", (MO_core_ ? "True" : "False"));
    outfile->Printf("    Hold Metric:             %11s\n", (hold_met_ ? "True" : "False"));
    outfile->Printf("    Metric Power:            %11.3f\n", mpower_);
    outfile->Printf("    Fitting Condition:       %11.0E\n", condition_);
    outfile->Printf("    Q Shell Max:             %11d\n", (int)Qshell_max_);
    outfile->Printf("\n\n");
}

void DFHelper::prepare_sparsity() {
    if (sparsity_prepared_) return;
    timer_on("DFH: sparsity prep");

    // => Initialize vectors <=
    std::vector<double> shell_max_vals(pshells_ * pshells_, 0.0);
    std::vector<double> fun_max_vals(nbf_ * nbf_, 0.0);
    schwarz_shell_mask_.resize(pshells_ * pshells_);
    schwarz_fun_index_.resize(nbf_ * nbf_);
    symm_ignored_columns_.resize(nbf_);
    symm_big_skips_.resize(nbf_ + 1);
    symm_small_skips_.resize(nbf_);
    small_skips_.resize(nbf_ + 1);
    big_skips_.resize(nbf_ + 1);

    // => Populate a vector of TwoBodyAOInt to make ERIs, one per screening thread. <=
    size_t screen_threads = (nthreads_ == 1 ? 1 : 2);  // TODO: Replace screen_threads with nthreads_?
    auto rifactory = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(screen_threads);
    eri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
#pragma omp parallel num_threads(screen_threads) if (nbf_ > 1000)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        if (rank) eri[rank] = std::shared_ptr<TwoBodyAOInt>(eri.front()->clone());
    }

    // => For each shell pair and basis pair, store the max (mn|mn)-type integral for screening. <=
    double max_val = 0.0;
#pragma omp parallel for num_threads(screen_threads) if (nbf_ > 1000) schedule(guided) reduction(max : max_val)
    for (size_t MU = 0; MU < pshells_; ++MU) {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        const auto& buffers = eri[rank]->buffers();
        size_t nmu = primary_->shell(MU).nfunction();
        for (size_t NU = 0; NU <= MU; ++NU) {
            size_t nnu = primary_->shell(NU).nfunction();
            eri[rank]->compute_shell(MU, NU, MU, NU);
            const auto *buffer = buffers[0];
            // Loop over basis functions inside shell pair
            for (size_t mu = 0; mu < nmu; ++mu) {
                size_t omu = primary_->shell(MU).function_index() + mu;
                for (size_t nu = 0; nu < nnu; ++nu) {
                    size_t onu = primary_->shell(NU).function_index() + nu;

                    // Find shell and function maximums
                    if (omu >= onu) {
                        size_t index = mu * (nnu * nmu * nnu + nnu) + nu * (nmu * nnu + 1);
                        double val = fabs(buffer[index]);
                        max_val = std::max(val, max_val);
                        if (shell_max_vals[MU * pshells_ + NU] <= val) {
                            shell_max_vals[MU * pshells_ + NU] = val;
                            shell_max_vals[NU * pshells_ + MU] = val;
                        }
                        if (fun_max_vals[omu * nbf_ + onu] <= val) {
                            fun_max_vals[omu * nbf_ + onu] = val;
                            fun_max_vals[onu * nbf_ + omu] = val;
                        }
                    }
                }
            }
        }
    }

    // => Prepare screening/indexing data <=
    double tolerance = cutoff_ * cutoff_ / max_val;

    // ==> Is this shell pair significant? <==
    for (size_t i = 0; i < pshells_ * pshells_; i++) schwarz_shell_mask_[i] = (shell_max_vals[i] >= tolerance);
    
    // ==> Is this basis function pair significant? Also start storing non-symmetric indexing. <==
    for (size_t i = 0, count = 0; i < nbf_; i++) {
        count = 0;
        for (size_t j = 0; j < nbf_; j++) {
            if (fun_max_vals[i * nbf_ + j] >= tolerance) {
                count++;
                schwarz_fun_index_[i * nbf_ + j] = count;
            } else
                schwarz_fun_index_[i * nbf_ + j] = 0;
        }
        small_skips_[i] = count;
    }

    // ==> Non-symmetric indexing. <==
    big_skips_[0] = 0;
    size_t coltots = 0;
    for (size_t j = 0; j < nbf_; j++) {
        size_t cols = small_skips_[j];
        size_t size = cols * naux_;
        coltots += cols;
        big_skips_[j + 1] = size + big_skips_[j];
    }
    small_skips_[nbf_] = coltots;

    // ==> Symmetric indexing. <==
    for (size_t i = 0; i < nbf_; i++) {
        size_t size = 0;
        size_t skip = 0;
        for (size_t j = 0; j < nbf_; j++) {
            if (schwarz_fun_index_[i * nbf_ + j]) {
                (j >= i ? size++ : skip++);
            }
        }
        symm_small_skips_[i] = size;
        symm_ignored_columns_[i] = skip;
    }

    symm_big_skips_[0] = 0;
    for (size_t i = 1; i < nbf_ + 1; i++) {
        symm_big_skips_[i] = symm_big_skips_[i - 1] + symm_small_skips_[i - 1] * naux_;
    }

    sparsity_prepared_ = true;
    timer_off("DFH: sparsity prep");
}

void DFHelper::prepare_AO() {
    // prepare eris
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    auto rifactory = std::make_shared<IntegralFactory>(aux_, zero, primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthreads_);
    eri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
#pragma omp parallel num_threads(nthreads_)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(eri.front()->clone());
    }

    // gather blocking info
    std::vector<std::pair<size_t, size_t>> psteps;
    std::pair<size_t, size_t> plargest = pshell_blocks_for_AO_build(memory_, 0, psteps);

    // declare largest necessary
    std::unique_ptr<double[]> M(new double[std::get<0>(plargest) / 2]);  // there was a factor of two built in
    std::unique_ptr<double[]> F(new double[std::get<0>(plargest) / 2]);
    std::unique_ptr<double[]> metric;
    double* Mp = M.get();
    double* Fp = F.get();

    // grab metric
    double* metp;
    if (!hold_met_) {
        metric = std::unique_ptr<double[]>(new double[naux_ * naux_]);
        metp = metric.get();
        std::string filename = return_metfile(mpower_);
        get_tensor_(std::get<0>(files_[filename]), metp, 0, naux_ - 1, 0, naux_ - 1);

    } else
        metp = metric_prep_core(mpower_);

    // prepare files
    AO_filename_maker(1);
    AO_filename_maker(2);
    std::string putf = AO_names_[1];
    std::string op = "ab";

    // Contract metric according to previously calculated scheme
    size_t count = 0;
    for (size_t i = 0; i < psteps.size(); i++) {
        // setup
        size_t start = std::get<0>(psteps[i]);
        size_t stop = std::get<1>(psteps[i]);
        size_t begin = pshell_aggs_[start];
        size_t end = pshell_aggs_[stop + 1] - 1;
        size_t block_size = end - begin + 1;
        size_t size = big_skips_[end + 1] - big_skips_[begin];

        // compute
        timer_on("DFH: Total Workflow");
        timer_on("DFH: AO Construction");
        compute_sparse_pQq_blocking_p(start, stop, Mp, eri);
        timer_off("DFH: AO Construction");

        // loop and contract
        timer_on("DFH: AO-Met. Contraction");
#pragma omp parallel for num_threads(nthreads_) schedule(guided)
        for (size_t j = 0; j < block_size; j++) {
            size_t mi = small_skips_[begin + j];
            size_t skips = big_skips_[begin + j] - big_skips_[begin];
            C_DGEMM('N', 'N', naux_, mi, naux_, 1.0, metp, naux_, &Mp[skips], mi, 0.0, &Fp[skips], mi);
        }
        timer_off("DFH: AO-Met. Contraction");
        timer_off("DFH: Total Workflow");

        // put
        put_tensor_AO(putf, Fp, size, count, op);
        count += size;
    }
}
void DFHelper::prepare_AO_wK() {
    // prepare eris
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    auto rifactory = std::make_shared<IntegralFactory>(aux_, zero, primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthreads_);
    eri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
#pragma omp parallel num_threads(nthreads_)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(eri.front()->clone());
    }

    // gather blocking info
    std::vector<std::pair<size_t, size_t>> psteps;
    std::pair<size_t, size_t> plargest = pshell_blocks_for_AO_build(memory_, 0, psteps);
}
void DFHelper::prepare_AO_core() {
    // get each thread an eri object
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    auto rifactory = std::make_shared<IntegralFactory>(aux_, zero, primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthreads_);
    eri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
#pragma omp parallel num_threads(nthreads_)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        if(rank) eri[rank] = std::shared_ptr<TwoBodyAOInt>(eri.front()->clone());
    }

    // determine blocking
    std::vector<std::pair<size_t, size_t>> psteps;
    std::pair<size_t, size_t> plargest = pshell_blocks_for_AO_build(memory_, 1, psteps);

    // allocate final AO vector
    if (direct_iaQ_) {
        Ppq_ = std::unique_ptr<double[]>(new double[naux_ * nbf_ * nbf_]);
    } else {
        Ppq_ = std::unique_ptr<double[]>(new double[big_skips_[nbf_]]);
    }

    double* ppq = Ppq_.get();

    // outfile->Printf("\n    ==> Begin AO Blocked Construction <==\n\n");
    if (direct_iaQ_ || direct_) {
        timer_on("DFH: AO Construction");
        if (direct_iaQ_) {
            compute_dense_Qpq_blocking_Q(0, Qshells_ - 1, &Ppq_[0], eri);
        } else {
            compute_sparse_pQq_blocking_p(0, pshells_ - 1, &Ppq_[0], eri);
        }
        timer_off("DFH: AO Construction");

    } else {
        // declare sparse buffer
        std::unique_ptr<double[]> Qpq(new double[std::get<0>(plargest)]);
        double* Mp = Qpq.get();
        std::unique_ptr<double[]> metric;
        double* metp;

        if (!hold_met_) {
            metric = std::unique_ptr<double[]>(new double[naux_ * naux_]);
            metp = metric.get();
            std::string filename = return_metfile(mpower_);
            get_tensor_(std::get<0>(files_[filename]), metp, 0, naux_ - 1, 0, naux_ - 1);
        } else
            metp = metric_prep_core(mpower_);

        for (size_t i = 0; i < psteps.size(); i++) {
            size_t start = std::get<0>(psteps[i]);
            size_t stop = std::get<1>(psteps[i]);
            size_t begin = pshell_aggs_[start];
            size_t end = pshell_aggs_[stop + 1] - 1;

            // compute
            timer_on("DFH: AO Construction");
            compute_sparse_pQq_blocking_p_symm(start, stop, Mp, eri);
            timer_off("DFH: AO Construction");

            // contract metric
            timer_on("DFH: AO-Met. Contraction");
            contract_metric_AO_core_symm(Mp, ppq, metp, begin, end);
            timer_off("DFH: AO-Met. Contraction");
        }
        // no more need for metrics
        if (hold_met_) metrics_.clear();
    }
    // outfile->Printf("\n    ==> End AO Blocked Construction <==");
}
void DFHelper::prepare_AO_wK_core() {
    // get each thread an eri object
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    auto rifactory = std::make_shared<IntegralFactory>(aux_, zero, primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthreads_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> weri(nthreads_);

    eri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    weri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_));
#pragma omp parallel num_threads(nthreads_)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        if (rank) {
            eri[rank] = std::shared_ptr<TwoBodyAOInt>(eri.front()->clone());
            weri[rank] = std::shared_ptr<TwoBodyAOInt>(weri.front()->clone());
        }
    }

    // use blocking as for prepare_AO_core
    std::vector<std::pair<size_t, size_t>> psteps;
    std::pair<size_t, size_t> plargest = pshell_blocks_for_AO_build(memory_, 1, psteps);

    wPpq_ = std::make_unique<double[]>(big_skips_[nbf_]);
    m1Ppq_ = std::make_unique<double[]>(big_skips_[nbf_]);

    if (!wcombine_) Ppq_ = std::make_unique<double[]>(big_skips_[nbf_]);


    double* wppq = wPpq_.get();
    double* ppq = nullptr;

    if (!wcombine_) ppq = Ppq_.get();

    double* m1ppq = m1Ppq_.get();

    std::unique_ptr<double[]> Qpq(new double[std::get<0>(plargest)]);
    std::unique_ptr<double[]> Qpq2(new double[std::get<0>(plargest)]);
    double* M1p = Qpq.get();
    double* M2p = Qpq2.get();
    std::unique_ptr<double[]> metric1;
    double* met1p = nullptr;
    std::unique_ptr<double[]> metric;
    double* metp = nullptr;

    if (!hold_met_) {
        metric1 = std::make_unique<double[]>(naux_ * naux_);
        met1p = metric1.get();
        std::string filename = return_metfile(wmpower_);
        get_tensor_(std::get<0>(files_[filename]), met1p, 0, naux_ - 1, 0, naux_ - 1);
        metric = std::make_unique<double[]>(naux_ * naux_);
        metp = metric.get();
        filename = return_metfile(mpower_);
        get_tensor_(std::get<0>(files_[filename]), metp, 0, naux_ - 1, 0, naux_ - 1);
    } else {
        met1p = metric_prep_core(wmpower_);
        metp = metric_prep_core(mpower_);
    }

    for (size_t i = 0; i < psteps.size(); i++) {
        size_t start = std::get<0>(psteps[i]);
        size_t stop = std::get<1>(psteps[i]);
        size_t begin = pshell_aggs_[start];
        size_t end = pshell_aggs_[stop + 1] - 1;

        if ( wcombine_ ) {
            // computes (Q|mn) and (Q|w|mn)
            timer_on("DFH: AO Construction");
            compute_sparse_pQq_blocking_p_symm_abw(start,stop, M1p, M2p, eri, weri);
            timer_off("DFH: AO Construction");
            // computes  [J^{-1.0}](Q|mn)
            timer_on("DFH: AO-Met. Contraction");
            contract_metric_AO_core_symm(M1p, m1ppq, met1p, begin, end);
            timer_off("DFH: AO-Met. Contraction");

            copy_upper_lower_wAO_core_symm(M2p, wppq, begin, end);
        } else {

            timer_on("DFH: AO Construction");
            compute_sparse_pQq_blocking_p_symm(start, stop, M1p, eri);
            timer_off("DFH: AO Construction");

            // contract full metric inverse computes  [J^{-1.0}](Q|mn)
            timer_on("DFH: AO-Met. Contraction");
            contract_metric_AO_core_symm(M1p, m1ppq, met1p, begin, end);
            timer_off("DFH: AO-Met. Contraction");
    
            // contract half metric inverse
            timer_on("DFH: AO-Met. Contraction");
            contract_metric_AO_core_symm(M1p, ppq, metp, begin, end);
            timer_off("DFH: AO-Met. Contraction");
    
            // compute (A|w|mn)
            timer_on("DFH: wAO Construction");
            compute_sparse_pQq_blocking_p_symm(start, stop, M1p, weri);
            timer_off("DFH: wAO Construction");
    
            timer_on("DFH: wAO Copy");
            copy_upper_lower_wAO_core_symm(M1p, wppq, begin, end);
            timer_off("DFH: wAO Copy");

        }
    }
    // no more need for metrics
    if (hold_met_) metrics_.clear();
}
std::pair<size_t, size_t> DFHelper::pshell_blocks_for_AO_build(const size_t mem, size_t symm,
                                                               std::vector<std::pair<size_t, size_t>>& b) {
    size_t full_3index = (symm ? big_skips_[nbf_] : 0);
    size_t constraint, end, begin, current, block_size, tmpbs, total, count, largest;
    block_size = tmpbs = total = count = largest = 0;
    for (size_t i = 0; i < pshells_; i++) {
        count++;
        begin = pshell_aggs_[i];
        end = pshell_aggs_[i + 1] - 1;
        tmpbs += end - begin + 1;

        if (symm) {
            // in-core symmetric
            // get current cost of this block of AOs and add it to the total
            // the second buffer is accounted for with full AO_core
            current = symm_big_skips_[end + 1] - symm_big_skips_[begin];
            if (do_wK_) {
                if (wcombine_) {
                    current *= 2; 
                }
                else {
                    current *= 3;
                } 
            }
            total += current;
        } else {
            // on-disk
            // get current cost of this block of AOs and add it to the total
            // count current twice, for both pre and post contracted buffers
            current = big_skips_[end + 1] - big_skips_[begin];
            if (do_wK_) {
                current *= 3;
            }
            total += 2 * current;
        }

        constraint = total;
        constraint += full_3index;
        constraint += (hold_met_ ? naux_ * naux_ : total);
        if (constraint > mem || i == pshells_ - 1) {
            if (count == 1 && i != pshells_ - 1) {
                std::stringstream error;
                error << "DFHelper: not enough memory for (p shell) AO blocking!"
                      << " required memory: " << constraint * 8 / (1024 * 1024 * 1024.0) << " [GiB].";
                throw PSIEXCEPTION(error.str().c_str());
            }
            if (constraint > mem) {
                total -= current;
                tmpbs -= end - begin + 1;
                b.push_back(std::make_pair(i - count + 1, i - 1));
                i--;
            } else if (i == pshells_ - 1)
                b.push_back(std::make_pair(i - count + 1, i));
            if (largest < total) {
                largest = total;
                block_size = tmpbs;
            }
            count = 0;
            total = 0;
            tmpbs = 0;
        }
    }
    // returns pair(largest buffer size, largest block size)
    return std::make_pair(largest, block_size);
}
std::pair<size_t, size_t> DFHelper::Qshell_blocks_for_transform(const size_t mem, size_t wtmp, size_t wfinal,
                                                                std::vector<std::pair<size_t, size_t>>& b) {
    size_t extra = (hold_met_ ? naux_ * naux_ : 0);
    size_t end, begin, current, block_size, tmpbs, total, count, largest;
    block_size = tmpbs = total = count = largest = 0;
    for (size_t i = 0; i < Qshells_; i++) {
        count++;
        begin = Qshell_aggs_[i];
        end = Qshell_aggs_[i + 1] - 1;
        tmpbs += end - begin + 1;

        if (direct_iaQ_) {
            // the direct_iaQ method does not use sparse storage
            current = (end - begin + 1) * nbf_ * nbf_;
            total += current;
            total = (AO_core_ ? naux_ * nbf_ * nbf_ : total);
        } else {
            current = (end - begin + 1) * small_skips_[nbf_];
            total += current;
            total = (AO_core_ ? big_skips_[nbf_] : total);
        }

        size_t constraint = total + (wtmp * nbf_ + 2 * wfinal) * tmpbs + extra;
        // AOs + worst half transformed + worst final
        if (constraint > mem || i == Qshells_ - 1) {
            if (count == 1 && i != Qshells_ - 1) {
                std::stringstream error;
                error << "DFHelper: not enough memory for transformation blocking!";
                throw PSIEXCEPTION(error.str().c_str());
            }
            if (constraint > mem) {
                if (!AO_core_) total -= current;
                tmpbs -= end - begin + 1;
                b.push_back(std::make_pair(i - count + 1, i - 1));
                i--;
            } else if (i == Qshells_ - 1)
                b.push_back(std::make_pair(i - count + 1, i));
            if (block_size < tmpbs) {
                block_size = tmpbs;
                largest = total;
            }
            count = 0;
            total = 0;
            tmpbs = 0;
        }
    }

    return std::make_pair(largest, block_size);
}
std::tuple<size_t, size_t> DFHelper::Qshell_blocks_for_JK_build(std::vector<std::pair<size_t, size_t>>& b,
                                                                size_t max_nocc, bool lr_symmetric) {
    // strategy here:
    // 1. depending on lr_symmetric, T2 can either be the same as T1 or
    // it can just be used as a Jtmp.
    // 2. T3 is always used, includes C_buffers

    // K tmps
    size_t T1 = nbf_ * max_nocc;
    size_t T2 = (lr_symmetric ? nbf_ * nbf_ : nbf_ * max_nocc);

    // C_buffers
    size_t T3 = std::max(nthreads_ * nbf_ * nbf_, nthreads_ * nbf_ * max_nocc);

    // total AO buffer size is max if core alg is used, otherwise init to 0
    size_t total_AO_buffer = (AO_core_ ? big_skips_[nbf_] : 0);

    size_t block_size = 0, largest = 0;
    for (size_t i = 0, tmpbs = 0, count = 1; i < Qshells_; i++, count++) {
        // get shell info
        size_t begin = Qshell_aggs_[i];
        size_t end = Qshell_aggs_[i + 1] - 1;

        // update AO buffer, block sizes
        size_t current = (end - begin + 1) * small_skips_[nbf_];
        total_AO_buffer += (AO_core_ ? 0 : current);
        tmpbs += end - begin + 1;

        // compute total memory used by aggregate block
        size_t constraint = total_AO_buffer + T1 * tmpbs + T3;
        constraint += (lr_symmetric ? T2 : T2 * tmpbs);

        if (constraint > memory_ || i == Qshells_ - 1) {
            if (count == 1 && i != Qshells_ - 1) {
                std::stringstream error;
                error << "DFHelper: not enough memory for JK blocking!";
                throw PSIEXCEPTION(error.str().c_str());
            }
            if (constraint > memory_) {
                total_AO_buffer -= current;
                tmpbs -= end - begin + 1;
                b.push_back(std::make_pair(i - count + 1, i - 1));
                i--;
            } else if (i == Qshells_ - 1) {
                b.push_back(std::make_pair(i - count + 1, i));
            }
            if (block_size < tmpbs) {
                largest = total_AO_buffer;
                block_size = tmpbs;
            }
            count = total_AO_buffer = tmpbs = 0;
        }
    }

    return std::make_tuple(largest, block_size);
}

FILE* DFHelper::stream_check(std::string filename, std::string op) {
    if (file_streams_.count(filename) == 0) {
        file_streams_[filename] = std::make_shared<Stream>(filename, op);
    }

    return file_streams_[filename]->get_stream(op);
}

DFHelper::StreamStruct::StreamStruct(std::string filename, std::string op, bool activate) {
    op_ = op;
    filename_ = filename;
    if (activate) {
        fp_ = fopen(filename.c_str(), op_.c_str());
        open_ = true;
    }
}

DFHelper::StreamStruct::StreamStruct() {}

DFHelper::StreamStruct::~StreamStruct() {
    fflush(fp_);
    fclose(fp_);
    std::remove(filename_.c_str());
}

FILE* DFHelper::StreamStruct::get_stream(std::string op) {
    if (op.compare(op_)) {
        change_stream(op);
    } else {
        if (!open_) {
            fp_ = fopen(filename_.c_str(), op_.c_str());
            open_ = true;
        }
    }

    return fp_;
}

void DFHelper::StreamStruct::change_stream(std::string op) {
    if (open_) {
        close_stream();
    }
    op_ = op;
    fp_ = fopen(filename_.c_str(), op_.c_str());
}

void DFHelper::StreamStruct::close_stream() {
    fflush(fp_);
    fclose(fp_);
}

void DFHelper::put_tensor(std::string file, double* b, std::pair<size_t, size_t> i0, std::pair<size_t, size_t> i1,
                          std::pair<size_t, size_t> i2, std::string op) {
    // collapse to 2D, assume file has form (i1 | i2 i3)
    size_t A2 = std::get<2>(sizes_[file]);

    size_t sta0 = std::get<0>(i0);
    size_t sto0 = std::get<1>(i0);
    size_t sta1 = std::get<0>(i1);
    size_t sto1 = std::get<1>(i1);
    size_t sta2 = std::get<0>(i2);
    size_t sto2 = std::get<1>(i2);

    size_t a0 = sto0 - sta0 + 1;
    size_t a1 = sto1 - sta1 + 1;
    size_t a2 = sto2 - sta2 + 1;

    // check contiguity (a2)
    if (A2 == a2) {
        put_tensor(file, b, sta0, sto0, a2 * sta1, a2 * (sto1 + 1) - 1, op);
    } else {  // loop (a0, a1)
        for (size_t j = 0; j < a0; j++) {
            for (size_t i = 0; i < a1; i++) {
                put_tensor(file, &b[j * (a1 * a2) + i * a2], sta0 + j, sta0 + j, (i + sta1) * A2 + sta2,
                           (i + sta1) * A2 + sta2 + a2 - 1, op);
            }
        }
    }
}
void DFHelper::put_tensor(std::string file, double* Mp, const size_t start1, const size_t stop1, const size_t start2,
                          const size_t stop2, std::string op) {
    size_t a0 = stop1 - start1 + 1;
    size_t a1 = stop2 - start2 + 1;
    size_t A0 = std::get<0>(sizes_[file]);
    size_t A1 = std::get<1>(sizes_[file]) * std::get<2>(sizes_[file]);
    size_t st = A1 - a1;

    // begin stream
    FILE* fp = stream_check(file, op);

    // adjust position
    fseek(fp, (start1 * A1 + start2) * sizeof(double), SEEK_SET);

    // is everything contiguous?
    if (st == 0) {
        size_t s = fwrite(&Mp[0], sizeof(double), a0 * a1, fp);
        if (!s) {
            std::stringstream error;
            error << "DFHelper:put_tensor: write error";
            throw PSIEXCEPTION(error.str().c_str());
        }
    } else {
        for (size_t i = start1; i < stop1; i++) {
            // write
            size_t s = fwrite(&Mp[i * a1], sizeof(double), a1, fp);
            if (!s) {
                std::stringstream error;
                error << "DFHelper:put_tensor: write error";
                throw PSIEXCEPTION(error.str().c_str());
            }
            // advance stream
            fseek(fp, st * sizeof(double), SEEK_CUR);
        }
        // manual last one
        size_t s = fwrite(&Mp[(a0 - 1) * a1], sizeof(double), a1, fp);
        if (!s) {
            std::stringstream error;
            error << "DFHelper:put_tensor: write error";
            throw PSIEXCEPTION(error.str().c_str());
        }
    }
}
void DFHelper::put_tensor_AO(std::string file, double* Mp, size_t size, size_t start, std::string op) {
    // begin stream
    FILE* fp = stream_check(file, op);

    // adjust position
    fseek(fp, start, SEEK_SET);

    // everything is contiguous
    size_t s = fwrite(&Mp[0], sizeof(double), size, fp);
    if (!s) {
        std::stringstream error;
        error << "DFHelper:put_tensor_AO: write error";
        throw PSIEXCEPTION(error.str().c_str());
    }
}
void DFHelper::get_tensor_AO(std::string file, double* Mp, size_t size, size_t start) {
    // begin stream
    FILE* fp = stream_check(file, "rb");

    // adjust position
    fseek(fp, start * sizeof(double), SEEK_SET);

    // everything is contiguous
    size_t s = fread(&Mp[0], sizeof(double), size, fp);
    if (!s) {
        std::stringstream error;
        error << "DFHelper:get_tensor_AO: read error";
        throw PSIEXCEPTION(error.str().c_str());
    }
}
void DFHelper::get_tensor_(std::string file, double* b, std::pair<size_t, size_t> i0, std::pair<size_t, size_t> i1,
                           std::pair<size_t, size_t> i2) {
    // has this integral been transposed?
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(file) != tsizes_.end() ? tsizes_[file] : sizes_[file]);

    // collapse to 2D, assume file has form (i1 | i2 i3)
    size_t A2 = std::get<2>(sizes);

    size_t sta0 = std::get<0>(i0);
    size_t sto0 = std::get<1>(i0);
    size_t sta1 = std::get<0>(i1);
    size_t sto1 = std::get<1>(i1);
    size_t sta2 = std::get<0>(i2);
    size_t sto2 = std::get<1>(i2);

    size_t a0 = sto0 - sta0 + 1;
    size_t a1 = sto1 - sta1 + 1;
    size_t a2 = sto2 - sta2 + 1;

    // check contiguity (a2)
    if (A2 == a2) {
        get_tensor_(file, b, sta0, sto0, a2 * sta1, a2 * (sto1 + 1) - 1);
    } else {  // loop (a0, a1)
        for (size_t j = 0; j < a0; j++) {
            for (size_t i = 0; i < a1; i++) {
                get_tensor_(file, &b[j * (a1 * a2) + i * a2], sta0 + j, sta0 + j, (i + sta1) * A2 + sta2,
                            (i + sta1) * A2 + sta2 + a2 - 1);
            }
        }
    }
}
void DFHelper::get_tensor_(std::string file, double* b, const size_t start1, const size_t stop1, const size_t start2,
                           const size_t stop2) {
    size_t a0 = stop1 - start1 + 1;
    size_t a1 = stop2 - start2 + 1;

    // has this integral been transposed?
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(file) != tsizes_.end() ? tsizes_[file] : sizes_[file]);

    size_t A0 = std::get<0>(sizes);
    size_t A1 = std::get<1>(sizes) * std::get<2>(sizes);
    size_t st = A1 - a1;

    // check stream
    FILE* fp = stream_check(file, "rb");

    // adjust position
    fseek(fp, (start1 * A1 + start2) * sizeof(double), SEEK_SET);

    // is everything contiguous?
    if (st == 0) {
        size_t s = fread(&b[0], sizeof(double), a0 * a1, fp);
        if (!s) {
            std::stringstream error;
            error << "DFHelper:get_tensor: read error";
            throw PSIEXCEPTION(error.str().c_str());
        }
    } else {
        for (size_t i = 0; i < a0 - 1; i++) {
            // read
            size_t s = fread(&b[i * a1], sizeof(double), a1, fp);
            if (!s) {
                std::stringstream error;
                error << "DFHelper:get_tensor: read error";
                throw PSIEXCEPTION(error.str().c_str());
            }
            // advance stream
            s = fseek(fp, st * sizeof(double), SEEK_CUR);
            if (s) {
                std::stringstream error;
                error << "DFHelper:get_tensor: read error";
                throw PSIEXCEPTION(error.str().c_str());
            }
        }
        // manual last one
        size_t s = fread(&b[(a0 - 1) * a1], sizeof(double), a1, fp);
        if (!s) {
            std::stringstream error;
            error << "DFHelper:get_tensor: read error";
            throw PSIEXCEPTION(error.str().c_str());
        }
    }
}

void DFHelper::compute_dense_Qpq_blocking_Q(const size_t start, const size_t stop, double* Mp,
                                            std::vector<std::shared_ptr<TwoBodyAOInt>> eri) {
    // Here, we compute dense AO integrals in the Qpq memory layout.
    // Sparsity and permutational symmetry are used in the computation,
    // but not in the resulting tensor.

    size_t begin = Qshell_aggs_[start];
    size_t end = Qshell_aggs_[stop + 1] - 1;
    size_t block_size = end - begin + 1;

    // stripe the buffer
    fill(Mp, block_size * nbf_ * nbf_, 0.0);

    // prepare eri buffers
    size_t nthread = nthreads_;
    if (eri.size() != nthreads_) nthread = eri.size();
    std::vector<const double*> buffer(nthread);
#pragma omp parallel num_threads(nthread)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        buffer[rank] = eri[rank]->buffer();
    }

#pragma omp parallel for schedule(guided) num_threads(nthreads_)
    for (size_t MU = 0; MU < pshells_; MU++) {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        size_t nummu = primary_->shell(MU).nfunction();
        for (size_t NU = 0; NU < pshells_; NU++) {
            size_t numnu = primary_->shell(NU).nfunction();
            if (!schwarz_shell_mask_[MU * pshells_ + NU]) {
                continue;
            }
            for (size_t Pshell = start; Pshell <= stop; Pshell++) {
                size_t PHI = aux_->shell(Pshell).function_index();
                size_t numP = aux_->shell(Pshell).nfunction();
                eri[rank]->compute_shell(Pshell, 0, MU, NU);
                buffer[rank] = eri[rank]->buffer();
                for (size_t mu = 0; mu < nummu; mu++) {
                    size_t omu = primary_->shell(MU).function_index() + mu;
                    for (size_t nu = 0; nu < numnu; nu++) {
                        size_t onu = primary_->shell(NU).function_index() + nu;
                        if (!schwarz_fun_index_[omu * nbf_ + onu]) {
                            continue;
                        }
                        for (size_t P = 0; P < numP; P++) {
                            Mp[(PHI + P - begin) * nbf_ * nbf_ + omu * nbf_ + onu] =
                                Mp[(PHI + P - begin) * nbf_ * nbf_ + onu * nbf_ + omu] =
                                    buffer[rank][P * nummu * numnu + mu * numnu + nu];
                        }
                    }
                }
            }
        }
    }
}

void DFHelper::compute_sparse_pQq_blocking_Q(const size_t start, const size_t stop, double* Mp,
                                             std::vector<std::shared_ptr<TwoBodyAOInt>> eri) {
    size_t begin = Qshell_aggs_[start];
    size_t end = Qshell_aggs_[stop + 1] - 1;
    size_t block_size = end - begin + 1;

    // prepare eri buffers
    size_t nthread = nthreads_;
    if (eri.size() != nthreads_) nthread = eri.size();

    std::vector<const double*> buffer(nthread);
#pragma omp parallel num_threads(nthread)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        buffer[rank] = eri[rank]->buffer();
    }

#pragma omp parallel for schedule(guided) num_threads(nthreads_)
    for (size_t MU = 0; MU < pshells_; MU++) {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        size_t nummu = primary_->shell(MU).nfunction();
        for (size_t NU = 0; NU < pshells_; NU++) {
            size_t numnu = primary_->shell(NU).nfunction();
            if (!schwarz_shell_mask_[MU * pshells_ + NU]) {
                continue;
            }
            for (size_t Pshell = start; Pshell <= stop; Pshell++) {
                size_t PHI = aux_->shell(Pshell).function_index();
                size_t numP = aux_->shell(Pshell).nfunction();
                eri[rank]->compute_shell(Pshell, 0, MU, NU);
                buffer[rank] = eri[rank]->buffer();
                for (size_t mu = 0; mu < nummu; mu++) {
                    size_t omu = primary_->shell(MU).function_index() + mu;
                    for (size_t nu = 0; nu < numnu; nu++) {
                        size_t onu = primary_->shell(NU).function_index() + nu;
                        if (!schwarz_fun_index_[omu * nbf_ + onu]) {
                            continue;
                        }
                        for (size_t P = 0; P < numP; P++) {
                            Mp[(big_skips_[omu] * block_size) / naux_ + (PHI + P - begin) * small_skips_[omu] +
                               schwarz_fun_index_[omu * nbf_ + onu] - 1] =
                                buffer[rank][P * nummu * numnu + mu * numnu + nu];
                        }
                    }
                }
            }
        }
    }
}
void DFHelper::compute_sparse_pQq_blocking_p(const size_t start, const size_t stop, double* Mp,
                                             std::vector<std::shared_ptr<TwoBodyAOInt>> eri) {
    size_t begin = pshell_aggs_[start];
    size_t end = pshell_aggs_[stop + 1] - 1;
    size_t block_size = end - begin + 1;
    size_t startind = big_skips_[begin];
    //    outfile->Printf("      MU shell: (%zu, %zu)", start, stop);
    //    outfile->Printf(", nbf_ index: (%zu, %zu), size: %zu\n", begin, end, block_size);

    // prepare eri buffers
    size_t nthread = nthreads_;
    if (eri.size() != nthreads_) nthread = eri.size();

    std::vector<const double*> buffer(nthread);
#pragma omp parallel num_threads(nthread)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        buffer[rank] = eri[rank]->buffer();
    }

#pragma omp parallel for schedule(guided) num_threads(nthread)
    for (size_t MU = start; MU <= stop; MU++) {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        size_t nummu = primary_->shell(MU).nfunction();
        for (size_t NU = 0; NU < pshells_; NU++) {
            size_t numnu = primary_->shell(NU).nfunction();
            if (!schwarz_shell_mask_[MU * pshells_ + NU]) {
                continue;
            }
            for (size_t Pshell = 0; Pshell < Qshells_; Pshell++) {
                size_t PHI = aux_->shell(Pshell).function_index();
                size_t numP = aux_->shell(Pshell).nfunction();
                eri[rank]->compute_shell(Pshell, 0, MU, NU);
                buffer[rank] = eri[rank]->buffer();
                for (size_t mu = 0; mu < nummu; mu++) {
                    size_t omu = primary_->shell(MU).function_index() + mu;
                    for (size_t nu = 0; nu < numnu; nu++) {
                        size_t onu = primary_->shell(NU).function_index() + nu;
                        if (!schwarz_fun_index_[omu * nbf_ + onu]) {
                            continue;
                        }
                        for (size_t P = 0; P < numP; P++) {
                            Mp[big_skips_[omu] - startind + (PHI + P) * small_skips_[omu] +
                               schwarz_fun_index_[omu * nbf_ + onu] - 1] =
                                buffer[rank][P * nummu * numnu + mu * numnu + nu];
                        }
                    }
                }
            }
        }
    }
}
void DFHelper::compute_sparse_pQq_blocking_p_symm(const size_t start, const size_t stop, double* Mp,
                                                  std::vector<std::shared_ptr<TwoBodyAOInt>> eri) {
    size_t begin = pshell_aggs_[start];
    size_t end = pshell_aggs_[stop + 1] - 1;
    size_t block_size = end - begin + 1;
    size_t startind = symm_big_skips_[begin];
    //    outfile->Printf("      MU shell: (%zu, %zu)", start, stop);
    //    outfile->Printf(", nbf_ index: (%zu, %zu), size: %zu\n", begin, end, block_size);

    // prepare eri buffers
    size_t nthread = nthreads_;
    if (eri.size() != nthreads_) nthread = eri.size();

    std::vector<const double*> buffer(nthread);
#pragma omp parallel num_threads(nthread)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
    }

// Block over p in pQq 3-index integrals
#pragma omp parallel for schedule(guided) num_threads(nthread)
    for (size_t MU = start; MU <= stop; MU++) {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        size_t nummu = primary_->shell(MU).nfunction();
        for (size_t NU = MU; NU < pshells_; NU++) {
            size_t numnu = primary_->shell(NU).nfunction();
            if (!schwarz_shell_mask_[MU * pshells_ + NU]) {
                continue;
            }

            // Loop over Auxiliary index
            for (size_t Pshell = 0; Pshell < Qshells_; Pshell++) {
                size_t PHI = aux_->shell(Pshell).function_index();
                size_t numP = aux_->shell(Pshell).nfunction();
                eri[rank]->compute_shell(Pshell, 0, MU, NU);
                buffer[rank] = eri[rank]->buffer();

                for (size_t mu = 0; mu < nummu; mu++) {
                    size_t omu = primary_->shell(MU).function_index() + mu;
                    for (size_t nu = 0; nu < numnu; nu++) {
                        size_t onu = primary_->shell(NU).function_index() + nu;

                        // Remove sieved integrals or lower triangular
                        if (!schwarz_fun_index_[omu * nbf_ + onu] || omu > onu) {
                            continue;
                        }

                        for (size_t P = 0; P < numP; P++) {
                            size_t jump = schwarz_fun_index_[omu * nbf_ + onu] - schwarz_fun_index_[omu * nbf_ + omu];
                            size_t ind1 = symm_big_skips_[omu] - startind + (PHI + P) * symm_small_skips_[omu] + jump;
                            Mp[ind1] = buffer[rank][P * nummu * numnu + mu * numnu + nu];
                        }
                    }
                }
            }
        }
    } // ends MU loop
}

void DFHelper::compute_sparse_pQq_blocking_p_symm_abw(const size_t start, const size_t stop, double* just_Mp, double* param_Mp,
                                                  std::vector<std::shared_ptr<TwoBodyAOInt>> eri,
                                                  std::vector<std::shared_ptr<TwoBodyAOInt>> weri
) {
    size_t begin = pshell_aggs_[start];
    size_t end = pshell_aggs_[stop + 1] - 1;
    size_t block_size = end - begin + 1;
    size_t startind = symm_big_skips_[begin];
    //    outfile->Printf("      MU shell: (%zu, %zu)", start, stop);
    //    outfile->Printf(", nbf_ index: (%zu, %zu), size: %zu\n", begin, end, block_size);

    // prepare eri buffers
    size_t nthread = nthreads_;
    if (eri.size() != nthreads_) nthread = eri.size();

    std::vector<const double*> buffer(nthread);
    std::vector<const double*> wbuffer(nthread);
#pragma omp parallel num_threads(nthread)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
    }


#pragma omp parallel for schedule(guided) num_threads(nthread)
    for (size_t MU = start; MU <= stop; MU++) {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        size_t nummu = primary_->shell(MU).nfunction();
        for (size_t NU = MU; NU < pshells_; NU++) {
            size_t numnu = primary_->shell(NU).nfunction();
            if (!schwarz_shell_mask_[MU * pshells_ + NU]) {
                continue;
            }

            // Loop over Auxiliary index
            for (size_t Pshell = 0; Pshell < Qshells_; Pshell++) {
                size_t PHI = aux_->shell(Pshell).function_index();
                size_t numP = aux_->shell(Pshell).nfunction();
                eri[rank]->compute_shell(Pshell, 0, MU, NU);
                weri[rank]->compute_shell(Pshell, 0, MU, NU);
                buffer[rank] = eri[rank]->buffers()[0];
                wbuffer[rank] = weri[rank]->buffers()[0];

                for (size_t mu = 0; mu < nummu; mu++) {
                    size_t omu = primary_->shell(MU).function_index() + mu;
                    for (size_t nu = 0; nu < numnu; nu++) {
                        size_t onu = primary_->shell(NU).function_index() + nu;

                        // Remove sieved integrals or lower triangular
                        if (!schwarz_fun_index_[omu * nbf_ + onu] || omu > onu) {
                            continue;
                        }

                        for (size_t P = 0; P < numP; P++) {
                            size_t jump = schwarz_fun_index_[omu * nbf_ + onu] - schwarz_fun_index_[omu * nbf_ + omu];
                            size_t ind1 = symm_big_skips_[omu] - startind + (PHI + P) * symm_small_skips_[omu] + jump;
                            just_Mp[ind1] = buffer[rank][P * nummu * numnu + mu * numnu + nu];
                            param_Mp[ind1] = omega_alpha_ * buffer[rank][P * nummu * numnu + mu * numnu + nu] + omega_beta_ * wbuffer[rank][P * nummu * numnu + mu * numnu + nu];
                        }
                    }
                }
            }
        }
    } // ends MU loop
}
void DFHelper::grab_AO(const size_t start, const size_t stop, double* Mp) {
    size_t begin = Qshell_aggs_[start];
    size_t end = Qshell_aggs_[stop + 1] - 1;
    size_t block_size = end - begin + 1;
    auto getf = AO_names_[1];

    // presumably not thread safe or inherently sequential, but could revisit
    for (size_t i = 0, sta = 0; i < nbf_; i++) {
        size_t size = block_size * small_skips_[i];
        size_t jump = begin * small_skips_[i];
        get_tensor_AO(getf, &Mp[sta], size, big_skips_[i] + jump);
        sta += size;
    }
}
void DFHelper::prepare_metric_core() {
    timer_on("DFH: metric construction");
    FittingMetric J(aux_, true);
    J.form_fitting_metric();
    metrics_[1.0] = J.get_metric();
    timer_off("DFH: metric construction");
}
double* DFHelper::metric_prep_core(double m_pow) {
    bool on = false;
    double power;
    for (auto& kv : metrics_) {
        if (!(std::fabs(m_pow - kv.first) > 1e-13)) {
            on = true;
            power = kv.first;
            break;
        }
    }
    if (!on) {
        power = m_pow;
        SharedMatrix J = metrics_[1.0];
        if ( fabs(m_pow - 1.0) < 1e-13 ) {
            return J->pointer()[0];
        } else {
            J->power(power, condition_);
        }
        metrics_[power] = J;
    }
    return metrics_[power]->pointer()[0];
}
void DFHelper::prepare_metric() {
    // construct metric
    FittingMetric J(aux_, true);
    J.form_fitting_metric();
    auto metric = J.get_metric();
    auto Mp = metric->pointer()[0];

    // create file
    std::string filename = "metric.1.0";
    filename_maker(filename, naux_, naux_, 1);
    metric_keys_.emplace_back(1.0, filename);
    // store
    std::string putf = std::get<0>(files_[filename]);
    put_tensor(putf, Mp, 0, naux_ - 1, 0, naux_ - 1, "wb");
}
std::string DFHelper::return_metfile(double m_pow) {
    bool on = 0;
    std::string key;
    for (size_t i = 0; i < metric_keys_.size() && !on; i++) {
        double power = std::get<0>(metric_keys_[i]);
        if (std::fabs(power - m_pow) < 1e-12) {
            key = std::get<1>(metric_keys_[i]);
            on = 1;
        }
    }

    if (!on) key = compute_metric(m_pow);
    return key;
}
std::string DFHelper::compute_metric(double m_pow) {
    // ensure J
    if (std::fabs(m_pow - 1.0) < 1e-13)
        prepare_metric();
    else {
        // get metric
        auto metric = std::make_shared<Matrix>("met", naux_, naux_);
        double* metp = metric->pointer()[0];
        std::string filename = return_metfile(1.0);

        // get and compute
        get_tensor_(std::get<0>(files_[filename]), metp, 0, naux_ - 1, 0, naux_ - 1);
        metric->power(m_pow, condition_);

        // make new file
        std::string name = "metric";
        name.append(".");
        name.append(std::to_string(m_pow));
        filename_maker(name, naux_, naux_, 1);
        metric_keys_.push_back(std::make_pair(m_pow, name));

        // store
        std::string putf = std::get<0>(files_[name]);
        put_tensor(putf, metp, 0, naux_ - 1, 0, naux_ - 1, "wb");
    }
    return return_metfile(m_pow);
}
double* DFHelper::metric_inverse_prep_core() {
    bool on = false;
    double power;
    for (auto& kv : metrics_) {
        if ((std::fabs(-1.0 - kv.first) < 1e-13)) {
            on = true;
            power = kv.first;
            break;
        }
    }
    if (!on) {
        power = -1.0;
        timer_on("DFH: metric power");
        SharedMatrix J = metrics_[1.0];
        J->invert();
        timer_off("DFH: metric power");
    }

    return metrics_[power]->pointer()[0];
}
void DFHelper::metric_contraction_blocking(std::vector<std::pair<size_t, size_t>>& steps, size_t blocking_index,
                                           size_t block_sizes, size_t total_mem, size_t memory_factor,
                                           size_t memory_bump) {
    for (size_t i = 0, count = 1; i < blocking_index; i++, count++) {
        if (total_mem < count * block_sizes || i == blocking_index - 1) {
            if (count == 1 && i != blocking_index - 1) {
                std::stringstream error;
                error << "DFHelper:contract_metric: not enough memory, ";
                error << "needs at least "
                      << ((count * block_sizes) * memory_factor + memory_bump) / (1024 * 1024 * 1024.0) * 8. << "[GiB]";
                throw PSIEXCEPTION(error.str().c_str());
            }
            if (total_mem < count * block_sizes) {
                steps.push_back(std::make_pair(i - count + 1, i - 1));
                i--;
            } else {
                steps.push_back(std::make_pair(i - count + 1, i));
            }
            count = 0;
        }
    }
}

void DFHelper::contract_metric_Qpq(std::string file, double* metp, double* Mp, double* Fp, const size_t total_mem) {
    std::string getf = std::get<0>(files_[file]);
    std::string putf = std::get<1>(files_[file]);

    size_t Q = std::get<0>(sizes_[getf]);
    size_t l = std::get<1>(sizes_[getf]);
    size_t r = std::get<2>(sizes_[getf]);

    std::string op = "wb";
    std::vector<std::pair<size_t, size_t>> steps;
    metric_contraction_blocking(steps, l, Q * r, total_mem, 2, naux_ * naux_);

    for (size_t i = 0; i < steps.size(); i++) {
        size_t begin = std::get<0>(steps[i]);
        size_t end = std::get<1>(steps[i]);
        size_t bs = end - begin + 1;

        get_tensor_(getf, Mp, 0, Q - 1, begin * r, (end + 1) * r - 1);
        timer_on("DFH: Total Workflow");
        C_DGEMM('T', 'N', bs * r, Q, Q, 1.0, Mp, bs * r, metp, Q, 0.0, Fp, Q);
        timer_off("DFH: Total Workflow");
        put_tensor(putf, Fp, begin, end, 0, r * Q - 1, op);
    }
}

void DFHelper::contract_metric(std::string file, double* metp, double* Mp, double* Fp, const size_t total_mem) {
    std::string getf = std::get<0>(files_[file]);
    std::string putf = std::get<1>(files_[file]);
    size_t a0 = std::get<0>(sizes_[getf]);
    size_t a1 = std::get<1>(sizes_[getf]);
    size_t a2 = std::get<2>(sizes_[getf]);

    std::string op = "wb";
    std::vector<std::pair<size_t, size_t>> steps;

    // contract in steps
    if (std::get<2>(transf_[file])) {
        // determine blocking
        // both pqQ and pQq formats block through p, which is index 0
        metric_contraction_blocking(steps, a0, a1 * a2, total_mem, 2, naux_ * naux_);

        // grab val, the inner contractions are different depending on the form
        size_t val = std::get<2>(transf_[file]);
        for (size_t i = 0; i < steps.size(); i++) {
            size_t begin = std::get<0>(steps[i]);
            size_t end = std::get<1>(steps[i]);
            size_t bs = end - begin + 1;

            get_tensor_(getf, Mp, begin, end, 0, a1 * a2 - 1);
            timer_on("DFH: Total Workflow");

            if (val == 2) {
                C_DGEMM('N', 'N', bs * a1, a2, a2, 1.0, Mp, a2, metp, a2, 0.0, Fp, a2);
            } else {
#pragma omp parallel for num_threads(nthreads_)
                for (size_t i = 0; i < bs; i++) {
                    C_DGEMM('N', 'N', a1, a2, a1, 1.0, metp, a1, &Mp[i * a1 * a2], a2, 0.0, &Fp[i * a1 * a2], a2);
                }
            }
            timer_off("DFH: Total Workflow");
            put_tensor(putf, Fp, begin, end, 0, a1 * a2 - 1, op);
        }

    } else {
        // determine blocking
        // the Qpq format blocks through p, which is index 1
        metric_contraction_blocking(steps, a1, a0 * a2, total_mem, 2, naux_ * naux_);

        for (size_t i = 0; i < steps.size(); i++) {
            size_t begin = std::get<0>(steps[i]);
            size_t end = std::get<1>(steps[i]);
            size_t bs = end - begin + 1;

            get_tensor_(getf, Mp, 0, a0 - 1, begin * a2, (end + 1) * a2 - 1);
            timer_on("DFH: Total Workflow");
            C_DGEMM('N', 'N', a0, bs * a2, a0, 1.0, metp, a0, Mp, bs * a2, 0.0, Fp, bs * a2);
            timer_off("DFH: Total Workflow");
            put_tensor(putf, Fp, 0, a0 - 1, begin * a2, (end + 1) * a2 - 1, op);
        }
    }
}

void DFHelper::contract_metric_AO_core(double* Qpq, double* metp) {
// loop and contract
#pragma omp parallel for num_threads(nthreads_) schedule(guided)
    for (size_t j = 0; j < nbf_; j++) {
        size_t mi = small_skips_[j];
        size_t skips = big_skips_[j];
        C_DGEMM('N', 'N', naux_, mi, naux_, 1.0, metp, naux_, &Qpq[skips], mi, 0.0, &Ppq_[skips], mi);
    }
}

void DFHelper::contract_metric_AO_core_symm(double* Qpq, double* Ppq, double* metp, size_t begin, size_t end) {
    // loop and contract
    size_t startind = symm_big_skips_[begin];
#pragma omp parallel for num_threads(nthreads_) schedule(guided)
    for (size_t j = begin; j <= end; j++) {
        size_t mi = symm_small_skips_[j];
        size_t si = small_skips_[j];
        size_t jump = symm_ignored_columns_[j];
        size_t skip1 = big_skips_[j];
        size_t skip2 = symm_big_skips_[j] - startind;
        C_DGEMM('N', 'N', naux_, mi, naux_, 1.0, metp, naux_, &Qpq[skip2], mi, 0.0, &Ppq[skip1 + jump], si);
    }
// copy upper-to-lower
#pragma omp parallel for num_threads(nthreads_) schedule(static)
    for (size_t omu = begin; omu <= end; omu++) {
        for (size_t Q = 0; Q < naux_; Q++) {
            for (size_t onu = omu + 1; onu < nbf_; onu++) {
                if (schwarz_fun_index_[omu * nbf_ + onu]) {
                    size_t ind1 = big_skips_[onu] + Q * small_skips_[onu] + schwarz_fun_index_[onu * nbf_ + omu] - 1;
                    size_t ind2 = big_skips_[omu] + Q * small_skips_[omu] + schwarz_fun_index_[omu * nbf_ + onu] - 1;
                    Ppq[ind1] = Ppq[ind2];
                }
            }
        }
    }
}
void DFHelper::copy_upper_lower_wAO_core_symm(double* Qpq, double* Ppq, size_t begin, size_t end) {
    // copy out of symm
    size_t startind = symm_big_skips_[begin];
    for (size_t j = begin; j <= end; j++) {
        size_t mi = symm_small_skips_[j];
        size_t si = small_skips_[j];
        size_t jump = symm_ignored_columns_[j];
        size_t skip1 = big_skips_[j];
        size_t skip2 = symm_big_skips_[j] - startind;
        for (size_t j2 = 0; j2 < naux_; j2++) {
            C_DCOPY(mi, &Qpq[skip2 + mi * j2], 1, &Ppq[skip1 + jump + si * j2], 1);
        }
    }
// copy upper-to-lower
#pragma omp parallel for num_threads(nthreads_) schedule(static)
    for (size_t omu = begin; omu <= end; omu++) {
        for (size_t Q = 0; Q < naux_; Q++) {
            for (size_t onu = omu + 1; onu < nbf_; onu++) {
                if (schwarz_fun_index_[omu * nbf_ + onu]) {
                    size_t ind1 = big_skips_[onu] + Q * small_skips_[onu] + schwarz_fun_index_[onu * nbf_ + omu] - 1;
                    size_t ind2 = big_skips_[omu] + Q * small_skips_[omu] + schwarz_fun_index_[omu * nbf_ + onu] - 1;
                    Ppq[ind1] = Ppq[ind2];
                }
            }
        }
    }
}

void DFHelper::add_space(std::string key, SharedMatrix M) {
    size_t a0 = M->rowspi()[0];
    size_t a1 = M->colspi()[0];

    if (!built_) {
        throw PSIEXCEPTION("DFHelper:add_space: call initialize() before adding spaces!");
    } else if (a0 != nbf_) {
        std::stringstream error;
        error << "DFHelper:add_space: illegal space (" << key << "), primary axis is not nbf_";
        throw PSIEXCEPTION(error.str().c_str());
    } else if (spaces_.find(key) != spaces_.end()) {
        if (a1 != std::get<1>(spaces_[key])) {
            std::stringstream error;
            error << "DFHelper:add_space: illegal space (" << key << "), new space has incorrect dimension!";
            throw PSIEXCEPTION(error.str().c_str());
        }
    }
    sorted_spaces_.push_back(std::make_pair(key, a1));
    spaces_[key] = std::make_tuple(M, a1);
}
void DFHelper::add_transformation(std::string name, std::string key1, std::string key2, std::string order) {
    if (spaces_.find(key1) == spaces_.end()) {
        std::stringstream error;
        error << "DFHelper:add_transformation: first space (" << key1 << "), is not in space list!";
        throw PSIEXCEPTION(error.str().c_str());
    } else if (spaces_.find(key2) == spaces_.end()) {
        std::stringstream error;
        error << "DFHelper:add_transformation: second space (" << key2 << "), is not in space list!";
        throw PSIEXCEPTION(error.str().c_str());
    }

    int op;
    if (!order.compare("Qpq")) {
        op = 0;
    } else if (!order.compare("pQq")) {
        op = 1;
    } else if (!order.compare("pqQ")) {
        op = 2;
    } else {
        throw PSIEXCEPTION("DF_Hepler:add_transformation: incorrect integral format, use 'Qpq', 'pQq', or 'pqQ'");
    }
    transf_[name] = std::make_tuple(key1, key2, op);

    size_t a1 = std::get<1>(spaces_[key1]);
    size_t a2 = std::get<1>(spaces_[key2]);
    filename_maker(name, naux_, a1, a2, op);
}
void DFHelper::clear_spaces() {
    // clear spaces
    spaces_.clear();
    sorted_spaces_.clear();
    order_.clear();
    bspace_.clear();
    strides_.clear();

    // no ordering
    ordered_ = false;
    transformed_ = false;
}

void DFHelper::clear_transformations() {
    // clears transformations

    transf_.clear();
    transf_core_.clear();
}

void DFHelper::clear_all() {
    // invokes destructors, eliminating all files.
    file_streams_.clear();

    // clears all info
    clear_spaces();
    files_.clear();
    sizes_.clear();
    tsizes_.clear();
    clear_transformations();
}

std::pair<size_t, size_t> DFHelper::identify_order() {
    // Identify order of transformations to use strategic intermediates
    std::sort(sorted_spaces_.begin(), sorted_spaces_.end(),
              [](const std::pair<std::string, size_t>& left, const std::pair<std::string, size_t>& right) {
                  return left.second < right.second;
              });

    // copy transf_ keys into a list of needs
    std::list<std::string> needs;
    for (auto const& itr : transf_) needs.push_back(itr.first);

    // construct best transformation order
    size_t largest = 0, maximum = 0, small, large, op;
    for (size_t i = 0; i < sorted_spaces_.size(); i++) {
        bool on = false;
        size_t st = 0;
        std::string str = sorted_spaces_[i].first;

        auto itr = needs.begin();
        while (itr != needs.end()) {
            op = 0;
            op = (!(std::get<0>(transf_[*itr]).compare(str)) ? 1 : op);
            op = (!(std::get<1>(transf_[*itr]).compare(str)) ? 2 : op);
            if (op != 0) {
                if (!on) {
                    bspace_.push_back(str);
                    on = true;
                }
                small = (op == 1 ? std::get<1>(spaces_[std::get<0>(transf_[*itr])])
                                 : std::get<1>(spaces_[std::get<1>(transf_[*itr])]));
                large = (op == 1 ? std::get<1>(spaces_[std::get<1>(transf_[*itr])])
                                 : std::get<1>(spaces_[std::get<0>(transf_[*itr])]));
                maximum = (maximum < small * large ? small * large : maximum);
                largest = (largest < small ? small : largest);
                order_.push_back(*itr);
                st++;
                itr = needs.erase(itr);
            } else
                itr++;
        }
        if (st > 0) {
            strides_.push_back(st);
        }
    }
    // print_order();
    ordered_ = true;
    return std::make_pair(largest, maximum);
}
void DFHelper::print_order() {
    size_t o = order_.size();
    size_t b = bspace_.size();
    outfile->Printf("\n     ==> DFHelper:--Begin Transformations Information <==\n\n");
    outfile->Printf("   Transformation order:\n");
    for (size_t i = 0; i < o; i++) {
        outfile->Printf("         %s: (%s, %s)\n", order_[i].c_str(), std::get<0>(transf_[order_[i]]).c_str(),
                        std::get<1>(transf_[order_[i]]).c_str());
    }
    outfile->Printf("\n    Best Spaces:\n");
    for (size_t i = 0; i < b; i++) {
        outfile->Printf("         (space: %s, size: %zu)\n", bspace_[i].c_str(), std::get<1>(spaces_[bspace_[i]]));
    }
    outfile->Printf("\n    Transformation strides: ");
    for (size_t i = 0; i < b; i++) {
        outfile->Printf("%zu", strides_[i]);
        if (i < b - 1) outfile->Printf(", ");
    }
    outfile->Printf("\n\n     ==> DFHelper:--End Transformations Information <==\n\n");
}

void DFHelper::transform() {
    if (debug_) {
        outfile->Printf("Entering DFHelper::transform\n");
    }

    timer_on("DFH: transform()");
    // outfile->Printf("\n     ==> DFHelper:--Begin Transformations <==\n\n");

    size_t nthreads = nthreads_;

    // reset tranposes (in case the transpose() function was called)
    tsizes_.clear();

    // get optimal path and info
    if (!ordered_) info_ = identify_order();
    size_t wtmp = std::get<0>(info_);
    size_t wfinal = std::get<1>(info_);

    // prep AO file stream if STORE + !AO_core_
    if (!direct_iaQ_ && !direct_ && !AO_core_) stream_check(AO_names_[1], "rb");

    // get Q blocking scheme
    std::vector<std::pair<size_t, size_t>> Qsteps;
    std::pair<size_t, size_t> Qlargest = Qshell_blocks_for_transform(memory_, wtmp, wfinal, Qsteps);
    size_t max_block = std::get<1>(Qlargest);

    // prepare eri and C buffers per thread
    size_t nthread = nthreads_;
    std::vector<std::vector<double>> C_buffers(nthreads_);
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    auto rifactory = std::make_shared<IntegralFactory>(aux_, zero, primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthread);
    eri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
#pragma omp parallel num_threads(nthreads_)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        std::vector<double> Cp(nbf_ * wtmp);
        C_buffers[rank] = Cp;
        if (rank) eri[rank] = std::shared_ptr<TwoBodyAOInt>(eri.front()->clone());
    }

    // allocate in-core transformed integrals if necessary
    if (MO_core_) {
        for (auto& kv : transf_) {
            size_t size = std::get<1>(spaces_[std::get<0>(kv.second)]) * std::get<1>(spaces_[std::get<1>(kv.second)]);
            transf_core_[kv.first] = std::unique_ptr<double[]>(new double[size * naux_]);
        }
    }

    // scope buffer declarations
    {
        // declare buffers
        std::unique_ptr<double[]> T(new double[max_block * nbf_ * wtmp]);
        std::unique_ptr<double[]> F(new double[max_block * wfinal]);
        std::unique_ptr<double[]> N;
        double* Tp = T.get();
        double* Fp = F.get();
        double* Np;
        if (!MO_core_) {
            N = std::unique_ptr<double[]>(new double[max_block * wfinal]);
            Np = N.get();
        }

        // AO buffer, allocate if not in-core, else point to in-core
        std::unique_ptr<double[]> M;
        double* Mp;
        if (!AO_core_) {
            M = std::unique_ptr<double[]>(new double[std::get<0>(Qlargest)]);
            Mp = M.get();
        } else {
            Mp = Ppq_.get();
        }

        // transform in steps, blocking over the auxiliary basis (Q blocks)
        for (size_t j = 0, bcount = 0, block_size; j < Qsteps.size(); j++, bcount += block_size) {
            // Qshell step info
            size_t start = std::get<0>(Qsteps[j]);
            size_t stop = std::get<1>(Qsteps[j]);
            size_t begin = Qshell_aggs_[start];
            size_t end = Qshell_aggs_[stop + 1] - 1;
            block_size = end - begin + 1;

            // print step info
            // outfile->Printf("      Qshell: (%zu, %zu)", start, stop);
            // outfile->Printf(", PHI: (%zu, %zu), size: %zu\n", begin, end, block_size);

            // get AO chunk according to directives
            if (AO_core_) {
                ;  // pass
            } else if (direct_iaQ_) {
                timer_on("DFH: Total Workflow");
                compute_dense_Qpq_blocking_Q(start, stop, Mp, eri);
                timer_off("DFH: Total Workflow");
            } else if (direct_) {
                timer_on("DFH: Total Workflow");
                compute_sparse_pQq_blocking_Q(start, stop, Mp, eri);
                timer_off("DFH: Total Workflow");
            } else {
                timer_on("DFH: Grabbing AOs");
                grab_AO(start, stop, Mp);
                timer_off("DFH: Grabbing AOs");
            }

            // stride through best spaces
            for (size_t i = 0, count = 0; i < bspace_.size(); count += strides_[i], i++) {
                // grab best space
                std::string bspace = bspace_[i];
                size_t bsize = std::get<1>(spaces_[bspace]);
                double* Bp = std::get<0>(spaces_[bspace])->pointer()[0];

                // index bump if AOs are in core
                size_t bump = (AO_core_ ? bcount * nbf_ * nbf_ : 0);

                // perform first contraction
                timer_on("DFH: Total Workflow");
                timer_on("DFH: Total Transform");
                timer_on("DFH: 1st Contraction");
                if (direct_iaQ_) {
                    // (qb)(Q|pq)->(Q|pb)
                    C_DGEMM('N', 'N', block_size * nbf_, bsize, nbf_, 1.0, &Mp[bump], nbf_, Bp, bsize, 0.0, Tp, bsize);
                } else {
                    // (bq)(p|Qq)->(p|Qb)
                    first_transform_pQq(bsize, bcount, block_size, Mp, Tp, Bp, C_buffers);
                }
                timer_off("DFH: 1st Contraction");
                timer_off("DFH: Total Transform");
                timer_off("DFH: Total Workflow");

                // to completion per transformation
                for (size_t k = 0; k < strides_[i]; k++) {
                    // get transformation info
                    auto transf_name = order_[count + k];
                    auto left = std::get<0>(transf_[transf_name]);
                    auto right = std::get<1>(transf_[transf_name]);
                    bool bleft = (bspace.compare(left) == 0 ? true : false);

                    // get worst space
                    std::tuple<SharedMatrix, size_t> I = (bleft ? spaces_[right] : spaces_[left]);
                    double* Wp = std::get<0>(I)->pointer()[0];
                    size_t wsize = std::get<1>(I);

                    // grab in-core pointer
                    if (direct_iaQ_ && MO_core_) {
                        Fp = transf_core_[order_[count + k]].get();
                    } else if (MO_core_) {
                        Np = transf_core_[order_[count + k]].get();
                    }

                    // perform final contraction
                    // (wp)(p|Qb)->(w|Qb)
                    timer_on("DFH: Total Workflow");
                    timer_on("DFH: Total Transform");
                    timer_on("DFH: 2nd Contraction");
                    if (direct_iaQ_) {
                        size_t bump = (MO_core_ ? begin * wsize * bsize : 0);
                        // (pw)(Q|pb)->(Q|bw)
                        if (bleft) {
#pragma omp parallel for num_threads(nthreads_)
                            for (size_t i = 0; i < block_size; i++) {
                                C_DGEMM('T', 'N', bsize, wsize, nbf_, 1.0, &Tp[i * nbf_ * bsize], bsize, Wp, wsize, 0.0,
                                        &Fp[bump + i * wsize * bsize], wsize);
                            }
                        } else {
// (pw)(Q|pb)->(Q|wb)
#pragma omp parallel for num_threads(nthreads_)
                            for (size_t i = 0; i < block_size; i++) {
                                C_DGEMM('T', 'N', wsize, bsize, nbf_, 1.0, Wp, wsize, &Tp[i * nbf_ * bsize], bsize, 0.0,
                                        &Fp[bump + i * wsize * bsize], bsize);
                            }
                        }
                    } else {
                        // (pw)(p|Qb)->(w|Qb)
                        C_DGEMM('T', 'N', wsize, block_size * bsize, nbf_, 1.0, Wp, wsize, Tp, block_size * bsize, 0.0,
                                Fp, block_size * bsize);
                    }
                    timer_off("DFH: 2nd Contraction");
                    timer_off("DFH: Total Transform");
                    timer_off("DFH: Total Workflow");

                    // put the transformations away
                    timer_on("DFH: MO to disk, " + transf_name);
                    if (direct_iaQ_) {
                        put_transformations_Qpq(begin, end, wsize, bsize, Fp, count + k, bleft);
                    } else {
                        put_transformations_pQq(begin, end, block_size, bcount, wsize, bsize, Np, Fp, count + k, bleft);
                    }
                    timer_off("DFH: MO to disk, " + transf_name);
                }
            }
        }
    }  // buffers destroyed with std housekeeping

    // outfile->Printf("\n     ==> DFHelper:--End Transformations (disk)<==\n\n");

    // transformations complete, time for metric contractions

    timer_on("DFH: Metric Contractions");

    // release in-core AO 
    if (AO_core_ && release_core_AO_before_metric_)  {
        Ppq_.reset();
    }

    if (direct_iaQ_ || direct_) {
        // prepare metric
        std::unique_ptr<double[]> metric;
        double* metp;
        if (!hold_met_) {
            metric = std::unique_ptr<double[]>(new double[naux_ * naux_]);
            metp = metric.get();
            std::string filename = return_metfile(mpower_);
            get_tensor_(std::get<0>(files_[filename]), metp, 0, naux_ - 1, 0, naux_ - 1);
        } else
            metp = metric_prep_core(mpower_);

        if (direct_iaQ_) {
            if (MO_core_) {
                std::unique_ptr<double[]> N(new double[naux_ * wfinal]);
                double* Np = N.get();

                for (auto& kv : transf_core_) {
                    size_t l = std::get<0>(sizes_[std::get<1>(files_[kv.first])]);
                    size_t r = std::get<1>(sizes_[std::get<1>(files_[kv.first])]);
                    size_t Q = std::get<2>(sizes_[std::get<1>(files_[kv.first])]);

                    double* Lp = kv.second.get();
                    C_DCOPY(l * r * Q, Lp, 1, Np, 1);

                    // (Q|ia) (PQ) -> (ia|Q)
                    C_DGEMM('T', 'N', l * r, Q, Q, 1.0, Np, l * r, metp, Q, 0.0, Lp, Q);
                }

            } else {
                // total size allowed, in doubles
                size_t rem_mem = memory_ - (AO_core_ && !release_core_AO_before_metric_ ? naux_ * nbf_ * nbf_ : 0);
                size_t total_mem =
                    (rem_mem > wfinal * naux_ * 2 + naux_ * naux_ ? wfinal * naux_ : (rem_mem - naux_ * naux_) / 2);

                std::unique_ptr<double[]> M(new double[total_mem]);
                std::unique_ptr<double[]> F(new double[total_mem]);
                double* Mp = M.get();
                double* Fp = F.get();
                for (std::vector<std::string>::iterator itr = order_.begin(); itr != order_.end(); itr++)
                    contract_metric_Qpq(*itr, metp, Mp, Fp, total_mem);
            }

        } else if (direct_) {
            if (!MO_core_) {
                // total size allowed, in doubles.
                // note that memory - naux_2 cannot be negative (handled in init)
                size_t rem_mem = memory_ - (AO_core_ && !release_core_AO_before_metric_ ? big_skips_[nbf_] : 0);
                size_t total_mem =
                    (rem_mem > wfinal * naux_ * 2 + naux_ * naux_ ? wfinal * naux_ : (rem_mem - naux_ * naux_) / 2);

                std::unique_ptr<double[]> M(new double[total_mem]);
                std::unique_ptr<double[]> F(new double[total_mem]);
                double* Mp = M.get();
                double* Fp = F.get();

                for (std::vector<std::string>::iterator itr = order_.begin(); itr != order_.end(); itr++)
                    contract_metric(*itr, metp, Mp, Fp, total_mem);

            } else {
                std::unique_ptr<double[]> N(new double[naux_ * wfinal]);
                double* Np = N.get();

                for (auto& kv : transf_core_) {
                    size_t a0 = std::get<0>(sizes_[std::get<1>(files_[kv.first])]);
                    size_t a1 = std::get<1>(sizes_[std::get<1>(files_[kv.first])]);
                    size_t a2 = std::get<2>(sizes_[std::get<1>(files_[kv.first])]);

                    double* Lp = kv.second.get();
                    C_DCOPY(a0 * a1 * a2, Lp, 1, Np, 1);

                    // the following differs depending on the form being outputted
                    // 0 : Qpq -- 1 : pQq -- 2 : pqQ

                    if (std::get<2>(transf_[kv.first]) == 2) {
                        C_DGEMM('N', 'N', a0 * a1, a2, a2, 1.0, Np, a2, metp, a2, 0.0, Lp, a2);
                    } else if (std::get<2>(transf_[kv.first]) == 0) {
                        C_DGEMM('N', 'N', a0, a1 * a2, a0, 1.0, metp, naux_, Np, a1 * a2, 0.0, Lp, a1 * a2);
                    } else {
#pragma omp parallel for num_threads(nthreads_)
                        for (size_t i = 0; i < a0; i++) {
                            C_DGEMM('N', 'N', a1, a2, a1, 1.0, metp, naux_, &Np[i * a1 * a2], a2, 0.0, &Lp[i * a1 * a2],
                                    a2);
                        }
                    }
                }
            }
        }
    }
    timer_off("DFH: Metric Contractions");
    timer_off("DFH: transform()");
    transformed_ = true;

    if (debug_) {
        outfile->Printf("Exiting DFHelper::transform\n");
    }
}

void DFHelper::first_transform_pQq(size_t bsize, size_t bcount, size_t block_size, double* Mp, double* Tp, double* Bp,
                                   std::vector<std::vector<double>>& C_buffers) {
// perform first contraction on pQq, thread over p.
#pragma omp parallel for schedule(guided) num_threads(nthreads_)
    for (size_t k = 0; k < nbf_; k++) {
        // truncate transformation matrix according to fun_mask
        size_t sp_size = small_skips_[k];
        size_t jump = (AO_core_ ? big_skips_[k] + bcount * sp_size : (big_skips_[k] * block_size) / naux_);

        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        for (size_t m = 0, sp_count = -1; m < nbf_; m++) {
            if (schwarz_fun_index_[k * nbf_ + m]) {
                sp_count++;
                C_DCOPY(bsize, &Bp[m * bsize], 1, &C_buffers[rank][sp_count * bsize], 1);
            }
        }

        // (Qm)(mb)->(Qb)
        C_DGEMM('N', 'N', block_size, bsize, sp_size, 1.0, &Mp[jump], sp_size, &C_buffers[rank][0], bsize, 0.0,
                &Tp[k * block_size * bsize], bsize);
    }
}

void DFHelper::put_transformations_Qpq(int begin, int end, int wsize, int bsize, double* Fp, int ind, bool bleft) {
    // incoming transformed integrals to this function are in a Qpq format.
    // if MO_core is on, do nothing
    // else, the buffers are put to disk.

    if (!MO_core_) {
        // "ab" is great here since we are actually appending.
        std::string putf = std::get<0>(files_[order_[ind]]);
        std::string op = "ab";

        if (bleft) {
            put_tensor(putf, Fp, std::make_pair(begin, end), std::make_pair(0, bsize - 1), std::make_pair(0, wsize - 1),
                       op);
        } else {
            put_tensor(putf, Fp, std::make_pair(begin, end), std::make_pair(0, wsize - 1), std::make_pair(0, bsize - 1),
                       op);
        }
    }
}

void DFHelper::put_transformations_pQq(int begin, int end, int rblock_size, int bcount, int wsize, int bsize,
                                       double* Np, double* Fp, int ind, bool bleft) {
    // incoming transformed integrals to this function are in a pQq format.
    // first, the integrals are tranposed to the desired format specified in add_transformation().
    // if MO_core is on, then the LHS buffers are final destinations.
    // else, the buffers are put to disk.

    // setup ~
    int lblock_size = rblock_size;
    std::string putf, op;
    if (!MO_core_) {
        putf = (!direct_ ? std::get<1>(files_[order_[ind]]) : std::get<0>(files_[order_[ind]]));
        op = "wb";
        bcount = 0;
    } else {
        lblock_size = naux_;
    }

    if (bleft) {
        // result is in pqQ format
        if (std::get<2>(transf_[order_[ind]]) == 2) {
// (w|Qb)->(bw|Q)
#pragma omp parallel for num_threads(nthreads_)
            for (size_t z = 0; z < wsize; z++) {
                for (size_t y = 0; y < bsize; y++) {
                    for (size_t x = 0; x < rblock_size; x++) {
                        Np[y * wsize * lblock_size + z * lblock_size + (bcount + x)] =
                            Fp[z * bsize * rblock_size + x * bsize + y];
                    }
                }
            }
            if (!MO_core_) {
                put_tensor(putf, Np, std::make_pair(0, bsize - 1), std::make_pair(0, wsize - 1),
                           std::make_pair(begin, end), op);
            }

            // result is in Qpq format
        } else if (std::get<2>(transf_[order_[ind]]) == 0) {
// (w|Qb)->(Q|bw)
#pragma omp parallel for num_threads(nthreads_)
            for (size_t x = 0; x < rblock_size; x++) {
                for (size_t z = 0; z < wsize; z++) {
                    for (size_t y = 0; y < bsize; y++) {
                        Np[(bcount + x) * bsize * wsize + y * wsize + z] = Fp[z * bsize * rblock_size + x * bsize + y];
                    }
                }
            }
            if (!MO_core_) {
                put_tensor(putf, Np, std::make_pair(begin, end), std::make_pair(0, bsize - 1),
                           std::make_pair(0, wsize - 1), op);
            }

            // result is in pQq format
        } else {
// (w|Qb)->(bQw)
#pragma omp parallel for num_threads(nthreads_)
            for (size_t x = 0; x < rblock_size; x++) {
                for (size_t y = 0; y < bsize; y++) {
                    for (size_t z = 0; z < wsize; z++) {
                        Np[y * lblock_size * wsize + (bcount + x) * wsize + z] =
                            Fp[z * bsize * rblock_size + x * bsize + y];
                    }
                }
            }
            if (!MO_core_) {
                put_tensor(putf, Np, std::make_pair(0, bsize - 1), std::make_pair(begin, end),
                           std::make_pair(0, wsize - 1), op);
            }
        }

    } else {
        // result is in pqQ format
        if (std::get<2>(transf_[order_[ind]]) == 2) {
// (w|Qb)->(wbQ)
#pragma omp parallel for num_threads(nthreads_)
            for (size_t z = 0; z < wsize; z++) {
                for (size_t x = 0; x < rblock_size; x++) {
                    for (size_t y = 0; y < bsize; y++) {
                        Np[z * lblock_size * bsize + y * lblock_size + (bcount + x)] =
                            Fp[z * rblock_size * bsize + x * bsize + y];
                    }
                }
            }
            if (!MO_core_) {
                put_tensor(putf, Np, std::make_pair(0, wsize - 1), std::make_pair(0, bsize - 1),
                           std::make_pair(begin, end), op);
            }

            // result is in Qpq format
        } else if (std::get<2>(transf_[order_[ind]]) == 0) {
// (w|Qb)->(Q|wb)
#pragma omp parallel for num_threads(nthreads_)
            for (size_t x = 0; x < rblock_size; x++) {
                for (size_t z = 0; z < wsize; z++) {
                    C_DCOPY(bsize, &Fp[z * rblock_size * bsize + x * bsize], 1,
                            &Np[(bcount + x) * wsize * bsize + z * bsize], 1);
                }
            }
            if (!MO_core_) {
                put_tensor(putf, Np, std::make_pair(begin, end), std::make_pair(0, wsize - 1),
                           std::make_pair(0, bsize - 1), op);
            }

            // result is in pQq format
        } else {
            // (w|Qb)
            if (!MO_core_) {
                put_tensor(putf, Fp, std::make_pair(0, wsize - 1), std::make_pair(begin, end),
                           std::make_pair(0, bsize - 1), op);
            } else {
// we have to copy over the buffer
#pragma omp parallel for num_threads(nthreads_)
                for (size_t x = 0; x < wsize; x++) {
                    for (size_t y = 0; y < rblock_size; y++) {
                        C_DCOPY(bsize, &Fp[x * rblock_size * bsize + y * bsize], 1,
                                &Np[x * lblock_size * bsize + (bcount + y) * bsize], 1);
                    }
                }
            }
        }
    }
}

// Fill using a pointer, be cautious of bounds!!
void DFHelper::fill_tensor(std::string name, double* b) {
    check_file_key(name);
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    fill_tensor(name, b, {0, std::get<0>(sizes)}, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DFHelper::fill_tensor(std::string name, double* b, std::vector<size_t> a1) {
    check_file_key(name);
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    fill_tensor(name, b, a1, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DFHelper::fill_tensor(std::string name, double* b, std::vector<size_t> a1, std::vector<size_t> a2) {
    check_file_key(name);
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    fill_tensor(name, b, a1, a2, {0, std::get<2>(sizes)});
}
void DFHelper::fill_tensor(std::string name, double* b, std::vector<size_t> a1, std::vector<size_t> a2,
                           std::vector<size_t> a3) {
    if (a1.size() != 2) {
        std::stringstream error;
        error << "DFHelper:fill_tensor:  axis 0 tensor indexing vector has " << a1.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (a2.size() != 2) {
        std::stringstream error;
        error << "DFHelper:fill_tensor:  axis 1 tensor indexing vector has " << a2.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (a3.size() != 2) {
        std::stringstream error;
        error << "DFHelper:fill_tensor:  axis 2 tensor indexing vector has " << a3.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }

    check_file_key(name);
    std::string filename = std::get<1>(files_[name]);

    // being pythonic ;)
    std::pair<size_t, size_t> i0 = std::make_pair(a1[0], a1[1] - 1);
    std::pair<size_t, size_t> i1 = std::make_pair(a2[0], a2[1] - 1);
    std::pair<size_t, size_t> i2 = std::make_pair(a3[0], a3[1] - 1);

    get_tensor_(filename, b, i0, i1, i2);
}

// Fill using a pre-allocated SharedMatrix
void DFHelper::fill_tensor(std::string name, SharedMatrix M) {
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    fill_tensor(name, M, {0, std::get<0>(sizes)}, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DFHelper::fill_tensor(std::string name, SharedMatrix M, std::vector<size_t> a1) {
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    fill_tensor(name, M, a1, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DFHelper::fill_tensor(std::string name, SharedMatrix M, std::vector<size_t> a1, std::vector<size_t> a2) {
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    fill_tensor(name, M, a1, a2, {0, std::get<2>(sizes)});
}
void DFHelper::fill_tensor(std::string name, SharedMatrix M, std::vector<size_t> t0, std::vector<size_t> t1,
                           std::vector<size_t> t2) {
    std::string filename = std::get<1>(files_[name]);
    // has this integral been transposed?
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    if (t0.size() != 2) {
        std::stringstream error;
        error << "DFHelper:fill_tensor:  axis 0 tensor indexing vector has " << t0.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (t1.size() != 2) {
        std::stringstream error;
        error << "DFHelper:fill_tensor:  axis 1 tensor indexing vector has " << t1.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (t2.size() != 2) {
        std::stringstream error;
        error << "DFHelper:fill_tensor:  axis 2 tensor indexing vector has " << t2.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }

    // be pythonic - adjust stops
    size_t sta0 = t0[0];
    size_t sto0 = t0[1] - 1;
    size_t sta1 = t1[0];
    size_t sto1 = t1[1] - 1;
    size_t sta2 = t2[0];
    size_t sto2 = t2[1] - 1;

    std::pair<size_t, size_t> i0 = std::make_pair(sta0, sto0);
    std::pair<size_t, size_t> i1 = std::make_pair(sta1, sto1);
    std::pair<size_t, size_t> i2 = std::make_pair(sta2, sto2);

    check_file_key(name);
    check_file_tuple(name, i0, i1, i2);
    check_matrix_size(name, M, i0, i1, i2);

    size_t A0 = (sto0 - sta0 + 1);
    size_t A1 = (sto1 - sta1 + 1);
    size_t A2 = (sto2 - sta2 + 1);

    double* Mp = M->pointer()[0];
    if (MO_core_) {
        size_t a0 = std::get<0>(sizes);
        size_t a1 = std::get<1>(sizes);
        size_t a2 = std::get<2>(sizes);

        double* Fp = transf_core_[name].get();
#pragma omp parallel for num_threads(nthreads_)
        for (size_t i = 0; i < A0; i++) {
            for (size_t j = 0; j < A1; j++) {
#pragma omp simd
                for (size_t k = 0; k < A2; k++) {
                    Mp[i * A1 * A2 + j * A2 + k] = Fp[(sta0 + i) * a1 * a2 + (sta1 + j) * a2 + (sta2 + k)];
                }
            }
        }
    } else {
        get_tensor_(filename, Mp, i0, i1, i2);
    }
    M->set_numpy_shape({(int)A0, (int)A1, (int)A2});
}

// Return a SharedMatrix
SharedMatrix DFHelper::get_tensor(std::string name) {
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    return get_tensor(name, {0, std::get<0>(sizes)}, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
SharedMatrix DFHelper::get_tensor(std::string name, std::vector<size_t> a1) {
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    return get_tensor(name, a1, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
SharedMatrix DFHelper::get_tensor(std::string name, std::vector<size_t> a1, std::vector<size_t> a2) {
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    return get_tensor(name, a1, a2, {0, std::get<2>(sizes)});
}
SharedMatrix DFHelper::get_tensor(std::string name, std::vector<size_t> t0, std::vector<size_t> t1,
                                  std::vector<size_t> t2) {
    // has this integral been transposed?
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    if (t0.size() != 2) {
        std::stringstream error;
        error << "DFHelper:fill_tensor:  axis 0 tensor indexing vector has " << t0.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (t1.size() != 2) {
        std::stringstream error;
        error << "DFHelper:fill_tensor:  axis 1 tensor indexing vector has " << t1.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (t2.size() != 2) {
        std::stringstream error;
        error << "DFHelper:fill_tensor:  axis 2 tensor indexing vector has " << t2.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }

    // be pythonic - adjust stops
    size_t sta0 = t0[0];
    size_t sto0 = t0[1] - 1;
    size_t sta1 = t1[0];
    size_t sto1 = t1[1] - 1;
    size_t sta2 = t2[0];
    size_t sto2 = t2[1] - 1;

    std::pair<size_t, size_t> i0 = std::make_pair(sta0, sto0);
    std::pair<size_t, size_t> i1 = std::make_pair(sta1, sto1);
    std::pair<size_t, size_t> i2 = std::make_pair(sta2, sto2);

    check_file_key(name);
    check_file_tuple(name, i0, i1, i2);

    size_t A0 = (sto0 - sta0 + 1);
    size_t A1 = (sto1 - sta1 + 1);
    size_t A2 = (sto2 - sta2 + 1);

    auto M = std::make_shared<Matrix>("M", A0, A1 * A2);
    double* Mp = M->pointer()[0];

    if (MO_core_) {
        size_t a0 = std::get<0>(sizes);
        size_t a1 = std::get<1>(sizes);
        size_t a2 = std::get<2>(sizes);

        double* Fp = transf_core_[name].get();
#pragma omp parallel for num_threads(nthreads_)
        for (size_t i = 0; i < A0; i++) {
            for (size_t j = 0; j < A1; j++) {
#pragma omp simd
                for (size_t k = 0; k < A2; k++) {
                    Mp[i * A1 * A2 + j * A2 + k] = Fp[(sta0 + i) * a1 * a2 + (sta1 + j) * a2 + (sta2 + k)];
                }
            }
        }
    } else {
        get_tensor_(filename, Mp, i0, i1, i2);
    }
    M->set_numpy_shape({(int)A0, (int)A1, (int)A2});
    return M;
}

// Add a disk tensor
void DFHelper::add_disk_tensor(std::string key, std::tuple<size_t, size_t, size_t> dimensions) {
    if (files_.count(key)) {
        std::stringstream error;
        error << "DFHelper:add_disk_tensor:  tensor already exists: (" << key << "!";
        throw PSIEXCEPTION(error.str().c_str());
    }

    filename_maker(key, std::get<0>(dimensions), std::get<1>(dimensions), std::get<2>(dimensions));
}

// Write to a disk tensor from Sharedmatrix
void DFHelper::write_disk_tensor(std::string key, SharedMatrix M) {
    check_file_key(key);
    std::string filename = std::get<1>(files_[key]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);
    write_disk_tensor(key, M, {0, std::get<0>(sizes)}, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DFHelper::write_disk_tensor(std::string key, SharedMatrix M, std::vector<size_t> a1) {
    check_file_key(key);
    std::string filename = std::get<1>(files_[key]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);
    write_disk_tensor(key, M, a1, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DFHelper::write_disk_tensor(std::string key, SharedMatrix M, std::vector<size_t> a1, std::vector<size_t> a2) {
    check_file_key(key);
    std::string filename = std::get<1>(files_[key]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);
    write_disk_tensor(key, M, a1, a2, {0, std::get<2>(sizes)});
}
void DFHelper::write_disk_tensor(std::string key, SharedMatrix M, std::vector<size_t> a0, std::vector<size_t> a1,
                                 std::vector<size_t> a2) {
    // being pythonic ;)
    std::pair<size_t, size_t> i0 = std::make_pair(a0[0], a0[1] - 1);
    std::pair<size_t, size_t> i1 = std::make_pair(a1[0], a1[1] - 1);
    std::pair<size_t, size_t> i2 = std::make_pair(a2[0], a2[1] - 1);

    check_file_key(key);
    check_file_tuple(key, i0, i1, i2);
    check_matrix_size(key, M, i0, i1, i2);

    // "wb" is the way to go. the stream will change when when the tensor is read,
    // but this should be extendible to back-and-forth read/writes.
    std::string op = "wb";
    put_tensor(std::get<1>(files_[key]), M->pointer()[0], i0, i1, i2, op);
}

// Write to a disk tensor from pointer, be careful!
void DFHelper::write_disk_tensor(std::string key, double* b) {
    check_file_key(key);
    std::string filename = std::get<1>(files_[key]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);
    write_disk_tensor(key, b, {0, std::get<0>(sizes)}, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DFHelper::write_disk_tensor(std::string key, double* b, std::vector<size_t> a0) {
    check_file_key(key);
    std::string filename = std::get<1>(files_[key]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);
    write_disk_tensor(key, b, a0, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DFHelper::write_disk_tensor(std::string key, double* b, std::vector<size_t> a0, std::vector<size_t> a1) {
    check_file_key(key);
    std::string filename = std::get<1>(files_[key]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);
    write_disk_tensor(key, b, a0, a1, {0, std::get<2>(sizes)});
}
void DFHelper::write_disk_tensor(std::string key, double* b, std::vector<size_t> a0, std::vector<size_t> a1,
                                 std::vector<size_t> a2) {
    // being pythonic ;)
    std::pair<size_t, size_t> i0 = std::make_pair(a0[0], a0[1] - 1);
    std::pair<size_t, size_t> i1 = std::make_pair(a1[0], a1[1] - 1);
    std::pair<size_t, size_t> i2 = std::make_pair(a2[0], a2[1] - 1);

    check_file_key(key);
    check_file_tuple(key, i0, i1, i2);

    // "wb" is the way to go. the stream will change when when the tensor is read,
    // but this should be extendible to back-and-forth read/writes.
    std::string op = "wb";
    put_tensor(std::get<1>(files_[key]), b, i0, i1, i2, op);
}

void DFHelper::check_file_key(std::string name) {
    if (files_.find(name) == files_.end()) {
        std::stringstream error;
        error << "DFHelper:get_tensor OR write_tensor: " << name << " not found.";
        throw PSIEXCEPTION(error.str().c_str());
    }
}
void DFHelper::check_matrix_size(std::string name, SharedMatrix M, std::pair<size_t, size_t> t0,
                                 std::pair<size_t, size_t> t1, std::pair<size_t, size_t> t2) {
    size_t A0 = std::get<1>(t0) - std::get<0>(t0) + 1;
    size_t A1 = (std::get<1>(t1) - std::get<0>(t1) + 1) * (std::get<1>(t2) - std::get<0>(t2) + 1);

    size_t a0 = M->rowspi()[0];
    size_t a1 = M->colspi()[0];

    if (A0 * A1 > a0 * a1) {
        std::stringstream error;
        error << "DFHelper:get_tensor: your matrix contridicts your tuple sizes when obtaining the (" << name
              << ") integral.  ";
        error << "you gave me a matrix of size: (" << a0 << "," << a1 << "), but tuple sizes give:(" << A0 << "," << A1
              << ")";
        throw PSIEXCEPTION(error.str().c_str());
    }
}
void DFHelper::check_file_tuple(std::string name, std::pair<size_t, size_t> t0, std::pair<size_t, size_t> t1,
                                std::pair<size_t, size_t> t2) {
    size_t sta0 = std::get<0>(t0);
    size_t sto0 = std::get<1>(t0);
    size_t sta1 = std::get<0>(t1);
    size_t sto1 = std::get<1>(t1);
    size_t sta2 = std::get<0>(t2);
    size_t sto2 = std::get<1>(t2);
    std::string filename = std::get<1>(files_[name]);

    // has this integral been transposed?
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    if (sta0 > sto0) {
        std::stringstream error;
        error << "when getting integral: (" << name << ")"
              << " your axis 0 tuple has a larger start index: " << sta0 << " than its stop index: " << sto0;
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (sta1 > sto1) {
        std::stringstream error;
        error << "when getting integral: (" << name << ")"
              << " your axis 1 tuple has a larger start index: " << sta1 << " than its stop index: " << sto1;
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (sta2 > sto2) {
        std::stringstream error;
        error << "when getting integral: (" << name << ")"
              << " your axis 2 tuple has a larger start index: " << sta2 << " than its stop index: " << sto2;
        throw PSIEXCEPTION(error.str().c_str());
    }
    size_t M0 = std::get<0>(sizes);
    if (sto0 > M0 - 1) {
        std::stringstream error;
        error << "your axis 0 tuple goes out of bounds when getting integral: " << name;
        error << ". you entered (" << sto0 << "), but bounds is (" << M0 - 1 << ").";
        throw PSIEXCEPTION(error.str().c_str());
    }
    size_t M1 = std::get<1>(sizes);
    if (sto1 > M1 - 1) {
        std::stringstream error;
        error << "your axis 1 tuple goes out of bounds when getting integral: " << name;
        error << ". you entered (" << sto1 << "), but bounds is (" << M1 - 1 << ").";
        throw PSIEXCEPTION(error.str().c_str());
    }
    size_t M2 = std::get<2>(sizes);
    if (sto2 > M2 - 1) {
        std::stringstream error;
        error << "your axis 2 tuple goes out of bounds when getting integral: " << name;
        error << ". you entered (" << sto2 << "), but bounds is (" << M2 - 1 << ").";
        throw PSIEXCEPTION(error.str().c_str());
    }
}
void DFHelper::transpose(std::string name, std::tuple<size_t, size_t, size_t> order) {
    if (!files_.count(name)) {
        std::stringstream error;
        error << "DFHelper::transpose(): cannot transpose input (" << name << "), name doe not exist!";
        throw PSIEXCEPTION(error.str().c_str());
    }

    (MO_core_ ? transpose_core(name, order) : transpose_disk(name, order));
}
void DFHelper::transpose_core(std::string name, std::tuple<size_t, size_t, size_t> order) {
    size_t a0 = std::get<0>(order);
    size_t a1 = std::get<1>(order);
    size_t a2 = std::get<2>(order);

    std::string filename = std::get<1>(files_[name]);
    size_t M0 = std::get<0>(sizes_[filename]);
    size_t M1 = std::get<1>(sizes_[filename]);
    size_t M2 = std::get<2>(sizes_[filename]);
    std::tuple<size_t, size_t, size_t> sizes;

    std::unique_ptr<double[]> M(new double[M0 * M1 * M2]);
    double* Mp = M.get();
    double* Fp = transf_core_[name].get();
    C_DCOPY(M0 * M1 * M2, Fp, 1, Mp, 1);

    bool on = false;
    if (a0 == 0) {
        if (a1 == 2) {
            sizes = std::make_tuple(M0, M2, M1);
            on = true;
        }
    } else if (a0 == 1) {
        if (a1 == 0) {
            sizes = std::make_tuple(M1, M0, M2);
            on = true;
        } else if (a1 == 2) {
            sizes = std::make_tuple(M1, M2, M0);
            on = true;
        }
    } else {
        if (a1 == 0) {
            sizes = std::make_tuple(M2, M0, M1);
            on = true;
        } else if (a1 == 1) {
            sizes = std::make_tuple(M2, M1, M0);
            on = true;
        }
    }

    if (!on) throw PSIEXCEPTION("you transposed all wrong!");

    if (a0 == 0) {
        if (a1 == 2) {  // (0|12) -> (0|21)
#pragma omp parallel for num_threads(nthreads_)
            for (size_t i = 0; i < M0; i++) {
                for (size_t j = 0; j < M1; j++) {
                    for (size_t k = 0; k < M2; k++) {
                        Fp[i * M1 * M2 + k * M1 + j] = Mp[i * M1 * M2 + j * M2 + k];
                    }
                }
            }
        }
    } else if (a0 == 1) {
        if (a1 == 0) {  // (0|12) -> (1|02)
#pragma omp parallel for num_threads(nthreads_)
            for (size_t i = 0; i < M0; i++) {
                for (size_t j = 0; j < M1; j++) {
#pragma omp simd
                    for (size_t k = 0; k < M2; k++) {
                        Fp[j * M0 * M2 + i * M2 + k] = Mp[i * M1 * M2 + j * M2 + k];
                    }
                }
            }
        } else if (a1 == 2) {  // (0|12) -> (1|20)
#pragma omp parallel for num_threads(nthreads_)
            for (size_t i = 0; i < M0; i++) {
                for (size_t j = 0; j < M1; j++) {
                    for (size_t k = 0; k < M2; k++) {
                        Fp[j * M0 * M2 + k * M0 + i] = Mp[i * M1 * M2 + j * M2 + k];
                    }
                }
            }
        }
    } else if (a0 == 2) {
        if (a1 == 0) {  // (0|12) -> (2|01)
#pragma omp parallel for num_threads(nthreads_)
            for (size_t i = 0; i < M0; i++) {
                for (size_t j = 0; j < M1; j++) {
                    for (size_t k = 0; k < M2; k++) {
                        Fp[k * M1 * M0 + i * M1 + j] = Mp[i * M1 * M2 + j * M2 + k];
                    }
                }
            }
        } else if (a1 == 1) {  // (0|12) -> (2|10)
#pragma omp parallel for num_threads(nthreads_)
            for (size_t i = 0; i < M0; i++) {
                for (size_t j = 0; j < M1; j++) {
                    for (size_t k = 0; k < M2; k++) {
                        Fp[k * M1 * M0 + j * M0 + i] = Mp[i * M1 * M2 + j * M2 + k];
                    }
                }
            }
        }
    }
    // keep tsizes_ separate and do not ovwrt sizes_ in case of STORE directive
    tsizes_[filename] = sizes;
}
void DFHelper::transpose_disk(std::string name, std::tuple<size_t, size_t, size_t> order) {
    size_t a0 = std::get<0>(order);
    size_t a1 = std::get<1>(order);
    size_t a2 = std::get<2>(order);

    // determine blocking
    std::string filename = std::get<1>(files_[name]);
    size_t M0 = std::get<0>(sizes_[filename]);
    size_t M1 = std::get<1>(sizes_[filename]);
    size_t M2 = std::get<2>(sizes_[filename]);

    size_t current = 0, count = 0, largest = 0;
    std::vector<std::pair<size_t, size_t>> steps;
    for (size_t i = 0; i < M0; i++) {
        current += M1 * M2;
        count++;
        if ((current * 2 > memory_) || (i == M0 - 1)) {  //
            if (count == 1 && i != M0 - 1) {
                std::stringstream error;
                error << "DFHelper:transpose_disk: not enough memory.";
                throw PSIEXCEPTION(error.str().c_str());
            }
            if (i == M0 - 1)
                steps.push_back(std::make_pair(i - count + 1, i));
            else {
                current -= M1 * M2;
                steps.push_back(std::make_pair(i - count + 1, i - 1));
                i--;
            }
            if (largest < current) largest = current;
            count = 0;
            current = 0;
        }
    }

    // declare
    std::unique_ptr<double[]> M(new double[largest]);
    std::unique_ptr<double[]> F(new double[largest]);
    double* Mp = M.get();
    double* Fp = F.get();
    std::tuple<size_t, size_t, size_t> sizes;

    bool on = false;
    if (a0 == 0) {
        if (a1 == 2) {
            sizes = std::make_tuple(M0, M2, M1);
            on = true;
        }
    } else if (a0 == 1) {
        if (a1 == 0) {
            sizes = std::make_tuple(M1, M0, M2);
            on = true;
        } else if (a1 == 2) {
            sizes = std::make_tuple(M1, M2, M0);
            on = true;
        }
    } else {
        if (a1 == 0) {
            sizes = std::make_tuple(M2, M0, M1);
            on = true;
        } else if (a1 == 1) {
            sizes = std::make_tuple(M2, M1, M0);
            on = true;
        }
    }

    if (!on) throw PSIEXCEPTION("you transposed all wrong!");

    std::string new_file = "newfilefortransposition";
    filename_maker(new_file, std::get<0>(sizes), std::get<1>(sizes), std::get<2>(sizes));
    std::string new_filename = std::get<1>(files_[new_file]);

    for (size_t m = 0; m < steps.size(); m++) {
        std::string op = (m ? "r+b" : "wb");
        size_t start = std::get<0>(steps[m]);
        size_t stop = std::get<1>(steps[m]);
        M0 = stop - start + 1;

        // grab
        get_tensor_(filename, Mp, start, stop, 0, M1 * M2 - 1);

        if (a0 == 0) {
            if (a1 == 2) {  // (0|12) -> (0|21)
#pragma omp parallel for num_threads(nthreads_)
                for (size_t i = 0; i < M0; i++) {
                    for (size_t j = 0; j < M1; j++) {
                        for (size_t k = 0; k < M2; k++) {
                            Fp[i * M1 * M2 + k * M1 + j] = Mp[i * M1 * M2 + j * M2 + k];
                        }
                    }
                }
                put_tensor(new_filename, Fp, std::make_pair(start, stop), std::make_pair(0, M2 - 1),
                           std::make_pair(0, M1 - 1), op);
            }
        } else if (a0 == 1) {
            if (a1 == 0) {  // (0|12) -> (1|02)
#pragma omp parallel for num_threads(nthreads_)
                for (size_t i = 0; i < M0; i++) {
                    for (size_t j = 0; j < M1; j++) {
#pragma omp simd
                        for (size_t k = 0; k < M2; k++) {
                            Fp[j * M0 * M2 + i * M2 + k] = Mp[i * M1 * M2 + j * M2 + k];
                        }
                    }
                }
                put_tensor(new_filename, Fp, std::make_pair(0, M1 - 1), std::make_pair(start, stop),
                           std::make_pair(0, M2 - 1), op);
            } else if (a1 == 2) {  // (0|12) -> (1|20)
#pragma omp parallel for num_threads(nthreads_)
                for (size_t i = 0; i < M0; i++) {
                    for (size_t j = 0; j < M1; j++) {
                        for (size_t k = 0; k < M2; k++) {
                            Fp[j * M0 * M2 + k * M0 + i] = Mp[i * M1 * M2 + j * M2 + k];
                        }
                    }
                }
                put_tensor(new_filename, Fp, std::make_pair(0, M1 - 1), std::make_pair(0, M2 - 1),
                           std::make_pair(start, stop), op);
            }
        } else if (a0 == 2) {
            if (a1 == 0) {  // (0|12) -> (2|01)
#pragma omp parallel for num_threads(nthreads_)
                for (size_t i = 0; i < M0; i++) {
                    for (size_t j = 0; j < M1; j++) {
                        for (size_t k = 0; k < M2; k++) {
                            Fp[k * M1 * M0 + i * M1 + j] = Mp[i * M1 * M2 + j * M2 + k];
                        }
                    }
                }
                put_tensor(new_filename, Fp, std::make_pair(0, M2 - 1), std::make_pair(start, stop),
                           std::make_pair(0, M1 - 1), op);
            } else if (a1 == 1) {  // (0|12) -> (2|10)
#pragma omp parallel for num_threads(nthreads_)
                for (size_t i = 0; i < M0; i++) {
                    for (size_t j = 0; j < M1; j++) {
                        for (size_t k = 0; k < M2; k++) {
                            Fp[k * M1 * M0 + j * M0 + i] = Mp[i * M1 * M2 + j * M2 + k];
                        }
                    }
                }
                put_tensor(new_filename, Fp, std::make_pair(0, M2 - 1), std::make_pair(0, M1 - 1),
                           std::make_pair(start, stop), op);
            }
        }
    }
    // better be careful
    remove(filename.c_str());
    rename(new_filename.c_str(), filename.c_str());
    file_streams_[filename] = file_streams_[new_filename];
    stream_check(filename, "rb");
    file_streams_.erase(new_filename);

    // keep tsizes_ separate and do not ovwrt sizes_ in case of STORE directive
    files_.erase(new_file);
    tsizes_[filename] = sizes;
}
size_t DFHelper::get_space_size(std::string name) {
    if (spaces_.find(name) == spaces_.end()) {
        std::stringstream error;
        error << "DFHelper:get_space_size: " << name << " not found.";
        throw PSIEXCEPTION(error.str().c_str());
    }
    return std::get<1>(spaces_[name]);
}
size_t DFHelper::get_tensor_size(std::string name) {
    if (transf_.find(name) == transf_.end()) {
        std::stringstream error;
        error << "DFHelper:get_tensor_size: " << name << " not found.";
        throw PSIEXCEPTION(error.str().c_str());
    }
    std::tuple<size_t, size_t, size_t> s = sizes_[std::get<1>(files_[name])];
    return std::get<0>(s) * std::get<1>(s) * std::get<2>(s);
}
std::tuple<size_t, size_t, size_t> DFHelper::get_tensor_shape(std::string name) {
    if (transf_.find(name) == transf_.end()) {
        std::stringstream error;
        error << "DFHelper:get_tensor_size: " << name << " not found.";
        throw PSIEXCEPTION(error.str().c_str());
    }
    return sizes_[std::get<1>(files_[name])];
}
void DFHelper::build_JK(std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright, std::vector<SharedMatrix> D,
                        std::vector<SharedMatrix> J, std::vector<SharedMatrix> K, std::vector<SharedMatrix> wK,
                        size_t max_nocc, bool do_J, bool do_K, bool do_wK, bool lr_symmetric) {
    if (debug_) {
        outfile->Printf("Entering DFHelper::build_JK\n");
    }
    bool store_k = do_K;
    if ( do_wK && wcombine_ ) { do_K = false; } 

    // This was an if-else statement. Presumably, we could manage J construction
    //   to more effectively manage memory, so I think that was what was going on.
    if (do_J || do_K) {
        timer_on("DFH: compute_JK()");
        compute_JK(Cleft, Cright, D, J, K, max_nocc, do_J, do_K, do_wK, lr_symmetric);
        timer_off("DFH: compute_JK()");
    }

    if (do_wK_) {
        timer_on("DFH: compute_wK()");
        compute_wK(Cleft, Cright, wK, max_nocc, do_J, do_K, do_wK);
        timer_off("DFH: compute_wK()");
    }

    if (store_k) { do_K = true;}

    if (debug_) {
        outfile->Printf("Exiting DFHelper::build_JK\n");
    }
}
void DFHelper::compute_JK(std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright,
                          std::vector<SharedMatrix> D, std::vector<SharedMatrix> J, std::vector<SharedMatrix> K,
                          size_t max_nocc, bool do_J, bool do_K, bool do_wK, bool lr_symmetric) {
    // outfile->Printf("\n     ==> DFHelper:--Begin J/K builds <==\n\n");
    // outfile->Printf("\n     ==> Using the %s directive with AO_CORE = %d <==\n\n", method_.c_str(), AO_core_);

    // size checks for C matrices occur in jk.cc
    // computing D occurs inside of jk.cc

    // determine buffer sizes and blocking scheme
    // would love to move this to initialize(), but
    // we would need to know max_nocc_ beforehand
    // two major advantages would arise from moving this
    // 1. we could predict the blocks used, write them to different files, and improve IO, otherwise
    // the strided disk reads for the AOs will result in a definite loss to DiskDFJK in the disk-bound realm
    // 2. we could allocate the buffers only once, instead of every time compute_JK() is called

    // Each element of Qsteps specifies the endpoints of a batch of auxiliary shells.
    // We'll treat all (PN|Q) for Q in this batch at once.
    std::vector<std::pair<size_t, size_t>> Qsteps;
    std::tuple<size_t, size_t> info = Qshell_blocks_for_JK_build(Qsteps, max_nocc, lr_symmetric);
    size_t tots = std::get<0>(info);
    size_t totsb = std::get<1>(info);

    // prep stream, blocking
    if (!direct_ && !AO_core_) stream_check(AO_names_[1], "rb");

    std::vector<std::vector<double>> C_buffers(nthreads_);

// prepare C buffers
#pragma omp parallel num_threads(nthreads_)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        C_buffers[rank] = std::vector<double>(nbf_ * std::max(max_nocc, nbf_));
    }

    // declare bufs
    std::unique_ptr<double[]> M;   // AOs
    std::unique_ptr<double[]> T1;  // Ktmp1
    std::unique_ptr<double[]> T2;  // Ktmp2

    // allocate first Ktmp
    size_t Ktmp_size = (!max_nocc ? totsb * 1 : totsb * max_nocc);
    Ktmp_size = std::max(Ktmp_size * nbf_, nthreads_ * naux_);  // max necessary
    Ktmp_size = std::max(Ktmp_size, nbf_ * nbf_); 
    T1 = std::make_unique<double[]>(Ktmp_size);
    double* T1p = T1.get();

    // if lr_symmetric, we can be more clever with mem usage. T2 is used for both the
    // second tmp in the K build, as well as the completed, pruned J build.
    if (lr_symmetric) {
        Ktmp_size = nbf_ * nbf_;  // size for pruned J build
    } else {
        Ktmp_size = std::max(nbf_ * nbf_, Ktmp_size);  // max necessary
    }

    // necessary for wcombine case for small molecules
    Ktmp_size = std::max(Ktmp_size, nbf_ * nbf_); 
    Ktmp_size = std::max(Ktmp_size, nthreads_ * naux_);

    T2 = std::make_unique<double[]>(Ktmp_size);
    double* T2p = T2.get();

    double* Mp;
    if (!AO_core_) {
        M = std::make_unique<double[]>(tots);
        Mp = M.get();
    } else
        if (!wcombine_) {Mp = Ppq_.get();}

    double* M1p;
    if (wcombine_) {
        M1p = m1Ppq_.get();
    }

    // Transform a single batch of integrals
    size_t bcount = 0;
    for (const auto& Qstep :Qsteps) {
        // Qshell step info
        auto start = std::get<0>(Qstep);
        auto stop = std::get<1>(Qstep);
        auto begin = Qshell_aggs_[start];
        auto end = Qshell_aggs_[stop + 1] - 1;
        auto block_size = end - begin + 1;

        // get AO chunk according to directive
        timer_on("DFH: Grabbing AOs");
        if (!AO_core_) {
            grab_AO(start, stop, Mp);
        }
        timer_off("DFH: Grabbing AOs");
        if (do_J) {
            timer_on("DFH: compute_J");
            if (wcombine_ && do_wK_) {
                compute_J_combined(D, J, M1p, T1p, T2p, C_buffers, bcount, block_size);
            } else {
                if (lr_symmetric) {
                    compute_J_symm(D, J, Mp, T1p, T2p, C_buffers, bcount, block_size);
                } else {
                    compute_J(D, J, Mp, T1p, T2p, C_buffers, bcount, block_size);
                }
            }
            timer_off("DFH: compute_J");
        }

        if (do_K) {
            timer_on("DFH: compute_K");
            compute_K(Cleft, Cright, K, T1p, T2p, Mp, bcount, block_size, C_buffers, lr_symmetric);
            timer_off("DFH: compute_K");
        }

        bcount += block_size;
    }
    // outfile->Printf("\n     ==> DFHelper:--End J/K Builds (disk)<==\n\n");
}
void DFHelper::compute_J_symm(std::vector<SharedMatrix> D, std::vector<SharedMatrix> J, double* Mp, double* T1p,
                              double* T2p, std::vector<std::vector<double>>& D_buffers, size_t bcount,
                              size_t block_size) {
    for (size_t i = 0; i < J.size(); i++) {
        // grab orbital spaces
        double* Dp = D[i]->pointer()[0];
        double* Jp = J[i]->pointer()[0];

        // initialize Tmp (pQ)
        fill(T1p, nthreads_ * naux_, 0.0);

#pragma omp parallel for schedule(guided) num_threads(nthreads_)
        for (size_t k = 0; k < nbf_; k++) {
            size_t si = small_skips_[k];
            size_t mi = symm_small_skips_[k];
            size_t skip = symm_ignored_columns_[k];
            size_t jump = (AO_core_ ? big_skips_[k] + bcount * si : (big_skips_[k] * block_size) / naux_);

            int rank = 0;
#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif

            for (size_t m = k, sp_count = -1; m < nbf_; m++) {
                if (schwarz_fun_index_[k * nbf_ + m]) {
                    sp_count++;
                    D_buffers[rank][sp_count] = (m == k ? Dp[nbf_ * k + m] : 2 * Dp[nbf_ * k + m]);
                }
            }

            // (Qm)(m) -> (Q)
            C_DGEMV('N', block_size, mi, 1.0, &Mp[jump + skip], si, &D_buffers[rank][0], 1, 1.0, &T1p[rank * naux_], 1);
        }

        // reduce
        for (size_t k = 1; k < nthreads_; k++) {
            for (size_t l = 0; l < naux_; l++) T1p[l] += T1p[k * naux_ + l];
        }

// complete pruned J
#pragma omp parallel for schedule(guided) num_threads(nthreads_)
        for (size_t k = 0; k < nbf_; k++) {
            size_t si = small_skips_[k];
            size_t mi = symm_small_skips_[k];
            size_t skip = symm_ignored_columns_[k];
            size_t jump = (AO_core_ ? big_skips_[k] + bcount * si : (big_skips_[k] * block_size) / naux_);
            C_DGEMV('T', block_size, mi, 1.0, &Mp[jump + skip], si, T1p, 1, 0.0, &T2p[k * nbf_], 1);
        }

        // unpack from sparse to dense
        for (size_t k = 0; k < nbf_; k++) {
            for (size_t m = k + 1, count = 0; m < nbf_; m++) {  // assumes diagonal exists to avoid if  FIXME
                if (schwarz_fun_index_[k * nbf_ + m]) {
                    count++;
                    Jp[k * nbf_ + m] += T2p[k * nbf_ + count];
                    Jp[m * nbf_ + k] += T2p[k * nbf_ + count];
                }
            }
        }
        for (size_t k = 0; k < nbf_; k++) Jp[k * nbf_ + k] += T2p[k * nbf_];
    }
}
void DFHelper::fill(double* b, size_t count, double value) {
#pragma omp parallel for simd num_threads(nthreads_) schedule(static)
    for (size_t i = 0; i < count; i++) {
        b[i] = value;
    }
}
void DFHelper::compute_J(const std::vector<SharedMatrix> D, std::vector<SharedMatrix> J, double* Mp, double* T1p, double* T2p,
                         std::vector<std::vector<double>>& D_buffers, size_t bcount, size_t block_size) {
    for (size_t i = 0; i < J.size(); i++) {
        // grab orbital spaces
        auto Dp = D[i]->pointer()[0];
        auto Jp = J[i]->pointer()[0];

        // initialize Tmp (pQ)
        // TODO: Make T1p a std::vector, so we don't need to know the length.
        fill(T1p, nthreads_ * naux_, 0.0);

#pragma omp parallel for schedule(guided) num_threads(nthreads_)
        for (size_t k = 0; k < nbf_; k++) {
            size_t sp_size = small_skips_[k];
            size_t jump = (AO_core_ ? big_skips_[k] + bcount * sp_size : (big_skips_[k] * block_size) / naux_);

            int rank = 0;
#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif

            for (size_t m = 0, sp_count = -1; m < nbf_; m++) {
                if (schwarz_fun_index_[k * nbf_ + m]) {
                    sp_count++;
                    D_buffers[rank][sp_count] = Dp[nbf_ * k + m];
                }
            }
            // (Qm)(m) -> (Q)
            C_DGEMV('N', block_size, sp_size, 1.0, &Mp[jump], sp_size, &D_buffers[rank][0], 1, 1.0, &T1p[rank * naux_],
                    1);
        }

        // reduce
        for (size_t k = 1; k < nthreads_; k++) {
            for (size_t l = 0; l < naux_; l++) T1p[l] += T1p[k * naux_ + l];
        }

// complete pruned J
#pragma omp parallel for schedule(guided) num_threads(nthreads_)
        for (size_t k = 0; k < nbf_; k++) {
            size_t sp_size = small_skips_[k];
            size_t jump = (AO_core_ ? big_skips_[k] + bcount * sp_size : (big_skips_[k] * block_size) / naux_);
            C_DGEMV('T', block_size, sp_size, 1.0, &Mp[jump], sp_size, T1p, 1, 0.0, &T2p[k * nbf_], 1);
        }

        // unpack from sparse to dense
        for (size_t k = 0; k < nbf_; k++) {
            for (size_t m = 0, count = -1; m < nbf_; m++) {
                if (schwarz_fun_index_[k * nbf_ + m]) {
                    count++;
                    Jp[k * nbf_ + m] += T2p[k * nbf_ + count];
                }
            }
        }
    }
}
void DFHelper::compute_J_combined(std::vector<SharedMatrix> D, std::vector<SharedMatrix> J, double* Mp, double* T1p, double* T2p, std::vector<std::vector<double>>& D_buffers, size_t bcount, size_t block_size) {
//    double* coul_metp = metric_prep_core(1.0);
    double* coul_metp;
    std::unique_ptr<double[]> metric;


    if (!hold_met_) {
        metric = std::make_unique<double[]>(naux_ * naux_);
        coul_metp = metric.get();
        std::string filename = return_metfile(1.0);
        get_tensor_(std::get<0>(files_[filename]), coul_metp, 0, naux_ - 1, 0, naux_ - 1);
    } else
        coul_metp = metric_prep_core(1.0);


    for (size_t i = 0; i < J.size(); i++) {
        double* Dp = D[i]->pointer()[0];
        double* Jp = J[i]->pointer()[0];
        memset(T1p, 0, naux_ *  nthreads_ * sizeof(double));
        memset(T2p, 0, naux_ * sizeof(double));

#pragma omp parallel for schedule(guided) num_threads(nthreads_)
        for (size_t k = 0; k < nbf_; k++) {
            size_t sp_size = small_skips_[k];
            size_t jump = (AO_core_ ? big_skips_[k] + bcount * sp_size : (big_skips_[k] * block_size) / naux_ );
            int rank = 0;
#ifdef _OPENMP 
            rank = omp_get_thread_num();
#endif
            for (size_t m = 0, sp_count = -1; m < nbf_; m++) {
                if (schwarz_fun_index_[k * nbf_ + m]) {
                    sp_count++;
                    D_buffers[rank][sp_count] = Dp[nbf_ * m + k];
                }
            }
            C_DGEMV('N', block_size, sp_size, 1.0, &Mp[jump], sp_size, &D_buffers[rank][0], 1, 1.0, &T1p[rank * naux_], 1);
        }
        //reduce
        for (size_t k = 1; k < nthreads_; k++) {
            for (size_t l = 0; l < naux_; l++) { T1p[l] += T1p[ k * naux_ + l ]; }
        }

        //metric contraction
        C_DGEMV('N', naux_, naux_, 1.0, coul_metp, naux_, T1p, 1, 0.0, T2p, 1);

        memset(T1p, 0, nbf_ * nbf_ * sizeof(double));
        // comute pruned J
#pragma omp parallel for schedule(guided) num_threads(nthreads_)
        for (size_t k = 0; k < nbf_; k++ ) {
            size_t sp_size = small_skips_[k];
            size_t jump = (AO_core_ ? big_skips_[k] + bcount * sp_size : (big_skips_[k] * block_size) / naux_ );
            C_DGEMV('T', block_size, sp_size, 1.0, &Mp[jump], sp_size, T2p, 1, 0.0, &T1p[k*nbf_], 1);
        }

        for (size_t k = 0; k < nbf_; k++) {
            for (size_t m = 0, count = -1; m < nbf_; m++) {
                if (schwarz_fun_index_[k * nbf_ + m]) {
                    count++;
                    Jp[k * nbf_ + m ] += T1p[k * nbf_ + count];
                }
            }
        }
    }
}
void DFHelper::compute_K(std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright, std::vector<SharedMatrix> K,
                         double* T1p, double* T2p, double* Mp, size_t bcount, size_t block_size,
                         std::vector<std::vector<double>>& C_buffers, bool lr_symmetric) {
    for (size_t i = 0; i < K.size(); i++) {
        size_t nocc = Cleft[i]->colspi()[0];
        if (!nocc) {
            continue;
        }

        double* Clp = Cleft[i]->pointer()[0];
        double* Crp = Cright[i]->pointer()[0];
        double* Kp = K[i]->pointer()[0];

        // compute first tmp
        first_transform_pQq(nocc, bcount, block_size, Mp, T1p, Clp, C_buffers);

        // compute second tmp
        if (lr_symmetric) {
            T2p = T1p;
        } else {
            first_transform_pQq(nocc, bcount, block_size, Mp, T2p, Crp, C_buffers);
        }

        // compute K
        C_DGEMM('N', 'T', nbf_, nbf_, nocc * block_size, 1.0, T1p, nocc * block_size, T2p, nocc * block_size, 1.0, Kp,
                nbf_);
    }
}
void DFHelper::compute_wK(std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright,
                          std::vector<SharedMatrix> wK, size_t max_nocc, bool do_J, bool do_K, bool do_wK) {
    std::vector<std::pair<size_t, size_t>> Qsteps;
    std::tuple<size_t, size_t> info = Qshell_blocks_for_JK_build(Qsteps, max_nocc, false);
    size_t tots = std::get<0>(info);
    size_t totsb = std::get<1>(info);

    double* wMp = wPpq_.get();
    double* M1p = m1Ppq_.get();
    std::vector<std::vector<double>> C_buffers(nthreads_);
// prepare C buffers
#pragma omp parallel num_threads(nthreads_)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        C_buffers[rank] = std::vector<double>(nbf_ * std::max(max_nocc, nbf_));
    }

    // declare bufs
    std::unique_ptr<double[]> T1;  // Ktmp1
    std::unique_ptr<double[]> T2;  // Ktmp2

    // allocate first Ktmp
    size_t Ktmp_size = (!max_nocc ? totsb * 1 : totsb * max_nocc);
    Ktmp_size = std::max(Ktmp_size * nbf_, nthreads_ * naux_);  // max necessary
    T1 = std::make_unique<double[]>(Ktmp_size);
    double* T1p = T1.get();
    T2 = std::make_unique<double[]>(Ktmp_size);
    double* T2p = T2.get();

    /* The rest of the file would usually be in compute_K */
    /* */
    for (size_t bind = 0, bcount = 0; bind < Qsteps.size(); bind++) {
        size_t start = std::get<0>(Qsteps[bind]);
        size_t stop = std::get<1>(Qsteps[bind]);
        size_t begin = Qshell_aggs_[start];
        size_t end = Qshell_aggs_[stop + 1] - 1;
        size_t block_size = end - begin + 1;
        for (size_t i = 0; i < wK.size(); i++) {
            size_t nocc = Cleft[i]->colspi()[0];
            if (!nocc) {
                continue;
            }
            double* Clp = Cleft[i]->pointer()[0];
            double* Crp = Cright[i]->pointer()[0];
            double* wKp = wK[i]->pointer()[0];

            // compute first tmp
            first_transform_pQq(nocc, bcount, block_size, M1p, T1p, Clp, C_buffers);
            // compute second tmp
            first_transform_pQq(nocc, bcount, block_size, wMp, T2p, Crp, C_buffers);

            // compute wK
            C_DGEMM('N', 'T', nbf_, nbf_, nocc * block_size, 1.0, T1p, nocc * block_size, T2p, nocc * block_size, 1.0,
                    wKp, nbf_);
        }
        bcount += block_size;
    }
}

}  // End namespaces
