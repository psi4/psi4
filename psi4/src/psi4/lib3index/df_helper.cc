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

#include "df_helper.h"

#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libfock/jk.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/sieve.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"

#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

DF_Helper::DF_Helper(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> aux) {
    primary_ = primary;
    aux_ = aux;
    nao_ = primary_->nbf();
    naux_ = aux_->nbf();
    wMO_ = nao_ / 2;
    prepare_blocking();
}
DF_Helper::~DF_Helper() {
    // close streams
    for (auto& kv : file_status_) fclose(kv.second.fp);

    // destroy all files
    std::string file1, file2;
    for (auto& kv : files_) {
        remove(std::get<0>(kv.second).c_str());
        remove(std::get<1>(kv.second).c_str());
    }
    for (auto& kv : AO_files_) remove(kv.second.c_str());
}
void DF_Helper::prepare_blocking() {
    Qshells_ = aux_->nshell();
    pshells_ = primary_->nshell();

    Qshell_aggs_.reserve(Qshells_ + 1);
    pshell_aggs_.reserve(pshells_ + 1);

    // Aux shell blocking
    Qshell_aggs_[0] = 0;
    for (size_t i = 0; i < Qshells_; i++) Qshell_aggs_[i + 1] = Qshell_aggs_[i] + aux_->shell(i).nfunction();

    // AO shell blocking
    pshell_aggs_[0] = 0;
    for (size_t i = 0; i < pshells_; i++) pshell_aggs_[i + 1] = pshell_aggs_[i] + primary_->shell(i).nfunction();
}
void DF_Helper::AO_filename_maker(size_t i) {
    std::string name = PSIOManager::shared_object()->get_default_path();
    name.append("dfh.AO");
    name.append(std::to_string(i));
    AO_names_.push_back(name);
    std::string file = name;
    file.append(".");
    file.append(std::to_string(getpid()));
    file.append(".");
    file.append(primary_->molecule()->name());
    file.append(".dat");
    AO_files_[name] = file;
}
void DF_Helper::filename_maker(std::string name, size_t a0, size_t a1, size_t a2, size_t op) {
    std::string pfile = "dfh.p";
    pfile.append(name);
    pfile.append(".");
    pfile.append(std::to_string(getpid()));
    pfile.append(".");
    pfile.append(primary_->molecule()->name());
    pfile.append(".dat");
    std::string file = pfile;
    file.erase(4, 1);

    std::string scratch = PSIOManager::shared_object()->get_default_path();
    std::string pfilename = scratch;
    pfilename.append(pfile);
    std::string filename = scratch;
    filename.append(file);

    std::tuple<std::string, std::string> files(pfilename.c_str(), filename.c_str());
    files_[name] = files;

    std::tuple<size_t, size_t, size_t> sizes;
    if (op) {
        sizes = std::make_tuple(a1, a2, a0);
    } else {
        sizes = std::make_tuple(a0, a1, a2);
    }

    sizes_[pfilename] = sizes;
    sizes_[filename] = sizes;
}
void DF_Helper::initialize() {
    
    timer_on("DFH: initialize()");
    
    // have the algorithm specified before init
    if (method_.compare("DIRECT") && method_.compare("STORE")) {
        std::stringstream error;
        error << "DF_Helper:initialize: specified method (" << method_ << ") is incorrect";
        throw PSIEXCEPTION(error.str().c_str());
    }
    direct_ = (!method_.compare("DIRECT") ? true : false);

    // if metric power is not zero, prepare it
    if (!(std::fabs(mpower_ - 0.0) < 1e-13)) (hold_met_ ? prepare_metric_core() : prepare_metric());

    // prepare sparsity masks
    timer_on("DFH: sparsity prep");
    prepare_sparsity();
    timer_off("DFH: sparsity prep");

    // shut off AO core if necessary
    if (memory_ * 0.8 < big_skips_[nao_]) {
        AO_core_ = false;
        outfile->Printf("\n    ==> DF_Helper memory announcement <==\n");
        std::stringstream ann;
        ann << "        Turning off in-core AOs.  DF_Helper needs at least (";
        ann << big_skips_[nao_] / 0.9 * 8 / 1e9 << "GB), but it only got (" << memory_ * 8 / 1e9 << "GB)\n";
        outfile->Printf(ann.str().c_str());
    }

    // prepare AOs for STORE directive
    if (AO_core_)
        prepare_AO_core();
    else if (!direct_)
        prepare_AO();

    built = true;
    timer_off("DFH: initialize()");
}
void DF_Helper::print_header() {
    outfile->Printf("\n    ==> DF_Helper <==\n");
    outfile->Printf("      nao           = %zu\n", nao_);
    outfile->Printf("      naux          = %zu\n", naux_);
    outfile->Printf("      Scwarz cutoff = %E\n", cutoff_);
    outfile->Printf("      mem (doubles) = %zu\n", memory_);
    outfile->Printf("      nthreads      = %zu\n", nthreads_);
    outfile->Printf("      Algorithm     = %s\n", method_.c_str());
    outfile->Printf("      AO_core?      = %d\n", (int)AO_core_);
    outfile->Printf("      MO_core?      = %d\n", (int)MO_core_);
    outfile->Printf("      hold metric   = %d\n", (int)hold_met_);
    outfile->Printf("      metric power  = %E\n", mpower_);
    outfile->Printf("  Fitting condition = %E\n", condition_);
    outfile->Printf("  Mask sparsity (%%) = (%f)\n", 100. * (1.0 - (double)small_skips_[nao_] / (double)(nao_ * nao_)));
    outfile->Printf("\n\n");
}
void DF_Helper::prepare_sparsity() {
    // prep info vectors
    std::vector<double> fun_prints(nao_ * nao_, 0.0);
    std::vector<double> shell_prints(pshells_ * pshells_, 0.0);
    schwarz_fun_mask_.reserve(nao_ * nao_);
    schwarz_shell_mask_.reserve(pshells_ * pshells_);
    small_skips_.reserve(nao_ + 1);
    big_skips_.reserve(nao_ + 1);

    // prepare eri buffers
    size_t nthreads = nthreads_;  // for now
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(primary_, primary_, primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthreads);
    std::vector<const double*> buffer(nthreads);

    int rank = 0;
#pragma omp parallel for private(rank) num_threads(nthreads) if (nao_ > 1000)
    for (size_t i = 0; i < nthreads; i++) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
        buffer[rank] = eri[rank]->buffer();
    }

    double val, max_val = 0.0;
    size_t MU, NU, mu, nu, omu, onu, nummu, numnu, index;
#pragma omp parallel for private(MU, NU, mu, nu, omu, onu, nummu, numnu, index, val, \
                                 rank) num_threads(nthreads_) if (nao_ > 1000) schedule(guided)
    for (MU = 0; MU < pshells_; ++MU) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        nummu = primary_->shell(MU).nfunction();
        for (NU = 0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU).nfunction();
            eri[rank]->compute_shell(MU, NU, MU, NU);
            for (mu = 0; mu < nummu; ++mu) {
                omu = primary_->shell(MU).function_index() + mu;
                for (nu = 0; nu < numnu; ++nu) {
                    onu = primary_->shell(NU).function_index() + nu;
                    if (omu >= onu) {
                        index = mu * (numnu * nummu * numnu + numnu) + nu * (nummu * numnu + 1);
                        val = fabs(buffer[rank][index]);
                        max_val = (max_val < val ? val : max_val);
                        if (shell_prints[MU * pshells_ + NU] <= val) {
                            shell_prints[MU * pshells_ + NU] = val;
                            shell_prints[NU * pshells_ + MU] = val;
                        }
                        if (fun_prints[omu * nao_ + onu] <= val) {
                            fun_prints[omu * nao_ + onu] = val;
                            fun_prints[onu * nao_ + omu] = val;
                        }
                    }
                }
            }
        }
    }
    tolerance_ = cutoff_ * cutoff_ / max_val;

    //#pragma omp parallel for simd num_threads(nthreads_) schedule(static)
    for (size_t i = 0; i < pshells_ * pshells_; i++) schwarz_shell_mask_[i] = (shell_prints[i] < tolerance_ ? 0 : 1);

    //#pragma omp parallel for private(count) num_threads(nthreads_)
    for (size_t i = 0, count = 0; i < nao_; i++) {
        count = 0;
        for (size_t j = 0; j < nao_; j++) {
            if (fun_prints[i * nao_ + j] >= tolerance_) {
                count++;
                schwarz_fun_mask_[i * nao_ + j] = count;
            } else
                schwarz_fun_mask_[i * nao_ + j] = 0;
        }
        small_skips_[i] = count;
    }

    // build indexing skips
    big_skips_[0] = 0;
    size_t coltots = 0;
    for (size_t j = 0; j < nao_; j++) {
        size_t cols = small_skips_[j];
        size_t size = cols * naux_;
        coltots += cols;
        big_skips_[j + 1] = size + big_skips_[j];
    }
    small_skips_[nao_] = coltots;

    symm_skips_.reserve(nao_);
    symm_sizes_.reserve(nao_);
    for (size_t i = 0; i < nao_; i++) {
        size_t size = 0;
        size_t skip = 0;
        for (size_t j = 0; j < nao_; j++)
            if (schwarz_fun_mask_[i * nao_ + j]) (j >= i ? size++ : skip++);
        symm_sizes_[i] = size;
        symm_skips_[i] = skip;
    }

    symm_agg_sizes_.reserve(nao_ + 1);
    symm_agg_sizes_[0] = 0;
    for (size_t i = 1; i < nao_ + 1; i++) symm_agg_sizes_[i] = symm_agg_sizes_[i - 1] + symm_sizes_[i - 1];
}
void DF_Helper::prepare_AO() {
    // prepare eris
    size_t rank = 0;
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(aux_, zero, primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthreads_);
#pragma omp parallel for schedule(static) num_threads(nthreads_) private(rank)
    for (size_t i = 0; i < nthreads_; i++) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    }

    // gather blocking info
    std::vector<std::pair<size_t, size_t>> psteps;
    std::pair<size_t, size_t> plargest = pshell_blocks_for_AO_build(memory_, 0, psteps);

    // declare largest necessary
    std::vector<double> M;
    std::vector<double> F;
    std::vector<double> metric;
    M.reserve(std::get<0>(plargest));
    F.reserve(std::get<0>(plargest));
    double* Mp = M.data();
    double* Fp = F.data();

    // grab metric
    double* metp;
    if (!hold_met_) {
        metric.reserve(naux_ * naux_);
        metp = metric.data();
        std::string filename = return_metfile(mpower_);
        get_tensor_(std::get<0>(files_[filename]), metp, 0, naux_ - 1, 0, naux_ - 1);

    } else
        metp = metric_prep_core(mpower_);

    // prepare files
    AO_filename_maker(1);
    AO_filename_maker(2);
    std::string putf = AO_files_[AO_names_[1]];
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
        compute_AO_p(start, stop, Mp, eri);
        timer_off("DFH: AO Construction");
        timer_on("DFH: AO-Met. Contraction");
        
        // loop and contract
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
void DF_Helper::prepare_AO_core() {
    // prepare eris
    size_t rank = 0;
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(aux_, zero, primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthreads_);
#pragma omp parallel for schedule(static) num_threads(nthreads_) private(rank)
    for (size_t i = 0; i < nthreads_; i++) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    }

    // determine blocking for symmetric scheme
    std::vector<std::pair<size_t, size_t>> psteps;
    std::pair<size_t, size_t> plargest = pshell_blocks_for_AO_build(memory_, 1, psteps);

    // allocate final AO vector
    Ppq_.reserve(big_skips_[nao_]);

    // outfile->Printf("\n    ==> Begin AO Blocked Construction <==\n\n");
    if (!direct_) {
        // declare sparse buffer
        std::vector<double> Qpq;
        Qpq.reserve(std::get<0>(plargest));
        double* Mp = Qpq.data();
        double* metp;
        std::vector<double> metric;

        if (!hold_met_) {
            metric.reserve(naux_ * naux_);
            metp = metric.data();
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
            compute_AO_p_symm(start, stop, Mp, eri);
            timer_off("DFH: AO Construction");

            // contract metric
            timer_on("DFH: AO-Met. Contraction");
            contract_metric_AO_core_symm(Mp, metp, begin, end);
            timer_off("DFH: AO-Met. Contraction");
        }
        // no more need for metrics
        if (hold_met_) metrics_.clear();
    } else {
        timer_on("DFH: AO Construction");
        compute_AO_p(0, pshells_ - 1, &Ppq_[0], eri);
        timer_off("DFH: AO Construction");
    }
    // outfile->Printf("\n    ==> End AO Blocked Construction <==");
}
std::pair<size_t, size_t> DF_Helper::pshell_blocks_for_AO_build(const size_t mem, size_t symm,
                                                                std::vector<std::pair<size_t, size_t>>& b) {
    size_t full_3index = big_skips_[nao_];
    size_t constraint, end, begin, current, block_size, tmpbs, total, count, largest;
    block_size = tmpbs = total = count = largest = 0;
    for (size_t i = 0; i < pshells_; i++) {
        count++;
        begin = pshell_aggs_[i];
        end = pshell_aggs_[i + 1] - 1;
        tmpbs += end - begin + 1;

        if (symm) {  // in-core symmetric
            current = symm_agg_sizes_[end + 1] - symm_agg_sizes_[begin];
            total += current * naux_;
        } else {  // on-disk
            current = big_skips_[end + 1] - big_skips_[begin];
            total += 2 * current;
        }

        constraint = (hold_met_ ? (total + naux_ * naux_) : total);
        constraint += (symm ? full_3index : 0);
        if (constraint > mem || i == pshells_ - 1) {
            if (count == 1 && i != pshells_ - 1) {
                std::stringstream error;
                error << "DF_Helper: not enough memory for (p shell) AO blocking!"
                      << " required memory: " << constraint * 8 / 1e9 << "GB.";
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
std::pair<size_t, size_t> DF_Helper::Qshell_blocks_for_transform(const size_t mem, size_t wtmp, size_t wfinal,
                                                                 std::vector<std::pair<size_t, size_t>>& b) {
    size_t extra = (hold_met_ ? naux_ * naux_ : 0);
    size_t end, begin, current, block_size, tmpbs, total, count, largest;
    block_size = tmpbs = total = count = largest = 0;
    for (size_t i = 0; i < Qshells_; i++) {
        count++;
        begin = Qshell_aggs_[i];
        end = Qshell_aggs_[i + 1] - 1;
        tmpbs += end - begin + 1;
        current = (end - begin + 1) * small_skips_[nao_];
        total += current;
        total = (AO_core_ ? big_skips_[nao_] : total);

        size_t constraint = total + (wtmp * nao_ + 2 * wfinal) * tmpbs + extra;
        // AOs + worst half transformed + worst final
        if (constraint > mem || i == Qshells_ - 1) {
            if (count == 1 && i != Qshells_ - 1) {
                std::stringstream error;
                error << "DF_Helper: not enough memory for transformation blocking!";
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
    // returns pair(largest buffer size, largest block size)
    return std::make_pair(largest, block_size);
}
std::tuple<size_t, size_t, size_t, size_t> DF_Helper::Qshell_blocks_for_JK_build(
    std::vector<std::pair<size_t, size_t>>& b, std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright) {
    // if wcleft != wcright something went wrong
    size_t wcleft = 0, wcright = 0;
    for (size_t i = 0; i < Cleft.size(); i++) {
        wcleft = (Cleft[i]->colspi()[0] > wcleft ? Cleft[i]->colspi()[0] : wcleft);
        wcright = (Cright[i]->colspi()[0] > wcright ? Cright[i]->colspi()[0] : wcright);
    }

    // K tmps
    size_t T1 = nao_ * wcleft;
    size_t T2 = (JK_hint_ ? nao_ * nao_ : nao_ * wcright);

    // C_buffers
    size_t T3 = nthreads_ * nao_ * nao_;

    size_t current, block_size = 0, tmpbs = 0, count = 0, largest = 0;
    size_t total = (AO_core_ ? big_skips_[nao_] : 0);
    for (size_t i = 0; i < Qshells_; i++) {
        count++;
        size_t begin = Qshell_aggs_[i];
        size_t end = Qshell_aggs_[i + 1] - 1;
        current = (end - begin + 1) * small_skips_[nao_];
        total += (AO_core_ ? 0 : current);
        tmpbs += end - begin + 1;

        size_t constraint = total + T1 * tmpbs + T3;
        constraint += (JK_hint_ ? T2 : T2 * tmpbs);
        if (constraint > memory_ || i == Qshells_ - 1) {
            if (count == 1 && i != Qshells_ - 1) {
                std::stringstream error;
                error << "DF_Helper: not enough memory for JK blocking!";
                throw PSIEXCEPTION(error.str().c_str());
            }
            if (constraint > memory_) {
                total -= current;
                tmpbs -= end - begin + 1;
                b.push_back(std::make_pair(i - count + 1, i - 1));
                i--;
            } else if (i == Qshells_ - 1)
                b.push_back(std::make_pair(i - count + 1, i));
            if (block_size < tmpbs) {
                largest = total;
                block_size = tmpbs;
            }
            count = 0;
            total = 0;
            tmpbs = 0;
        }
    }
    // returns tuple(largest buffer size, largest block size, wcleft, wcright)
    return std::make_tuple(largest, block_size, wcleft, wcright);
}
FILE* DF_Helper::stream_check(std::string filename, std::string op) {
    // timer_on("stream checks             ");
    if (file_status_.find(filename) == file_status_.end()) {
        // outfile->Printf("file not found: %s\n", filename.c_str());
        file_status_[filename].op = op;
        file_status_[filename].fp = fopen(filename.c_str(), op.c_str());
    } else {
        if (op.compare(file_status_[filename].op)) {
            // timer_on("stream change           ");
            // outfile->Printf("changing......: %s\n", filename.c_str());
            file_status_[filename].op = op;
            fflush(file_status_[filename].fp);
            fclose(file_status_[filename].fp);
            file_status_[filename].fp = fopen(filename.c_str(), op.c_str());
            // timer_off("stream change           ");
        }
    }
    // timer_off("stream checks             ");
    return file_status_[filename].fp;
}
void DF_Helper::put_tensor(std::string file, double* b, std::pair<size_t, size_t> i0, std::pair<size_t, size_t> i1,
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
void DF_Helper::put_tensor(std::string file, double* Mp, const size_t start1, const size_t stop1, const size_t start2,
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
            error << "DF_Helper:put_tensor: write error";
            throw PSIEXCEPTION(error.str().c_str());
        }
    } else {
        for (size_t i = start1; i < stop1; i++) {
            // write
            size_t s = fwrite(&Mp[i * a1], sizeof(double), a1, fp);
            if (!s) {
                std::stringstream error;
                error << "DF_Helper:put_tensor: write error";
                throw PSIEXCEPTION(error.str().c_str());
            }
            // advance stream
            fseek(fp, st * sizeof(double), SEEK_CUR);
        }
        // manual last one
        size_t s = fwrite(&Mp[(a0 - 1) * a1], sizeof(double), a1, fp);
        if (!s) {
            std::stringstream error;
            error << "DF_Helper:put_tensor: write error";
            throw PSIEXCEPTION(error.str().c_str());
        }
    }
}
void DF_Helper::put_tensor_AO(std::string file, double* Mp, size_t size, size_t start, std::string op) {
    // begin stream
    FILE* fp = stream_check(file, op);

    // adjust position
    fseek(fp, start, SEEK_SET);

    // everything is contiguous
    size_t s = fwrite(&Mp[0], sizeof(double), size, fp);
    if (!s) {
        std::stringstream error;
        error << "DF_Helper:put_tensor_AO: write error";
        throw PSIEXCEPTION(error.str().c_str());
    }
}
void DF_Helper::get_tensor_AO(std::string file, double* Mp, size_t size, size_t start) {
    // begin stream
    FILE* fp = stream_check(file, "rb");

    // adjust position
    fseek(fp, start * sizeof(double), SEEK_SET);

    // everything is contiguous
    size_t s = fread(&Mp[0], sizeof(double), size, fp);
    if (!s) {
        std::stringstream error;
        error << "DF_Helper:get_tensor_AO: read error";
        throw PSIEXCEPTION(error.str().c_str());
    }
}
void DF_Helper::get_tensor_(std::string file, double* b, std::pair<size_t, size_t> i0, std::pair<size_t, size_t> i1,
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
void DF_Helper::get_tensor_(std::string file, double* b, const size_t start1, const size_t stop1, const size_t start2,
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
            error << "DF_Helper:get_tensor: read error";
            throw PSIEXCEPTION(error.str().c_str());
        }
    } else {
        for (size_t i = 0; i < a0 - 1; i++) {
            // read
            size_t s = fread(&b[i * a1], sizeof(double), a1, fp);
            if (!s) {
                std::stringstream error;
                error << "DF_Helper:get_tensor: read error";
                throw PSIEXCEPTION(error.str().c_str());
            }
            // advance stream
            s = fseek(fp, st * sizeof(double), SEEK_CUR);
            if (s) {
                std::stringstream error;
                error << "DF_Helper:get_tensor: read error";
                throw PSIEXCEPTION(error.str().c_str());
            }
        }
        // manual last one
        size_t s = fread(&b[(a0 - 1) * a1], sizeof(double), a1, fp);
        if (!s) {
            std::stringstream error;
            error << "DF_Helper:get_tensor: read error";
            throw PSIEXCEPTION(error.str().c_str());
        }
    }
}
void DF_Helper::compute_AO_Q(const size_t start, const size_t stop, double* Mp,
                             std::vector<std::shared_ptr<TwoBodyAOInt>> eri) {
    size_t begin = Qshell_aggs_[start];
    size_t end = Qshell_aggs_[stop + 1] - 1;
    size_t block_size = end - begin + 1;

    // prepare eri buffers
    size_t nthread = nthreads_;
    if (eri.size() != nthreads_) nthread = eri.size();

    int rank = 0;
    std::vector<const double*> buffer(nthread);
#pragma omp parallel private(rank) num_threads(nthread)
    {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        buffer[rank] = eri[rank]->buffer();
    }

    size_t MU, nummu, NU, numnu, Pshell, numP, mu, omu, nu, onu, P, PHI;
#pragma omp parallel for private(numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, \
                                 rank) schedule(guided) num_threads(nthreads_)
    for (MU = 0; MU < pshells_; MU++) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        nummu = primary_->shell(MU).nfunction();
        for (NU = 0; NU < pshells_; NU++) {
            numnu = primary_->shell(NU).nfunction();
            if (!schwarz_shell_mask_[MU * pshells_ + NU]) {
                continue;
            }
            for (Pshell = start; Pshell <= stop; Pshell++) {
                PHI = aux_->shell(Pshell).function_index();
                numP = aux_->shell(Pshell).nfunction();
                eri[rank]->compute_shell(Pshell, 0, MU, NU);
                for (mu = 0; mu < nummu; mu++) {
                    omu = primary_->shell(MU).function_index() + mu;
                    for (nu = 0; nu < numnu; nu++) {
                        onu = primary_->shell(NU).function_index() + nu;
                        if (!schwarz_fun_mask_[omu * nao_ + onu]) {
                            continue;
                        }
                        for (P = 0; P < numP; P++) {
                            Mp[(big_skips_[omu] * block_size) / naux_ + (PHI + P - begin) * small_skips_[omu] +
                               schwarz_fun_mask_[omu * nao_ + onu] - 1] =
                                buffer[rank][P * nummu * numnu + mu * numnu + nu];
                        }
                    }
                }
            }
        }
    }
}
void DF_Helper::compute_AO_p(const size_t start, const size_t stop, double* Mp,
                             std::vector<std::shared_ptr<TwoBodyAOInt>> eri) {
    size_t begin = pshell_aggs_[start];
    size_t end = pshell_aggs_[stop + 1] - 1;
    size_t block_size = end - begin + 1;
    size_t startind = big_skips_[begin];
    //    outfile->Printf("      MU shell: (%zu, %zu)", start, stop);
    //    outfile->Printf(", nao index: (%zu, %zu), size: %zu\n", begin, end, block_size);

    // prepare eri buffers
    size_t nthread = nthreads_;
    if (eri.size() != nthreads_) nthread = eri.size();

    int rank = 0;
    std::vector<const double*> buffer(nthread);
#pragma omp parallel for private(rank) num_threads(nthread) schedule(static)
    for (size_t i = 0; i < nthread; i++) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        buffer[rank] = eri[rank]->buffer();
    }

    size_t MU, nummu, NU, numnu, Pshell, numP, mu, omu, nu, onu, P, PHI;
#pragma omp parallel for private(numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, \
                                 rank) schedule(guided) num_threads(nthread)
    for (MU = start; MU <= stop; MU++) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        nummu = primary_->shell(MU).nfunction();
        for (NU = 0; NU < pshells_; NU++) {
            numnu = primary_->shell(NU).nfunction();
            if (!schwarz_shell_mask_[MU * pshells_ + NU]) {
                continue;
            }
            for (Pshell = 0; Pshell < Qshells_; Pshell++) {
                PHI = aux_->shell(Pshell).function_index();
                numP = aux_->shell(Pshell).nfunction();
                eri[rank]->compute_shell(Pshell, 0, MU, NU);
                for (mu = 0; mu < nummu; mu++) {
                    omu = primary_->shell(MU).function_index() + mu;
                    for (nu = 0; nu < numnu; nu++) {
                        onu = primary_->shell(NU).function_index() + nu;
                        if (!schwarz_fun_mask_[omu * nao_ + onu]) {
                            continue;
                        }
                        for (P = 0; P < numP; P++) {
                            Mp[big_skips_[omu] - startind + (PHI + P) * small_skips_[omu] +
                               schwarz_fun_mask_[omu * nao_ + onu] - 1] =
                                buffer[rank][P * nummu * numnu + mu * numnu + nu];
                        }
                    }
                }
            }
        }
    }
}
void DF_Helper::compute_AO_p_symm(const size_t start, const size_t stop, double* Mp,
                                  std::vector<std::shared_ptr<TwoBodyAOInt>> eri) {
    size_t begin = pshell_aggs_[start];
    size_t end = pshell_aggs_[stop + 1] - 1;
    size_t block_size = end - begin + 1;
    size_t startind = symm_agg_sizes_[begin] * naux_;
    //    outfile->Printf("      MU shell: (%zu, %zu)", start, stop);
    //    outfile->Printf(", nao index: (%zu, %zu), size: %zu\n", begin, end, block_size);

    // prepare eri buffers
    size_t nthread = nthreads_;
    if (eri.size() != nthreads_) nthread = eri.size();

    int rank = 0;
    std::vector<const double*> buffer(nthread);
#pragma omp parallel for private(rank) num_threads(nthread) schedule(static)
    for (size_t i = 0; i < nthread; i++) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        buffer[rank] = eri[rank]->buffer();
    }

    size_t MU, nummu, NU, numnu, Pshell, numP, mu, omu, nu, onu, P, PHI;
#pragma omp parallel for private(numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, \
                                 rank) schedule(guided) num_threads(nthread)
    for (MU = start; MU <= stop; MU++) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        nummu = primary_->shell(MU).nfunction();
        for (NU = MU; NU < pshells_; NU++) {
            numnu = primary_->shell(NU).nfunction();
            if (!schwarz_shell_mask_[MU * pshells_ + NU]) {
                continue;
            }
            for (Pshell = 0; Pshell < Qshells_; Pshell++) {
                PHI = aux_->shell(Pshell).function_index();
                numP = aux_->shell(Pshell).nfunction();
                eri[rank]->compute_shell(Pshell, 0, MU, NU);
                for (mu = 0; mu < nummu; mu++) {
                    omu = primary_->shell(MU).function_index() + mu;
                    for (nu = 0; nu < numnu; nu++) {
                        onu = primary_->shell(NU).function_index() + nu;
                        if (!schwarz_fun_mask_[omu * nao_ + onu] || omu > onu) {
                            continue;
                        }
                        for (P = 0; P < numP; P++) {
                            size_t jump = schwarz_fun_mask_[omu * nao_ + onu] - schwarz_fun_mask_[omu * nao_ + omu];
                            size_t ind1 = symm_agg_sizes_[omu] * naux_ - startind + (PHI + P) * symm_sizes_[omu] + jump;
                            Mp[ind1] = buffer[rank][P * nummu * numnu + mu * numnu + nu];
                        }
                    }
                }
            }
        }
    }
}
void DF_Helper::grab_AO(const size_t start, const size_t stop, double* Mp) {
    size_t begin = Qshell_aggs_[start];
    size_t end = Qshell_aggs_[stop + 1] - 1;
    size_t block_size = end - begin + 1;
    std::string getf = AO_files_[AO_names_[1]];

    // presumably not thread safe or inherently sequential, but could revisit
    for (size_t i = 0, sta = 0; i < nao_; i++) {
        size_t size = block_size * small_skips_[i];
        size_t jump = begin * small_skips_[i];
        get_tensor_AO(getf, &Mp[sta], size, big_skips_[i] + jump);
        sta += size;
    }
}
void DF_Helper::prepare_metric_core() {
    timer_on("DFH: metric contrsuction");
    std::shared_ptr<FittingMetric> Jinv(new FittingMetric(aux_, true));
    Jinv->form_fitting_metric();
    metrics_[1.0] = Jinv->get_metric();
    timer_off("DFH: metric contrsuction");
}
double* DF_Helper::metric_prep_core(double pow) {
    bool on = false;
    double power;
    for (auto& kv : metrics_) {
        if (!(std::fabs(pow - kv.first) > 1e-13)) {
            on = true;
            power = kv.first;
            break;
        }
    }
    if (!on) {
        power = pow;
        timer_on("DFH: metric power");
        SharedMatrix J = metrics_[1.0];
        J->power(power, condition_);
        metrics_[power] = J;
        timer_off("DFH: metric power");
    }
    return metrics_[power]->pointer()[0];
}
void DF_Helper::prepare_metric() {
    // construct metric
    std::shared_ptr<FittingMetric> Jinv(new FittingMetric(aux_, true));
    Jinv->form_fitting_metric();
    SharedMatrix metric = Jinv->get_metric();
    double* Mp = metric->pointer()[0];

    // create file
    std::string filename = "metric";
    filename.append(".");
    filename.append(std::to_string(1.0));
    filename_maker(filename, naux_, naux_, 1);
    metric_keys_.push_back(std::make_pair(1.0, filename));

    // store
    std::string putf = std::get<0>(files_[filename]);
    put_tensor(putf, Mp, 0, naux_ - 1, 0, naux_ - 1, "wb");
}
std::string DF_Helper::return_metfile(double pow) {
    bool on = 0;
    std::string key;
    for (size_t i = 0; i < metric_keys_.size() && !on; i++) {
        double pos = std::get<0>(metric_keys_[i]);
        if (std::fabs(pos - pow) < 1e-12) {
            key = std::get<1>(metric_keys_[i]);
            on = 1;
        }
    }

    if (!on) key = compute_metric(pow);
    return key;
}
std::string DF_Helper::compute_metric(double pow) {
    // ensure J
    if (std::fabs(pow - 1.0) < 1e-13)
        prepare_metric();
    else {
        // get metric
        SharedMatrix metric(new Matrix("met", naux_, naux_));
        double* metp = metric->pointer()[0];
        std::string filename = return_metfile(1.0);

        // get and compute
        get_tensor_(std::get<0>(files_[filename]), metp, 0, naux_ - 1, 0, naux_ - 1);
        metric->power(pow, condition_);

        // make new file
        std::string name = "metric";
        name.append(".");
        name.append(std::to_string(pow));
        filename_maker(name, naux_, naux_, 1);
        metric_keys_.push_back(std::make_pair(pow, name));

        // store
        std::string putf = std::get<0>(files_[name]);
        put_tensor(putf, metp, 0, naux_ - 1, 0, naux_ - 1, "wb");
    }
    return return_metfile(pow);
}
void DF_Helper::contract_metric(std::string file, double* metp, double* Mp, double* Fp, const size_t tots) {
    std::string getf = std::get<0>(files_[file]);
    std::string putf = std::get<1>(files_[file]);

    size_t a0 = std::get<0>(sizes_[getf]);
    size_t a1 = std::get<1>(sizes_[getf]);
    size_t a2 = std::get<2>(sizes_[getf]);
    size_t count = 0;
    size_t maxim = 0;

    std::string op = "wb";

    // contract in steps
    if (std::get<2>(transf_[file])) {
        // determine blocking
        std::vector<std::pair<size_t, size_t>> steps;
        for (size_t i = 0; i < a0; i++) {
            count++;
            if (tots < count * a1 * a2 || i == a0 - 1) {
                if (count == 1 && i != a0 - 1) {
                    std::stringstream error;
                    error << "DF_Helper:contract_metric: not enough memory.";
                    throw PSIEXCEPTION(error.str().c_str());
                }
                if (tots < count * a1 * a2) {
                    maxim = (maxim < count ? count - 1 : maxim);
                    steps.push_back(std::make_pair(i - count + 1, i - 1));
                    i--;
                } else {
                    maxim = (maxim < count ? count : maxim);
                    steps.push_back(std::make_pair(i - count + 1, i));
                }
                count = 0;
            }
        }

        for (size_t i = 0; i < steps.size(); i++) {
            size_t begin = std::get<0>(steps[i]);
            size_t end = std::get<1>(steps[i]);
            size_t bs = end - begin + 1;

            get_tensor_(getf, Mp, begin, end, 0, a1 * a2 - 1);
            timer_on("DFH: Total Workflow");
            C_DGEMM('N', 'N', bs * a1, a2, a2, 1.0, Mp, a2, metp, a2, 0.0, Fp, a2);
            timer_off("DFH: Total Workflow");
            put_tensor(putf, Fp, begin, end, 0, a1 * a2 - 1, op);
        }
    } else {
        // determine blocking
        std::vector<std::pair<size_t, size_t>> steps;
        for (size_t i = 0; i < a1; i++) {
            count++;
            if (tots < count * a0 * a2 || i == a1 - 1) {
                if (count == 1 && i != a1 - 1) {
                    std::stringstream error;
                    error << "DF_Helper:contract_metric: not enough memory.";
                    throw PSIEXCEPTION(error.str().c_str());
                }
                if (tots < count * a0 * a2) {
                    maxim = (maxim < count ? count - 1 : maxim);
                    steps.push_back(std::make_pair(i - count + 1, i - 1));
                    i--;
                } else {
                    maxim = (maxim < count ? count : maxim);
                    steps.push_back(std::make_pair(i - count + 1, i));
                }
                count = 0;
            }
        }

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

void DF_Helper::contract_metric_AO_core(double* Qpq, double* metp) {
// loop and contract
#pragma omp parallel for num_threads(nthreads_) schedule(guided)
    for (size_t j = 0; j < nao_; j++) {
        size_t mi = small_skips_[j];
        size_t skips = big_skips_[j];
        C_DGEMM('N', 'N', naux_, mi, naux_, 1.0, metp, naux_, &Qpq[skips], mi, 0.0, &Ppq_[skips], mi);
    }
}
void DF_Helper::contract_metric_AO_core_symm(double* Qpq, double* metp, size_t begin, size_t end) {
    // loop and contract
    size_t startind = symm_agg_sizes_[begin] * naux_;
#pragma omp parallel for num_threads(nthreads_) schedule(guided)
    for (size_t j = begin; j <= end; j++) {
        size_t mi = symm_sizes_[j];
        size_t si = small_skips_[j];
        size_t jump = symm_skips_[j];
        size_t skip1 = big_skips_[j];
        size_t skip2 = symm_agg_sizes_[j] * naux_ - startind;
        C_DGEMM('N', 'N', naux_, mi, naux_, 1.0, metp, naux_, &Qpq[skip2], mi, 0.0, &Ppq_[skip1 + jump], si);
    }
    // copy upper-to-lower
    double* Ppq = Ppq_.data();
#pragma omp parallel for num_threads(nthreads_) schedule(static, nao_ / nthreads_)
    for (size_t omu = begin; omu <= end; omu++) {
        for (size_t Q = 0; Q < naux_; Q++) {
            for (size_t onu = omu + 1; onu < nao_; onu++) {
                if (schwarz_fun_mask_[omu * nao_ + onu]) {
                    size_t ind1 = big_skips_[onu] + Q * small_skips_[onu] + schwarz_fun_mask_[onu * nao_ + omu] - 1;
                    size_t ind2 = big_skips_[omu] + Q * small_skips_[omu] + schwarz_fun_mask_[omu * nao_ + onu] - 1;
                    Ppq[ind1] = Ppq[ind2];
                }
            }
        }
    }
}
void DF_Helper::add_space(std::string key, SharedMatrix M) {
    size_t a0 = M->rowspi()[0];
    size_t a1 = M->colspi()[0];

    if (!built) {
        throw PSIEXCEPTION("DF_Helper:add_space: call initialize() before adding spaces!");
    } else if (a0 != nao_) {
        std::stringstream error;
        error << "DF_Helper:add_space: illegal space (" << key << "), primary axis is not nao";
        throw PSIEXCEPTION(error.str().c_str());
    } else if (spaces_.find(key) != spaces_.end()) {
        if (a1 != std::get<1>(spaces_[key])) {
            std::stringstream error;
            error << "DF_Helper:add_space: illegal space (" << key << "), new space has incorrect dimension!";
            throw PSIEXCEPTION(error.str().c_str());
        }
    }
    //    else if(wMO_ < a1 && !direct_){
    //            std::stringstream error;
    //            error << "DF_Helper:add_space: illegal space ("<< key <<"), new space is larger than the "   <<
    //            "worst MO size -(" << wMO_ << "<" << a1 << ")- specified at the time when initialize() was " <<
    //            "called, use set_MO_hint() before calling initialize()";
    //            throw PSIEXCEPTION(error.str().c_str());
    //    }
    //    else
    sorted_spaces_.push_back(std::make_pair(key, a1));

    spaces_[key] = std::make_tuple(M, a1);
}
void DF_Helper::add_transformation(std::string name, std::string key1, std::string key2, std::string order) {
    if (spaces_.find(key1) == spaces_.end()) {
        std::stringstream error;
        error << "DF_Helper:add_transformation: first space (" << key1 << "), is not in space list!";
        throw PSIEXCEPTION(error.str().c_str());
    } else if (spaces_.find(key2) == spaces_.end()) {
        std::stringstream error;
        error << "DF_Helper:add_transformation: second space (" << key2 << "), is not in space list!";
        throw PSIEXCEPTION(error.str().c_str());
    }

    int op;
    if (!order.compare("Qpq"))
        op = 0;
    else if (!order.compare("pqQ"))
        op = 1;
    else {
        throw PSIEXCEPTION("Matt doesnt do exceptions.");
    }

    size_t a1 = std::get<1>(spaces_[key1]);
    size_t a2 = std::get<1>(spaces_[key2]);
    if (op)
        filename_maker(name, naux_, a1, a2, 1);
    else
        filename_maker(name, naux_, a1, a2);

    transf_[name] = std::make_tuple(key1, key2, op);
}
void DF_Helper::clear_spaces() {
    // clear spaces
    spaces_.clear();
    sorted_spaces_.clear();
    order_.clear();
    bspace_.clear();
    strides_.clear();

    // no ordering
    ordered_ = false;
}
void DF_Helper::clear_all() {
    // outfile->Printf("\n clearing everything (spaces, integrals)! \n");

    clear_spaces();

    // clear transformations
    std::string left;
    std::string right;

    for (auto& kv : transf_) {
        left = std::get<0>(files_[kv.first]);
        right = std::get<1>(files_[kv.first]);

        // close streams
        if (!MO_core_) {
            if (direct_) fclose(file_status_[left].fp);
            fclose(file_status_[right].fp);
        }

        // delete stream holders
        if (!MO_core_) {
            if (direct_) file_status_.erase(left);
            file_status_.erase(right);
        }

        // delete file holders, sizes
        sizes_.erase(std::get<0>(files_[kv.first]));
        sizes_.erase(std::get<1>(files_[kv.first]));
        files_.erase(kv.first);
    }

    // tranposes too
    tsizes_.clear();

    // what transformations?
    if (MO_core_) transf_core_.clear();
    transf_.clear();
    transformed_ = false;
}
std::pair<size_t, size_t> DF_Helper::identify_order() {
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
        std::list<std::string>::iterator itr, end;
        for (itr = needs.begin(), end = needs.end(); itr != end; ++itr) {
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
                needs.erase(itr);
                itr--;
            }
        }
        if (st > 0) {
            strides_.push_back(st);
        }
    }
    // print_order();
    ordered_ = true;
    return std::make_pair(largest, maximum);
}
void DF_Helper::print_order() {
    size_t o = order_.size();
    size_t b = bspace_.size();
    outfile->Printf("\n     ==> DF_Helper:--Begin Transformations Information <==\n\n");
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
    outfile->Printf("\n\n     ==> DF_Helper:--End Transformations Information <==\n\n");
}
void DF_Helper::transform() {
    timer_on("DFH: transform()");
    (MO_core_ ? transform_core() : transform_disk());
    timer_off("DFH: transform()");
    transformed_ = true;
}
void DF_Helper::transform_core() {
    // outfile->Printf("\n     ==> DF_Helper:--Begin Transformations (core)<==\n\n");

    size_t nthreads = nthreads_;
    size_t naux = naux_;
    size_t nao = nao_;
    int rank = 0;

    // reset tranposes
    tsizes_.clear();

    // gather transformation info
    if (!ordered_) info_ = identify_order();
    size_t wtmp = std::get<0>(info_);
    size_t wfinal = std::get<1>(info_);

    // prep stream
    if (!direct_ && !AO_core_) stream_check(AO_files_[AO_names_[1]], "rb");

    // blocking
    std::vector<std::pair<size_t, size_t>> Qsteps;
    std::pair<size_t, size_t> Qlargest = Qshell_blocks_for_transform(memory_, wtmp, wfinal, Qsteps);
    size_t max_block = std::get<1>(Qlargest);

    // prepare eri buffers
    std::vector<std::vector<double>> C_buffers(nthreads);
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(aux_, zero, primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthreads);
#pragma omp parallel private(rank) num_threads(nthreads)
    {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        std::vector<double> Cp(nao * wtmp);
        C_buffers[rank] = Cp;
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    }

    // stripe transformed integrals
    for (auto& kv : transf_) {
        size_t size = std::get<1>(spaces_[std::get<0>(kv.second)]) * std::get<1>(spaces_[std::get<1>(kv.second)]);
        transf_core_[kv.first].reserve(size * naux);
    }

    // scoped buffer declarations
    {
        // declare bufs
        std::vector<double> T;
        T.reserve(max_block * nao * wtmp);
        std::vector<double> F;
        F.reserve(max_block * wfinal);
        std::vector<double> N;
        N.reserve(max_block * wfinal);
        double* Tp = T.data();
        double* Fp = F.data();
        double* Np = N.data();

        // declare bufs
        std::vector<double> M;  // AOs
        double* Mp;
        if (!AO_core_) {
            M.reserve(std::get<0>(Qlargest));
            Mp = M.data();
        } else
            Mp = Ppq_.data();

        // transform in steps
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
            timer_on("DFH: Grabbing AOs");
            if (AO_core_)
                Mp = Ppq_.data();
            else if (direct_)
                compute_AO_Q(start, stop, Mp, eri);
            else
                grab_AO(start, stop, Mp);
            timer_off("DFH: Grabbing AOs");

            // stride through best spaces
            for (size_t i = 0, count = 0; i < bspace_.size(); count += strides_[i], i++) {
                // grab best space
                std::string bspace = bspace_[i];
                size_t bsize = std::get<1>(spaces_[bspace]);
                double* Bpt = std::get<0>(spaces_[bspace])->pointer()[0];

                // form tmp, thread over spM (nao)
                timer_on("DFH: First Contraction");
#pragma omp parallel for firstprivate(nao, naux, bsize, block_size) private(rank) schedule(guided) num_threads(nthreads)
                for (size_t k = 0; k < nao_; k++) {
                    // truncate transformation matrix according to fun_mask
                    size_t sp_size = small_skips_[k];
                    size_t jump = (AO_core_ ? big_skips_[k] + bcount * sp_size : (big_skips_[k] * block_size) / naux_);

#ifdef _OPENMP
                    rank = omp_get_thread_num();
#endif
                    for (size_t m = 0, sp_count = -1; m < nao_; m++) {
                        if (schwarz_fun_mask_[k * nao_ + m]) {
                            sp_count++;
                            C_DCOPY(bsize, &Bpt[m * bsize], 1, &C_buffers[rank][sp_count * bsize], 1);
                        }
                    }

                    // (Qm)(mb)->(Qb)
                    C_DGEMM('N', 'N', block_size, bsize, sp_size, 1.0, &Mp[jump], sp_size, &C_buffers[rank][0], bsize,
                            0.0, &Tp[k * block_size * bsize], bsize);
                }
                timer_off("DFH: First Contraction");

                // to completion per transformation
                for (size_t k = 0; k < strides_[i]; k++) {
                    std::string left = std::get<0>(transf_[order_[count + k]]);
                    std::string right = std::get<1>(transf_[order_[count + k]]);
                    std::tuple<SharedMatrix, size_t> I = (bspace.compare(left) == 0 ? spaces_[right] : spaces_[left]);

                    size_t wsize = std::get<1>(I);
                    double* Wp = std::get<0>(I)->pointer()[0];

                    // (wp)(p|Qb)->(w|Qb)
                    timer_on("DFH: Second Contraction");
                    C_DGEMM('T', 'N', wsize, block_size * bsize, nao_, 1.0, Wp, wsize, Tp, block_size * bsize, 0.0, Fp,
                            block_size * bsize);
                    timer_off("DFH: Second Contraction");

                    // Transpose in memory for convenient formats
                    size_t st1 = bsize * wsize;
                    size_t st2 = bsize * block_size;
                    size_t st3 = wsize * block_size;

                    // final in-core MO pointer
                    double* Lp = transf_core_[order_[count + k]].data();

                    timer_on("DFH: (w|Qb)->(bw|Q)");
                    if (bspace.compare(left) == 0) {  // (w|Qb)->(Q|bw)

                        if (std::get<2>(transf_[order_[count + k]])) {  // (w|Qb)->(bw|Q)

#pragma omp parallel for num_threads(nthreads_)
                            for (size_t x = 0; x < block_size; x++) {
                                for (size_t y = 0; y < bsize; y++) {
                                    for (size_t z = 0; z < wsize; z++) {
                                        Lp[y * wsize * naux + z * naux + (bcount + x)] = Fp[z * st2 + x * bsize + y];
                                    }
                                }
                            }

                        } else {  // (w|Qb)->(Q|bw)

#pragma omp parallel for num_threads(nthreads_)
                            for (size_t x = 0; x < block_size; x++) {
                                for (size_t y = 0; y < bsize; y++) {
                                    for (size_t z = 0; z < wsize; z++) {
                                        Lp[(bcount + x) * st1 + y * wsize + z] = Fp[z * st2 + x * bsize + y];
                                    }
                                }
                            }
                        }

                    } else {
                        if (std::get<2>(transf_[order_[count + k]])) {  // // (w|Qb)->(wbQ)

#pragma omp parallel for num_threads(nthreads_)
                            for (size_t x = 0; x < block_size; x++) {
                                for (size_t z = 0; z < wsize; z++) {
                                    for (size_t y = 0; y < bsize; y++) {
                                        Lp[z * naux * bsize + y * naux + (bcount + x)] = Fp[z * st2 + x * bsize + y];
                                    }
                                }
                            }

                        } else {  // (w|Qb)->(Q|wb)

#pragma omp parallel for num_threads(nthreads_)
                            for (size_t x = 0; x < block_size; x++) {
                                for (size_t z = 0; z < wsize; z++) {
#pragma omp simd
                                    for (size_t y = 0; y < bsize; y++) {
                                        Lp[(bcount + x) * st1 + z * bsize + y] = Fp[z * st2 + x * bsize + y];
                                    }
                                }
                            }
                        }
                    }
                    timer_off("DFH: (w|Qb)->(bw|Q)");
                }
            }
        }
    }  // buffers destroyed with std housekeeping

    if (direct_) {
        double* metp;
        std::vector<double> metric;
        if (!hold_met_) {
            metric.reserve(naux_ * naux_);
            metp = metric.data();
            std::string filename = return_metfile(mpower_);
            get_tensor_(std::get<0>(files_[filename]), metp, 0, naux_ - 1, 0, naux_ - 1);
        } else
            metp = metric_prep_core(mpower_);

        std::vector<double> N;
        N.reserve(naux * wfinal);
        double* Np = N.data();

        for (auto& kv : transf_core_) {
            size_t a0 = std::get<0>(sizes_[std::get<1>(files_[kv.first])]);
            size_t a1 = std::get<1>(sizes_[std::get<1>(files_[kv.first])]);
            size_t a2 = std::get<2>(sizes_[std::get<1>(files_[kv.first])]);

            double* Lp = kv.second.data();
            C_DCOPY(a0 * a1 * a2, Lp, 1, Np, 1);

            if (std::get<2>(transf_[kv.first]))
                C_DGEMM('N', 'N', a0 * a1, a2, a2, 1.0, Np, a2, metp, a2, 0.0, Lp, a2);
            else
                C_DGEMM('N', 'N', a0, a1 * a2, a0, 1.0, metp, naux, Np, a1 * a2, 0.0, Lp, a1 * a2);
        }
    }
    // outfile->Printf("\n     ==> DF_Helper:--End Transformations (core)<==\n\n");
}
void DF_Helper::transform_disk() {
    // outfile->Printf("\n     ==> DF_Helper:--Begin Transformations (disk)<==\n\n");

    size_t nthreads = nthreads_;
    size_t naux = naux_;
    size_t nao = nao_;
    int rank = 0;

    // reset tranposes
    tsizes_.clear();

    // gather transformation info
    if (!ordered_) info_ = identify_order();
    size_t wtmp = std::get<0>(info_);
    size_t wfinal = std::get<1>(info_);

    // prep stream
    if (!direct_ && !AO_core_) stream_check(AO_files_[AO_names_[1]], "rb");

    // blocking
    std::vector<std::pair<size_t, size_t>> Qsteps;
    std::pair<size_t, size_t> Qlargest = Qshell_blocks_for_transform(memory_, wtmp, wfinal, Qsteps);
    size_t max_block = std::get<1>(Qlargest);

    // prepare eri buffers
    size_t nthread = nthreads_;
    std::vector<std::vector<double>> C_buffers(nthreads_);
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(aux_, zero, primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthread);
#pragma omp parallel private(rank) num_threads(nthreads_)
    {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        std::vector<double> Cp(nao * wtmp);
        C_buffers[rank] = Cp;
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    }

    // scope buffer declarations
    {
        // declare bufs
        std::vector<double> T;
        T.reserve(max_block * nao * wtmp);
        std::vector<double> F;
        F.reserve(max_block * wfinal);
        std::vector<double> N;
        N.reserve(max_block * wfinal);
        double* Tp = T.data();
        double* Fp = F.data();
        double* Np = N.data();

        // declare bufs
        std::vector<double> M;  // AOs
        double* Mp;
        if (!AO_core_) {
            M.reserve(std::get<0>(Qlargest));
            Mp = M.data();
        } else
            Mp = Ppq_.data();

        // transform in steps
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
            timer_on("DFH: Grabbing AOs");
            if (AO_core_)
                Mp = Ppq_.data();
            else if (direct_) {
                timer_on("DFH: Total Workflow");
                compute_AO_Q(start, stop, Mp, eri);
                timer_off("DFH: Total Workflow");
            } else
                grab_AO(start, stop, Mp);
            timer_off("DFH: Grabbing AOs");

            // stride through best spaces
            for (size_t i = 0, count = 0; i < bspace_.size(); count += strides_[i], i++) {
                // grab best space
                std::string bspace = bspace_[i];
                size_t bsize = std::get<1>(spaces_[bspace]);
                double* Bpt = std::get<0>(spaces_[bspace])->pointer()[0];

                timer_on("DFH: Total Workflow");
                timer_on("DFH: Total Transform");
                timer_on("DFH: First Contraction");
// form temp, thread over spM (nao)
#pragma omp parallel for firstprivate(nao, naux, bsize, \
                                      block_size) private(rank) schedule(guided) num_threads(nthreads_)
                for (size_t k = 0; k < nao_; k++) {
                    // truncate transformation matrix according to fun_mask
                    size_t sp_size = small_skips_[k];
                    size_t jump = (AO_core_ ? big_skips_[k] + bcount * sp_size : (big_skips_[k] * block_size) / naux_);

#ifdef _OPENMP
                    rank = omp_get_thread_num();
#endif
                    for (size_t m = 0, sp_count = -1; m < nao_; m++) {
                        if (schwarz_fun_mask_[k * nao_ + m]) {
                            sp_count++;
                            C_DCOPY(bsize, &Bpt[m * bsize], 1, &C_buffers[rank][sp_count * bsize], 1);
                        }
                    }

                    // (Qm)(mb)->(Qb)
                    C_DGEMM('N', 'N', block_size, bsize, sp_size, 1.0, &Mp[jump], sp_size, &C_buffers[rank][0], bsize,
                            0.0, &Tp[k * block_size * bsize], bsize);
                }
                timer_off("DFH: First Contraction");
                timer_off("DFH: Total Transform");
                timer_off("DFH: Total Workflow");

                // to completion per transformation
                for (size_t k = 0; k < strides_[i]; k++) {
                    std::string left = std::get<0>(transf_[order_[count + k]]);
                    std::string right = std::get<1>(transf_[order_[count + k]]);
                    std::tuple<SharedMatrix, size_t> I = (bspace.compare(left) == 0 ? spaces_[right] : spaces_[left]);

                    size_t wsize = std::get<1>(I);
                    double* Wp = std::get<0>(I)->pointer()[0];

                    // (wp)(p|Qb)->(w|Qb)
                    timer_on("DFH: Total Workflow");
                    timer_on("DFH: Total Transform");
                    timer_on("DFH: Second Contraction");
                    C_DGEMM('T', 'N', wsize, block_size * bsize, nao_, 1.0, Wp, wsize, Tp, block_size * bsize, 0.0, Fp,
                            block_size * bsize);
                    timer_off("DFH: Second Contraction");
                    timer_off("DFH: Total Transform");
                    timer_off("DFH: Total Workflow");

                    // setup putf
                    std::string putf =
                        (!direct_ ? std::get<1>(files_[order_[count + k]]) : std::get<0>(files_[order_[count + k]]));

                    std::string op = "ab";

                    // Transpose in memory for convenient formats
                    size_t st1 = bsize * wsize;
                    size_t st2 = bsize * block_size;
                    size_t st3 = wsize * block_size;

                    timer_on("DFH: (w|Qb)->(bw|Q)");
                    if (bspace.compare(left) == 0) {
                        if (std::get<2>(transf_[order_[count + k]])) {  // (w|Qb)->(bw|Q)

#pragma omp parallel for num_threads(nthreads_)
                            for (size_t x = 0; x < block_size; x++) {
                                for (size_t y = 0; y < bsize; y++) {
                                    for (size_t z = 0; z < wsize; z++) {
                                        Np[y * st3 + z * block_size + x] = Fp[z * st2 + x * bsize + y];
                                    }
                                }
                            }
                            timer_off("DFH: (w|Qb)->(bw|Q)");
                            timer_on("DFH: MO to disk");
                            put_tensor(putf, Np, std::make_pair(0, bsize - 1), std::make_pair(0, wsize - 1),
                                       std::make_pair(begin, end), "wb");

                        } else {  // (w|Qb)->(Q|bw)

#pragma omp parallel for num_threads(nthreads_)
                            for (size_t x = 0; x < block_size; x++) {
                                for (size_t y = 0; y < bsize; y++) {
                                    for (size_t z = 0; z < wsize; z++) {
                                        Np[x * st1 + y * wsize + z] = Fp[z * st2 + x * bsize + y];
                                    }
                                }
                            }
                            timer_off("DFH: (w|Qb)->(bw|Q)");
                            timer_on("DFH: MO to disk");
                            put_tensor(putf, Np, std::make_pair(begin, end), std::make_pair(0, bsize - 1),
                                       std::make_pair(0, wsize - 1), op);
                        }

                    } else {
                        if (std::get<2>(transf_[order_[count + k]])) {  // // (w|Qb)->(wbQ)

#pragma omp parallel for num_threads(nthreads_)
                            for (size_t x = 0; x < block_size; x++) {
                                for (size_t z = 0; z < wsize; z++) {
                                    for (size_t y = 0; y < bsize; y++) {
                                        Np[z * st2 + y * block_size + x] = Fp[z * st2 + x * bsize + y];
                                    }
                                }
                            }
                            timer_off("DFH: (w|Qb)->(bw|Q)");
                            timer_on("DFH: MO to disk");
                            put_tensor(putf, Np, std::make_pair(0, wsize - 1), std::make_pair(0, bsize - 1),
                                       std::make_pair(begin, end), "wb");

                        } else {  // (w|Qb)->(Q|wb)

#pragma omp parallel for num_threads(nthreads_)
                            for (size_t x = 0; x < block_size; x++) {
                                for (size_t z = 0; z < wsize; z++) {
#pragma omp simd
                                    for (size_t y = 0; y < bsize; y++) {
                                        Np[x * st1 + z * bsize + y] = Fp[z * st2 + x * bsize + y];
                                    }
                                }
                            }
                            timer_off("DFH: (w|Qb)->(bw|Q)");
                            timer_on("DFH: MO to disk");
                            put_tensor(putf, Np, std::make_pair(begin, end), std::make_pair(0, wsize - 1),
                                       std::make_pair(0, bsize - 1), op);
                        }
                    }
                    timer_off("DFH: MO to disk");
                }
            }
        }
    }
    // outfile->Printf("\n     ==> DF_Helper:--End Transformations (disk)<==\n\n");

    if (direct_) {
        // total size allowed, in doubles
        size_t total_mem =
            (memory_ > wfinal * naux * 2 + naux_ * naux_ ? wfinal * naux : (memory_ - naux_ * naux_) / 2);

        std::vector<double> M;
        std::vector<double> F;
        M.reserve(total_mem);
        F.reserve(total_mem);
        double* Mp = M.data();
        double* Fp = F.data();
        double* metp;

        if (hold_met_) {
            metp = metric_prep_core(mpower_);
            for (std::vector<std::string>::iterator itr = order_.begin(); itr != order_.end(); itr++)
                contract_metric(*itr, metp, Mp, Fp, total_mem);

        } else {
            std::vector<double> metric;
            metric.reserve(naux * naux);
            metp = metric.data();

            std::string mfilename = return_metfile(mpower_);
            get_tensor_(std::get<0>(files_[mfilename]), metp, 0, naux_ - 1, 0, naux_ - 1);

            for (std::vector<std::string>::iterator itr = order_.begin(); itr != order_.end(); itr++)
                contract_metric(*itr, metp, Mp, Fp, total_mem);
        }
    }
}

// Fill using a pointer, be cautious of bounds!!
void DF_Helper::fill_tensor(std::string name, double* b) {
    check_file_key(name);
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    fill_tensor(name, b, {0, std::get<0>(sizes)}, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DF_Helper::fill_tensor(std::string name, double* b, std::vector<size_t> a1) {
    check_file_key(name);
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    fill_tensor(name, b, a1, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DF_Helper::fill_tensor(std::string name, double* b, std::vector<size_t> a1, std::vector<size_t> a2) {
    check_file_key(name);
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    fill_tensor(name, b, a1, a2, {0, std::get<2>(sizes)});
}
void DF_Helper::fill_tensor(std::string name, double* b, std::vector<size_t> a1, std::vector<size_t> a2,
                            std::vector<size_t> a3) {
    if (a1.size() != 2) {
        std::stringstream error;
        error << "DF_Helper:fill_tensor:  axis 0 tensor indexing vector has " << a1.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (a2.size() != 2) {
        std::stringstream error;
        error << "DF_Helper:fill_tensor:  axis 1 tensor indexing vector has " << a2.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (a3.size() != 2) {
        std::stringstream error;
        error << "DF_Helper:fill_tensor:  axis 2 tensor indexing vector has " << a3.size() << " elements!";
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
void DF_Helper::fill_tensor(std::string name, SharedMatrix M) {
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    fill_tensor(name, M, {0, std::get<0>(sizes)}, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DF_Helper::fill_tensor(std::string name, SharedMatrix M, std::vector<size_t> a1) {
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    fill_tensor(name, M, a1, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DF_Helper::fill_tensor(std::string name, SharedMatrix M, std::vector<size_t> a1, std::vector<size_t> a2) {
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    fill_tensor(name, M, a1, a2, {0, std::get<2>(sizes)});
}
void DF_Helper::fill_tensor(std::string name, SharedMatrix M, std::vector<size_t> t0, std::vector<size_t> t1,
                            std::vector<size_t> t2) {
    std::string filename = std::get<1>(files_[name]);
    // has this integral been transposed?
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    if (t0.size() != 2) {
        std::stringstream error;
        error << "DF_Helper:fill_tensor:  axis 0 tensor indexing vector has " << t0.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (t1.size() != 2) {
        std::stringstream error;
        error << "DF_Helper:fill_tensor:  axis 1 tensor indexing vector has " << t1.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (t2.size() != 2) {
        std::stringstream error;
        error << "DF_Helper:fill_tensor:  axis 2 tensor indexing vector has " << t2.size() << " elements!";
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

        double* Fp = transf_core_[name].data();
#pragma omp parallel num_threads(nthreads_)
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
SharedMatrix DF_Helper::get_tensor(std::string name) {
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    return get_tensor(name, {0, std::get<0>(sizes)}, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
SharedMatrix DF_Helper::get_tensor(std::string name, std::vector<size_t> a1) {
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    return get_tensor(name, a1, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
SharedMatrix DF_Helper::get_tensor(std::string name, std::vector<size_t> a1, std::vector<size_t> a2) {
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    return get_tensor(name, a1, a2, {0, std::get<2>(sizes)});
}
SharedMatrix DF_Helper::get_tensor(std::string name, std::vector<size_t> t0, std::vector<size_t> t1,
                                   std::vector<size_t> t2) {
    // has this integral been transposed?
    std::string filename = std::get<1>(files_[name]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);

    if (t0.size() != 2) {
        std::stringstream error;
        error << "DF_Helper:fill_tensor:  axis 0 tensor indexing vector has " << t0.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (t1.size() != 2) {
        std::stringstream error;
        error << "DF_Helper:fill_tensor:  axis 1 tensor indexing vector has " << t1.size() << " elements!";
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (t2.size() != 2) {
        std::stringstream error;
        error << "DF_Helper:fill_tensor:  axis 2 tensor indexing vector has " << t2.size() << " elements!";
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

    SharedMatrix M(new Matrix("M", A0, A1 * A2));
    double* Mp = M->pointer()[0];

    if (MO_core_) {
        size_t a0 = std::get<0>(sizes);
        size_t a1 = std::get<1>(sizes);
        size_t a2 = std::get<2>(sizes);

        double* Fp = transf_core_[name].data();
#pragma omp parallel num_threads(nthreads_)
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
void DF_Helper::add_disk_tensor(std::string key, std::tuple<size_t, size_t, size_t> dimensions) {
    if (files_.count(key)) {
        std::stringstream error;
        error << "DF_Helper:add_disk_tensor:  tensor already exists: (" << key << "!";
        throw PSIEXCEPTION(error.str().c_str());
    }

    filename_maker(key, std::get<0>(dimensions), std::get<1>(dimensions), std::get<2>(dimensions));
}

// Write to a disk tensor from Sharedmatrix
void DF_Helper::write_disk_tensor(std::string key, SharedMatrix M) {
    check_file_key(key);
    std::string filename = std::get<1>(files_[key]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);
    write_disk_tensor(key, M, {0, std::get<0>(sizes)}, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DF_Helper::write_disk_tensor(std::string key, SharedMatrix M, std::vector<size_t> a1) {
    check_file_key(key);
    std::string filename = std::get<1>(files_[key]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);
    write_disk_tensor(key, M, a1, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DF_Helper::write_disk_tensor(std::string key, SharedMatrix M, std::vector<size_t> a1, std::vector<size_t> a2) {
    check_file_key(key);
    std::string filename = std::get<1>(files_[key]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);
    write_disk_tensor(key, M, a1, a2, {0, std::get<2>(sizes)});
}
void DF_Helper::write_disk_tensor(std::string key, SharedMatrix M, std::vector<size_t> a0, std::vector<size_t> a1,
                                  std::vector<size_t> a2) {
    // being pythonic ;)
    std::pair<size_t, size_t> i0 = std::make_pair(a0[0], a0[1] - 1);
    std::pair<size_t, size_t> i1 = std::make_pair(a1[0], a1[1] - 1);
    std::pair<size_t, size_t> i2 = std::make_pair(a2[0], a2[1] - 1);

    check_file_key(key);
    check_file_tuple(key, i0, i1, i2);
    check_matrix_size(key, M, i0, i1, i2);

    // you write over transformed integrals or you write new disk tensors
    std::string op = (transf_.count(key) ? "r+b" : "ab");
    put_tensor(std::get<1>(files_[key]), M->pointer()[0], i0, i1, i2, op);
}

// Write to a disk tensor from pointer, be careful!
void DF_Helper::write_disk_tensor(std::string key, double* b) {
    check_file_key(key);
    std::string filename = std::get<1>(files_[key]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);
    write_disk_tensor(key, b, {0, std::get<0>(sizes)}, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DF_Helper::write_disk_tensor(std::string key, double* b, std::vector<size_t> a0) {
    check_file_key(key);
    std::string filename = std::get<1>(files_[key]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);
    write_disk_tensor(key, b, a0, {0, std::get<1>(sizes)}, {0, std::get<2>(sizes)});
}
void DF_Helper::write_disk_tensor(std::string key, double* b, std::vector<size_t> a0, std::vector<size_t> a1) {
    check_file_key(key);
    std::string filename = std::get<1>(files_[key]);
    std::tuple<size_t, size_t, size_t> sizes;
    sizes = (tsizes_.find(filename) != tsizes_.end() ? tsizes_[filename] : sizes_[filename]);
    write_disk_tensor(key, b, a0, a1, {0, std::get<2>(sizes)});
}
void DF_Helper::write_disk_tensor(std::string key, double* b, std::vector<size_t> a0, std::vector<size_t> a1,
                                  std::vector<size_t> a2) {
    // being pythonic ;)
    std::pair<size_t, size_t> i0 = std::make_pair(a0[0], a0[1] - 1);
    std::pair<size_t, size_t> i1 = std::make_pair(a1[0], a1[1] - 1);
    std::pair<size_t, size_t> i2 = std::make_pair(a2[0], a2[1] - 1);

    check_file_key(key);
    check_file_tuple(key, i0, i1, i2);

    // you write over transformed integrals or you write new disk tensors
    std::string op = (transf_.count(key) ? "r+b" : "wb");
    put_tensor(std::get<1>(files_[key]), b, i0, i1, i2, op);
}

void DF_Helper::check_file_key(std::string name) {
    if (files_.find(name) == files_.end()) {
        std::stringstream error;
        error << "DF_Helper:get_tensor OR write_tensor: " << name << " not found.";
        throw PSIEXCEPTION(error.str().c_str());
    }
}
void DF_Helper::check_matrix_size(std::string name, SharedMatrix M, std::pair<size_t, size_t> t0,
                                  std::pair<size_t, size_t> t1, std::pair<size_t, size_t> t2) {
    size_t A0 = std::get<1>(t0) - std::get<0>(t0) + 1;
    size_t A1 = (std::get<1>(t1) - std::get<0>(t1) + 1) * (std::get<1>(t2) - std::get<0>(t2) + 1);

    size_t a0 = M->rowspi()[0];
    size_t a1 = M->colspi()[0];

    if (A0 * A1 > a0 * a1) {
        std::stringstream error;
        error << "DF_Helper:get_tensor: your matrix contridicts your tuple sizes when obtaining the (" << name
              << ") integral.  ";
        error << "you gave me a matrix of size: (" << a0 << "," << a1 << "), but tuple sizes give:(" << A0 << "," << A1
              << ")";
        throw PSIEXCEPTION(error.str().c_str());
    }
}
void DF_Helper::check_file_tuple(std::string name, std::pair<size_t, size_t> t0, std::pair<size_t, size_t> t1,
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
void DF_Helper::transpose(std::string name, std::tuple<size_t, size_t, size_t> order) {
    if (!files_.count(name)) {
        std::stringstream error;
        error << "DF_Helper::transpose(): cannot transpose input (" << name << "), name doe not exist!";
        throw PSIEXCEPTION(error.str().c_str());
    }

    (MO_core_ ? transpose_core(name, order) : transpose_disk(name, order));
}
void DF_Helper::transpose_core(std::string name, std::tuple<size_t, size_t, size_t> order) {
    size_t a0 = std::get<0>(order);
    size_t a1 = std::get<1>(order);
    size_t a2 = std::get<2>(order);

    std::string filename = std::get<1>(files_[name]);
    size_t M0 = std::get<0>(sizes_[filename]);
    size_t M1 = std::get<1>(sizes_[filename]);
    size_t M2 = std::get<2>(sizes_[filename]);
    std::tuple<size_t, size_t, size_t> sizes;

    std::vector<double> M;
    M.reserve(M0 * M1 * M2);
    double* Mp = M.data();
    double* Fp = transf_core_[name].data();
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
#pragma omp parallel num_threads(nthreads_)
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
#pragma omp parallel num_threads(nthreads_)
            for (size_t i = 0; i < M0; i++) {
                for (size_t j = 0; j < M1; j++) {
#pragma omp simd
                    for (size_t k = 0; k < M2; k++) {
                        Fp[j * M0 * M2 + i * M2 + k] = Mp[i * M1 * M2 + j * M2 + k];
                    }
                }
            }
        } else if (a1 == 2) {  // (0|12) -> (1|20)
#pragma omp parallel num_threads(nthreads_)
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
#pragma omp parallel num_threads(nthreads_)
            for (size_t i = 0; i < M0; i++) {
                for (size_t j = 0; j < M1; j++) {
                    for (size_t k = 0; k < M2; k++) {
                        Fp[k * M1 * M0 + i * M1 + j] = Mp[i * M1 * M2 + j * M2 + k];
                    }
                }
            }
        } else if (a1 == 1) {  // (0|12) -> (2|10)
#pragma omp parallel num_threads(nthreads_)
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
void DF_Helper::transpose_disk(std::string name, std::tuple<size_t, size_t, size_t> order) {
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
                error << "DF_Helper:transpose_disk: not enough memory.";
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
    std::vector<double> M;
    std::vector<double> F;
    M.reserve(largest);
    F.reserve(largest);
    double* Mp = M.data();
    double* Fp = F.data();
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
#pragma omp parallel num_threads(nthreads_)
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
#pragma omp parallel num_threads(nthreads_)
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
#pragma omp parallel num_threads(nthreads_)
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
#pragma omp parallel num_threads(nthreads_)
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
#pragma omp parallel num_threads(nthreads_)
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
    file_status_[filename] = file_status_[new_filename];
    stream_check(filename, "rb");
    file_status_.erase(new_filename);

    // keep tsizes_ separate and do not ovwrt sizes_ in case of STORE directive
    files_.erase(new_file);
    tsizes_[filename] = sizes;
}
size_t DF_Helper::get_space_size(std::string name) {
    if (spaces_.find(name) == spaces_.end()) {
        std::stringstream error;
        error << "DF_Helper:get_space_size: " << name << " not found.";
        throw PSIEXCEPTION(error.str().c_str());
    }
    return std::get<1>(spaces_[name]);
}
size_t DF_Helper::get_tensor_size(std::string name) {
    if (transf_.find(name) == transf_.end()) {
        std::stringstream error;
        error << "DF_Helper:get_tensor_size: " << name << " not found.";
        throw PSIEXCEPTION(error.str().c_str());
    }
    std::tuple<size_t, size_t, size_t> s = sizes_[std::get<1>(files_[name])];
    return std::get<0>(s) * std::get<1>(s) * std::get<2>(s);
}
std::tuple<size_t, size_t, size_t> DF_Helper::get_tensor_shape(std::string name) {
    if (transf_.find(name) == transf_.end()) {
        std::stringstream error;
        error << "DF_Helper:get_tensor_size: " << name << " not found.";
        throw PSIEXCEPTION(error.str().c_str());
    }
    return sizes_[std::get<1>(files_[name])];
}
void DF_Helper::build_JK(std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright, std::vector<SharedMatrix> J,
                         std::vector<SharedMatrix> K) {
    timer_on("DFH: build_JK()");
    compute_JK(Cleft, Cright, J, K);
    timer_off("DFH: build_JK()");
}
void DF_Helper::compute_JK(std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright,
                           std::vector<SharedMatrix> J, std::vector<SharedMatrix> K) {
    // outfile->Printf("\n     ==> DF_Helper:--Begin J/K builds <==\n\n");
    // outfile->Printf("\n     ==> Using the %s directive with AO_CORE = %d <==\n\n", method_.c_str(), AO_core_);

    // size checks
    if (Cleft.size() != Cright.size()) {
        std::stringstream error;
        error << "DF_Helper:compute_D - Cleft size (" << Cleft.size() << ") is not equal to Cright size ("
              << Cright.size() << ")";
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (Cleft.size() != J.size()) {
        std::stringstream error;
        error << "DF_Helper:compute_D - Cleft size (" << Cleft.size() << ") is not equal to J size (" << J.size()
              << ")";
        throw PSIEXCEPTION(error.str().c_str());
    }
    if (Cleft.size() != K.size()) {
        std::stringstream error;
        error << "DF_Helper:compute_D - Cleft size (" << Cleft.size() << ") is not equal to K size (" << K.size()
              << ")";
        throw PSIEXCEPTION(error.str().c_str());
    }

    size_t naux = naux_;
    size_t nao = nao_;

    // determine largest buffers needed
    std::vector<std::pair<size_t, size_t>> Qsteps;
    std::tuple<size_t, size_t, size_t, size_t> info = Qshell_blocks_for_JK_build(Qsteps, Cleft, Cright);
    size_t tots = std::get<0>(info);
    size_t totsb = std::get<1>(info);
    size_t wcleft = std::get<2>(info);
    size_t wcright = std::get<3>(info);

    // prep stream, blocking
    if (!direct_ && !AO_core_) stream_check(AO_files_[AO_names_[1]], "rb");

    int rank = 0;
    std::vector<std::vector<double>> C_buffers(nthreads_);
    std::vector<double*> C_bufsp;
    C_bufsp.reserve(nthreads_);

    // prepare eri buffers
    size_t nthread = nthreads_;
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(aux_, zero, primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthreads_);

    // manual cache alignment..(poor std vectors..) (assumes cache size is 64bytes)
    size_t rem = nao * nao;
    size_t align_size = (rem % 8 ? rem + (8 - rem % 8) : rem);

#pragma omp parallel for private(rank) num_threads(nthreads_) schedule(static)
    for (size_t i = 0; i < nthreads_; i++) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        std::vector<double> Cp(align_size);
        C_buffers[rank] = Cp;
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    }

    // declare bufs
    std::vector<double> M;   // AOs
    std::vector<double> T1;  // Tmps
    std::vector<double> T2;

    // cache alignment (is this over-zealous?)
    rem = totsb * wcleft;
    align_size = (rem % 8 ? rem + (8 - rem % 8) : rem);
    T1.reserve(nao * align_size);
    if (JK_hint_) align_size = (nao % 8 ? nao + (8 - nao % 8) : nao);
    T2.reserve(nao * align_size);

    double* T1p = T1.data();
    double* T2p = T2.data();
    double* Mp;

    if (!AO_core_) {
        M.reserve(tots);
        Mp = M.data();
    } else
        Mp = Ppq_.data();

    std::vector<SharedMatrix> D;
    compute_D(D, Cleft, Cright);

    // transform in steps (blocks of Q)
    for (size_t j = 0, bcount = 0; j < Qsteps.size(); j++) {
        // Qshell step info
        size_t start = std::get<0>(Qsteps[j]);
        size_t stop = std::get<1>(Qsteps[j]);
        size_t begin = Qshell_aggs_[start];
        size_t end = Qshell_aggs_[stop + 1] - 1;
        size_t block_size = end - begin + 1;

        // print step info
        // outfile->Printf("      Qshell: (%zu, %zu)", start, stop);
        // outfile->Printf(", PHI: (%zu, %zu), size: %zu\n", begin, end, block_size);

        // get AO chunk according to directive
        timer_on("DFH: Grabbing AOs");
        if (direct_) {
            compute_AO_Q(start, stop, Mp, eri);
        } else if (!AO_core_)
            grab_AO(start, stop, Mp);
        timer_off("DFH: Grabbing AOs");

        timer_on("DFH: compute_J");
        compute_J_symm(D, J, Mp, T1p, T2p, C_buffers, bcount, block_size);
        timer_off("DFH: compute_J");

        timer_on("DFH: compute_K");
        compute_K(Cleft, Cright, K, T1p, T2p, Mp, bcount, block_size, C_buffers, J, D);
        timer_off("DFH: compute_K");

        bcount += block_size;
    }
    // outfile->Printf("\n     ==> DF_Helper:--End J/K Builds (disk)<==\n\n");
}
void DF_Helper::compute_D(std::vector<SharedMatrix>& D, std::vector<SharedMatrix> Cleft,
                          std::vector<SharedMatrix> Cright) {
    for (size_t i = 0; i < Cleft.size(); i++) {
        std::stringstream s;
        s << "D " << i << " (SO)";
        D.push_back(SharedMatrix(new Matrix(s.str(), nao_, nao_)));
    }

    for (size_t i = 0; i < Cleft.size(); i++) {
        size_t cleft = Cleft[i]->colspi()[0];
        size_t cright = Cright[i]->colspi()[0];
        if (cleft != cright) {
            std::stringstream error;
            error << "DF_Helper:compute_D - Cleft[" << i << "] has space size (" << cleft
                  << "), which is not equal to Cright[" << i << "] space size (" << cright << ")";
            throw PSIEXCEPTION(error.str().c_str());
        }
        double* Clp = Cleft[i]->pointer()[0];
        double* Crp = Cright[i]->pointer()[0];
        double* Dp = D[i]->pointer()[0];
        C_DGEMM('N', 'T', nao_, nao_, cleft, 1.0, Clp, cleft, Crp, cleft, 0.0, Dp, nao_);
    }
}
void DF_Helper::compute_J_symm(std::vector<SharedMatrix> D, std::vector<SharedMatrix> J, double* Mp, double* T1p,
                               double* T2p, std::vector<std::vector<double>> D_buffers, size_t bcount,
                               size_t block_size) {
    size_t nao = nao_;
    size_t naux = naux_;
    int rank = 0;

    for (size_t i = 0; i < J.size(); i++) {
        // grab orbital spaces
        double* Dp = D[i]->pointer()[0];
        double* Jp = J[i]->pointer()[0];

        // cache alignment
        size_t align_size = (naux % 8 ? naux + 8 - naux % 8 : naux);

// initialize Tmp (pQ)
#pragma omp parallel for simd num_threads(nthreads_)
        for (size_t k = 0; k < nthreads_ * align_size; k++) T1p[k] = 0.0;

#pragma omp parallel for private(rank) schedule(guided) num_threads(nthreads_)
        for (size_t k = 0; k < nao; k++) {
            size_t si = small_skips_[k];
            size_t mi = symm_sizes_[k];
            size_t skip = symm_skips_[k];
            size_t jump = (AO_core_ ? big_skips_[k] + bcount * si : (big_skips_[k] * block_size) / naux);

#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif

            for (size_t m = k, sp_count = -1; m < nao; m++) {
                if (schwarz_fun_mask_[k * nao + m]) {
                    sp_count++;
                    D_buffers[rank][sp_count] = (m == k ? Dp[nao * k + m] : 2 * Dp[nao * k + m]);
                }
            }

            // (Qm)(m) -> (Q)
            C_DGEMV('N', block_size, mi, 1.0, &Mp[jump + skip], si, &D_buffers[rank][0], 1, 1.0,
                    &T1p[rank * align_size], 1);
        }
        
        // reduce
        for (size_t k = 1; k < nthreads_; k++) {
            for (size_t l = 0; l < naux; l++) T1p[l] += T1p[k * align_size + l];
        }

        // align cache to reduce false sharing
        align_size = (nao % 8 ? nao + (8 - nao % 8) : nao);

// complete pruned J
#pragma omp parallel for schedule(guided) num_threads(nthreads_)
        for (size_t k = 0; k < nao; k++) {
            size_t si = small_skips_[k];
            size_t mi = symm_sizes_[k];
            size_t skip = symm_skips_[k];
            size_t jump = (AO_core_ ? big_skips_[k] + bcount * si : (big_skips_[k] * block_size) / naux);
            C_DGEMV('T', block_size, mi, 1.0, &Mp[jump + skip], si, T1p, 1, 0.0, &T2p[k * align_size], 1);
        }
        
        // unpack from sparse to dense
        for (size_t k = 0; k < nao; k++) {
            for (size_t m = k + 1, count = 0; m < nao; m++) {  // assumes diagonal exists to avoid if  FIXME
                if (schwarz_fun_mask_[k * nao + m]) {
                    count++;
                    Jp[k * nao + m] += T2p[k * align_size + count];
                    Jp[m * nao + k] += T2p[k * align_size + count];
                }
            }
        }
        for (size_t k = 0; k < nao; k++) Jp[k * nao + k] += T2p[k * align_size];
    }
}
void DF_Helper::compute_J(std::vector<SharedMatrix> D, std::vector<SharedMatrix> J, double* Mp, double* T1p,
                          double* T2p, std::vector<std::vector<double>> D_buffers, size_t bcount, size_t block_size) {
    size_t nao = nao_;
    size_t naux = naux_;
    int rank = 0;

    for (size_t i = 0; i < J.size(); i++) {
        // grab orbital spaces
        double* Dp = D[i]->pointer()[0];
        double* Jp = J[i]->pointer()[0];

        // cache alignment
        size_t align_size = (naux % 8 ? naux + 8 - naux % 8 : naux);

        // initialize Tmp (pQ)
#pragma omp parallel for simd num_threads(nthreads_)
        for (size_t k = 0; k < nthreads_ * align_size; k++) T1p[k] = 0.0;

#pragma omp parallel for firstprivate(nao, naux, block_size) private(rank) schedule(guided) num_threads(nthreads_)
        for (size_t k = 0; k < nao; k++) {
            size_t sp_size = small_skips_[k];
            size_t jump = (AO_core_ ? big_skips_[k] + bcount * sp_size : (big_skips_[k] * block_size) / naux);

#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif

            for (size_t m = 0, sp_count = -1; m < nao; m++) {
                if (schwarz_fun_mask_[k * nao + m]) {
                    sp_count++;
                    D_buffers[rank][sp_count] = Dp[nao * k + m];
                }
            }
            // (Qm)(m) -> (Q)
            C_DGEMV('N', block_size, sp_size, 1.0, &Mp[jump], sp_size, &D_buffers[rank][0], 1, 1.0,
                    &T1p[rank * align_size], 1);
        }
        // reduce
        for (size_t k = 1; k < nthreads_; k++) {
            for (size_t l = 0; l < naux; l++) T1p[l] += T1p[k * align_size + l];
        }

        // align cache to reduce false sharing
        align_size = (nao % 8 ? nao + (8 - nao % 8) : nao);

        // complete pruned J
#pragma omp parallel for schedule(guided) num_threads(nthreads_)
        for (size_t k = 0; k < nao; k++) {
            size_t sp_size = small_skips_[k];
            size_t jump = (AO_core_ ? big_skips_[k] + bcount * sp_size : (big_skips_[k] * block_size) / naux);
            C_DGEMV('T', block_size, sp_size, 1.0, &Mp[jump], sp_size, T1p, 1, 0.0, &T2p[k * align_size], 1);
        }
        // unpack from sparse to dense
        for (size_t k = 0; k < nao; k++) {
            for (size_t m = 0, count = -1; m < nao; m++) {
                if (schwarz_fun_mask_[k * nao + m]) {
                    count++;
                    Jp[k * nao + m] += T2p[k * align_size + count];
                }
            }
        }
    }
}
void DF_Helper::compute_K(std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright,
                          std::vector<SharedMatrix> K, double* Tp, double* T2p, double* Mp, size_t bcount,
                          size_t block_size, std::vector<std::vector<double>> C_buffers, std::vector<SharedMatrix> J,
                          std::vector<SharedMatrix> D) {
    size_t nao = nao_;
    size_t naux = naux_;
    int rank = 0;

    for (size_t i = 0; i < K.size(); i++) {
        // grab orbital spaces
        double* Clp = Cleft[i]->pointer()[0];
        double* Crp = Cright[i]->pointer()[0];
        size_t cleft = Cleft[i]->colspi()[0];
        size_t cright = Cright[i]->colspi()[0];
        double* Kp = K[i]->pointer()[0];

        // cache alignment (is this over-zealous?)
        size_t rem = block_size * cleft;
        size_t align_size = (rem % 8 ? rem + (8 - rem % 8) : rem);

// form temp, thread over spM (nao)
#pragma omp parallel for private(rank) schedule(guided) num_threads(nthreads_)
        for (size_t k = 0; k < nao_; k++) {
            size_t sp_size = small_skips_[k];
            size_t jump = (AO_core_ ? big_skips_[k] + bcount * sp_size : (big_skips_[k] * block_size) / naux_);

#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif

            for (size_t m = 0, sp_count = -1; m < nao_; m++) {
                if (schwarz_fun_mask_[k * nao_ + m]) {
                    sp_count++;
                    C_DCOPY(cleft, &Clp[m * cleft], 1, &C_buffers[rank][sp_count * cleft], 1);
                }
            }
            // (Qm)(mb)->(Qb)
            C_DGEMM('N', 'N', block_size, cleft, sp_size, 1.0, &Mp[jump], sp_size, &C_buffers[rank][0], cleft, 0.0,
                    &Tp[k * align_size], cleft);
        }
        if (!JK_hint_ && (Cleft[i] != Cright[i])) {
#pragma omp parallel for private(rank) schedule(guided) num_threads(nthreads_)
            for (size_t k = 0; k < nao_; k++) {
                size_t sp_size = small_skips_[k];
                size_t jump = (AO_core_ ? big_skips_[k] + bcount * sp_size : (big_skips_[k] * block_size) / naux_);

#ifdef _OPENMP
                rank = omp_get_thread_num();
#endif

                for (size_t m = 0, sp_count = -1; m < nao_; m++) {
                    if (schwarz_fun_mask_[k * nao_ + m]) {
                        sp_count++;
                        C_DCOPY(cright, &Crp[m * cright], 1, &C_buffers[rank][sp_count * cright], 1);
                    }
                }
                // (Qm)(mb)->(Qb)
                C_DGEMM('N', 'N', block_size, cright, sp_size, 1.0, &Mp[jump], sp_size, &C_buffers[rank][0], cright,
                        0.0, &T2p[k * align_size], cright);
            }
            C_DGEMM('N', 'T', nao, nao, cleft * block_size, 1.0, Tp, align_size, T2p, align_size, 1.0, Kp, nao);
        } else {
            C_DGEMM('N', 'T', nao, nao, cleft * block_size, 1.0, Tp, align_size, Tp, align_size, 1.0, Kp, nao);
        }
    }
}
}  // End namespaces
