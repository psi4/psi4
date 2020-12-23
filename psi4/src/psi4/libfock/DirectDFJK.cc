/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#include "psi4/libfock/jk.h"

#include "psi4/libmints/siminteri.h"
#include "psi4/lib3index/3index.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/matrix.h"
#include "psi4/pragma.h"
#include "psi4/libqt/qt.h"

#include <cstdlib>
#include <memory>
#include <vector>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"
#include <sstream>
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {

DirectDFJK::DirectDFJK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary)
    : JK(primary), auxiliary_(auxiliary) {
    // Not much happens in common init, just that the number of elements in the
    //  basis sets gets saved as private member data in the DirectDFJK object.
    //  However, the JK() constructor gets called as well. Not to fear!
    //  In our case, all that does is set a couple booleans.
    common_init();
}

DirectDFJK::~DirectDFJK() {}

// It's hard to say what does and doesn't belong in this function. I am
// open to suggestions as to what should change. JSOB
void DirectDFJK::common_init() {
    nbf_ = primary_->nbf();
    naux_ = auxiliary_->nbf();

    p_shells_ = primary_->nshell();
    Q_shells_ = auxiliary_->nshell();
}

// Used as testing machinery. It's nice to have a public function
// to be able to look at private member data from pyside.
// I plan on nuking this as soon as I get this project done,
// don't you worry. JSOB

// syntax copied and adapted from DiskDFJK.
// I decided to forgo copying from MemDFJK because
// that code seems more specialized.
void DirectDFJK::print_header() const {
    if (print_) {
        outfile->Printf("==> DirectDFJK: Density-Fitted J/K Matrices <==\n\n");

        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("Don't interpret these results.\nWe don't have w integrals in DirectDFJK\n");
        outfile->Printf("    OpenMP threads:    %11d\n", omp_nthread_);
        outfile->Printf("    Integrals threads: %11d\n", df_ints_num_threads_);
        outfile->Printf("    Memory [MiB]:      %11ld\n", (memory_ * 8L) / (1024L * 1024L));
        outfile->Printf("    Algorithm:         %11s\n", "Direct");
        outfile->Printf("    Schwarz Cutoff:    %11.0E\n", cutoff_);
        outfile->Printf("    Fitting Condition: %11.0E\n\n", condition_);

        outfile->Printf("   => Auxiliary Basis Set <=\n\n");
        auxiliary_->print_by_level("outfile", print_);
    }
}

// In the long term, I'm not sure what belongs here and what
// belongs in common_init() JSOB

// I think this should call the memory estimate and determine the number of blocks.
// Memory will be allocated locally by the compute_JK() method.

void DirectDFJK::preiterations() {
    // cutoff_ = 1e-12;
    // set_df_ints_num_threads(omp_nthread_);
    sparsity_prep_pQq();

    // prepare blocks tells us at which indices we start and stop each block.
    // I'm making it a separate function because I've seen evidence that
    // DFHelper does some resizing, and I want to be able to match that.
    // dirty
    biggest_shell_ = 0;
    for (int i = 0; i < primary_->nshell(); i++) {
        if (primary_->shell(i).nfunction() > biggest_shell_) {
            biggest_shell_ = primary_->shell(i).nfunction();
        }
    }
    // end dirty

    // prepares the coulomb metric. Of course, our calculations depend on the
    //  inverse of the coulomb metric and the (-0.5) power of the metric, so
    //  we will prepare those immediately afterwards
    // prepare_metric();
    get_met();

    prepare_metric_power(-0.5);

    // for what it's worth, we should have the density matrix
    // prepare_D_symm();
}

void DirectDFJK::sparsity_prep_pQq() {
    int procs = 1;
#ifdef _OPENMP
    procs = omp_nthread_;
#endif

    // variables to hold data
    double global_max_int = 0.0;

    std::unique_ptr<double[]> shel_maxes(new double[p_shells_ * p_shells_]);
    std::unique_ptr<double[]> schwarz_test_vals(new double[nbf_ * nbf_]);

    for (size_t shel_sq_iter = 0; shel_sq_iter < p_shells_ * p_shells_; shel_sq_iter++) shel_maxes[shel_sq_iter] = 0.0;

    for (size_t func_sq_iter = 0; func_sq_iter < nbf_ * nbf_; func_sq_iter++) schwarz_test_vals[func_sq_iter] = 0.0;

    // AO machinery
    auto rifactory = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);

    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(procs);
    std::vector<const double*> buffer(procs);

#pragma omp parallel num_threads(procs)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
        buffer[rank] = eri[rank]->buffer();
    }

    // Sparsity member data
    schwarz_func_ints_.resize(nbf_);

    std::vector<size_t> dense_its(nbf_);

    // AO shell metadata
    size_t Shell_M_count;
    size_t Shell_M_start;
    size_t Shell_N_count;
    size_t Shell_N_start;
    double val = 0.0;

//#pragma omp parallel for schedule(guided) num_threads(procs) private(Shell_M_count, Shell_M_start, Shell_N_count, Shell_N_start, val) reduction(max : global_max_int)
    for (size_t Shell_M = 0; Shell_M < p_shells_; Shell_M++) {
        int rank = 0;
//#ifdef _OPENMP
//        rank = omp_get_thread_num();
//#endif
        Shell_M_start = primary_->shell(Shell_M).function_index();
        Shell_M_count = primary_->shell(Shell_M).nfunction();
        for (size_t Shell_N = 0; Shell_N < p_shells_; Shell_N++) {
            Shell_N_start = primary_->shell(Shell_N).function_index();
            Shell_N_count = primary_->shell(Shell_N).nfunction();
            eri[rank]->compute_shell(Shell_M, Shell_N, Shell_M, Shell_N);
            for (size_t func_m = 0; func_m < Shell_M_count; func_m++) {
                for (size_t func_n = 0; func_n < Shell_N_count; func_n++) {
                    schwarz_test_vals[(Shell_M_start + func_m) * nbf_ + Shell_N_start + func_n] =
                        val = fabs(
                            buffer[rank][func_m * Shell_N_count * Shell_M_count * Shell_N_count +
                                         func_n * Shell_M_count * Shell_N_count + func_m * Shell_N_count + func_n]);
                    if (global_max_int < val) global_max_int = val;

                    if (shel_maxes[Shell_M * p_shells_ + Shell_N] < val) {
                        shel_maxes[Shell_M * p_shells_ + Shell_N]
                            = val;
                    }
                }
            }
        }
    }

    tolerance_ = cutoff_ * cutoff_ / global_max_int;
    schwarz_shell_mask_pQq_.resize(p_shells_);
    schwarz_func_starts_pQq_.resize(p_shells_);
    schwarz_dense_funcs_.resize(p_shells_);
    size_t count;

    for (size_t shell_it_out = 0U; shell_it_out < p_shells_; shell_it_out++) {
        count = 0U;
        for (size_t shell_it_in = 0U; shell_it_in < p_shells_; shell_it_in++) {
            if (shel_maxes[shell_it_out * p_shells_ + shell_it_in] > tolerance_) {
                schwarz_shell_mask_pQq_[shell_it_out].push_back(shell_it_in);
                schwarz_func_starts_pQq_[shell_it_out].push_back(count);
                count += primary_->shell(shell_it_in).nfunction();
            }
        }
        schwarz_dense_funcs_[shell_it_out] = count;
    }
}

// calculates a block of AO's then conducts the proper contractions to get them
//   where they need to go. I don't think it makes sense to have separate J and
//   K function calls as J and K are never fully built in a single function.
/* Functions are named for the memory layout of the 3-index eri tensors. 
 * Multiple functions are kept to save future programmers from having to write
 * Pilot implementations
 * In most cases, only the pQp_sparse function will run. This is due to
 * the fact that it is the only one that takes advantage of Schwarz Screeining 
 * and it performs best in testing. The BB variable controls whether or not
 * the Coulomb Metric Contraction comes before the orbital coefficient
 * contraction ("K1" in DiskDJK) */
void DirectDFJK::compute_JK() {
    BB_ = false;

    if (pQq_) { /* No need for separate blocking and direct functions:
                   Both are part of the same general case encoded      */
        if (ao_sparse_) {
            JK_build_pQq_sparse(); 
            } else {
            JK_build_pQq();
        }

    } else {
        prepare_Q_blocks();
        if (num_blocks_ == 1) {
            build_jk_CC_Qpq_direct(); /* No code for recomputing eri
                                      terms that don't fit in memory */
        } else {
            build_jk_CC_Qpq_blocks(); /* Has code for recomputing eri
                                     terms that don't fit in memory */
        }
    }
}

void DirectDFJK::postiterations() {}

// Determines the maximum amount of memory that DirectDFJK will be using
//   at any one time. To be sure, the memory used to construct AO blocks
//   will only be allocated during the compute_JK() method, so it will
//   be available during other procedures.
// Uses blocking and sparsity information to set total_needs which is
//   returned by memory_estimate.
void DirectDFJK::our_needs() {
    free_memory_ = memory_;
    size_t blocks;

    size_t our_needs = 0;
    double charges_f = 0.0;
    size_t charges_z;
    std::shared_ptr<Molecule> mol_ptr = primary_->molecule();

    for (int atom_iter = 0; atom_iter < mol_ptr->natom(); atom_iter++) charges_f += mol_ptr->fZ(atom_iter);

    /* The fact that this code is necessary is a headache.
     * The actual dimension of the SCF orbital coefficient matrix
     * is unknown when this funciton is called, so half the number 
     * of electrons in the molecule is used as a proxy. A separate
     * function is used during the SCF iterations to determine
     * the number of blocks over which AO ERI's need to be calculated */
    charges_f = charges_f / 2.0;

    charges_z = static_cast<size_t>(charges_f);
    ++charges_z;

    // D
    our_needs += nbf_ * nbf_;
    free_memory_ -= nbf_ * nbf_;

    // C_a
    our_needs += nbf_ * charges_z;
    free_memory_ -= nbf_ * charges_z;
    // C_b
    our_needs += nbf_ * charges_z;
    free_memory_ -= nbf_ * charges_z;
    // K
    our_needs += nbf_ * nbf_;
    free_memory_ -= nbf_ * nbf_;
    // J
    our_needs += nbf_ * nbf_;
    free_memory_ -= nbf_ * nbf_;
    //[ J^{-1} ]
    our_needs += naux_ * naux_;
    free_memory_ -= naux_ * naux_;
    //[ J^{-\frac{1}{2} ]
    our_needs += naux_ * naux_;
    free_memory_ -= naux_ * naux_;
    blocks = (2 * (nbf_ * nbf_ * naux_)) / free_memory_ + 1;
    ABX_block_size_ = ((naux_ / blocks) + 1) * nbf_ * nbf_;

    total_needs_ = our_needs;  // + 2*ABX_block_size_;
}
size_t DirectDFJK::memory_estimate() {
    if (Qpq_) {
        prepare_Q_blocks();
    }

    our_needs();
    return total_needs_;
}

void DirectDFJK::prepare_metric_power(double power) {
    SharedMatrix my_met;
    auto coul_met = std::make_shared<FittingMetric>(auxiliary_, true);
    coul_met->form_fitting_metric();
    my_met = coul_met->get_metric();
    if (!(1e-13 > fabs(1.0 - power))) {
        my_met->power(power, condition_);
    }

    met_cols_.push_back(my_met->ncol());
    met_rows_.push_back(my_met->nrow());
    metric_.push_back(my_met);
    met_powers_.push_back(power);

    while (power < 0.0) {
        power = power + 1.0;
    }
}

double* DirectDFJK::get_metric_power(double power) {
    bool pow_on = false;
    double ret_power;
    size_t met_ind;
    for (size_t met_iter = 0; met_iter < metric_.size(); met_iter++) {
        if (std::fabs(power - met_powers_[met_iter]) < condition_) {
            met_ind = met_iter;
            pow_on = true;
        }
    }
    if (!pow_on) {
        met_ind = metric_.size();
        SharedMatrix my_met;
        auto coul_met = std::make_shared<FittingMetric>(auxiliary_, true);
        coul_met->form_fitting_metric();
        my_met = coul_met->get_metric();
        my_met->power(power, condition_);
        metric_.push_back(my_met);
        met_powers_.push_back(power);
    }
    return metric_[met_ind]->pointer()[0];
}

// Computes the LU decomposition of the coulomb metric and stores it in
//   member vector CMPQ_LU_.
//   pivots are storec in vector PERMUTE_
void DirectDFJK::get_met() {
    auto coul_met = std::make_shared<FittingMetric>(auxiliary_, true);
    int info;

    coul_met->form_fitting_metric();
    SharedMatrix my_met;
    my_met = coul_met->get_metric();
    PERMUTE_.resize(naux_);
    CMPQ_LU_.resize(naux_ * naux_);

    CMPQ_inv_.resize(naux_ * naux_);

//    double ZERO[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//    double* zerop = &ZERO[0];
    std::vector<double> work_p(naux_);
    //double* work_p =  &WORK_P[0];///.get();  // (double*)malloc(sizeof(double) * naux_);
    // C_DCOPY( naux_, zerop, 0, work_p, 1);

    memset(&work_p[0], 0, naux_ * 8);

    double* cmpq_lup = &CMPQ_LU_.front();
    double* cmpq_invp = &CMPQ_inv_.front();
    double* metp = my_met->pointer()[0];
    int* pert = &PERMUTE_.front();

    info = C_DGETRF(naux_, naux_, metp, naux_, pert);
    C_DCOPY(naux_ * naux_, metp, 1, cmpq_lup, 1);

    C_DCOPY(naux_ * naux_, cmpq_lup, 1, cmpq_invp, 1);

    C_DGETRI(naux_, cmpq_invp, naux_, pert, &work_p[0], naux_);

//    if (work_p) {
//        free(work_p);
//    }
}

// called after sparsity_prep_pQq() 
// is called in the memory estimate because the maximum permissable
// block size is needed for an effective memory estimate, and determining
// that is equivalent to estimating the total amount of memory
// Determines the number of blocks and
void DirectDFJK::prepare_p_blocks() {
    Block_funcs_.resize(0);
    Shell_starts_.resize(0);
    Shell_stops_.resize(0);
    k_disps_.resize(0);

    size_t current_costs = 0;
    // block max functions
    size_t bmf = 0;

    size_t block_size = 0;
    size_t shell_funcs;
    size_t biggest_shell = 0;

    // one of two things. Either nelectron/2 + 1 or C_left_ao->ncol()
    size_t charges_z;

    /* Used because psi won't tell us the number of occupied orbitals until 
     * the actual scf iterations, so we need to use this as a stand in.
     * I'm not a huge fan either. */
    if (C_left_ao_.size() == 0) {
        // atomic charges float
        double charges_f;

        /*	for (size_t func_it = 0; func_it < schwarz_func_ints_.size(); func_it++) {
                        sparse_fpf_ += schwarz_func_ints_[func_it];
                }	*/

        // We need to figure out the biggest shell in the primary basis set
        // to allocate memory for AO construction
        for (size_t i = 0; i < primary_->nshell(); i++) {
            if (primary_->shell(i).nfunction() > biggest_shell) {
                biggest_shell = primary_->shell(i).nfunction();
            }
        }

        std::shared_ptr<Molecule> mol_ptr = primary_->molecule();

        for (int atom_iter = 0; atom_iter < mol_ptr->natom(); atom_iter++) charges_f += mol_ptr->fZ(atom_iter);

        charges_f = charges_f / 2.0;

        charges_z = static_cast<size_t>(charges_f);
        ++charges_z;
    } else {
        charges_z = C_left_ao_[0]->ncol();
    }
    //
    // D_ao_ and D_
    current_costs += 3 * nbf_ * nbf_;
    // Coulomb Metrics
    current_costs += 2 * naux_ * naux_;
    current_costs += 4 * nbf_ * nbf_;
    // adding in costs for U and A
    current_costs += biggest_shell * naux_ * charges_z;
    current_costs += biggest_shell * naux_ * charges_z;
    current_costs += omp_nthread_ * 5 * nbf_;

    if (memory_ - current_costs > memory_) {
        std::stringstream error_message;
        error_message << "DirectDFJK: not enough memory for overhead! Overhead requires "
                      << current_costs * 8 / (1024 * 1024) << " [MiB].";
        throw PSIEXCEPTION(error_message.str().c_str());
    }

    Shell_starts_.push_back(0);
    for (size_t shell_iter = 0; shell_iter < p_shells_; shell_iter++) {
        shell_funcs = primary_->shell(shell_iter).nfunction();
        // The 5 is where we're getting our fudge factor we're getting charges_z
        //   to stand in for n_occ, and there should only be a factor of 2 included for
        //   the number of tensors with dimension occ. The extra 3 is a massive amount of
        //   memory that should cover the overhead from the reast of the calculation.
        //   charges_z also includes a fudge factor
        if ((block_size + shell_funcs + 2) * naux_ * (2 * charges_z) > memory_ - current_costs) {
            Block_funcs_.push_back(block_size);
            block_size = shell_funcs;
            Shell_stops_.push_back(shell_iter - 1);
            Shell_starts_.push_back(shell_iter);
        } else {
            block_size += shell_funcs;
            if (block_size > bmf) {
                bmf = block_size;
            }
        }
    }
    Shell_stops_.push_back(p_shells_ - 1);

    Block_funcs_.push_back(block_size);

    biggest_block_ = bmf * naux_ * nbf_;
    num_blocks_ = Block_funcs_.size();

    k_disps_.resize(num_blocks_);
    size_t col_disp = 0U;
    size_t row_disp = 0U;

    for (size_t row_it = 0; row_it < num_blocks_; row_it++) {
        k_disps_[row_it].resize(num_blocks_);
        for (size_t col_it = 0; col_it < num_blocks_; col_it++) {
            k_disps_[row_it][col_it] = nbf_ * row_disp + col_disp;
            col_disp += Block_funcs_[col_it];
        }
        row_disp += Block_funcs_[row_it];
        col_disp = 0U;
    }
}

// decides where all the blocks will stop and start
// initializes Shell_starts_ and Shell_stops_
// initializes num_blocks_
void DirectDFJK::prepare_Q_blocks() {
    //	std::vector<size_t> starts;//Shell_starts_.get();
    //	std::vector<size_t> stops;//stops_.get();
    size_t current_costs = 0;
    double charges_f = 0.0;
    size_t charges_z = 0;

    std::shared_ptr<Molecule> mol_ptr = primary_->molecule();

    for (int atom_iter = 0; atom_iter < mol_ptr->natom(); atom_iter++) charges_f += mol_ptr->fZ(atom_iter);

    charges_f = charges_f / 2.0;

    charges_z = static_cast<size_t>(charges_f);
    ++charges_z;
    current_costs += 3 * nbf_ * nbf_;
    current_costs += 2 * naux_ * naux_;
    current_costs += 4 * nbf_ * nbf_;
    current_costs += 4 * nbf_ * nbf_;

    size_t block_max = (free_memory_) / 2;
    size_t block_size = 0;
    size_t shell_funcs;

    Shell_stops_.clear();
    Shell_starts_.clear();
    Block_funcs_.clear();

    Shell_starts_.push_back(0);
    biggest_block_ = 0;
    for (size_t shell_iter = 0; shell_iter < Q_shells_; shell_iter++) {
        shell_funcs = auxiliary_->shell(shell_iter).nfunction();
        if ((shell_funcs + block_size + 3) * nbf_ * (nbf_ + 2 * charges_z) > memory_ - current_costs) {
            Block_funcs_.push_back(block_size);
            block_size = shell_funcs;
            Shell_stops_.push_back(shell_iter - 1);  //[block_iter] = shell_iter - 1;
            Shell_starts_.push_back(shell_iter);     //[block_iter] = shell_iter;
        } else {
            block_size += shell_funcs;
            if (block_size * nbf_ * nbf_ > biggest_block_) {
                biggest_block_ = block_size * nbf_ * nbf_;
            }
        }
    }

    Block_funcs_.push_back(block_size);

    Shell_stops_.push_back(Q_shells_ - 1);

    num_blocks_ = Shell_starts_.size();
}

void DirectDFJK::compute_dense_AO_block_p_pQq(size_t shell, double* ao_block,
                                              std::vector<std::shared_ptr<TwoBodyAOInt>> eri) {
    int procs = 1;
#ifdef _OPENMP
    procs = omp_nthread_;
#endif

    std::vector<const double*> buffer(procs);

#pragma omp parallel num_threads(procs)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        buffer[rank] = eri[rank]->buffer();
    }

    int rank = 0;
    size_t shell_count = primary_->shell(shell).nfunction();
    size_t ShellP_count;
    size_t ShellP_start;
    size_t ShellN_count;
    size_t ShellN_start;
    size_t Shell_ind_0;
    size_t Buff_ind_0;
    // the value put in mat_size is nbf*naux
    size_t mat_size = nbf_ * naux_;

// In Matt's Code and Rob's code (and, indeed, some of mine), we see 3 loops.
//   however, in this code, we only need two because we only have one shell
//   to compute in the slowest running index, so we're just going to write
//   2 loops to accomodate that.
#pragma omp parallel for schedule(guided) num_threads(procs) private(ShellP_count, ShellP_start, ShellN_count, ShellN_start, Shell_ind_0, Buff_ind_0, rank)
    for (size_t ShellP = 0; ShellP < Q_shells_; ShellP++) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        ShellP_count = auxiliary_->shell(ShellP).nfunction();
        ShellP_start = auxiliary_->shell(ShellP).function_index();
        for (size_t ShellN = 0; ShellN < p_shells_; ShellN++) {
            ShellN_count = primary_->shell(ShellN).nfunction();
            ShellN_start = primary_->shell(ShellN).function_index();
            Buff_ind_0 = ShellN_count * shell_count;
            eri[rank]->compute_shell(ShellP, 0, shell, ShellN);
            int i = 0;
            for (size_t func_m = 0; func_m < shell_count; func_m++) {
                for (size_t func_p = 0; func_p < ShellP_count; func_p++) {
                    for (size_t func_n = 0; func_n < ShellN_count; func_n++) {
                        ao_block[func_m * nbf_ * naux_ + (ShellP_start + func_p) * nbf_ + ShellN_start + func_n] =
                            buffer[rank][func_p * Buff_ind_0 + func_m * ShellN_count + func_n];
                    }
                }
            }
        }
    }
}

void DirectDFJK::compute_sparse_AO_block_p_pQq(size_t shell, double* ao_block,
                                               std::vector<std::shared_ptr<TwoBodyAOInt>> eri) {
    int procs = 1;
#ifdef _OPENMP
    procs = omp_nthread_;
#endif

    std::vector<const double*> buffer(procs);
#pragma omp parallel num_threads(procs)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        buffer[rank] = eri[rank]->buffer();
    }

    // ERI indexing Variables
    // Incides are done with romanization of greek letters \Mu -> M \Nu ->N
    size_t ShellP_start;
    size_t ShellP_count;
    size_t Shell_count = primary_->shell(shell).nfunction();
    size_t ShellN_start;
    size_t ShellN_count;
    size_t Shell_ind_0;
    size_t func_starts_sum;  // schwarz_func_starts_pQq_

// If we try to save something at &ao_block + ShellP_start *nbv_^2, we
//   will get a segfault in all but the 0th block. disp short for displacemtnt

// We're building A_\mu P \nu

// Loop over auxes
#pragma omp parallel for schedule(guided) num_threads(procs) private(ShellP_start, ShellP_count, ShellN_start, ShellN_count, Shell_ind_0)
    for (size_t ShellP = 0; ShellP < Q_shells_; ShellP++) {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        ShellP_start = auxiliary_->shell(ShellP).function_index();
        ShellP_count = auxiliary_->shell(ShellP).nfunction();
        // Inner Loop over prims
        for (size_t ShellN = 0; ShellN < schwarz_shell_mask_pQq_[shell].size(); ShellN++) {
            ShellN_start = schwarz_func_starts_pQq_[shell][ShellN];
            // primary_->shell(schwarz_shell_mask_pQq_[shell][ShellN]).function_index();
            ShellN_count = primary_->shell(schwarz_shell_mask_pQq_[shell][ShellN]).nfunction();
            Shell_ind_0 = Shell_count * ShellN_count;
            // Calculate AO's
            eri[rank]->compute_shell(ShellP, 0, shell, schwarz_shell_mask_pQq_[shell][ShellN]);
            // repack AO's
            for (size_t intp = 0; intp < ShellP_count; intp++) {
                for (size_t intm = 0; intm < Shell_count; intm++) {
                    for (size_t intn = 0; intn < ShellN_count; intn++) {
                        ao_block[intm * schwarz_dense_funcs_[shell] * naux_ +
                                 (intp + ShellP_start) * schwarz_dense_funcs_[shell] + ShellN_start + intn] =
                            // ao_block[ schwarz_func_starts_pQq_[intn] + schwarz_func_ints_[intn] * intp + intm]
                            //                        		=
                            buffer[rank][intp * Shell_ind_0 + intm * ShellN_count + intn];
                    }  // inner loop over prim ints
                }      // outer loop over prim ints
            }          // Loop over aux ints
        }              // Inner Loop over prim shells
    }                  // Loop over aux shells
}

// basic block AO function to get myself off the ground.
// this might not be the optimal blocking scheme to calculate AO's
// I'll feel just fine nuking this and testing replacements as the
// project progresses. JSOB

// A note on functionality: it might seem odd that we're passing member data to
//   our own function, but we don't want this function to think that hard.
//   Passing the ao block might also seem odd, but this will be the biggest
//   single piece of memory used in the scf calculation, so we don't want to
//   reallocate it all the time, and we can't keep it around as member data
//   because then it will live forever, and that's too much memory to keep
//   around. The memory in &ao_block will live in a unique pointer in the
//   module that wraps this function, and in this function, we'll just worry
//   about filling it with what belongs inside.
void DirectDFJK::compute_AO_block_Qpq(size_t start_Q, size_t stop_Q, double* ao_block,
                                      std::vector<std::shared_ptr<TwoBodyAOInt>> eri) {
    int procs = 1;
#ifdef _OPENMP
    procs = omp_nthread_;
#endif

    std::vector<const double*> buffer(procs);

#pragma omp parallel num_threads(procs)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        buffer[rank] = eri[rank]->buffer();
    }

    // ERI indexing Variables
    // Incides are done with romanization of greek letters \Mu -> M \Nu ->N
    size_t ShellP_start;
    size_t ShellP_count;
    size_t ShellM_start;
    size_t ShellM_count;
    size_t ShellN_start;
    size_t ShellN_count;
    size_t Shell_ind_0;

    // If we try to save something at &ao_block + ShellP_start *nbv_^2, we
    //   will get a segfault in all but the 0th block. disp short for displacemtnt
    size_t ShellP_disp = auxiliary_->shell(start_Q).function_index();

    size_t nbf_squared = nbf_ * nbf_;

// We're building A^P_{\mu \nu}
// Loop over auxes
// timer_on("DDF AO_CONST");
#pragma omp parallel for schedule(guided) num_threads(procs) private(ShellP_start, ShellP_count, ShellM_start, ShellM_count, ShellN_start, ShellN_count, Shell_ind_0)
    for (size_t ShellP = start_Q; ShellP <= stop_Q; ShellP++) {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        ShellP_start = auxiliary_->shell(ShellP).function_index();
        ShellP_count = auxiliary_->shell(ShellP).nfunction();
        // Outer Loop over prims
        for (size_t ShellM = 0; ShellM < p_shells_; ShellM++) {
            ShellM_start = primary_->shell(ShellM).function_index();
            ShellM_count = primary_->shell(ShellM).nfunction();
            // Inner Loop over prims
            for (size_t ShellN = ShellM; ShellN < p_shells_; ShellN++) {
                ShellN_start = primary_->shell(ShellN).function_index();
                ShellN_count = primary_->shell(ShellN).nfunction();
                Shell_ind_0 = ShellM_count * ShellN_count;
                // Calculate AO's
                eri[rank]->compute_shell(ShellP, 0, ShellM, ShellN);
                // repack AO's
                for (size_t intp = 0; intp < ShellP_count; intp++) {
                    for (size_t intm = 0; intm < ShellM_count; intm++) {
                        for (size_t intn = 0; intn < ShellN_count; intn++) {
                            ao_block[(ShellP_start - ShellP_disp + intp) * nbf_squared + (ShellM_start + intm) * nbf_ +
                                     ShellN_start + intn] =
                                ao_block[(ShellP_start - ShellP_disp + intp) * nbf_squared +
                                         (ShellN_start + intn) * nbf_ + ShellM_start + intm] =

                                    buffer[rank][intp * Shell_ind_0 + intm * ShellN_count + intn];

                        }  // inner loop over prim ints
                    }      // outer loop over prim ints
                }          // Loop over aux ints
            }              // Inner Loop over prim shells
        }                  // Outer Loop over prim shells
    }                      // Loop over aux shells
                           // timer_off("DDF AO_CONST");
}

void DirectDFJK::build_jk_CC_Qpq_direct() {
printf("Direct\n");

    int procs = 1;
#ifdef _OPENMP
    procs = omp_nthread_;
#endif

    std::shared_ptr<BasisSet> zero_ = BasisSet::zero_ao_basis_set();
    auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, zero_, primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(procs);

#pragma omp parallel num_threads(procs)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    }

    x_size_ = C_left_ao_[0]->ncol() * nbf_ * naux_;

    double* c = C_left_ao_[0]->pointer()[0];
    double* d = D_ao_[0]->pointer()[0];
    double* met_m_1_0 = get_metric_power(-1.0);
    double* met_m_0_5 = get_metric_power(-0.5);

    /* allocation for eri terms */
    std::unique_ptr<double[]> AMN(new double[nbf_ * nbf_ * naux_]);
    //std::unique_ptr<double[]> B(new double[nbf_ * nbf_ * naux_]);
    std::unique_ptr<double[]> AMI(new double[x_size_]);
    std::unique_ptr<double[]> QMI(new double[x_size_]);
    std::unique_ptr<double[]> V(new double[naux_]);
    std::unique_ptr<double[]> PHI(new double[naux_]);

    /* AO ERI terms (A|mn) */
    double* amn = AMN.get();
    //double* b = B.get();
    /* Contracted ERI terms (A|mi) = C_{ni}(A|mn) */
    double* ami = AMI.get();
    /* Contracted ERI terms (Q|mi) = [J^{-1/2}]_{mi} */
    double* qmi = QMI.get();
    double* v = V.get();
    double* phi = PHI.get();

    double* j = J_ao_[0]->pointer()[0];
    double* k = K_ao_[0]->pointer()[0];

    memset(j, 0, nbf_ * nbf_ * sizeof(double));
    memset(k, 0, nbf_ * nbf_ * sizeof(double));
    memset(v, 0, naux_ * sizeof(double) );
    memset(qmi, 0, x_size_ * sizeof(double) );

    // x_slice_.push_back(nbf_*C_left_ao_[0]->ncol());

    size_t nbf_squared = nbf_ * nbf_;

//    for (size_t x_iter = 0; x_iter < x_size_; x_iter++) {
//        qmi[x_iter] = 0.0;
//    }

    timer_on("DDF Qpq AO Construction");
    compute_AO_block_Qpq(Shell_starts_[0], Shell_stops_[0], amn, eri);
    timer_off("DDF Qpq AO Construction");

    // timer_on("Tensor Contractions");

 /* // this is the BB formulation of J: J = BBD
  * //  keep it here to save someone writing their own pilate implementations
    timer_on("DDF Qpq Metric Contraction");
            C_DGEMM( 'N', 'N', naux_, nbf_*nbf_, naux_, 1.0, met_m_0_5, naux_, a, nbf_*nbf_, 0.0, b, nbf_*nbf_);
    timer_off("DDF Qpq Metric Contraction");

    timer_on("DDF Qpq J1");
            C_DGEMV( 'N', naux_, nbf_*nbf_, 1.0, b, nbf_*nbf_, d, 1, 0.0, v, 1 );
    timer_off("DDF Qpq J1");

    timer_on("DDF Qpq J2");
            C_DGEMV( 'T', naux_, nbf_*nbf_, 1.0, b, nbf_*nbf_, v, 1, 0.0, j, 1 );
    timer_off("DDF Qpq J2");

    //changed a to b. changed a to b. changed a to b. changed a to b.
    timer_on("DDF Qpq K1");
            C_DGEMM( 'N', 'N', (int) (nbf_*naux_), C_left_ao_[0]->ncol(), (int) nbf_, 1.0, b, (int) nbf_, c,
    C_left_ao_[0]->ncol(), 0.0, x, C_left_ao_[0]->ncol());
    timer_off("DDF Qpq K1");
    //changed a to b. changed a to b. changed a to b. changed a to b.
    */

    /* cmpq_invp
       c oulomb
       m etric
       pq indexing variables
       _
       inv erse
       p ointer */
    double* cmpq_invp = &CMPQ_inv_.front();
    int* pert = &PERMUTE_.front();

    // density contraction V^A = (A|mn)D_{mn}
    timer_on("DDF Qpq J1");
    C_DGEMV('N', (int)naux_, (int)nbf_squared, 1.0, amn, (int)nbf_squared, d, 1, 0.0, v, 1);
    timer_off("DDF Qpq J1");

    //metric inverse for coulomb matrix contraction PHI /gets \cmpq V
    // P^{Q} = (A|Q) V^A
    C_DGEMV('T', (int)naux_, (int)naux_, 1.0, met_m_1_0, (int)naux_, v, 1, 0.0, phi, 1);
    // C_DGETRS( 'N', naux_, 1, cmpq_invp, naux_, pert, v, naux_);

    //            'N', naux_, 1, metp, naux_, pert, coulomb_vector, naux_

    // C_DCOPY(naux_, v, 1, phi, 1);

    timer_on("DDF Qpq J2");
    C_DGEMV('T', (int)naux_, (int)nbf_squared, 1.0, amn, (int)nbf_squared, phi, 1, 0.0, j, 1);
    timer_off("DDF Qpq J2");

    // (A|mi) = (A|mn)C_{ni}
    timer_on("DDF Qpq K1");
    C_DGEMM('N', 'N', (int)(nbf_ * naux_), C_left_ao_[0]->ncol(), (int)nbf_, 1.0, amn, (int)nbf_, c,
            C_left_ao_[0]->ncol(), 0.0, ami, C_left_ao_[0]->ncol());
    timer_off("DDF Qpq K1");

    // (Q|mi) = (Q|A)(A|mi)
    timer_on("DDF Qpq Metric Contraction");
    C_DGEMM('N', 'N', (int)naux_, ((int)nbf_) * C_left_ao_[0]->ncol(), (int)naux_, 1.0, met_m_0_5, (int)naux_, ami,
            ((int)nbf_) * C_left_ao_[0]->ncol(), 0.0, qmi, ((int)nbf_) * C_left_ao_[0]->ncol());
    timer_off("DDF Qpq Metric Contraction");

    //	Accumulate_K_c_is_c(Block_funcs_[0], x_slice_[0], static_cast<size_t>(C_left_ao_[0]->ncol()), k, x);
    /* Keep to avoid rewriting a pilot implementation. */
    // timer_on("DDF Qpq K2");
    //	C_DGEMM( 'N', 'T', (int) nbf_, (int) nbf_, ((int) naux_) * C_left_ao_[0]->ncol(), 1.0, x, ((int) naux_) *
    //C_left_ao_[0]->ncol(), x, ((int) naux_) * C_left_ao_[0]->ncol()  , 0.0, k, (int) nbf_);
    // timer_off("DDF Qpq K2");


    // K_{mn} = (Q|mi)(Q|ni)
    timer_on("DDF Qpq K2");
    for (size_t R_it = 0; R_it < naux_; R_it++) {
        C_DGEMM('N', 'T', nbf_, nbf_, C_left_ao_[0]->ncol(), 1.0, qmi + R_it * nbf_ * C_left_ao_[0]->ncol(),
                C_left_ao_[0]->ncol(), qmi + R_it * nbf_ * C_left_ao_[0]->ncol(), C_left_ao_[0]->ncol(), 1.0, k, nbf_);
    }
    timer_off("DDF Qpq K2");

    //	J_ao_[0]->save("/theoryfs2/ds/obrien/Debug/Psi4/directdfjk_J.txt", false, false, true);
    //	K_ao_[0]->save("/theoryfs2/ds/obrien/Debug/Psi4/directdfjk_K.txt", false, false, true);

    //	C_left_ao_[0]->save("/theoryfs2/ds/obrien/Debug/Psi4/C_directdfjk.txt", false, false, true);
    // timer_off("Tensor Contractions");
    // x_slice_.erase(x_slice_.begin(), x_slice_.end() );
}

/* same as the above function, but it works in the case where more 
 * eri's must be calculated than fit in memory */
void DirectDFJK::build_jk_CC_Qpq_blocks() {
prepare_Q_blocks();
    int procs = 1;
#ifdef _OPENMP
    procs = omp_nthread_;
#endif

    std::shared_ptr<BasisSet> zero_ = BasisSet::zero_ao_basis_set();
    auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, zero_, primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(procs);

#pragma omp parallel num_threads(procs)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    }

    x_size_ = (biggest_block_ / nbf_) * C_left_ao_[0]->ncol();

    double* c = C_left_ao_[0]->pointer()[0];
    double* d = D_ao_[0]->pointer()[0];
    double* met_m_1_0 = get_metric_power(-1.0);
    double* met_m_0_5 = get_metric_power(-0.5);

    /* Storage allocation for eri terms */
    /* AO ERI terms (A|mn) */
    std::unique_ptr<double[]> AMN(new double[biggest_block_]);
    /* Contracted ERI terms (A|mi) = (A|mn)C_{ni} */
    std::unique_ptr<double[]> AMI(new double[x_size_]);
    /* Doubly Contracted ERI terms (Q|mi) = [J^{-1/2}]_{QA}(A|mi) */
    std::unique_ptr<double[]> QMI(new double[x_size_]);
    /* V^A = (A|mn)D_{mn} */
    std::unique_ptr<double[]> V(new double[naux_]);

    /* Pointers to access eri terms*/
    double* amn = AMN.get();
    double* ami = AMI.get();
    double* qmi = QMI.get();
    double* v = V.get();

    double* j = J_ao_[0]->pointer(0)[0];
    double* k = K_ao_[0]->pointer(0)[0];

    memset(j, 0, nbf_ * nbf_ * sizeof(double));
    memset(k, 0, nbf_ * nbf_ * sizeof(double));
    memset(v, 0, naux_ * sizeof(double) );
    memset(qmi, 0, x_size_ * sizeof(double) );

    x_slice_.push_back(nbf_ * C_left_ao_[0]->ncol());

    // Inner Block offset. Used for pointer math.
    size_t in_block_off{0};
    size_t ou_block_off{0};

    size_t nbf_squared = nbf_ * nbf_;

    double* metp = &CMPQ_LU_.front();
    int* pert = &PERMUTE_.front();

    for (size_t block_iter_ou = 0; block_iter_ou < num_blocks_; block_iter_ou++) {
        // there should be a blas call for this
        for (size_t x_iter = 0; x_iter < x_size_; x_iter++) {
            qmi[x_iter] = 0.0;
        }
        in_block_off = 0;
        if (block_iter_ou == 1) {
            // Metric contraction for Coulomb matrix
            C_DGETRS('N', naux_, 1, metp, naux_, pert, v, naux_);
        }
        for (size_t block_iter_in = 0; block_iter_in < num_blocks_; block_iter_in++) {
            timer_on("DDF Qpq AO Construction");
            //construct terms (A|mn)
            compute_AO_block_Qpq(Shell_starts_[block_iter_in], Shell_stops_[block_iter_in], amn, eri);
            timer_off("DDF Qpq AO Construction");
            if (block_iter_ou == 0) {
                timer_on("DDF Qpq J1");
                // V^A = (A|mn) D_{mn}
                C_DGEMV('N', Block_funcs_[block_iter_in], nbf_squared, 1.0, amn, nbf_squared, d, 1, 0.0, v + in_block_off,
                        1);
                timer_off("DDF Qpq J1");
            }
            if (block_iter_ou == 1) {
                timer_on("DDF Qpq J2");
                // J_{mn} = P^Q (Q|mn)
                C_DGEMV('T', Block_funcs_[block_iter_in], nbf_squared, 1.0, amn, nbf_squared, v + in_block_off, 1, 1.0,
                        j, 1);
                timer_off("DDF Qpq J2");
            }
            // timer_on("Tensor Contractions");
            timer_on("DDF Qpq K1");
            // (A|mi) = (A|mn) C_{ni}
            C_DGEMM('N', 'N', Block_funcs_[block_iter_in] * nbf_, C_left_ao_[0]->ncol(), nbf_, 1.0, amn, nbf_, c,
                    C_left_ao_[0]->ncol(), 0.0, ami, C_left_ao_[0]->ncol());
            timer_off("DDF Qpq K1");
            timer_on("DDF Qpq Metric Contraction");
            // (Q|mi) =  [J^{-1/2}]_{QA}(A|mi)
            C_DGEMM('T', 'N', nbf_ * C_left_ao_[0]->ncol(), Block_funcs_[block_iter_ou], Block_funcs_[block_iter_in],
                    1.0, ami, nbf_ * C_left_ao_[0]->ncol(), met_m_0_5 + (in_block_off * naux_) + ou_block_off, naux_, 1.0,
                    qmi, Block_funcs_[block_iter_ou]);
            timer_off("DDF Qpq Metric Contraction");
            in_block_off += Block_funcs_[block_iter_in];
            // timer_off("Tensor Contractions");
        }
        // timer_on("Tensor Contractions");
        timer_on("DDF Qpq K2");
        // K_{mn} = (Q|mi)(Q|ni)
        C_DGEMM('N', 'T', nbf_, nbf_, Block_funcs_[block_iter_ou] * C_left_ao_[0]->ncol(), 1.0, qmi,
                Block_funcs_[block_iter_ou] * C_left_ao_[0]->ncol(), qmi,
                Block_funcs_[block_iter_ou] * C_left_ao_[0]->ncol(), 1.0, k, nbf_);
        timer_off("DDF Qpq K2");
        // timer_off("Tensor Contractions");
        ou_block_off += Block_funcs_[block_iter_ou];
        //		x_block = x_block + Block_funcs[block_iter_ou] * x_slice_[0];
    }
    x_slice_.clear();
}

/* Similar to the function below except that schwarz screening is not included
 * This code is useful for benchmarking as there is a significant amount
 * of overhead involved with spatial sparsity, and it may be useful to measure
 * the cost and loss of accuracy on a given architecture */
void DirectDFJK::JK_build_pQq() {
    // In principle, this function call should be in preiterations or in
    //   memory estimator. However, it depends on knowing the memory_
    //   value from input which precludes its calling in either of those
    //   places.
    prepare_p_blocks();

    int procs = 1;
#ifdef _OPENMP
    procs = omp_nthread_;
#endif

    std::shared_ptr<BasisSet> zero_ = BasisSet::zero_ao_basis_set();
    auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, zero_, primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(procs);

#pragma omp parallel num_threads(procs)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    }

    int nocc = C_left_ao_[0]->ncol();

    double* j = J_ao_[0]->pointer(0)[0];
    double* k = K_ao_[0]->pointer(0)[0];

    /* Allocation of memory for eri terms */
    /* Eri terms (A|mn) */
    std::unique_ptr<double[]> AMN(new double[biggest_shell_ * naux_ * nbf_]);
    /* Eri terms (A|mi) = (A|mn)C_{ni} */
    std::unique_ptr<double[]> AMI(new double[1U]);
    if (BB_) { /* As mentioned earlier, Coulomb Metric ordering. BB_ is 
                * false by default */
        AMI.reset(new double[biggest_shell_ * naux_ * nbf_]);
    } else {
        AMI.reset(new double[biggest_shell_ * naux_ * nocc]);
    }
    /* Look at comments in DirectDFJK::pQp_sparse to see the reason why there
     * are two of these */
    /* Doubly contracted terms (Q|mi) = [J^{-1/2}]_{QA}(A|mi) 
     * The reason for having two blocks of memory is described in 
     * DirectDFJK::pQp_sparse() and in JSOB thesis (see header) */
    std::unique_ptr<double[]> QMI_NEW(new double[biggest_block_ / nbf_ * nocc]);
    std::unique_ptr<double[]> QMI_OLD(new double[biggest_block_ / nbf_ * nocc]);
    std::unique_ptr<double[]> V(new double[naux_]);

    /* Access pointers to eri arrays */
    double* amn = AMN.get();
    double* ami = AMI.get();
    double* qmi_new = QMI_NEW.get();
    double* qmi_old = QMI_OLD.get();
    double* qmi_hold;
    double* v = V.get();

    memset(j, 0, nbf_ * nbf_ * 8);
    memset(k, 0, nbf_ * nbf_ * 8);
    memset(v, 0, naux_ * 8);

    // for (size_t i = 0; i < naux_; i++) { v[i] = 0.0;}
    char first_char = (num_blocks_ == 1 ? 'B' : 'V');

    size_t qmi_old_ind = 0U;
    size_t qmi_new_ind;

    QMI_Block(first_char, true, 0, amn, qmi_old, ami, v, eri);
    timer_on("DDF pQq K2");
    C_DGEMM('N', 'T', Block_funcs_[0], Block_funcs_[0], naux_ * nocc, 1.0, qmi_old, naux_ * nocc, qmi_old, naux_ * nocc, 1.0, k,
            nbf_);
    timer_off("DDF pQq K2");
    /* x stands for (Q|mi) = qmi to shorten the iterating variable names */
    /* explanation of loops in the pQp_sparse codes. */
    for (size_t xo_iter = 0; xo_iter < num_blocks_ - 1; xo_iter++) {
        for (size_t xn_iter = 1; xn_iter < num_blocks_ - xo_iter; xn_iter++) {
            if (xo_iter == 0 && xn_iter != num_blocks_ - 1) {
                QMI_Block('V', true, xn_iter, amn, qmi_new, ami, v, eri);
                qmi_new_ind = xn_iter;
                timer_on("DDF pQq K2");
                C_DGEMM('N', 'T', Block_funcs_[qmi_new_ind], Block_funcs_[qmi_new_ind], naux_ * nocc, 1.0, qmi_new, naux_ * nocc, qmi_new,
                        naux_ * nocc, 1.0, k + k_disps_[qmi_new_ind][qmi_new_ind], nbf_);
                timer_off("DDF pQq K2");
                timer_on("DDF pQq K2");
                C_DGEMM('N', 'T', Block_funcs_[qmi_old_ind], Block_funcs_[qmi_new_ind], naux_ * nocc, 1.0, qmi_old, naux_ * nocc, qmi_new,
                        naux_ * nocc, 1.0, k + k_disps_[qmi_old_ind][qmi_new_ind], nbf_);
                timer_off("DDF pQq K2");
            } else if (xo_iter == 0 && xn_iter == num_blocks_ - 1) {
                QMI_Block('B', true, xn_iter, amn, qmi_new, ami, v, eri);
                qmi_new_ind = xn_iter;
                timer_on("DDF pQq K2");
                C_DGEMM('N', 'T', Block_funcs_[qmi_new_ind], Block_funcs_[qmi_new_ind], naux_ * nocc, 1.0, qmi_new, naux_ * nocc, qmi_new,
                        naux_ * nocc, 1.0, k + k_disps_[qmi_new_ind][qmi_new_ind], nbf_);
                timer_off("DDF pQq K2");
                timer_on("DDF pQq K2");
                C_DGEMM('N', 'T', Block_funcs_[qmi_old_ind], Block_funcs_[qmi_new_ind], naux_ * nocc, 1.0, qmi_old, naux_ * nocc, qmi_new,
                        naux_ * nocc, 1.0, k + k_disps_[qmi_old_ind][qmi_new_ind], nbf_);
                timer_off("DDF pQq K2");
            } else if (xo_iter == 1) {
                QMI_Block('P', true, xn_iter, amn, qmi_new, ami, v, eri);
                qmi_new_ind = xn_iter;
                timer_on("DDF pQq K2");
                C_DGEMM('N', 'T', Block_funcs_[qmi_new_ind], Block_funcs_[qmi_old_ind], naux_ * nocc, 1.0, qmi_new, naux_ * nocc, qmi_old,
                        naux_ * nocc, 1.0, k + k_disps_[qmi_new_ind][qmi_old_ind], nbf_);
                timer_off("DDF pQq K2");
            } else {
                QMI_Block('N', true, xn_iter, amn, qmi_new, ami, nullptr, eri);
                qmi_new_ind = xn_iter;
                timer_on("DDF pQq K2");
                C_DGEMM('N', 'T', Block_funcs_[qmi_new_ind], Block_funcs_[qmi_old_ind], naux_ * nocc, 1.0, qmi_new, naux_ * nocc, qmi_old,
                        naux_ * nocc, 0.0, k + k_disps_[qmi_new_ind][qmi_old_ind], nbf_);
                timer_off("DDF pQq K2");
            }
        }
        qmi_hold = qmi_new;
        qmi_new = qmi_old;  // XO.get();
        qmi_old = qmi_hold;  // XN.get();
        qmi_old_ind = qmi_new_ind;
    }
    for (size_t kf_i = 0; kf_i < nbf_; kf_i++) {
        for (size_t kf_j = 0; kf_j < kf_i; kf_j++) {
            k[kf_i * nbf_ + kf_j] = k[kf_j * nbf_ + kf_i];
        }
    }

    if (first_char != 'B') {
        QMI_Block('P', false, 0, amn, nullptr, ami, v, eri);
    }
}

void DirectDFJK::JK_build_pQq_sparse() {
    // In principle, this function call should be in preiterations or in
    //   memory estimator. However, it depends on knowing the memory_
    //   value from input which precludes its calling in either of those
    //   places.
    prepare_p_blocks();

    int procs = 1;
#ifdef _OPENMP
    procs = omp_nthread_;
#endif

    std::shared_ptr<BasisSet> zero_ = BasisSet::zero_ao_basis_set();
    auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, zero_, primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(procs);

#pragma omp parallel num_threads(procs)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    }

    int nocc = C_left_ao_[0]->ncol();

    double* j = J_ao_[0]->pointer(0)[0];
    double* k = K_ao_[0]->pointer(0)[0];

    /* Choices of variable letters below from 
       doi.org/10.1039/B204199P */
    /* Allocation for holding uncontracted eri terms (A|mn) */
    std::unique_ptr<double[]> AMN(new double[biggest_shell_ * naux_ * nbf_]);
    /* Allocation for holding contracted eri terms (A|mi) = C_{ n i}(A|m n) */
    std::unique_ptr<double[]> AMI(new double[biggest_shell_ * naux_ * nocc]);
    /* Allocation for holding one block of doubly contracted integrals
       (Q|mi) = C_{mi}(A|mn)[J^{-1/2}]_{AQ} */
    std::unique_ptr<double[]> QMI_NEW(new double[biggest_block_ / nbf_ * nocc]);
    /* Allocation for retaining another block of doubly contracted integrals
       (Q|mi) = C_{mi}(A|mn)[J^{-1/2}]_{AQ} */
    std::unique_ptr<double[]> QMI_OLD(new double[biggest_block_ / nbf_ * nocc]);
    /* Allocation for holding terms for constructing the coulomb matrix */
    std::unique_ptr<double[]> V(new double[naux_]);

    /* Pointer to uncontracted eri terms (A|mn) */
    double* amn = AMN.get();
    /* Pointer to contracted eri terms (A|mi) = C_{ n i}(A|m n) */
    double* ami = AMI.get();
    /* Pointer to one block of doubly contracted integrals
       (Q|mi) = C_{mi}(A|mn)[J^{-1/2}]_{AQ} */
    double* qmi_new = QMI_NEW.get();
    /* Pointer to another block of doubly contracted integrals
       (Q|mi) = C_{mi}(A|mn)[J^{-1/2}]_{AQ} */
    double* qmi_old = QMI_OLD.get();
    /* Pointer used for swapping addresses */
    double* qmi_hold;
    double* v = V.get();

    /* containers for temporary storage of data shaped for sparsity */
    std::unique_ptr<double[]> P_C(new double[nbf_ * nocc]);
    std::unique_ptr<double[]> P_D(new double[nbf_]);
    std::unique_ptr<double[]> P_J(new double[nbf_]);

    /* pointers to temperary storage for data shaped for sparsity */
    double* pruned_c = P_C.get();
    double* pruned_d = P_D.get();
    double* pruned_j = P_J.get();

    /* zeroing out terms that will be accumulated into */
    memset(j, 0, nbf_ * nbf_ * 8);
    memset(k, 0, nbf_ * nbf_ * 8);
    memset(v, 0, naux_ * 8);
    char first_char = (num_blocks_ == 1 ? 'B' : 'V');

    /* variables to track which block of eri's each pointer holds*/
    size_t qmi_old_ind = 0U;
    size_t qmi_new_ind;

    outfile->Printf("num_blocks_ is %zu\n", num_blocks_);

    /* See comments below regarding X_Block */
    QMI_Block_sparse(first_char, true, 0UL, pruned_c, pruned_d, amn, qmi_old, ami, v, pruned_j, eri);
    timer_on("DDF pQq big K DGEMM");
    C_DGEMM('N', 'T', Block_funcs_[0], Block_funcs_[0], naux_ * nocc, 1.0, qmi_old, naux_ * nocc, qmi_old, naux_ * nocc, 1.0, k,
            nbf_);
    timer_off("DDF pQq big K DGEMM");
    /* Loop over blocks of exchange matrix construction. The blocking scheme 
       is over blocks (A| [m]_i n) where [m]_i is the i^th block over the 
       primary basis set. The exchange matrix will be constructed as 
       
       old_iter = 0
       x o o o  before loop
       x x o o  new_iter = 1
       x o x o  new_iter = 2
       x o o x  new_iter = 3

       old_iter = 1          old_iter = 2
       x o o o               x o o o
       x x o o               x x o o
       x o x o               x x x o
       x x x x               x x x x
         1 2<----new_iter      1<------new_iter

       final: copy lower to upper
       Legend x: constructed o: empty

       Details can be found by requesting
       Joseph Senan O'Brien thesis from the Georgia Tech Library
       In that document (Q|mi) <--> X^Q_{\mu a}
    */
    for (size_t old_iter = 0; old_iter < num_blocks_ - 1; old_iter++) {
        for (size_t new_iter = 1; new_iter < num_blocks_ - old_iter; new_iter++) {
            if (old_iter == 0 && new_iter != num_blocks_ - 1) {
    /* QMI_Block: function  to: 
       construct integrals (A|mn)
       contract integrals  (A|mi)  = C_{ni}(A|mn)
       contract integrals  V^{A}  += (Q|mn) D_{mn}
       solve               F^P     = (P|Q)^{-1}V^A
       contract integrals  (Q|mi)  = [J^{-1/2}]_{QA} (A|mi) 
       contract integrals  J_{mn} += (A|mn)F^A */
                QMI_Block_sparse('V', true, new_iter, pruned_c, pruned_d, amn, qmi_new, ami, v, pruned_j, eri);
                qmi_new_ind = new_iter;
                timer_on("DDF pQq big K DGEMM");
    /* All D_GEMM's in this function are to calculate
       K_{[m]_j [n]_k} = (Q|[m]_j i)(Q|[n]_k i)  */
                C_DGEMM('N', 'T', Block_funcs_[qmi_new_ind], Block_funcs_[qmi_new_ind], naux_ * nocc, 1.0, qmi_new, naux_ * nocc, qmi_new,
                        naux_ * nocc, 1.0, k + k_disps_[qmi_new_ind][qmi_new_ind], nbf_);
                C_DGEMM('N', 'T', Block_funcs_[qmi_old_ind], Block_funcs_[qmi_new_ind], naux_ * nocc, 1.0, qmi_old, naux_ * nocc, qmi_new,
                        naux_ * nocc, 1.0, k + k_disps_[qmi_old_ind][qmi_new_ind], nbf_);
                timer_off("DDF pQq big K DGEMM");
            } else if (old_iter == 0 && new_iter == num_blocks_ - 1) {
                QMI_Block_sparse('B', true, new_iter, pruned_c, pruned_d, amn, qmi_new, ami, v, pruned_j, eri);
                qmi_new_ind = new_iter;
                timer_on("DDF pQq big K DGEMM");
                C_DGEMM('N', 'T', Block_funcs_[qmi_new_ind], Block_funcs_[qmi_new_ind], naux_ * nocc, 1.0, qmi_new, naux_ * nocc, qmi_new,
                        naux_ * nocc, 1.0, k + k_disps_[qmi_new_ind][qmi_new_ind], nbf_);
                C_DGEMM('N', 'T', Block_funcs_[qmi_old_ind], Block_funcs_[qmi_new_ind], naux_ * nocc, 1.0, qmi_old, naux_ * nocc, qmi_new,
                        naux_ * nocc, 1.0, k + k_disps_[qmi_old_ind][qmi_new_ind], nbf_);
                timer_off("DDF pQq big K DGEMM");
            } else if (old_iter == 1) {
                QMI_Block_sparse('P', true, new_iter, pruned_c, pruned_d, amn, qmi_new, ami, v, pruned_j, eri);
                qmi_new_ind = new_iter;
                timer_on("DDF pQq big K DGEMM");
                C_DGEMM('N', 'T', Block_funcs_[qmi_new_ind], Block_funcs_[qmi_old_ind], naux_ * nocc, 1.0, qmi_new, naux_ * nocc, qmi_old,
                        naux_ * nocc, 1.0, k + k_disps_[qmi_new_ind][qmi_old_ind], nbf_);
                timer_off("DDF pQq big K DGEMM");
            } else {
                QMI_Block_sparse('N', true, new_iter, pruned_c, pruned_d, amn, qmi_new, ami, nullptr, nullptr, eri);
                qmi_new_ind = new_iter;
                timer_on("DDF pQq big K DGEMM");
                C_DGEMM('N', 'T', Block_funcs_[qmi_new_ind], Block_funcs_[qmi_old_ind], naux_ * nocc, 1.0, qmi_new, naux_ * nocc, qmi_old,
                        naux_ * nocc, 0.0, k + k_disps_[qmi_new_ind][qmi_old_ind], nbf_);
                timer_off("DDF pQq big K DGEMM");
            }
        }
        qmi_hold = qmi_new;
        qmi_new = qmi_old;  
        qmi_old = qmi_hold; 
        qmi_old_ind = qmi_new_ind;
    }

    for (size_t kf_i = 0; kf_i < nbf_; kf_i++) {
        for (size_t kf_j = 0; kf_j < kf_i; kf_j++) {
            k[kf_i * nbf_ + kf_j] = k[kf_j * nbf_ + kf_i];
        }
    }

    if (first_char != 'B') {
        QMI_Block_sparse('P', false, 0, pruned_c, pruned_d, amn, nullptr, ami, v, pruned_j, eri);
    }
}

// Function that produces tensors for the final contraction
//    of the exchange matrix build. The issue is that it may or may not
//    have to construct one of various terms for the coulomb matrix.
//    We handle this with a switch.
// coul_work \in { 'V', 'P', 'B', 'N' }
// 'V' means we compute a vector to be contracted against the coulomb Metric.
void DirectDFJK::QMI_Block(char coul_work, bool compute_k, size_t block, double* amn, double* qmi, double* ami,
                         double* coulomb_vector, std::vector<std::shared_ptr<TwoBodyAOInt>> eri) {
    double* cmpq_invp = &CMPQ_inv_.front();
    size_t nocc = C_left_ao_[0]->ncol();
    double* c = C_left_ao_[0]->pointer()[0];
    double* met_m_0_5 = get_metric_power(-0.5);
    double* j = J_ao_[0]->pointer()[0];
    double* d = D_ao_[0]->pointer()[0];
/* Explanation of code blocks in QMI_Block_sparse */
    switch (coul_work) {
        case 'N':
            for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
timer_on("DDF pQq dense AO eri");
                compute_dense_AO_block_p_pQq(shell_iter, amn, eri);
timer_off("DDF pQq dense AO eri");
                if (compute_k) {
timer_on("DDF pQq small K DGEMM");
                    C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction() * naux_, nocc, nbf_, 1.0, amn, nbf_,
                            c, nocc, 0.0, ami, nocc);
timer_off("DDF pQq small K DGEMM");
                    for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
timer_on("DDF pQq Metric Contraction");
                        C_DGEMM('N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, ami + func_it * naux_ * nocc, nocc,
                                0.0, qmi +
                                         (primary_->shell(shell_iter).function_index() -
                                          primary_->shell(Shell_starts_[block]).function_index() + func_it) *
                                             naux_ * nocc,
                                nocc);
timer_off("DDF pQq Metric Contraction");
                    }
                }
            }
            break;
        case 'V':
            for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
                // compute ao blocks
timer_on("DDF pQq dense AO eri");
                compute_dense_AO_block_p_pQq(shell_iter, amn, eri);
timer_off("DDF pQq dense AO eri");
                for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
                    // Form V for Coulomb Matrix construction
timer_on("DDF pQq J orbitals");
                    C_DGEMV('N', (int)naux_, (int)nbf_, 1.0, amn + func_it * naux_ * nbf_, nbf_,
                            d + (primary_->shell(shell_iter).function_index() + func_it) * nbf_, 1, 1.0, coulomb_vector,
                            1);
timer_off("DDF pQq J orbitals");
                }
                if (compute_k) {
timer_on("DDF pQq small K DGEMM");
                    // Form U for Exchange Matrix construction
                    C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction() * naux_, nocc, nbf_, 1.0, amn, nbf_,
                            c, nocc, 0.0, ami, nocc);
timer_off("DDF pQq small K DGEMM");
                    // Contract this u into the corresponding portion of x
                    for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
timer_on("DDF pQq Metric Contraction");
                        C_DGEMM('N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, ami + func_it * naux_ * nocc, nocc,
                                0.0, qmi +
                                         (primary_->shell(shell_iter).function_index() -
                                          primary_->shell(Shell_starts_[block]).function_index() + func_it) *
                                             naux_ * nocc,
                                nocc);
timer_off("DDF pQq Metric Contraction");
                    }
                }
            }
            break;
        case 'P':
            for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
timer_on("DDF pQq dense AO eri");
                compute_dense_AO_block_p_pQq(shell_iter, amn, eri);
timer_off("DDF pQq dense AO eri");
                for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
                    timer_on("DDF pQq J2");
timer_on("DDF pQq J vector");
                    C_DGEMV('T', (int)naux_, (int)nbf_, 1.0, amn + func_it * naux_ * nbf_, nbf_, coulomb_vector, 1,
                            0.0, j + nbf_ * (primary_->shell(shell_iter).function_index() + func_it), 1);
timer_off("DDF pQq J vector");
                    timer_off("DDF pQq J2");
                }
                if (compute_k) {
                    // timer_on("DDF pQq small K DGEMM");
                    timer_on("DDF pQq K1");
timer_on("DDF pQq small K DGEMM");
                    C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction() * naux_, nocc, nbf_, 1.0, amn, nbf_,
                            c, nocc, 0.0, ami, nocc);
timer_off("DDF pQq small K DGEMM");
                    timer_off("DDF pQq K1");
                    for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
timer_on("DDF pQq Metric Contraction");
                        C_DGEMM('N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, qmi + func_it * naux_ * nocc, nocc,
                                0.0, qmi +
                                         (primary_->shell(shell_iter).function_index() -
                                          primary_->shell(Shell_starts_[block]).function_index() + func_it) *
                                             naux_ * nocc,
                                nocc);
timer_off("DDF pQq Metric Contraction");
                    }
                    // timer_off("DDF pQq small K DGEMM");
                }
            }
            break;
        case 'B':
            double* metp = &CMPQ_LU_.front();
            int* pert = &PERMUTE_.front();
            for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
timer_on("DDF pQq dense AO eri");
                compute_dense_AO_block_p_pQq(shell_iter, amn, eri);
timer_off("DDF pQq dense AO eri");
timer_on("DDF pQq small K DGEMM");
                C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction() * naux_, nocc, nbf_, 1.0, amn, nbf_, c,
                        nocc, 0.0, ami, nocc);
timer_off("DDF pQq small K DGEMM");
                for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
timer_on("DDF pQq Metric Contraction");
                    C_DGEMM('N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, ami + func_it * naux_ * nocc, nocc, 0.0,
                            qmi +
                                (primary_->shell(shell_iter).function_index() -
                                 primary_->shell(Shell_starts_[block]).function_index() + func_it) *
                                    naux_ * nocc,
                            nocc);
timer_off("DDF pQq Metric Contraction");
timer_on("DDF pQq J orbitals");
                    C_DGEMV('N', naux_, nbf_, 1.0, amn + func_it * naux_ * nbf_, nbf_,
                            d + (primary_->shell(shell_iter).function_index() + func_it) * nbf_, 1, 1.0, coulomb_vector,
                            1);
timer_on("DDF pQq J orbitals");
                }
                // timer_off("DDF pQq small K DGEMM");
            }

            timer_on("DDF pQq Metric Contraction");
            C_DGETRS('N', naux_, 1, metp, naux_, pert, coulomb_vector, naux_);
            timer_off("DDF pQq Metric Contraction");
            // double* my_vec = (double*) malloc(naux_ * sizeof(double));
            // C_DGEMV('N', naux_, naux_, 1.0, cmpq_invp, naux_, coulomb_vector, naux_, 0.0, my_vec, 1 );
            // C_DCOPY(naux_, my_vec, 1, coulomb_vector, 1);
            // free(my_vec);

            for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
timer_on("DDF pQq dense AO eri");
                compute_dense_AO_block_p_pQq(shell_iter, amn, eri);
timer_on("DDF pQq dense AO eri");
                for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
timer_on("DDF pQq J vector");
                    C_DGEMV('T', naux_, nbf_, 1.0, amn + func_it * naux_ * nbf_, nbf_, coulomb_vector, 1, 0.0,
                            j + (primary_->shell(shell_iter).function_index() + func_it) * nbf_, 1);
timer_off("DDF pQq J vector");
                }
            }
            break;
    }
}

void DirectDFJK::QMI_Block_sparse(char coul_work, bool compute_k, size_t block, double* pruned_c, double* pruned_d,
                                double* ao_block, double* x, double* u, double* coulomb_vector, double* pruned_j,
                                std::vector<std::shared_ptr<TwoBodyAOInt>> eri) {
    size_t nocc = C_left_ao_[0]->ncol();
    double* c = C_left_ao_[0]->pointer()[0];
    double* met_m_0_5 = get_metric_power(-0.5);
    double* j = J_ao_[0]->pointer()[0];
    double* d = D_ao_[0]->pointer()[0];
/*
 We're going to let the timers take the place of comments to match
 The labels in my thesis
DDF pQq sparse AO eri : (A|mn)
DDF pQq small K DGEMM : (A|mi) = (A|mn)C_{ni}
DDF pQq Metric Contraction : (Q|mi) = [J^{-1/2}]_{AQ} (A|mi)
DDF pQq J orbitals : V^A = (A|[m]_i n) D_{[m]_i n}
DDF pQq J vector   : J_{[m]_i n} = 

Not all timers are done in all cases

This function manages construction of doubly contracted eri terms
(Q|mi) = C_{ni} (A|mi) [J^{-1/2}]_{Q A}

as well as construction of the Coulomb matrix

Cases distinguished by what work is done for the coulomb matrix in each of them
N is NO coulomb work
V is accumulation (reduction) into the vector V as V += (A|[m]_i n)D_{[m]_i n}
P is formation of a block of the Coulomb Matrix as 
          J_{[m]_i n } = (A|[m]i n)F^Q
B involves the work in the V and P cases as well as 
    F^P = (P|A)V^A
*/
    switch (coul_work) {
        case 'N':
            for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
timer_on("DDF pQq sparse AO eri");
                compute_sparse_AO_block_p_pQq(shell_iter, ao_block, eri);
timer_off("DDF pQq sparse AO eri");
                if (compute_k) {
                    prune_c(shell_iter, nocc, pruned_c, c);
timer_on("DDF pQq small K DGEMM");
                    C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction() * naux_, nocc,
                            schwarz_dense_funcs_[shell_iter], 1.0, ao_block, schwarz_dense_funcs_[shell_iter], pruned_c,
                            nocc, 0.0, u, nocc);
timer_off("DDF pQq small K DGEMM");
                    for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
timer_on("DDF pQq Metric Contraction");
                        C_DGEMM('N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, u + func_it * naux_ * nocc, nocc,
                                0.0, x +
                                         (primary_->shell(shell_iter).function_index() -
                                          primary_->shell(Shell_starts_[block]).function_index() + func_it) *
                                             naux_ * nocc,
                                nocc);
timer_off("DDF pQq Metric Contraction");
                    }
                }
            }
            break;
        case 'V':
            for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
timer_on("DDF pQq sparse AO eri");
                compute_sparse_AO_block_p_pQq(shell_iter, ao_block, eri);
timer_off("DDF pQq sparse AO eri");
                for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
                    prune_d(shell_iter, pruned_d, d + nbf_ * (func_it + primary_->shell(shell_iter).function_index()));
timer_on("DDF pQq J orbitals");
                    C_DGEMV('N', (int)naux_, (int)schwarz_dense_funcs_[shell_iter], 1.0,
                            ao_block + func_it * naux_ * schwarz_dense_funcs_[shell_iter],
                            schwarz_dense_funcs_[shell_iter], pruned_d, 1, 1.0, coulomb_vector, 1);
timer_off("DDF pQq J orbitals");
                }
                if (compute_k) {
                    prune_c(shell_iter, nocc, pruned_c, c);
timer_on("DDF pQq small K DGEMM");
                    C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction() * naux_, nocc,
                            schwarz_dense_funcs_[shell_iter], 1.0, ao_block, schwarz_dense_funcs_[shell_iter], pruned_c,
                            nocc, 0.0, u, nocc);
timer_off("DDF pQq small K DGEMM");
                    for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
timer_on("DDF pQq Metric Contraction");
                        C_DGEMM('N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, u + func_it * naux_ * nocc, nocc,
                                0.0, x +
                                         (primary_->shell(shell_iter).function_index() -
                                          primary_->shell(Shell_starts_[block]).function_index() + func_it) *
                                             naux_ * nocc,
                                nocc);
timer_off("DDF pQq Metric Contraction");
                    }
                }
            }
            break;
        case 'P':
            for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
timer_on("DDF pQq sparse AO eri");
                compute_sparse_AO_block_p_pQq(shell_iter, ao_block, eri);
timer_off("DDF pQq sparse AO eri");
                for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
timer_on("DDF pQq J vector");
                    C_DGEMV('T', (int)naux_, (int)schwarz_dense_funcs_[shell_iter], 1.0,
                            ao_block + func_it * naux_ * schwarz_dense_funcs_[shell_iter],
                            schwarz_dense_funcs_[shell_iter], coulomb_vector, 1, 0.0, pruned_j, 1);
                    unprune_J(shell_iter, j + nbf_ * (primary_->shell(shell_iter).function_index() + func_it),
                              pruned_j);
timer_off("DDF pQq J vector");
                }
                if (compute_k) {
                    prune_c(shell_iter, nocc, pruned_c, c);
timer_on("DDF pQq small K DGEMM");
                    C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction() * naux_, nocc,
                            schwarz_dense_funcs_[shell_iter], 1.0, ao_block, schwarz_dense_funcs_[shell_iter], pruned_c,
                            nocc, 0.0, u, nocc);
timer_off("DDF pQq small K DGEMM");
                    for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
timer_on("DDF pQq Metric Contraction");
                        C_DGEMM('N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, u + func_it * naux_ * nocc, nocc,
                                0.0, x +
                                         (primary_->shell(shell_iter).function_index() -
                                          primary_->shell(Shell_starts_[block]).function_index() + func_it) *
                                             naux_ * nocc,
                                nocc);
timer_off("DDF pQq Metric Contraction");
                    }
                }
            }
            break;
        case 'B':
            double* metp = &CMPQ_LU_.front();
            int* pert = &PERMUTE_.front();
            for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
timer_on("DDF pQq sparse AO eri");
                compute_sparse_AO_block_p_pQq(shell_iter, ao_block, eri);
timer_off("DDF pQq sparse AO eri");
                if (compute_k) {
                    prune_c(shell_iter, nocc, pruned_c, c);
timer_on("DDF pQq small K DGEMM");
                    C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction() * naux_, nocc,
                            schwarz_dense_funcs_[shell_iter], 1.0, ao_block, schwarz_dense_funcs_[shell_iter], pruned_c,
                            nocc, 0.0, u, nocc);
timer_off("DDF pQq small K DGEMM");
                    for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
timer_on("DDF pQq Metric Contraction");
                        C_DGEMM('N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, u + func_it * naux_ * nocc, nocc,
                                0.0, x +
                                         (primary_->shell(shell_iter).function_index() -
                                          primary_->shell(Shell_starts_[block]).function_index() + func_it) *
                                             naux_ * nocc,
                                nocc);
timer_off("DDF pQq Metric Contraction");
                    }
                }
                for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
                    prune_d(shell_iter, pruned_d, d + nbf_ * (func_it + primary_->shell(shell_iter).function_index()));
timer_on("DDF pQq J orbitals");
                    C_DGEMV('N', naux_, schwarz_dense_funcs_[shell_iter], 1.0,
                            ao_block + func_it * naux_ * schwarz_dense_funcs_[shell_iter],
                            schwarz_dense_funcs_[shell_iter], pruned_d, 1, 1.0, coulomb_vector, 1);
timer_off("DDF pQq J orbitals");
                }
            }

            C_DGETRS('N', naux_, 1, metp, naux_, pert, coulomb_vector, naux_);
            for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
timer_on("DDF pQq sparse AO eri");
                compute_sparse_AO_block_p_pQq(shell_iter, ao_block, eri);
timer_off("DDF pQq sparse AO eri");
                for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
timer_on("DDF pQq J vector");
                    C_DGEMV('T', (int)naux_, (int)schwarz_dense_funcs_[shell_iter], 1.0,
                            ao_block + func_it * naux_ * schwarz_dense_funcs_[shell_iter],
                            schwarz_dense_funcs_[shell_iter], coulomb_vector, 1, 0.0, pruned_j, 1);
                    unprune_J(shell_iter, j + nbf_ * (primary_->shell(shell_iter).function_index() + func_it),
                              pruned_j);
timer_off("DDF pQq J vector");
                }
            }
            break;
    }
}

/* Function to transform the orbital coefficient matrix for sparsity. 
 * Schwarz screening effectively makes a set of primary basis
 * functions that are kept for each member of the primary basis set.
 * Only eri's from these functions are stored, so we have to make a
 * coefficient matrix that only has these rows.*/
void DirectDFJK::prune_c(size_t& mu, size_t nocc, double* pruned_c, double* raw_c) {
    int procs = 1;
#ifdef _OPENMP
    procs = omp_nthread_;
#endif
    size_t pru_add;
    size_t raw_add;

#pragma omp parallel for schedule(dynamic) num_threads(procs)
    for (size_t shell_it = 0U; shell_it < schwarz_shell_mask_pQq_[mu].size(); shell_it++) {
        for (size_t in_shell = 0U; in_shell < primary_->shell(schwarz_shell_mask_pQq_[mu][shell_it]).nfunction();
             in_shell++) {
            for (size_t mol_it = 0U; mol_it < nocc; mol_it++) {
                pruned_c[(schwarz_func_starts_pQq_[mu][shell_it] + in_shell) * nocc + mol_it] =
                    raw_c[(primary_->shell(schwarz_shell_mask_pQq_[mu][shell_it]).function_index() + in_shell) * nocc +
                          mol_it];
            }
        }
    }
}

/* Function to transform the scf density matrix for sparsity. 
 * Schwarz screening effectively makes a set of primary basis
 * functions that are kept for each member of the primary basis set.
 * Only eri's from these functions are stored, so we have to make a
 * density matrix that only has these columns.*/
void DirectDFJK::prune_d(size_t& mu, double* pruned_d, double* raw_d) {
    int procs = 1;
#ifdef _OPENMP
    procs = omp_nthread_;
#endif

    for (size_t shell_it = 0U; shell_it < schwarz_shell_mask_pQq_[mu].size(); shell_it++) {
        for (size_t in_shell = 0U; in_shell < primary_->shell(schwarz_shell_mask_pQq_[mu][shell_it]).nfunction();
             in_shell++) {
            pruned_d[schwarz_func_starts_pQq_[mu][shell_it] + in_shell] =
                raw_d[primary_->shell(schwarz_shell_mask_pQq_[mu][shell_it]).function_index() + in_shell];
        }
    }
}

/* Function to transform the scf density matrix for sparsity. 
 * Schwarz screening effectively makes a set of primary basis
 * functions that are kept for each member of the primary basis set.
 * Only eri's from these functions are stored, so we have to make a
 * density matrix that only has these columns.*/
// You should only pass in one row of J in at a time. We will trust
//    the wrapping function to take care of this.
// mu is a shell index. we need it to get at sparsity information.
void DirectDFJK::unprune_J(size_t& mu, double* raw_j, double* pruned_j) {
    size_t address = 0;
    for (size_t shell_iter = 0; shell_iter < schwarz_shell_mask_pQq_[mu].size(); shell_iter++) {
        for (size_t in_shell = 0; in_shell < primary_->shell(schwarz_shell_mask_pQq_[mu][shell_iter]).nfunction();
             in_shell++) {
            raw_j[primary_->shell(schwarz_shell_mask_pQq_[mu][shell_iter]).function_index() + in_shell] =
                pruned_j[schwarz_func_starts_pQq_[mu][shell_iter] + in_shell];
        }
    }
}

}  // namespace psi
