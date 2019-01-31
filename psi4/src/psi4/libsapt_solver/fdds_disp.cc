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

#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libsapt_solver/fdds_disp.h"
#include "psi4/libmints/3coverlap.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/lib3index/dfhelper.h"

#include <iomanip>

// OMP
#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

namespace sapt {

FDDS_Dispersion::FDDS_Dispersion(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                                 std::map<std::string, SharedMatrix> matrix_cache,
                                 std::map<std::string, SharedVector> vector_cache)
    : primary_(primary), auxiliary_(auxiliary), matrix_cache_(matrix_cache), vector_cache_(vector_cache) {
    Options& options = Process::environment.options;

    // ==> Check incoming cache <==
    std::vector<std::string> matrix_cache_check = {"Cocc_A", "Cvir_A", "Cocc_B", "Cvir_B"};
    for (auto key : matrix_cache_check) {
        if (matrix_cache_.find(key) == matrix_cache_.end()) {
            outfile->Printf("FDDS_Dispersion: Missing matrix_cache key %s\n", key.c_str());
            throw PSIEXCEPTION("FDDS_Dispersion: Missing values in the matrix_cache!");
        }
    }

    std::vector<std::string> vector_cache_check = {"eps_occ_A", "eps_vir_A", "eps_occ_B", "eps_vir_B"};
    for (auto key : vector_cache_check) {
        if (vector_cache_.find(key) == vector_cache_.end()) {
            outfile->Printf("FDDS_Dispersion: Missing vector_cache key %s\n", key.c_str());
            throw PSIEXCEPTION("FDDS_Dispersion: Missing values in the vector_cache!");
        }
    }

    // ==> Form Metric <==

    size_t naux = auxiliary_->nbf();
    metric_ = std::make_shared<Matrix>("Inv Coulomb Metric", naux, naux);

    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    size_t nthread = 1;
#ifdef _OPENMP
    nthread = omp_get_max_threads();
#endif

    // ==> (P|Q) Metric <==
    IntegralFactory metric_factory(auxiliary_, zero, auxiliary_, zero);

    std::vector<std::shared_ptr<TwoBodyAOInt>> metric_ints(nthread);
    std::vector<const double*> metric_buff(nthread);
    for (size_t thread = 0; thread < nthread; thread++) {
        metric_ints[thread] = std::shared_ptr<TwoBodyAOInt>(metric_factory.eri());
        metric_buff[thread] = metric_ints[thread]->buffer();
    }

    double** metricp = metric_->pointer();

#pragma omp parallel for schedule(dynamic) num_threads(nthread)
    for (size_t MU = 0; MU < auxiliary_->nshell(); ++MU) {
        size_t nummu = auxiliary_->shell(MU).nfunction();

        size_t thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // Triangular
        for (size_t NU = 0; NU <= MU; ++NU) {
            size_t numnu = auxiliary_->shell(NU).nfunction();

            metric_ints[thread]->compute_shell(MU, 0, NU, 0);

            size_t index = 0;
            // #pragma simd collapse(2)
            for (size_t mu = 0; mu < nummu; ++mu) {
                size_t omu = auxiliary_->shell(MU).function_index() + mu;

                for (size_t nu = 0; nu < numnu; ++nu, ++index) {
                    size_t onu = auxiliary_->shell(NU).function_index() + nu;

                    metricp[omu][onu] = metricp[onu][omu] = metric_buff[thread][index];
                }
            }
        }
    }

    metric_ints.clear();
    metric_buff.clear();

    metric_inv_ = metric_->clone();
    metric_inv_->power(-1.0, 1.e-12);

    // ==> Form Aux overlap <==

    IntegralFactory factory(auxiliary_);
    std::shared_ptr<OneBodyAOInt> overlap(factory.ao_overlap());
    aux_overlap_ = std::make_shared<Matrix>("Auxiliary Overlap", naux, naux);
    overlap->compute(aux_overlap_);

    // ==> Form 3-index object <==

    // Build C Stack
    std::vector<SharedMatrix> Cstack_vec;
    Cstack_vec.push_back(matrix_cache_["Cocc_A"]);
    Cstack_vec.push_back(matrix_cache_["Cvir_A"]);
    Cstack_vec.push_back(matrix_cache_["Cocc_B"]);
    Cstack_vec.push_back(matrix_cache_["Cvir_B"]);

    size_t doubles = Process::environment.get_memory() * 0.8 / sizeof(double);
    size_t max_MO = 0;
    for (auto& mat : Cstack_vec) max_MO = std::max(max_MO, (size_t)mat->ncol());

    // Build DFHelper
    dfh_ = std::make_shared<DFHelper>(primary_, auxiliary_);
    dfh_->set_memory(doubles);
    dfh_->set_method("DIRECT_iaQ");
    dfh_->set_nthreads(nthread);
    dfh_->set_metric_pow(0.0);
    dfh_->initialize();

    // Define spaces
    dfh_->add_space("i", Cstack_vec[0]);
    dfh_->add_space("a", Cstack_vec[1]);
    dfh_->add_space("j", Cstack_vec[2]);
    dfh_->add_space("b", Cstack_vec[3]);

    // add transformations
    dfh_->add_transformation("iaQ", "i", "a", "pqQ");
    dfh_->add_transformation("jbQ", "j", "b", "pqQ");

    // transform
    dfh_->transform();
}

FDDS_Dispersion::~FDDS_Dispersion() {}

std::vector<SharedMatrix> FDDS_Dispersion::project_densities(std::vector<SharedMatrix> densities) {
    // Perform the contraction
    // (PQS) (S|R)^-1 (R|pq) Dpq -> PQ

    // ==> Contract (R|pq) Dpq <== //

    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    size_t nthread = 1;
#ifdef _OPENMP
    nthread = omp_get_max_threads();
#endif

    // Build integral threads
    IntegralFactory df_factory(auxiliary_, zero, primary_, primary_);

    std::vector<std::shared_ptr<TwoBodyAOInt>> df_ints(nthread);
    std::vector<const double*> df_buff(nthread);
    for (size_t thread = 0; thread < nthread; thread++) {
        df_ints[thread] = std::shared_ptr<TwoBodyAOInt>(df_factory.eri());
        df_buff[thread] = df_ints[thread]->buffer();
    }

    size_t naux = auxiliary_->nbf();
    size_t nbf = primary_->nbf();
    size_t nbf2 = nbf * nbf;

    // Pack the DF pairs
    std::vector<std::pair<size_t, size_t>> df_pairs;
    for (size_t M = 0; M < primary_->nshell(); M++) {
        for (size_t N = 0; N <= M; N++) {
            df_pairs.push_back(std::pair<size_t, size_t>(M, N));
        }
    }

    // Check on memory real quick
    size_t doubles = Process::environment.get_memory() * 0.8 / sizeof(double);
    size_t mem_size = nbf2 * auxiliary_->max_nprimitive() * nthread;
    if (mem_size > doubles) {
        std::stringstream message;
        double mem_gb = ((double)(mem_size) / 0.8 * sizeof(double));
        message << "FDDS Dispersion requires at least nbf^2 * max_ang * nthread of memory." << std::endl;
        message << "       After taxes this is " << std::setprecision(2) << mem_gb << " GB of memory.";
        throw PSIEXCEPTION(message.str());
    }

    // Build result and temp vectors
    std::vector<SharedVector> aux_dens;
    for (size_t i = 0; i < densities.size(); i++) {
        aux_dens.push_back(std::make_shared<Vector>(naux));
    }

    std::vector<SharedMatrix> collapse_temp;
    for (size_t i = 0; i < nthread; i++) {
        collapse_temp.push_back(std::make_shared<Matrix>(auxiliary_->max_function_per_shell(), nbf2));
    }

// Do the contraction
#pragma omp parallel for schedule(dynamic) num_threads(nthread)
    for (size_t Rshell = 0; Rshell < auxiliary_->nshell(); Rshell++) {
        size_t thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        collapse_temp[thread]->zero();
        double** tempp = collapse_temp[thread]->pointer();

        size_t num_r = auxiliary_->shell(Rshell).nfunction();
        size_t index_r = auxiliary_->shell(Rshell).function_index();

        // Loop over our PQ shells
        for (auto PQshell : df_pairs) {
            size_t Pshell = PQshell.first;
            size_t Qshell = PQshell.second;

            df_ints[thread]->compute_shell(Rshell, 0, Pshell, Qshell);

            size_t num_p = primary_->shell(Pshell).nfunction();
            size_t index_p = primary_->shell(Pshell).function_index();

            size_t num_q = primary_->shell(Qshell).nfunction();
            size_t index_q = primary_->shell(Qshell).function_index();

            size_t index = 0;
            for (size_t r = 0; r < num_r; r++) {
                for (size_t p = index_p; p < index_p + num_p; p++) {
                    for (size_t q = index_q; q < index_q + num_q; q++) {
                        tempp[r][p * nbf + q] = tempp[r][q * nbf + p] = df_buff[thread][index++];
                    }
                }
            }
        }

        // Stitch it together
        for (size_t i = 0; i < densities.size(); i++) {
            C_DGEMV('N', num_r, nbf2, 1.0, tempp[0], nbf2, densities[i]->pointer()[0], 1, 0.0,
                    (aux_dens[i]->pointer() + index_r), 1);
        }

    }  // End Rshell

    // Clear a few temps
    collapse_temp.clear();
    df_buff.clear();
    df_ints.clear();
    df_pairs.clear();

    // ==> Contract (R|pq) Dpq <== //
    std::vector<SharedVector> aux_dens_inv;
    for (size_t i = 0; i < densities.size(); i++) {
        aux_dens_inv.push_back(std::make_shared<Vector>(naux));
        aux_dens_inv[i]->gemv(false, 1.0, metric_inv_.get(), aux_dens[i].get(), 0.0);
    }

    // ==> Contract (PQS) S -> PQ <== //
    std::vector<SphericalTransform> trans;
    for (size_t i = 0; i <= auxiliary_->max_am(); i++) {
        trans.push_back(SphericalTransform(i));
    }
    std::vector<std::shared_ptr<ThreeCenterOverlapInt>> aux_ints(nthread);
    std::vector<const double*> aux_buff(nthread);

    for (size_t i = 0; i < nthread; i++) {
        aux_ints[i] = std::shared_ptr<ThreeCenterOverlapInt>(
            new ThreeCenterOverlapInt(trans, auxiliary_, auxiliary_, auxiliary_));
        aux_buff[i] = aux_ints[i]->buffer();
    }

    // Pack the Aux pairs
    std::vector<std::pair<size_t, size_t>> aux_pairs;
    for (size_t M = 0; M < auxiliary_->nshell(); M++) {
        for (size_t N = 0; N <= M; N++) {
            aux_pairs.push_back(std::pair<size_t, size_t>(M, N));
        }
    }

    // Build result and temp vectors
    std::vector<SharedMatrix> ret;
    for (size_t i = 0; i < densities.size(); i++) {
        ret.push_back(std::make_shared<Matrix>(naux, naux));
    }

    size_t max_func = auxiliary_->max_function_per_shell();
    for (size_t i = 0; i < nthread; i++) {
        collapse_temp.push_back(std::make_shared<Matrix>(max_func * max_func, naux));
    }

#pragma omp parallel for schedule(dynamic) num_threads(nthread)
    for (size_t PQi = 0; PQi < aux_pairs.size(); PQi++) {
        size_t thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        size_t Pshell = aux_pairs[PQi].first;
        size_t Qshell = aux_pairs[PQi].second;

        size_t num_p = auxiliary_->shell(Pshell).nfunction();
        size_t index_p = auxiliary_->shell(Pshell).function_index();

        size_t num_q = auxiliary_->shell(Qshell).nfunction();
        size_t index_q = auxiliary_->shell(Qshell).function_index();

        double** tempp = collapse_temp[thread]->pointer();

        // Build aux temp
        for (size_t Rshell = 0; Rshell < auxiliary_->nshell(); Rshell++) {
            size_t num_r = auxiliary_->shell(Rshell).nfunction();
            size_t index_r = auxiliary_->shell(Rshell).function_index();

            aux_ints[thread]->compute_shell(Pshell, Qshell, Rshell);

            size_t index = 0;
            for (size_t p = 0; p < num_p; p++) {
                for (size_t q = 0; q < num_q; q++) {
                    for (size_t r = index_r; r < num_r + index_r; r++) {
                        tempp[p * num_q + q][r] = aux_buff[thread][index++];
                    }
                }
            }
        }

        // Contract back to aux
        for (size_t i = 0; i < densities.size(); i++) {
            double** retp = ret[i]->pointer();
            double* dens = aux_dens_inv[i]->pointer();
            for (size_t p = 0; p < num_p; p++) {
                size_t abs_p = index_p + p;
                for (size_t q = 0; q < num_q; q++) {
                    size_t abs_q = index_q + q;

                    retp[abs_p][abs_q] = retp[abs_q][abs_p] = 2.0 * C_DDOT(naux, tempp[p * num_q + q], 1, dens, 1);
                }
            }
        }
    }

    return ret;
}

SharedMatrix FDDS_Dispersion::form_unc_amplitude(std::string monomer, double omega) {
    // ==> Configuration <==
    SharedVector eps_occ, eps_vir;
    std::string ovQ_tensor_name;

    if (monomer == "A") {
        eps_occ = vector_cache_["eps_occ_A"];
        eps_vir = vector_cache_["eps_vir_A"];
        ovQ_tensor_name = "iaQ";
    } else if (monomer == "B") {
        eps_occ = vector_cache_["eps_occ_B"];
        eps_vir = vector_cache_["eps_vir_B"];
        ovQ_tensor_name = "jbQ";

    } else {
        throw PSIEXCEPTION("FDDS_Dispersion::form_unc_amplitude: Monomer must be A or B!");
    }

    // Sizes
    size_t nocc = eps_occ->dim(0);
    size_t nvir = eps_vir->dim(0);
    size_t naux = auxiliary_->nbf();

    // Check on memory real quick
    size_t doubles = Process::environment.get_memory() * 0.8 / sizeof(double);
    size_t mem_size = 2 * naux * nvir + naux * naux + nvir * nocc;
    if (mem_size > doubles) {
        std::stringstream message;
        double mem_gb = ((double)(mem_size) / 0.8 * sizeof(double));
        message << "FDDS Dispersion requires at least naux * nvir + naux * naux of memory." << std::endl;
        message << "       After taxes this is " << std::setprecision(2) << mem_gb << " GB of memory.";
        throw PSIEXCEPTION(message.str());
    }

    // ==> Uncoupled Amplitudes <==
    auto amp = std::make_shared<Matrix>(nocc, nvir);

    double** ampp = amp->pointer();
    double* eoccp = eps_occ->pointer();
    double* evirp = eps_vir->pointer();

#pragma omp parallel for
    for (size_t i = 0; i < nocc; i++) {
        for (size_t a = 0; a < nvir; a++) {
            double val = -1.0 * (eoccp[i] - evirp[a]);
            double tmp = 4.0 * val / (val * val + omega * omega);
            // Lets see how stable this is, should be fine
            if (tmp < 1.e-14) {
                ampp[i][a] = 0.0;
            } else {
                ampp[i][a] = std::pow(tmp, 0.5);
            }
        }
    }

    // amp->print();

    // ==> Contract <==

    size_t dmem = doubles - naux * naux - nvir * nocc;
    size_t bsize = dmem / (naux * nvir);
    if (bsize > nocc) {
        bsize = nocc;
    }
    size_t nblocks = 1 + ((nocc - 1) / bsize);

    // printf("dmem:    %zu\n", dmem);
    // printf("bsize:   %zu\n", bsize);

    auto ret = std::make_shared<Matrix>("UNC Amplitude", naux, naux);
    auto tmp = std::make_shared<Matrix>("iaQ tmp", bsize * nvir, naux);

    double** tmpp = tmp->pointer();
    double** retp = ret->pointer();

    size_t rstat, osize;
    for (size_t block = 0, bcount = 0; block < nblocks; block++) {
        // printf("Block %zu\n", block);
        if (((block + 1) * bsize) > nocc) {
            osize = nocc - block * bsize;
            tmp->zero();
        } else {
            osize = bsize;
        }

        dfh_->fill_tensor(ovQ_tensor_name, tmp, {bcount, bcount + osize});
        size_t shift_i = block * bsize;

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < osize; i++) {
            for (size_t a = 0; a < nvir; a++) {
                double val = ampp[i + shift_i][a];
#pragma omp simd
                for (size_t Q = 0; Q < naux; Q++) {
                    tmpp[i * nvir + a][Q] *= val;
                }
            }
        }

        ret->gemm(true, false, 1.0, tmp, tmp, 1.0);
        bcount += osize;
    }

    return ret;
}
}  // namespace sapt
}  // namespace psi
