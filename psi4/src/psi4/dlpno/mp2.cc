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

#include "mp2.h"
#include "sparse.h"

#include "psi4/lib3index/3index.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/local.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/orthog.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libqt/qt.h"

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {
namespace dlpno {

DLPNOMP2::DLPNOMP2(SharedWavefunction ref_wfn, Options& options) : Wavefunction(options) {
    shallow_copy(ref_wfn);
    reference_wavefunction_ = ref_wfn;

    common_init();
}
DLPNOMP2::~DLPNOMP2() {}
void DLPNOMP2::common_init() {
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");

    T_CUT_PNO_ = options_.get_double("T_CUT_PNO");
    T_CUT_DO_ = options_.get_double("T_CUT_DO");

    // did the user manually change expert level options?
    bool T_CUT_PNO_changed = options_["T_CUT_PNO"].has_changed();
    bool T_CUT_DO_changed = options_["T_CUT_DO"].has_changed();

    // if not, values are determined by the user-friendly "PNO_CONVERGENCE"
    if (options_.get_str("PNO_CONVERGENCE") == "LOOSE") {
        if (!T_CUT_PNO_changed) T_CUT_PNO_ = 1e-7;
        if (!T_CUT_DO_changed) T_CUT_DO_ = 2e-2;
    } else if (options_.get_str("PNO_CONVERGENCE") == "NORMAL") {
        if (!T_CUT_PNO_changed) T_CUT_PNO_ = 1e-8;
        if (!T_CUT_DO_changed) T_CUT_DO_ = 1e-2;
    } else if (options_.get_str("PNO_CONVERGENCE") == "TIGHT") {
        if (!T_CUT_PNO_changed) T_CUT_PNO_ = 1e-9;
        if (!T_CUT_DO_changed) T_CUT_DO_ = 5e-3;
    }

    name_ = "DLPNO-MP2";
    module_ = "dlpno";

    variables_["SCF TOTAL ENERGY"] = reference_wavefunction_->energy();

    ribasis_ = get_basisset("DF_BASIS_MP2");
}

/* Utility function for making C_DGESV calls
 *
 * C_DGESV solves AX=B for X, given symmetric NxN matrix A and NxM matrix B
 * B is expected in fortran layout, which complicates the call when (M > 1)
 * The workaround used here is to switch the layout of B before and after the call
 */
void C_DGESV_wrapper(SharedMatrix A, SharedMatrix B) {
    int N = B->rowspi(0);
    int M = B->colspi(0);
    if (N == 0 || M == 0) return;

    // create a copy of B in fortran ordering
    std::vector<double> B_fortran(N * M, 0.0);
    for (int n = 0; n < N; n++) {
        for (int m = 0; m < M; m++) {
            B_fortran[m * N + n] = B->get(n, m);
        }
    }

    // make the C_DGESV call, solving AX=B for X
    std::vector<int> ipiv(N);
    int errcode = C_DGESV(N, M, A->pointer()[0], N, ipiv.data(), B_fortran.data(), N);

    // copy the fortran-ordered X into the original matrix, reverting to C-ordering
    for (int n = 0; n < N; n++) {
        for (int m = 0; m < M; m++) {
            B->set(n, m, B_fortran[m * N + n]);
        }
    }
}

/*
 * In order to use DIIS to accelerate convergence of the MP2 amplitudes,
 * we need to "flatten" the sparse data structure (list of Matrix objects)
 * into a single Matrix
 */
SharedVector flatten_mats(const std::vector<SharedMatrix>& mat_list) {
    size_t total_size = 0;
    for (SharedMatrix mat : mat_list) {
        total_size += mat->size();
    }

    auto flat = std::make_shared<Vector>("flattened matrix list", total_size);
    double* flatp = flat->pointer();

    size_t flat_ind = 0;
    for (SharedMatrix mat : mat_list) {
        if (mat->size() == 0) continue;
        ::memcpy(&flatp[flat_ind], mat->pointer()[0], sizeof(double) * mat->size());
        flat_ind += mat->size();
    }

    return flat;
}

/* This function is a complement to flatten_mats(). A flattened Matrix is
 * copied into a list of Matrix objects.
 */
void copy_flat_mats(SharedVector flat, std::vector<SharedMatrix>& mat_list) {
    double* flatp = flat->pointer();
    size_t flat_ind = 0;
    for (SharedMatrix mat : mat_list) {
        if (mat->size() == 0) continue;
        ::memcpy(mat->pointer()[0], &flatp[flat_ind], sizeof(double) * mat->size());
        flat_ind += mat->size();
    }
}

/* Args: orthonormal orbitals C (ao x mo) and fock matrix F (ao x ao)
 * Return: transformation matrix X (mo x mo) and energy vector e (mo)
 *
 * CX are canonical orbitals (i.e. F(CX) = e(CX))
 */
std::pair<SharedMatrix, SharedVector> canonicalizer(SharedMatrix C, SharedMatrix F) {
    SharedMatrix X = std::make_shared<Matrix>("eigenvectors", C->colspi(0), C->colspi(0));
    SharedVector e = std::make_shared<Vector>("eigenvalues", C->colspi(0));

    auto temp = linalg::triplet(C, F, C, true, false, false);
    temp->diagonalize(X, e, descending);

    return std::make_pair(X, e);
}

/* Args: overlap matrix S (mo x mo), fock matrix F (mo x mo)
 * Return: transformation matrix X (mo x mo_new) and energy vector e (mo_new)
 *
 * CX are orthonormal orbitals (i.e. S_ao (CX) = CX) and also canonical (i.e. F_ao (CX) = e(CX))
 * linear dependencies are removed with keyword S_CUT, so (mo_new <= mo)
 */
std::pair<SharedMatrix, SharedVector> DLPNOMP2::orthocanonicalizer(SharedMatrix S, SharedMatrix F) {
    BasisSetOrthogonalization orthog(BasisSetOrthogonalization::PartialCholesky, S, 0.0, options_.get_double("S_CUT"),
                                     0);
    auto X = orthog.basis_to_orthog_basis();

    int nmo_initial = X->rowspi(0);
    int nmo_final = X->colspi(0);

    auto U = std::make_shared<Matrix>("eigenvectors", nmo_final, nmo_final);
    auto e = std::make_shared<Vector>("eigenvalues", nmo_final);

    auto F_orth = linalg::triplet(X, F, X, true, false, false);
    F_orth->diagonalize(U, e, descending);

    X = linalg::doublet(X, U, false, false);

    return std::make_pair(X, e);
}

void DLPNOMP2::compute_overlap_ints() {
    int nbf = basisset_->nbf();
    int naocc = C_lmo_->colspi(0);
    int npao = C_pao_->colspi(0);  // same as nbf

    timer_on("Construct Grid");
    auto grid = DFTGrid(molecule_, basisset_, options_);
    timer_off("Construct Grid");

    size_t nthread = 1;
#ifdef _OPENMP
    nthread = omp_get_max_threads();
#endif
    std::vector<std::shared_ptr<BasisFunctions>> point_funcs(nthread);
    std::vector<Matrix> DOI_ij_temps(nthread);
    std::vector<Matrix> DOI_iu_temps(nthread);
    for (size_t thread = 0; thread < nthread; thread++) {
        point_funcs[thread] = std::make_shared<BasisFunctions>(basisset_, grid.max_points(), nbf);
        DOI_ij_temps[thread] = Matrix("(i,j) Differential Overlap Integrals", naocc, naocc);
        DOI_iu_temps[thread] = Matrix("(i,u) Differential Overlap Integrals", naocc, nbf);
    }

    timer_on("Integration");
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t Q = 0; Q < grid.blocks().size(); Q++) {
        size_t thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        std::shared_ptr<BlockOPoints> block = grid.blocks()[Q];
        int nbf_block = block->local_nbf();
        int npoints_block = block->npoints();

        // compute values of each basis function at each point in this block
        point_funcs[thread]->compute_functions(block);

        // the values we just computed (max_points x max_functions)
        auto point_values = point_funcs[thread]->basis_values()["PHI"];

        std::vector<int> bf_map = block->functions_local_to_global();

        // resize point_values buffer to size of this block
        auto point_values_trim =
            std::make_shared<Matrix>("DFTGrid PHI Buffer", npoints_block, nbf_block);  // points x bf_block
        for (size_t p = 0; p < npoints_block; p++) {
            for (size_t k = 0; k < nbf_block; k++) {
                point_values_trim->set(p, k, point_values->get(p, k));
            }
        }

        auto C_lmo_slice = submatrix_rows(*C_lmo_, bf_map);  // bf_block x naocc
        auto C_pao_slice = submatrix_rows(*C_pao_, bf_map);  // bf_block x npao

        // value of mo at each point squared
        C_lmo_slice = linalg::doublet(point_values_trim, C_lmo_slice, false, false);  // points x naocc
        C_pao_slice = linalg::doublet(point_values_trim, C_pao_slice, false, false);  // points x npao

        for (size_t p = 0; p < npoints_block; p++) {
            for (size_t i = 0; i < naocc; ++i) {
                C_lmo_slice->set(p, i, pow(C_lmo_slice->get(p, i), 2));
            }
            for (size_t u = 0; u < npao; ++u) {
                C_pao_slice->set(p, u, pow(C_pao_slice->get(p, u), 2));
            }
        }

        auto C_lmo_slice_w = std::make_shared<Matrix>(C_lmo_slice);  // points x naocc
        for (size_t p = 0; p < npoints_block; p++) {
            C_lmo_slice_w->scale_row(0, p, block->w()[p]);
        }

        DOI_ij_temps[thread].add(linalg::doublet(C_lmo_slice_w, C_lmo_slice, true, false));  // naocc x naocc
        DOI_iu_temps[thread].add(linalg::doublet(C_lmo_slice_w, C_pao_slice, true, false));  // naocc x npao
    }
    timer_off("Integration");

    DOI_ij_ = std::make_shared<Matrix>("(i,j) Differential Overlap Integrals", naocc, naocc);
    DOI_iu_ = std::make_shared<Matrix>("(i,u) Differential Overlap Integrals", naocc, nbf);

    for (size_t thread = 0; thread < nthread; thread++) {
        DOI_ij_->add(DOI_ij_temps[thread]);
        DOI_iu_->add(DOI_iu_temps[thread]);
    }

    DOI_ij_->sqrt_this();
    DOI_iu_->sqrt_this();
}

void DLPNOMP2::compute_dipole_ints() {
    int natom = molecule_->natom();
    int naocc = C_lmo_->colspi(0);
    int nbf = C_lmo_->rowspi(0);

    const auto ao_dipole = MintsHelper(basisset_, options_).ao_dipole();

    auto lmo_lmo_dipx = linalg::triplet(C_lmo_, ao_dipole[0], C_lmo_, true, false, false);
    auto lmo_lmo_dipy = linalg::triplet(C_lmo_, ao_dipole[1], C_lmo_, true, false, false);
    auto lmo_lmo_dipz = linalg::triplet(C_lmo_, ao_dipole[2], C_lmo_, true, false, false);

    auto lmo_pao_dipx = linalg::triplet(C_lmo_, ao_dipole[0], C_pao_, true, false, false);
    auto lmo_pao_dipy = linalg::triplet(C_lmo_, ao_dipole[1], C_pao_, true, false, false);
    auto lmo_pao_dipz = linalg::triplet(C_lmo_, ao_dipole[2], C_pao_, true, false, false);

    // < i | dipole | i >
    std::vector<Vector3> R_i;

    // < i | dipole | u >
    std::vector<std::vector<Vector3>> lmo_pao_dr(naocc);

    // e_u
    std::vector<SharedVector> lmo_pao_e(naocc);

    for (size_t i = 0; i < naocc; ++i) {
        R_i.push_back(Vector3(lmo_lmo_dipx->get(i, i), lmo_lmo_dipy->get(i, i), lmo_lmo_dipz->get(i, i)));
    }

    double T_CUT_DO_PRE = options_.get_double("T_CUT_DO_PRE");
    for (size_t i = 0; i < naocc; ++i) {
        std::vector<int> pao_inds;
        for (size_t u = 0; u < nbf; u++) {
            if (fabs(DOI_iu_->get(i, u)) > T_CUT_DO_PRE) {
                pao_inds.push_back(u);
            }
        }
        pao_inds = contract_lists(pao_inds, atom_to_bf_);

        auto C_pao_i = submatrix_cols(*C_pao_, pao_inds);
        auto S_pao_i = submatrix_rows_and_cols(*S_pao_, pao_inds, pao_inds);
        auto F_pao_i = submatrix_rows_and_cols(*F_pao_, pao_inds, pao_inds);

        SharedMatrix X_pao_i;
        SharedVector e_pao_i;
        std::tie(X_pao_i, e_pao_i) = orthocanonicalizer(S_pao_i, F_pao_i);

        auto lmo_pao_dipx_i = submatrix_rows_and_cols(*lmo_pao_dipx, {(int)i}, pao_inds);
        auto lmo_pao_dipy_i = submatrix_rows_and_cols(*lmo_pao_dipy, {(int)i}, pao_inds);
        auto lmo_pao_dipz_i = submatrix_rows_and_cols(*lmo_pao_dipz, {(int)i}, pao_inds);

        lmo_pao_dipx_i = linalg::doublet(lmo_pao_dipx_i, X_pao_i);
        lmo_pao_dipy_i = linalg::doublet(lmo_pao_dipy_i, X_pao_i);
        lmo_pao_dipz_i = linalg::doublet(lmo_pao_dipz_i, X_pao_i);

        int npao_i = X_pao_i->colspi(0);

        for (size_t u = 0; u < npao_i; u++) {
            lmo_pao_dr[i].push_back(
                Vector3(lmo_pao_dipx_i->get(0, u), lmo_pao_dipy_i->get(0, u), lmo_pao_dipz_i->get(0, u)));
        }
        lmo_pao_e[i] = e_pao_i;
    }

    dipole_pair_e_ = std::make_shared<Matrix>("Dipole SC MP2 Energies", naocc, naocc);
    dipole_pair_e_bound_ = std::make_shared<Matrix>("Parallel Dipole SC MP2 Energies", naocc, naocc);

    for (size_t i = 0; i < naocc; ++i) {
        for (size_t j = i + 1; j < naocc; ++j) {
            auto R_ij = R_i[i] - R_i[j];
            auto Rh_ij = R_ij / R_ij.norm();

            double dipole_pair_e_temp = 0.0;
            double dipole_pair_e_bound_temp = 0.0;

            for (int u = 0; u < lmo_pao_dr[i].size(); u++) {
                for (int v = 0; v < lmo_pao_dr[j].size(); v++) {
                    Vector3 iu = lmo_pao_dr[i][u];
                    Vector3 jv = lmo_pao_dr[j][v];

                    double num_actual = iu.dot(jv) - 3 * (iu.dot(Rh_ij) * jv.dot(Rh_ij));
                    num_actual *= num_actual;

                    double num_linear = -2 * iu.dot(jv);
                    num_linear *= num_linear;

                    double denom =
                        (lmo_pao_e[i]->get(u) + lmo_pao_e[j]->get(v)) - (F_lmo_->get(i, i) + F_lmo_->get(j, j));

                    dipole_pair_e_temp += (num_actual / denom);
                    dipole_pair_e_bound_temp += (num_linear / denom);
                }
            }

            dipole_pair_e_temp *= (-4 * pow(R_ij.norm(), -6));
            dipole_pair_e_bound_temp *= (-4 * pow(R_ij.norm(), -6));

            dipole_pair_e_->set(i, j, dipole_pair_e_temp);
            dipole_pair_e_->set(j, i, dipole_pair_e_temp);

            dipole_pair_e_bound_->set(i, j, dipole_pair_e_bound_temp);
            dipole_pair_e_bound_->set(j, i, dipole_pair_e_bound_temp);
        }
    }
}

void DLPNOMP2::prep_sparsity() {
    int natom = molecule_->natom();
    int nbf = basisset_->nbf();
    int nshell = basisset_->nshell();
    int naux = ribasis_->nbf();
    int naocc = nalpha_ - nfrzc();
    int npao = C_pao_->colspi(0);  // same as nbf

    auto bf_to_atom = std::vector<int>(nbf);
    auto ribf_to_atom = std::vector<int>(naux);

    for (size_t i = 0; i < nbf; ++i) {
        bf_to_atom[i] = basisset_->function_to_center(i);
    }

    for (size_t i = 0; i < naux; ++i) {
        ribf_to_atom[i] = ribasis_->function_to_center(i);
    }

    outfile->Printf("  ==> Forming Local MO Domains <==\n");

    // map from LMO to local DF domain (aux basis functions)
    // locality determined via mulliken charges

    lmo_to_ribfs_.resize(naocc);
    lmo_to_riatoms_.resize(naocc);

    for (size_t i = 0; i < naocc; ++i) {
        // atomic mulliken populations for this orbital
        std::vector<double> mkn_pop(natom, 0.0);

        auto P_i = reference_wavefunction_->S()->clone();

        for (size_t u = 0; u < nbf; u++) {
            P_i->scale_row(0, u, C_lmo_->get(u, i));
            P_i->scale_column(0, u, C_lmo_->get(u, i));
        }

        for (size_t u = 0; u < nbf; u++) {
            int centerU = basisset_->function_to_center(u);
            double p_uu = P_i->get(u, u);

            for (size_t v = 0; v < nbf; v++) {
                int centerV = basisset_->function_to_center(v);
                double p_vv = P_i->get(v, v);

                // off-diag pops (p_uv) split between u and v prop to diag pops
                double p_uv = P_i->get(u, v);
                mkn_pop[centerU] += p_uv * ((p_uu) / (p_uu + p_vv));
                mkn_pop[centerV] += p_uv * ((p_vv) / (p_uu + p_vv));
            }
        }

        // if non-zero mulliken pop on atom, include atom in the LMO's fitting domain
        for (size_t a = 0; a < natom; a++) {
            if (fabs(mkn_pop[a]) > options_.get_double("T_CUT_MKN")) {
                lmo_to_riatoms_[i].push_back(a);

                // each atom's aux orbitals are all-or-nothing for each LMO
                for (int u : atom_to_ribf_[a]) {
                    lmo_to_ribfs_[i].push_back(u);
                }
            }
        }
    }

    // map from LMO to local virtual domain (PAOs)
    // locality determined via differential overlap integrals

    lmo_to_paos_.resize(naocc);
    lmo_to_paoatoms_.resize(naocc);
    for (size_t i = 0; i < naocc; ++i) {
        // PAO domains determined by differential overlap integral
        std::vector<int> lmo_to_paos_temp;
        for (size_t u = 0; u < nbf; ++u) {
            if (fabs(DOI_iu_->get(i, u)) > T_CUT_DO_) {
                lmo_to_paos_temp.push_back(u);
            }
        }

        // if any PAO on an atom is in the list, we take all of the PAOs on that atom
        lmo_to_paos_[i] = contract_lists(lmo_to_paos_temp, atom_to_bf_);

        // contains the same information as previous map
        lmo_to_paoatoms_[i] = block_list(lmo_to_paos_[i], bf_to_atom);
    }

    // map from LMO to local occupied domain (other LMOs)
    // locality determined via differential overlap integrals
    //   and also approximated pair energies from dipole integrals

    i_j_to_ij_.resize(naocc);
    de_dipole_ = 0.0;

    for (size_t i = 0, ij = 0; i < naocc; i++) {
        for (size_t j = 0; j < naocc; j++) {
            bool overlap_big = (DOI_ij_->get(i, j) > options_.get_double("T_CUT_DO_ij"));
            bool energy_big = (fabs(dipole_pair_e_bound_->get(i, j)) > options_.get_double("T_CUT_PRE"));

            if (overlap_big || energy_big) {
                i_j_to_ij_[i].push_back(ij);
                ij_to_i_j_.push_back(std::make_pair(i, j));
                ij++;
            } else {
                de_dipole_ += dipole_pair_e_->get(i, j);
                i_j_to_ij_[i].push_back(-1);
            }
        }
    }

    int n_lmo_pairs = ij_to_i_j_.size();

    for (size_t ij = 0; ij < n_lmo_pairs; ++ij) {
        size_t i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        ij_to_ji_.push_back(i_j_to_ij_[j][i]);
    }

    print_aux_domains();
    print_pao_domains();
    print_lmo_domains();

    outfile->Printf("\n  ==> Merging LMO Domains into LMO Pair Domains <==\n");

    // map from (LMO, LMO) pair to local auxiliary and virtual domains
    // LMO pair domains are the union of LMO domains

    lmopair_to_paos_.resize(n_lmo_pairs);
    lmopair_to_paoatoms_.resize(n_lmo_pairs);
    lmopair_to_ribfs_.resize(n_lmo_pairs);
    lmopair_to_riatoms_.resize(n_lmo_pairs);

#pragma omp parallel for
    for (size_t ij = 0; ij < n_lmo_pairs; ++ij) {
        size_t i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        lmopair_to_paos_[ij] = merge_lists(lmo_to_paos_[i], lmo_to_paos_[j]);
        lmopair_to_paoatoms_[ij] = merge_lists(lmo_to_paoatoms_[i], lmo_to_paoatoms_[j]);

        lmopair_to_ribfs_[ij] = merge_lists(lmo_to_ribfs_[i], lmo_to_ribfs_[j]);
        lmopair_to_riatoms_[ij] = merge_lists(lmo_to_riatoms_[i], lmo_to_riatoms_[j]);
    }

    // Create a list of lmos that "interact" with a lmo_pair by differential overlap
    lmopair_to_lmos_.resize(n_lmo_pairs);
    lmopair_to_lmos_dense_.resize(n_lmo_pairs);

#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        lmopair_to_lmos_dense_[ij] = std::vector<int>(naocc, -1);

        int m_ij = 0;
        for (int m = 0; m < naocc; ++m) {
            int im = i_j_to_ij_[i][m];
            int jm = i_j_to_ij_[j][m];

            if (im != -1 && jm != -1) {
                lmopair_to_lmos_[ij].push_back(m);
                lmopair_to_lmos_dense_[ij][m] = m_ij;
                m_ij++;
            }
        }
    }

    print_aux_pair_domains();
    print_pao_pair_domains();

    // => Coefficient Sparsity <= //

    // which basis functions (on which atoms) contribute to each local MO?
    lmo_to_bfs_.resize(naocc);
    lmo_to_atoms_.resize(naocc);

    for (int i = 0; i < naocc; ++i) {
        for (int bf_ind = 0; bf_ind < nbf; ++bf_ind) {
            if (fabs(C_lmo_->get(bf_ind, i)) > options_.get_double("T_CUT_CLMO")) {
                lmo_to_bfs_[i].push_back(bf_ind);
            }
        }
        lmo_to_atoms_[i] = block_list(lmo_to_bfs_[i], bf_to_atom);
    }

    // which basis functions (on which atoms) contribute to each projected AO?
    pao_to_bfs_.resize(nbf);
    pao_to_atoms_.resize(nbf);

    for (int u = 0; u < nbf; ++u) {
        for (int bf_ind = 0; bf_ind < nbf; ++bf_ind) {
            if (fabs(C_pao_->get(bf_ind, u)) > options_.get_double("T_CUT_CPAO")) {
                pao_to_bfs_[u].push_back(bf_ind);
            }
        }
        pao_to_atoms_[u] = block_list(pao_to_bfs_[u], bf_to_atom);
    }

    // determine maps to extended LMO domains, which are the union of an LMO's domain with domains
    //   of all interacting LMOs

    lmo_to_riatoms_ext_ = extend_maps(lmo_to_riatoms_, ij_to_i_j_);
    riatom_to_lmos_ext_ = invert_map(lmo_to_riatoms_ext_, natom);
    riatom_to_paos_ext_ = chain_maps(riatom_to_lmos_ext_, lmo_to_paos_);

    // We'll use these maps to screen the the local MO transform (first index):
    //   (mn|Q) * C_mi -> (in|Q)
    riatom_to_atoms1_ = chain_maps(riatom_to_lmos_ext_, lmo_to_atoms_);
    riatom_to_shells1_ = chain_maps(riatom_to_atoms1_, atom_to_shell_);
    riatom_to_bfs1_ = chain_maps(riatom_to_atoms1_, atom_to_bf_);

    // We'll use these maps to screen the projected AO transform (second index):
    //   (mn|Q) * C_nu -> (mu|Q)
    riatom_to_atoms2_ = chain_maps(riatom_to_lmos_ext_, chain_maps(lmo_to_paos_, pao_to_atoms_));
    riatom_to_shells2_ = chain_maps(riatom_to_atoms2_, atom_to_shell_);
    riatom_to_bfs2_ = chain_maps(riatom_to_atoms2_, atom_to_bf_);

    // Need dense versions of previous maps for quick lookup

    // riatom_to_lmos_ext_dense_[riatom][lmo] is the index of lmo in riatom_to_lmos_ext_[riatom]
    //   (if present), else -1
    riatom_to_lmos_ext_dense_.resize(natom);
    // riatom_to_paos_ext_dense_[riatom][pao] is the index of pao in riatom_to_paos_ext_[riatom]
    //   (if present), else -1
    riatom_to_paos_ext_dense_.resize(natom);

    // riatom_to_atoms1_dense_(1,2)[riatom][a] is true if the orbitals basis functions of atom A
    //   are needed for the (LMO,PAO) transform
    riatom_to_atoms1_dense_.resize(natom);
    riatom_to_atoms2_dense_.resize(natom);

    for (int a_ri = 0; a_ri < natom; a_ri++) {
        riatom_to_lmos_ext_dense_[a_ri] = std::vector<int>(naocc, -1);
        riatom_to_paos_ext_dense_[a_ri] = std::vector<int>(npao, -1);
        riatom_to_atoms1_dense_[a_ri] = std::vector<bool>(natom, false);
        riatom_to_atoms2_dense_[a_ri] = std::vector<bool>(natom, false);

        for (int i_ind = 0; i_ind < riatom_to_lmos_ext_[a_ri].size(); i_ind++) {
            int i = riatom_to_lmos_ext_[a_ri][i_ind];
            riatom_to_lmos_ext_dense_[a_ri][i] = i_ind;
        }
        for (int u_ind = 0; u_ind < riatom_to_paos_ext_[a_ri].size(); u_ind++) {
            int u = riatom_to_paos_ext_[a_ri][u_ind];
            riatom_to_paos_ext_dense_[a_ri][u] = u_ind;
        }
        for (int a_bf : riatom_to_atoms1_[a_ri]) {
            riatom_to_atoms1_dense_[a_ri][a_bf] = true;
        }
        for (int a_bf : riatom_to_atoms2_[a_ri]) {
            riatom_to_atoms2_dense_[a_ri][a_bf] = true;
        }
    }

    lmopair_lmo_to_riatom_lmo_.resize(n_lmo_pairs);
    lmopair_pao_to_riatom_pao_.resize(n_lmo_pairs);
#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int naux_ij = lmopair_to_ribfs_[ij].size();
        lmopair_lmo_to_riatom_lmo_[ij].resize(naux_ij);
        lmopair_pao_to_riatom_pao_[ij].resize(naux_ij);

        for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
            int q = lmopair_to_ribfs_[ij][q_ij];
            int centerq = ribasis_->function_to_center(q);

            int nlmo_ij = lmopair_to_lmos_[ij].size();
            int npao_ij = lmopair_to_paos_[ij].size();
            lmopair_lmo_to_riatom_lmo_[ij][q_ij].resize(nlmo_ij);
            lmopair_pao_to_riatom_pao_[ij][q_ij].resize(npao_ij);

            for (int m_ij = 0; m_ij < nlmo_ij; m_ij++) {
                int m = lmopair_to_lmos_[ij][m_ij];
                int m_sparse = riatom_to_lmos_ext_dense_[centerq][m];
                lmopair_lmo_to_riatom_lmo_[ij][q_ij][m_ij] = m_sparse;
            }

            for (int a_ij = 0; a_ij < npao_ij; a_ij++) {
                int a = lmopair_to_paos_[ij][a_ij];
                int a_sparse = riatom_to_paos_ext_dense_[centerq][a];
                lmopair_pao_to_riatom_pao_[ij][q_ij][a_ij] = a_sparse;
            }
        }
    }
}

void DLPNOMP2::compute_df_ints() {
    timer_on("(mn|K) -> (ia|K)");

    int nbf = basisset_->nbf();
    int naux = ribasis_->nbf();

    size_t nthread = 1;
#ifdef _OPENMP
    nthread = omp_get_max_threads();
#endif

    std::shared_ptr<IntegralFactory> factory =
        std::make_shared<IntegralFactory>(ribasis_, BasisSet::zero_ao_basis_set(), basisset_, basisset_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eris(nthread);

    eris[0] = std::shared_ptr<TwoBodyAOInt>(factory->eri());
    for (size_t thread = 1; thread < nthread; thread++) {
        eris[thread] = std::shared_ptr<TwoBodyAOInt>(eris.front()->clone());
    }

    outfile->Printf("\n  ==> Transforming 3-Index Integrals to LMO/PAO basis <==\n");

    print_integral_sparsity();

    auto SC_lmo =
        linalg::doublet(reference_wavefunction_->S(), C_lmo_, false, false);  // intermediate for coefficient fitting

    qia_.resize(naux);

#pragma omp parallel for schedule(dynamic, 1)
    for (int Q = 0; Q < ribasis_->nshell(); Q++) {
        int nq = ribasis_->shell(Q).nfunction();
        int qstart = ribasis_->shell(Q).function_index();
        int centerQ = ribasis_->shell_to_center(Q);

        size_t thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        // sparse lists of non-screened basis functions
        auto bf_map1 = riatom_to_bfs1_[centerQ];
        auto bf_map2 = riatom_to_bfs2_[centerQ];

        // inverse map, from global (non-screened) bf-index to Q-specific (screened) index
        std::vector<int> bf_map1_inv(nbf, -1);
        std::vector<int> bf_map2_inv(nbf, -1);
        for (int m_ind = 0; m_ind < bf_map1.size(); m_ind++) {
            bf_map1_inv[bf_map1[m_ind]] = m_ind;
        }
        for (int n_ind = 0; n_ind < bf_map2.size(); n_ind++) {
            bf_map2_inv[bf_map2[n_ind]] = n_ind;
        }

        for (size_t q = 0; q < nq; q++) {
            qia_[qstart + q] = std::make_shared<Matrix>("(mn|Q)", bf_map1.size(), bf_map2.size());
        }

        for (int M : riatom_to_shells1_[centerQ]) {
            int nm = basisset_->shell(M).nfunction();
            int mstart = basisset_->shell(M).function_index();
            int centerM = basisset_->shell_to_center(M);

            for (int N : riatom_to_shells2_[centerQ]) {
                int nn = basisset_->shell(N).nfunction();
                int nstart = basisset_->shell(N).function_index();
                int centerN = basisset_->shell_to_center(N);

                // is (N in the list of M's) and (M in the list of N's)?
                bool MN_symmetry =
                    (riatom_to_atoms1_dense_[centerQ][centerN] && riatom_to_atoms2_dense_[centerQ][centerM]);

                // if so, we want to exploit (MN|Q) <-> (NM|Q) symmetry
                if (N < M && MN_symmetry) continue;

                eris[thread]->compute_shell(Q, 0, M, N);
                const double* buffer = eris[thread]->buffer();

                for (int q = 0, index = 0; q < nq; q++) {
                    for (int m = 0; m < nm; m++) {
                        for (int n = 0; n < nn; n++, index++) {
                            qia_[qstart + q]->set(bf_map1_inv[mstart + m], bf_map2_inv[nstart + n], buffer[index]);
                        }
                    }
                }

                // (MN|Q) <-> (NM|Q) symmetry
                if (N > M && MN_symmetry) {
                    for (int q = 0, index = 0; q < nq; q++) {
                        for (int m = 0; m < nm; m++) {
                            for (int n = 0; n < nn; n++, index++) {
                                qia_[qstart + q]->set(bf_map1_inv[nstart + n], bf_map2_inv[mstart + m], buffer[index]);
                            }
                        }
                    }
                }

            }  // N loop
        }      // M loop

        auto C_pao_slice = submatrix_rows_and_cols(*C_pao_, riatom_to_bfs2_[centerQ], riatom_to_paos_ext_[centerQ]);

        //// Here we'll refit the coefficients of C_lmo_slice to minimize residual from unscreened orbitals
        //// This lets us get away with agressive coefficient screening
        //// Boughton and Pulay 1992 JCC, Equation 3

        // Solve for C_lmo_slice such that S[local,local] @ C_lmo_slice ~= S[local,all] @ C_lmo_
        auto C_lmo_slice =
            submatrix_rows_and_cols(*SC_lmo, riatom_to_bfs1_[centerQ], riatom_to_lmos_ext_[centerQ]);
        auto S_aa =
            submatrix_rows_and_cols(*reference_wavefunction_->S(), riatom_to_bfs1_[centerQ], riatom_to_bfs1_[centerQ]);
        C_DGESV_wrapper(S_aa, C_lmo_slice);

        // (mn|Q) C_mi C_nu -> (iu|Q)
        for (size_t q = 0; q < nq; q++) {
            qia_[qstart + q] = linalg::triplet(C_lmo_slice, qia_[qstart + q], C_pao_slice, true, false, false);
        }

    }  // Q loop

    timer_off("(mn|K) -> (ia|K)");

    timer_on("(K|L)");

    // Compute the full metric, don't invert
    auto metric = FittingMetric(ribasis_, true);
    metric.form_fitting_metric();
    full_metric_ = std::make_shared<Matrix>(metric.get_metric());

    timer_off("(K|L)");
}

void DLPNOMP2::pno_transform() {
    int nbf = basisset_->nbf();
    int n_lmo_pairs = ij_to_i_j_.size();

    outfile->Printf("\n  ==> Forming Pair Natural Orbitals <==\n");

    K_iajb_.resize(n_lmo_pairs);   // exchange operators (i.e. (ia|jb) integrals)
    T_iajb_.resize(n_lmo_pairs);   // amplitudes
    Tt_iajb_.resize(n_lmo_pairs);  // antisymmetrized amplitudes
    X_pno_.resize(n_lmo_pairs);    // global PAOs -> canonical PNOs
    e_pno_.resize(n_lmo_pairs);    // PNO orbital energies

    n_pno_.resize(n_lmo_pairs);   // number of pnos
    de_pno_.resize(n_lmo_pairs);  // PNO truncation error
    de_pno_os_.resize(n_lmo_pairs);  // opposite-spin contributions to de_pno_
    de_pno_ss_.resize(n_lmo_pairs);  // same-spin contributions to de_pno_

    D_ij_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij];

        if (i > j) continue;

        //                                                   //
        // ==> Assemble (ia|jb) for pair ij in PAO basis <== //
        //                                                   //

        // number of PAOs in the pair domain (before removing linear dependencies)
        int npao_ij = lmopair_to_paos_[ij].size();  // X_pao_ij->rowspi(0);

        // number of auxiliary basis in the domain
        int naux_ij = lmopair_to_ribfs_[ij].size();

        auto i_qa = std::make_shared<Matrix>("Three-index Integrals", naux_ij, npao_ij);
        auto j_qa = std::make_shared<Matrix>("Three-index Integrals", naux_ij, npao_ij);

        for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
            int q = lmopair_to_ribfs_[ij][q_ij];
            int centerq = ribasis_->function_to_center(q);
            for (int a_ij = 0; a_ij < npao_ij; a_ij++) {
                int a = lmopair_to_paos_[ij][a_ij];
                i_qa->set(q_ij, a_ij, qia_[q]->get(riatom_to_lmos_ext_dense_[centerq][i], riatom_to_paos_ext_dense_[centerq][a]));
                j_qa->set(q_ij, a_ij, qia_[q]->get(riatom_to_lmos_ext_dense_[centerq][j], riatom_to_paos_ext_dense_[centerq][a]));
            }
        }

        auto A_solve = submatrix_rows_and_cols(*full_metric_, lmopair_to_ribfs_[ij], lmopair_to_ribfs_[ij]);
        C_DGESV_wrapper(A_solve, i_qa);

        auto K_pao_ij = linalg::doublet(i_qa, j_qa, true, false);

        //                                      //
        // ==> Canonicalize PAOs of pair ij <== //
        //                                      //

        auto S_pao_ij = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[ij]);
        auto F_pao_ij = submatrix_rows_and_cols(*F_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[ij]);

        SharedMatrix X_pao_ij;  // canonical transformation of this domain's PAOs to
        SharedVector e_pao_ij;  // energies of the canonical PAOs
        std::tie(X_pao_ij, e_pao_ij) = orthocanonicalizer(S_pao_ij, F_pao_ij);

        // S_pao_ij = linalg::triplet(X_pao_ij, S_pao_ij, X_pao_ij, true, false, false);
        F_pao_ij = linalg::triplet(X_pao_ij, F_pao_ij, X_pao_ij, true, false, false);
        K_pao_ij = linalg::triplet(X_pao_ij, K_pao_ij, X_pao_ij, true, false, false);

        // number of PAOs in the domain after removing linear dependencies
        int npao_can_ij = X_pao_ij->colspi(0);
        auto T_pao_ij = K_pao_ij->clone();
        for (int a = 0; a < npao_can_ij; ++a) {
            for (int b = 0; b < npao_can_ij; ++b) {
                T_pao_ij->set(a, b, T_pao_ij->get(a, b) /
                                        (-e_pao_ij->get(b) + -e_pao_ij->get(a) + F_lmo_->get(i, i) + F_lmo_->get(j, j)));
            }
        }

        //                                           //
        // ==> Canonical PAOs  to Canonical PNOs <== //
        //                                           //

        // PNOs defined in (DOI: 10.1063/1.3086717), EQ 17 through EQ 24

        size_t nvir_ij = K_pao_ij->rowspi(0);

        auto Tt_pao_ij = T_pao_ij->clone();
        Tt_pao_ij->scale(2.0);
        Tt_pao_ij->subtract(T_pao_ij->transpose());

        // mp2 energy of this LMO pair before transformation to PNOs
        double e_ij_initial = K_pao_ij->vector_dot(Tt_pao_ij);
        double e_ij_os_initial = K_pao_ij->vector_dot(T_pao_ij);
        double e_ij_ss_initial = e_ij_initial - e_ij_os_initial;

        // Construct pair density from amplitudes
        auto D_ij = linalg::doublet(Tt_pao_ij, T_pao_ij, false, true);
        D_ij->add(linalg::doublet(Tt_pao_ij, T_pao_ij, true, false));
        auto D_ij_copy = D_ij->clone();

        // Diagonalization of pair density gives PNOs (in basis of the LMO's virtual domain) and PNO occ numbers
        auto X_pno_ij = std::make_shared<Matrix>("eigenvectors", nvir_ij, nvir_ij);
        Vector pno_occ("eigenvalues", nvir_ij);
        D_ij->diagonalize(*X_pno_ij, pno_occ, descending);

        int nvir_ij_final = 0;
        for (size_t a = 0; a < nvir_ij; ++a) {
            if (fabs(pno_occ.get(a)) >= T_CUT_PNO_) {
                nvir_ij_final++;
            }
        }

        Dimension zero(1);
        Dimension dim_final(1);
        dim_final.fill(nvir_ij_final);

        // This transformation gives orbitals that are orthonormal but not canonical
        X_pno_ij = X_pno_ij->get_block({zero, X_pno_ij->rowspi()}, {zero, dim_final});
        pno_occ = pno_occ.get_block({zero, dim_final});

        SharedMatrix pno_canon;
        SharedVector e_pno_ij;
        std::tie(pno_canon, e_pno_ij) = canonicalizer(X_pno_ij, F_pao_ij);

        // This transformation gives orbitals that are orthonormal and canonical
        X_pno_ij = linalg::doublet(X_pno_ij, pno_canon, false, false);

        auto K_pno_ij = linalg::triplet(X_pno_ij, K_pao_ij, X_pno_ij, true, false, false);
        auto T_pno_ij = linalg::triplet(X_pno_ij, T_pao_ij, X_pno_ij, true, false, false);
        auto Tt_pno_ij = linalg::triplet(X_pno_ij, Tt_pao_ij, X_pno_ij, true, false, false);

        // mp2 energy of this LMO pair after transformation to PNOs and truncation
        double e_ij_trunc = K_pno_ij->vector_dot(Tt_pno_ij);
        double e_ij_os_trunc = K_pno_ij->vector_dot(T_pno_ij);
        double e_ij_ss_trunc = e_ij_trunc - e_ij_os_trunc;

        // truncation error
        double de_pno_ij = e_ij_initial - e_ij_trunc;
        double de_pno_ij_os = e_ij_os_initial - e_ij_os_trunc;
        double de_pno_ij_ss = e_ij_ss_initial - e_ij_ss_trunc;

        X_pno_ij = linalg::doublet(X_pao_ij, X_pno_ij, false, false);

        K_iajb_[ij] = K_pno_ij;
        T_iajb_[ij] = T_pno_ij;
        X_pno_[ij] = X_pno_ij;
        e_pno_[ij] = e_pno_ij;
        n_pno_[ij] = X_pno_ij->colspi(0);
        de_pno_[ij] = de_pno_ij;
        de_pno_os_[ij] = de_pno_ij_os;
        de_pno_ss_[ij] = de_pno_ij_ss;

        D_ij_[ij] = D_ij_copy;

        // account for symmetry
        if (i < j) {
            K_iajb_[ji] = K_iajb_[ij]->transpose();
            T_iajb_[ji] = T_iajb_[ij]->transpose();
            X_pno_[ji] = X_pno_[ij];
            e_pno_[ji] = e_pno_[ij];
            n_pno_[ji] = n_pno_[ij];
            de_pno_[ji] = de_pno_ij;
            de_pno_os_[ji] = de_pno_ij_os;
            de_pno_ss_[ji] = de_pno_ij_os;

            D_ij_[ji] = D_ij_[ij];
        }
    }

    int pno_count_total = 0, pno_count_min = nbf, pno_count_max = 0;
    de_pno_total_ = 0.0, de_pno_total_os_ = 0.0, de_pno_total_ss_ = 0.0;
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        pno_count_total += n_pno_[ij];
        pno_count_min = std::min(pno_count_min, n_pno_[ij]);
        pno_count_max = std::max(pno_count_max, n_pno_[ij]);
        de_pno_total_ += de_pno_[ij];
        de_pno_total_os_ += de_pno_os_[ij];
        de_pno_total_ss_ += de_pno_ss_[ij];
    }

    outfile->Printf("  \n");
    outfile->Printf("    Natural Orbitals per Local MO pair:\n");
    outfile->Printf("      Avg: %3d NOs \n", pno_count_total / n_lmo_pairs);
    outfile->Printf("      Min: %3d NOs \n", pno_count_min);
    outfile->Printf("      Max: %3d NOs \n", pno_count_max);
    outfile->Printf("  \n");
    outfile->Printf("    PNO truncation energy = %.12f\n", de_pno_total_);

#pragma omp parallel for schedule(static, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        Tt_iajb_[ij] = T_iajb_[ij]->clone();
        Tt_iajb_[ij]->scale(2.0);
        Tt_iajb_[ij]->subtract(T_iajb_[ij]->transpose());
    }
}

void DLPNOMP2::compute_pno_overlaps() {
    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    S_pno_ij_kj_.resize(n_lmo_pairs);
    S_pno_ij_ik_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij];

        if (n_pno_[ij] == 0) continue;
        if (i < j) continue;

        S_pno_ij_kj_[ij] = std::vector<SharedMatrix>(naocc);
        S_pno_ij_ik_[ij] = std::vector<SharedMatrix>(naocc);

        for (int k = 0; k < naocc; ++k) {
            int kj = i_j_to_ij_[k][j];
            if (kj != -1 && n_pno_[kj] > 0) {
                S_pno_ij_kj_[ij][k] = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[kj]);
                S_pno_ij_kj_[ij][k] = linalg::triplet(X_pno_[ij], S_pno_ij_kj_[ij][k], X_pno_[kj], true, false, false);
            }

            int ik = i_j_to_ij_[i][k];
            if (ik != -1 && n_pno_[ik] > 0) {
                S_pno_ij_ik_[ij][k] = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[ik]);
                S_pno_ij_ik_[ij][k] = linalg::triplet(X_pno_[ij], S_pno_ij_ik_[ij][k], X_pno_[ik], true, false, false);
            }
        }

        if (i > j) {
            S_pno_ij_kj_[ji] = S_pno_ij_ik_[ij];
            S_pno_ij_ik_[ji] = S_pno_ij_kj_[ij];
        }
    }
}

void DLPNOMP2::lmp2_iterations() {
    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    outfile->Printf("\n  ==> Local MP2 <==\n\n");
    outfile->Printf("    E_CONVERGENCE = %.2e\n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE = %.2e\n\n", options_.get_double("R_CONVERGENCE"));
    outfile->Printf("                     Corr. Energy    Delta E     Max R\n");

    std::vector<SharedMatrix> R_iajb(n_lmo_pairs);

    int iteration = 0, max_iteration = options_.get_int("DLPNO_MAXITER");
    double e_curr = 0.0, e_prev = 0.0, r_curr = 0.0;
    bool e_converged = false, r_converged = false;
    DIISManager diis(options_.get_int("DIIS_MAX_VECS"), "LMP2 DIIS", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::InCore);

    while (!(e_converged && r_converged)) {
        // RMS of residual per LMO pair, for assessing convergence
        std::vector<double> R_iajb_rms(n_lmo_pairs, 0.0);

        // Calculate residuals from current amplitudes
#pragma omp parallel for schedule(static, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];

            R_iajb[ij] = std::make_shared<Matrix>("Residual", n_pno_[ij], n_pno_[ij]);

            if (n_pno_[ij] == 0) continue;

            for (int a = 0; a < n_pno_[ij]; ++a) {
                for (int b = 0; b < n_pno_[ij]; ++b) {
                    R_iajb[ij]->set(a, b,
                                    K_iajb_[ij]->get(a, b) +
                                        (e_pno_[ij]->get(a) + e_pno_[ij]->get(b) - F_lmo_->get(i, i) - F_lmo_->get(j, j)) *
                                            T_iajb_[ij]->get(a, b));
                }
            }

            for (int k = 0; k < naocc; ++k) {
                int kj = i_j_to_ij_[k][j];
                int ik = i_j_to_ij_[i][k];

                if (kj != -1 && i != k && fabs(F_lmo_->get(i, k)) > options_.get_double("F_CUT") && n_pno_[kj] > 0) {
                    auto temp =
                        linalg::triplet(S_pno_ij_kj_[ij][k], T_iajb_[kj], S_pno_ij_kj_[ij][k], false, false, true);
                    temp->scale(-1.0 * F_lmo_->get(i, k));
                    R_iajb[ij]->add(temp);
                }
                if (ik != -1 && j != k && fabs(F_lmo_->get(k, j)) > options_.get_double("F_CUT") && n_pno_[ik] > 0) {
                    auto temp =
                        linalg::triplet(S_pno_ij_ik_[ij][k], T_iajb_[ik], S_pno_ij_ik_[ij][k], false, false, true);
                    temp->scale(-1.0 * F_lmo_->get(k, j));
                    R_iajb[ij]->add(temp);
                }
            }

            R_iajb_rms[ij] = R_iajb[ij]->rms();
        }

        // evaluate convergence using current amplitudes and residuals
        e_prev = e_curr;
        e_curr = compute_iteration_energy(R_iajb);
        r_curr = *max_element(R_iajb_rms.begin(), R_iajb_rms.end());

        r_converged = (fabs(r_curr) < options_.get_double("R_CONVERGENCE"));
        e_converged = (fabs(e_curr - e_prev) < options_.get_double("E_CONVERGENCE"));

// use residuals to get next amplitudes
#pragma omp parallel for schedule(static, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];
            for (int a = 0; a < n_pno_[ij]; ++a) {
                for (int b = 0; b < n_pno_[ij]; ++b) {
                    T_iajb_[ij]->add(a, b, -R_iajb[ij]->get(a, b) / ((e_pno_[ij]->get(a) + e_pno_[ij]->get(b)) -
                                                                    (F_lmo_->get(i, i) + F_lmo_->get(j, j))));
                }
            }
        }

        // DIIS extrapolation
        auto T_iajb_flat = flatten_mats(T_iajb_);
        auto R_iajb_flat = flatten_mats(R_iajb);

        if (iteration == 0) {
            diis.set_error_vector_size(R_iajb_flat.get());
            diis.set_vector_size(T_iajb_flat.get());
        }

        diis.add_entry(R_iajb_flat.get(), T_iajb_flat.get());
        diis.extrapolate(T_iajb_flat.get());

        copy_flat_mats(T_iajb_flat, T_iajb_);

#pragma omp parallel for schedule(static, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            Tt_iajb_[ij]->copy(T_iajb_[ij]);
            Tt_iajb_[ij]->scale(2.0);
            Tt_iajb_[ij]->subtract(T_iajb_[ij]->transpose());
        }

        outfile->Printf("  @LMP2 iter %3d: %16.12f %10.3e %10.3e\n", iteration, e_curr, e_curr - e_prev, r_curr);

        iteration++;

        if(iteration > max_iteration) {
            throw PSIEXCEPTION("Maximum DLPNO iterations exceeded.");
        }
    }

    e_lmp2_ = e_curr;

    e_lmp2_os_ = 0.0;
    for (int ij = 0; ij < T_iajb_.size(); ++ij) {
        e_lmp2_os_ += K_iajb_[ij]->vector_dot(T_iajb_[ij]);
    }
    e_lmp2_ss_ = e_curr - e_lmp2_os_;

}

double DLPNOMP2::compute_iteration_energy(const std::vector<SharedMatrix> &R_iajb) {
    double e_mp2 = 0.0;
    for (int ij = 0; ij < Tt_iajb_.size(); ++ij) {
        e_mp2 += K_iajb_[ij]->vector_dot(Tt_iajb_[ij]);
        e_mp2 += R_iajb[ij]->vector_dot(Tt_iajb_[ij]);
    }
    return e_mp2;
}

void DLPNOMP2::setup_orbitals() {
    int natom = molecule_->natom();
    int nbf = basisset_->nbf();
    int nshell = basisset_->nshell();
    int naux = ribasis_->nbf();
    int nshellri = ribasis_->nshell();
    int naocc = nalpha_ - nfrzc();

    auto C_occ = reference_wavefunction_->Ca_subset("AO", "OCC");

    timer_on("Local MOs");
    // Localize active occupied orbitals
    if (options_.get_str("DLPNO_LOCAL_ORBITALS") == "BOYS") {
        BoysLocalizer localizer = BoysLocalizer(basisset_, reference_wavefunction_->Ca_subset("AO", "ACTIVE_OCC"));
        localizer.set_convergence(options_.get_double("LOCAL_CONVERGENCE"));
        localizer.set_maxiter(options_.get_int("LOCAL_MAXITER"));
        localizer.localize();
        C_lmo_ = localizer.L();
    } else if (options_.get_str("DLPNO_LOCAL_ORBITALS") == "PIPEK_MEZEY") {
        PMLocalizer localizer = PMLocalizer(basisset_, reference_wavefunction_->Ca_subset("AO", "ACTIVE_OCC"));
        localizer.set_convergence(options_.get_double("LOCAL_CONVERGENCE"));
        localizer.set_maxiter(options_.get_int("LOCAL_MAXITER"));
        localizer.localize();
        C_lmo_ = localizer.L();
    } else {
        throw PSIEXCEPTION("Invalid option for DLPNO_LOCAL_ORBITALS");
    }
    timer_off("Local MOs");

    F_lmo_ = linalg::triplet(C_lmo_, reference_wavefunction_->Fa(), C_lmo_, true, false, false);

    timer_on("Projected AOs");

    // Form projected atomic orbitals by removing occupied space from the basis
    C_pao_ = std::make_shared<Matrix>("Projected Atomic Orbitals", nbf, nbf);
    C_pao_->identity();
    C_pao_->subtract(linalg::triplet(C_occ, C_occ, reference_wavefunction_->S(), false, true, false));
    S_pao_ = linalg::triplet(C_pao_, reference_wavefunction_->S(), C_pao_, true, false, false);

    // normalize PAOs
    for (size_t i = 0; i < C_pao_->colspi(0); ++i) {
        C_pao_->scale_column(0, i, pow(S_pao_->get(i, i), -0.5));
    }
    S_pao_ = linalg::triplet(C_pao_, reference_wavefunction_->S(), C_pao_, true, false, false);
    F_pao_ = linalg::triplet(C_pao_, reference_wavefunction_->Fa(), C_pao_, true, false, false);

    timer_off("Projected AOs");

    // map from atomic center to orbital/aux basis function/shell index

    atom_to_bf_.resize(natom);
    atom_to_ribf_.resize(natom);
    atom_to_shell_.resize(natom);
    atom_to_rishell_.resize(natom);

    for (size_t i = 0; i < nbf; ++i) {
        atom_to_bf_[basisset_->function_to_center(i)].push_back(i);
    }

    for (size_t i = 0; i < naux; ++i) {
        atom_to_ribf_[ribasis_->function_to_center(i)].push_back(i);
    }

    for (size_t s = 0; s < nshell; s++) {
        atom_to_shell_[basisset_->shell_to_center(s)].push_back(s);
    }

    for (size_t s = 0; s < nshellri; s++) {
        atom_to_rishell_[ribasis_->shell_to_center(s)].push_back(s);
    }
}

double DLPNOMP2::compute_energy() {
    timer_on("DLPNO-MP2");

    print_header();

    timer_on("Setup Orbitals");
    setup_orbitals();
    timer_off("Setup Orbitals");

    timer_on("Overlap Ints");
    compute_overlap_ints();
    timer_off("Overlap Ints");

    timer_on("Dipole Ints");
    compute_dipole_ints();
    timer_off("Dipole Ints");

    timer_on("Sparsity");
    prep_sparsity();
    timer_off("Sparsity");

    timer_on("DF Ints");
    compute_df_ints();
    timer_off("DF Ints");

    timer_on("PNO Transform");
    pno_transform();
    timer_off("PNO Transform");

    timer_on("PNO Overlaps");
    compute_pno_overlaps();
    timer_off("PNO Overlaps");

    timer_on("LMP2");
    lmp2_iterations();
    timer_off("LMP2");

    print_results();

    timer_off("DLPNO-MP2");

    double e_scf = reference_wavefunction_->energy();
    double e_mp2_corr = e_lmp2_ + de_dipole_ + de_pno_total_;
    double e_mp2_corr_os = e_lmp2_os_ + 0.5 * de_dipole_ + de_pno_total_os_;
    double e_mp2_corr_ss = e_lmp2_ss_ + 0.5 * de_dipole_ + de_pno_total_ss_;
    double e_mp2_total = e_scf + e_mp2_corr;

    double sss = options_.get_double("MP2_SS_SCALE");
    double oss = options_.get_double("MP2_OS_SCALE");

    set_scalar_variable("MP2 CORRELATION ENERGY", e_mp2_corr);
    set_scalar_variable("CURRENT CORRELATION ENERGY", e_mp2_corr);
    set_scalar_variable("MP2 TOTAL ENERGY", e_mp2_total);
    set_scalar_variable("CURRENT ENERGY", e_mp2_total);

    set_scalar_variable("MP2 SINGLES ENERGY", 0.0);
    set_scalar_variable("MP2 DOUBLES ENERGY", e_mp2_corr);
    set_scalar_variable("MP2 SAME-SPIN CORRELATION ENERGY", e_mp2_corr_ss);
    set_scalar_variable("MP2 OPPOSITE-SPIN CORRELATION ENERGY", e_mp2_corr_os);
    set_scalar_variable("SCS-MP2 CORRELATION ENERGY", 6.0/5.0 * e_mp2_corr_os +
                                                      1.0/3.0 * e_mp2_corr_ss);
    set_scalar_variable("SCS-MP2 TOTAL ENERGY", e_scf +
                                                6.0/5.0 * e_mp2_corr_os +
                                                1.0/3.0 * e_mp2_corr_ss);
    set_scalar_variable("CUSTOM SCS-MP2 CORRELATION ENERGY", oss * e_mp2_corr_os +
                                                             sss * e_mp2_corr_ss);
    set_scalar_variable("CUSTOM SCS-MP2 TOTAL ENERGY", e_scf +
                                                       oss * e_mp2_corr_os +
                                                       sss * e_mp2_corr_ss);

    if (options_.get_bool("DLPNO_CCSD")) {
        timer_on("DLPNO CCSD");

        timer_on("DLPNO CCSD : Integrals");

        // Compute Integrals
        estimate_memory();
        compute_qij();
        compute_qab();
        compute_cc_ints();

        qij_.clear();
        qia_.clear();
        if (virtual_storage_ != DIRECT) qab_.clear();

        timer_off("DLPNO CCSD : Integrals");

        compute_cc_pno_overlaps();

        // Run DLPNO-CCSD
        lccsd_iterations();

        double e_ccsd_corr = e_lccsd_ + de_dipole_ + de_pno_total_;
        double e_ccsd_total = e_scf + e_ccsd_corr;

        set_scalar_variable("CCSD CORRELATION ENERGY", e_ccsd_corr);
        set_scalar_variable("CURRENT CORRELATION ENERGY", e_ccsd_corr);
        set_scalar_variable("CCSD TOTAL ENERGY", e_ccsd_total);
        set_scalar_variable("CURRENT ENERGY", e_ccsd_total);

        print_ccsd_results();

        timer_off("DLPNO CCSD");

        if (options_.get_bool("DLPNO_CCSD_T")) {

            timer_on("DLPNO CCSD(T)");

            tno_transform();
            compute_pno_tno_overlaps();
            compute_tno_overlaps();
            compute_W_iajbkc();
            lccsd_t_iterations();

            double e_ccsd_t_corr = e_lccsd_t_ + de_dipole_ + de_pno_total_;
            double e_ccsd_t_total = e_scf + e_ccsd_t_corr;

            set_scalar_variable("CCSD(T) CORRELATION ENERGY", e_ccsd_t_corr);
            set_scalar_variable("CURRENT CORRELATION ENERGY", e_ccsd_t_corr);
            set_scalar_variable("CCSD TOTAL ENERGY", e_ccsd_t_total);
            set_scalar_variable("CURRENT ENERGY", e_ccsd_t_total);

            timer_off("DLPNO CCSD(T)");

            return e_ccsd_t_total;
        }

        return e_ccsd_total;
    }

    return e_mp2_total;
}

void DLPNOMP2::print_header() {
    outfile->Printf("   --------------------------------------------\n");
    outfile->Printf("                     DLPNO-MP2                 \n");
    outfile->Printf("                   by Zach Glick               \n");
    outfile->Printf("   --------------------------------------------\n\n");
    outfile->Printf("  DLPNO convergence set to %s.\n\n", options_.get_str("PNO_CONVERGENCE").c_str());
    outfile->Printf("  Detailed DLPNO thresholds and cutoffs:\n");
    outfile->Printf("    T_CUT_DO     = %6.3e \n", T_CUT_DO_);
    outfile->Printf("    T_CUT_PNO    = %6.3e \n", T_CUT_PNO_);
    outfile->Printf("    T_CUT_DO_ij  = %6.3e \n", options_.get_double("T_CUT_DO_ij"));
    outfile->Printf("    T_CUT_PRE    = %6.3e \n", options_.get_double("T_CUT_PRE"));
    outfile->Printf("    T_CUT_DO_PRE = %6.3e \n", options_.get_double("T_CUT_DO_PRE"));
    outfile->Printf("    T_CUT_MKN    = %6.3e \n", options_.get_double("T_CUT_MKN"));
    outfile->Printf("    T_CUT_CLMO   = %6.3e \n", options_.get_double("T_CUT_CLMO"));
    outfile->Printf("    T_CUT_CPAO   = %6.3e \n", options_.get_double("T_CUT_CPAO"));
    outfile->Printf("    S_CUT        = %6.3e \n", options_.get_double("S_CUT"));
    outfile->Printf("    F_CUT        = %6.3e \n", options_.get_double("F_CUT"));
    outfile->Printf("\n");
}

void DLPNOMP2::print_aux_domains() {
    size_t total_atoms = 0, min_atoms = lmo_to_riatoms_[0].size(), max_atoms = 0;
    for (const auto &atom_list : lmo_to_riatoms_) {
        total_atoms += atom_list.size();
        min_atoms = std::min(min_atoms, atom_list.size());
        max_atoms = std::max(max_atoms, atom_list.size());
    }

    size_t total_bfs = 0, min_bfs = lmo_to_ribfs_[0].size(), max_bfs = 0;
    for (const auto &bf_list : lmo_to_ribfs_) {
        total_bfs += bf_list.size();
        min_bfs = std::min(min_bfs, bf_list.size());
        max_bfs = std::max(max_bfs, bf_list.size());
    }

    size_t naocc = lmo_to_ribfs_.size();
    outfile->Printf("\n");
    outfile->Printf("    Auxiliary BFs per Local MO:\n");
    outfile->Printf("      Average = %4zu AUX BFs (%zu atoms)\n", total_bfs / naocc, total_atoms / naocc);
    outfile->Printf("      Min     = %4zu AUX BFs (%zu atoms)\n", min_bfs, min_atoms);
    outfile->Printf("      Max     = %4zu AUX BFs (%zu atoms)\n", max_bfs, max_atoms);
}

void DLPNOMP2::print_pao_domains() {
    size_t total_atoms = 0, min_atoms = lmo_to_paoatoms_[0].size(), max_atoms = 0;
    for (const auto &atom_list : lmo_to_paoatoms_) {
        total_atoms += atom_list.size();
        min_atoms = std::min(min_atoms, atom_list.size());
        max_atoms = std::max(max_atoms, atom_list.size());
    }

    size_t total_paos = 0, min_paos = lmo_to_paos_[0].size(), max_paos = 0;
    for (const auto &pao_list : lmo_to_paos_) {
        total_paos += pao_list.size();
        min_paos = std::min(min_paos, pao_list.size());
        max_paos = std::max(max_paos, pao_list.size());
    }

    size_t naocc = lmo_to_paos_.size();
    outfile->Printf("  \n");
    outfile->Printf("    Projected AOs per Local MO:\n");
    outfile->Printf("      Average = %4zu PAOs (%zu atoms)\n", total_paos / naocc, total_atoms / naocc);
    outfile->Printf("      Min     = %4zu PAOs (%zu atoms)\n", min_paos, min_atoms);
    outfile->Printf("      Max     = %4zu PAOs (%zu atoms)\n", max_paos, max_atoms);
}

void DLPNOMP2::print_lmo_domains() {
    int naocc = i_j_to_ij_.size();

    int exclude_pairs_overlap = 0;
    int exclude_pairs_energy = 0;
    int total_lmos = 0, min_lmos = naocc, max_lmos = 0;
    for (size_t i = 0; i < naocc; i++) {
        int lmos = 0;
        for (size_t j = 0; j < naocc; ++j) {
            if (i_j_to_ij_[i][j] != -1) {
                lmos += 1;
            }
            bool overlap_big = (DOI_ij_->get(i, j) > options_.get_double("T_CUT_DO_ij"));
            bool energy_big = (fabs(dipole_pair_e_bound_->get(i, j)) > options_.get_double("T_CUT_PRE"));
            if (!overlap_big) exclude_pairs_overlap++;
            if (!energy_big) exclude_pairs_energy++;
        }
        total_lmos += lmos;
        min_lmos = std::min(min_lmos, lmos);
        max_lmos = std::max(max_lmos, lmos);
    }

    outfile->Printf("\n");
    outfile->Printf("    Local MOs per Local MO:\n");
    outfile->Printf("      Average = %4d LMOs\n", total_lmos / naocc);
    outfile->Printf("      Min     = %4d LMOs\n", min_lmos);
    outfile->Printf("      Max     = %4d LMOs\n", max_lmos);
    outfile->Printf(" \n");
    outfile->Printf("    Screened %d of %d LMO pairs (%.2f %%)\n", naocc * naocc - total_lmos, naocc * naocc,
                    100.0 - (total_lmos * 100.0) / (naocc * naocc));
    outfile->Printf("             %d pairs met overlap criteria\n", exclude_pairs_overlap);
    outfile->Printf("             %d pairs met energy criteria\n", exclude_pairs_energy);
    outfile->Printf(" \n");
    outfile->Printf("    Screened LMO pair energy =  %.12f \n", de_dipole_);
}

void DLPNOMP2::print_aux_pair_domains() {
    int n_lmo_pairs = lmopair_to_ribfs_.size();
    int min_domain_ri = lmopair_to_ribfs_[0].size(), max_domain_ri = 0, total_domain_ri = 0;
    int min_domain_ri_atom = lmopair_to_riatoms_[0].size(), max_domain_ri_atom = 0, total_domain_ri_atom = 0;
    for (size_t ij = 0; ij < n_lmo_pairs; ij++) {
        int pair_domain_size_ri = lmopair_to_ribfs_[ij].size();
        int pair_domain_size_ri_atom = lmopair_to_riatoms_[ij].size();

        total_domain_ri += pair_domain_size_ri;
        total_domain_ri_atom += pair_domain_size_ri_atom;

        min_domain_ri = std::min(min_domain_ri, pair_domain_size_ri);
        min_domain_ri_atom = std::min(min_domain_ri_atom, pair_domain_size_ri_atom);

        max_domain_ri = std::max(max_domain_ri, pair_domain_size_ri);
        max_domain_ri_atom = std::max(max_domain_ri_atom, pair_domain_size_ri_atom);
    }

    outfile->Printf("  \n");
    outfile->Printf("    Auxiliary BFs per Local MO pair:\n");
    outfile->Printf("      Average = %4d AUX BFs (%d atoms)\n", total_domain_ri / n_lmo_pairs,
                    total_domain_ri_atom / n_lmo_pairs);
    outfile->Printf("      Min     = %4d AUX BFs (%d atoms)\n", min_domain_ri, min_domain_ri_atom);
    outfile->Printf("      Max     = %4d AUX BFs (%d atoms)\n", max_domain_ri, max_domain_ri_atom);
}

void DLPNOMP2::print_pao_pair_domains() {
    int n_lmo_pairs = lmopair_to_paos_.size();
    int min_domain_pao = lmopair_to_paos_[0].size(), max_domain_pao = 0, total_domain_pao = 0;
    int min_domain_atom = lmopair_to_paoatoms_[0].size(), max_domain_atom = 0, total_domain_atom = 0;
    for (size_t ij = 0; ij < n_lmo_pairs; ij++) {
        int pair_domain_size_pao = lmopair_to_paos_[ij].size();
        int pair_domain_size_atom = lmopair_to_paoatoms_[ij].size();

        total_domain_pao += pair_domain_size_pao;
        total_domain_atom += pair_domain_size_atom;

        min_domain_pao = std::min(min_domain_pao, pair_domain_size_pao);
        min_domain_atom = std::min(min_domain_atom, pair_domain_size_atom);

        max_domain_pao = std::max(max_domain_pao, pair_domain_size_pao);
        max_domain_atom = std::max(max_domain_atom, pair_domain_size_atom);
    }

    outfile->Printf("  \n");
    outfile->Printf("    Projected AOs per Local MO pair:\n");
    outfile->Printf("      Average = %4d PAOs (%d atoms)\n", total_domain_pao / n_lmo_pairs,
                    total_domain_atom / n_lmo_pairs);
    outfile->Printf("      Min     = %4d PAOs (%d atoms)\n", min_domain_pao, min_domain_atom);
    outfile->Printf("      Max     = %4d PAOs (%d atoms)\n", max_domain_pao, max_domain_atom);
}

void DLPNOMP2::print_integral_sparsity() {
    // statistics for number of (MN|K) shell triplets we need to compute

    int nbf = basisset_->nbf();
    int nshell = basisset_->nshell();
    int naux = ribasis_->nbf();
    int naocc = nalpha_ - nfrzc();

    size_t triplets = 0;       // computed (MN|K) triplets with no screening
    size_t triplets_lmo = 0;   // computed (MN|K) triplets with only LMO screening
    size_t triplets_pao = 0;   // computed (MN|K) triplets with only PAO screening
    size_t triplets_both = 0;  // computed (MN|K) triplets with LMO and PAO screening

    for (size_t atom = 0; atom < riatom_to_shells1_.size(); atom++) {
        size_t nshellri_atom = atom_to_rishell_[atom].size();
        triplets += nshell * nshell * nshellri_atom;
        triplets_lmo += riatom_to_shells1_[atom].size() * nshell * nshellri_atom;
        triplets_pao += nshell * riatom_to_shells2_[atom].size() * nshellri_atom;
        triplets_both += riatom_to_shells1_[atom].size() * riatom_to_shells2_[atom].size() * nshellri_atom;
    }
    size_t screened_total = triplets - triplets_both;
    size_t screened_lmo = triplets - triplets_lmo;
    size_t screened_pao = triplets - triplets_pao;

    // statistics for the number of (iu|Q) integrals we're left with after the transformation

    size_t total_integrals = (size_t)naocc * nbf * naux;
    size_t actual_integrals = 0;

    for (size_t atom = 0; atom < riatom_to_shells1_.size(); atom++) {
        actual_integrals +=
            riatom_to_lmos_ext_[atom].size() * riatom_to_paos_ext_[atom].size() * atom_to_ribf_[atom].size();
    }

    // number of doubles * (2^3 bytes / double) * (1 GiB / 2^30 bytes)
    double total_memory = total_integrals * pow(2.0, -27);
    double actual_memory = actual_integrals * pow(2.0, -27);
    double screened_memory = total_memory - actual_memory;

    outfile->Printf("\n");
    outfile->Printf("    Coefficient sparsity in AO -> LMO transform: %6.2f %% \n", screened_lmo * 100.0 / triplets);
    outfile->Printf("    Coefficient sparsity in AO -> PAO transform: %6.2f %% \n", screened_pao * 100.0 / triplets);
    outfile->Printf("    Coefficient sparsity in combined transforms: %6.2f %% \n", screened_total * 100.0 / triplets);
    outfile->Printf("\n");
    outfile->Printf("    Storing transformed LMO/PAO integrals in sparse format.\n");
    outfile->Printf("    Required memory: %.3f GiB (%.2f %% reduction from dense format) \n", actual_memory,
                    screened_memory * 100.0 / total_memory);
}

void DLPNOMP2::print_results() {
    outfile->Printf("  \n");
    outfile->Printf("  Total DLPNO-MP2 Correlation Energy: %16.12f \n", e_lmp2_ + de_pno_total_ + de_dipole_);
    outfile->Printf("    MP2 Correlation Energy:           %16.12f \n", e_lmp2_);
    outfile->Printf("    LMO Truncation Correction:        %16.12f \n", de_dipole_);
    outfile->Printf("    PNO Truncation Correction:        %16.12f \n", de_pno_total_);
}

void DLPNOMP2::print_ccsd_results() {
    outfile->Printf("  \n");
    outfile->Printf("  Total DLPNO-CCSD Correlation Energy: %16.12f \n", e_lccsd_ + de_pno_total_ + de_dipole_);
    outfile->Printf("    CCSD Correlation Energy:           %16.12f \n", e_lccsd_);
    outfile->Printf("    LMO Truncation Correction:         %16.12f \n", de_dipole_);
    outfile->Printf("    PNO Truncation Correction:         %16.12f \n", de_pno_total_);
    outfile->Printf("    Get Diagonalized!!! --Andy Jiang, February 2023 \n\n\n");
    outfile->Printf("  @Total DLPNO-CCSD Energy: %16.12f \n", variables_["SCF TOTAL ENERGY"] + e_lccsd_ + de_pno_total_ + de_dipole_);
}

void DLPNOMP2::store_information() {
    // Store LMO information
    lmo_matrices_["C_LMO"] = C_lmo_;
    lmo_matrices_["F_LMO"] = F_lmo_;

    // Store PAO information
    pao_matrices_["C_PAO"] = C_pao_;
    pao_matrices_["F_PAO"] = F_pao_;
    pao_matrices_["S_PAO"] = S_pao_;

    // Store PNO information
    pno_matrices_["K_IAJB"] = K_iajb_;
    pno_matrices_["T_IAJB"] = T_iajb_;
    pno_matrices_["Tt_IAJB"] = Tt_iajb_;
    pno_matrices_["X_PNO"] = X_pno_;

    // Stores SparseMap information
    sparse_maps_["ATOM_TO_BF"] = atom_to_bf_;
    sparse_maps_["ATOM_TO_RIBF"] = atom_to_ribf_;
    sparse_maps_["ATOM_TO_SHELL"] = atom_to_shell_;
    sparse_maps_["ATOM_TO_RISHELL"] = atom_to_rishell_;

    sparse_maps_["LMO_TO_BF"] = lmo_to_bfs_;
    sparse_maps_["LMO_TO_ATOM"] = lmo_to_atoms_;
    sparse_maps_["PAO_TO_BF"] = pao_to_bfs_;
    sparse_maps_["PAO_TO_ATOM"] = pao_to_atoms_;

    sparse_maps_["LMO_TO_RIBF"] = lmo_to_ribfs_;
    sparse_maps_["LMO_TO_RIATOM"] = lmo_to_riatoms_;
    sparse_maps_["LMO_TO_PAO"] = lmo_to_paos_;
    sparse_maps_["LMO_TO_PAOATOM"] = lmo_to_paoatoms_;

    sparse_maps_["LMOPAIR_TO_RIBF"] = lmopair_to_ribfs_;
    sparse_maps_["LMOPAIR_TO_RIATOM"] = lmopair_to_riatoms_;
    sparse_maps_["LMOPAIR_TO_PAO"] = lmopair_to_paos_;
    sparse_maps_["LMOPAIR_TO_PAOATOM"] = lmopair_to_paoatoms_;
    sparse_maps_["LMOPAIR_TO_LMO"] = lmopair_to_lmos_;

    sparse_maps_["LMO_TO_RIATOM_EXT"] = lmo_to_riatoms_ext_;
    sparse_maps_["RIATOM_TO_LMO_EXT"] = riatom_to_lmos_ext_;
    sparse_maps_["RIATOM_TO_PAO_EXT"] = riatom_to_paos_ext_;
    sparse_maps_["RIATOM_TO_ATOM_1"] = riatom_to_atoms1_;
    sparse_maps_["RIATOM_TO_SHELL_1"] = riatom_to_shells1_;
    sparse_maps_["RIATOM_TO_BF_1"] = riatom_to_bfs1_;
    sparse_maps_["RIATOM_TO_ATOM_2"] = riatom_to_atoms2_;
    sparse_maps_["RIATOM_TO_SHELL_2"] = riatom_to_shells2_;
    sparse_maps_["RIATOM_TO_BF_2"] = riatom_to_bfs2_;

    sparse_maps_["RIATOM_TO_LMO_EXT_DENSE"] = riatom_to_lmos_ext_dense_;
    sparse_maps_["RIATOM_TO_PAO_EXT_DENSE"] = riatom_to_paos_ext_dense_;
    sparse_maps_["LMOPAIR_TO_LMO_DENSE"] = lmopair_to_lmos_dense_;
}

void DLPNOMP2::compute_qij() {
    timer_on("(mn|K)->(ij|K)");

    outfile->Printf("   Computing qij...\n\n");

    int nbf = basisset_->nbf();
    int naux = ribasis_->nbf();

    size_t nthread = 1;
#ifdef _OPENMP
    nthread = omp_get_max_threads();
#endif

    std::shared_ptr<IntegralFactory> factory =
        std::make_shared<IntegralFactory>(ribasis_, BasisSet::zero_ao_basis_set(), basisset_, basisset_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eris(nthread);

    eris[0] = std::shared_ptr<TwoBodyAOInt>(factory->eri());
    for (size_t thread = 1; thread < nthread; thread++) {
        eris[thread] = std::shared_ptr<TwoBodyAOInt>(eris.front()->clone());
    }

    auto SC_lmo = linalg::doublet(reference_wavefunction_->S(), C_lmo_, false, false);

    qij_.resize(naux);

    // LMO-LMO DF ints
#pragma omp parallel for schedule(dynamic, 1)
    for (int Q = 0; Q < ribasis_->nshell(); Q++) {
        int nq = ribasis_->shell(Q).nfunction();
        int qstart = ribasis_->shell(Q).function_index();
        int centerQ = ribasis_->shell_to_center(Q);

        size_t thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // Sparse lists for LMO/LMO
        auto bf_map_lmo = riatom_to_bfs1_[centerQ];

        // inverse map, from global basis function index to auxiliary specific index
        std::vector<int> bf_map_lmo_inv(nbf, -1);
        for (int m_ind = 0; m_ind < bf_map_lmo.size(); m_ind++) {
            bf_map_lmo_inv[bf_map_lmo[m_ind]] = m_ind;
        }

        for (size_t q = 0; q < nq; q++) {
            qij_[qstart + q] = std::make_shared<Matrix>("(mn|Q)", bf_map_lmo.size(), bf_map_lmo.size());
        }

        for (int M : riatom_to_shells1_[centerQ]) {
            int nm = basisset_->shell(M).nfunction();
            int mstart = basisset_->shell(M).function_index();
            int centerM = basisset_->shell_to_center(M);

            for (int N : riatom_to_shells1_[centerQ]) {
                int nn = basisset_->shell(N).nfunction();
                int nstart = basisset_->shell(N).function_index();
                int centerN = basisset_->shell_to_center(N);

                // TODO: Permutational Symmetry
                eris[thread]->compute_shell(Q, 0, M, N);
                const double* buffer = eris[thread]->buffer();

                for (int q = 0; q < nq; q++) {
                    for (int m = 0; m < nm; m++) {
                        for (int n = 0; n < nn; n++) {
                            int index_m = bf_map_lmo_inv[mstart + m];
                            int index_n = bf_map_lmo_inv[nstart + n];
                            qij_[qstart + q]->set(index_m, index_n, *(buffer));
                            buffer++;
                        }
                    }
                }

            } // N loop
        } // M loop

        auto C_lmo_slice = submatrix_rows_and_cols(*SC_lmo, riatom_to_bfs1_[centerQ], riatom_to_lmos_ext_[centerQ]);
        auto S_aa = submatrix_rows_and_cols(*reference_wavefunction_->S(), riatom_to_bfs1_[centerQ], riatom_to_bfs1_[centerQ]);
        C_DGESV_wrapper(S_aa, C_lmo_slice);

        // (mn|Q) C_mi C_nj ->(ij|Q)
        for (size_t q = 0; q < nq; q++) {
            qij_[qstart+q] = linalg::triplet(C_lmo_slice, qij_[qstart+q], C_lmo_slice, true, false, false);
        }
    }

    timer_off("(mn|K)->(ij|K)");
}

void DLPNOMP2::compute_qab() {
    timer_on("(mn|K)->(ab|K)");

    outfile->Printf("   Computing qab...\n\n");

    int nbf = basisset_->nbf();
    int naux = ribasis_->nbf();

    size_t nthread = 1;
#ifdef _OPENMP
    nthread = omp_get_max_threads();
#endif

    std::shared_ptr<IntegralFactory> factory =
        std::make_shared<IntegralFactory>(ribasis_, BasisSet::zero_ao_basis_set(), basisset_, basisset_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eris(nthread);

    eris[0] = std::shared_ptr<TwoBodyAOInt>(factory->eri());
    for (size_t thread = 1; thread < nthread; thread++) {
        eris[thread] = std::shared_ptr<TwoBodyAOInt>(eris.front()->clone());
    }

    qab_.resize(naux);

    // PAO-PAO DF ints
#pragma omp parallel for schedule(dynamic, 1)
    for (int Q = 0; Q < ribasis_->nshell(); Q++) {
        int nq = ribasis_->shell(Q).nfunction();
        int qstart = ribasis_->shell(Q).function_index();
        int centerQ = ribasis_->shell_to_center(Q);

        size_t thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // Sparse lists for PAO/PAO
        auto bf_map_pao = riatom_to_bfs2_[centerQ];

        // inverse map, from global basis function index to auxiliary specific index
        std::vector<int> bf_map_pao_inv(nbf, -1);
        for (int m_ind = 0; m_ind < bf_map_pao.size(); m_ind++) {
            bf_map_pao_inv[bf_map_pao[m_ind]] = m_ind;
        }

        for (size_t q = 0; q < nq; q++) {
            qab_[qstart + q] = std::make_shared<Matrix>("(mn|Q)", bf_map_pao.size(), bf_map_pao.size());
        }

        for (int M : riatom_to_shells2_[centerQ]) {
            int nm = basisset_->shell(M).nfunction();
            int mstart = basisset_->shell(M).function_index();
            int centerM = basisset_->shell_to_center(M);

            for (int N : riatom_to_shells2_[centerQ]) {
                int nn = basisset_->shell(N).nfunction();
                int nstart = basisset_->shell(N).function_index();
                int centerN = basisset_->shell_to_center(N);

                // TODO: Permutational Symmetry
                eris[thread]->compute_shell(Q, 0, M, N);
                const double* buffer = eris[thread]->buffer();

                for (int q = 0; q < nq; q++) {
                    for (int m = 0; m < nm; m++) {
                        for (int n = 0; n < nn; n++) {
                            int index_m = bf_map_pao_inv[mstart + m];
                            int index_n = bf_map_pao_inv[nstart + n];
                            qab_[qstart + q]->set(index_m, index_n, *(buffer));
                            buffer++;
                        }
                    }
                }

            }  // N loop
        } // M loop

        auto C_pao_slice = submatrix_rows_and_cols(*C_pao_, riatom_to_bfs2_[centerQ], riatom_to_paos_ext_[centerQ]);

        // (mn|Q) C_mi C_nj ->(ij|Q)
        for (size_t q = 0; q < nq; q++) {
            qab_[qstart + q] = linalg::triplet(C_pao_slice, qab_[qstart + q], C_pao_slice, true, false, false);
        }
    }

    timer_off("(mn|K)->(ab|K)");
}

void DLPNOMP2::estimate_memory() {
    outfile->Printf("   => DLPNO-CCSD Memory Estimate <= \n\n");

    int nbf = basisset_->nbf();
    int naux = ribasis_->nbf();
    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();
    int npao = C_pao_->colspi(0);

    size_t qoo_lmo = 0L;
    size_t qov_pao = 0L;
    size_t qvv_pao = 0L;
    for (int q = 0; q < naux; q++) {
        int centerq = ribasis_->function_to_center(q);
        int nlmo_q = riatom_to_lmos_ext_[centerq].size();
        int npao_q = riatom_to_paos_ext_[centerq].size();

        qoo_lmo += nlmo_q * nlmo_q;
        qov_pao += nlmo_q * npao_q;
        qvv_pao += npao_q * npao_q;
    }

    size_t s_oooo = 0L;
    size_t oooo = 0L;
    size_t ooov = 0L;
    size_t oovv = 0L;
    size_t ovvv = 0L;
    size_t qvv = 0L;
    size_t vvvv = 0L;

    for (int ij = 0; ij < n_lmo_pairs; ij++) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        int naux_ij = lmopair_to_ribfs_[ij].size();
        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int npno_ij = n_pno_[ij];

        oooo += nlmo_ij * nlmo_ij;
        ooov += nlmo_ij * npno_ij;
        oovv += 4 * npno_ij * npno_ij;
        ovvv += npno_ij * npno_ij * npno_ij;
        if (i <= j) qvv += naux_ij * npno_ij * npno_ij;
        if (i <= j) vvvv += npno_ij * npno_ij * npno_ij * npno_ij;

        if (i >= j) {
            for (int mn_ij = 0; mn_ij < nlmo_ij * nlmo_ij; mn_ij++) {
                int m_ij = mn_ij / nlmo_ij, n_ij = mn_ij % nlmo_ij;
                int m = lmopair_to_lmos_[ij][m_ij], n = lmopair_to_lmos_[ij][n_ij];
                int mn = i_j_to_ij_[m][n];

                int npno_mn = n_pno_[mn];
                if (m_ij >= n_ij) {
                    s_oooo += npno_ij * npno_mn;
                }
            }
        }
    }

    size_t direct_memory = (s_oooo + qoo_lmo + qov_pao + qvv_pao + oooo + ooov + oovv + ovvv) * sizeof(double);
    size_t qvv_memory = direct_memory + qvv * sizeof(double);
    size_t vvvv_memory = direct_memory + vvvv * sizeof(double);

    outfile->Printf("    Memory Required to Store Each Integral Type:\n");
    outfile->Printf("    S_pno   [ij/kl]   : %8.4f [GiB]\n", s_oooo * sizeof(double) / (1024 * 1024 * 1024.0));
    outfile->Printf("    (q|oo)  [LMO/LMO] : %8.4f [GiB]\n", qoo_lmo * sizeof(double) / (1024 * 1024 * 1024.0));
    outfile->Printf("    (q|ov)  [LMO/PAO] : %8.4f [GiB]\n", qov_pao * sizeof(double) / (1024 * 1024 * 1024.0));
    outfile->Printf("    (q|vv)  [PAO/PAO] : %8.4f [GiB]\n", qvv_pao * sizeof(double) / (1024 * 1024 * 1024.0));
    outfile->Printf("    (oo|oo) [Pair ij] : %8.4f [GiB]\n", oooo * sizeof(double) / (1024 * 1024 * 1024.0));
    outfile->Printf("    (oo|ov) [Pair ij] : %8.4f [GiB]\n", ooov * sizeof(double) / (1024 * 1024 * 1024.0));
    outfile->Printf("    (oo|vv) [Pair ij] : %8.4f [GiB]\n", oovv * sizeof(double) / (1024 * 1024 * 1024.0));
    outfile->Printf("    (ov|vv) [Pair ij] : %8.4f [GiB]\n", ovvv * sizeof(double) / (1024 * 1024 * 1024.0));
    outfile->Printf("    (q|vv)  [Pair ij] : %8.4f [GiB]\n", qvv * sizeof(double) / (1024 * 1024 * 1024.0));
    outfile->Printf("    (vv|vv) [Pair ij] : %8.4f [GiB]\n", vvvv * sizeof(double) / (1024 * 1024 * 1024.0));
    outfile->Printf("    DIRECT Mem Req    : %8.4f [GiB]\n", direct_memory / (1024 * 1024 * 1024.0));
    outfile->Printf("    DF-STORE Mem Req  : %8.4f [GiB]\n", qvv_memory / (1024 * 1024 * 1024.0));
    outfile->Printf("    STORE Mem Req     : %8.4f [GiB]\n", vvvv_memory / (1024 * 1024 * 1024.0));
    outfile->Printf("    Memory Given      : %8.4f [GiB]\n", memory_ / (1024 * 1024 * 1024.0));
    outfile->Printf("    Using 80% of Mem  : %8.4f [GiB]\n\n", 0.8 * memory_ / (1024 * 1024 * 1024.0));

    /*
    if (vvvv_memory < 0.8 * memory_) {
        outfile->Printf("   Storing 4-virtual integrals [HIGH MEMORY]...\n\n");
        virtual_storage_ = STORE;
    } else 
    */
    
    if (qvv_memory < 0.8 * memory_) {
        outfile->Printf("   Storing DF virtual/virtual integrals [MED MEMORY]...\n\n");
        virtual_storage_ = DF_STORE;
    } else if (direct_memory < 0.8 * memory_) {
        outfile->Printf("   Computing 4-virtual integrals as needed [LOW MEMORY]...\n\n");
        virtual_storage_ = DIRECT;
    } else {
        throw PSIEXCEPTION("Not enough memory given to complete DLPNO-CCSD computation!");
    }
}

void DLPNOMP2::compute_cc_pno_overlaps() {
    timer_on("Compute CC PNO overlaps");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    S_pno_ij_mn_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        int npno_ij = n_pno_[ij];
        int nlmo_ij = lmopair_to_lmos_[ij].size();

        if (npno_ij == 0 || i < j) continue;

        S_pno_ij_mn_[ij].resize(nlmo_ij * nlmo_ij);

        for (int mn_ij = 0; mn_ij < nlmo_ij * nlmo_ij; mn_ij++) {
            int m_ij = mn_ij / nlmo_ij, n_ij = mn_ij % nlmo_ij;
            int m = lmopair_to_lmos_[ij][m_ij], n = lmopair_to_lmos_[ij][n_ij];
            int mn = i_j_to_ij_[m][n];

            int npno_mn = n_pno_[mn];
            if (npno_mn == 0 || m_ij < n_ij) continue;

            S_pno_ij_mn_[ij][mn_ij] = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[mn]);
            S_pno_ij_mn_[ij][mn_ij] = linalg::triplet(X_pno_[ij], S_pno_ij_mn_[ij][mn_ij], X_pno_[mn], true, false, false);
        }
    }

    timer_off("Compute CC PNO overlaps");
}

inline SharedMatrix DLPNOMP2::S_PNO(int ij, int mn) {
    int i, j, m, n;
    std::tie(i, j) = ij_to_i_j_[ij];
    std::tie(m, n) = ij_to_i_j_[mn];
    
    if (i == m) {
        return S_pno_ij_ik_[ij][n];
    } else if (i == n) {
        return S_pno_ij_ik_[ij][m];
    } else if (j == m) {
        return S_pno_ij_kj_[ij][n];
    } else if (j == n) {
        return S_pno_ij_kj_[ij][m];
    } else {
        int m_ij = lmopair_to_lmos_dense_[ij][m], n_ij = lmopair_to_lmos_dense_[ij][n];
        if (m_ij == -1 || n_ij == -1) {
            outfile->Printf("Invalid PNO Pairs (%d, %d) and (%d, %d)\n", i, j, m, n);
            throw PSIEXCEPTION("Invalid PNO pairs!");
        }

        int nlmo_ij = lmopair_to_lmos_[ij].size();

        int mn_ij; 
        if (m_ij < n_ij) {
            mn_ij = n_ij * nlmo_ij + m_ij;
        } else {
            mn_ij = m_ij * nlmo_ij + n_ij;
        }

        if (i < j) {
            int ji = ij_to_ji_[ij];
            return S_pno_ij_mn_[ji][mn_ij];
        } else {
            return S_pno_ij_mn_[ij][mn_ij];
        }
    }
}

void DLPNOMP2::compute_cc_ints() {
    timer_on("Compute CC Ints");

    outfile->Printf("   Computing CC integrals...\n\n");

    int nbf = basisset_->nbf();
    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();
    int npao = C_pao_->colspi(0);

    K_mnij_.resize(n_lmo_pairs);
    K_mbij_.resize(n_lmo_pairs);
    J_ijab_.resize(n_lmo_pairs);
    K_maef_.resize(n_lmo_pairs);
    if (virtual_storage_ == STORE) K_abef_.resize(n_lmo_pairs);
    if (virtual_storage_ == DF_STORE) Qab_ij_.resize(n_lmo_pairs);
    L_iajb_.resize(n_lmo_pairs);
    Lt_iajb_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij];

        // number of PNOs in the pair domain
        int npno_ij = n_pno_[ij];
        if (npno_ij == 0) continue;

        // number of LMOs in the pair domain
        int nlmo_ij = lmopair_to_lmos_[ij].size();
        // number of PAOs in the pair domain (before removing linear dependencies)
        int npao_ij = lmopair_to_paos_[ij].size();
        // number of auxiliary functions in the pair domain
        int naux_ij = lmopair_to_ribfs_[ij].size();

        auto q_pair = std::make_shared<Matrix>(naux_ij, 1);

        auto q_io = std::make_shared<Matrix>(naux_ij, nlmo_ij);
        auto q_jo = std::make_shared<Matrix>(naux_ij, nlmo_ij);

        auto q_jv = std::make_shared<Matrix>(naux_ij, npno_ij);
        auto q_vv = std::make_shared<Matrix>(naux_ij, npno_ij * npno_ij);

        for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
            int q = lmopair_to_ribfs_[ij][q_ij];
            int centerq = ribasis_->function_to_center(q);

            int i_sparse = riatom_to_lmos_ext_dense_[centerq][i];
            int j_sparse = riatom_to_lmos_ext_dense_[centerq][j];
            std::vector<int> i_slice(1, i_sparse);
            std::vector<int> j_slice(1, j_sparse);

            q_pair->set(q_ij, 0, (*qij_[q])(i_sparse, j_sparse));
            
            auto q_io_tmp = submatrix_rows_and_cols(*qij_[q], i_slice, 
                                lmopair_lmo_to_riatom_lmo_[ij][q_ij]);
            C_DCOPY(nlmo_ij, &(*q_io_tmp)(0,0), 1, &(*q_io)(q_ij, 0), 1);

            auto q_jo_tmp = submatrix_rows_and_cols(*qij_[q], j_slice, 
                                lmopair_lmo_to_riatom_lmo_[ij][q_ij]);
            C_DCOPY(nlmo_ij, &(*q_jo_tmp)(0,0), 1, &(*q_jo)(q_ij, 0), 1);

            auto q_jv_tmp = submatrix_rows_and_cols(*qia_[q], j_slice,
                                lmopair_pao_to_riatom_pao_[ij][q_ij]);
            q_jv_tmp = linalg::doublet(q_jv_tmp, X_pno_[ij], false, false);
            C_DCOPY(npno_ij, &(*q_jv_tmp)(0,0), 1, &(*q_jv)(q_ij, 0), 1);

            auto q_vv_tmp = submatrix_rows_and_cols(*qab_[q], lmopair_pao_to_riatom_pao_[ij][q_ij],
                                lmopair_pao_to_riatom_pao_[ij][q_ij]);
            q_vv_tmp = linalg::triplet(X_pno_[ij], q_vv_tmp, X_pno_[ij], true, false, false);
            C_DCOPY(npno_ij * npno_ij, &(*q_vv_tmp)(0,0), 1, &(*q_vv)(q_ij, 0), 1);
        }

        auto A_solve = submatrix_rows_and_cols(*full_metric_, lmopair_to_ribfs_[ij], lmopair_to_ribfs_[ij]);
        A_solve->power(0.5, 1.0e-14);

        C_DGESV_wrapper(A_solve->clone(), q_pair);
        C_DGESV_wrapper(A_solve->clone(), q_io);
        C_DGESV_wrapper(A_solve->clone(), q_jo);
        C_DGESV_wrapper(A_solve->clone(), q_jv);
        C_DGESV_wrapper(A_solve, q_vv);

        K_mnij_[ij] = linalg::doublet(q_io, q_jo, true, false);
        K_mbij_[ij] = linalg::doublet(q_io, q_jv, true, false);
        J_ijab_[ij] = linalg::doublet(q_pair, q_vv, true, false);
        J_ijab_[ij]->reshape(npno_ij, npno_ij);

        K_maef_[ji] = linalg::doublet(q_jv, q_vv, true, false);
        auto K_maef_tmp = K_maef_[ji]->clone();
        for (int a_ij = 0; a_ij < npno_ij; a_ij++) {
            for (int e_ij = 0; e_ij < npno_ij; e_ij++) {
                for (int f_ij = 0; f_ij < npno_ij; f_ij++) {
                    (*K_maef_[ji])(a_ij, e_ij * npno_ij + f_ij) = 
                                    (*K_maef_tmp)(e_ij, a_ij * npno_ij + f_ij);
                }
            }
        }
        
        if (virtual_storage_ == STORE && i <= j) {
            K_abef_[ij] = linalg::doublet(q_vv, q_vv, true, false);
            auto K_abef_tmp = K_abef_[ij]->clone();
            for (int a_ij = 0; a_ij < npno_ij; a_ij++) {
                for (int b_ij = 0; b_ij < npno_ij; b_ij++) {
                    for (int e_ij = 0; e_ij < npno_ij; e_ij++) {
                        for (int f_ij = 0; f_ij < npno_ij; f_ij++) {
                            (*K_abef_[ij])(a_ij * npno_ij + b_ij, e_ij * npno_ij + f_ij) = 
                                            (*K_abef_tmp)(a_ij * npno_ij + e_ij, b_ij * npno_ij + f_ij);
                        }
                    }
                }
            }
        }

        if (virtual_storage_ == DF_STORE && i <= j) {
            Qab_ij_[ij].resize(naux_ij);
            for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
                Qab_ij_[ij][q_ij] = std::make_shared<Matrix>(npno_ij, npno_ij);
                C_DCOPY(npno_ij * npno_ij, &(*q_vv)(q_ij, 0), 1, &(*Qab_ij_[ij][q_ij])(0,0), 1);
            }
        }

        // L_iajb
        L_iajb_[ij] = K_iajb_[ij]->clone();
        L_iajb_[ij]->scale(2.0);
        L_iajb_[ij]->subtract(K_iajb_[ij]->transpose());

        // Lt_iajb
        Lt_iajb_[ij] = K_iajb_[ij]->clone();
        Lt_iajb_[ij]->scale(2.0);
        Lt_iajb_[ij]->subtract(J_ijab_[ij]);
    }

    timer_off("Compute CC Ints");
}

void DLPNOMP2::tno_transform() {
    timer_on("TNO transform");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();
    int npao = C_pao_->colspi(0);

    int ijk = 0;
    for (int ij = 0; ij < n_lmo_pairs; ij++) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        for (int k_ij : lmopair_to_lmos_[ij]) {
            int k = lmopair_to_lmos_[ij][k_ij];
            ijk_to_i_j_k_.push_back(std::make_tuple(i, j, k));
            i_j_k_to_ijk_[i * naocc * naocc + j * naocc + k] = ijk;
            ijk++;
        }
    }

    int n_lmo_triplets = ijk_to_i_j_k_.size();
    lmotriplet_to_paos_.resize(n_lmo_triplets);
    X_tno_.resize(n_lmo_triplets);
    e_tno_.resize(n_lmo_triplets);
    n_tno_.resize(n_lmo_triplets);

    std::vector<SharedMatrix> D_ij(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        if (i > j) continue;

        D_ij[ij] = linalg::doublet(Tt_iajb_[ij], T_iajb_[ij], false, true);
        D_ij[ij]->add(linalg::doublet(Tt_iajb_[ij], T_iajb_[ij], true, false));

        D_ij[ij] = linalg::triplet(X_pno_[ij], D_ij[ij], X_pno_[ij], false, false, true);

        if (i < j) {
            int ji = ij_to_ji_[ij];
            D_ij[ji] = D_ij[ij]->clone();
        }
    }

    std::vector<std::vector<int>> global_pao_to_pao_ij(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        global_pao_to_pao_ij[ij] = std::vector<int>(npao, -1);

        for (int u_ij = 0; u_ij < lmopair_to_paos_[ij].size(); u_ij++) {
            int u = lmopair_to_paos_[ij][u_ij];
            global_pao_to_pao_ij[ij][u] = u_ij;
        }
    }

#pragma omp parallel for schedule(static, 1)
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        int i, j, k;
        std::tie(i, j, k) = ijk_to_i_j_k_[ijk];
        int ij = i_j_to_ij_[i][j], jk = i_j_to_ij_[j][k], ik = i_j_to_ij_[i][k];

        lmotriplet_to_paos_[ijk] = merge_lists(lmopair_to_paos_[ij], lmo_to_paos_[k]);

        if (i > j || i > k || j > k) continue;

        // number of PAOs in the triplet domain (before removing linear dependencies)
        int npao_ijk = lmotriplet_to_paos_[ijk].size();

        // Form the triplet density (from pair densities in redundant basis)
        auto D_ijk = std::make_shared<Matrix>("D_ijk", npao_ijk, npao_ijk);
        D_ijk->zero();

        for (int u_ijk = 0; u_ijk < lmotriplet_to_paos_[ijk].size(); u_ijk++) {
            int u = lmotriplet_to_paos_[ijk][u_ijk];
            int u_ij = global_pao_to_pao_ij[ij][u], u_jk = global_pao_to_pao_ij[jk][u], u_ik = global_pao_to_pao_ij[ik][u];
            
            for (int v_ijk = 0; v_ijk < lmotriplet_to_paos_[ijk].size(); v_ijk++) {
                int v = lmotriplet_to_paos_[ijk][v_ijk];
                int v_ij = global_pao_to_pao_ij[ij][v], v_jk = global_pao_to_pao_ij[jk][v], v_ik = global_pao_to_pao_ij[ik][v];

                if (u_ij != -1 && v_ij != -1) (*D_ijk)(u_ijk, v_ijk) += (*D_ij[ij])(u_ij, v_ij);
                if (u_jk != -1 && v_jk != -1) (*D_ijk)(u_ijk, v_ijk) += (*D_ij[jk])(u_jk, v_jk);
                if (u_ik != -1 && v_ik != -1) (*D_ijk)(u_ijk, v_ijk) += (*D_ij[ik])(u_ik, v_ik);
            }
        }
        D_ijk->scale(1.0 / 3.0);

        // Canonicalize PAOs of triplet ijk
        auto S_pao_ijk = submatrix_rows_and_cols(*S_pao_, lmotriplet_to_paos_[ijk], lmotriplet_to_paos_[ijk]);
        auto F_pao_ijk = submatrix_rows_and_cols(*F_pao_, lmotriplet_to_paos_[ijk], lmotriplet_to_paos_[ijk]);

        SharedMatrix X_pao_ijk;
        SharedVector e_pao_ijk;
        std::tie(X_pao_ijk, e_pao_ijk) = orthocanonicalizer(S_pao_ijk, F_pao_ijk);

        F_pao_ijk = linalg::triplet(X_pao_ijk, F_pao_ijk, X_pao_ijk, true, false, false);
        D_ijk = linalg::triplet(X_pao_ijk, D_ijk, X_pao_ijk, true, false, false);

        size_t nvir_ijk = F_pao_ijk->rowspi(0);

        // Diagonalization of triplet density gives TNOs (in basis of LMO's virtual domain)
        // as well as TNO occ numbers
        auto X_tno_ijk = std::make_shared<Matrix>("eigenvectors", nvir_ijk, nvir_ijk);
        Vector tno_occ("eigenvalues", nvir_ijk);
        D_ijk->diagonalize(*X_tno_ijk, tno_occ, descending);

        int nvir_ijk_final = 0;
        for (size_t a = 0; a < nvir_ijk; ++a) {
            // TODO: Introduce a T_CUT_TNO parameter
            if (fabs(tno_occ.get(a)) >= T_CUT_PNO_) {
                nvir_ijk_final++;
            }
        }

        Dimension zero(1);
        Dimension dim_final(1);
        dim_final.fill(nvir_ijk_final);

        // This transformation gives orbitals that are orthonormal but not canonical
        X_tno_ijk = X_tno_ijk->get_block({zero, X_tno_ijk->rowspi()}, {zero, dim_final});
        tno_occ = tno_occ.get_block({zero, dim_final});

        SharedMatrix tno_canon;
        SharedVector e_tno_ijk;
        std::tie(tno_canon, e_tno_ijk) = canonicalizer(X_tno_ijk, F_pao_ijk);

        X_tno_ijk = linalg::doublet(X_tno_ijk, tno_canon, false, false);
        X_tno_ijk = linalg::doublet(X_pao_ijk, X_tno_ijk, false, false);

        X_tno_[ijk] = X_tno_ijk;
        e_tno_[ijk] = e_tno_ijk;
        n_tno_[ijk] = X_tno_ijk->colspi(0);

        // account for symmetry
        if (i != j || j != k || i != k) {
            int ikj = i_j_k_to_ijk_[i * naocc * naocc + k * naocc + j];
            int jik = i_j_k_to_ijk_[j * naocc * naocc + i * naocc + k];
            int jki = i_j_k_to_ijk_[j * naocc * naocc + k * naocc + i];
            int kij = i_j_k_to_ijk_[k * naocc * naocc + i * naocc + j];
            int kji = i_j_k_to_ijk_[k * naocc * naocc + j * naocc + i];

            std::vector<int> perms{ikj, jik, jki, kij, kji};
            for (int perm : perms) {
                X_tno_[perm] = X_tno_ijk;
                e_tno_[perm] = e_tno_ijk;
                n_tno_[perm] = X_tno_ijk->colspi(0);
            }
        }
    }

    timer_off("TNO transform");
}

void DLPNOMP2::compute_pno_tno_overlaps() {
    timer_on("PNO/TNO overlaps");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    S_pno_tno_ij_ilm_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(static, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        if (n_pno_[ij] == 0) continue;

        S_pno_tno_ij_ilm_[ij].resize(n_lmo_pairs);

        for (int lm = 0; lm < n_lmo_pairs; ++lm) {
            int l, m;
            std::tie(l, m) = ij_to_i_j_[lm];
            int ilm_dense = i * naocc * naocc + l * naocc + m;
            if (i_j_k_to_ijk_.count(ilm_dense)) {
                int ilm = i_j_k_to_ijk_[ilm_dense];

                if (n_tno_[ilm] == 0) continue;

                S_pno_tno_ij_ilm_[ij][lm] = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], \
                                                                    lmotriplet_to_paos_[ilm]);
                S_pno_tno_ij_ilm_[ij][lm] = linalg::triplet(X_pno_[ij], S_pno_tno_ij_ilm_[ij][lm], X_tno_[ilm], 
                                                            true, false, false);
            }
        }
    }

    timer_off("PNO/TNO overlaps");
}

void DLPNOMP2::compute_tno_overlaps() {

    timer_on("TNO overlaps");

    int n_lmo_triplets = ijk_to_i_j_k_.size();
    int naocc = nalpha_ - nfrzc();

    S_tno_ijk_ljk_.resize(n_lmo_triplets);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        int i, j, k;
        std::tie(i, j, k) = ijk_to_i_j_k_[ijk];

        int ntno_ijk = n_tno_[ijk];
        if (ntno_ijk == 0) continue;

        S_tno_ijk_ljk_[ijk].resize(naocc);

        for (int l = 0; l < naocc; l++) {
            int ljk_dense = l * naocc * naocc + j * naocc + k;
            if (!i_j_k_to_ijk_.count(ljk_dense)) continue;

            int ljk = i_j_k_to_ijk_[ljk_dense];
            if (n_tno_[ljk] == 0) continue;

            auto S_pao_ijk_ljk = submatrix_rows_and_cols(*S_pao_, lmotriplet_to_paos_[ijk], lmotriplet_to_paos_[ljk]);
            S_tno_ijk_ljk_[ijk][l] = linalg::triplet(X_tno_[ijk], S_pao_ijk_ljk, X_tno_[ljk], true, false, false);
        }
    }

    timer_off("TNO overlaps");
}

void DLPNOMP2::compute_W_iajbkc() {
    timer_on("Compute W_iajbkc");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_triplets = ijk_to_i_j_k_.size();

    std::vector<SharedMatrix> W_temp(n_lmo_triplets);

#pragma omp parallel for schedule(dynamic)
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        int i, j, k;
        std::tie(i, j, k) = ijk_to_i_j_k_[ijk];

        int ii = i_j_to_ij_[i][i];
        int ij = i_j_to_ij_[i][j], jk = i_j_to_ij_[j][k], ik = i_j_to_ij_[i][k];
        int kj = ij_to_ji_[jk];

        int ntno_ijk = n_tno_[ijk];

        if (ntno_ijk == 0) continue;
        
        W_temp[ijk] = std::make_shared<Matrix>(ntno_ijk, ntno_ijk * ntno_ijk);

        auto T_kj = linalg::doublet(T_iajb_[kj], S_PNO(kj, ii));
        auto K_temp1 = linalg::doublet(T_kj, K_maef_[ii]);
        K_temp1 = linalg::doublet(S_pno_tno_ij_ilm_[kj][ij], K_temp1, true, false);

        for (int c_ijk = 0; c_ijk < ntno_ijk; c_ijk++) {
            std::vector<int> c_ijk_slice(1, c_ijk);
            auto K_temp1a = submatrix_rows(*K_temp1, c_ijk_slice);
            K_temp1a->reshape(ntno_ijk, ntno_ijk);

            K_temp1a = linalg::triplet(S_pno_tno_ij_ilm_[ii][kj], K_temp1a, S_pno_tno_ij_ilm_[ii][kj], true, false, false);

            for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
                for (int b_ijk = 0; b_ijk < ntno_ijk; b_ijk++) {
                    (*W_temp[ijk])(a_ijk, b_ijk * ntno_ijk + c_ijk) = (*K_temp1a)(a_ijk, b_ijk);
                }
            }
        }

        for (int l_jk = 0; l_jk < lmopair_to_lmos_[jk].size(); l_jk++) {
            int l = lmopair_to_lmos_[jk][l_jk];
            int il = i_j_to_ij_[i][l];

            if (il == -1 || n_pno_[il] == 0) continue;

            std::vector<int> l_jk_slice(1, l_jk);
            auto K_temp2 = linalg::doublet(submatrix_rows(*K_mbij_[jk], l_jk_slice), S_pno_tno_ij_ilm_[jk][ik]);

            auto T_il = linalg::doublet(S_pno_tno_ij_ilm_[il][jk], T_iajb_[il], true, false);
            T_il = linalg::doublet(T_il, S_pno_tno_ij_ilm_[il][jk], false, false);

            for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
                for (int b_ijk = 0; b_ijk < ntno_ijk; b_ijk++) {
                    for (int c_ijk = 0; c_ijk < ntno_ijk; c_ijk++) {
                        (*W_temp[ijk])(a_ijk, b_ijk * ntno_ijk + c_ijk) -= (*T_il)(a_ijk, b_ijk) * (*K_temp2)(0, c_ijk);
                    }
                }
            }

            //C_DGER(ntno_ijk * ntno_ijk, ntno_ijk, -1.0, &(*T_il)(0,0), 1, &(*K_temp2)(0,0), 1, &(*W_temp[ijk])(0,0), ntno_ijk * ntno_ijk);
        }
    }

    W_iajbkc_.resize(n_lmo_triplets);
#pragma omp parallel for schedule(dynamic)
    for (int ijk = 0; ijk < n_lmo_triplets; ijk++) {
        int i, j, k;
        std::tie(i, j, k) = ijk_to_i_j_k_[ijk];

        int ntno_ijk = n_tno_[ijk];
        W_iajbkc_[ijk] = std::make_shared<Matrix>(ntno_ijk, ntno_ijk * ntno_ijk);

        int ikj = i_j_k_to_ijk_[i * naocc * naocc + k * naocc + j];
        int jik = i_j_k_to_ijk_[j * naocc * naocc + i * naocc + k];
        int jki = i_j_k_to_ijk_[j * naocc * naocc + k * naocc + i];
        int kij = i_j_k_to_ijk_[k * naocc * naocc + i * naocc + j];
        int kji = i_j_k_to_ijk_[k * naocc * naocc + j * naocc + i];

        for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
            for (int b_ijk = 0; b_ijk < ntno_ijk; b_ijk++) {
                for (int c_ijk = 0; c_ijk < ntno_ijk; c_ijk++) {
                    (*W_iajbkc_[ijk])(a_ijk, b_ijk * ntno_ijk + c_ijk) = (*W_temp[ijk])(a_ijk, b_ijk * ntno_ijk + c_ijk) +
                        (*W_temp[ikj])(a_ijk, c_ijk * ntno_ijk + b_ijk) + (*W_temp[jik])(b_ijk, a_ijk * ntno_ijk + c_ijk) +
                        (*W_temp[jki])(b_ijk, c_ijk * ntno_ijk + a_ijk) + (*W_temp[kij])(c_ijk, a_ijk * ntno_ijk + b_ijk) + 
                        (*W_temp[kji])(c_ijk, b_ijk * ntno_ijk + a_ijk);
                }
            }
        }

    }

    timer_off("Compute W_iajbkc");
}

double DLPNOMP2::compute_t_energy() {
    timer_on("Compute (T) Energy");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_triplets = ijk_to_i_j_k_.size();

    double E_T = 0.0;

#pragma omp parallel for schedule(dynamic) reduction(+ : E_T)
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        int i, j, k;
        std::tie(i, j, k) = ijk_to_i_j_k_[ijk];

        int ij = i_j_to_ij_[i][j], jk = i_j_to_ij_[j][k], ik = i_j_to_ij_[i][k];
        int ii = i_j_to_ij_[i][i], jj = i_j_to_ij_[j][j], kk = i_j_to_ij_[k][k];

        int ntno_ijk = n_tno_[ijk];

        auto V_ijk = W_iajbkc_[ijk]->clone();

        auto T_temp1 = linalg::doublet(T_ia_[i], S_pno_tno_ij_ilm_[ii][jk], true, false);
        auto K_temp1 = linalg::triplet(S_pno_tno_ij_ilm_[jk][ik], K_iajb_[jk], S_pno_tno_ij_ilm_[jk][ik], true, false, false);

        auto T_temp2 = linalg::doublet(T_ia_[j], S_pno_tno_ij_ilm_[jj][ik], true, false);
        auto K_temp2 = linalg::triplet(S_pno_tno_ij_ilm_[ik][jk], K_iajb_[ik], S_pno_tno_ij_ilm_[ik][jk], true, false, false);

        auto T_temp3 = linalg::doublet(T_ia_[k], S_pno_tno_ij_ilm_[kk][ij], true, false);
        auto K_temp3 = linalg::triplet(S_pno_tno_ij_ilm_[ij][jk], K_iajb_[ij], S_pno_tno_ij_ilm_[ij][jk], true, false, false);

        for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
            for (int b_ijk = 0; b_ijk < ntno_ijk; b_ijk++) {
                for (int c_ijk = 0; c_ijk < ntno_ijk; c_ijk++) {
                    (*V_ijk)(a_ijk, b_ijk * ntno_ijk + c_ijk) += (*T_temp1)(0, a_ijk) * (*K_temp1)(b_ijk, c_ijk) +
                            (*T_temp2)(0, b_ijk) * (*K_temp2)(a_ijk, c_ijk) + (*T_temp3)(0, c_ijk) * (*K_temp3)(a_ijk, b_ijk);
                }
            }
        }

        auto Vt_ijk = V_ijk->clone();
        Vt_ijk->scale(4.0);
        for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
            for (int b_ijk = 0; b_ijk < ntno_ijk; b_ijk++) {
                for (int c_ijk = 0; c_ijk < ntno_ijk; c_ijk++) {
                    (*Vt_ijk)(a_ijk, b_ijk * ntno_ijk + c_ijk) += -6.0 * (*V_ijk)(c_ijk, b_ijk * ntno_ijk + a_ijk) +
                            2.0 * (*V_ijk)(c_ijk, a_ijk * ntno_ijk + b_ijk);
                }
            }
        }

        E_T += Vt_ijk->vector_dot(T_iajbkc_[ijk]) / 3.0;
    }

    timer_off("Compute (T) Energy");

    return E_T;
}

void DLPNOMP2::lccsd_t_iterations() {
    timer_on("LCCSD(T) Iterations");

    outfile->Printf("\n  ==> Local CCSD(T) <==\n\n");
    outfile->Printf("    E_CONVERGENCE = %.2e\n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE = %.2e\n\n", options_.get_double("R_CONVERGENCE"));
    outfile->Printf("                         Corr. Energy    Delta E     Max R\n");

    // => Initialize Triples Amplitude <= //

    int naocc = nalpha_ - nfrzc();
    int n_lmo_triplets = ijk_to_i_j_k_.size();
    T_iajbkc_.resize(n_lmo_triplets);

#pragma omp parallel for schedule(dynamic)
    for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        int ntno_ijk = n_tno_[ijk];

        T_iajbkc_[ijk] = std::make_shared<Matrix>(ntno_ijk, ntno_ijk * ntno_ijk);
    }

    // => Initialize Triples Residuals <= //

    std::vector<SharedMatrix> R_iajbkc(n_lmo_triplets);

    int iteration = 0, max_iteration = options_.get_int("DLPNO_MAXITER");
    double e_curr = 0.0, e_prev = 0.0, r_curr = 0.0;
    bool e_converged = false, r_converged = false;

    DIISManager diis(options_.get_int("DIIS_MAX_VECS"), "LCCSD(T) DIIS", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::InCore);

    while (!(e_converged && r_converged)) {
        // RMS of residual per single LMO, for assesing convergence
        std::vector<double> R_iajbkc_rms(n_lmo_triplets, 0.0);

#pragma omp parallel for schedule(dynamic)
        for (int ijk = 0; ijk < n_lmo_triplets; ++ijk) {
            int i, j, k;
            std::tie(i, j, k) = ijk_to_i_j_k_[ijk];

            int ntno_ijk = n_tno_[ijk];

            if (ntno_ijk == 0) continue;

            R_iajbkc[ijk] = W_iajbkc_[ijk]->clone();
            for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
                for (int b_ijk = 0; b_ijk < ntno_ijk; b_ijk++) {
                    for (int c_ijk = 0; c_ijk < ntno_ijk; c_ijk++) {
                        (*R_iajbkc[ijk])(a_ijk, b_ijk * ntno_ijk + c_ijk) += (*T_iajbkc_[ijk])(a_ijk, b_ijk * ntno_ijk + c_ijk) * 
                                (e_tno_[ijk]->get(a_ijk) + e_tno_[ijk]->get(b_ijk) + e_tno_[ijk]->get(c_ijk) 
                                - (*F_lmo_)(i, i) - (*F_lmo_)(j, j) - (*F_lmo_)(k, k));
                    }
                }
            }

            int kij = i_j_k_to_ijk_[k * naocc * naocc + i * naocc + j];
            int jik = i_j_k_to_ijk_[j * naocc * naocc + i * naocc + k];

            for (int l = 0; l < naocc; l++) {
                int ijl_dense = i * naocc * naocc + j * naocc + l;
                if (l != k && i_j_k_to_ijk_.count(ijl_dense)) {
                    int ijl = i_j_k_to_ijk_[ijl_dense];
                    if (n_tno_[ijl] == 0) continue;

                    auto T_temp1 = linalg::doublet(S_tno_ijk_ljk_[kij][l], T_iajbkc_[ijl]);
                    for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
                        std::vector<int> a_ijk_slice(1, a_ijk);
                        auto T_temp2 = submatrix_rows(*T_temp1, a_ijk_slice);
                        T_temp2->reshape(ntno_ijk, ntno_ijk);
                        T_temp2 = linalg::triplet(S_tno_ijk_ljk_[kij][l], T_temp2, S_tno_ijk_ljk_[kij][l], false, false, true);

                        C_DAXPY(ntno_ijk * ntno_ijk, -(*F_lmo_)(l,k), &(*T_temp2)(0,0), 1, &(*R_iajbkc[ijk])(a_ijk, 0), 1);
                    }
                }

                int ilk_dense = i * naocc * naocc + l * naocc + k;
                if (l != j && i_j_k_to_ijk_.count(ilk_dense)) {
                    int ilk = i_j_k_to_ijk_[ilk_dense];
                    if (n_tno_[ilk] == 0) continue;

                    auto T_temp1 = linalg::doublet(S_tno_ijk_ljk_[jik][l], T_iajbkc_[ilk]);
                    for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
                        std::vector<int> a_ijk_slice(1, a_ijk);
                        auto T_temp2 = submatrix_rows(*T_temp1, a_ijk_slice);
                        T_temp2->reshape(ntno_ijk, ntno_ijk);
                        T_temp2 = linalg::triplet(S_tno_ijk_ljk_[jik][l], T_temp2, S_tno_ijk_ljk_[jik][l], false, false, true);

                        C_DAXPY(ntno_ijk * ntno_ijk, -(*F_lmo_)(l,j), &(*T_temp2)(0,0), 1, &(*R_iajbkc[ijk])(a_ijk, 0), 1);
                    }
                }

                int ljk_dense = l * naocc * naocc + j * naocc + k;
                if (l != i && i_j_k_to_ijk_.count(ljk_dense)) {
                    int ljk = i_j_k_to_ijk_[ljk_dense];
                    if (n_tno_[ljk] == 0) continue;

                    auto T_temp1 = linalg::doublet(S_tno_ijk_ljk_[ijk][l], T_iajbkc_[ljk]);
                    for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
                        std::vector<int> a_ijk_slice(1, a_ijk);
                        auto T_temp2 = submatrix_rows(*T_temp1, a_ijk_slice);
                        T_temp2->reshape(ntno_ijk, ntno_ijk);
                        T_temp2 = linalg::triplet(S_tno_ijk_ljk_[ijk][l], T_temp2, S_tno_ijk_ljk_[ijk][l], false, false, true);

                        C_DAXPY(ntno_ijk * ntno_ijk, -(*F_lmo_)(l,i), &(*T_temp2)(0,0), 1, &(*R_iajbkc[ijk])(a_ijk, 0), 1);
                    }
                }
            }

            // => Update T2 Amplitudes <= //
            for (int a_ijk = 0; a_ijk < ntno_ijk; a_ijk++) {
                for (int b_ijk = 0; b_ijk < ntno_ijk; b_ijk++) {
                    for (int c_ijk = 0; c_ijk < ntno_ijk; c_ijk++) {
                        (*T_iajbkc_[ijk])(a_ijk, b_ijk * ntno_ijk + c_ijk) += -(*R_iajbkc[ijk])(a_ijk, b_ijk * ntno_ijk + c_ijk) /
                                (e_tno_[ijk]->get(a_ijk) + e_tno_[ijk]->get(b_ijk) + e_tno_[ijk]->get(c_ijk) 
                                - (*F_lmo_)(i, i) - (*F_lmo_)(j, j) - (*F_lmo_)(k, k));
                    }
                }
            }

            R_iajbkc_rms[ijk] = R_iajbkc[ijk]->rms();
        }

        // => DIIS Extrapolation <= //
        auto T_iajbkc_flat = flatten_mats(T_iajbkc_);
        auto R_iajbkc_flat = flatten_mats(R_iajbkc);

        if (iteration == 0) {
            diis.set_error_vector_size(R_iajbkc_flat.get());
            diis.set_vector_size(T_iajbkc_flat.get());
        }

        diis.add_entry(R_iajbkc_flat.get(), T_iajbkc_flat.get());
        diis.extrapolate(T_iajbkc_flat.get());

        copy_flat_mats(T_iajbkc_flat, T_iajbkc_);

        // evaluate convergence
        e_prev = e_curr;
        // Compute LCCSD(T) energy
        e_curr = e_lccsd_ + compute_t_energy();

        double r_curr = *max_element(R_iajbkc_rms.begin(), R_iajbkc_rms.end());

        r_converged = fabs(r_curr) < options_.get_double("R_CONVERGENCE");
        e_converged = fabs(e_curr - e_prev) < options_.get_double("E_CONVERGENCE");

        outfile->Printf("  @LCCSD(T) iter %3d: %16.12f %10.3e %10.3e\n", iteration, e_curr, e_curr - e_prev, r_curr);

        iteration++;

        if (iteration > max_iteration) {
            throw PSIEXCEPTION("Maximum DLPNO iterations exceeded.");
        }
    }

    e_lccsd_t_ = e_curr;

    timer_off("LCCSD(T) Iterations");
}

SharedMatrix DLPNOMP2::get_lmo_matrix(std::string key) {
    if (!lmo_matrices_.size()) store_information();
    return lmo_matrices_[key];
}

SharedMatrix DLPNOMP2::get_pao_matrix(std::string key) {
    if (!pao_matrices_.size()) store_information();
    return pao_matrices_[key];
}

std::vector<SharedMatrix> DLPNOMP2::get_pno_matrix(std::string key) {
    if (!pno_matrices_.size()) store_information();
    return pno_matrices_[key];
}

SparseMap DLPNOMP2::get_sparse_map(std::string key) {
    if (!sparse_maps_.size()) store_information();
    return sparse_maps_[key];
}

std::vector<SharedMatrix> DLPNOMP2::get_qia() {
    return qia_;
}

std::vector<SharedMatrix> DLPNOMP2::get_qij() {
    if (qij_.empty()) compute_qij();
    return qij_;
}

std::vector<SharedMatrix> DLPNOMP2::get_qab() {
    if (qab_.empty()) compute_qab();
    return qab_;
}

SharedMatrix DLPNOMP2::compute_Fmi(const std::vector<SharedMatrix>& tau_tilde) {
    timer_on("Compute Fmi");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    // Equation 40, Term 1
    auto Fmi = F_lmo_->clone();

#pragma omp parallel for schedule(dynamic, 1)
    for (int mi = 0; mi < n_lmo_pairs; ++mi) {
        int m, i;
        std::tie(m, i) = ij_to_i_j_[mi];

        if (n_pno_[mi] == 0) continue;

        for (int n_mi = 0; n_mi < lmopair_to_lmos_[mi].size(); n_mi++) {
            int n = lmopair_to_lmos_[mi][n_mi];
            int mn = i_j_to_ij_[m][n], nn = i_j_to_ij_[n][n], in = i_j_to_ij_[i][n];
            int nm = ij_to_ji_[mn];

            if (n_pno_[mn] == 0 || n_pno_[in] == 0) continue;

            int i_mn = lmopair_to_lmos_dense_[mn][i];

            // Equation 40, Term 3
            std::vector<int> i_mn_slice(1, i_mn);
            auto l_mn_temp = submatrix_rows(*K_mbij_[mn], i_mn_slice);
            l_mn_temp->scale(2.0);
            l_mn_temp->subtract(submatrix_rows(*K_mbij_[nm], i_mn_slice));
            SharedMatrix S_nn_mn = S_pno_ij_kj_[nn][m];
            l_mn_temp = linalg::doublet(S_nn_mn, l_mn_temp, false, true);
            (*Fmi)(m,i) += l_mn_temp->vector_dot(T_ia_[n]);

            // Equation 40, Term 4
            auto S_in_mn = S_pno_ij_kj_[in][m];
            l_mn_temp = linalg::triplet(S_in_mn, L_iajb_[mn], S_in_mn, false, false, true);
            (*Fmi)(m,i) += l_mn_temp->vector_dot(tau_tilde[in]);
        }
    }

    timer_off("Compute Fmi");

    return Fmi;
}

std::vector<SharedMatrix> DLPNOMP2::compute_Fbe(const std::vector<SharedMatrix>& tau_tilde) {
    timer_on("Compute Fbe");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    std::vector<SharedMatrix> Fbe(n_lmo_pairs);

#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        int npno_ij = n_pno_[ij];
        if (npno_ij == 0) continue;

        // Equation 39, Term 1
        Fbe[ij] = std::make_shared<Matrix>("Fbe", npno_ij, npno_ij);
        Fbe[ij]->zero();
        Fbe[ij]->set_diagonal(e_pno_[ij]);

        // See to it, brothers and sisters, that none of you has a sinful, 
        // unbelieving heart that turns away from the living God (Hebrews 3:12)

        for (int m_ij = 0; m_ij < lmopair_to_lmos_[ij].size(); m_ij++) {
            int m = lmopair_to_lmos_[ij][m_ij];
            int im = i_j_to_ij_[i][m], mm = i_j_to_ij_[m][m], mj = i_j_to_ij_[m][j];
            int npno_mm = n_pno_[mm], npno_mj = n_pno_[mj];

            if (n_pno_[mj] != 0) {
                auto S_ij_mj = S_pno_ij_kj_[ij][m];
                auto S_mm_mj = S_pno_ij_ik_[mm][j];
                auto T_m_temp = linalg::doublet(S_mm_mj, T_ia_[m], true, false);

                SharedMatrix Fbe_mj = std::make_shared<Matrix>("Fbe_mj", npno_mj, npno_mj);
                Fbe_mj->zero();

                for (int a_mj = 0; a_mj < npno_mj; ++a_mj) {
                    std::vector<int> a_mj_slice(1, a_mj);
                    SharedMatrix K_maef_slice = submatrix_rows(*K_maef_[mj], a_mj_slice);
                    K_maef_slice->reshape(npno_mj, npno_mj);

                    SharedMatrix Fbe_temp1 = linalg::doublet(K_maef_slice, T_m_temp, true, false);
                    C_DAXPY(npno_mj, 2.0, &(*Fbe_temp1)(0,0), 1, &(*Fbe_mj)(a_mj, 0), 1);

                    SharedMatrix Fbe_temp2 = linalg::doublet(K_maef_slice, T_m_temp, false, false);
                    C_DAXPY(npno_mj, -1.0, &(*Fbe_temp2)(0,0), 1, &(*Fbe_mj)(a_mj, 0), 1);
                }

                Fbe[ij]->add(linalg::triplet(S_ij_mj, Fbe_mj, S_ij_mj, false, false, true));
            }

            for (int n_ij = 0; n_ij < lmopair_to_lmos_[ij].size(); n_ij++) {
                int n = lmopair_to_lmos_[ij][n_ij];
                int mn = i_j_to_ij_[m][n];

                if (mn != -1 && n_pno_[mn] != 0) {
                    SharedMatrix S_ij_mn = S_PNO(ij, mn);
                    auto tau_L_temp = linalg::triplet(tau_tilde[mn], L_iajb_[mn], S_ij_mn, false, true, true);
                    Fbe[ij]->subtract(linalg::doublet(S_ij_mn, tau_L_temp, false, false));
                }
            }
        }

        // But encourage one another daily, as long as it is called "Today",
        // so that none of you may be hardened by sin's deceitfulness (Hebrews 3:13)

    }

    timer_off("Compute Fbe");

    return Fbe;
}

std::vector<SharedMatrix> DLPNOMP2::compute_Fme() {
    timer_on("Compute Fme");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    std::vector<SharedMatrix> Fme(n_lmo_pairs);

#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int npno_ij = n_pno_[ij];

        if (npno_ij == 0) continue;

        Fme[ij] = std::make_shared<Matrix>("Fme", nlmo_ij, npno_ij);
        Fme[ij]->zero();

        for (int m_ij = 0; m_ij < lmopair_to_lmos_[ij].size(); m_ij++) {
            int m = lmopair_to_lmos_[ij][m_ij];
            int mm = i_j_to_ij_[m][m], mj = i_j_to_ij_[m][j];

            for (int n_ij = 0; n_ij < lmopair_to_lmos_[ij].size(); n_ij++) {
                int n = lmopair_to_lmos_[ij][n_ij];
                int mn = i_j_to_ij_[m][n], nn = i_j_to_ij_[n][n];

                if (mn != -1 && n_pno_[mn] != 0) {
                    SharedMatrix S_mn_nn = S_pno_ij_kj_[mn][n];
                    SharedMatrix T_n_temp = linalg::doublet(S_mn_nn, T_ia_[n], false, false);

                    SharedMatrix S_mn_ij = S_PNO(mn, ij);
                    auto F_me_temp = linalg::triplet(S_mn_ij, L_iajb_[mn], T_n_temp, true, false, false);
                    C_DAXPY(npno_ij, 1.0, &(*F_me_temp)(0,0), 1, &(*Fme[ij])(m_ij, 0), 1);
                }
            }
        }
    }

    timer_off("Compute Fme");

    return Fme;
}

std::vector<SharedMatrix> DLPNOMP2::compute_Wmnij(const std::vector<SharedMatrix>& tau) {
    timer_on("Compute Wmnij");
    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    std::vector<SharedMatrix> Wmnij(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int mn = 0; mn < n_lmo_pairs; ++mn) {
        int m, n;
        std::tie(m, n) = ij_to_i_j_[mn];
        int nm = ij_to_ji_[mn];

        int npno_mn = n_pno_[mn];
        if (npno_mn == 0) continue;

        Wmnij[mn] = K_mnij_[mn]->clone();

        int nlmo_mn = lmopair_to_lmos_[mn].size();
        SharedMatrix T_i_mn = std::make_shared<Matrix>(nlmo_mn, npno_mn);
        T_i_mn->zero();

        for (int i_mn = 0; i_mn < lmopair_to_lmos_[mn].size(); i_mn++) {
            int i = lmopair_to_lmos_[mn][i_mn];
            int ii = i_j_to_ij_[i][i], im = i_j_to_ij_[i][m], in = i_j_to_ij_[i][n];

            auto S_ii_mn = S_PNO(ii, mn);
            auto T_temp = linalg::doublet(S_ii_mn, T_ia_[i], true, false);
            C_DCOPY(npno_mn, &(*T_temp)(0,0), 1, &(*T_i_mn)(i_mn, 0), 1);
        }

        Wmnij[mn]->add(linalg::doublet(K_mbij_[mn], T_i_mn, false, true));
        Wmnij[mn]->add(linalg::doublet(T_i_mn, K_mbij_[nm], false, true));

        for (int ij_mn = 0; ij_mn < nlmo_mn * nlmo_mn; ++ij_mn) {
            int i_mn = ij_mn / nlmo_mn, j_mn = ij_mn % nlmo_mn;
            int i = lmopair_to_lmos_[mn][i_mn], j = lmopair_to_lmos_[mn][j_mn];
            int ij = i_j_to_ij_[i][j];
            
            if (ij != -1 && n_pno_[ij] != 0) {
                auto S_mn_ij = S_PNO(mn, ij);
                auto K_temp = linalg::triplet(S_mn_ij, K_iajb_[mn], S_mn_ij, true, false, false);
                (*Wmnij[mn])(i_mn, j_mn) += K_temp->vector_dot(tau[ij]);
            }
        }
    }

    timer_off("Compute Wmnij");

    return Wmnij;
}

std::vector<SharedMatrix> DLPNOMP2::compute_Wmbej(const std::vector<SharedMatrix>& tau_bar) {
    timer_on("Compute Wmbej");
    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    std::vector<SharedMatrix> Wmbej(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int mj = 0; mj < n_lmo_pairs; ++mj) {
        int m, j;
        std::tie(m, j) = ij_to_i_j_[mj];
        int mm = i_j_to_ij_[m][m], jj = i_j_to_ij_[j][j], jm = ij_to_ji_[mj];
        int npno_mm = n_pno_[mm], npno_mj = n_pno_[mj];

        if (npno_mj == 0) continue;

        Wmbej[mj] = K_iajb_[jm]->clone();
        
        SharedMatrix S_mm_mj = S_pno_ij_ik_[mm][j];
        SharedMatrix S_mj_jj = S_pno_ij_kj_[mj][j];
        SharedMatrix tia_temp = linalg::doublet(S_mj_jj, T_ia_[j], false, false);

        SharedMatrix K_temp1 = std::make_shared<Matrix>(npno_mj, npno_mj);
        K_temp1->zero();

        /*
        for (int q_mm = 0; q_mm < lmopair_to_ribfs_[mm].size(); q_mm++) {
            SharedMatrix q_e = Qma_mm_[m][q_mm];
            SharedMatrix q_bf = Qab_mm_[m][q_mm];

            SharedMatrix w_temp = linalg::doublet(q_bf, tia_temp);
            C_DGER(npno_mm, npno_mm, 1.0, &(*w_temp)(0,0), 1, &(*q_e)(0,0), 1, &(*K_temp1)(0,0), npno_mm);
        }
        */
        for (int b_mj = 0; b_mj < npno_mj; b_mj++) {
            std::vector<int> b_mj_slice(1, b_mj);
            SharedMatrix K_vv = submatrix_rows(*K_maef_[mj], b_mj_slice);
            K_vv->reshape(npno_mj, npno_mj);
            SharedMatrix K_temp2 = linalg::doublet(K_vv, tia_temp, false, false);
            C_DAXPY(npno_mj, 1.0, &(*K_temp2)(0,0), 1, &(*K_temp1)(b_mj,0), 1);
        }
        Wmbej[mj]->add(K_temp1);

        for (int n_mj = 0; n_mj < lmopair_to_lmos_[mj].size(); n_mj++) {
            int n = lmopair_to_lmos_[mj][n_mj];
            int mn = i_j_to_ij_[m][n], nn = i_j_to_ij_[n][n], jn = i_j_to_ij_[j][n];
            int nm = ij_to_ji_[mn], nj = ij_to_ji_[jn];

            SharedMatrix S_nn_mj = S_PNO(nn, mj);
            SharedMatrix t_n_temp = linalg::doublet(S_nn_mj, T_ia_[n], true, false);
            C_DGER(npno_mj, npno_mj, -1.0, &(*t_n_temp)(0, 0), 1, &(*K_mbij_[jm])(n_mj, 0), 1, &(*Wmbej[mj])(0, 0), npno_mj);

            if (n_pno_[mn] == 0 || n_pno_[nj] == 0) continue;

            SharedMatrix tau_temp = linalg::triplet(S_pno_ij_kj_[mn][j], tau_bar[jn], \
                                        S_pno_ij_ik_[jn][m], false, false, false);
            SharedMatrix K_mn_temp = linalg::doublet(S_pno_ij_ik_[mj][n], K_iajb_[mn], false, false);
            Wmbej[mj]->subtract(linalg::doublet(tau_temp, K_mn_temp, true, true));

            SharedMatrix T_nj_temp = linalg::triplet(S_pno_ij_kj_[mn][j], T_iajb_[nj], \
                                        S_pno_ij_ik_[jn][m], false, false, false);
            SharedMatrix L_mn_temp = linalg::doublet(S_pno_ij_ik_[mj][n], L_iajb_[mn], false, false);
            SharedMatrix TL_temp = linalg::doublet(T_nj_temp, L_mn_temp, true, true);
            TL_temp->scale(0.5);
            Wmbej[mj]->add(TL_temp);
        }
    }

    timer_off("Compute Wmbej");

    return Wmbej;
}

std::vector<SharedMatrix> DLPNOMP2::compute_Wmbje(const std::vector<SharedMatrix>& tau_bar) {

    timer_on("Compute Wmbje");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    std::vector<SharedMatrix> Wmbje(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int mj = 0; mj < n_lmo_pairs; ++mj) {
        int m, j;
        std::tie(m, j) = ij_to_i_j_[mj];
        int mm = i_j_to_ij_[m][m], jj = i_j_to_ij_[j][j], jm = ij_to_ji_[mj];
        int npno_mm = n_pno_[mm], npno_mj = n_pno_[mj];

        if (npno_mj == 0) continue;

        Wmbje[mj] = J_ijab_[mj]->clone();
        Wmbje[mj]->scale(-1.0);
        
        SharedMatrix S_mm_mj = S_pno_ij_ik_[mm][j];
        SharedMatrix S_mj_jj = S_pno_ij_kj_[mj][j];
        SharedMatrix tia_temp = linalg::doublet(S_mj_jj, T_ia_[j], false, false);

        SharedMatrix K_temp1 = std::make_shared<Matrix>(npno_mj, npno_mj);
        K_temp1->zero();

        /*
        for (int q_mm = 0; q_mm < lmopair_to_ribfs_[mm].size(); q_mm++) {
            SharedMatrix q_f = Qma_mm_[m][q_mm];
            SharedMatrix q_be = Qab_mm_[m][q_mm];

            SharedMatrix w_temp = q_be->clone();
            w_temp->scale(q_f->vector_dot(tia_temp));
            K_temp1->add(w_temp);
        }
        */
        for (int b_mj = 0; b_mj < npno_mj; b_mj++) {
            std::vector<int> b_mj_slice(1, b_mj);
            SharedMatrix K_vv = submatrix_rows(*K_maef_[mj], b_mj_slice);
            K_vv->reshape(npno_mj, npno_mj);
            SharedMatrix K_temp2 = linalg::doublet(K_vv, tia_temp, true, false);
            C_DAXPY(npno_mj, 1.0, &(*K_temp2)(0,0), 1, &(*K_temp1)(b_mj,0), 1);
        }

        Wmbje[mj]->subtract(K_temp1);

        for (int n_mj = 0; n_mj < lmopair_to_lmos_[mj].size(); n_mj++) {
            int n = lmopair_to_lmos_[mj][n_mj];
            int mn = i_j_to_ij_[m][n], nn = i_j_to_ij_[n][n], jn = i_j_to_ij_[j][n];
            int nm = ij_to_ji_[mn], nj = ij_to_ji_[jn];

            SharedMatrix S_nn_mj = S_PNO(nn, mj);
            SharedMatrix t_n_temp = linalg::doublet(S_nn_mj, T_ia_[n], true, false);

            if (n_pno_[nj] == 0 || n_pno_[mn] == 0) continue;

            int j_mn = lmopair_to_lmos_dense_[mn][j];
            std::vector<int> j_mn_slice(1, j_mn);
            SharedMatrix K_mn_temp = linalg::doublet(submatrix_rows(*K_mbij_[mn], j_mn_slice), S_pno_ij_ik_[mn][j])->transpose();
            C_DGER(npno_mj, npno_mj, 1.0, &(*t_n_temp)(0,0), 1, &(*K_mn_temp)(0, 0), 1, &(*Wmbje[mj])(0, 0), npno_mj);

            SharedMatrix tau_temp = linalg::triplet(S_pno_ij_kj_[mn][j], tau_bar[jn], \
                                        S_pno_ij_ik_[jn][m], false, false, false);
            K_mn_temp = linalg::doublet(K_iajb_[mn], S_pno_ij_ik_[mn][j], false, false);
            Wmbje[mj]->add(linalg::doublet(tau_temp, K_mn_temp, true, false));
        }
    }

    timer_off("Compute Wmbje");

    return Wmbje;
}

void DLPNOMP2::lccsd_iterations() {

    timer_on("LCCSD iterations");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    outfile->Printf("\n  ==> Local CCSD <==\n\n");
    outfile->Printf("    E_CONVERGENCE = %.2e\n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE = %.2e\n\n", options_.get_double("R_CONVERGENCE"));
    outfile->Printf("                     Corr. Energy    Delta E     Max R1     Max R2\n");

    // => Initialize Singles Amplitudes <= //

    T_ia_.resize(naocc);
#pragma omp parallel for
    for (int i = 0; i < naocc; ++i) {
        int ii = i_j_to_ij_[i][i];
        T_ia_[i] = std::make_shared<Matrix>(n_pno_[ii], 1);
        T_ia_[i]->zero();
    }

    // => Initialize Dressed Doubles Amplitudes <= //

    std::vector<SharedMatrix> tau(n_lmo_pairs);
    std::vector<SharedMatrix> tau_tilde(n_lmo_pairs);
    std::vector<SharedMatrix> tau_bar(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        if (n_pno_[ij] == 0) continue;

        tau[ij] = T_iajb_[ij]->clone();
        tau_tilde[ij] = T_iajb_[ij]->clone();
        tau_bar[ij] = T_iajb_[ij]->clone();
        tau_bar[ij]->scale(0.5);
    }

    // => Initialize Residuals <= //

    std::vector<SharedMatrix> R_ia(naocc);
    std::vector<SharedMatrix> Rn_iajb(n_lmo_pairs);
    std::vector<SharedMatrix> R_iajb(n_lmo_pairs);

    int iteration = 0, max_iteration = options_.get_int("DLPNO_MAXITER");
    double e_curr = 0.0, e_prev = 0.0, r1_curr = 0.0, r2_curr = 0.0;
    bool e_converged = false, r_converged = false;

    DIISManager diis1(options_.get_int("DIIS_MAX_VECS"), "LCCSD DIIS (T1)", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::InCore);
    DIISManager diis2(options_.get_int("DIIS_MAX_VECS"), "LCCSD DIIS (T2)", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::InCore);

    while (!(e_converged && r_converged)) {
        // RMS of residual per single LMO, for assesing convergence
        std::vector<double> R_ia_rms(naocc, 0.0);
        // RMS of residual per LMO pair, for assessing convergence
        std::vector<double> R_iajb_rms(n_lmo_pairs, 0.0);

        // Build one-particle intermediates
        auto Fmi = compute_Fmi(tau_tilde);
        auto Fbe = compute_Fbe(tau_tilde);
        auto Fme = compute_Fme();

        // Build two-particle W intermediates
        auto Wmnij = compute_Wmnij(tau);
        auto Wmbej = compute_Wmbej(tau_bar);
        auto Wmbje = compute_Wmbje(tau_bar);

        // Calculate singles residuals from current amplitudes
#pragma omp parallel for
        for (int i = 0; i < naocc; i++) {
            int ii = i_j_to_ij_[i][i];
            int npno_ii = n_pno_[ii];

            // Madriaga Eq. 34, Term 2
            R_ia[i] = linalg::doublet(Fbe[ii], T_ia_[i], false, false);
            /*
            for (int a_ii = 0; a_ii < npno_ii; a_ii++) {
                (*R_ia[i])(a_ii, 0) += e_pno_[ii]->get(a_ii) * (*T_ia_[i])(a_ii, 0);
            }
            */

            for (int m = 0; m < naocc; m++) {
                int im = i_j_to_ij_[i][m], mm = i_j_to_ij_[m][m], mi = i_j_to_ij_[m][i];
                int npno_mm = n_pno_[mm], npno_mi = n_pno_[mi];

                if (im == -1 || n_pno_[im] == 0) continue;

                // Madriaga Eq. 34, Term 5
                auto S_mm_im = S_pno_ij_kj_[mm][i];
                auto temp_t1 = linalg::doublet(S_mm_im, T_ia_[m], true, false);
                auto S_im_ii = S_pno_ij_ik_[im][i];
                R_ia[i]->add(linalg::triplet(S_im_ii, Lt_iajb_[im], temp_t1, true, false, false));

                // Madriaga Eq. 34, Term 3
                auto S_mm_ii = S_PNO(mm, ii);
                auto T_m_temp = linalg::doublet(S_mm_ii, T_ia_[m], true, false);
                T_m_temp->scale(Fmi->get(m,i));
                R_ia[i]->subtract(T_m_temp);

                // Madriaga Eq. 34, Term 4
                int m_im = lmopair_to_lmos_dense_[im][m];
                std::vector<int> m_im_slice(1, m_im);
                auto Fe_im = submatrix_rows(*Fme[im], m_im_slice);
                R_ia[i]->add(linalg::triplet(S_im_ii, Tt_iajb_[im], Fe_im, true, false, true));

                // Madriaga Eq. 34, Term 6
                auto T_mi_mm = Tt_iajb_[mi]->clone();
                T_mi_mm->reshape(npno_mi * npno_mi, 1);
                // (npno_ii, npno_mm) (npno_mm, npno_mm * npno_mm) (npno_mm * npno_mm, 1)
                R_ia[i]->add(linalg::triplet(S_im_ii, K_maef_[mi], T_mi_mm, true, false, false));
                /*
                for (int q_mm = 0; q_mm < lmopair_to_ribfs_[mm].size(); q_mm++) {
                    SharedMatrix q_e = Qma_mm_[m][q_mm];
                    SharedMatrix q_af = Qab_mm_[m][q_mm];

                    SharedMatrix r1_temp = linalg::triplet(q_af, T_mi_mm, q_e, false, true, false);
                    R_ia[i]->add(linalg::doublet(S_mm_ii, r1_temp, true, false));
                }
                */

                // Madriaga Eq. 34, Term 7
                for (int n_im = 0; n_im < lmopair_to_lmos_[im].size(); n_im++) {
                    int n = lmopair_to_lmos_[im][n_im];
                    int mn = i_j_to_ij_[m][n], nm = i_j_to_ij_[n][m];
                    if (n_pno_[mn] == 0) continue;

                    int i_mn = lmopair_to_lmos_dense_[mn][i];
                    std::vector<int> i_mn_slice(1, i_mn);
                    auto K_temp = submatrix_rows(*K_mbij_[mn], i_mn_slice);
                    K_temp->scale(2.0);
                    K_temp->subtract(submatrix_rows(*K_mbij_[nm], i_mn_slice));

                    auto S_ii_mn = S_PNO(ii, mn);
                    R_ia[i]->subtract(linalg::triplet(S_ii_mn, T_iajb_[mn], K_temp, false, false, true));
                }
            }
            R_ia_rms[i] = R_ia[i]->rms();
        }

        // Calculate residuals from current amplitudes
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];
            int ii = i_j_to_ij_[i][i], jj = i_j_to_ij_[j][j], ji = ij_to_ji_[ij];
            int npno_ij = n_pno_[ij], npno_ii = n_pno_[ii], npno_jj = n_pno_[jj];

            if (npno_ij == 0) continue;

            // Buffers for R2 (Save memory)
            SharedMatrix r2_temp;

            // Madriaga Eq. 35, Term 1
            Rn_iajb[ij] = K_iajb_[ij]->clone();
            Rn_iajb[ij]->scale(0.5);

            // Madriaga Eq. 35, Term 2a
            Rn_iajb[ij]->add(linalg::doublet(T_iajb_[ij], Fbe[ij], false, true));
            /*
            for (int a_ij = 0; a_ij < npno_ij; a_ij++) {
                for (int b_ij = 0; b_ij < npno_ij; b_ij++) {
                    (*Rn_iajb[ij])(a_ij, b_ij) += e_pno_[ij]->get(b_ij) * (*T_iajb_[ij])(a_ij, b_ij);
                }
            }
            */

            // Madriaga Eq. 35, Term 5
            if (virtual_storage_ == STORE) { // High Memory Algorithm
                r2_temp = tau[ij]->clone();
                r2_temp->reshape(npno_ij * npno_ij, 1);
                if (i > j) r2_temp = linalg::doublet(K_abef_[ji], r2_temp);
                else r2_temp = linalg::doublet(K_abef_[ij], r2_temp);
                r2_temp->reshape(npno_ij, npno_ij);
                r2_temp->scale(0.5);
                Rn_iajb[ij]->add(r2_temp);
            } else if (virtual_storage_ == DF_STORE) { // Intermediate Memory Algorithm
                for (int q_ij = 0; q_ij < lmopair_to_ribfs_[ij].size(); q_ij++) {
                    if (i > j) r2_temp = linalg::triplet(Qab_ij_[ji][q_ij], tau[ij], Qab_ij_[ji][q_ij]);
                    else r2_temp = linalg::triplet(Qab_ij_[ij][q_ij], tau[ij], Qab_ij_[ij][q_ij]);
                    r2_temp->scale(0.5);
                    Rn_iajb[ij]->add(r2_temp);
                }
            } else { // Low Memory Algorithm
                int naux_ij = lmopair_to_ribfs_[ij].size();
                auto q_ab = std::make_shared<Matrix>(naux_ij, npno_ij * npno_ij);
                for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
                    int q = lmopair_to_ribfs_[ij][q_ij];
                    int centerq = ribasis_->function_to_center(q);

                    auto ab_temp = submatrix_rows_and_cols(*qab_[q], lmopair_pao_to_riatom_pao_[ij][q_ij], 
                                                        lmopair_pao_to_riatom_pao_[ij][q_ij]);
                    ab_temp = linalg::triplet(X_pno_[ij], ab_temp, X_pno_[ij], true, false, false);
                    C_DCOPY(npno_ij * npno_ij, &(*ab_temp)(0,0), 1, &(*q_ab)(q_ij, 0), 1);
                }
                auto A_solve = submatrix_rows_and_cols(*full_metric_, lmopair_to_ribfs_[ij], lmopair_to_ribfs_[ij]);
                A_solve->power(0.5, 1.0e-14);
                C_DGESV_wrapper(A_solve, q_ab);

                for (int q_ij = 0; q_ij < lmopair_to_ribfs_[ij].size(); q_ij++) {
                    std::vector<int> q_ij_slice(1, q_ij);
                    r2_temp = submatrix_rows(*q_ab, q_ij_slice);
                    r2_temp->reshape(npno_ij, npno_ij);
                    r2_temp = linalg::triplet(r2_temp, tau[ij], r2_temp);
                    r2_temp->scale(0.5);
                    Rn_iajb[ij]->add(r2_temp);
                }
                q_ab = nullptr;
            }

            // Madriaga Eq. 35, Term 12
            auto S_ij_ii = S_pno_ij_ik_[ij][i];
            r2_temp = linalg::doublet(S_ij_ii, T_ia_[i], false, false);
            /*
            SharedMatrix r_jj_temp = std::make_shared<Matrix>(npno_jj, npno_jj);
            r_jj_temp->zero();
            for (int q_jj = 0; q_jj < lmopair_to_ribfs_[jj].size(); q_jj++) {
                SharedMatrix q_b = Qma_mm_[j][q_jj];
                SharedMatrix q_ae = Qab_mm_[j][q_jj];

                SharedMatrix q_a = linalg::doublet(q_ae, r2_temp);
                C_DGER(npno_jj, npno_jj, 1.0, &(*q_a)(0,0), 1, &(*q_b)(0,0), 1, &(*r_jj_temp)(0,0), npno_jj);
            }
            */
            r2_temp = linalg::doublet(r2_temp, K_maef_[ji], true, false);
            r2_temp->reshape(npno_ij, npno_ij);
            auto eye = std::make_shared<Matrix>(npno_ij, npno_ij);
            eye->identity();
            Rn_iajb[ij]->add(linalg::doublet(r2_temp, eye, true, false));

            for (int m_ij = 0; m_ij < lmopair_to_lmos_[ij].size(); m_ij++) {
                int m = lmopair_to_lmos_[ij][m_ij];
                int im = i_j_to_ij_[i][m], mm = i_j_to_ij_[m][m], mj = i_j_to_ij_[m][j];
                int mi = ij_to_ji_[im];
                int npno_mj = n_pno_[mj];

                // Shared Intermediates
                auto S_ij_im = S_pno_ij_ik_[ij][m];
                auto S_ij_mj = S_pno_ij_kj_[ij][m];
                auto S_ij_mm = S_PNO(ij, mm);
                auto temp_t1 = linalg::doublet(S_ij_mm, T_ia_[m], false, false);

                // Madriaga Eq. 35, Term 2b
                std::vector<int> m_ij_slice(1, m_ij);
                r2_temp = submatrix_rows(*Fme[ij], m_ij_slice);
                r2_temp = linalg::doublet(T_iajb_[ij], r2_temp, false, true);
                C_DGER(npno_ij, npno_ij, -0.5, &(*r2_temp)(0,0), 1, &(*temp_t1)(0,0), 1, &(*Rn_iajb[ij])(0,0), npno_ij);
                
                if (n_pno_[mj] != 0) {
                    // Madriaga Eq. 35, Term 6 (Zmbij term)
                    r2_temp = linalg::triplet(S_ij_mj, tau[ij], S_ij_mj, true, false, false);
                    r2_temp->reshape(npno_mj * npno_mj, 1);
                    r2_temp = linalg::doublet(K_maef_[mj], r2_temp, false, false);
                    r2_temp = linalg::doublet(S_ij_mj, r2_temp, false, false);
                    C_DGER(npno_ij, npno_ij, -1.0, &(*temp_t1)(0,0), 1, &(*r2_temp)(0,0), 1, &(*Rn_iajb[ij])(0,0), npno_ij);

                    // Madriaga Eq. 35, Term 10
                    auto S_ii_mj = S_PNO(ii, mj);
                    auto T_i_temp = linalg::doublet(S_ii_mj, T_ia_[i], true, false);
                    r2_temp = linalg::triplet(T_i_temp, K_iajb_[mj], S_ij_mj, true, false, true);
                    C_DGER(npno_ij, npno_ij, -1.0, &(*temp_t1)(0,0), 1, &(*r2_temp)(0,0), 1, &(*Rn_iajb[ij])(0,0), npno_ij);

                    // Madriaga Eq. 35, Term 11
                    r2_temp = linalg::triplet(S_ij_mj, J_ijab_[mj], T_i_temp, false, false, false);
                    C_DGER(npno_ij, npno_ij, -1.0, &(*r2_temp)(0,0), 1, &(*temp_t1)(0,0), 1, &(*Rn_iajb[ij])(0,0), npno_ij);
                }

                if (n_pno_[im] != 0) {
                    // Madriaga Eq. 35, Term 3
                    auto S_im_ij = S_pno_ij_ik_[im][j];
                    r2_temp = linalg::triplet(S_im_ij, T_iajb_[im], S_im_ij, true, false, false);
                    int m_jj = lmopair_to_lmos_dense_[jj][m];
                    std::vector<int> m_jj_slice(1, m_jj);
                    double tf_dot = T_ia_[j]->vector_dot(submatrix_rows(*Fme[jj], m_jj_slice)->transpose());
                    r2_temp->scale((*Fmi)(m,j) + 0.5 * tf_dot);
                    Rn_iajb[ij]->subtract(r2_temp);
                }

                // Madriaga Eq. 35, Term 13
                r2_temp = submatrix_rows(*K_mbij_[ij], m_ij_slice)->transpose();
                C_DGER(npno_ij, npno_ij, -1.0, &(*temp_t1)(0,0), 1, &(*r2_temp)(0,0), 1, &(*Rn_iajb[ij])(0,0), npno_ij);

                if (n_pno_[mi] != 0 && n_pno_[mj] != 0) {
                    auto S_im_ij = S_pno_ij_ik_[im][j];
                    auto S_mj_mi = S_pno_ij_ik_[mj][i];
                    auto Wmbej_mj = linalg::triplet(S_ij_mj, Wmbej[mj], S_mj_mi, false, false, false);
                    auto Wmbje_mj = linalg::triplet(S_ij_mj, Wmbje[mj], S_mj_mi, false, false, false);
                    auto Wmbje_mi = linalg::triplet(S_ij_im, Wmbje[mi], S_mj_mi, false, false, true);

                    // Madriaga Eq. 35, Term 7
                    r2_temp = T_iajb_[im]->clone();
                    r2_temp->subtract(T_iajb_[im]->transpose());
                    Rn_iajb[ij]->add(linalg::triplet(S_im_ij, r2_temp, Wmbej_mj, true, false, true));

                    // Madriaga Eq. 35, Term 8
                    r2_temp = Wmbej_mj->clone();
                    r2_temp->add(Wmbje_mj);
                    Rn_iajb[ij]->add(linalg::triplet(S_im_ij, T_iajb_[im], r2_temp, true, false, true));

                    // Madriaga Eq. 35, Term 9
                    Rn_iajb[ij]->add(linalg::triplet(S_ij_mj, T_iajb_[mj], Wmbje_mi, false, false, true));
                }

                // Madriaga Eq. 35, Term 4
                for (int n_ij = 0; n_ij < lmopair_to_lmos_[ij].size(); n_ij++) {
                    int n = lmopair_to_lmos_[ij][n_ij];
                    int mn = i_j_to_ij_[m][n];

                    if (mn == -1 || n_pno_[mn] == 0) continue;

                    auto S_ij_mn = S_PNO(ij, mn);

                    int i_mn = lmopair_to_lmos_dense_[mn][i], j_mn = lmopair_to_lmos_dense_[mn][j];
                    r2_temp = linalg::triplet(S_ij_mn, tau[mn], S_ij_mn, false, false, true);
                    r2_temp->scale(0.5 * Wmnij[mn]->get(i_mn, j_mn));
                    Rn_iajb[ij]->add(r2_temp);
                }
            }
        }

        // Compute Doubles Residual from Non-Symmetrized Doubles Residual
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            int ji = ij_to_ji_[ij];
            R_iajb[ij] = std::make_shared<Matrix>("R_iajb", n_pno_[ij], n_pno_[ij]);

            if (n_pno_[ij] == 0) continue;

            R_iajb[ij] = Rn_iajb[ij]->clone();
            R_iajb[ij]->add(Rn_iajb[ji]->transpose());

            R_iajb_rms[ij] = R_iajb[ij]->rms();
        }

        // Update Singles Amplitude
#pragma omp parallel for
        for (int i = 0; i < naocc; ++i) {
            int ii = i_j_to_ij_[i][i];
            for (int a_ij = 0; a_ij < n_pno_[ii]; ++a_ij) {
                (*T_ia_[i])(a_ij, 0) -= (*R_ia[i])(a_ij, 0) / (e_pno_[ii]->get(a_ij) - F_lmo_->get(i,i));
            }
        }

        // Update Doubles Amplitude
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];
            int ii = i_j_to_ij_[i][i], jj = i_j_to_ij_[j][j];

            if (n_pno_[ij] == 0) continue;

            for (int a_ij = 0; a_ij < n_pno_[ij]; ++a_ij) {
                for (int b_ij = 0; b_ij < n_pno_[ij]; ++b_ij) {
                    (*T_iajb_[ij])(a_ij, b_ij) -= (*R_iajb[ij])(a_ij, b_ij) / 
                                (e_pno_[ij]->get(a_ij) + e_pno_[ij]->get(b_ij) - F_lmo_->get(i,i) - F_lmo_->get(j,j));
                }
            }
        }

        // DIIS Extrapolation
        auto T_ia_flat = flatten_mats(T_ia_);
        auto R_ia_flat = flatten_mats(R_ia);

        auto T_iajb_flat = flatten_mats(T_iajb_);
        auto R_iajb_flat = flatten_mats(R_iajb);

        if (iteration == 0) {
            diis1.set_error_vector_size(R_ia_flat.get());
            diis1.set_vector_size(T_ia_flat.get());

            diis2.set_error_vector_size(R_iajb_flat.get());
            diis2.set_vector_size(T_iajb_flat.get());
        }

        diis1.add_entry(R_ia_flat.get(), T_ia_flat.get());
        diis1.extrapolate(T_ia_flat.get());

        diis2.add_entry(R_iajb_flat.get(), T_iajb_flat.get());
        diis2.extrapolate(T_iajb_flat.get());

        copy_flat_mats(T_ia_flat, T_ia_);
        copy_flat_mats(T_iajb_flat, T_iajb_);

        // Update Special Doubles Amplitudes
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ij++) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];
            int ii = i_j_to_ij_[i][i], jj = i_j_to_ij_[j][j];

            if (n_pno_[ij] == 0) continue;

            Tt_iajb_[ij] = T_iajb_[ij]->clone();
            Tt_iajb_[ij]->scale(2.0);
            Tt_iajb_[ij]->subtract(T_iajb_[ij]->transpose());

            auto S_ii_ij = S_pno_ij_ik_[ii][j];
            auto S_jj_ij = S_pno_ij_kj_[jj][i];
            auto tia_temp = linalg::doublet(S_ii_ij, T_ia_[i], true, false);
            auto tjb_temp = linalg::doublet(S_jj_ij, T_ia_[j], true, false);

            for (int a_ij = 0; a_ij < n_pno_[ij]; ++a_ij) {
                for (int b_ij = 0; b_ij < n_pno_[ij]; ++b_ij) {
                    double t1_cont = tia_temp->get(a_ij, 0) * tjb_temp->get(b_ij, 0);
                    double t2_cont = T_iajb_[ij]->get(a_ij, b_ij);

                    tau[ij]->set(a_ij, b_ij, t2_cont + t1_cont);
                    tau_tilde[ij]->set(a_ij, b_ij, t2_cont + 0.5 * t1_cont);
                    tau_bar[ij]->set(a_ij, b_ij, 0.5 * t2_cont + t1_cont);
                }
            }
        }

        // evaluate convergence using current amplitudes and residuals
        e_prev = e_curr;
        // Compute LCCSD energy
        e_curr = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : e_curr)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            if (n_pno_[ij] == 0) continue;

            e_curr += tau[ij]->vector_dot(L_iajb_[ij]);
        }
        double r_curr1 = *max_element(R_ia_rms.begin(), R_ia_rms.end());
        double r_curr2 = *max_element(R_iajb_rms.begin(), R_iajb_rms.end());

        r_converged = (fabs(r_curr1) < options_.get_double("R_CONVERGENCE"));
        r_converged &= (fabs(r_curr2) < options_.get_double("R_CONVERGENCE"));
        e_converged = (fabs(e_curr - e_prev) < options_.get_double("E_CONVERGENCE"));

        outfile->Printf("  @LCCSD iter %3d: %16.12f %10.3e %10.3e %10.3e\n", iteration, e_curr, e_curr - e_prev, r_curr1, r_curr2);

        iteration++;

        if (iteration > max_iteration) {
            throw PSIEXCEPTION("Maximum DLPNO iterations exceeded.");
        }
    }

    e_lccsd_ = e_curr;

    timer_off("LCCSD iterations");
}

void DLPNOMP2::compute_triples_info() {
    tno_transform();
    compute_pno_tno_overlaps();
}

}  // namespace dlpno
}  // namespace psi
