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

#include "mp2.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/psi4-dec.h"
#include "psi4/physconst.h"
#include "psi4/psifiles.h"

#include "psi4/lib3index/3index.h"
#include "psi4/libfock/apps.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libfock/v.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/extern.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/local.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace psi {
namespace dlpnomp2 {

DLPNOMP2::DLPNOMP2(SharedWavefunction ref_wfn, Options& options) : Wavefunction(options) {
    shallow_copy(ref_wfn);
    reference_wavefunction_ = ref_wfn;

    common_init();
}
DLPNOMP2::~DLPNOMP2() {}
void DLPNOMP2::common_init() {
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");

    // if (options_.get_str("REFERENCE") == "ROHF" || options_.get_str("REFERENCE") == "CUHF")
    //    throw PSIEXCEPTION("SemiCanonical transform does not work at the moment");
    // reference_wavefunction_->semicanonicalize();

    // copy(reference_wavefunction_);
    name_ = "DLPNO-MP2";
    module_ = "dlpnomp2";

    //variables_["MP2 SINGLES ENERGY"] = 0.0;
    //variables_["MP2 DOUBLES ENERGY"] = 0.0;
    //variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = 0.0;
    //variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = 0.0;
    variables_["SCF TOTAL ENERGY"] = reference_wavefunction_->energy();

    //sss_ = options_.get_double("MP2_SS_SCALE");
    //oss_ = options_.get_double("MP2_OS_SCALE");

    ribasis_ = get_basisset("DF_BASIS_MP2");

}


void print_aux_domains(const std::vector<std::vector<int>> &lmo_to_ribfs, const std::vector<std::vector<int>> &lmo_to_riatoms) {

    size_t total_atoms = 0, min_atoms = lmo_to_riatoms[0].size(), max_atoms = 0;
    for (auto atom_list : lmo_to_riatoms) {
        total_atoms += atom_list.size();
        min_atoms = std::min(min_atoms, atom_list.size());
        max_atoms = std::max(max_atoms, atom_list.size());
    }

    size_t total_bfs = 0, min_bfs = lmo_to_ribfs[0].size(), max_bfs = 0;
    for (auto bf_list : lmo_to_ribfs) {
        total_bfs += bf_list.size();
        min_bfs = std::min(min_bfs, bf_list.size());
        max_bfs = std::max(max_bfs, bf_list.size());
    }

    size_t naocc = lmo_to_ribfs.size();
    outfile->Printf("\n");
    outfile->Printf("    Auxiliary BFs per Local MO:\n");
    outfile->Printf("      Average = %4d AUX BFs (%d atoms)\n", total_bfs / naocc, total_atoms / naocc);
    outfile->Printf("      Min     = %4d AUX BFs (%d atoms)\n", min_bfs, min_atoms);
    outfile->Printf("      Max     = %4d AUX BFs (%d atoms)\n", max_bfs, max_atoms);

}

void print_aux_pair_domains(const std::vector<std::vector<int>> &lmopair_to_ribfs, const std::vector<std::vector<int>> &lmopair_to_riatoms) {

    int n_lmo_pairs = lmopair_to_ribfs.size();
    int min_domain_ri = lmopair_to_ribfs[0].size(), max_domain_ri = 0, total_domain_ri = 0;
    int min_domain_ri_atom = lmopair_to_riatoms[0].size(), max_domain_ri_atom = 0, total_domain_ri_atom = 0;
    for(size_t ij = 0; ij < n_lmo_pairs; ij++) {
        int pair_domain_size_ri = lmopair_to_ribfs[ij].size();
        int pair_domain_size_ri_atom = lmopair_to_riatoms[ij].size();

        total_domain_ri += pair_domain_size_ri;
        total_domain_ri_atom += pair_domain_size_ri_atom;

        min_domain_ri = std::min(min_domain_ri, pair_domain_size_ri);
        min_domain_ri_atom = std::min(min_domain_ri_atom, pair_domain_size_ri_atom);

        max_domain_ri = std::max(max_domain_ri, pair_domain_size_ri);
        max_domain_ri_atom = std::max(max_domain_ri_atom, pair_domain_size_ri_atom);
    }

    outfile->Printf("  \n");
    outfile->Printf("    Auxiliary BFs per Local MO pair:\n");
    outfile->Printf("      Average = %4d AUX BFs (%d atoms)\n", total_domain_ri / n_lmo_pairs, total_domain_ri_atom / n_lmo_pairs);
    outfile->Printf("      Min     = %4d AUX BFs (%d atoms)\n", min_domain_ri, min_domain_ri_atom);
    outfile->Printf("      Max     = %4d AUX BFs (%d atoms)\n", max_domain_ri, max_domain_ri_atom);

}

void print_pao_domains(const std::vector<std::vector<int>> &lmo_to_paos, const std::vector<std::vector<int>> &lmo_to_paoatoms) {

    size_t total_atoms = 0, min_atoms = lmo_to_paoatoms[0].size(), max_atoms = 0;
    for (auto atom_list : lmo_to_paoatoms) {
        total_atoms += atom_list.size();
        min_atoms = std::min(min_atoms, atom_list.size());
        max_atoms = std::max(max_atoms, atom_list.size());
    }

    size_t total_paos = 0, min_paos = lmo_to_paos[0].size(), max_paos = 0;
    for (auto pao_list : lmo_to_paos) {
        total_paos += pao_list.size();
        min_paos = std::min(min_paos, pao_list.size());
        max_paos = std::max(max_paos, pao_list.size());
    }

    size_t naocc = lmo_to_paos.size();
    outfile->Printf("  \n");
    outfile->Printf("    Projected AOs per Local MO:\n");
    outfile->Printf("      Average = %4d PAOs (%d atoms)\n", total_paos / naocc, total_atoms / naocc);
    outfile->Printf("      Min     = %4d PAOs (%d atoms)\n", min_paos, min_atoms);
    outfile->Printf("      Max     = %4d PAOs (%d atoms)\n", max_paos, max_atoms);

}

void print_pao_pair_domains(const std::vector<std::vector<int>> &lmopair_to_paos, const std::vector<std::vector<int>> &lmopair_to_paoatoms) {

    //int n_lmo_pairs = lmopair_to_paoatoms.size();
    //int min_domain_pre = X_pao[0]->rowspi(0), max_domain_pre = 0, total_domain_pre = 0;
    //int min_domain = X_pao[0]->colspi(0), max_domain = 0, total_domain = 0;
    //for(size_t ij = 0; ij < n_lmo_pairs; ij++) {
    //    int pair_domain_pre_size = X_pao[ij]->rowspi(0);
    //    int pair_domain_size = X_pao[ij]->colspi(0);

    //    total_domain_pre += pair_domain_pre_size;
    //    total_domain += pair_domain_size;

    //    min_domain_pre = std::min(min_domain_pre, pair_domain_pre_size);
    //    min_domain = std::min(min_domain, pair_domain_size);

    //    max_domain_pre = std::max(max_domain_pre, pair_domain_pre_size);
    //    max_domain = std::max(max_domain, pair_domain_size);
    //}

    //outfile->Printf("  \n");
    //outfile->Printf("    Projected AOs per Local MO pair: (Before removing lindeps)\n");
    //outfile->Printf("      Average = (%4d -> %4d) PAOs\n", total_domain_pre / n_lmo_pairs, total_domain / n_lmo_pairs);
    //outfile->Printf("      Min     = (%4d -> %4d) PAOs\n", min_domain_pre, min_domain);
    //outfile->Printf("      Max     = (%4d -> %4d) PAOs\n", max_domain_pre, max_domain);


    int n_lmo_pairs = lmopair_to_paos.size();
    int min_domain_pao = lmopair_to_paos[0].size(), max_domain_pao = 0, total_domain_pao = 0;
    int min_domain_atom = lmopair_to_paoatoms[0].size(), max_domain_atom = 0, total_domain_atom = 0;
    for(size_t ij = 0; ij < n_lmo_pairs; ij++) {
        int pair_domain_size_pao = lmopair_to_paos[ij].size();
        int pair_domain_size_atom = lmopair_to_paoatoms[ij].size();

        total_domain_pao += pair_domain_size_pao;
        total_domain_atom += pair_domain_size_atom;

        min_domain_pao = std::min(min_domain_pao, pair_domain_size_pao);
        min_domain_atom = std::min(min_domain_atom, pair_domain_size_atom);

        max_domain_pao = std::max(max_domain_pao, pair_domain_size_pao);
        max_domain_atom = std::max(max_domain_atom, pair_domain_size_atom);
    }

    outfile->Printf("  \n");
    outfile->Printf("    Projected AOs per Local MO pair (possibly linearly dependent):\n");
    outfile->Printf("      Average = %4d AUX BFs (%d atoms)\n", total_domain_pao / n_lmo_pairs, total_domain_atom / n_lmo_pairs);
    outfile->Printf("      Min     = %4d AUX BFs (%d atoms)\n", min_domain_pao, min_domain_atom);
    outfile->Printf("      Max     = %4d AUX BFs (%d atoms)\n", max_domain_pao, max_domain_atom);

}

void print_lmo_domains(const std::vector<std::vector<int>> &i_j_to_ij, int exclude_pairs_overlap, int exclude_pairs_energy, double de_dipole) {

    int naocc = i_j_to_ij.size();
    int total_lmos = 0, min_lmos = naocc, max_lmos = 0;
    for(size_t i = 0; i < naocc; i++) {
        int lmos = 0;
        for(size_t j = 0; j < naocc; ++j) {
            if(i_j_to_ij[i][j] != -1) {
                lmos += 1;
            }
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
    outfile->Printf("    Screened %d of %d LMO pairs (%.2f %%)\n", naocc * naocc - total_lmos, naocc * naocc,  100.0 - (total_lmos * 100.0) / (naocc * naocc) );
    outfile->Printf("             %d pairs met overlap criteria\n", exclude_pairs_overlap);
    outfile->Printf("             %d pairs met energy criteria\n", exclude_pairs_energy);
    outfile->Printf(" \n");
    outfile->Printf("    Screened LMO pair energy =  %.12f \n", de_dipole);
}

void print_integral_sparsity(const std::vector<std::vector<int>> &riatom_to_shells1,
                             const std::vector<std::vector<int>> &riatom_to_shells2,
                             const std::vector<std::vector<int>> &atom_to_rishell,
                             const std::vector<std::vector<int>> &riatom_to_lmos_ext,
                             const std::vector<std::vector<int>> &riatom_to_paos_ext,
                             const std::vector<std::vector<int>> &atom_to_ribf,
                             int nshell,
                             int naocc,
                             int nbf,
                             int naux) {

    // statistics for number of (MN|K) shell triplets we need to compute

    size_t triplets = 0; // computed (MN|K) triplets with no screening
    size_t triplets_lmo = 0;  // computed (MN|K) triplets with only LMO screening
    size_t triplets_pao = 0;  // computed (MN|K) triplets with only PAO screening
    size_t triplets_both = 0; // computed (MN|K) triplets with LMO and PAO screening

    for(size_t atom = 0; atom < riatom_to_shells1.size(); atom++) {
        size_t nshellri_atom = atom_to_rishell[atom].size();
        triplets += nshell * nshell * nshellri_atom;
        triplets_lmo += riatom_to_shells1[atom].size() * nshell * nshellri_atom;
        triplets_pao += nshell * riatom_to_shells2[atom].size() * nshellri_atom;
        triplets_both += riatom_to_shells1[atom].size() * riatom_to_shells2[atom].size() * nshellri_atom;
    }
    size_t screened_total = triplets - triplets_both;
    size_t screened_lmo = triplets - triplets_lmo;
    size_t screened_pao = triplets - triplets_pao;

    // statistics for the number of (iu|Q) integrals we're left with after the transformation

    size_t total_integrals = (size_t) naocc * nbf * naux;
    size_t actual_integrals = 0;

    for(size_t atom = 0; atom < riatom_to_shells1.size(); atom++) {
        actual_integrals += riatom_to_lmos_ext[atom].size() * riatom_to_paos_ext[atom].size() * atom_to_ribf[atom].size();
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
    outfile->Printf("    Required memory: %.3f GiB (%.2f %% reduction from dense format) \n", actual_memory, screened_memory * 100.0 / total_memory);

}

std::vector<int> merge_lists(const std::vector<int> &l1, const std::vector<int> &l2) {
    /* Args: sorted lists l1 and l2
     * Returns: sorted union of l1 and l2
     */

    std::vector<int> l12;

    int i1 = 0, i2 = 0;
    while(i1 < l1.size() || i2 < l2.size()) {
        if(i1 == l1.size()) {
            l12.push_back(l2[i2]);
            ++i2;
        } else if(i2 == l2.size()) {
            l12.push_back(l1[i1]);
            ++i1;
        } else if(l1[i1] == l2[i2]) {
            l12.push_back(l1[i1]);
            ++i1;
            ++i2;
        } else if(l1[i1] < l2[i2]) {
            l12.push_back(l1[i1]);
            ++i1;
        } else {
            l12.push_back(l2[i2]);
            ++i2;
        }
    }

    return l12;

}


std::vector<int> contract_lists(const std::vector<int> &y, const std::vector<std::vector<int>> &A_to_y) {
    /* Args: sorted list of values y, sparse map from A to y (assume sorted)
     * Returns: sorted list yA
     *
     * For all a, every value in A_to_y[a] is included in yA if at least one is present in y
     */

    // TODO: runtime is proportional to A_to_y size (system size, O(N))
    // could maybe reduce to &y size (domain size, O(1)), probably doesn't matter

    std::vector<int> yA;

    for(int a = 0, y_ind = 0; a < A_to_y.size(); ++a) {

        bool is_a = false;
        for(auto y_val : A_to_y[a]) {
            if (y_ind < y.size() && y[y_ind] == y_val) {
                y_ind++;
                is_a = true;
            }
        }

        if(is_a) {
            for(auto y_val : A_to_y[a]) {
                yA.push_back(y_val);
            }
        }

    }

    return yA;

}


std::vector<int> block_list(const std::vector<int> &x_list, const std::vector<int> &x_to_y_map) {
    /* Args: x is a list of values (sorted), y is a map from values of x to values of y
     * Returns: a list of y values
     *
     * Multiple values in x may map to the same value in y (i.e. x is a list of bf, y is atoms)
     */

    std::vector<int> y_list;

    for(int x_val : x_list) {
        int y_val = x_to_y_map[x_val];
        if(y_list.size() == 0) {
            y_list.push_back(y_val);
        } else if(y_list[y_list.size() - 1] != y_val) {
            y_list.push_back(y_val);
        }
    }

    return y_list;

}



std::vector<std::vector<int>> invert_map(const std::vector<std::vector<int>> &x_to_y, int ny) {
    /* Args: 
     * Returns: 
     *
     *
     */

    int nx = x_to_y.size();
    std::vector<std::vector<int>> y_to_x(ny);

    for(int x = 0; x < nx; x++) {
        for(auto y : x_to_y[x]) {
            y_to_x[y].push_back(x);
        }
    }

    return y_to_x;

}



std::vector<std::vector<int>> chain_maps(const std::vector<std::vector<int>> &x_to_y, const std::vector<std::vector<int>> &y_to_z) {
    /* Args: 
     * Returns: 
     *
     *
     */

    int nx = x_to_y.size();
    std::vector<std::vector<int>> x_to_z(nx);

    for(int x = 0; x < nx; x++) {
        for(auto y : x_to_y[x]) {
            for(auto z : y_to_z[y]) {
                //if(x_to_z[x].size() == 0) {
                    x_to_z[x].push_back(z);
                //} else if(x_to_z[x][x_to_z[x].size() - 1] != z) {
                    x_to_z[x].push_back(z);
                //}
            }
        }
        std::sort(x_to_z[x].begin(), x_to_z[x].end());
        x_to_z[x].erase(std::unique(x_to_z[x].begin(), x_to_z[x].end()), x_to_z[x].end());

        for(auto z: x_to_z[x]) {
            //outfile->Printf(" %d", z);
        }
    }

    return x_to_z;

}



std::vector<std::vector<int>> extend_maps(const std::vector<std::vector<int>> &i_to_y, const std::vector<std::pair<int,int>> &ij_to_i_j) {
    /* Args: 
     * Returns: 
     *
     *
     */

    int ni = i_to_y.size();
    std::vector<std::vector<int>> iext_to_y(ni);

    for(auto pair : ij_to_i_j) {
        size_t i, j;
        std::tie(i,j) = pair;
        iext_to_y[i] = merge_lists(iext_to_y[i], i_to_y[j]);
    }

    return iext_to_y;

}


void normalize_cols(SharedMatrix C, SharedMatrix S) {

    // could do this way more efficiently (we only need the diagonals)
    SharedMatrix C_norm = linalg::triplet(C, S, C, true, false, false);
    for(size_t i = 0; i < C->colspi(0); ++i) {
        C->scale_column(0, i, pow(C_norm->get(i,i), -0.5));
    }
    return;
}

std::pair<SharedMatrix, SharedVector> get_canonicalizer(SharedMatrix C, SharedMatrix F) {
    /* Args: orthonormal orbitals C (ao x mo) and fock matrix F (ao x ao)
     * Return: canonical transformation matrix U (mo x mo) and energy vector e (mo) 
     *
     * U and e are defined as: FCU = eCU (CU is canonical)
     */

    SharedMatrix U = std::make_shared<Matrix>("eigenvectors", C->colspi(0), C->colspi(0));
    SharedVector e = std::make_shared<Vector>("eigenvalues", C->colspi(0));

    auto temp = linalg::triplet(C, F, C, true, false, false);
    temp->diagonalize(U, e, descending);

    return std::make_pair(U, e);

}

SharedMatrix get_orthogonalizer(SharedMatrix C, SharedMatrix S, Options &options) {
    /* Args: normalized orbitals C (ao x mo) and overlap matrix S (ao x ao)
     * Return: orthogonalization matrix X (mo x mo_new)
     *
     * X is defined as: SCX = CX (CX is orthogonal)
     * linear dependencies are removed with keyword S_CUT, so (mo_new <= mo)
     */

    int nbf = S->rowspi(0);
    int nmo_initial = C->colspi(0);
    int nmo_final = C->colspi(0);

    SharedMatrix X = std::make_shared<Matrix>("eigenvectors", C->colspi(0), C->colspi(0));
    SharedVector n = std::make_shared<Vector>("eigenvalues", C->colspi(0));

    auto S_mo = linalg::triplet(C, S, C, true, false, false);
    S_mo->diagonalize(X, n, descending);

    for(size_t i = 0; i < nmo_initial; ++i) {
        if (fabs(n->get(i)) < options.get_double("S_CUT")) {
            nmo_final -= 1;
        }
    }

    //outfile->Printf("  orthogonalizer: %d to %d \n", nmo_initial, nmo_final);

    Dimension zero = Dimension(1);
    Dimension dim_final = Dimension(1);
    dim_final.fill(nmo_final);

    // Can I get a block of X without having to make these dimension objects?
    X = X->get_block({zero, X->rowspi()}, {zero, dim_final});
    n = n->get_block({zero, dim_final});

    S_mo = linalg::triplet(X, S_mo, X, true, false, false);

    // Can/should I do this scaling earlier? Can/should I also scale eigenvalues?
    for(size_t i = 0; i < nmo_final; ++i) {
        X->scale_column(0, i, pow(S_mo->get(i,i), -0.5));
    }

    return X;

}


std::pair<SharedMatrix, SharedVector> get_orthocanonicalizer(SharedMatrix Smo, SharedMatrix Fmo, Options &options) {
    /* Args: normalized orbitals C (ao x mo), overlap matrix S (ao x ao), fock matrix F (ao x ao)
     * Return: transformation matrix X (mo x mo_new) and energy vector e (mo_new)
     *
     * X is defined as: SCX = CX (CX is orthogonal) AND FCX = eCX (CX is canonical, w/ energies e)
     * linear dependencies are removed with keyword S_CUT, so (mo_new <= mo)
     */

    int nmo_initial = Smo->colspi(0);
    int nmo_final = Smo->colspi(0);

    SharedMatrix X = std::make_shared<Matrix>("eigenvectors", nmo_initial, nmo_initial);
    SharedVector n = std::make_shared<Vector>("eigenvalues", nmo_initial);

    Smo->diagonalize(X, n, descending);
    //n->print_out();

    for(size_t i = 0; i < nmo_initial; ++i) {
        if (fabs(n->get(i)) < options.get_double("S_CUT")) {
            nmo_final -= 1;
        }
    }

    //outfile->Printf("  orthogonalizer: %d to %d \n", nmo_initial, nmo_final);

    Dimension zero = Dimension(1);
    Dimension dim_final = Dimension(1);
    dim_final.fill(nmo_final);

    // Can I get a block of X without having to make these dimension objects?
    X = X->get_block({zero, X->rowspi()}, {zero, dim_final});
    n = n->get_block({zero, dim_final});

    auto Smo_orth = linalg::triplet(X, Smo, X, true, false, false);

    // Can/should I do this scaling earlier? Can/should I also scale eigenvalues?
    for(size_t i = 0; i < nmo_final; ++i) {
        X->scale_column(0, i, pow(Smo_orth->get(i,i), -0.5));
    }

    SharedMatrix U = std::make_shared<Matrix>("eigenvectors", nmo_final, nmo_final);
    SharedVector e = std::make_shared<Vector>("eigenvalues", nmo_final);

    auto Fmo_orth = linalg::triplet(X, Fmo, X, true, false, false);
    Fmo_orth->diagonalize(U, e, descending);

    X = linalg::doublet(X, U, false, false);

    return std::make_pair(X, e);

}

SharedMatrix get_rows(SharedMatrix mat, const std::vector<int> &row_inds) {

    SharedMatrix mat_new = std::make_shared<Matrix>("blah", row_inds.size(), mat->colspi(0));
    for(int r_new = 0; r_new < row_inds.size(); r_new++) {
        int r_old = row_inds[r_new];
        for(int c = 0; c < mat->colspi(0); c++) {
            mat_new->set(r_new, c, mat->get(r_old, c));
        }
    }
    return mat_new;
}

SharedMatrix get_cols(SharedMatrix mat, const std::vector<int> &col_inds) {

    SharedMatrix mat_new = std::make_shared<Matrix>("blah", mat->rowspi(0), col_inds.size());
    for(int r = 0; r < mat->rowspi(0); r++) {
        for(int c_new = 0; c_new < col_inds.size(); c_new++) {
            int c_old = col_inds[c_new];
            mat_new->set(r, c_new, mat->get(r, c_old));
        }
    }
    return mat_new;
}

SharedMatrix get_rows_and_cols(SharedMatrix mat, const std::vector<int> &row_inds, const std::vector<int> &col_inds) {

    SharedMatrix mat_new = std::make_shared<Matrix>("blah", row_inds.size(), col_inds.size());
    for(int r_new = 0; r_new < row_inds.size(); r_new++) {
        int r_old = row_inds[r_new];
        for(int c_new = 0; c_new < col_inds.size(); c_new++) {
            int c_old = col_inds[c_new];
            mat_new->set(r_new, c_new, mat->get(r_old, c_old));
        }
    }
    return mat_new;
}


std::pair<SharedMatrix,SharedMatrix> calculate_dois(SharedWavefunction ref_wfn, Options& options, SharedMatrix C_lmo, SharedMatrix C_pao)
{
    timer_on("Differential Overlap Integrals");

    std::shared_ptr<Molecule> mol = ref_wfn->molecule();
    std::shared_ptr<BasisSet> basis = ref_wfn->basisset();

    int nbf = basis->nbf();
    int naocc = C_lmo->colspi(0);
    int npao = C_pao->colspi(0); // same as nbf

    SharedMatrix DOI_ij = SharedMatrix(new Matrix("(i,i) Differential Overlap Integral", naocc, naocc));
    SharedMatrix DOI_iu = SharedMatrix(new Matrix("(i,u) Differential Overlap Integral", naocc, nbf));

    timer_on("Construct Grid");
    std::shared_ptr<DFTGrid> grid = std::make_shared<DFTGrid>(mol, basis, options);

    size_t nthread = 1;
#ifdef _OPENMP
    nthread = omp_get_max_threads();
#endif
    std::vector<std::shared_ptr<BasisFunctions>> point_funcs(nthread);
    std::vector<SharedMatrix> DOI_ij_temps(nthread);
    std::vector<SharedMatrix> DOI_iu_temps(nthread);
    for(size_t thread = 0; thread < nthread; thread++) {
        point_funcs[thread] = std::make_shared<BasisFunctions>(basis, grid->max_points(), basis->nbf());
        DOI_ij_temps[thread] = SharedMatrix(new Matrix("(i,i) Differential Overlap Integral", naocc, naocc));
        DOI_iu_temps[thread] = SharedMatrix(new Matrix("(i,u) Differential Overlap Integral", naocc, nbf));
    }
    
    timer_off("Construct Grid");

    // TODO: parallelize integration
    timer_on("Integration");
#pragma omp parallel for schedule(static, 1)    
    for (size_t Q = 0; Q < grid->blocks().size(); Q++) {

        size_t thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        std::shared_ptr<BlockOPoints> block = grid->blocks()[Q];
        int nbf_block = block->local_nbf();
        int npoints_block = block->npoints();

        // compute values of each basis function at each point in this block
        point_funcs[thread]->compute_functions(block);

        // the values we just computed (max_points x max_functions)
        SharedMatrix point_values = point_funcs[thread]->basis_values()["PHI"]; 
         
        std::vector<int> bf_map = block->functions_local_to_global();

        // resize point_values buffer to size of this block
        SharedMatrix point_values_trim = std::make_shared<Matrix>("blah", npoints_block, nbf_block); // points x bf_block
        for (size_t p = 0; p < npoints_block; p++) {
            for(size_t k = 0; k < nbf_block; k++) {
                point_values_trim->set(p,k,point_values->get(p,k));
            }
        }

        SharedMatrix C_lmo_slice = get_rows(C_lmo, bf_map); //bf_block x naocc
        SharedMatrix C_pao_slice = get_rows(C_pao, bf_map); //bf_block x npao

        // value of mo at each point squared
        C_lmo_slice = linalg::doublet(point_values_trim, C_lmo_slice, false, false); // points x naocc
        C_pao_slice = linalg::doublet(point_values_trim, C_pao_slice, false, false); // points x npao

        for (size_t p = 0; p < npoints_block; p++) {
            for (size_t i = 0; i < naocc; ++i) {
                C_lmo_slice->set(p, i, pow(C_lmo_slice->get(p, i), 2));
            }
            for (size_t u = 0; u < npao; ++u) {
                C_pao_slice->set(p, u, pow(C_pao_slice->get(p, u), 2));
            }
        }

        SharedMatrix C_lmo_slice_w = std::make_shared<Matrix>(C_lmo_slice); // points x naocc
        for (size_t p = 0; p < npoints_block; p++) {
            C_lmo_slice_w->scale_row(0, p, block->w()[p]);
        }

        DOI_ij_temps[thread]->add(linalg::doublet(C_lmo_slice_w, C_lmo_slice, true, false)); // naocc x naocc
        DOI_iu_temps[thread]->add(linalg::doublet(C_lmo_slice_w, C_pao_slice, true, false)); // naocc x npao

    }
    timer_off("Integration");

    for(size_t thread = 0; thread < nthread; thread++) {
        DOI_ij->add(DOI_ij_temps[thread]);
        DOI_iu->add(DOI_iu_temps[thread]);
    }

    for(size_t i = 0; i < naocc; i++) {
        for(size_t j = 0; j < naocc; j++) {
            DOI_ij->set(i,j,sqrt(DOI_ij->get(i,j)));
        }
        for(size_t u = 0; u < nbf; u++) {
            DOI_iu->set(i,u,sqrt(DOI_iu->get(i,u)));
        }
    }

    timer_off("Differential Overlap Integrals");

    return std::make_pair(DOI_ij, DOI_iu);
}


std::tuple<SharedMatrix, SharedMatrix, SharedMatrix> dipole_approx(SharedWavefunction ref_wfn, 
                                                                   Options& options, 
                                                                   SharedMatrix DOI_iu,
                                                                   SharedMatrix C_lmo, 
                                                                   SharedVector e_lmo,
                                                                   SharedMatrix C_pao,
                                                                   SharedMatrix S_pao,
                                                                   SharedMatrix F_pao)
{

    timer_on("Pair Energy Estimation");
    std::shared_ptr<Molecule> mol = ref_wfn->molecule();
    std::shared_ptr<BasisSet> basis = ref_wfn->basisset();

    int natom = mol->natom();
    int naocc = C_lmo->colspi(0);
    int nbf = C_lmo->rowspi(0);

    timer_on("<m|r|n>");
    std::shared_ptr<MintsHelper> mints = std::make_shared<MintsHelper>(basis, options);
    std::vector<SharedMatrix> ao_dipole = mints->ao_dipole();
    timer_off("<m|r|n>");

    std::vector<Vector3> R_i;

    // probably a BLAS-ier way to do this
    for (size_t i = 0; i < naocc; ++i) {
        double rx = 0, ry = 0, rz = 0;
        for (size_t u = 0; u < nbf; ++u) {
            for (size_t v = 0; v < nbf; ++v) {
                rx += C_lmo->get(u,i) * ao_dipole[0]->get(u,v) * C_lmo->get(v,i);
                ry += C_lmo->get(u,i) * ao_dipole[1]->get(u,v) * C_lmo->get(v,i);
                rz += C_lmo->get(u,i) * ao_dipole[2]->get(u,v) * C_lmo->get(v,i);
            }
        }
        R_i.push_back(Vector3(rx, ry, rz));
    }

    SharedMatrix e_actual = std::make_shared<Matrix>("Dipole SC MP2 Energies", naocc, naocc);
    SharedMatrix e_linear = std::make_shared<Matrix>("Parallel Dipole SC MP2 Energies", naocc, naocc);

    timer_on("<i|r|n>");
    SharedMatrix lmo_bf_dx = linalg::doublet(C_lmo, ao_dipole[0], true, false);
    SharedMatrix lmo_bf_dy = linalg::doublet(C_lmo, ao_dipole[1], true, false);
    SharedMatrix lmo_bf_dz = linalg::doublet(C_lmo, ao_dipole[2], true, false);
    timer_off("<i|r|n>");

    std::vector<std::vector<Vector3>> lmo_pao_dr(naocc);
    std::vector<SharedVector> lmo_pao_e(naocc);

    std::vector<std::vector<int>> atom_to_bf(natom);
    for(size_t i = 0; i < nbf; ++i) {
        atom_to_bf[basis->function_to_center(i)].push_back(i);
    }

    timer_on("<i|r|u>");
    for (size_t i = 0; i < naocc; ++i) {

        std::vector<int> pao_inds;
        for(size_t u = 0; u < nbf; u++) {
            if(fabs(DOI_iu->get(i, u)) > options.get_double("T_CUT_DO_PRE")) {
                pao_inds.push_back(u);
            }
        }
        pao_inds = contract_lists(pao_inds, atom_to_bf);;
        
        SharedMatrix C_pao_i = get_cols(C_pao, pao_inds);
        SharedMatrix S_pao_i = get_rows_and_cols(S_pao, pao_inds, pao_inds);
        SharedMatrix F_pao_i = get_rows_and_cols(F_pao, pao_inds, pao_inds);

        SharedMatrix X_pao_i;
        SharedVector e_pao_i;
        std::tie(X_pao_i, e_pao_i) = get_orthocanonicalizer(S_pao_i, F_pao_i, options);
        C_pao_i = linalg::doublet(C_pao_i, X_pao_i, false, false); // now in a nonredundant basis

        int npao_i_new = X_pao_i->colspi(0);
        //outfile->Printf("  LMO %d has domain size of %d / %d AOs (%d / %d VIR) \n", i, npao_i, nbf, npao_i_new,nbf-naocc);

        for(size_t u = 0; u < npao_i_new; u++) {
            double dx_iu = 0.0;
            double dy_iu = 0.0;
            double dz_iu = 0.0;
            for(size_t v = 0; v < nbf; v++) {
                dx_iu += lmo_bf_dx->get(i,v) * C_pao_i->get(v,u);
                dy_iu += lmo_bf_dy->get(i,v) * C_pao_i->get(v,u);
                dz_iu += lmo_bf_dz->get(i,v) * C_pao_i->get(v,u);
            }
            lmo_pao_dr[i].push_back(Vector3(dx_iu, dy_iu, dz_iu));
        }
        lmo_pao_e[i] = e_pao_i;

    }
    timer_off("<i|r|u>");

    SharedMatrix distances = std::make_shared<Matrix>("blah", naocc, naocc);
    timer_on("e_ij");
    for (size_t i = 0; i < naocc; ++i) {
        for (size_t j = i+1; j < naocc; ++j) {

            Vector3 R_ij = R_i[i] - R_i[j];
            Vector3 Rh_ij = R_ij / R_ij.norm();
            distances->set(i, j, R_ij.norm());
            distances->set(j, i, R_ij.norm());

            double e_actual_temp = 0.0;
            double e_linear_temp = 0.0;

            for (int u = 0; u < lmo_pao_dr[i].size(); u++) {
                for (int v = 0; v < lmo_pao_dr[j].size(); v++) {

                    Vector3 iu = lmo_pao_dr[i][u];
                    Vector3 jv = lmo_pao_dr[j][v];

                    double num_actual = iu.dot(jv) - 3 * (iu.dot(Rh_ij) * jv.dot(Rh_ij)); 
                    num_actual *= num_actual;

                    double num_linear = -2 * iu.dot(jv);
                    num_linear *= num_linear;

                    double denom = (lmo_pao_e[i]->get(u) + lmo_pao_e[j]->get(v)) - (e_lmo->get(i) + e_lmo->get(j));

                    if (denom < 0.0) {
                        outfile->Printf("  ERROR: denom=%f \n", denom);
                    }

                    e_actual_temp += (num_actual / denom);
                    e_linear_temp += (num_linear / denom);
              }
            }

            e_actual_temp *= (-4 * pow(R_ij.norm(), -6));
            e_linear_temp *= (-4 * pow(R_ij.norm(), -6));

            e_actual->set(i, j, e_actual_temp);
            e_actual->set(j, i, e_actual_temp);

            e_linear->set(i, j, e_linear_temp);
            e_linear->set(j, i, e_linear_temp);

        }
    }
    timer_off("e_ij");

    timer_off("Pair Energy Estimation");
                      
    return std::make_tuple(e_actual, e_linear, distances);
}



std::tuple<SharedMatrix, SharedMatrix, SharedMatrix, SharedVector, double> pno(SharedMatrix F, SharedMatrix K_ij, SharedVector e_vir_ij, double ei, double ej, Options &options) {

    size_t nvir_ij = K_ij->rowspi(0);

    // Build amplitudes for this pair of LMOs
    SharedMatrix T_ij = std::make_shared<Matrix>("blah", nvir_ij, nvir_ij);
    for(size_t a = 0; a < nvir_ij; a++) {
        for(size_t b = 0; b < nvir_ij; b++) {
            T_ij->set(a, b,  -1.0 * K_ij->get(a,b) / (e_vir_ij->get(a) + e_vir_ij->get(b) + -ei + -ej));
        }
    }

    SharedMatrix Tt_ij = T_ij->clone();
    Tt_ij->scale(2.0);
    Tt_ij->subtract(T_ij->transpose());

    // mp2 energy of this LMO pair before transformation to PNOs
    double e_mp2_vir = K_ij->vector_dot(Tt_ij);
    
    // Construct pair density from amplitudes
    SharedMatrix D_ij = linalg::doublet(Tt_ij, T_ij, false, true);
    D_ij->add(linalg::doublet(Tt_ij, T_ij, true, false));

    // Diagonalization of pair density gives PNOs (in basis of the LMO's virtual domain) and PNO occ numbers
    SharedMatrix vir_to_pno = std::make_shared<Matrix>("eigenvectors", nvir_ij, nvir_ij);
    SharedVector pno_occ = std::make_shared<Vector>("eigenvalues", nvir_ij);
    D_ij->diagonalize(vir_to_pno, pno_occ, descending);

    int nvir_ij_final= 0;
    for(size_t a = 0; a < nvir_ij; ++a) {
        if (fabs(pno_occ->get(a)) >= options.get_double("T_CUT_PNO")) { // TODO: can this be negative? Or is that just a precision thing?
            nvir_ij_final++;
        }
    }
    //outfile->Printf("  %d PNOs for pair (reduced from %d) \n", nvir_ij_final, nvir_ij);

    Dimension zero = Dimension(1);
    Dimension dim_final = Dimension(1);
    dim_final.fill(nvir_ij_final);

    // This transformation gives orbitals that are orthonormal but not canonical
    vir_to_pno = vir_to_pno->get_block({zero, vir_to_pno->rowspi()}, {zero, dim_final});
    pno_occ = pno_occ->get_block({zero, dim_final});

    SharedMatrix pno_canon;
    SharedVector e_ij_pno;
    std::tie(pno_canon, e_ij_pno) = get_canonicalizer(vir_to_pno, F);

    // This transformation gives orbitals that are orthonormal and canonical
    vir_to_pno = linalg::doublet(vir_to_pno, pno_canon, false, false);

    SharedMatrix K_ij_pno = linalg::triplet(vir_to_pno, K_ij, vir_to_pno, true, false, false);
    SharedMatrix T_ij_pno = linalg::triplet(vir_to_pno, T_ij, vir_to_pno, true, false, false);

    SharedMatrix Tt_ij_pno = T_ij_pno->clone();
    Tt_ij_pno->scale(2.0);
    Tt_ij_pno->subtract(T_ij_pno->transpose());

    // mp2 energy of this LMO pair after transformation to PNOs
    double e_mp2_pno = K_ij_pno->vector_dot(Tt_ij_pno);

    //outfile->Printf("  PNO truncation energy error of %16.12f (from %16.12f to %16.12f) \n", e_mp2_vir - e_mp2_pno, e_mp2_vir, e_mp2_pno);

    return std::make_tuple(K_ij_pno, T_ij_pno, vir_to_pno, e_ij_pno, e_mp2_vir - e_mp2_pno);

}

double mp2_energy(std::vector<SharedMatrix> &Tt_iajb, std::vector<SharedMatrix> &K_iajb) {
    double e_mp2 = 0.0;
    for (int ij = 0; ij < Tt_iajb.size(); ++ij) {
        e_mp2 += K_iajb[ij]->vector_dot(Tt_iajb[ij]);
    }
    return e_mp2;
}


double DLPNOMP2::compute_energy() {

    timer_on("DLPNO");

    print_header();

    int natom = molecule_->natom();
    int nbf = basisset_->nbf();
    int nshell = basisset_->nshell();
    int naux = ribasis_->nbf();
    int nshellri = ribasis_->nshell();

    int nocc = nalpha_;
    int nvir = nbf - nocc;

    int nfocc = basisset_->n_frozen_core(options_.get_str("FREEZE_CORE"), molecule_);
    int naocc = nocc - nfocc;

    outfile->Printf("  Debug Info:\n");
    outfile->Printf("    NATOM = %d\n", natom);
    outfile->Printf("    NFOCC = %d\n", nfocc);
    outfile->Printf("    NAOCC = %d\n", naocc);
    outfile->Printf("    NAVIR = %d\n", nvir);
    outfile->Printf("    NBF   = %d\n", nbf);
    outfile->Printf("    NAUX  = %d\n", naux);
    outfile->Printf("\n");

    std::vector<std::vector<int>> atom_to_bf(natom);
    std::vector<std::vector<int>> shell_to_bf(nshell);
    std::vector<int> bf_to_atom(nbf);
    for(size_t i = 0; i < nbf; ++i) {
        atom_to_bf[basisset_->function_to_center(i)].push_back(i);
        shell_to_bf[basisset_->function_to_shell(i)].push_back(i);
        bf_to_atom[i] = basisset_->function_to_center(i);
    }

    std::vector<std::vector<int>> atom_to_ribf(natom);
    std::vector<std::vector<int>> shell_to_ribf(nshellri);
    std::vector<int> ribf_to_atom(naux);
    for(size_t i = 0; i < naux; ++i) {
        atom_to_ribf[ribasis_->function_to_center(i)].push_back(i);
        shell_to_ribf[ribasis_->function_to_shell(i)].push_back(i);
        ribf_to_atom[i] = ribasis_->function_to_center(i);
    }

    std::vector<std::vector<int>> atom_to_shell(natom);
    for(size_t s = 0; s < nshell; s++) {
        atom_to_shell[basisset_->shell_to_center(s)].push_back(s);
    }

    std::vector<std::vector<int>> atom_to_rishell(natom);
    for(size_t s = 0; s < nshellri; s++) {
        atom_to_rishell[ribasis_->shell_to_center(s)].push_back(s);
    }

    SharedMatrix C_vir = reference_wavefunction_->Ca_subset("AO", "VIR");
    SharedVector e_vir = reference_wavefunction_->epsilon_a_subset("AO", "VIR");

    SharedMatrix C_occ = reference_wavefunction_->Ca_subset("AO", "OCC");
    SharedVector e_occ = reference_wavefunction_->epsilon_a_subset("AO", "OCC");


    // Localize active occupied orbitals
    auto localizer = std::make_shared<BoysLocalizer>(BoysLocalizer(basisset_, reference_wavefunction_->Ca_subset("AO", "ACTIVE_OCC")));
    localizer->set_convergence(options_.get_double("LOCAL_CONVERGENCE"));
    localizer->set_maxiter(options_.get_int("LOCAL_MAXITER"));
    
    localizer->localize();
    SharedMatrix C_lmo = localizer->L();

    SharedMatrix F_lmo = linalg::triplet(C_lmo, reference_wavefunction_->Fa_subset("AO"), C_lmo, true, false, false);
    SharedMatrix S_lmo = linalg::triplet(C_lmo, reference_wavefunction_->S(), C_lmo, true, false, false);
    SharedMatrix H_lmo = linalg::triplet(C_lmo, reference_wavefunction_->H(), C_lmo, true, false, false);
    SharedVector e_lmo = std::make_shared<Vector>(naocc); // LMOs aren't canonical, so not technically orbital energies
    for (size_t i=0; i < naocc; ++i) {
        e_lmo->set(i, F_lmo->get(i, i));
    }
    SharedMatrix SC_lmo = linalg::doublet(reference_wavefunction_->S(), C_lmo, false, false); // intermediate for coefficient fitting

    outfile->Printf("  ==> Forming Local MO Domains <==\n");


    // Form and normalize projected atomic orbitals
    // Note that we project out all (frozen + active) occupied orbitals
    SharedMatrix C_pao = std::make_shared<Matrix>("Projected Atomic Orbitals", nbf, nbf);
    C_pao->identity();
    C_pao->subtract(linalg::triplet(C_occ, C_occ, reference_wavefunction_->S(), false, true, false));
    normalize_cols(C_pao, reference_wavefunction_->S());
    SharedMatrix S_pao = linalg::triplet(C_pao, reference_wavefunction_->S(), C_pao, true, false, false);
    SharedMatrix F_pao = linalg::triplet(C_pao, reference_wavefunction_->Fa(), C_pao, true, false, false);


    timer_on("LMO Domains");

    // Calculate differential overlap integrals for (LMO, LMO) and (LMO, PAO) pairs
    SharedMatrix DOI_ij, DOI_iu;
    std::tie(DOI_ij, DOI_iu) = calculate_dois(reference_wavefunction_, options_, C_lmo, C_pao);
 
    // Calculate approximate dipole energies for (LMO, LMO) pairs
    SharedMatrix e_actual, e_linear, lmo_distances;
    std::tie(e_actual, e_linear, lmo_distances) = dipole_approx(reference_wavefunction_, options_, DOI_iu, C_lmo, e_lmo, C_pao, S_pao, F_pao);


    //                                               //
    // =>    (LMO -> Auxiliary Basis)  Domains    <= //
    //                                               //

    // sparse map from local MOs to auxiliary basis functions in the LMO's domain
    std::vector<std::vector<int>> lmo_to_ribfs(naocc);

    // same  as above, but map to atom indices instead of bf indices
    std::vector<std::vector<int>> lmo_to_riatoms(naocc);

    for (size_t i=0; i < naocc; ++i) {

        // atomic mulliken populations for this orbital
        std::vector<double> mkn_pop(natom, 0.0);

        for(size_t u=0; u < nbf; u++){
            int centerU = basisset_->function_to_center(u);
            double p_uu = reference_wavefunction_->S()->get(u, u) * pow(C_lmo->get(u, i), 2);

            for(size_t v=0; v < nbf; v++){
                int centerV = basisset_->function_to_center(u);
                double p_vv = reference_wavefunction_->S()->get(v, v) * pow(C_lmo->get(v, i), 2);

                // off-diag pops (p_uv) split between u and v prop to diag pops
                double p_uv = reference_wavefunction_->S()->get(u, v) * C_lmo->get(u, i) * C_lmo->get(v, i);
                mkn_pop[centerU] += p_uv * ((p_uu) / (p_uu + p_vv));
                mkn_pop[centerV] += p_uv * ((p_vv) / (p_uu + p_vv));

            }
        }

        // if non-zero mulliken pop on atom, include atom in the MO's fitting domain
        for(size_t a = 0; a < natom; a++) {
            if (fabs(mkn_pop[a]) > options_.get_double("T_CUT_MKN")) {
                lmo_to_riatoms[i].push_back(a);
                // atom's aux orbitals are all-or-nothing
                for(int u : atom_to_ribf[a]) { 
                    lmo_to_ribfs[i].push_back(u);
                }
            }
        }
    }

    print_aux_domains(lmo_to_ribfs, lmo_to_riatoms);


    //                                  //
    // =>    (LMO -> PAO) Domains    <= //
    //                                  //

    // sparse map from LMO [i] to list of PAOs [u] in the domain of [i].
    std::vector<std::vector<int>> lmo_to_paos(naocc);

    // same  as above, but map to atom indices instead of pao indices
    std::vector<std::vector<int>> lmo_to_paoatoms(naocc);

    for(size_t i = 0; i < naocc; ++i) {

        // PAO domains determined by differential overlap integral
        std::vector<int> lmo_to_paos_temp;
        for(size_t u = 0; u < nbf; ++u) {
            if(fabs(DOI_iu->get(i,u)) > options_.get_double("T_CUT_DO")) {
                lmo_to_paos_temp.push_back(u);
            }
        }

        // if any PAO on an atom is in the list, we take all of the PAOs on that atom
        lmo_to_paos[i] = contract_lists(lmo_to_paos_temp, atom_to_bf);

        // equivalent to previous map
        lmo_to_paoatoms[i] = block_list(lmo_to_paos[i], bf_to_atom);

    }

    print_pao_domains(lmo_to_paos, lmo_to_paoatoms);


    //                                  //
    // =>    (LMO -> LMO) Domains    <= //
    //                                  //

    // map from LMO pair index (ij) to both LMO indices (i, j)
    std::vector<std::pair<int,int>> ij_to_i_j;

    // map from LMO indices (i, j) to LMO pair index (ij)
    std::vector<std::vector<int>> i_j_to_ij(naocc);

    int exclude_pairs_overlap = 0;
    int exclude_pairs_energy = 0;
    double de_dipole = 0.0;

    for(size_t i = 0, ij = 0; i < naocc; i++) {
        for(size_t j = 0; j < naocc; j++) {

            bool overlap_big = (DOI_ij->get(i,j) > options_.get_double("T_CUT_DO_ij"));
            if (!overlap_big) exclude_pairs_overlap++;

            bool energy_big = (fabs(e_linear->get(i,j)) > options_.get_double("T_CUT_PRE"));
            if (!energy_big) exclude_pairs_energy++;

            if(overlap_big || energy_big) {
                i_j_to_ij[i].push_back(ij);
                ij_to_i_j.push_back(std::make_pair(i,j));
                ij++;
            } else {
                de_dipole += e_actual->get(i,j); // TODO: scale this number?
                i_j_to_ij[i].push_back(-1);
            }
        }
    }

    int n_lmo_pairs = ij_to_i_j.size();

    // map from LMO pair index (ij) to LMO pair index (ji)
    std::vector<int> ij_to_ji;

    for(size_t ij = 0; ij < n_lmo_pairs; ++ij) {
        size_t i, j;
        std::tie(i,j) = ij_to_i_j[ij];
        ij_to_ji.push_back(i_j_to_ij[j][i]);
    }

    print_lmo_domains(i_j_to_ij, exclude_pairs_overlap, exclude_pairs_energy, de_dipole);


    //                          //
    // =>    Pair Domains    <= //
    //                          //

    timer_on("LMO Pair Domains");

    outfile->Printf("\n  ==> Merging LMO Domains into LMO Pair Domains <==\n");

    // sparse map from LMO pair ij to all PAO's u in the pair domain of i and j.
    std::vector<std::vector<int>> lmopair_to_paos(n_lmo_pairs);
    std::vector<std::vector<int>> lmopair_to_paoatoms(n_lmo_pairs);

    // sparse map from LMO pair ij to all AUXBF u (or atom a) in the pair domain of i and j.
    std::vector<std::vector<int>> lmopair_to_ribfs(n_lmo_pairs);
    std::vector<std::vector<int>> lmopair_to_riatoms(n_lmo_pairs);

#pragma omp parallel for
    for(size_t ij = 0; ij < n_lmo_pairs; ++ij) {
        size_t i, j;
        std::tie(i,j) = ij_to_i_j[ij];

        lmopair_to_paos[ij] = merge_lists(lmo_to_paos[i], lmo_to_paos[j]);
        lmopair_to_paoatoms[ij] = merge_lists(lmo_to_paoatoms[i], lmo_to_paoatoms[j]);

        lmopair_to_ribfs[ij] = merge_lists(lmo_to_ribfs[i], lmo_to_ribfs[j]);
        lmopair_to_riatoms[ij] = merge_lists(lmo_to_riatoms[i], lmo_to_riatoms[j]);
    }

    print_aux_pair_domains(lmopair_to_ribfs, lmopair_to_riatoms);
    print_pao_pair_domains(lmopair_to_paos, lmopair_to_paoatoms);

    timer_off("LMO Pair Domains");










    timer_on("Coefficient Sparsity");

    // which basis functions (on which atoms) contribute to each local MO?
    std::vector<std::vector<int>> lmo_to_bfs(naocc);
    std::vector<std::vector<int>> lmo_to_atoms(naocc);

    for(int i = 0; i < naocc; ++i) {
        for(int bf_ind = 0; bf_ind < nbf; ++bf_ind) {
            if(fabs(C_lmo->get(bf_ind, i)) > options_.get_double("T_CUT_CLMO")) {
                lmo_to_bfs[i].push_back(bf_ind);
            }
        }
        lmo_to_atoms[i] = block_list(lmo_to_bfs[i], bf_to_atom);
        //outfile->Printf("  LMO %d -> BF:  %d / %d bf    %d / %d atoms\n", i, lmo_to_bfs[i].size(), nbf, lmo_to_atoms[i].size(), natom);
    }

    // which basis functions (on which atoms) contribute to each projected AO?
    std::vector<std::vector<int>> pao_to_bfs(nbf);
    std::vector<std::vector<int>> pao_to_atoms(nbf);

    for(int u = 0; u < nbf; ++u) {
        for(int bf_ind = 0; bf_ind < nbf; ++bf_ind) {
            if(fabs(C_pao->get(bf_ind, u)) > options_.get_double("T_CUT_CPAO")) {
                pao_to_bfs[u].push_back(bf_ind);
            }
        }
        pao_to_atoms[u] = block_list(pao_to_bfs[u], bf_to_atom);
        //outfile->Printf("  PAO %d -> BF:  %d / %d bf    %d / %d atoms\n", u, pao_to_bfs[u].size(), nbf, pao_to_atoms[u].size(), natom);
    }

    // extended fitting map (include fitting domains of all local MOs in your domain)
    std::vector<std::vector<int>> lmo_to_riatoms_ext = extend_maps(lmo_to_riatoms, ij_to_i_j);

    // target of the integral transformation
    std::vector<std::vector<int>> riatom_to_lmos_ext = invert_map(lmo_to_riatoms_ext, natom);
    std::vector<std::vector<int>> riatom_to_paos_ext = chain_maps(riatom_to_lmos_ext, lmo_to_paos);
    for(size_t a=0; a < natom; a++) {
        outfile->Printf("%d / %d \n", riatom_to_paos_ext[a].size(), nbf);
    }

    // We'll use these maps to screen the local MO transform: (mn|Q) * C_mi -> (in|Q)
    std::vector<std::vector<int>> riatom_to_atoms1 = chain_maps(riatom_to_lmos_ext, lmo_to_atoms);
    std::vector<std::vector<int>> riatom_to_shells1 = chain_maps(riatom_to_atoms1, atom_to_shell);
    std::vector<std::vector<int>> riatom_to_bfs1 = chain_maps(riatom_to_atoms1, atom_to_bf);

    // We'll use these maps to screen the projected AO transform: (mn|Q) * C_nu -> (mu|Q) 
    std::vector<std::vector<int>> riatom_to_atoms2 = chain_maps(riatom_to_lmos_ext, chain_maps(lmo_to_paos, pao_to_atoms));
    std::vector<std::vector<int>> riatom_to_shells2 = chain_maps(riatom_to_atoms2, atom_to_shell);
    std::vector<std::vector<int>> riatom_to_bfs2 = chain_maps(riatom_to_atoms2, atom_to_bf);

    // Need dense versions of previous maps for quick lookup
    // arr[riatom][lmo] is the index of lmo in riatom_to_lmos_ext[riatom] (if present), else -1
    std::vector<std::vector<int>> riatom_to_lmos_ext_dense(natom, std::vector<int>(naocc, -1));
    std::vector<std::vector<bool>> riatom_to_atoms1_dense(natom, std::vector<bool>(natom, false));
    std::vector<std::vector<bool>> riatom_to_atoms2_dense(natom, std::vector<bool>(natom, false));

    for(int a_ri = 0; a_ri < natom; a_ri++) {
        for(int a_bf : riatom_to_atoms1[a_ri]) {
            riatom_to_atoms1_dense[a_ri][a_bf] = true;
        }
        for(int a_bf : riatom_to_atoms2[a_ri]) {
            riatom_to_atoms2_dense[a_ri][a_bf] = true;
        }
        for(int i_ind = 0; i_ind < riatom_to_lmos_ext[a_ri].size(); i_ind++) {
            int i = riatom_to_lmos_ext[a_ri][i_ind];
            riatom_to_lmos_ext_dense[a_ri][i] = i_ind;
        }
    }

    timer_off("Coefficient Sparsity");

    timer_on("AO ints: (mn|K) -> (ia|K)");
    
    size_t nthread = 1;
#ifdef _OPENMP
    nthread = omp_get_max_threads();
#endif

    std::shared_ptr<IntegralFactory> factory = std::make_shared<IntegralFactory>(ribasis_, BasisSet::zero_ao_basis_set(), basisset_, basisset_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eris(nthread);

    for (size_t thread = 0; thread < nthread; thread++) {
        eris[thread] = std::shared_ptr<TwoBodyAOInt>(factory->eri());
    }

    outfile->Printf("\n  ==> Transforming 3-Index Integrals to LMO/PAO basis <==\n");

    print_integral_sparsity(riatom_to_shells1,
                            riatom_to_shells2, 
                            atom_to_rishell,
                            riatom_to_lmos_ext,
                            riatom_to_paos_ext,
                            atom_to_ribf,
                            nshell,
                            naocc,
                            nbf,
                            naux);


    std::vector<SharedMatrix> qia(naux);

#pragma omp parallel for schedule(static, 1)    
    for (int Q = 0; Q < ribasis_->nshell(); Q++) {
        int nq = ribasis_->shell(Q).nfunction();
        int qstart = ribasis_->shell(Q).function_index();
        int centerQ = ribasis_->shell_to_center(Q);

        size_t thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        // sparse lists of non-screened basis functions
        auto bf_map1 = riatom_to_bfs1[centerQ];
        auto bf_map2 = riatom_to_bfs2[centerQ];

        // inverse map, from global (non-screened) bf-index to Q-specific (screened) index
        std::vector<int> bf_map1_inv(nbf, -1);
        std::vector<int> bf_map2_inv(nbf, -1);
        for(int m_ind = 0; m_ind < bf_map1.size(); m_ind++) {
            bf_map1_inv[bf_map1[m_ind]] = m_ind;
        }
        for(int n_ind = 0; n_ind < bf_map2.size(); n_ind++) {
            bf_map2_inv[bf_map2[n_ind]] = n_ind;
        }
        
        for(size_t q = 0; q < nq; q++) {
            qia[qstart+q] = std::make_shared<Matrix>("(mn|Q)", bf_map1.size(), bf_map2.size());
        }

        for (int M : riatom_to_shells1[centerQ]) {
            int nm = basisset_->shell(M).nfunction();
            int mstart = basisset_->shell(M).function_index();
            int centerM = basisset_->shell_to_center(M);

            for (int N : riatom_to_shells2[centerQ]) {
                int nn = basisset_->shell(N).nfunction();
                int nstart = basisset_->shell(N).function_index();
                int centerN = basisset_->shell_to_center(N);

                // is (N in the list of M's) and (M in the list of N's)?
                bool MN_symmetry = (riatom_to_atoms1_dense[centerQ][centerN] && riatom_to_atoms2_dense[centerQ][centerM]);

                // if so, we want to exploit (MN|Q) <-> (NM|Q) symmetry
                if(N < M && MN_symmetry) continue;

                eris[thread]->compute_shell(Q, 0, M, N);
                const double* buffer = eris[thread]->buffer();

                for (int q = 0, index = 0; q < nq; q++) {
                    for (int m = 0; m < nm; m++) {
                        for (int n = 0; n < nn; n++, index++) {
                            qia[qstart+q]->set(bf_map1_inv[mstart+m], bf_map2_inv[nstart+n], buffer[index]);
                        }
                    }
                }

                // (MN|Q) <-> (NM|Q) symmetry
                if(N > M && MN_symmetry) {
                    for (int q = 0, index = 0; q < nq; q++) {
                        for (int m = 0; m < nm; m++) {
                            for (int n = 0; n < nn; n++, index++) {
                                qia[qstart+q]->set(bf_map1_inv[nstart+n], bf_map2_inv[mstart+m], buffer[index]);
                            }
                        }
                    }
                }

            } // N loop
        } // M loop

        SharedMatrix C_lmo_slice = get_rows_and_cols(C_lmo, riatom_to_bfs1[centerQ], riatom_to_lmos_ext[centerQ]);
        SharedMatrix C_pao_slice = get_rows(C_pao, riatom_to_bfs2[centerQ]); // TODO: PAO slices

        //// Here we'll refit the coefficients of C_lmo_slice to minimize residual from unscreened orbitals
        //// This lets us get away with agressive coefficient screening
        //// Boughton and Pulay 1992 JCC, Equation 3

        // Solve for C_lmo_slice such that S[local,local] @ C_lmo_slice ~= S[local,all] @ C_lmo
        SharedMatrix SC_lmo_slice = get_rows_and_cols(SC_lmo, riatom_to_bfs1[centerQ], riatom_to_lmos_ext[centerQ]);
        SharedMatrix S_aa = get_rows_and_cols(reference_wavefunction_->S(), riatom_to_bfs1[centerQ], riatom_to_bfs1[centerQ]);

        int nbf_solve = riatom_to_bfs1[centerQ].size();
        int nmo_solve = riatom_to_lmos_ext[centerQ].size();

        // DGESV expects data in fortran ordering, so we have to be more explicit with memory
        // TODO: Can I use this function w/o doing this stupid copy?
        double* C_lmo_slice_fit = new double[nbf_solve * nmo_solve];
        for(int ind1 = 0; ind1 < nbf_solve; ind1++) {
            for(int ind2 = 0; ind2 < nmo_solve; ind2++) {
                C_lmo_slice_fit[ind2 * nbf_solve + ind1] = SC_lmo_slice->get(ind1, ind2);
            }
        }

        if(nmo_solve > 0) {
            int *ipiv = init_int_array(nbf_solve);
            int errcode = C_DGESV(nbf_solve, nmo_solve, S_aa->pointer()[0], nbf_solve, ipiv, C_lmo_slice_fit, nbf_solve);
            if(errcode != 0) {
                outfile->Printf("C_DGESV Error: %d\n", errcode);
            }
            free(ipiv);
        }

        // reverse the fortran ordering
        for(int ind1 = 0; ind1 < nbf_solve; ind1++) {
            for(int ind2 = 0; ind2 < nmo_solve; ind2++) {
                C_lmo_slice->set(ind1, ind2, C_lmo_slice_fit[ind2 * nbf_solve + ind1]);
            }
        }

        delete[] C_lmo_slice_fit;

        // (mn|Q) C_mi C_nu -> (iu|Q)
        for(size_t q = 0; q < nq; q++) {
            qia[qstart+q] = linalg::triplet(C_lmo_slice, qia[qstart+q], C_pao_slice, true, false, false);
        }

    } // Q loop

    timer_off("AO ints: (mn|K) -> (ia|K)");

    timer_on("Full Metric: (K|L)");

    // Compute the full metric, don't invert
    std::shared_ptr<FittingMetric> metric = std::make_shared<FittingMetric>(ribasis_, true);
    metric->form_fitting_metric();
    SharedMatrix full_metric = std::make_shared<Matrix>(metric->get_metric());

    timer_off("Full Metric: (K|L)");

    timer_on("PNO Transformation");

    outfile->Printf("\n  ==> Forming Pair Natural Orbitals <==\n");

    std::vector<SharedMatrix> K_iajb(n_lmo_pairs);  // exchange operators (i.e. (ia|jb) integrals)
    std::vector<SharedMatrix> T_iajb(n_lmo_pairs);  // amplitudes
    std::vector<SharedMatrix> Tt_iajb(n_lmo_pairs); // antisymmetrized amplitudes
    std::vector<SharedMatrix> X_pno(n_lmo_pairs);   // global PAOs -> canonical PNOs
    std::vector<SharedVector> e_pno(n_lmo_pairs);   // PNO orbital energies

    std::vector<int> n_pno(n_lmo_pairs, 0);         // number of pnos
    std::vector<double> de_pno(n_lmo_pairs, 0.0);   // PNO truncation error

#pragma omp parallel for schedule(static, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {

        int i, j;
        std::tie(i, j) = ij_to_i_j[ij];
        int ji = ij_to_ji[ij];

        if(i > j) continue;

        //                                       //
        // ==> Assemble (ia|jb) in PAO basis <== //
        //                                       //

        // number of PAOs in the pair domain (before removing linear dependencies)
        int npao_ij = lmopair_to_paos[ij].size();//X_pao_ij->rowspi(0);

        // number of auxiliary basis in the domain
        int naux_ij = lmopair_to_ribfs[ij].size();

        SharedMatrix i_qa = std::make_shared<Matrix>("blah", naux_ij, npao_ij);
        SharedMatrix j_qa = std::make_shared<Matrix>("blah", naux_ij, npao_ij);

        for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
            int q = lmopair_to_ribfs[ij][q_ij];
            int centerq = ribasis_->function_to_center(q);
            for (int a_ij = 0; a_ij < npao_ij; a_ij++) {
                int a = lmopair_to_paos[ij][a_ij];
                i_qa->set(q_ij, a_ij, qia[q]->get(riatom_to_lmos_ext_dense[centerq][i], a));  
                j_qa->set(q_ij, a_ij, qia[q]->get(riatom_to_lmos_ext_dense[centerq][j], a));  
            }
        }

        int N_solve = naux_ij;
        int M_solve = npao_ij;

        SharedMatrix A_solve = get_rows_and_cols(full_metric, lmopair_to_ribfs[ij], lmopair_to_ribfs[ij]);

        // DGESV expects data in fortran ordering, so we have to be more explicit with memory
        // TODO: Can I use this function w/o doing this stupid copy?
        double* B_solve = new double[N_solve * M_solve];
        for(int indn = 0; indn < N_solve; indn++) {
            for(int indm = 0; indm < M_solve; indm++) {
                B_solve[indm * N_solve + indn] = i_qa->get(indn, indm);
            }
        }

        int *ipiv = init_int_array(N_solve);
        int errcode = C_DGESV(N_solve, M_solve, A_solve->pointer()[0], N_solve, ipiv, B_solve, N_solve);
        if(errcode != 0) {
            outfile->Printf("C_DGESV Error: %d\n", errcode);
        }
        free(ipiv);

        // reverse the fortran ordering
        for(int indn = 0; indn < N_solve; indn++) {
            for(int indm = 0; indm < M_solve; indm++) {
                i_qa->set(indn, indm, B_solve[indm * N_solve + indn]);
            }
        }

        delete[] B_solve;
        SharedMatrix K_pao_ij = linalg::doublet(i_qa, j_qa, true, false);

        //                                       //
        // ==> Canonicalize Pair Domain PAOs <== //
        //                                       //

        SharedMatrix S_pao_ij = get_rows_and_cols(S_pao, lmopair_to_paos[ij], lmopair_to_paos[ij]);
        SharedMatrix F_pao_ij = get_rows_and_cols(F_pao, lmopair_to_paos[ij], lmopair_to_paos[ij]);

        SharedMatrix X_pao_ij; // canonical transformation of this domain's PAOs to
        SharedVector e_pao_ij; // energies of the canonical PAOs
        std::tie(X_pao_ij, e_pao_ij) = get_orthocanonicalizer(S_pao_ij, F_pao_ij, options_);

        //S_pao_ij = linalg::triplet(X_pao_ij, S_pao_ij, X_pao_ij, true, false, false);
        F_pao_ij = linalg::triplet(X_pao_ij, F_pao_ij, X_pao_ij, true, false, false);
        K_pao_ij = linalg::triplet(X_pao_ij, K_pao_ij, X_pao_ij, true, false, false);

        // number of PAOs in the domain after removing linear dependencies
        int npao_can_ij = X_pao_ij->colspi(0);
        SharedMatrix T_pao_ij = K_pao_ij->clone();
        for (int a = 0; a < npao_can_ij; ++a) {
            for (int b = 0; b < npao_can_ij; ++b) {
                T_pao_ij->set(a, b, T_pao_ij->get(a,b) / (-e_pao_ij->get(b) + -e_pao_ij->get(a) + e_lmo->get(i) + e_lmo->get(j)));
            }
        }

        //                                           //
        // ==> Canonical PAOs  to Canonical PNOs <== //
        //                                           //

        SharedVector e_pno_ij;
        SharedMatrix K_pno_ij, T_pno_ij, X_pno_ij;
        double de_pno_ij;
        std::tie(K_pno_ij, T_pno_ij, X_pno_ij, e_pno_ij, de_pno_ij) = pno(F_pao_ij, K_pao_ij, e_pao_ij, e_lmo->get(i), e_lmo->get(j), options_);
        
        X_pno_ij = linalg::doublet(X_pao_ij, X_pno_ij, false, false);

        K_iajb[ij] = K_pno_ij;
        T_iajb[ij] = T_pno_ij;
        X_pno[ij] = X_pno_ij;
        e_pno[ij] = e_pno_ij;
        n_pno[ij] = X_pno_ij->colspi(0);
        de_pno[ij] = de_pno_ij;

        // account for symmetry
        if (i < j) {
            K_iajb[ji] = K_iajb[ij]->transpose();
            T_iajb[ji] = T_iajb[ij]->transpose();
            X_pno[ji] = X_pno[ij];
            e_pno[ji] = e_pno[ij];
            n_pno[ji] = n_pno[ij];
            de_pno[ji] = de_pno_ij;
        }

    }

    int pno_count_total = 0, pno_count_min = nbf, pno_count_max = 0;
    double de_pno_total = 0.0;
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        pno_count_total += n_pno[ij];
        pno_count_min = std::min(pno_count_min, n_pno[ij]);
        pno_count_max = std::max(pno_count_max, n_pno[ij]);
        de_pno_total += de_pno[ij];
    }

    outfile->Printf("  \n");
    outfile->Printf("    Natural Orbitals per Local MO pair:\n");
    outfile->Printf("      Avg: %3d NOs \n", pno_count_total / n_lmo_pairs);
    outfile->Printf("      Min: %3d NOs \n", pno_count_min);
    outfile->Printf("      Max: %3d NOs \n", pno_count_max);
    outfile->Printf("  \n");
    outfile->Printf("    PNO truncation energy = %.12f\n", de_pno_total);

#pragma omp parallel for schedule(static, 1)    
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        Tt_iajb[ij] = T_iajb[ij]->clone();
        Tt_iajb[ij]->scale(2.0);
        Tt_iajb[ij]->subtract(T_iajb[ij]->transpose());
    }

    timer_off("PNO Transformation");

    timer_on("LMP2 Iterations");

    //outfile->Printf("\n  ==> Precalculating PNO Overlaps <==\n");

    timer_on("PNO Overlaps");
    std::vector<std::vector<SharedMatrix>> pno_overlaps1(n_lmo_pairs, std::vector<SharedMatrix>(naocc));
    std::vector<std::vector<SharedMatrix>> pno_overlaps2(n_lmo_pairs, std::vector<SharedMatrix>(naocc));

#pragma omp parallel for schedule(static, 1)    
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {

        int i, j;
        std::tie(i, j) = ij_to_i_j[ij];
        int ji = ij_to_ji[ij];

        if(n_pno[ij] == 0) continue;
        if(i < j) continue;
        
        for (int k = 0; k < naocc; ++k) {

            int kj = i_j_to_ij[k][j];
            if (kj != -1 && i != k && fabs(F_lmo->get(i,k)) > options_.get_double("F_CUT") && n_pno[kj] > 0) {
                pno_overlaps1[ij][k] = get_rows_and_cols(S_pao, lmopair_to_paos[ij], lmopair_to_paos[kj]);
                pno_overlaps1[ij][k] = linalg::triplet(X_pno[ij], pno_overlaps1[ij][k], X_pno[kj], true, false, false);

            }

            int ik = i_j_to_ij[i][k];
            if (ik != -1 && j != k && fabs(F_lmo->get(k,j)) > options_.get_double("F_CUT") && n_pno[ik] > 0) {
                pno_overlaps2[ij][k] = get_rows_and_cols(S_pao, lmopair_to_paos[ij], lmopair_to_paos[ik]);
                pno_overlaps2[ij][k] = linalg::triplet(X_pno[ij], pno_overlaps2[ij][k], X_pno[ik], true, false, false);

            }

            if(i > j) {
                pno_overlaps1[ji][k] = pno_overlaps2[ij][k];
                pno_overlaps2[ji][k] = pno_overlaps1[ij][k];
            }

        }

    }
    timer_off("PNO Overlaps");

    outfile->Printf("\n  ==> Local MP2 <==\n\n");
    outfile->Printf("    E_CONVERGENCE = %.2e\n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE = %.2e\n", options_.get_double("R_CONVERGENCE"));

    // computed every iteration
    std::vector<double> iter_energies;
    std::vector<double> iter_delta_energies;
    std::vector<double> iter_max_residuals;

    iter_energies.push_back(mp2_energy(Tt_iajb, K_iajb));
    std::vector<double> pair_energies_mp2;
    iter_delta_energies.push_back(iter_energies[0]);
    iter_max_residuals.push_back(1000.0);

    outfile->Printf("\n");
    outfile->Printf("                     Corr. Energy    Delta E     Max R");
    outfile->Printf("\n");
    outfile->Printf("  @LMP2 iter  SC: %16.12f ---------- ----------\n", iter_energies[0]);

    int iteration = 0;
    while(fabs(iter_delta_energies[iteration]) > options_.get_double("E_CONVERGENCE") || iter_max_residuals[iteration] > options_.get_double("R_CONVERGENCE")) {

        std::vector<SharedMatrix> R_iajb(n_lmo_pairs);
        std::vector<double> R_iajb_rms(n_lmo_pairs, 0.0);

        timer_on("Main Loop");
#pragma omp parallel for schedule(static, 1)    
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {

            int i, j;
            std::tie(i, j) = ij_to_i_j[ij];

            R_iajb[ij] = std::make_shared<Matrix>("Residual", n_pno[ij], n_pno[ij]);

            if(n_pno[ij] == 0) continue;

            for (int a = 0; a < n_pno[ij]; ++a) {
                for (int b = 0; b < n_pno[ij]; ++b) {
                    R_iajb[ij]->set(a, b, K_iajb[ij]->get(a,b) + (e_pno[ij]->get(a) + e_pno[ij]->get(b) - e_lmo->get(i) - e_lmo->get(j)) * T_iajb[ij]->get(a,b));
                }
            }

            for (int k = 0; k < naocc; ++k) {
                int kj = i_j_to_ij[k][j];
                int ik = i_j_to_ij[i][k];

                if (kj != -1 && i != k && fabs(F_lmo->get(i,k)) > options_.get_double("F_CUT") && n_pno[kj] > 0) {
                    SharedMatrix temp = linalg::triplet(pno_overlaps1[ij][k], T_iajb[kj], pno_overlaps1[ij][k], false, false, true);
                    temp->scale(-1.0 * F_lmo->get(i,k));
                    R_iajb[ij]->add(temp);
                }
                if (ik != -1 && j != k && fabs(F_lmo->get(k,j)) > options_.get_double("F_CUT") && n_pno[ik] > 0) {
                    SharedMatrix temp = linalg::triplet(pno_overlaps2[ij][k], T_iajb[ik], pno_overlaps2[ij][k], false, false, true);
                    temp->scale(-1.0 * F_lmo->get(k,j));
                    R_iajb[ij]->add(temp);
                }
            }

            R_iajb_rms[ij] = R_iajb[ij]->rms();

        }

#pragma omp parallel for schedule(static, 1)    
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {

            int i, j;
            std::tie(i, j) = ij_to_i_j[ij];
            for (int a = 0; a < n_pno[ij]; ++a) {
                for (int b = 0; b < n_pno[ij]; ++b) {
                    T_iajb[ij]->add(a, b, -R_iajb[ij]->get(a,b) / ((e_pno[ij]->get(a) + e_pno[ij]->get(b)) - (e_lmo->get(i) + e_lmo->get(j))));
                }
            }

            Tt_iajb[ij]->zero();
            Tt_iajb[ij]->add(T_iajb[ij]);
            Tt_iajb[ij]->scale(2.0);
            Tt_iajb[ij]->subtract(T_iajb[ij]->transpose());
        }

        timer_off("Main Loop");

        iteration++;

        // the max residual over all ij pairs
        double max_R_rms = *max_element(R_iajb_rms.begin(), R_iajb_rms.end());
        double e_curr = mp2_energy(Tt_iajb, K_iajb);
        double e_prev = iter_energies[iter_energies.size() - 1];

        iter_max_residuals.push_back(max_R_rms);
        iter_energies.push_back(e_curr);
        iter_delta_energies.push_back(e_curr - e_prev);

        outfile->Printf("  @LMP2 iter %3d: %16.12f %10.3e %10.3e\n", iteration, e_curr, e_curr - e_prev, max_R_rms);
    }

    timer_off("LMP2 Iterations");

    outfile->Printf("  \n"); 
    outfile->Printf("  Total DLPNO-MP2 Correlation Energy: %16.12f \n", iter_energies[iter_energies.size() - 1] + de_pno_total + de_dipole);
    outfile->Printf("    MP2 Correlation Energy:           %16.12f \n", iter_energies[iter_energies.size() - 1]);
    outfile->Printf("    LMO Truncation Correction:        %16.12f \n", de_dipole);
    outfile->Printf("    PNO Truncation Correction:        %16.12f \n", de_pno_total);

    // TODO: put vars in the right place
    //ref_wfn->set_scalar_variable("DLPNO-MP2 CORRELATION ENERGY", iter_energies[iter_energies.size() - 1] + de_pno_total + de_dipole);

    //ref_wfn->set_array_variable("LMO FOCK MATRIX", F_lmo);
    //ref_wfn->set_array_variable("LMO OVERLAP MATRIX", S_lmo);
    //ref_wfn->set_array_variable("LMO CORE MATRIX", H_lmo);
    //ref_wfn->set_array_variable("LMO DISTANCE MATRIX", lmo_distances);

    //ref_wfn->set_array_variable("MP2 PAIR ENERGIES", pair_energies_mat_mp2);
    //ref_wfn->set_array_variable("MP2_SC PAIR ENERGIES", pair_energies_mat_sc);
    //ref_wfn->set_array_variable("MP2_DIPOLE PAIR ENERGIES", e_actual);

    timer_off("DLPNO");

    return 0.0;
}

void DLPNOMP2::print_header() {

    outfile->Printf("   --------------------------------------------\n");
    outfile->Printf("                     DLPNO-MP2                 \n");
    outfile->Printf("                   by Zach Glick               \n");
    outfile->Printf("   --------------------------------------------\n");
    outfile->Printf("\n");

    outfile->Printf("    T_CUT_DO_ij  = %6.3e \n", options_.get_double("T_CUT_DO_ij"));
    outfile->Printf("    T_CUT_PRE    = %6.3e \n", options_.get_double("T_CUT_PRE"));
    outfile->Printf("    T_CUT_DO_PRE = %6.3e \n", options_.get_double("T_CUT_DO_PRE"));
    outfile->Printf("    T_CUT_DO     = %6.3e \n", options_.get_double("T_CUT_DO"));
    outfile->Printf("    T_CUT_MKN    = %6.3e \n", options_.get_double("T_CUT_MKN"));
    outfile->Printf("    T_CUT_CLMO   = %6.3e \n", options_.get_double("T_CUT_CLMO"));
    outfile->Printf("    T_CUT_CPAO   = %6.3e \n", options_.get_double("T_CUT_CPAO"));
    outfile->Printf("    S_CUT        = %6.3e \n", options_.get_double("S_CUT"));
    outfile->Printf("    T_CUT_PNO    = %6.3e \n", options_.get_double("T_CUT_PNO"));
    outfile->Printf("    F_CUT        = %6.3e \n", options_.get_double("F_CUT"));
    outfile->Printf("\n");

}

}  // namespace dlpnomp2
}  // namespace psi
