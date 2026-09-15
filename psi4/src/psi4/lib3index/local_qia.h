/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#ifndef three_index_local_qia_h_
#define three_index_local_qia_h_

#include "psi4/pragma.h"

#include <memory>
#include <utility>
#include <vector>

namespace psi {

class BasisSet;
class Matrix;
typedef std::shared_ptr<Matrix> SharedMatrix;

/// Three-index integrals transformed into per-atom orbital domains, built only
/// where a local correlation method will read them.
///
/// The quantity is the same one a density-fitted method half-transforms out of
/// the AO integrals,
///
///     (Q|iu) = sum_{mn} C_lmo(m,i) (mn|Q) C_pao(n,u),
///
/// with no fitting metric applied, because a local method fits per domain
/// rather than globally. What differs from building it through DFHelper or
/// MintsHelper is that the caller says up front which (i, u) pairs it will read
/// against the auxiliary functions of each atom, and nothing outside that
/// demand is ever computed. For a local method the demand is a vanishing
/// fraction of the full tensor at large system size, and the AO integrals it
/// implies are a vanishing fraction of the screened AO integrals, so this is
/// where a linear-scaling method's three-index phase wants to live.
///
/// Two screenings are in force and they are not the same kind of thing:
///
///   - the shell-pair tolerance (set_ints_tolerance) drops AO integrals that
///     are numerically negligible. It is controlled: the error goes to zero
///     with the tolerance, and at zero the result is the exact transform.
///   - the domain restriction drops contributions from basis functions outside
///     the extended domain of the requested orbitals. That is an approximation
///     even at exact arithmetic, mitigated but not removed by the
///     Boughton-Pulay refit. Its size is controlled by the coefficient
///     tolerances.
///
/// With both coefficient tolerances negative, the integral tolerance zero and
/// the refit off, the domains become the full basis and the result is the exact
/// transform of the full AO integrals - which is how a caller can check its
/// indexing against a dense reference.
///
/// Usage mirrors DFHelper: construct, set, declare the demand, compute.
///
///     LocalQiaBuilder builder(primary, aux);
///     builder.set_nthreads(nthread);
///     builder.set_domains(atom_to_lmos, atom_to_paos);
///     auto blocks = builder.compute(C_lmo, C_pao, S);
///
/// This is a standalone builder rather than a method on DFHelper because
/// DFHelper's subject is the fitted B tensors of a full screened AO build, and
/// this never forms that build: the two share an integral factory and nothing
/// else. It is deliberately not coupled to psi4's own DLPNO code either, though
/// the algorithm is that of DLPNO::compute_qia, so that a caller supplying its
/// own domains does not inherit DLPNO's option handling or its member state.
class PSI_API LocalQiaBuilder {
   public:
    LocalQiaBuilder(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> aux);

    /// Threads for the auxiliary-shell loop. Defaults to the process thread
    /// count, so a caller that forgets this still gets the machine it asked
    /// psi4 for.
    void set_nthreads(int nthread);
    int nthreads() const { return nthread_; }

    /// Shell-pair screening on the AO integrals, in the same units as DLPNO's
    /// DLPNO_AO_INTS_TOL: a shell pair is skipped when
    /// max_Q (Q|Q) * (MN|MN)^1/2 < tol^2. Zero or negative disables it.
    void set_ints_tolerance(double tol);
    double ints_tolerance() const { return ints_tolerance_; }

    /// A basis function enters an atom's LMO-side domain when some requested
    /// LMO has |C_lmo| above this on any function of the same atom. DLPNO's
    /// T_CUT_CLMO. Negative admits every basis function.
    void set_coef_tolerance_lmo(double tol);
    double coef_tolerance_lmo() const { return coef_tol_lmo_; }

    /// The same for the PAO-side domain. DLPNO's T_CUT_CPAO.
    void set_coef_tolerance_pao(double tol);
    double coef_tolerance_pao() const { return coef_tol_pao_; }

    /// Refit the LMO coefficients onto each restricted basis-function domain
    /// before transforming, Boughton and Pulay 1992 JCC eq. 3. On by default:
    /// it is what lets the coefficient screening be aggressive without the
    /// dropped tail costing accuracy. Turning it off makes the transform a
    /// plain slice of C_lmo, which is what an exactness check wants.
    void set_bp_refit(bool on);
    bool bp_refit() const { return bp_refit_; }

    /// The demand: for each ATOM of the auxiliary basis, which LMO columns and
    /// which PAO columns will be read against the auxiliary functions on that
    /// atom. Atom indices follow aux->molecule(), and both vectors must have
    /// one entry per atom; an atom given empty lists is skipped entirely and
    /// answers with a null block.
    ///
    /// The ordering of each list is the caller's and is preserved in the
    /// output, so a caller can scatter with it directly.
    void set_domains(const std::vector<std::vector<int>>& atom_to_lmos,
                     const std::vector<std::vector<int>>& atom_to_paos);

    /// Build the declared blocks.
    ///
    /// C_lmo is (nbf, nlmo) and C_pao is (nbf, npao) in the primary basis. S is
    /// the AO overlap, read only when the Boughton-Pulay refit is on; it may be
    /// null otherwise.
    ///
    /// Returns one matrix per auxiliary atom, of shape
    /// (naux on that atom, |lmos| * |paos|), the second axis row-major in
    /// (lmo, pao) with the caller's own list ordering. Atoms with an empty
    /// domain, and atoms carrying no auxiliary functions, answer with nullptr.
    std::vector<SharedMatrix> compute(SharedMatrix C_lmo, SharedMatrix C_pao, SharedMatrix S);

    /// The [first, last) auxiliary function index range of each atom, so a
    /// caller can scatter each block into a full-size tensor. Empty atoms give
    /// an empty range.
    std::vector<std::pair<int, int>> atom_aux_ranges() const { return atom_aux_range_; }

    /// How many AO shell pairs the last compute() actually evaluated, against
    /// how many the domains alone would have asked for. The honest measure of
    /// what the screening bought, and zero-zero before the first compute().
    std::pair<size_t, size_t> last_shell_pair_counts() const {
        return {shell_pairs_computed_, shell_pairs_offered_};
    }

   private:
    /// Max of the auxiliary metric diagonal over each auxiliary shell, the
    /// (Q|Q) factor of the screening estimate. Built on first use and only when
    /// screening is on.
    void build_metric_shell_diag();

    /// Per-atom basis-function and shell domains implied by the declared
    /// orbital domains and the coefficient tolerances.
    void build_bf_domains(SharedMatrix C_lmo, SharedMatrix C_pao);

    std::shared_ptr<BasisSet> primary_;
    std::shared_ptr<BasisSet> aux_;

    int natom_;
    int nbf_;
    int naux_;

    int nthread_;
    double ints_tolerance_;
    double coef_tol_lmo_;
    double coef_tol_pao_;
    bool bp_refit_;

    bool domains_set_;
    std::vector<std::vector<int>> atom_to_lmos_;
    std::vector<std::vector<int>> atom_to_paos_;

    /// atom_to_bf_[A] / atom_to_shell_[A]: the primary basis functions and
    /// shells centered on atom A. Static, built in the constructor.
    std::vector<std::vector<int>> atom_to_bf_;
    std::vector<std::vector<int>> atom_to_shell_;
    std::vector<std::pair<int, int>> atom_aux_range_;

    /// Per auxiliary atom, the LMO-side (1) and PAO-side (2) domains: which
    /// primary atoms contribute, and the shells and basis functions those
    /// imply. Rebuilt by each compute() because they depend on the
    /// coefficients.
    std::vector<std::vector<char>> atom_in_domain1_;
    std::vector<std::vector<char>> atom_in_domain2_;
    std::vector<std::vector<int>> domain_shells1_;
    std::vector<std::vector<int>> domain_shells2_;
    std::vector<std::vector<int>> domain_bfs1_;
    std::vector<std::vector<int>> domain_bfs2_;

    std::vector<double> metric_shell_diag_;

    size_t shell_pairs_computed_;
    size_t shell_pairs_offered_;
};

}  // namespace psi

#endif
