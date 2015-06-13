/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */


/**!
This is a generalized second-order SCF module capable of computing the
optimal orbital rotations at each macroiteration for RHF, RASSCF, and CASSCF.

Derived from the purple book
Chapters 10.8.8,

Indices:
mu, nu, sigma, ro - General AO indices
mnopqrs           - General MO indices
ijkl              - Occupied indices
tuvw              - Active indices
abcd              - Virtual indices

 - DGA Smith, 22 May, 2015
**/

#ifndef SOSCF_H
#define SOSCF_H

#include <libmints/typedefs.h>
#include <map>
#include <boost/tuple/tuple.hpp>

namespace psi {

/// Forward declare
class JK;
class DFERI;

class SORHF {

public:
    /**
     * Initialize the SORHF object
     * @param jk object to use
     */
    SORHF(boost::shared_ptr<JK> jk);

    ~SORHF(void);

    /**
     * Rotate the current orbitals for a given rotation matrix.
     * @param  x The [o, v] non-redundant orbital rotation parameters.
     * @return   The rotated orbitals.
     */
    SharedMatrix Ck(SharedMatrix x);

    /**
     * Update to the current macroiteration.
     * @param Cocc Occupied orbitals
     * @param Cvir Virtual orbitals
     * @param Fock SO Fock Matrix
     */
    void update(SharedMatrix Cocc, SharedMatrix Cvir, SharedMatrix Fock);

    /**
     * Returns the hessian times a trial vector or, in effect, the rotated inactive Fock matrix.
     * IF_k = IF_mp K_np + IF_pn K_mp + (4 G_mnip - G_mpin - G_npim) K_ip
     * @param  x  The [o, v] matrix of non-redundant orbital rotations.
     * @return Hx The [o, v] block of the rotated Fock matrix.
     */
    SharedMatrix Hk(SharedMatrix x);

    /**
     * Solves the set of linear equations Hx = gradient using CG. In this particular case only
     * orbital-orbital contributions are considered.
     * @return x The [o, v] matrix of non-redundant orbital rotation parameters.
     */
    SharedMatrix solve(int max_iter=5, double conv=1.e-10, bool print=true);

protected:

    /// Parameters
    size_t nirrep_;
    size_t nocc_;
    Dimension noccpi_;
    size_t nvir_;
    Dimension nvirpi_;
    size_t nso_;
    Dimension nsopi_;

    /// Global JK object
    boost::shared_ptr<JK> jk_;

    /// Map of matrices
    std::map<std::string, SharedMatrix > matrices_;

}; // SORHF class

class SOMCSCF {

public:

    /**
     * Initialize the SOMCSCF object
     * @param jk      JK object to use.
     * @param H       Core hamiltonian in the SO basis.
     * @param casscf  Is this a CAS calculation? (ignore active-active rotations)
     */
    // SOMCSCF(boost::shared_ptr<JK> jk, SharedMatrix H, bool casscf);
    SOMCSCF(boost::shared_ptr<JK> jk, boost::shared_ptr<DFERI> df, SharedMatrix AOTOSO,
            SharedMatrix H, bool casscf);

    ~SOMCSCF(void);

    /**
     * DGAS -- I dont think we actually need this.
     * Sets the frozen core. Here frozen core are orbitals that do not rotate.
     * @param Cfzc The frozen core orbitals
     */
    void set_frozen_core(SharedMatrix Cfzc);

    /**
     * DGAS -- I dont think we actually need this.
     * Sets the virtual core. Here virtual core are orbitals that do not rotate.
     * @param Cfzv The frozen core orbitals
     */
    void set_frozen_virtual(SharedMatrix Cfzv);

    /**
     * Updates all internal variables to the new reference frame and builds
     * non-rotated intermediates.
     * @param Cocc Current doubly occupied orbitals
     * @param Cocc Current active orbitals
     * @param Cocc Current virtual orbitals
     * @param OPDM Current active one-particle density matrix
     * @param TPDM Current active two-particle density matrix (symmetrized, dense)
    */
    void update(SharedMatrix Cocc, SharedMatrix Cact, SharedMatrix Cvir,
                SharedMatrix OPDM, SharedMatrix TPDM);
    /**
     * Returns the approximate diagonal hessian.
     * @return Hdiag The [oa, av] block of the diagonal hessian.
     */
    SharedMatrix H_approx_diag();

    /**
     * Returns the hessian times a trial vector or, in effect, the rotated inactive Fock matrix.
     * IF_k = IF_mp K_np + IF_pn K_mp + (4 G_mnip - G_mpin - G_npim) K_ip
     * @param  x  The [oa, av] matrix of non-redundant orbital rotations.
     * @return Hx The [oa, av] block of the rotated Fock matrix.
     */
    SharedMatrix Hk(SharedMatrix x);

    /**
     * Solves the set of linear equations Hx = gradient using CG.
     * @return x The [oa, av] matrix of non-redundant orbital rotation parameters.
     */
    SharedMatrix solve(int max_iter=5, double conv=1.e-10, bool print=true);


protected:

    /// Parameters
    bool casscf_;

    // Frozen core orbitals in the sense that they cannot rotate
    size_t nfzc_;
    Dimension nfzcpi_;
    bool freeze_core_;

    size_t nocc_;
    Dimension noccpi_;
    size_t nact_;
    Dimension nactpi_;
    size_t nvir_;
    Dimension nvirpi_;

    // Frozen virtual orbitals in the sense that they cannot rotate
    size_t nfzv_;
    Dimension nfzvpi_;
    bool freeze_virtual_;

    // General info
    size_t nirrep_;
    size_t nmo_;
    Dimension nmopi_;
    size_t nso_;
    size_t nao_;
    Dimension nsopi_;

    // Non-redunant rotatations
    Dimension noapi_;
    Dimension navpi_;

    // AOTOSO
    SharedMatrix AOTOSO_;

    /// Global JK object
    boost::shared_ptr<JK> jk_;

    /// Global libtrans object
    /// boost::shared_ptr<IntegralTransform> ints_;

    /// Global DFERI object
    boost::shared_ptr<DFERI> dferi_;

    /// Map of matrices
    std::map<std::string, SharedMatrix > matrices_;

    /// Zeros out the active-active part of a trial vector
    void zero_act(SharedMatrix vector);


}; // SOMCSCF class

} // Namespace psi


#endif
