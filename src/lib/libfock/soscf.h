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

class JK;
class DFTensor;

class SORHF {

public:
    /**
     * Initialize the SORHF object
     * @param jk object to compute
     */
    SORHF(boost::shared_ptr<JK> jk);

    ~SORHF(void);

    /**
     * Update to the current macroiteration
     * @param Cocc Occupied orbitals
     * @param Cvir Virtual orbitals
     * @param Fock SO Fock Matrix
     */
    void update(SharedMatrix Cocc, SharedMatrix Cvir, SharedMatrix Fock);

    /**
     * Returns the hessian times a trial vector or, in effect, the rotated inactive Fock matrix
     * IF_k = IF_mp K_np + IF_pn K_mp + (4 G_mnip - G_mpin - G_npim) K_ip
     * @param  x  The [o, v] matrix of non-redundant orbital rotations.
     * @return Hx The [o, v] block of the rotated Fock matrix
     */
    SharedMatrix Hx(SharedMatrix x);

    /**
     * Solves the set of linear equations Hx = gradient using CG. In this particular case only
     * orbital-orbital contributions are considered.
     * @return x The [o, v] matrix of non-redundant orbital rotation parameters.
     */
    SharedMatrix solve(int max_iter=5, double conv=1.e-10, bool print=true);

protected:

    /// Parameters
    size_t nocc_;
    size_t nvir_;
    size_t nao_;

    /// Global JK object
    boost::shared_ptr<JK> jk_;

    /// Map of matrices
    std::map<std::string, SharedMatrix > matrices_;

}; // SOSCF class

} // Namespace psi

// class SOSCF {

// public:

//     /**
//      * What kind of wavefunction is this?
//      *
//      */
//     enum WaveType {RHF, RASSCF, CASSCF};

//     /**
//      * Initialize the SOSCF object
//      * @param nocc Number of doubly occupied orbitals
//      * @param nact Number of variably occupied orbitals
//      * @param nvir Number of virtual orbitals
//      * @param cas  Is this a CAS calculation? (ignore active-active rotations)
//      */
//     SOSCF(boost::shared_ptr<JK> jk);

//     /**
//      * RHF Constructor
//      */
//     //SOSCF(size_t wave_type, size_t nocc, size_t nfzv);

//     ~SOSCF(void);


//     *
//      * Updates all internal variables to the new reference frame and builds
//      * non-rotated intermediates. The orbitals should enter ordered like
//      * fzc, occ, act, vir, fzv
//      * @param C    Current orbitals
//      * @param OPDM Current active one-particle density matrix
//      * @param TPDM Current active two-particle density matrix (symmetrized, dense)

//     //void update(SharedMatrix C, SharedMatrix OPDM, SharedMatrix TPDM);

//     /**
//      * Same as for above, but for RHF
//      */
//     void update(SharedMatrix Cocc, SharedMatrix Cvir, SharedMatrix Fock);

//     /**
//      * Constructs the full antisymmetrix U matrix
//      * @param  x The occupied-active by active-virtual non-redundant orbital rotations
//      * @return   U
//      */
//     //SharedMatrix build_U_from_x(SharedMatrix x);

//     /**
//      * The inactive Fock matrix
//      * IFock = h_mn + \sum_i (2 g_mnii - g_miin)
//      * @return  IFock
//      */
//     //SharedMatrix IFock(void);

//     /**
//      * Returns the rotated inactive Fock matrix
//      * IFock_k = IF_mp K_np + IF_pn K_mp + (4 G_mnip - G_mpin - G_npim) K_ip
//      * @param  x The occupied-active by active-virtual non-redundant orbital rotations
//      * @return   IFock_k
//      */
//     //SharedMatrix IFock_k(SharedMatrix x);

#endif
