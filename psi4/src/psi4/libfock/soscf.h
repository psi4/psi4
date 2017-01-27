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
 * @END LICENSE
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

#include "psi4/libmints/typedefs.h"
#include <map>

namespace psi {

/// Forward declare
class JK;
class DFERI;
class IntegralTransform;
class PSIO;

class SOMCSCF {

public:

    /**
     * Initialize the SOMCSCF object
     * @param jk      JK object to use.
     * @param AOTOSO  AOTOSO object
     * @param H       Core hamiltonian in the SO basis.
     */
    // SOMCSCF(std::shared_ptr<JK> jk, SharedMatrix H, bool casscf);
    SOMCSCF(std::shared_ptr<JK> jk, SharedMatrix AOTOSO,
            SharedMatrix H);

    virtual ~SOMCSCF(void);

    /**
     * Sets the ras spaces, rotations will not happen inside of a space
     * @param ras_spaces Standard vector of arbitrary length with Dimension objects
     */
    void set_ras(std::vector<Dimension> ras_spaces);

    /**
     * Sets the frozen core orbitals, these orbitals do not rotate.
     * @param Cfzc The frozen core orbitals
     */
    void set_frozen_orbitals(SharedMatrix Cfzc);

    /**
     * Sets the AO based IFock matrix. It should be noted that SOMCSCF will never compute
     * the IFock matrix again, therefore this should be updated on every update() call.
     * Frozen core orbitals should be included
     * @param IFock The AO based IFock matrix
     */
    void set_AO_IFock(SharedMatrix IFock);
    /**
     * Forms the rotation matrix exp(U).
     * @param  x The [oa, av] non-redundant orbital rotation parameters.
     * @return   The rotation matrix.
     */
    SharedMatrix form_rotation_matrix(SharedMatrix x, size_t order = 2);

    /**
     * Rotate the current orbitals for a given rotation matrix.
     * @param  C The orbital matrix to rotate.
     * @param  x The [oa, av] non-redundant orbital rotation parameters.
     * @return   The rotated orbitals.
     */
    SharedMatrix Ck(SharedMatrix C, SharedMatrix x);

    /**
     * Computes the RHF energy for a given C matrix
     * @param  C The desired orbitals
     * @return   The RHF energy
     */
    double rhf_energy(SharedMatrix C);

    /**
     * Updates all internal variables to the new reference frame and builds
     * non-rotated intermediates.
     * @param Cocc Current doubly occupied orbitals
     * @param Cact Current active orbitals
     * @param Cvir Current virtual orbitals
     * @param OPDM Current active one-particle density matrix
     * @param TPDM Current active two-particle density matrix (symmetrized, dense)
    */
    void update(SharedMatrix Cocc, SharedMatrix Cact, SharedMatrix Cvir,
                SharedMatrix OPDM, SharedMatrix TPDM);

    /**
     * Sets the occupied Fock Matrix
     * @param occ_Fock The current doubly occupied Fock matrix
     */
    void set_occ_fock(SharedMatrix occ_Fock);

    /**
     * Returns the approximate diagonal hessian.
     * @return Hdiag The [oa, av] block of the diagonal hessian.
     */
    SharedMatrix H_approx_diag();

    /**
     * Returns the hessian times a trial vector or, in effect, the rotated inactive Fock matrix.
     * @param  x  The [oa, av] matrix of non-redundant orbital rotations.
     * @return Hx The [oa, av] block of the rotated Fock matrix.
     */
    SharedMatrix Hk(SharedMatrix x);

    /**
     * Uses the approximate H diagonal hessian for an update.
     * @return x         The [oa, av] matrix of non-redundant orbital rotation parameters.
     */
    SharedMatrix approx_solve();

    /**
     * Solves the set of linear equations Hx = gradient using CG.
     * @return x The [oa, av] matrix of non-redundant orbital rotation parameters.
     */
    SharedMatrix solve(int max_iter=5, double conv=1.e-10, bool print=true);

    /**
     * @return gradient Returns the MO gradient.
     */
    SharedMatrix gradient();

    /**
     * Computes the Q matrix Q_pw = (pt|uv) \Gamma_{tuvw}
     * @param  TPDM Dense nact*nact by nact*nact symmetrized TPDM
     * @return      The symmetry blocked Q matrix
     */
    virtual SharedMatrix compute_Q(SharedMatrix TPDM);


    /**
     * Computes the Qk matrix Q_pw = (pt|uv)^k \Gamma_{tuvw} with rotated integrals
     * @param  TPDM Dense nact*nact by nact*nact symmetrized TPDM
     * @param  U    The full rotation matrix U
     * @param  Uact The active portion of the matrix U
     * @return      The symmetry blocked Q matrix
     */
    virtual SharedMatrix compute_Qk(SharedMatrix TPDM, SharedMatrix U, SharedMatrix Uact);

    /**
     * Computes the Q matrix AFock_pq = (pq|uv) - 0.5 (pu|qv) \gamma_{uv}
     * @param  OPDM Dense nact*nact by nact*nact symmetrized OPDM
     * @return      The symmetry AFock matrix
     */
    SharedMatrix compute_AFock(SharedMatrix OPDM);

    /**
     * @return gradient_rms Returns the RMS of the gradient.
     */
    double gradient_rms();

    /**
     * Zeros out the redundant rotations
     * @param vector Zero redundant rotations
     */
    void zero_redundant(SharedMatrix vector);

    double current_total_energy() { return (energy_drc_ + energy_ci_); }
    double current_docc_energy() { return energy_drc_; }
    double current_ci_energy() { return energy_ci_; }
    SharedMatrix current_AFock() { return matrices_["AFock"]; }
    SharedMatrix current_IFock() { return matrices_["IFock"]; }
    virtual void set_eri_tensors(SharedMatrix, SharedMatrix)
    {
    }

protected:

    /// Parameters
    bool casscf_;
    bool has_fzc_;
    bool compute_IFock_;

    /// Doubles
    double energy_drc_;
    double energy_ci_;

    /// Orbital info
    size_t nocc_;
    size_t nact_;
    size_t nvir_;
    Dimension noccpi_;
    Dimension nactpi_;
    Dimension nvirpi_;

    size_t nirrep_;
    size_t nmo_;
    size_t nso_;
    size_t nao_;
    Dimension nmopi_;
    Dimension nsopi_;

    // Non-redunant rotatations
    Dimension noapi_;
    Dimension navpi_;

    /// Integral objects
    std::shared_ptr<JK> jk_;

    /// Map of matrices
    std::map<std::string, SharedMatrix > matrices_;

    /// RAS arrays
    std::vector<Dimension> ras_spaces_;
    void check_ras();

    /// Zero's out redundant rotations, will chose act or ras.
    void zero_act(SharedMatrix vector);
    void zero_ras(SharedMatrix vector);

    // Transform the integrals
    virtual void transform(bool approx_only);

    // Grab actMO (dense)
    virtual void set_act_MO();


}; // SOMCSCF class


/**
 * Class DFSOMCSCF
 *
 * Density fitted second-order MCSCF
 */
class DFSOMCSCF : public SOMCSCF {

public:
    /**
     * Initialize the DF SOMCSCF object
     * @param jk      JK object to use.
     * @param df      DFERI object to use.
     * @param H       Core hamiltonian in the SO basis.
     */
    DFSOMCSCF(std::shared_ptr<JK> jk, std::shared_ptr<DFERI> df, SharedMatrix AOTOSO,
            SharedMatrix H);

    virtual ~DFSOMCSCF();

protected:

    std::shared_ptr<DFERI> dferi_;
    virtual void transform(bool approx_only);
    virtual void set_act_MO();
    virtual SharedMatrix compute_Q(SharedMatrix TPDM);
    virtual SharedMatrix compute_Qk(SharedMatrix TPDM, SharedMatrix U, SharedMatrix Uact);

}; // DFSOMCSCF class

/**
 * Class DiskSOMCSCF
 *
 * Disk-based fitted second-order MCSCF
 */
class DiskSOMCSCF : public SOMCSCF {

public:
    /**
     * Initialize the DF SOMCSCF object.
     * @param jk      JK object to use.
     * @param H       Core hamiltonian in the SO basis.
     */
    DiskSOMCSCF(std::shared_ptr<JK> jk, std::shared_ptr<IntegralTransform> ints, SharedMatrix AOTOSO, SharedMatrix H);

    virtual ~DiskSOMCSCF();

protected:

    std::shared_ptr<IntegralTransform> ints_;
    std::shared_ptr<PSIO>  psio_;
    virtual void transform(bool approx_only);
    virtual void set_act_MO();
    virtual SharedMatrix compute_Q(SharedMatrix TPDM);
    virtual SharedMatrix compute_Qk(SharedMatrix TPDM, SharedMatrix U, SharedMatrix Uact);

}; // DiskSOMCSCF class

/**
 * Class IncoreSOMCSCF
 *
 * Second-order MCSCF using inc ore tensors
 * Note: set_eri_tensors should be called _before_ update.
 */
class IncoreSOMCSCF : public SOMCSCF {

public:
    /**
     * Initialize the DF SOMCSCF object.
     * @param jk      JK object to use.
     * @param H       Core hamiltonian in the SO basis.
     */
    IncoreSOMCSCF(std::shared_ptr<JK> jk, SharedMatrix AOTOSO, SharedMatrix H);

    virtual ~IncoreSOMCSCF();
    virtual void set_eri_tensors(SharedMatrix aaaa, SharedMatrix aaar);

protected:

    virtual void set_act_MO();
    virtual SharedMatrix compute_Q(SharedMatrix TPDM);
    virtual SharedMatrix compute_Qk(SharedMatrix TPDM, SharedMatrix U, SharedMatrix Uact);

    bool eri_tensor_set_;
    SharedMatrix mo_aaaa_;
    SharedMatrix mo_aaar_;

};  ///Incore SOMCSCF


} // Namespace psi


#endif
