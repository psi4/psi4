/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <libmints/typedefs.h>
#include <map>
#include <boost/tuple/tuple.hpp>

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
     * @param H       Core hamiltonian in the SO basis.
     * @param casscf  Is this a CAS calculation? (ignore active-active rotations)
     */
    // SOMCSCF(boost::shared_ptr<JK> jk, SharedMatrix H, bool casscf);
    SOMCSCF(boost::shared_ptr<JK> jk, SharedMatrix AOTOSO,
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
     * @return gradient_rms Returns the RMS of the gradient.
     */
    double gradient_rms();

protected:

    /// Parameters
    bool casscf_;
    bool has_fzc_;

    /// Doubles
    double efzc_;
    double edrc_;

    /// Orbital info
    size_t nocc_;
    Dimension noccpi_;
    size_t nact_;
    Dimension nactpi_;
    size_t nvir_;
    Dimension nvirpi_;

    size_t nirrep_;
    size_t nmo_;
    Dimension nmopi_;
    size_t nso_;
    size_t nao_;
    Dimension nsopi_;

    // Non-redunant rotatations
    Dimension noapi_;
    Dimension navpi_;

    /// Integral objects
    boost::shared_ptr<JK> jk_;

    /// Map of matrices
    std::map<std::string, SharedMatrix > matrices_;

    /// RAS arrays
    std::vector<Dimension> ras_spaces_;
    void check_ras();

    /// Zero's out redundant rotations, will chose act or ras.
    void zero_redundant(SharedMatrix vector);
    void zero_act(SharedMatrix vector);
    void zero_ras(SharedMatrix vector);

    // Transform the integrals
    virtual void transform(bool approx_only);

    // Grab actMO (dense)
    virtual void set_act_MO();

    // Build the Q matrices
    virtual void compute_Q();
    virtual void compute_Qk(SharedMatrix U, SharedMatrix Uact);


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
     * @param df      AOTOSO object to use.
     * @param H       Core hamiltonian in the SO basis.
     */
    DFSOMCSCF(boost::shared_ptr<JK> jk, boost::shared_ptr<DFERI> df, SharedMatrix AOTOSO,
            SharedMatrix H);

    virtual ~DFSOMCSCF();

protected:

    boost::shared_ptr<DFERI> dferi_;
    virtual void transform(bool approx_only);
    virtual void set_act_MO();
    virtual void compute_Q();
    virtual void compute_Qk(SharedMatrix U, SharedMatrix Uact);

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
    DiskSOMCSCF(boost::shared_ptr<JK> jk, boost::shared_ptr<IntegralTransform> ints, SharedMatrix AOTOSO, SharedMatrix H);

    virtual ~DiskSOMCSCF();

protected:

    boost::shared_ptr<IntegralTransform> ints_;
    boost::shared_ptr<PSIO>  psio_;
    virtual void transform(bool approx_only);
    virtual void set_act_MO();
    virtual void compute_Q();
    virtual void compute_Qk(SharedMatrix U, SharedMatrix Uact);

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
    IncoreSOMCSCF(boost::shared_ptr<JK> jk, SharedMatrix AOTOSO, SharedMatrix H);

    virtual ~IncoreSOMCSCF();

protected:

    virtual void set_act_MO();
    virtual void compute_Q();
    virtual void compute_Qk(SharedMatrix U, SharedMatrix Uact);

    void set_eri_tensors(SharedMatrix aaaa, SharedMatrix aaar);
    bool eri_tensor_set_;
    SharedMatrix mo_aaaa_;
    SharedMatrix mo_aaar_;

}; // DFSOMCSCF class


} // Namespace psi


#endif