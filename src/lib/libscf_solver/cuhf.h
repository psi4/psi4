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

#ifndef __math_test_cuhf_h__
#define __math_test_cuhf_h__

#include <libpsio/psio.hpp>
#include "hf.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi { namespace scf {

/*

    Constrained Unrestricted Hartree-Fock

    Reference: T. Tsuchimochi and G.E. Scuseria, J. Chem. Phys. 133,
               141102 (2010)

    This is an alternative formulation of ROHF as a contrained UHF. A
    Lagrangian constraint is placed on the usual UHF procedure to remove
    the spin contamination. The result is an ROHF energy and semicanonical
    ROHF orbitals. The need to pick coupling coefficients is removed.
    Koopmans' theorem is valid for CUHF orbital energes.

    It is claimed that CUHF does not suffer from the convergence problems
    of certain ROHF implementations (not sure how PSI's ROHF code does).
    CUHF retains the UHF-like trait that Ca != Cb. Also, the converged CUHF
    wavefunction yields the correct value for <S^2>, however, this is only
    true at convergence. It is possible that this increased flexibility
    improves convergence.

    -- EGH, August 15th, 2011

    TODO:

    Probably can't handle NSO != NMO right now, should either fix this code
    or the transform functions from Matrix.

    Using the UHF form for the Lagrangian, this is probably correct, but
    should be checked.

*/

class CUHF : public HF {
protected:
    SharedMatrix Dt_, Dt_old_;
    SharedMatrix Da_old_, Db_old_;
    SharedMatrix J_, Ka_, Kb_;
    // Contributions to the Fock matrix from charge and spin density
    SharedMatrix Fp_, Fm_;
    // Charge denisty and natural orbitals (eigenvectors of charge density)
    SharedMatrix Dp_, Cno_, Cno_temp_;
    // Natural orbital occupations
    SharedVector No_;

    void form_initialF();
    void form_C();
    void form_D();
    double compute_initial_E();
    virtual bool stability_analysis();
    virtual double compute_E();

    virtual void form_G();
    virtual void form_F();

    virtual void compute_orbital_gradient(bool save_diis);

    bool diis();

    bool test_convergency();
    void save_information();
    void compute_spin_contamination();

    void common_init();

    void save_density_and_energy();

    virtual void damp_update();

    // Finalize memory/files
    virtual void finalize();

public:
    CUHF(SharedWavefunction ref_wfn, Options& options, boost::shared_ptr<PSIO> psio);
    virtual ~CUHF();
};

}}

#endif
