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

#ifndef __math_test_uhf_h__
#define __math_test_uhf_h__

#include "psi4/libpsio/psio.hpp"
#include "hf.h"

namespace psi {
class Matrix;
class Vector;
namespace scf {

class UHF : public HF {
protected:
    SharedMatrix Dt_, Dt_old_;
    SharedMatrix Da_old_, Db_old_;
    SharedMatrix Ga_, Gb_, J_, Ka_, Kb_, wKa_, wKb_;

    void form_initialF();
    void form_C();
    void form_V();
    void form_D();
    double compute_initial_E();
    virtual double compute_E();
    virtual bool stability_analysis();
    bool stability_analysis_pk();

    virtual void form_G();
    virtual void form_F();

    virtual void compute_orbital_gradient(bool save_diis);
    bool diis();

    bool test_convergency();
    void save_information();

    void common_init();

    void save_density_and_energy();

    // Finalize memory/files
    virtual void finalize();

    // Scaling factor for orbital rotation
    double step_scale_;
    // Increment to explore different scaling factors
    double step_increment_;
    // Stability eigenvalue, for doing smart eigenvector following
    double stab_val;

    // Compute UHF NOs
    void compute_nos();

    // Damp down the density update
    virtual void damp_update();

    // Second-order convergence code
    void Hx(SharedMatrix x_a, SharedMatrix IFock_a, SharedMatrix Cocc_a,
            SharedMatrix Cvir_a, SharedMatrix ret_a,
            SharedMatrix x_b, SharedMatrix IFock_b, SharedMatrix Cocc_b,
            SharedMatrix Cvir_b, SharedMatrix ret_b);
    virtual int soscf_update(void);

public:
    UHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional);
    UHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional,
        Options& options, std::shared_ptr<PSIO> psio);
    virtual ~UHF();

    virtual bool same_a_b_orbs() const { return false; }
    virtual bool same_a_b_dens() const { return false; }
};

}}

#endif
