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

#ifndef __rohf_psi_h__
#define __rohf_psi_h__

#include <vector>
#include "psi4/libpsio/psio.hpp"
#include "hf.h"

namespace psi { namespace scf {

class ROHF : public HF {
protected:
    SharedMatrix moFeff_;
    SharedMatrix soFeff_;
    SharedMatrix Dt_;
    SharedMatrix Da_old_;
    SharedMatrix Db_old_;
    SharedMatrix Dt_old_;
    SharedMatrix Ct_;
    SharedMatrix Ga_;
    SharedMatrix Gb_;
    SharedMatrix Ka_;
    SharedMatrix Kb_;
    SharedMatrix moFa_;
    SharedMatrix moFb_;

    void form_initialF();
    void form_initial_C();
    void form_C();
    void form_D();
    double compute_initial_E();
    double compute_E();
    virtual bool stability_analysis();
    virtual void prepare_canonical_orthogonalization();
    void semicanonicalize();

    void form_G();
    void form_F();

    virtual void compute_orbital_gradient(bool save_diis);
    bool diis();

    bool test_convergency();

    void save_information();
    // Finalize memory/files
    virtual void finalize();

    void save_density_and_energy();
    void format_guess();

    // Second-order convergence code
    void Hx(SharedMatrix x, SharedMatrix ret);
    virtual int soscf_update(void);

    /** Applies damping to the density update */
    virtual void damp_update();

    void common_init();
public:
    ROHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional);
    ROHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional,
         Options& options, std::shared_ptr<PSIO> psio);
    virtual ~ROHF();

    SharedMatrix moFeff() const {return moFeff_; }
    SharedMatrix moFa() const {return moFa_; }
    SharedMatrix moFb() const {return moFb_; }

};

}}

#endif
