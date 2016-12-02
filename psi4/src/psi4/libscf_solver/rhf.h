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

#ifndef RHF_H
#define RHF_H

#include "psi4/libpsio/psio.hpp"
#include "hf.h"

namespace psi {

class TwoBodySOInt;
class PSIO;
class Chkpt;
class Matrix;
class Vector;

namespace scf {

class RHF : public HF {
protected:

    // Temporary matrices
    SharedMatrix D_;
    SharedMatrix Dold_;
    SharedMatrix G_;
    SharedMatrix J_;
    SharedMatrix K_;
    SharedMatrix wK_;

    void form_C();
    void form_D();
    virtual void damp_update();
    double compute_initial_E();
    virtual double compute_E();
    virtual bool stability_analysis();

    virtual void form_F();
    virtual void form_G();
    virtual void form_V();
    virtual void compute_orbital_gradient(bool save_fock);

    bool diis();

    bool test_convergency();

    void common_init();

    // Finalize memory/files
    virtual void finalize();

    void save_density_and_energy();

    // Second-order convergence code
    virtual int soscf_update(void);

public:
    RHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional);
    RHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional,
        Options& options, std::shared_ptr<PSIO> psio);
    virtual ~RHF();

    virtual SharedMatrix Da() const;

    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }

    /// Hessian-vector computers and solvers
    virtual std::vector<SharedMatrix> onel_Hx(std::vector<SharedMatrix> x);
    virtual std::vector<SharedMatrix> twoel_Hx(std::vector<SharedMatrix> x, bool combine = true,
                                               std::string return_basis = "MO");
    virtual std::vector<SharedMatrix> cphf_Hx(std::vector<SharedMatrix> x);
    virtual std::vector<SharedMatrix> cphf_solve(std::vector<SharedMatrix> x_vec,
                                                 double conv_tol = 1.e-4, int max_iter = 10,
                                                 int print_lvl = 1);
};

}}

#endif
