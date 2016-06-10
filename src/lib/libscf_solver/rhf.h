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

#ifndef RHF_H
#define RHF_H

#include <libpsio/psio.hpp>
#include "hf.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class TwoBodySOInt;
class PSIO;
class Chkpt;
class Matrix;
class Vector;

namespace scf {

class RHF : public HF {
protected:
    SharedMatrix D_;
    SharedMatrix Dold_;
    SharedMatrix G_;
    SharedMatrix J_;
    SharedMatrix K_;


    void form_C();
    void form_D();
    virtual void damp_update();
    double compute_initial_E();
    virtual double compute_E();
    virtual bool stability_analysis();

    virtual void form_F();
    virtual void form_G();
    virtual void compute_orbital_gradient(bool save_fock);

    bool diis();

    bool test_convergency();
    void save_information();

    void common_init();

    // Finalize memory/files
    virtual void finalize();

    void save_density_and_energy();

    // Second-order convergence code
    void Hx(SharedMatrix x, SharedMatrix IFock, SharedMatrix Cocc, SharedMatrix Cvir, SharedMatrix ret);
    virtual int soscf_update(void);

public:
    RHF(SharedWavefunction ref_wfn, Options& options, boost::shared_ptr<PSIO> psio);
    virtual ~RHF();


    virtual SharedMatrix Da() const;

    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }
};

}}

#endif