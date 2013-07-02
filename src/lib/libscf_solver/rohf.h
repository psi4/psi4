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

#ifndef __rohf_psi_h__
#define __rohf_psi_h__

#include <vector>
#include <libpsio/psio.hpp>
#include "hf.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi { namespace scf {

class ROHF : public HF {
protected:
    SharedMatrix Feff_;
    SharedMatrix soFeff_;
    SharedMatrix Dt_old_;
    SharedMatrix Dt_;
    SharedMatrix Ct_;
    SharedMatrix Ga_;
    SharedMatrix Gb_;
    SharedMatrix Ka_;
    SharedMatrix Kb_;
    SharedMatrix moFa_;
    SharedMatrix moFb_;

    /// Before semicanonicalize is called, this is true, but it becomes false
    bool restricted_;
    void form_initialF();
    void form_initial_C();
    void form_C();
    void form_D();
    double compute_initial_E();
    double compute_E();
    virtual void stability_analysis();
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

    void common_init();
public:
    ROHF(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    ROHF(Options& options, boost::shared_ptr<PSIO> psio);
    virtual ~ROHF();
    virtual bool same_a_b_orbs() const { return restricted_; }
    virtual bool same_a_b_dens() const { return false; }
};

}}

#endif
