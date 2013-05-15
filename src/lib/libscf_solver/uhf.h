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

#ifndef __math_test_uhf_h__
#define __math_test_uhf_h__

#include <libpsio/psio.hpp>
#include "hf.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {
class Matrix;
class Vector;
namespace scf {

class UHF : public HF {
protected:
    SharedMatrix Dt_, Dtold_;
    SharedMatrix Ga_, Gb_, J_, Ka_, Kb_;

    void form_initialF();
    void form_C();
    void form_D();
    double compute_initial_E();
    virtual double compute_E();
    virtual void stability_analysis();

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

public:
    UHF(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    UHF(Options& options, boost::shared_ptr<PSIO> psio);
    virtual ~UHF();

    virtual bool same_a_b_orbs() const { return false; }
    virtual bool same_a_b_dens() const { return false; }
};

}}

#endif
