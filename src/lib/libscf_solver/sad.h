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

#ifndef LIBSCF_SAD_H
#define LIBSCF_SAD_H

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class BasisSet;
class Molecule;
class Matrix;

namespace scf {

class SADGuess {

protected:

    int print_;
    int debug_;

    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<BasisSet> basis_;
    SharedMatrix AO2SO_;

    int nalpha_;
    int nbeta_;

    Options& options_;

    SharedMatrix Da_;
    SharedMatrix Db_;
    SharedMatrix Ca_;
    SharedMatrix Cb_;

    void common_init();

    SharedMatrix form_D_AO();
    void form_gradient(int norbs, SharedMatrix grad, SharedMatrix F, SharedMatrix D,
                      SharedMatrix S, SharedMatrix X);
    void get_uhf_atomic_density(boost::shared_ptr<BasisSet> atomic_basis,
                                int n_electrons, int multiplicity, SharedMatrix D);
    void form_C_and_D(int nocc, int norbs, SharedMatrix X, SharedMatrix F,
                                  SharedMatrix C, SharedMatrix Cocc, SharedVector occ,
                                  SharedMatrix D);

    void form_D();
    void form_C();

public:

    SADGuess(boost::shared_ptr<BasisSet> basis, int nalpha, int nbeta, Options& options);
    virtual ~SADGuess();

    void compute_guess();

    SharedMatrix Da() const { return Da_; }
    SharedMatrix Db() const { return Db_; }
    SharedMatrix Ca() const { return Ca_; }
    SharedMatrix Cb() const { return Cb_; }

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }

};

}}

#endif
