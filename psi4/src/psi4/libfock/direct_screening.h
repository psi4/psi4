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

#ifndef __psi_libfock_direct_screening__
#define __psi_libfock_direct_screening__

#include <cfloat>

#include "psi4/psi4-dec.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"

#include "psi4/libpsio/psio.hpp"

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"

#include "psi4/libfock/apps.h"
#include "psi4/libfock/jk.h"

#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/sieve.h"


namespace psi {

  class DirectScreening {

  protected:

    std::shared_ptr<BasisSet> basis_;

    bool do_J_;
    bool do_K_;
    bool do_wK_;

    //std::vector<SharedMatrix >& D_;
    std::vector<SharedMatrix> D_;

    // the coulomb and exchange matrices, to be filled in and passed back out
    std::vector<SharedMatrix> J_;
    std::vector<SharedMatrix> K_;
    //std::vector<SharedMatrix> wK_;

    // ERIs
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri_;
    /// Integral factory (must be retained for Spherical Transforms)
    std::shared_ptr<IntegralFactory> factory_;
    /// ERI Sieve
    std::shared_ptr<ERISieve> sieve_;

    ///////////////////// functions ////////////////////////////


  public:

    DirectScreening(std::shared_ptr<BasisSet> basis_in,
         std::vector<SharedMatrix>& density_in);

    ~DirectScreening();

    void Compute();

    // Assuming this always gets called before Compute()
    void Update(const std::vector<SharedMatrix>& D_new);

    std::vector<SharedMatrix>& J();

    std::vector<SharedMatrix>& K();

    void set_do_J(bool do_it);

    void set_do_K(bool do_it);

    void print_header() const;



  }; // class

}


#endif /* defined(____direct_screening__) */
