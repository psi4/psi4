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

#ifndef _jk_independent_h
#define _jk_independent_h

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

  template <class JDriver, class KDriver>
  class JKIndependent : public JK
  {

  protected:

    JDriver j_driver_;
    KDriver k_driver_;

    // true if we need to call j_driver_ for J computation and k_driver_ for K.
    // if false, then we assume j_driver does everything we need it to
    bool do_separately_;

    /// Do we need to backtransform to C1 under the hood?
    virtual bool C1() const { return allow_desymmetrization_; }
    /// Setup integrals, files, etc
    virtual void preiterations();
    /// Compute J/K for current C/D
    virtual void compute_JK();
    /// Delete integrals, files, etc
    virtual void postiterations();

    /// Common initialization
    void common_init();

  public:
    // => Constructors < = //

    JKIndependent(std::shared_ptr<BasisSet> primary, bool do_separately);
    /// Destructor
    virtual ~JKIndependent();

    // => Accessors <= //

    /**
     * Print header information regarding JK
     * type on output file
     */
    virtual void print_header() const;

  }; // class

} // namespace

// definitions
#include "jk_independent_impl.h"

#endif


