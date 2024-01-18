/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#ifndef GauXC_FUNCTIONAL_H
#define GauXC_FUNCTIONAL_H

#include "psi4/libfunctional/functional.h"
#include "psi4/libmints/typedefs.h"

#include <map>

#include <gauxc/types.hpp>

/**
 * GauXC functional wrapper
 **/
namespace psi {

class GauXCFunctional : public Functional {
    // Wrapper to the GauXC library
    
      std::unique_ptr<GauXC::functional_type> func_;

   public:
    GauXCFunctional(std::string xc_name, bool unpolarized);
    ~GauXCFunctional() override;

    void compute_functional(const std::map<std::string, SharedVector>& in,
                            const std::map<std::string, SharedVector>& out, int npoints, int deriv) override;

    // Clones a *polarized*, complete functional. Used, e.g., in spin-symmetry-
    // breaking eigenvectors of the MO hessian or linear response eigenproblem.
    std::shared_ptr<Functional> build_polarized() override;

};
}  // namespace psi

#endif
