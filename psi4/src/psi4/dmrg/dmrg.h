/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#ifndef _PSI4_SRC_BIN_DMRG_DMRG_H_
#define _PSI4_SRC_BIN_DMRG_DMRG_H_

#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/wavefunction.h"

namespace psi {
namespace dmrg {
class DMRGSolver : public Wavefunction {
   public:
    DMRGSolver(SharedWavefunction ref_wfn, Options& options);
    double compute_energy() override;

    /// Returns the occupation vectors
    std::shared_ptr<Vector> occupation_a() const;
    std::shared_ptr<Vector> occupation_b() const { return occupation_a(); };
};
}  // namespace dmrg
}  // namespace psi
#endif  // Header guard
