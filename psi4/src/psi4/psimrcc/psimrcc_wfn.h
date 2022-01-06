/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#ifndef _psi_src_bin_psimrcc_wfn_h_
#define _psi_src_bin_psimrcc_wfn_h_

#include "psi4/libmints/wavefunction.h"

#include "psi4/libqt/qt.h"
#define START_TIMER(a) timer_on((a));
#define END_TIMER(a) timer_off((a));

namespace psi {

class MOInfo;

namespace psimrcc {

class CCBLAS;

class PSIMRCCWfn : public Wavefunction {
   public:
    PSIMRCCWfn(SharedWavefunction ref_wfn, Options &options);
    double compute_energy() override;

    // Methods
    const std::shared_ptr<MOInfo> moinfo() const { return moinfo_; }
    const std::shared_ptr<CCBLAS> blas() const { return blas_; }
    // Estimate the free memory.
    size_t free_memory_;

   protected:
    // Class members
    std::shared_ptr<MOInfo> moinfo_;
    std::shared_ptr<CCBLAS> blas_;

    // Methods
    void active_space_warning() const;
};
}
}
#endif  // _psi_src_bin_psimrcc_wfn_h_
