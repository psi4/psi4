/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#include "psi4/libmints/wavefunction.h"

namespace psi {
namespace dlpnomp2 {

// References to DiStasio are to J Comput Chem 28: 839–856, 2007; DOI: 10.1002/jcc.20604

class DLPNOMP2 : public Wavefunction {
   protected:
    // Auxiliary basis
    std::shared_ptr<BasisSet> ribasis_;

    SharedMatrix Cfocc_;
    SharedMatrix Caocc_;
    SharedMatrix Cavir_;
    SharedMatrix Cfvir_;

    SharedVector eps_focc_;
    SharedVector eps_aocc_;
    SharedVector eps_avir_;
    SharedVector eps_fvir_;

    void common_init();

    void print_header();

   public:
    DLPNOMP2(SharedWavefunction ref_wfn, Options& options);
    ~DLPNOMP2() override;

    double compute_energy();
};

}  // namespace dlpnomp2
}  // namespace psi
