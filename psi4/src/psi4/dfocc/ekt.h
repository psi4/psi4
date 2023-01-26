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

#ifndef _dfocc_ekt_h_
#define _dfocc_ekt_h_

#include "psi4/libpsio/psio.h"
#include "psi4/libmints/typedefs.h"
#include "tensors.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

class Ektip;
typedef std::shared_ptr<Ektip> SharedEktip;

class Ektip {
   private:
    std::string name_;

    int nocc_;
    int nvir_;
    int norb_;
    int plevel_;

    double cutoff_;
    double scale_;
    double scale2_;

    // 1D
    SharedTensor1d diagG1_;
    SharedTensor1d ps_vec_;
    SharedTensor1d eorb_;
    SharedTensor1d eocc_;
    SharedTensor1d ps_occ_;

    // 2D
    SharedTensor2d GF_;
    SharedTensor2d GFt_;
    SharedTensor2d GF_copy_;
    SharedTensor2d GFp_;
    SharedTensor2d G1_;
    SharedTensor2d G1t_;
    SharedTensor2d G1_copy_;
    SharedTensor2d G1half_;
    SharedTensor2d Uvec_;
    SharedTensor2d Uvecp_;
    SharedTensor2d temp_;
    SharedTensor2d PS_;
    SharedTensor2d GCt_;

    // ektip
    void compute_ektip();

   public:
    Ektip(std::string name, int nocc, int norb, const SharedTensor2d& GFock, const SharedTensor2d& Gamma,
          double scale_gf, double scale_ps);
    ~Ektip();

    // Return occupied energies
    SharedTensor1d eocc() const { return eocc_; }

    // Return all energies
    SharedTensor1d eorb() const { return eorb_; }

    // Return occupied pole strengths
    SharedTensor1d psocc() const { return ps_occ_; }

    // Return all pole strengths
    SharedTensor1d ps() const { return ps_vec_; }
};

}  // namespace dfoccwave
}  // namespace psi
#endif  // _dfocc_ekt_h_
