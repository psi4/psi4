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

/*
 * deriv.h
 *
 *  Created on: Feb 24, 2009
 *      Author: jturney
 */

#ifndef _psi_src_lib_libmints_deriv_h_
#define _psi_src_lib_libmints_deriv_h_

#include <vector>
#include "matrix.h"
#include "psi4/libmints/cdsalclist.h"

namespace psi {

class Wavefunction;
class IntegralFactory;
class BasisSet;
class SOBasisSet;
class Molecule;
class CdSalcList;

class Deriv
{
    const std::shared_ptr<Wavefunction> wfn_;
    std::shared_ptr<IntegralFactory> integral_;
    std::shared_ptr<BasisSet> basis_;
    std::shared_ptr<SOBasisSet> sobasis_;
    std::shared_ptr<MatrixFactory> factory_;
    std::shared_ptr<Molecule> molecule_;

    CdSalcList cdsalcs_;

    SharedMatrix P_2_;
    SharedMatrix W_2_;
    SharedMatrix SCF_D_;

    int natom_;
    bool tpdm_presorted_;
    bool deriv_density_backtransformed_;
    bool ignore_reference_;

    std::vector<SharedMatrix> dH_;
    std::vector<SharedMatrix> dS_;

    // Results go here.
    /// One-electron contribution to the gradient
    SharedMatrix opdm_contr_;
    /// Reference one-electron contribution to the gradient
    SharedMatrix opdm_ref_contr_;
    /// Overlap contribution to the gradient
    SharedMatrix x_contr_;
    /// Reference overlap contribution to the gradient
    SharedMatrix x_ref_contr_;
    /// Two-electron contribution to the gradient
    SharedMatrix tpdm_contr_;
    /// Reference two-electron contribution to the gradient
    SharedMatrix tpdm_ref_contr_;
    /// Final gradient
    SharedMatrix gradient_;
public:
    /*!
     * Constructor for derivative object.
     *
     * \param wave Wavefunction object to compute derivative for
     * \param needed_irreps by default do A1 derivatives
     * \param project_out_translations remove translations from the CdSalcs
     * \param project_out_rotations remove rotations from the CdSalcs
     */
    Deriv(const std::shared_ptr<Wavefunction>& wave,
          char needed_irreps=0x1,
          bool project_out_translations=false,
          bool project_out_rotations=false);

    // Is the TPDM already presorted? Default: False
    void set_tpdm_presorted(bool val) { tpdm_presorted_ = val; }

    // Ignore reference contributions to the gradient? Default: False
    void set_ignore_reference(bool val) { ignore_reference_ = val; }

    // Is the deriv_density already backtransformed? Default: False
    void set_deriv_density_backtransformed(bool val) { deriv_density_backtransformed_ = val; }

    SharedMatrix compute();

    const SharedMatrix& one_electron() {
        return opdm_contr_;
    }

    const SharedMatrix& lagrangian() {
        return x_contr_;
    }

    const SharedMatrix& two_body() {
        return tpdm_contr_;
    }

    const SharedMatrix& gradient() {
        return gradient_;
    }
};

}

#endif /* DERIV_H_ */
