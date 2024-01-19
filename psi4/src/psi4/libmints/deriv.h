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

// Enum used to specify the type of derivative computation
// Default:     Use internal logic
// Correlated:  Correlated methods that write RDMs and Lagrangian to disk
enum class DerivCalcType { Default,
    Correlated };

class PSI_API Deriv {
    const std::shared_ptr<Wavefunction> wfn_;
    std::shared_ptr<IntegralFactory> integral_;
    std::shared_ptr<BasisSet> basis_;
    std::shared_ptr<SOBasisSet> sobasis_;
    std::shared_ptr<MatrixFactory> factory_;
    std::shared_ptr<Molecule> molecule_;

    CdSalcList cdsalcs_;

    int natom_;
    bool tpdm_presorted_;
    bool deriv_density_backtransformed_;
    bool ignore_reference_;

    // Results go here.
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
    Deriv(const std::shared_ptr<Wavefunction>& wave, char needed_irreps = 0x1, bool project_out_translations = false,
          bool project_out_rotations = false);

    // Is the TPDM already presorted? Default: False
    void set_tpdm_presorted(bool val) { tpdm_presorted_ = val; }

    // Ignore reference contributions to the gradient? Default: False
    void set_ignore_reference(bool val) { ignore_reference_ = val; }

    // Is the deriv_density already backtransformed? Default: False
    void set_deriv_density_backtransformed(bool val) { deriv_density_backtransformed_ = val; }

    SharedMatrix compute(DerivCalcType deriv_calc_type = DerivCalcType::Default);

    /*!
     * Computes gradients assuming density-fitted two electron integrals.
     * All density matrix intermediates are assumed stored on the wavefunction or written to disk.
     * That is the responsibility of the caller code. See deriv.cc for expected format.
     * Generating integral derivatives and contracting is the job of this code.
     * \param ref_aux_name Name of reference auxiliary basis set, e.g., DF_BASIS_SCF
     * \param cor_aux_name Name of correlated auxiliary basis set, e.g., DF_BASIS_CC
     */
    SharedMatrix compute_df(const std::string& ref_aux_name, const std::string& cor_aux_name);

    const SharedMatrix& two_body() { return tpdm_contr_; }

    const SharedMatrix& gradient() { return gradient_; }
};

}  // namespace psi

#endif /* DERIV_H_ */
