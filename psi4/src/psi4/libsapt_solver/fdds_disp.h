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

#ifndef FDDS_DISP_H
#define FDDS_DISP_H

#include "psi4/libmints/typedefs.h"

namespace psi {

class BasisSet;
class DFHelper;

namespace sapt {

class FDDS_Dispersion {
   protected:
    // BasisSets
    std::shared_ptr<BasisSet> primary_;
    std::shared_ptr<BasisSet> auxiliary_;

    // Coulomb metric inverse
    SharedMatrix metric_inv_;

    // Coulomb metric
    SharedMatrix metric_;

    // Auxiliary overlap matrix
    SharedMatrix aux_overlap_;

    // DFHelper object
    std::shared_ptr<DFHelper> dfh_;

    // Cache map
    std::map<std::string, SharedMatrix> matrix_cache_;
    std::map<std::string, SharedVector> vector_cache_;

   public:
    /**
     * Constructs the FDDS_Dispersion object.
     * @param primary   The primary basis
     * @param auxiliary The auxiliary basis
     * @param cache     A data cache containing "Cocc_A", "Cvir_A", "eps_occ_A", "eps_vir_A", "Cocc_B",
     * "Cvir_B", "eps_occ_B", "eps_vir_B" quantities
     */
    FDDS_Dispersion(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                    std::map<std::string, SharedMatrix> matrix_cache, std::map<std::string, SharedVector> vector_cache);

    ~FDDS_Dispersion();

    /**
     * Projects the densities from primary AO space to auxiliary AO space
     * @param  dens     Vector of densities to transform
     * @return axu_dens Vector of transformed densities
     */
    std::vector<SharedMatrix> project_densities(std::vector<SharedMatrix> dens);

    /**
     * Forms the uncoupled amplitude Qia,Eia,Pia->PQ
     * @param  monomer Monomer "A" or "B"
     * @param  omega   Time dependent value
     * @return         "PQ" amplitude tensor
     */
    SharedMatrix form_unc_amplitude(std::string monomer, double omega);

    /**
     * Returns the metric matrix
     * @return Metric
     */
    SharedMatrix metric() { return metric_; }

    /**
     * Returns the metric_inv matrix
     * @return Metric
     */
    SharedMatrix metric_inv() { return metric_inv_; }

    /**
     * Returns the auxiliary overlap matrix
     * @return Overlap
     */
    SharedMatrix aux_overlap() { return aux_overlap_; }

};  // End FDDS_Dispersion
}  // namespace sapt
}  // namespace psi

#endif
