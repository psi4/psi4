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

    // Coulomb metric -0.5 power
    SharedMatrix metric_half_inv_;

    // Auxiliary overlap matrix
    SharedMatrix aux_overlap_;

    // DFHelper object
    std::shared_ptr<DFHelper> dfh_;

    // Cache map
    std::map<std::string, SharedMatrix> matrix_cache_;
    std::map<std::string, SharedVector> vector_cache_;

    // Is hybrid functional? 
    bool is_hybrid_;

    // QR factorization result
    SharedMatrix R_A_, R_B_;

   public:
    /**
     * Constructs the FDDS_Dispersion object.
     * @param primary   The primary basis
     * @param auxiliary The auxiliary basis
     * @param cache     A data cache containing "Cocc_A", "Cvir_A", "eps_occ_A", "eps_vir_A", "Cocc_B",
     * "Cvir_B", "eps_occ_B", "eps_vir_B" quantities
     * @param is_hybrid Flag of hybrid functional
     */
    FDDS_Dispersion(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                    std::map<std::string, SharedMatrix> matrix_cache, std::map<std::string, SharedVector> vector_cache,
                    bool is_hybrid);

    ~FDDS_Dispersion();

    /**
     * Projects the densities from primary AO space to auxiliary AO space
     * @param  dens     Vector of densities to transform
     * @return aux_dens Vector of transformed densities
     */
    std::vector<SharedMatrix> project_densities(std::vector<SharedMatrix> dens);

    /**
     * Forms the uncoupled amplitude Qia,Eia,Pia->PQ and projector Qia,Pia->PQ
     * @param  monomer Monomer "A" or "B"
     * @param  omega   Time dependent value
     * @return         "PQ" amplitude 
     */
    SharedMatrix form_unc_amplitude(std::string monomer, double omega);

    /**
     * Forms the uncoupled amplitude and other PQ matrices in hybrid FDDS dispersion
     * @param  monomer Monomer "A" or "B"
     * @param  omega   Time dependent value
     * @return         Dictionary of PQ matrics 
     * @return ret["amp"]  PQ uncoupled amplitude
     * @return ret["K1LD"]  K1(lambda * d) = (P|ar) LDar (ar|X+Y|Q)
     * @return ret["K2LD"]  K2(lambda * d) = (P|ar) LDar (ar|X-Y|Q)
     * @return ret["K2L"]  K2(lambda) = (P|ar) Lar (ar|X-Y|Q)
     * @return ret["K21L"]  K21(lambda) = (P|X-Y|ar) Lar (ar|X+Y|Q)
     */
    std::map<std::string, SharedMatrix> form_aux_matrices(std::string monomer, double omega);

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

    /**
     * Returns R for QR decomposition of (ar|Q)
     * @return R_A
     */
    SharedMatrix R_A() { return R_A_; }

    /**
     * Returns R for QR decomposition of (bs|Q)
     * @return R_B
     */
    SharedMatrix R_B() { return R_B_; }

    /**
     * Forms X-type 3-index exchange integral tensor (ar|X|Q) = (ar'|a'r)(a'r'|Q)
     * @param monomer Monomer "A" or "B"
     */
    void form_X(std::string monomer);

    /**
     * Forms Y-type 3-index exchange integral tensor (ar|Y|Q) = (aa'|rr')(a'r'|Q)
     * @param monomer Monomer "A" or "B"
     */
    void form_Y(std::string monomer);

    /**
     * Performs QR factorization and store (R^t)^-1 into R_A_ or R_B_
     * @param monomer Monomer "A" or "B"
     * @return (R^t)^-1
     */
    SharedMatrix QR(std::string monomer);

    SharedMatrix get_tensor_pqQ(std::string name, std::tuple<size_t, size_t, size_t> dimensions);
    void print_tensor_pqQ(std::string tensor_name, std::string file_name, std::tuple<size_t, size_t, size_t> dimensions);

};  // End FDDS_Dispersion
}  // namespace sapt
}  // namespace psi

#endif
