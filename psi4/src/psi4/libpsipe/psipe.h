/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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
 
#ifndef PSIPE_H
#define PSIPE_H
#ifdef USING_cppe

#include <string>
#include <memory>

#include "psi4/libmints/dimension.h"
#include "psi4/libmints/typedefs.h"

#include <cppe/core/multipole.hh>
#include <cppe/core/molecule.hh>
#include <cppe/core/cppe_state.hh>
#include <cppe/core/pe_options.hh>


namespace psi {
  class BasisSet;
  class OneBodyAOInt;
  
  class PeIntegralHelper {
      public:
          PeIntegralHelper(std::shared_ptr<BasisSet> basisset) : basisset_(basisset) {}
          /*! \brief computes the potential integrals at a site through k-th order
           *   
           */
          SharedMatrix compute_multipole_potential_integrals(Vector3 site, int order, std::vector<double>& moments);
          
          /*! \brief computes the field integrals at a site
           *
           */
          SharedMatrix compute_field_integrals(Vector3 site, arma::vec moments);
          
          /*! \brief computes the field at a site using the density matrix
           *
           */
          Vector compute_field(Vector3 site, const SharedMatrix &D);
          
          
      private:
          std::shared_ptr<BasisSet> basisset_;
  };

  class PeState {

  public:
    enum CalcType { total, electronic_only };
    PeState() = default;
    PeState(libcppe::PeOptions options, std::shared_ptr<BasisSet> basisset);
    ~PeState() {}
    
    /*! \brief Compute PE energy and Fock matrix contribution
     *  \param[in] D density matrix
     *  \param[in] type (total Fock contribution or electronic contribution only)
     */
    std::pair<double, SharedMatrix> compute_pe_contribution(const SharedMatrix &D, CalcType type = CalcType::total,
                                                            bool subtract_scf_density = false);

    /* \brief prints the summary table of PE energy contributions to the Psi4 output file
    */
    void print_energy_summary();

    // void calculate_fock_contribution(arma::mat& Ptot, const libqints::multi_array<double>& out, double* energy);
    // void calculate_excited_state_energy_correction(double* p_exc, double* energy, bool is_tdm);

  private:
      std::shared_ptr<BasisSet> basisset_;
      std::vector<libcppe::Potential> potentials_;
      libcppe::CppeState cppe_state_;
      PeIntegralHelper int_helper_;
      
      SharedMatrix V_es_;
      SharedMatrix D_scf_;
      
      int iteration = 0;
  };


} // namespace psi

#endif // USING_cppe
#endif // PSIPE_H