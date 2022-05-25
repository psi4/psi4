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

#ifndef COMPOSITE_H
#define COMPOSITE_H

#include "psi4/libfock/jk.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/potential.h"
#include "psi4/libmints/twobody.h"
#include "psi4/lib3index/dfhelper.h"

#include <unordered_set>

namespace psi {

/**
 * Class DirectDFJ
 *
 * Builds the J matrix using a direct density-fitted algorithm
 */
class DirectDFJ : public SplitJKBase {
  protected:
   /// Auxiliary basis set
   std::shared_ptr<BasisSet> auxiliary_;
   /// Coulomb Metric
   SharedMatrix J_metric_;
   /// Numerical cutoff for ERI screening
   double cutoff_;

   /// Form J_PQ
   void build_metric();

   /// Builds the integrals for the DirectDFJ class
   void build_ints() override;

  public:
   /**
    * @brief Construct a new DirectDFJ object
    * 
    * @param primary The primary basis set used in DirectDFJ
    * @param auxiliary The auxiliary basis set used in DirectDFJ
    * @param options The options object
    */
   DirectDFJ(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options);

   /**
    * @author Zach Glick, Georgia Tech, May 2022
    *
    * @brief Builds the J matrix according to the DirectDFJ algorithm, described in [Weigand:2002:4285]_
    * doi: 10.1039/b204199p
    * 
    * @param D The list of AO density matrixes to contract to form the J matrix (1 for RHF, 2 for UHF/ROHF)
    * @param J The list of AO J matrices to build (Same size as D)
    */
   void build_G_component(const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J) override;

   /**
    * @brief Prints information regarding Direct-DF-J run
    * 
    */
   void print_header() override;
};

/**
 * Class LinK
 *
 * Builds the K matrix using Ochsenfeld's LinK algorithm
 */
class LinK : public SplitJKBase {
  protected:
   /// ERI Screening Cutoff
   double cutoff_;
   /// Density-based Sparsity Screening Cutoff for LinK
   double linK_ints_cutoff_;

   /// Builds the integrals for the LinK class
   void build_ints() override;

  public:
   /**
    * @brief Construct a new LinK object
    * 
    * @param primary The primary basis set used in LinK
    * @param options The options object
    */
   LinK(std::shared_ptr<BasisSet> primary, Options& options);

   /**
    * @author Andy Jiang, Georgia Tech, March 2022
    *
    * @brief Builds the K matrix according to the LinK algorithm, described in [Ochsenfeld:1998:1663]_
    * doi: 10.1063/1.476741
    * 
    * @param D The list of AO density matrixes to contract to form the K matrix (1 for RHF, 2 for UHF/ROHF)
    * @param K The list of AO K matrices to build (Same size as D)
    */
   void build_G_component(const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& K) override;

   /**
    * @brief Prints information regarding LinK run
    * 
    */
   void print_header() override;

};

/**
 * Class COSK
 *
 * Builds the K matrix using a seminumerical algorithm
 */
class COSK : public SplitJKBase {
  protected:
   /// Small DFTGrid for initial SCF iterations
   std::shared_ptr<DFTGrid> grid_init_;
   /// Large DFTGrid for the final SCF iteration
   std::shared_ptr<DFTGrid> grid_final_;
   /// Overlap fitting metric for grid_initial_
   SharedMatrix Q_init_;
   /// Overlap fitting metric for grid_final_
   SharedMatrix Q_final_;

   /// Builds grids and integrals for the COSK class
   void build_ints() override;

  public:
   /**
    * @brief Construct a new COSK object
    * 
    * @param primary The primary basis set used in COSK
    * @param options The options object
    */
   COSK(std::shared_ptr<BasisSet> primary, Options& options);

   /**
    * @author Zach Glick, Georgia Tech, May 2022
    *
    * @brief Builds the K matrix according to the COSK algorithm, described in [Neese:2009:98]_
    * doi: 10.1063/1.476741
    * 
    * @param D The list of AO density matrixes to contract to form the K matrix (1 for RHF, 2 for UHF/ROHF)
    * @param K The list of AO K matrices to build (Same size as D)
    */
   void build_G_component(const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& K) override;

   /**
    * @brief Prints information regarding COSK run
    */
   void print_header() override;

};

}

#endif