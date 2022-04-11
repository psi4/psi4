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

class DirectDFJ : public JBase {
  protected:
   /// The auxiliary basis set used in the DF algorithm
   std::shared_ptr<BasisSet> auxiliary_;
   /// The metrix used in the DF algorithm J_PQ = (P|Q)
   SharedMatrix Jmet_;
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
    * @author Andy Jiang, Georgia Tech, April 2022
    *
    * @brief Builds the J matrix according to the DirectDFJ algorithm, described in [Weigand:2002:4285]_
    * doi: 10.1039/b204199p
    * 
    * @param D The list of AO density matrixes to contract to form the J matrix (1 for RHF, 2 for UHF/ROHF)
    * @param J The list of AO J matrices to build (Same size as D)
    */
   void build_J(const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J) override;

};

class LinK : public KBase {
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
   void build_K(const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& K) override;

};

}

#endif