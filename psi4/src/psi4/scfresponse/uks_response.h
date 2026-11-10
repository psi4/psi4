/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#ifndef UKS_RESPONSE_H
#define UKS_RESPONSE_H

// for dft
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"


#include "ks_response.h" 

namespace psi {namespace scf_response{ 

class UKSResponseSolver: public KSResponseSolver {

  public:

    UKSResponseSolver(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_);

    ~UKSResponseSolver();

  protected:

    std::vector<std::vector<double>> first_order_response(
        std::vector<std::shared_ptr<Matrix>> op_a,
        std::vector<std::shared_ptr<Matrix>> op_b,
        double omega
    );

    void compute_polarizability(
        std::vector<double>X,
        std::vector<double>Y,
        double omega
    );

    void compute_hyperpolarizability(std::vector<std::vector<double>>amps_wx,
        std::vector<std::vector<double>>amps_wy,
        std::vector<std::vector<double>>amps_wz,
        std::string type, double omega
    );

    void build_Au_Bu(int N, int L, double *u, double *ABu);

};

}} // End namespaces

#endif 
