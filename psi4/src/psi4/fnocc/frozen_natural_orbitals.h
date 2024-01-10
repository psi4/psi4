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

#ifndef FROZENNO_H
#define FROZENNO_H
#include "psi4/psi4-dec.h"
#include "psi4/libmints/wavefunction.h"

namespace psi {
namespace fnocc {

// base class
class FrozenNO : public Wavefunction {
   public:
    FrozenNO(std::shared_ptr<Wavefunction> wfn, Options& options);
    ~FrozenNO() override;

    double compute_energy() override;
    void ComputeNaturalOrbitals();

   protected:
    // mp2 energy in full basis
    double emp2;
    long int nso, nmo, ndocc, nvirt, nfzc, nfzv, ndoccact, nvirt_no;

    void common_init();
};

class PSI_API DFFrozenNO : public FrozenNO {
   public:
    DFFrozenNO(std::shared_ptr<Wavefunction> wfn, Options& options);
    ~DFFrozenNO() override;

    double compute_energy() override;

    /// computes MP2 natural orbitals
    void ComputeNaturalOrbitals();

    /// generates 3-index integrals and writes them to disk
    void ThreeIndexIntegrals();

    /// generates 4-index eri's from 3-index integrals
    void FourIndexIntegrals();

   protected:
    void ModifyCa(const std::vector<double>& Dab);
    void ModifyCa_occ(double* Dij);
    void BuildFock(long int nQ, double* Qso, double* F);
    void TransformQ(long int nQ, double* Qso);
};
}
}

#endif
