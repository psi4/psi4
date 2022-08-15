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

#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "scf.h"

namespace psi {
namespace mcscf {

extern MemoryManager* memory_manager;

double SCF::energy(int cycle, double old_energy) {
    double electronic_energy = 0.0;

    T = H;
    T += Fc;
    electronic_energy += dot(Dc, T);

    if (reference == rohf) {
        T = H;
        T.scale(0.5);
        T += Fo;
        electronic_energy += dot(Do, T);
    }

    total_energy = electronic_energy + moinfo_scf->get_nuclear_energy();

    if (reference == tcscf) {
        // Compute the CI gradient
        norm_ci_grad = 0.0;
        for (int I = 0; I < nci; ++I) {
            ci_grad[I] = 0.0;
            for (int J = 0; J < nci; ++J) {
                ci_grad[I] += H_tcscf[I][J] * ci[J];
            }
            ci_grad[I] -= old_energy * ci[I];
            norm_ci_grad += std::fabs(ci_grad[I]);
        }

        double* eigenvalues;
        double** eigenvectors;
        allocate1(double, eigenvalues, nci);
        allocate2(double, eigenvectors, nci, nci);

        if (DSYEV_ascending(nci, H_tcscf, eigenvalues, eigenvectors) != 0){
            outfile->Printf("DSYEV failed in mcscf::SCF::energy()");
            throw PsiException("DSYEV failed in mcscf::SCF::energy()", __FILE__, __LINE__);
        }

        total_energy = eigenvalues[root];

        if (std::fabs(old_energy - total_energy) < 1.0e-5) {
            for (int I = 0; I < nci; ++I) ci[I] = eigenvectors[I][root];
        }
        release1(eigenvalues);
        release2(eigenvectors);
    }

    return (total_energy);
}

}  // namespace mcscf
}  // namespace psi
