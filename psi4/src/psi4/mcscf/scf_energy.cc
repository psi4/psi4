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

#include <iostream>
#include <cmath>

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "scf.h"
#include "sblock_matrix.h"

extern FILE* outfile;

namespace psi{ namespace mcscf{

extern MemoryManager* memory_manager;

double SCF::energy(int cycle,double old_energy)
{
  double electronic_energy = 0.0;

  T  = H;
  T += Fc;
  electronic_energy += dot(Dc,T);

  if(reference == rohf){
    T  = H;
    T.scale(0.5);
    T += Fo;
    electronic_energy += dot(Do,T);
  }

  total_energy = electronic_energy + moinfo_scf->get_nuclear_energy();

  if(reference == tcscf){
//     SBlockMatrix Dtc_sum("Dtc sum",nirreps,sopi,sopi);
//
//     // Compute diagonal elements of H
//     for(int I = 0 ; I < nci; ++I){
//       Dtc_sum  = Dc;
//       Dtc_sum += Dtc[I];
//       construct_G(Dtc_sum,G,"PK");
//       T  = H;
//       T.scale(2.0);
//       T += G;
//       H_tcscf[I][I] = dot(Dtc_sum,T) + moinfo_scf->get_nuclear_energy();
//     }
//
//     // Compute off-diagonal elements of H
//     for(int I = 0 ; I < nci; ++I){
//       for(int J = I + 1; J < nci; ++J){
//         construct_G(Dtc[I],G,"K");
//         H_tcscf[I][J] = H_tcscf[J][I] = - dot(Dtc[J],G);
//       }
//     }

//     outfile->Printf("\n  Hamiltonian");
//     for(int I = 0 ; I < nci; ++I){
//       outfile->Printf("\n    ");
//       for(int J = 0 ; J < nci; ++J)
//         outfile->Printf(" %11.8f ",H_tcscf[I][J]);
//     }


    // Compute the CI gradient
    norm_ci_grad = 0.0;
    for(int I = 0 ; I < nci; ++I){
      ci_grad[I] = 0.0;
      for(int J = 0; J < nci; ++J){
        ci_grad[I] += H_tcscf[I][J] * ci[J];
      }
      ci_grad[I] -= old_energy * ci[I];
      norm_ci_grad += fabs(ci_grad[I]);
    }

    double*  eigenvalues;
    double** eigenvectors;
    allocate1(double,eigenvalues,nci);
    allocate2(double,eigenvectors,nci,nci);

    sq_rsp(nci,nci,H_tcscf,eigenvalues,1,eigenvectors,1.0e-14);

    total_energy = eigenvalues[root];

    if(fabs(old_energy - total_energy) < 1.0e-5){
    for(int I = 0 ; I < nci; ++I)
      ci[I] = eigenvectors[I][root];
    }
    release1(eigenvalues);
    release2(eigenvectors);
  }

  return(total_energy);
}

}} /* End Namespaces */
