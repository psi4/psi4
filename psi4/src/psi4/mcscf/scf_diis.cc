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
#include <cstdio>

#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "scf.h"
#include "algebra_interface.h"

extern FILE* outfile;

namespace psi{ namespace mcscf{

extern MemoryManager* memory_manager;

void SCF::diis(int cycle)
{
  // Transform Feff from MO to orthogonal AO
  SBlockMatrix CFeffC("CFeffC",nirreps,sopi,sopi);

  transform(Feff_t,CFeffC,C_T);
  transform(CFeffC,Feff_oAO,S_sqrt);

  diis_F[current_diis] = Feff_oAO;

  // Build and transform the error vector from MO to orthogonal AO
  SBlockMatrix CeC("CeC",nirreps,sopi,sopi);
  SBlockMatrix e_oAO("e_oAO",nirreps,sopi,sopi);

  e = Feff_t;
  e->zero_diagonal();
  transform(e,CeC,C_T);
  transform(CeC,e_oAO,S_sqrt);

  diis_e[current_diis] = e_oAO;

  if(reference == tcscf){
    for(int I = 0 ; I < nci; ++I){
      diis_ci[I][current_diis] = ci[I];
    }
  }

  outfile->Printf(" S");

  if(cycle >= ndiis){
    outfile->Printf("/E");
    int  matrix_size = ndiis + 1;

    double** diis_B;
    double*  diis_A;
    allocate1(double,diis_A,matrix_size);
    allocate2(double,diis_B,matrix_size,matrix_size);

    // Zero A and B
    for(int i=0;i<ndiis;i++){
      diis_A[i]=0.0;
      diis_B[i][ndiis]=diis_B[ndiis][i]=-1.0;
      for(int j=0;j<ndiis;j++)
        diis_B[i][j]=0.0;
    }
    diis_B[ndiis][ndiis]=0.0;
    diis_A[ndiis]=-1.0;

    // Build the diis_B matrix
    for(int i=0;i<ndiis;i++){
      for(int j=i;j<ndiis;j++){
        diis_B[i][j] += dot(diis_e[i],diis_e[j]);
        diis_B[j][i] = diis_B[i][j];
      }
    }

    // Solve B x = A
    int* IPIV = new int[matrix_size];
    int nrhs = 1;
    int info = 0;
    F_DGESV(&matrix_size, &nrhs, &(diis_B[0][0]),&matrix_size, &(IPIV[0]), &(diis_A[0]),&matrix_size, &info);
    delete[] IPIV;

    // Update F = sum diis_F(i) * A(i);
    if(!info){
      Feff_oAO->zero();
      for(int i=0;i<ndiis;i++){
        e = diis_F[i];
        e->scale(diis_A[i]);
        Feff_oAO += e;
      }

      if(reference == tcscf && options_.get_bool("CI_DIIS") ){
        for(int I = 0 ; I < nci; ++I){
          ci[I] = 0.0;
          for(int i=0; i < ndiis;i++){
            ci[I] += diis_ci[I][i] * diis_A[i];
          }
        }
      }

    }else{
      outfile->Printf(" (singularities found)");
    }

    release1(diis_A);
    release2(diis_B);
  }
  current_diis++;
  if(current_diis == ndiis)
    current_diis  = 0;
}

}} /* End Namespaces */
