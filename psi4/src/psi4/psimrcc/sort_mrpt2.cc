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

/**
 *  @file ccsort_mrpt2.cpp
 *  @ingroup (PSIMRCC)
*/

#include "psi4/libmoinfo/libmoinfo.h"

#include "blas.h"
#include "matrix.h"
#include "sort.h"
#include "transform.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

using namespace std;

void CCSort::build_integrals_mrpt2(IntegralTransform *ints)
{
  trans->read_integrals_mrpt2(ints);
  frozen_core_energy_mrpt2();
  allocate_and_sort_integrals_mrpt2();
  trans->free_memory();
  allocate_amplitudes_mrpt2();
}

void CCSort::frozen_core_energy_mrpt2()
{
  // One electron contribution to the frozen core energy from each irrep
  efzc=0.0;
  for(int i=0;i<nfzc;i++){
    int ii=frozen_core[i];
    efzc+=2.0*trans->oei(ii,ii);
  }
  // Two electron contribution to the frozen core energy
  for(int i=0;i<nfzc;i++){
    for(int j=0;j<nfzc;j++){
      int ii=frozen_core[i];
      int jj=frozen_core[j];
      efzc+=2.0*trans->tei_mrpt2(ii,ii,jj,jj);
      efzc-=trans->tei_mrpt2(ii,jj,ii,jj);
    }
  }
}

void CCSort::allocate_and_sort_integrals_mrpt2()
{
  // Sort the TEI for CC computations
  MatrixMap matrix_map = blas->get_MatrixMap();
  for(MatrixMap::iterator iter = matrix_map.begin(); iter!=matrix_map.end(); ++iter){
    if(iter->second->is_integral() || iter->second->is_fock()){
      iter->second->allocate_memory();
      form_fock_mrpt2(iter);
      form_two_electron_integrals_mrpt2(iter);
    }
  }
}

void CCSort::allocate_amplitudes_mrpt2()
{
  // Sort the TEI for CC computations
  MatrixMap matrix_map = blas->get_MatrixMap();
  for(MatrixMap::iterator iter = matrix_map.begin(); iter!=matrix_map.end(); ++iter){
    if(!(iter->second->is_integral() || iter->second->is_fock() )){
      iter->second->allocate_memory();
    }
  }
}

void CCSort::form_fock_mrpt2(MatrixMap::iterator& iter)
{
  CCMatrix* Matrix = iter->second;
  if(Matrix->is_fock()){
    string label     = Matrix->get_label();
    double*** matrix = Matrix->get_matrix();
    short* pq = new short[2];
    const intvec& oa2p = moinfo->get_occ_to_mo();

    bool alpha = true;
    if((label.find("O")!=string::npos) || (label.find("V")!=string::npos) || (label.find("A")!=string::npos) || (label.find("F")!=string::npos)) // NB This was missing the last bit, this might be a problem
      alpha = false;

    // N.B. Never introduce Matrices/Vectors with O or V in the name before you compute the Fock matrix elements
    vector<int> aocc = moinfo->get_aocc(Matrix->get_reference(),AllRefs);
    vector<int> bocc = moinfo->get_bocc(Matrix->get_reference(),AllRefs);
    for(int h=0;h<moinfo->get_nirreps();h++){
      for(int i = 0;i<Matrix->get_left_pairpi(h);i++){
        for(int j = 0;j<Matrix->get_right_pairpi(h);j++){
          // Find p and q from the pairs
          Matrix->get_two_indices_pitzer(pq,h,i,j);
          // Add the h(p,q) contribution
          matrix[h][i][j]=trans->oei(pq[0],pq[1]);

          // Add the core contribution//
          for(int k=0;k<nfzc;k++){
            int kk=frozen_core[k];
            matrix[h][i][j]+=add_fock_two_mrpt2(pq[0],pq[1],kk,true);
            matrix[h][i][j]+=add_fock_two_mrpt2(pq[0],pq[1],kk,false);
          }
          for(int k=0;k<aocc.size();k++){
            int kk=oa2p[aocc[k]];
            if(alpha)
              matrix[h][i][j]+=add_fock_two_mrpt2(pq[0],pq[1],kk,true);
            else
              matrix[h][i][j]+=add_fock_two_mrpt2(pq[0],pq[1],kk,false);
          }
          for(int k=0;k<bocc.size();k++){
            int kk=oa2p[bocc[k]];
            if(!alpha)
              matrix[h][i][j]+=add_fock_two_mrpt2(pq[0],pq[1],kk,true);
            else
              matrix[h][i][j]+=add_fock_two_mrpt2(pq[0],pq[1],kk,false);
          }
        }
      }
    }
    delete[] pq;
  }
//   CCMatrix* Matrix = iter->second;
//   if(Matrix->is_fock()){
//     string label     = Matrix->get_label();
//     double*** matrix = Matrix->get_matrix();
//     short* pq = new short[2];
//     int* oa2p = moinfo->get_occ_to_pitzer();
//
//     bool alpha = true;
//     if((label.find("O")!=string::npos) || (label.find("V")!=string::npos))
//       alpha = false;
//
//     // N.B. Never introduce Matrices/Vectors with O or V in the name before you compute the Fock matrix elements
//     vector<int> aocc = moinfo->get_aocc("a",Matrix->get_reference());
//     vector<int> bocc = moinfo->get_bocc("a",Matrix->get_reference());
//
//     for(int n=0;n<moinfo->get_nirreps();n++)
//       for(int i = 0;i<Matrix->get_left_pairpi(n);i++)
//         for(int j = 0;j<Matrix->get_right_pairpi(n);j++){
//           // Find p and q from the pairs
//           Matrix->get_two_indices_pitzer(pq,n,i,j);
//           // Add the h(p,q) contribution
//           matrix[n][i][j]=trans->oei(pq[0],pq[1]);
//           // Add the core contribution//
//           for(int k=0;k<nfzc;k++){
//             int kk=frozen_core[k];
//             matrix[n][i][j]+=add_fock_two_mrpt2(pq[0],pq[1],kk,true);
//             matrix[n][i][j]+=add_fock_two_mrpt2(pq[0],pq[1],kk,false);
//           }
//           for(int k=0;k<aocc.size();k++){
//             int kk=oa2p[aocc[k]];
//             if(alpha)
//               matrix[n][i][j]+=add_fock_two_mrpt2(pq[0],pq[1],kk,true);
//             else
//               matrix[n][i][j]+=add_fock_two_mrpt2(pq[0],pq[1],kk,false);
//           }
//           for(int k=0;k<bocc.size();k++){
//             int kk=oa2p[bocc[k]];
//             if(!alpha)
//               matrix[n][i][j]+=add_fock_two_mrpt2(pq[0],pq[1],kk,true);
//             else
//               matrix[n][i][j]+=add_fock_two_mrpt2(pq[0],pq[1],kk,false);
//           }
//         }
//     delete[] pq;
//   }
}

double CCSort::add_fock_two_mrpt2(int p, int q, int k, bool exchange)
{
  // Add the (pq|kk) contribution
  double term = trans->tei_mrpt2(p,q,k,k);
  // Add the -(pk|qk) contribution
  if(exchange)
    term -= trans->tei_mrpt2(p,k,q,k);
  return(term);
}

void CCSort::form_two_electron_integrals_mrpt2(MatrixMap::iterator& iter)
{
  CCMatrix* Matrix = iter->second;
  if(Matrix->is_integral()){
    short*      pqrs = new short[4];
    double*** matrix = Matrix->get_matrix();
    bool antisymmetric = Matrix->is_antisymmetric();
    if(Matrix->is_chemist()){
      for(int n=0;n<moinfo->get_nirreps();n++)
        for(int i = 0;i<Matrix->get_left_pairpi(n);i++)
          for(int j = 0;j<Matrix->get_right_pairpi(n);j++){
            Matrix->get_four_indices_pitzer(pqrs,n,i,j);
            // From (pq|rs) = <pr|qs> we define
            // (pq:rs) = <pr:qs> = (pq|rs) - (ps|qr)

            // Add the +<pr|qs> = (pq|rs) contribution
            matrix[n][i][j] += trans->tei_mrpt2(pqrs[0],pqrs[1],pqrs[2],pqrs[3]);

            // Add the -<pq|sr> = -(ps|qr) contribution
            if(antisymmetric)
              matrix[n][i][j] -= trans->tei_mrpt2(pqrs[0],pqrs[3],pqrs[1],pqrs[2]);
          }
    }else{
      for(int n=0;n<moinfo->get_nirreps();n++)
        for(int i = 0;i<Matrix->get_left_pairpi(n);i++)
          for(int j = 0;j<Matrix->get_right_pairpi(n);j++){
            Matrix->get_four_indices_pitzer(pqrs,n,i,j);
            // Add the +<pq|rs> = (pr|qs) contribution
            matrix[n][i][j] += trans->tei_mrpt2(pqrs[0],pqrs[2],pqrs[1],pqrs[3]);

            // Add the -<pq|sr> = -(ps|qr) contribution
            if(antisymmetric)
              matrix[n][i][j] -= trans->tei_mrpt2(pqrs[0],pqrs[3],pqrs[1],pqrs[2]);
          }
    }
    delete[] pqrs;
  }
}

}} /* End Namespaces */
