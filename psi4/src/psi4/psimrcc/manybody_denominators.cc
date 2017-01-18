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
 *  @file manybody_denominators.cc
 *  @ingroup (PSIMRCC)
 *  @brief The base class for all the many-body methods
*/

#include <iostream>
#include <cmath>
#include <algorithm>

#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libmoinfo/libmoinfo.h"

#include "algebra_interface.h"
#include "blas.h"
#include "debugging.h"
#include "manybody.h"
#include "matrix.h"
#include "sort.h"


namespace psi{

    namespace psimrcc{
    extern MOInfo *moinfo;

using namespace std;

/**
 * Generates the MP denominators
 * \f[ \Delta_{ij...}^{ab...}(\mu) = f_{ii}(\mu) + f_{jj}(\mu) + ... - f_{aa}(\mu) - f_{bb}(\mu) - ... \f]
 * where the excitations that are not allowed in reference \f$ \mu \f$ are set to a large value (see huge)
 */
void CCManyBody::generate_denominators()
{
  START_TIMER(1,"Generating Denominators");

  bool keep_denominators_in_core = false;
  if(options_.get_str("CORR_WFN")=="PT2")
    keep_denominators_in_core = true;

  MatrixMap& matrix_map = blas->get_MatrixMap();
  MatrixMap::iterator end_iter = matrix_map.end();
  for(MatrixMap::iterator iter=matrix_map.begin();iter!=end_iter;++iter){
    string str = iter->first;

    if(str.find("d1")!=string::npos){

      // Load the temporary matrix
      CCMatTmp MatTmp  = blas->get_MatTmp(str,( keep_denominators_in_core ? none : dump));
      double*** matrix = MatTmp->get_matrix();

      // Get the reference number
      int reference = MatTmp->get_reference();

      vector<string> spl_str = split_indices(str);
      string ind;
      for(size_t i = 0; i < spl_str.size(); ++i){
        ind += spl_str[i];
      }

      std::vector<int> aocc     = moinfo->get_aocc(reference,AllRefs);
      std::vector<int> bocc     = moinfo->get_bocc(reference,AllRefs);
      std::vector<int> avir     = moinfo->get_avir(reference,AllRefs);
      std::vector<int> bvir     = moinfo->get_bvir(reference,AllRefs);


      // Build the is_arrays for reference ref
      std::vector<bool> is_aocc(moinfo->get_nocc(),false);
      std::vector<bool> is_bocc(moinfo->get_nocc(),false);
      std::vector<bool> is_avir(moinfo->get_nvir(),false);
      std::vector<bool> is_bvir(moinfo->get_nvir(),false);
      std::vector<bool> is_frzv(moinfo->get_nfvir(),true);
      std::vector<bool> is_element[2];

      for(size_t i = 0; i < aocc.size(); i++) is_aocc[aocc[i]]=true;
      for(size_t i = 0; i < bocc.size(); i++) is_bocc[bocc[i]]=true;
      for(size_t i = 0; i < avir.size(); i++) is_avir[avir[i]]=true;
      for(size_t i = 0; i < bvir.size(); i++) is_bvir[bvir[i]]=true;


      //      // Read the Fock matrices
      CCMatTmp f_oo_Matrix = blas->get_MatTmp("fock[oo]",reference,( keep_denominators_in_core ? none : dump));
      CCMatTmp f_OO_Matrix = blas->get_MatTmp("fock[OO]",reference,( keep_denominators_in_core ? none : dump));
      CCMatTmp f_vv_Matrix = blas->get_MatTmp("fock[vv]",reference,( keep_denominators_in_core ? none : dump));
      CCMatTmp f_VV_Matrix = blas->get_MatTmp("fock[VV]",reference,( keep_denominators_in_core ? none : dump));
      CCMatTmp f_ff_Matrix = blas->get_MatTmp("fock[ff]",reference,( keep_denominators_in_core ? none : dump));
      CCMatTmp f_FF_Matrix = blas->get_MatTmp("fock[FF]",reference,( keep_denominators_in_core ? none : dump));

      vector<CCMatrix*> f_Matrix;
      int k = 0;
      for(size_t l = 0; l < ind.size(); ++l){
        switch(ind[l]){
        case 'o':

          f_Matrix.push_back(f_oo_Matrix.get_CCMatrix());
          is_element[k] = is_aocc;
          k++;
          break;
        case 'O':
          f_Matrix.push_back(f_OO_Matrix.get_CCMatrix());
          is_element[k] = is_bocc;
          k++;
          break;
        case 'v':
          f_Matrix.push_back(f_vv_Matrix.get_CCMatrix());
          is_element[k] = is_avir;
          k++;
          break;
        case 'V':
          f_Matrix.push_back(f_VV_Matrix.get_CCMatrix());
          is_element[k] = is_bvir;
          k++;
          break;
        case 'f':
          f_Matrix.push_back(f_ff_Matrix.get_CCMatrix());
          is_element[k] = is_frzv;
          k++;
          break;
        case 'F':
          f_Matrix.push_back(f_FF_Matrix.get_CCMatrix());
          is_element[k] = is_frzv;
          k++;
          break;
        default:
          break;
        }
      }

      // N.B. Never introduce Matrices/Vectors with O or V in the name before you compute the Fock matrix elements

      short* ia = new short[2];
      for(int n=0;n<moinfo->get_nirreps();n++)
        for(size_t i = 0;i<MatTmp->get_left_pairpi(n);i++)
          for(size_t j = 0;j<MatTmp->get_right_pairpi(n);j++){
            // Set the denomiator to huge by default
            matrix[n][i][j]=huge;

            MatTmp->get_two_indices(ia,n,i,j);

            // Build the denominator
            if(is_element[0][ia[0]] && is_element[1][ia[1]]){
              matrix[n][i][j]  = f_Matrix[0]->get_two_address_element(ia[0],ia[0]);
              matrix[n][i][j] -= f_Matrix[1]->get_two_address_element(ia[1],ia[1]);
            }

          } // End loop over n,i,j
      delete[] ia;
    } // End "if d1"

    if(str.find("d2")!=string::npos){
      // Load the temporary matrix
      CCMatTmp MatTmp  = blas->get_MatTmp(str,( keep_denominators_in_core ? none : dump));
      double*** matrix = MatTmp->get_matrix();

      // Get the reference number
      int reference = MatTmp->get_reference();

      vector<string> spl_str = split_indices(str);
      string ind;
      for(size_t i = 0; i < spl_str.size(); ++i){
        ind += spl_str[i];
      }

      std::vector<int> aocc     = moinfo->get_aocc(reference,AllRefs);
      std::vector<int> bocc     = moinfo->get_bocc(reference,AllRefs);
      std::vector<int> avir     = moinfo->get_avir(reference,AllRefs);
      std::vector<int> bvir     = moinfo->get_bvir(reference,AllRefs);


      // Build the is_arrays for reference ref
      std::vector<bool> is_aocc(moinfo->get_nocc(),false);
      std::vector<bool> is_bocc(moinfo->get_nocc(),false);
      std::vector<bool> is_avir(moinfo->get_nvir(),false);
      std::vector<bool> is_bvir(moinfo->get_nvir(),false);
      std::vector<bool> is_frzv(moinfo->get_nfvir(),true);
      std::vector<bool> is_element[4];

      for(size_t i = 0; i < aocc.size(); i++) is_aocc[aocc[i]]=true;
      for(size_t i = 0; i < bocc.size(); i++) is_bocc[bocc[i]]=true;
      for(size_t i = 0; i < avir.size(); i++) is_avir[avir[i]]=true;
      for(size_t i = 0; i < bvir.size(); i++) is_bvir[bvir[i]]=true;


      //      // Read the Fock matrices
      CCMatTmp f_oo_Matrix = blas->get_MatTmp("fock[oo]",reference,( keep_denominators_in_core ? none : dump));
      CCMatTmp f_OO_Matrix = blas->get_MatTmp("fock[OO]",reference,( keep_denominators_in_core ? none : dump));
      CCMatTmp f_vv_Matrix = blas->get_MatTmp("fock[vv]",reference,( keep_denominators_in_core ? none : dump));
      CCMatTmp f_VV_Matrix = blas->get_MatTmp("fock[VV]",reference,( keep_denominators_in_core ? none : dump));
      CCMatTmp f_ff_Matrix = blas->get_MatTmp("fock[ff]",reference,( keep_denominators_in_core ? none : dump));
      CCMatTmp f_FF_Matrix = blas->get_MatTmp("fock[FF]",reference,( keep_denominators_in_core ? none : dump));

      vector<CCMatrix*> f_Matrix;
      int k = 0;
      for(size_t l = 0; l < ind.size(); ++l){
        switch(ind[l]){
        case 'o':

          f_Matrix.push_back(f_oo_Matrix.get_CCMatrix());
          is_element[k] = is_aocc;
          k++;
          break;
        case 'O':
          f_Matrix.push_back(f_OO_Matrix.get_CCMatrix());
          is_element[k] = is_bocc;
          k++;
          break;
        case 'v':
          f_Matrix.push_back(f_vv_Matrix.get_CCMatrix());
          is_element[k] = is_avir;
          k++;
          break;
        case 'V':
          f_Matrix.push_back(f_VV_Matrix.get_CCMatrix());
          is_element[k] = is_bvir;
          k++;
          break;
        case 'f':
          f_Matrix.push_back(f_ff_Matrix.get_CCMatrix());
          is_element[k] = is_frzv;
          k++;
          break;
        case 'F':
          f_Matrix.push_back(f_FF_Matrix.get_CCMatrix());
          is_element[k] = is_frzv;
          k++;
          break;
        default:
          break;
        }
      }

      // N.B. Never introduce Matrices/Vectors with O or V in the name before you compute the Fock matrix elements

      short* ijab = new short[4];
      for(int n=0;n<moinfo->get_nirreps();n++)
        for(size_t i = 0;i<MatTmp->get_left_pairpi(n);i++)
          for(size_t j = 0;j<MatTmp->get_right_pairpi(n);j++){
            // Set the denomiator to huge by default
            matrix[n][i][j]=huge;

            MatTmp->get_four_indices(ijab,n,i,j);

            // Build the denominator
            if(is_element[0][ijab[0]] && is_element[1][ijab[1]] && is_element[2][ijab[2]] && is_element[3][ijab[3]]){
              matrix[n][i][j]  = f_Matrix[0]->get_two_address_element(ijab[0],ijab[0]);
              matrix[n][i][j] += f_Matrix[1]->get_two_address_element(ijab[1],ijab[1]);
              matrix[n][i][j] -= f_Matrix[2]->get_two_address_element(ijab[2],ijab[2]);
              matrix[n][i][j] -= f_Matrix[3]->get_two_address_element(ijab[3],ijab[3]);
            }

          } // End loop over n,i,j
      delete[] ijab;
    } // End "if d2"

  } // End for each reference
  END_TIMER(1);
}

}} /* End Namespaces */
