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
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libciomr/libciomr.h"

#include "blas.h"
#include "index.h"
#include "matrix.h"

namespace psi{

    namespace psimrcc{
    extern MOInfo *moinfo;

using namespace std;

void CCOperation::sort()
{
  sort(B_Matrix->get_left(),B_Matrix->get_right(),B_Matrix->get_matrix(),factor);
}

void CCOperation::sort(CCIndex* T_left,CCIndex* T_right,double*** T_matrix,double constant)
{
  typedef std::vector<size_t>  Size_tVec;
  Timer sort_timer;
  // Setup reindexing_array
  // This assumes that the reindexing starts from 1 !!! This can cost you an headache
  if(reindexing.size()>6)
    throw PSIEXCEPTION("CCOperation::sort() doesn't support more that six indices");
  short* reindexing_array = new short[6];
  for(size_t i = 0; i< reindexing.size(); i++)
    reindexing_array[i]= to_integer(reindexing.substr(i, 1))-1;

  CCIndex* A_left    = A_Matrix->get_left();
  CCIndex* A_right   = A_Matrix->get_right();
  double*** A_matrix = A_Matrix->get_matrix();

  int A_left_nelements = A_left->get_nelements();
  int T_left_nelements = T_left->get_nelements();

  Size_tVec&    A_left_pairpi   = A_left->get_pairpi();
  Size_tVec&    A_right_pairpi  = A_right->get_pairpi();
  Size_tVec&    A_left_first    = A_left->get_first();
  Size_tVec&    A_right_first   = A_right->get_first();
  short**       A_left_tuples   = A_left->get_tuples();
  short**       A_right_tuples  = A_right->get_tuples();
  Size_tVec&    T_left_first    = T_left->get_first();
  Size_tVec&    T_right_first   = T_right->get_first();
  Size_tVec&    T_left_pairpi   = T_left->get_pairpi();
  Size_tVec&    T_right_pairpi  = T_right->get_pairpi();
  short**       T_left_tuples   = T_left->get_tuples();
  short**       T_right_tuples  = T_right->get_tuples();

  // Zero the target matrix (A) if the assignment operator requires it

  if(assignment=="=" || assignment==">="){
    for(int h=0;h<moinfo->get_nirreps();++h){
      size_t block_size = A_Matrix->get_block_sizepi(h);
      if(block_size>0){
        double* A_block = &(A_matrix[h][0][0]);
        for(size_t i = 0;i<block_size;i++)
          A_block[i]=0.0;
      }
    }
  }


  // This routine performs the reindexing of CCMatrix objects
  // for matrices that have the same number of elements or
  // when size(A) < size(T)

  if(assignment=="=" || assignment=="+="){
    if(reindexing.size()==2){
      int b_irrep;
      size_t b_left,b_right;
      short* pq = new short[2];
      // A[x][x] <- T[x][x]
      if((A_left_nelements == 1) && (T_left_nelements == 1)){
        int*    T_left_one_index_to_irrep    = T_left->get_one_index_to_irrep();
        size_t* T_left_one_index_to_tuple    = T_left->get_one_index_to_tuple_rel_index();
        size_t* T_right_one_index_to_tuple   = T_right->get_one_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<A_left_pairpi[n];i++){
            pq[0]=A_left_tuples[i+A_left_first[n]][0];
            for(size_t j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pq[1]=A_right_tuples[j+A_right_first[n]][0];
              b_irrep = T_left_one_index_to_irrep[pq[reindexing_array[0]]];
              b_left  = T_left_one_index_to_tuple[pq[reindexing_array[0]]];
              b_right = T_right_one_index_to_tuple[pq[reindexing_array[1]]];
              A_matrix[n][i][j]+=constant*T_matrix[b_irrep][b_left][b_right];
            }
          }
        }
      }
      // A[x][x] <- T[xx][]
      if((A_left_nelements == 1) && (T_left_nelements == 2)){
        int**    T_left_two_index_to_irrep = T_left->get_two_index_to_irrep();
        size_t** T_left_two_index_to_tuple = T_left->get_two_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<A_left_pairpi[n];i++){
            pq[0]=A_left_tuples[i+A_left_first[n]][0];
            for(size_t j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pq[1]=A_right_tuples[j+A_right_first[n]][0];
              b_irrep = T_left_two_index_to_irrep[pq[reindexing_array[0]]][pq[reindexing_array[1]]];
              b_left  = T_left_two_index_to_tuple[pq[reindexing_array[0]]][pq[reindexing_array[1]]];
              A_matrix[n][i][j]+=constant*T_matrix[b_irrep][b_left][0];
            }
          }
        }
      }
      // A[xx][] <- T[x][x]
      if((A_left_nelements == 2) && (T_left_nelements == 1)){
        int*    T_left_one_index_to_irrep    = T_left->get_one_index_to_irrep();
        size_t* T_left_one_index_to_tuple    = T_left->get_one_index_to_tuple_rel_index();
        size_t* T_right_one_index_to_tuple   = T_right->get_one_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<A_left_pairpi[n];i++){
            pq[0]=A_left_tuples[i+A_left_first[n]][0];
            pq[1]=A_left_tuples[i+A_left_first[n]][1];
            for(size_t j = 0;j<A_right_pairpi[n];j++){
              b_irrep = T_left_one_index_to_irrep[pq[reindexing_array[0]]];
              b_left  = T_left_one_index_to_tuple[pq[reindexing_array[0]]];
              b_right = T_right_one_index_to_tuple[pq[reindexing_array[1]]];
              A_matrix[n][i][0]+=constant*T_matrix[b_irrep][b_left][b_right];
            }
          }
        }
      }
      // A[xx][] <- T[xx][]
      if((A_left_nelements == 2) && (T_left_nelements == 2)){
        int**    T_left_two_index_to_irrep    = T_left->get_two_index_to_irrep();
        size_t** T_left_two_index_to_tuple    = T_left->get_two_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<A_left_pairpi[n];i++){
            pq[0]=A_left_tuples[i+A_left_first[n]][0];
            pq[1]=A_left_tuples[i+A_left_first[n]][1];
            for(size_t j = 0;j<A_right_pairpi[n];j++){
              b_irrep = T_left_two_index_to_irrep[pq[reindexing_array[0]]][pq[reindexing_array[1]]];
              b_left  = T_left_two_index_to_tuple[pq[reindexing_array[0]]][pq[reindexing_array[1]]];
              A_matrix[n][i][0]+=constant*T_matrix[b_irrep][b_left][0];
            }
          }
        }
      }
      delete[] pq;
    }else if(reindexing.size()==4){
      int b_irrep;
      size_t b_left,b_right;
      short* pqrs = new short[4];
      // A[x][xxx] <- B[x][xxx]
      if((A_left_nelements == 1) && (T_left_nelements == 1)){
        int*      T_left_one_index_to_irrep    = T_left->get_one_index_to_irrep();
        size_t*   T_left_one_index_to_tuple    = T_left->get_one_index_to_tuple_rel_index();
        size_t*** T_right_three_index_to_tuple = T_right->get_three_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<A_left_pairpi[n];i++){
            //Check that this should be handled by this thread
            pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
            for(size_t j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[1]=A_right_tuples[j+A_right_first[n]][0];
              pqrs[2]=A_right_tuples[j+A_right_first[n]][1];
              pqrs[3]=A_right_tuples[j+A_right_first[n]][2];
              b_irrep = T_left_one_index_to_irrep[pqrs[reindexing_array[0]]];
              b_left  = T_left_one_index_to_tuple[pqrs[reindexing_array[0]]];
              b_right = T_right_three_index_to_tuple[pqrs[reindexing_array[1]]]
                                                    [pqrs[reindexing_array[2]]]
                                                    [pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*T_matrix[b_irrep][b_left][b_right];
            }
          }
        }
      }
      // A[xx][xx] <- B[x][xxx]
      if((A_left_nelements == 2) && (T_left_nelements == 1)){
        int*      T_left_one_index_to_irrep    = T_left->get_one_index_to_irrep();
        size_t*   T_left_one_index_to_tuple    = T_left->get_one_index_to_tuple_rel_index();
        size_t*** T_right_three_index_to_tuple = T_right->get_three_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<A_left_pairpi[n];i++){
            //Check that this should be handled by this thread
            pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
            pqrs[1]=A_left_tuples[i+A_left_first[n]][1];
            for(size_t j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[2]=A_right_tuples[j+A_right_first[n]][0];
              pqrs[3]=A_right_tuples[j+A_right_first[n]][1];
              b_irrep = T_left_one_index_to_irrep[pqrs[reindexing_array[0]]];
              b_left  = T_left_one_index_to_tuple[pqrs[reindexing_array[0]]];
              b_right = T_right_three_index_to_tuple[pqrs[reindexing_array[1]]]
                                                    [pqrs[reindexing_array[2]]]
                                                    [pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*T_matrix[b_irrep][b_left][b_right];
            }
          }
        }
      }

      // A[xxx][x] <- B[x][xxx]
      if((A_left_nelements == 3) && (T_left_nelements == 1)){
        int*      T_left_one_index_to_irrep    = T_left->get_one_index_to_irrep();
        size_t*   T_left_one_index_to_tuple    = T_left->get_one_index_to_tuple_rel_index();
        size_t*** T_right_three_index_to_tuple = T_right->get_three_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<A_left_pairpi[n];i++){
            //Check that this should be handled by this thread
            pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
            pqrs[1]=A_left_tuples[i+A_left_first[n]][1];
            pqrs[2]=A_left_tuples[i+A_left_first[n]][2];
            for(size_t j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[3]=A_right_tuples[j+A_right_first[n]][0];
              b_irrep = T_left_one_index_to_irrep[pqrs[reindexing_array[0]]];
              b_left  = T_left_one_index_to_tuple[pqrs[reindexing_array[0]]];
              b_right = T_right_three_index_to_tuple[pqrs[reindexing_array[1]]]
                                                    [pqrs[reindexing_array[2]]]
                                                    [pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*T_matrix[b_irrep][b_left][b_right];
            }
          }
        }
      }

      // A[x][xxx] <- B[xx][xx]
      if((A_left_nelements == 1) && (T_left_nelements == 2)){
        int**    T_left_two_index_to_irrep    = T_left->get_two_index_to_irrep();
        size_t** T_left_two_index_to_tuple    = T_left->get_two_index_to_tuple_rel_index();
        size_t** T_right_two_index_to_tuple   = T_right->get_two_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<A_left_pairpi[n];i++){
            //Check that this should be handled by this thread
            pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
            for(size_t j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[1]=A_right_tuples[j+A_right_first[n]][0];
              pqrs[2]=A_right_tuples[j+A_right_first[n]][1];
              pqrs[3]=A_right_tuples[j+A_right_first[n]][2];
              b_irrep = T_left_two_index_to_irrep[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              b_left  = T_left_two_index_to_tuple[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              b_right = T_right_two_index_to_tuple[pqrs[reindexing_array[2]]][pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*T_matrix[b_irrep][b_left][b_right];
            }
          }
        }
      }

      // A[xx][xx] <- B[xx][xx]
      if((A_left_nelements == 2) && (T_left_nelements == 2)){
        int** T_left_two_index_to_irrep    = T_left->get_two_index_to_irrep();
        size_t** T_left_two_index_to_tuple    = T_left->get_two_index_to_tuple_rel_index();
        size_t** T_right_two_index_to_tuple   = T_right->get_two_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<A_left_pairpi[n];i++){
            //Check that this should be handled by this thread
            pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
            pqrs[1]=A_left_tuples[i+A_left_first[n]][1];
            for(size_t j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[2]=A_right_tuples[j+A_right_first[n]][0];
              pqrs[3]=A_right_tuples[j+A_right_first[n]][1];
              b_irrep = T_left_two_index_to_irrep[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              b_left  = T_left_two_index_to_tuple[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              b_right = T_right_two_index_to_tuple[pqrs[reindexing_array[2]]][pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*T_matrix[b_irrep][b_left][b_right];
            }
          }
        }
      }

      // A[xxx][x] <- B[xx][xx]
      if((A_left_nelements == 3) && (T_left_nelements == 2)){
        int** T_left_two_index_to_irrep    = T_left->get_two_index_to_irrep();
        size_t** T_left_two_index_to_tuple    = T_left->get_two_index_to_tuple_rel_index();
        size_t** T_right_two_index_to_tuple   = T_right->get_two_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<A_left_pairpi[n];i++){
            //Check that this should be handled by this thread
            pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
            pqrs[1]=A_left_tuples[i+A_left_first[n]][1];
            pqrs[2]=A_left_tuples[i+A_left_first[n]][2];
            for(size_t j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[3]=A_right_tuples[j+A_right_first[n]][0];
              b_irrep = T_left_two_index_to_irrep[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              b_left  = T_left_two_index_to_tuple[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              b_right = T_right_two_index_to_tuple[pqrs[reindexing_array[2]]][pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*T_matrix[b_irrep][b_left][b_right];
            }
          }
        }
      }

      // A[x][xxx] <- B[xxx][x]
      if((A_left_nelements == 1) && (T_left_nelements == 3)){
        size_t*** T_left_three_index_to_tuple  = T_left->get_three_index_to_tuple_rel_index();
        int*   T_right_one_index_to_irrep   = T_right->get_one_index_to_irrep();
        size_t*   T_right_one_index_to_tuple   = T_right->get_one_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<A_left_pairpi[n];i++){
            //Check that this should be handled by this thread
            pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
            for(size_t j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[1]=A_right_tuples[j+A_right_first[n]][0];
              pqrs[2]=A_right_tuples[j+A_right_first[n]][1];
              pqrs[3]=A_right_tuples[j+A_right_first[n]][2];
              b_irrep = T_right_one_index_to_irrep[pqrs[reindexing_array[3]]];
              b_left  = T_left_three_index_to_tuple[pqrs[reindexing_array[0]]]
                                                    [pqrs[reindexing_array[1]]]
                                                    [pqrs[reindexing_array[2]]];
              b_right = T_right_one_index_to_tuple[pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*T_matrix[b_irrep][b_left][b_right];
            }
          }
        }
      }

      // A[xx][xx] <- B[xxx][x]
      if((A_left_nelements == 2) && (T_left_nelements == 3)){
        size_t*** T_left_three_index_to_tuple  = T_left->get_three_index_to_tuple_rel_index();
        int*   T_right_one_index_to_irrep   = T_right->get_one_index_to_irrep();
        size_t*   T_right_one_index_to_tuple   = T_right->get_one_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<A_left_pairpi[n];i++){
            //Check that this should be handled by this thread
            pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
            pqrs[1]=A_left_tuples[i+A_left_first[n]][1];
            for(size_t j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[2]=A_right_tuples[j+A_right_first[n]][0];
              pqrs[3]=A_right_tuples[j+A_right_first[n]][1];
              b_irrep = T_right_one_index_to_irrep[pqrs[reindexing_array[3]]];
              b_left  = T_left_three_index_to_tuple[pqrs[reindexing_array[0]]]
                                                    [pqrs[reindexing_array[1]]]
                                                    [pqrs[reindexing_array[2]]];
              b_right = T_right_one_index_to_tuple[pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*T_matrix[b_irrep][b_left][b_right];
            }
          }
        }
      }

      // A[xxx][x] <- B[xxx][x]
      if((A_left_nelements == 3) && (T_left_nelements == 3)){
        size_t*** T_left_three_index_to_tuple  = T_left->get_three_index_to_tuple_rel_index();
        int*   T_right_one_index_to_irrep   = T_right->get_one_index_to_irrep();
        size_t*   T_right_one_index_to_tuple   = T_right->get_one_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<A_left_pairpi[n];i++){
            //Check that this should be handled by this thread
            pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
            pqrs[1]=A_left_tuples[i+A_left_first[n]][1];
            pqrs[2]=A_left_tuples[i+A_left_first[n]][2];
            for(size_t j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[3]=A_right_tuples[j+A_right_first[n]][0];
              b_irrep = T_right_one_index_to_irrep[pqrs[reindexing_array[3]]];
              b_left  = T_left_three_index_to_tuple[pqrs[reindexing_array[0]]]
                                                    [pqrs[reindexing_array[1]]]
                                                    [pqrs[reindexing_array[2]]];
              b_right = T_right_one_index_to_tuple[pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*T_matrix[b_irrep][b_left][b_right];
            }
          }
        }
      }
      delete[] pqrs;
    } else if(reindexing.size()==6){
      // This is a fast 6-index sorting algorithm used by the active-space computations
      int b_irrep;
      size_t b_left,b_right;
      short* pqrstu = new short[6];
      // A[xxx][xxx] <- T[xxx][xxx]
      if((A_left_nelements == 3) && (T_left_nelements == 3)){
        int***    T_left_three_index_to_irrep    = T_left->get_three_index_to_irrep();
        size_t*** T_left_three_index_to_tuple    = T_left->get_three_index_to_tuple_rel_index();
        size_t*** T_right_three_index_to_tuple   = T_right->get_three_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<A_left_pairpi[n];i++){
            // Get the pqr indices
            pqrstu[0]=A_left_tuples[i+A_left_first[n]][0];
            pqrstu[1]=A_left_tuples[i+A_left_first[n]][1];
            pqrstu[2]=A_left_tuples[i+A_left_first[n]][2];
            for(size_t j = 0;j<A_right_pairpi[n];j++){
              // Get the stu indices
              pqrstu[3]=A_right_tuples[j+A_right_first[n]][0];
              pqrstu[4]=A_right_tuples[j+A_right_first[n]][1];
              pqrstu[5]=A_right_tuples[j+A_right_first[n]][2];

              b_irrep = T_left_three_index_to_irrep[pqrstu[reindexing_array[0]]]
                                                   [pqrstu[reindexing_array[1]]]
                                                   [pqrstu[reindexing_array[2]]];

              b_left  = T_left_three_index_to_tuple[pqrstu[reindexing_array[0]]]
                                                   [pqrstu[reindexing_array[1]]]
                                                   [pqrstu[reindexing_array[2]]];

              b_right = T_right_three_index_to_tuple[pqrstu[reindexing_array[3]]]
                                                    [pqrstu[reindexing_array[4]]]
                                                    [pqrstu[reindexing_array[5]]];
              A_matrix[n][i][j]+=constant*T_matrix[b_irrep][b_left][b_right];
            }
          }
        }
      }
      delete[] pqrstu;
    }
  }


  //
  // A >= B for size(A) > size(B)

  if(assignment==">=" || assignment=="+>="){
    int     a_irrep;
    size_t  a_left,a_right;
    short   swap_temp;
    double  sym_constant = constant * A_Matrix->get_symmetry();
    if(reindexing.size()==4){
      short* pqrs = new short[4];
      // A[x][xxx] <- B[x][xxx]
      if((A_left_nelements == 1) && (T_left_nelements == 1)){
        throw FeatureNotImplemented("PSIMRCC", "A[x][xxx] <- B[x][xxx] with expansion",__FILE__,__LINE__);
      }

      // A[xx][xx] <- B[x][xxx]
      if((A_left_nelements == 2) && (T_left_nelements == 1)){
        throw FeatureNotImplemented("PSIMRCC","A[xx][xx] <- B[x][xxx] with expansion",__FILE__,__LINE__);
      }

      // A[xxx][x] <- B[x][xxx]
      if((A_left_nelements == 3) && (T_left_nelements == 1)){
        throw FeatureNotImplemented("PSIMRCC","A[xxx][x] <- B[x][xxx] with expansion",__FILE__,__LINE__);
      }

      // A[x][xxx] <- B[xx][xx]
      if((A_left_nelements == 1) && (T_left_nelements == 2)){
        throw FeatureNotImplemented("PSIMRCC","A[x][xxx] <- B[xx][xx] with expansion",__FILE__,__LINE__);
      }

      // A[xx][xx] <- B[xx][xx]
      if((A_left_nelements == 2) && (T_left_nelements == 2)){
        int** A_left_two_index_to_irrep    = A_left->get_two_index_to_irrep();
        size_t** A_left_two_index_to_tuple    = A_left->get_two_index_to_tuple_rel_index();
        size_t** A_right_two_index_to_tuple   = A_right->get_two_index_to_tuple_rel_index();
        for(int n=0;n<moinfo->get_nirreps();n++){
          for(size_t i = 0;i<T_left_pairpi[n];i++){
            pqrs[0]=T_left_tuples[i+T_left_first[n]][0];
            pqrs[1]=T_left_tuples[i+T_left_first[n]][1];
            for(size_t j = 0;j<T_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[2]=T_right_tuples[j+T_right_first[n]][0];
              pqrs[3]=T_right_tuples[j+T_right_first[n]][1];
              a_irrep = A_left_two_index_to_irrep[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              a_left  = A_left_two_index_to_tuple[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              a_right = A_right_two_index_to_tuple[pqrs[reindexing_array[2]]][pqrs[reindexing_array[3]]];
              A_matrix[a_irrep][a_left][a_right]+=constant*T_matrix[n][i][j];
              // Add the expandend term with the correct symmetry
              swap_temp =pqrs[2];
              pqrs[2] = pqrs[3];
              pqrs[3] = swap_temp;
              a_irrep = A_left_two_index_to_irrep[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              a_left  = A_left_two_index_to_tuple[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              a_right = A_right_two_index_to_tuple[pqrs[reindexing_array[2]]][pqrs[reindexing_array[3]]];
              A_matrix[a_irrep][a_left][a_right]+=sym_constant*T_matrix[n][i][j];
            }
          }
        }
      }

      // A[xxx][x] <- B[xx][xx]
      if((A_left_nelements == 3) && (T_left_nelements == 2)){
        throw FeatureNotImplemented("PSIMRCC","A[xxx][x] <- B[xx][xx] with expansion",__FILE__,__LINE__);
      }

      // A[x][xxx] <- B[xxx][x]
      if((A_left_nelements == 1) && (T_left_nelements == 3)){
        throw FeatureNotImplemented("PSIMRCC","A[x][xxx] <- B[xxx][x] with expansion",__FILE__,__LINE__);
      }

      // A[xx][xx] <- B[xxx][x]
      if((A_left_nelements == 2) && (T_left_nelements == 3)){
        throw FeatureNotImplemented("PSIMRCC","A[xx][xx] <- B[xxx][x] with expansion",__FILE__,__LINE__);
      }

      // A[xxx][x] <- B[xxx][x]
      if((A_left_nelements == 3) && (T_left_nelements == 3)){
        throw FeatureNotImplemented("PSIMRCC","A[xxx][x] <- B[xxx][x] with expansion",__FILE__,__LINE__);
      }

      delete[] pqrs;
    }
  }

  delete[] reindexing_array;
  sort_timing += sort_timer.get();
}


}} /* End Namespaces */
