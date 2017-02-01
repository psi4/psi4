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

#include <cstdlib>
#include <cstdio>
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "blas.h"
#include "index.h"
#include "matrix.h"


namespace psi{

    namespace psimrcc{

using namespace std;

bool CCOperation::compatible_dot()
{
  //
  //  Here we are doing the following check:
  //            --------------
  //           |              |
  //   A = B[B_l][B_r] . C[C_l][C_r]
  //               |               |
  //                ---------------
  //
  bool same = false;

  int A_left   = A_Matrix->get_left()->get_ntuples();
  int A_right  = A_Matrix->get_right()->get_ntuples();
  int B_left   = B_Matrix->get_left()->get_ntuples();
  int B_right  = B_Matrix->get_right()->get_ntuples();
  int C_left   = C_Matrix->get_left()->get_ntuples();
  int C_right  = C_Matrix->get_right()->get_ntuples();

  if((A_left==1) && (B_left==C_left) && (A_right==1) && (B_right==C_right))
    same = true;
  if(!same){
    outfile->Printf("\n\nSolve couldn't perform the operation ");
    print_operation();

    exit(1);
  }
  return(same);
}

bool CCOperation::compatible_element_by_element()
{
  //
  //  Here we are doing the following check:
  //            ---------------------------
  //           |             |             |
  //   A[A_l][A_r] = B[B_l][B_r] / C[C_l][C_r]
  //      |             |             |
  //       ---------------------------
  //
  bool same = false;

  int A_left   = A_Matrix->get_left()->get_ntuples();
  int A_right  = A_Matrix->get_right()->get_ntuples();
  int B_left   = B_Matrix->get_left()->get_ntuples();
  int B_right  = B_Matrix->get_right()->get_ntuples();

  // First case: We are comparing A and B
  if(C_Matrix==NULL){
    if((A_left==B_left) && (A_right==B_right))
      same = true;
  }else{
  // Second case: We are comparing A, B, and C
    int C_left   = C_Matrix->get_left()->get_ntuples();
    int C_right  = C_Matrix->get_right()->get_ntuples();
    if((A_left==B_left) && (B_left==C_left) && (A_right==B_right) && (B_right==C_right))
      same = true;
    if((B_left!=C_left) || (B_right!=C_right)){
      outfile->Printf("\n\nSolve couldn't perform the operation ");
      print_operation();

      exit(1);
    }
  }
  return(same);
}

bool CCOperation::compatible_contract()
{
  //
  //  Here we are doing the following check:
  //            -----------------------------
  //           |                             |
  //   A[A_l][A_r] = B[B_i][B_c] 2@1 C[C_c][C_i]
  //      |             |    |          |
  //       -------------      ----------
  //
  bool same = false;

  int A_left   = A_Matrix->get_left()->get_ntuples();
  int A_right  = A_Matrix->get_right()->get_ntuples();
  int B_contracted;
  int B_index;
  int C_contracted;
  int C_index;

  if(operation[0]=='1'){
    B_contracted   = B_Matrix->get_left()->get_ntuples();
    B_index        = B_Matrix->get_right()->get_ntuples();
  } else{
    B_contracted   = B_Matrix->get_right()->get_ntuples();
    B_index        = B_Matrix->get_left()->get_ntuples();
  }
  if(operation[2]=='1'){
    C_contracted   = C_Matrix->get_left()->get_ntuples();
    C_index        = C_Matrix->get_right()->get_ntuples();
  } else{
    C_contracted   = C_Matrix->get_right()->get_ntuples();
    C_index        = C_Matrix->get_left()->get_ntuples();
  }

  if((B_contracted==C_contracted) && (A_left==B_index) && (A_right==C_index))
    same = true;
  if(B_contracted!=C_contracted){
    outfile->Printf("\n\nSolve couldn't perform the operation ");
    print_operation();

    exit(1);
  }
  return(same);
}

}} /* End Namespaces */
