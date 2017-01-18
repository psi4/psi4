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

/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include <cstdlib>
#include <cstdio>

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "algebra_interface.h"
#include "blas.h"
#include "debugging.h"
#include "matrix.h"
#include "operation.h"

namespace psi{

    namespace psimrcc{
    extern MOInfo *moinfo;

using namespace std;

/**
 * This is the core of the CCBLAS class. Computes the expression in the CCOperation class.
 */
void CCOperation::compute()
{
  // Here we distinguis between all the possible cases
  DEBUGGING(2,
    outfile->Printf("\nPerforming ");
    print_operation();
  );

  Timer numerical_timer;
  // (1) Assigment of a number
  //     Expression of the type A = - 1/2
  if(operation=="add_factor")
    add_numerical_factor();
  numerical_timing += numerical_timer.get();

  Timer dot_timer;
  // (2) Dot Product
  //     operation = .
  if(operation==".")
    dot_product();
  dot_timing += dot_timer.get();

  Timer contract_timer;
  // (2) Contraction
  //     operation = i@j
  if(operation.substr(1,1)=="@")
    contract();
  contract_timing += contract_timer.get();

  Timer plus_timer;
  // (4) Add a matrix
  //     operation = plus
  if(operation=="plus")
     element_by_element_addition();
  plus_timing += plus_timer.get();

  Timer tensor_timer;
  // (5) Tensor Product of two matrices
  //     operation = X
  if(operation=="X")
    tensor_product();
  tensor_timing += tensor_timer.get();

  Timer product_timer;
  // (6) Element by element product
  //     operation = *
  if(operation=="*")
    element_by_element_product();
  product_timing += product_timer.get();

  Timer division_timer;
  // (7) Element by element division
  //     operation = /
  if(operation=="/")
    element_by_element_division();
  division_timing += division_timer.get();

  // (8) Zero two diagonal
  //     operation = "zero_two_diagonal"
  if(operation=="zero_two_diagonal")
    zero_two_diagonal();
}


/**
 * Add a number to each element of a matrix
 */
void CCOperation::add_numerical_factor()
{
  for(int h=0;h<moinfo->get_nirreps();++h){
    CCMatIrTmp AMatIrTmp = blas->get_MatIrTmp(A_Matrix,h,none);
    check_and_zero_target_block(h);
    AMatIrTmp->add_numerical_factor(factor,h);
  }
}

void CCOperation::dot_product()
{
  if(compatible_dot()){
    double dot_product = 0.0;
    for(int h=0;h<moinfo->get_nirreps();++h){
      CCMatIrTmp BMatIrTmp = blas->get_MatIrTmp(B_Matrix,h,none);
      CCMatIrTmp CMatIrTmp = blas->get_MatIrTmp(C_Matrix,h,none);
      dot_product += CCMatrix::dot_product(BMatIrTmp.get_CCMatrix(),CMatIrTmp.get_CCMatrix(),h);
    }
    CCMatTmp AMatTmp = blas->get_MatTmp(A_Matrix,none);
    if(assignment=="=" || assignment==">=")
      AMatTmp->set_scalar(dot_product*factor);
    else
      AMatTmp->add_scalar(dot_product*factor);
  }else
    fail_to_compute();
}

void CCOperation::element_by_element_product()
{
  if(compatible_element_by_element()){
    for(int h=0;h<moinfo->get_nirreps();++h){
      CCMatIrTmp AMatIrTmp = blas->get_MatIrTmp(A_Matrix,h,none);
      check_and_zero_target_block(h);
      CCMatIrTmp BMatIrTmp = blas->get_MatIrTmp(B_Matrix,h,none);
      CCMatIrTmp CMatIrTmp = blas->get_MatIrTmp(C_Matrix,h,none);
      AMatIrTmp->element_by_element_product(factor,BMatIrTmp.get_CCMatrix(),CMatIrTmp.get_CCMatrix(),h);
    }
  }else
    fail_to_compute();
}

void CCOperation::element_by_element_division()
{
  if(compatible_element_by_element()){
    for(int h=0;h<moinfo->get_nirreps();++h){
      CCMatIrTmp AMatIrTmp = blas->get_MatIrTmp(A_Matrix,h,none);
      check_and_zero_target_block(h);
      CCMatIrTmp BMatIrTmp = blas->get_MatIrTmp(B_Matrix,h,none);
      CCMatIrTmp CMatIrTmp = blas->get_MatIrTmp(C_Matrix,h,none);
      AMatIrTmp->element_by_element_division(factor,BMatIrTmp.get_CCMatrix(),CMatIrTmp.get_CCMatrix(),h);
    }
  }else
    fail_to_compute();
}

void CCOperation::element_by_element_addition()
{
  if(compatible_element_by_element() && (reindexing.size()==0)){
    DEBUGGING(4,
      outfile->Printf("\n...same indexing for the target and the output of this operation");
    );
    for(int h=0;h<moinfo->get_nirreps();++h){
      CCMatIrTmp AMatIrTmp = blas->get_MatIrTmp(A_Matrix,h,none);
      check_and_zero_target_block(h);
      CCMatIrTmp BMatIrTmp = blas->get_MatIrTmp(B_Matrix,h,none);
      AMatIrTmp->element_by_element_addition(factor,BMatIrTmp.get_CCMatrix(),h);
    }
  }else if(reindexing.size()!=0){
    DEBUGGING(4,
      outfile->Printf("\n...different indexing for the target and the output of this operation");
    );
    CCMatTmp AMatTmp = blas->get_MatTmp(A_Matrix,none);
    check_and_zero_target();
    CCMatTmp BMatTmp = blas->get_MatTmp(B_Matrix,none);
    sort();
  }else
    fail_to_compute();
}

void CCOperation::tensor_product()
{
  DEBUGGING(4,
    outfile->Printf("\n...different indexing for the target and the output of this operation"););
  if(reindexing.size()==0)
    reindexing = "1234";
  // Perform this for all the matrix at once
  CCMatTmp AMatTmp = blas->get_MatTmp(A_Matrix,none);
  check_and_zero_target();
  CCMatTmp BMatTmp = blas->get_MatTmp(B_Matrix,none);
  CCMatTmp CMatTmp = blas->get_MatTmp(C_Matrix,none);
  AMatTmp->tensor_product(reindexing,factor,BMatTmp.get_CCMatrix(),CMatTmp.get_CCMatrix());
}


void CCOperation::contract()
{
  if(compatible_contract() && (reindexing.size()==0)){
    // Same indexing contract, let BLAS directly handle this (although we guide it)
    DEBUGGING(4,
      outfile->Printf("\n...same indexing for the target and the output of this operation"););
  }else{
    // Different indexing contract, we work this by hand at first
    DEBUGGING(4,
      outfile->Printf("\n...different indexing for the target and the output of this operation"););
  }
  setup_contractions();
}

void CCOperation::zero_two_diagonal()
{
  A_Matrix->zero_two_diagonal();
}

void CCOperation::check_and_zero_target()
{
  if(assignment=="=" || assignment==">=")
    zero_target();
}

void CCOperation::check_and_zero_target_block(int h)
{
  if(assignment=="=" || assignment==">=")
    zero_target_block(h);
}

void CCOperation::zero_target()
{
  for(int h=0;h<moinfo->get_nirreps();++h)
    A_Matrix->zero_matrix();
}

void CCOperation::zero_target_block(int h)
{
  Timer zero_timer;
  A_Matrix->zero_matrix_block(h);
  zero_timing += zero_timer.get();
}

void CCOperation::fail_to_compute()
{
  outfile->Printf("\n\nSolve couldn't perform the operation ");
  print_operation();

  exit(EXIT_FAILURE);
}

//   Timer zero_timer;
//   // Parse the assignment for "= >=" and in this case zero A_Matrix
//   if(assignment=="=" || assignment==">="){
//     DEBUGGING(4,
//       outfile->Printf("\n...zero the target Matrix");
//
//   }
//   zero_timing += zero_timer.get();

// void CCBlas::compute(std::vector<double>& factors,std::vector<int>& types,std::vector<std::string>& operations,std::string& reindexing)
// {
//   double factor;
//   int A_type,B_type,C_type;
//   int group_index = 1;
//   for(int group=0;group<factors.size();group++){
//     // perform types[group*2+1] operations[group+1] types[group*2+2]
//     factor    = factors[group];
//     A_type    = types[0];
//     B_type    = types[group_index++];
//     string operation = operations[group+1];
//     if(operation!="plus")
//       C_type    = types[group_index++];
//
//     if(operation!="plus"){
//       DEBUGGING(4,
//
//         factor,matrix_labels(B_type).c_str(),operation.c_str(),matrix_labels(C_type).c_str(),matrix_labels(A_type).c_str());
//     }else{
//       DEBUGGING(4,
//       outfile->Printf("\n\nPerforming %lf %s -> %s",
//         factor,matrix_labels(B_type).c_str(),matrix_labels(A_type).c_str());
//     }
//
//     // What is the size of the first operation? This determines whether
//     // we are incrementing the existing array or assigning a new value
//     if((group==0) && (operations[0]=="=")){
//       DEBUGGING(4,
//         outfile->Printf("\n...zero the target Matrix");
//       zero_matrix(types[0]);
//     }
//
//
//   }
// }

}} /* End Namespaces */
