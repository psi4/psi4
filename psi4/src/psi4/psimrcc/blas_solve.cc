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

#include <cstdio>
#include "psi4/libmoinfo/libmoinfo.h"

#include "blas.h"
#include "debugging.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

using namespace std;

/**
 * Read and compute an expression
 * @param cstr
 */
void CCBLAS::solve(const char* cstr)
{
  string str(cstr);
  solve(str);
}

/**
 * Read and compute an expression
 * @param str
 */
void CCBLAS::solve(string str)
{
  append(str);
  compute();
}

/**
 * Read and store expressions without computing them
 * @param cstr
 */
void CCBLAS::append(const char* cstr)
{
  string str(cstr);
  append(str);
}

/**
 * Read and store expressions without computing them
 * @param str
 */
void CCBLAS::append(string str)
{
  // Main driver for solving expressions
  int noperations_added = 0;
  DEBUGGING(5,
    outfile->Printf("\n\nYou have requested the following operation :\n\t\"%s\"",str.c_str());
    outfile->Printf("\n\nCCBLAS::append() has parsed the following:");
  )
  vector<string> names = moinfo->get_matrix_names(str);
  for(int n=0;n<names.size();n++){
    noperations_added+=parse(names[n]);
  }
}

/**
 * Flush the operation deque
 */
// void CCBLAS::compute()
// {
//   while(!operations.empty()){
//     // Read the element
//     CCOperation& op = operations.front();
//     // Compute the operation
//     op.compute();
//     // Eliminate the element
//     operations.pop_front();
//   }
// }

/**
 * Flush the operation deque in a memory smart way!
 */
void CCBLAS::compute()
{

//   outfile->Printf("\n\nsmart_compute::size of current deque = %d",operations.size());
//
//   outfile->Printf("\nsmart_compute::content of the deque:");
//   for(OpDeque::iterator it = operations.begin();it!=operations.end();++it){
//     outfile->Printf("\n");
//     it->print();
//   }

//   // Create a map with all the target matrices and how many times they appear
//   map<string,int> target_count;
//   for(OpDeque::iterator it = operations.begin();it!=operations.end();++it){
//     if(it->get_A_Matrix()!=NULL)
//       target_count[it->get_A_Matrix()->get_label()]++;
//   }
//   // Create a map with all the source matrices and how many times they appear
//   map<string,int> source_count;
//   for(OpDeque::iterator it = operations.begin();it!=operations.end();++it){
//     if(it->get_B_Matrix()!=NULL)
//       source_count[it->get_B_Matrix()->get_label()]++;
//     if(it->get_C_Matrix()!=NULL)
//       source_count[it->get_C_Matrix()->get_label()]++;
//   }
//
//   // Create a map with all the intermediates defined as matrices that appear a source and target
//     map<string,int> intermediates_count;
//   for(map<string,int>::iterator it=matrix_count.begin();it!=matrix_count.end();++it){
//     map<string,int>::iterator target_it = target_count.find(it->first);
//     map<string,int>::iterator source_it = source_count.find(it->first);
//     if( (target_it!=target_count.end()) && (source_it!=source_count.end()))
//       intermediates_count[source_it->first]=source_it->second;
//   }
//
//   // Print the map for debugging purposes
//   outfile->Printf("\n\nsmart_compute::printing the matrix_count map");
//   for(map<string,int>::iterator it=matrix_count.begin();it!=matrix_count.end();++it){
//     outfile->Printf("\n %s(%d)",it->first.c_str(),it->second);
//   }
//
//   // Print the map for debugging purposes
//   outfile->Printf("\n\nsmart_compute::printing the target_count map");
//   for(map<string,int>::iterator it=target_count.begin();it!=target_count.end();++it){
//     outfile->Printf("\n %s(%d)",it->first.c_str(),it->second);
//   }
//
//   // Print the map for debugging purposes
//   outfile->Printf("\n\nsmart_compute::printing the source_count map");
//   for(map<string,int>::iterator it=source_count.begin();it!=source_count.end();++it){
//     outfile->Printf("\n %s(%d)",it->first.c_str(),it->second);
//   }
//
//   // Print the map for debugging purposes
//   outfile->Printf("\n\nsmart_compute::printing the intermediates_count map");
//   for(map<string,int>::iterator it=intermediates_count.begin();it!=intermediates_count.end();++it){
//     outfile->Printf("\n %s(%d)",it->first.c_str(),it->second);
//   }

//   map<string,bool> matrix_on_disk;
//   map<string,bool> matrix_in_core;
//   for(map<string,int>::iterator it=matrix_count.begin();it!=matrix_count.end();++it){
//     if(it->first[0]=='<'
//     matrix_on_disk[it->first]=false;
//     matrix_on_disk[it->first]=false;
//   }

//   for(OpDeque::iterator it = operations.begin();it!=operations.end();++it){
//
//   }
//   static int memory_for_solve = 0;
//
//   for(OpDeque::iterator it = operations.begin();it!=operations.end();++it){
//     int nop = distance(it-operations.begin());
//     outfile->Printf("\nI will process operation #%3d",nop);
//     int memA=0;
//     if(it->get_A_Matrix()!=NULL){
//       memA=it->get_A_Matrix()->get_memory();
//     }
//   }

  // Create a map with all the required matrices and how many times they appear

  for(OpDeque::iterator it = operations.begin();it!=operations.end();++it){
    if(it->get_A_Matrix()!=NULL){
      matrices_in_deque[it->get_A_Matrix()]++;
      matrices_in_deque_target[it->get_A_Matrix()]++;
    }
    if(it->get_B_Matrix()!=NULL){
      matrices_in_deque[it->get_B_Matrix()]++;
      matrices_in_deque_source[it->get_B_Matrix()]++;
    }
    if(it->get_C_Matrix()!=NULL){
      matrices_in_deque[it->get_C_Matrix()]++;
      matrices_in_deque_source[it->get_C_Matrix()]++;
    }
  }
  while(!operations.empty()){
    // Read the element
    CCOperation& op = operations.front();
    // Compute the operation
    op.compute();

    // Decrease the counters for the matrices to be processed
    if(op.get_A_Matrix()!=NULL){
      matrices_in_deque[op.get_A_Matrix()]--;
      matrices_in_deque_target[op.get_A_Matrix()]--;
    }
    if(op.get_B_Matrix()!=NULL){
      matrices_in_deque[op.get_B_Matrix()]--;
      matrices_in_deque_source[op.get_B_Matrix()]--;
    }
    if(op.get_C_Matrix()!=NULL){
      matrices_in_deque[op.get_C_Matrix()]--;
      matrices_in_deque_source[op.get_C_Matrix()]--;
    }
    // Eliminate the element
    operations.pop_front();
  }
}

/**
 * store a zero_two_diagonal operation without executing it
 * @param cstr
 */
void CCBLAS::append_zero_two_diagonal(const char* cstr)
{
  string str(cstr);
  // To zero diagonals of things like "Fae[v][v]{u}"
  vector<string> names = moinfo->get_matrix_names(str);
  for(int n=0;n<names.size();n++){
    CCMatrix* Matrix = get_Matrix(names[n]);
    CCOperation op(0.0,"","","zero_two_diagonal",Matrix,NULL,NULL,work[0],buffer[0]);
    operations.push_back(op);
  }
}

/**
 * store a zero_two_diagonal operation and executing it
 * @param str
 */
void CCBLAS::solve_zero_two_diagonal(const char* cstr)
{
  append_zero_two_diagonal(cstr);
  compute();
}


}} /* End Namespaces */
