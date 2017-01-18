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

#include "psi4/libmoinfo/libmoinfo.h"
#include <cstdio>
#include "blas.h"
#include "debugging.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include <algorithm>

namespace psi{

    namespace psimrcc{

typedef std::vector<std::string>            strvec;
typedef std::vector<std::pair<int,int> >    intpairvec;

bool is_number(const std::string& str);
double get_number(const std::string& str);
bool is_operation(const std::string& str);

// Parsing Algorithm
int CCBLAS::parse(std::string& str){

  int noperations_added = 0;
  // Split the operation
  strvec split_str = split(str);

  // Define the allowed operations
  strvec allowed_assignments = split("= += >= +>=");

  // Store the A Matrix
  CCMatrix* A_Matrix = get_Matrix(split_str[0],str);
  // Check the assigment operator
  std::string assignment(split_str[1]);
  strvec::iterator strveciter(find(allowed_assignments.begin(),allowed_assignments.end(),assignment));
  if(strveciter==allowed_assignments.end())
    outfile->Printf("\n\nCCBLAS::parse() %s is not a proper assigment\n\nin the expression:\n\t%s\n\n",assignment.c_str(),str.c_str());

  // Eliminate the first two strings and store the rest of the terms
  strvec::iterator iter = split_str.begin();
  iter=split_str.erase(iter,iter+2);
  std::string   reindexing;
  std::string   operation;
  CCMatrix*     B_Matrix;
  CCMatrix*     C_Matrix;

  while(iter!=split_str.end()){
    double factor=1.0;
    C_Matrix = NULL;
    B_Matrix = NULL;
    // Read the reindexing
    if(iter->find("#")!=std::string::npos){
      reindexing = *iter;
      reindexing = reindexing.substr(1,reindexing.size()-2);
      ++iter;
    }
    // Read a product of numerical factors
    while((iter!=split_str.end()) && get_factor(*iter,factor))
      ++iter;
    operation="add_factor";
    // Read the first matrix
    if(iter!=split_str.end()){
      B_Matrix = get_Matrix(*iter,str);
      ++iter;
      // And eventually an operation and a matrix
      if(iter!=split_str.end()){
        if(is_operation(*iter)){
          // Check operation
          operation = *iter;
          ++iter;
          if(iter!=split_str.end()){
            C_Matrix = get_Matrix(*iter,str);
            ++iter;
          }
        }else{
          operation="plus";
        }
      }else{ // When you reach the end
        operation="plus";
      }
    }
    if(noperations_added && assignment[0]!='+')
      assignment= "+" + assignment;
    CCOperation op(factor,assignment,reindexing,operation,A_Matrix,B_Matrix,C_Matrix,work[0],buffer[0]);
    operations.push_back(op);
    noperations_added++;
    DEBUGGING(5,
      op.print();
    )
  }
  return(noperations_added);
}

bool CCBLAS::get_factor(const std::string& str,double& factor)
{
  if(is_number(str)){
    factor *= get_number(str);
    return true;
  }
  if(str=="-"){
    factor *= -1.0;
    return true;
  }
  if(str=="+"){
    factor *= 1.0;
    return true;
  }
  if(str.substr(0,6)=="factor"){
    factor = get_scalar(str);
    return true;
  }
  return false;
}

bool is_number(const std::string& str){
  // Is the string str a number?
  static const std::string numbers = "1234567890.-+/e";
  bool numeric = true;
  for(int i=0;i<str.size();i++)
      if(numbers.find(str[i])==std::string::npos) numeric = false;
  // In case we have only a symbol (ex. "+")
  if((str.size()==1) && !isdigit(str[0])) numeric = false;
  return(numeric);
}

double get_number(const std::string& str){
  double value = 1.0;
  bool fraction = false;
  size_t fraction_sign;
  for(int i=0;i<str.size();i++)
    if(str[i]=='/'){
      fraction = true;
      fraction_sign = i;
    }
  if(fraction){
    std::string numerator   = str.substr(0,fraction_sign);
    std::string denominator = str.substr(fraction_sign+1,str.size()-fraction_sign-1);
    std::string unsigned_numerator = find_and_replace(numerator,"-","");
    if(unsigned_numerator.size() * denominator.size()!=0){
      value=to_double(numerator)/to_double(denominator);
    }else{
      outfile->Printf("\n\nSolve couldn't parse the numerical factor %s\n\n",str.c_str());
      outfile->Printf("\n\nCritical Breakdown of the Program. Blame the programmers!!!\n\n");

      exit(1);
    }

  }else{
    value=to_double(str);
  }
  return(value);
}

bool is_operation(const std::string& str){
  // Is the string str a number?
  strvec allowed_operations = split(". @ / * X");
  bool operation = false;
  for(int i=0;i<allowed_operations.size();i++)
      if(str.find(allowed_operations[i])!=std::string::npos) operation = true;
  return(operation);
}


// void split_index(const std::string& str,std::string& t2_label, int*& t2_indices,std::string& first_t1_label, int*& first_t1_indices,std::string& second_t1_label, int*& second_t1_indices)
// {
//   size_t opening = str.find_first_of("[");
//   string index_str = str.substr(opening,8);
//
//   size_t opening_ref = str.find_first_of("{");
//   size_t closing_ref = str.find_first_of("}");
//   string ref_str     = str.substr(opening_ref,closing_ref-opening_ref+1);
//
//   intpairvec pairs;
//   std::vector<char> labels;
//   int index=0;
//   for(int i=0;i<index_str.size();i++){
//     if(index_str[i]=='o' || index_str[i]=='O' || index_str[i]=='v' || index_str[i]=='V') labels.push_back(index_str[i]);
//
//     if(index_str[i]=='o') pairs.push_back(make_pair(0,index++));
//     if(index_str[i]=='O') pairs.push_back(make_pair(1,index++));
//     if(index_str[i]=='v') pairs.push_back(make_pair(2,index++));
//     if(index_str[i]=='V') pairs.push_back(make_pair(3,index++));
//   }
//   sort(pairs.begin(),pairs.end());
//   t2_label = "[";
//   for(int i=0;i<index;i++){
//     t2_indices[i]=pairs[i].second;
//     t2_label += labels[t2_indices[i]];
//     if(i==1)
//       t2_label += "][";
//   }
//   t2_label += "]" + ref_str;
//
//   first_t1_indices[0]=pairs[0].second;
//   first_t1_indices[1]=pairs[2].second;
//   second_t1_indices[0]=pairs[1].second;
//   second_t1_indices[1]=pairs[3].second;
//   first_t1_label = "[";
//   second_t1_label = "[";
//   for(int i=0;i<2;i++){
//     first_t1_label += labels[first_t1_indices[i]];
//     second_t1_label += labels[second_t1_indices[i]];
//     if(i==0){
//       first_t1_label += "][";
//       second_t1_label += "][";
//     }
//   }
//   first_t1_label += "]" + ref_str;
//   second_t1_label += "]" + ref_str;
//
// //   if(moinfo->get_debug()>7){
// //     printf("\n\nInput string : %s",str.c_str());
// //     printf("\n\nIndex string : %s + %s + %s",t2_label.c_str(),first_t1_label.c_str(),second_t1_label.c_str());
// //   }
// }




}} /* End Namespaces */
