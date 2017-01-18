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

#ifndef _psi_src_bin_psimrcc_ccoperation_h
#define _psi_src_bin_psimrcc_ccoperation_h

#include <string>

namespace psi{ namespace psimrcc{

class CCIndex;
class CCMatrix;

class CCOperation{
  // Used to define operations of the type:
  // A [+ -]= ## factor B [*  / @] C
  public:
    CCOperation(double in_factor,std::string in_assignment,
          std::string in_reindexing,std::string in_operation,
          CCMatrix* in_A_Matrix, CCMatrix* in_B_Matrix, CCMatrix* in_C_Matrix,double* work,double* buffer);
    ~CCOperation();
    double      get_factor()    {return(factor);}
    std::string get_assignment(){return(assignment);}
    std::string get_reindexing(){return(reindexing);}
    std::string get_operation() {return(operation);}
    CCMatrix*   get_A_Matrix()  {return(A_Matrix);}
    CCMatrix*   get_B_Matrix()  {return(B_Matrix);}
    CCMatrix*   get_C_Matrix()  {return(C_Matrix);}
    void        print();
    void        print_operation();
    void        compute();
    static void print_timing();
  private:
  //            Variable        Syntax (p,q,r,s=integers)
    double      factor;     // + - +p/q -p/q
    std::string assignment; // = += >= +>=
    std::string reindexing; // ## #pq# #pqrs#
    std::string operation;  // . @ / * X plus
    static double* out_of_core_buffer;
    static double* local_work;
    CCMatrix*   A_Matrix;
    CCMatrix*   B_Matrix;
    CCMatrix*   C_Matrix;
  private:
  // Check that an operation can be performed
    bool compatible_dot();
    bool compatible_contract();
    bool compatible_element_by_element();
  //
    static double zero_timing;
    static double numerical_timing;
    static double contract_timing;
    static double tensor_timing;
    static double dot_timing;
    static double plus_timing;
    static double product_timing;
    static double division_timing;
    static double sort_timing;
    static double PartA_timing;
    static double PartB_timing;
    static double PartC_timing;
    void fail_to_compute();
    void add_numerical_factor();
    void contract();
    void dot_product();
    void element_by_element_product();
    void element_by_element_division();
    void element_by_element_addition();
    void different_index_contract(int A_type,double factor,int B_type,std::string& operation,int C_type,std::string& reindexing);
    void reindex(double constant);
    void *reindex_thread(void * my_id);
    void *reindex_thread_right_expand(void * my_id);
    void tensor_product();
    void setup_contractions();
    void contract_in_core(double** A_matrix,double** B_matrix,double** C_matrix,bool B_on_disk,bool C_on_disk,int rows_A,int rows_B,int rows_C,int cols_A,int cols_B,int cols_C,int offset);
    void sort(CCIndex* T_left,CCIndex* T_right,double*** T_matrix,double constant);
    void sort();
    void check_and_zero_target();
    void check_and_zero_target_block(int h);
    void zero_target();
    void zero_target_block(int h);
    void zero_two_diagonal();
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_ccoperation_h