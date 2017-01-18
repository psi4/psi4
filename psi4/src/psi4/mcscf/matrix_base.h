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

#ifndef _psi_src_lib_libmemtrix_matrix_base_h_
#define _psi_src_lib_libmemtrix_matrix_base_h_

#include <cstring> // for size_t

#include "psi4/libpsi4util/memory_manager.h"
#include "vector_base.h"

namespace psi{ namespace mcscf{

class MatrixBase
{
public:
  MatrixBase();
  MatrixBase(size_t rows, size_t cols);
  ~MatrixBase();

  //Inlines
  size_t  get_rows()                  {return(rows_);}
  size_t  get_cols()                  {return(cols_);}
  size_t  get_elements()              {return(elements_);}
  void    set(size_t i, size_t j, double value) {matrix_[i][j]  = value;}
  void    add(size_t i, size_t j, double value) {matrix_[i][j] += value;}
  double  get(size_t i, size_t j)               {return(matrix_[i][j]);}
  double** get_matrix()                   {return(matrix_);}

  void    multiply(bool transpose_A, bool transpose_B, MatrixBase* A, MatrixBase* B);
  void    diagonalize(MatrixBase* eigenmatrix, VectorBase* eigenvalues);
  void    print();
  void    zero();
  void    zero_diagonal();
  void    scale(double factor);
  void    transpose();
  MatrixBase& operator+=(const MatrixBase& rhs);
  MatrixBase& operator-=(const MatrixBase& rhs);
  friend double dot(MatrixBase* A, MatrixBase* B);
private:
  // Matrix size
  size_t  rows_;
  size_t  cols_;
  size_t  elements_;

  // Matrix data
  double** matrix_;
};

}}

#endif // _psi_src_lib_libmemtrix_matrix_base_h_

/*

  int       get_nirreps()                          {return(nirreps);}
  int       get_nirreps() const                    {return(nirreps);}
  double*** get_matrix()                           {return(matrix);}
  double**  get_block(int h)                       {return(matrix[h]);}
  const double**  get_block(int h) const           {return((const double**)matrix[h]);}
  double    get(int h, int i, int j)               {return(matrix[h][i][j]);}
  void      set(int h, int i, int j, double value) {matrix[h][i][j]  = value;}
  void      add(int h, int i, int j, double value) {matrix[h][i][j] += value;}
  int       get_block_size(int h)                  {return(block_size[h]);}
  int       get_block_size(int h) const            {return(block_size[h]);}
  void      diagonalize(SBlockMatrix* eigenvectors, double* eigenvalues);
  void      DGEMM(bool transpose_A, bool transpose_B, SBlockMatrix* A, SBlockMatrix* B);
  void      minus(SBlockMatrix* B);

  friend double  operator^(const SBlockMatrix& rhs,const SBlockMatrix& lhs);
  SBlockMatrix& operator=(const SBlockMatrix& rhs);
  SBlockMatrix& operator+=(const SBlockMatrix& lhs);
*/
