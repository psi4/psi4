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

#ifndef _psi_src_bin_psimrcc_heff_h_
#define _psi_src_bin_psimrcc_heff_h_

#include <vector>

namespace psi{ namespace psimrcc{

class Hamiltonian
{
public:
  // Constructor and destructor
  Hamiltonian();
  ~Hamiltonian();

  double get_eigenvalue()               const {return eigenvalue;}
  double get_matrix(int mu,int nu)      const {return matrix[mu][nu];}
  double get_left_eigenvector(int mu)   const {return left_eigenvector[mu];}
  double get_right_eigenvector(int mu)  const {return right_eigenvector[mu];}
  double get_zeroth_order_eigenvector(int mu)  const {return zeroth_order_eigenvector[mu];}

  double expectation_value();
  double diagonalize(int root = 0);

  double trace();

  void add_matrix(int mu,int nu,double value) {matrix[mu][nu] += value;}

  void set_eigenvalue(double eigenvalue_) {eigenvalue = eigenvalue_;}
  void set_matrix(double** M,int ndets);
  void set_zeroth_order_eigenvector(double* v,int ndets);
  void set_left_eigenvector(double* v,int ndets);
  void set_right_eigenvector(double* v,int ndets);
  void print();
  void print_matrix();

private:
  void startup();
  void cleanup();

  int ndets;
  double eigenvalue;
  std::vector<std::vector<double> > matrix;
  std::vector<double> right_eigenvector;
  std::vector<double> left_eigenvector;
  std::vector<double> zeroth_order_eigenvector;
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_heff_h_