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

#include "heff.h"

#include <algorithm>
#include <functional>
#include <utility>

namespace psi{

    namespace psimrcc{
    extern MOInfo *moinfo;

Hamiltonian::Hamiltonian()
{
  startup();
}

Hamiltonian::~Hamiltonian()
{
  cleanup();
}

void Hamiltonian::startup()
{
}

void Hamiltonian::cleanup()
{
}

void Hamiltonian::print_matrix()
{
  if(ndets < 8){
    outfile->Printf("\n\n  Hamiltonian Matrix\n");
    for(int mu = 0; mu < ndets; ++mu){
      outfile->Printf("\n  ");
      for(int nu = 0; nu < ndets; ++nu)
        outfile->Printf(" %22.15f",matrix[mu][nu]);
    }
  }
}

void Hamiltonian::print()
{
  print_matrix();

  std::vector<std::pair<double,int> > eigenvector_index_pair;
  for(int mu = 0; mu < ndets; ++mu){
    eigenvector_index_pair.push_back(std::make_pair(right_eigenvector[mu] * right_eigenvector[mu],mu));
  }
  std::sort(eigenvector_index_pair.begin(),eigenvector_index_pair.end(),std::greater<std::pair<double,int> >());
  int max_size_list = std::min(10,static_cast<int>(eigenvector_index_pair.size()));
  outfile->Printf("\n\n  Most important determinants in the wave function");
  outfile->Printf("\n\n  determinant  eigenvector   eigenvector^2\n");
  for(int i = 0; i < max_size_list; ++i){
    outfile->Printf("\n  %11d   %9.6f    %9.6f  %s",eigenvector_index_pair[i].second
                     ,right_eigenvector[eigenvector_index_pair[i].second]
                     ,eigenvector_index_pair[i].first
                     ,moinfo->get_determinant_label(eigenvector_index_pair[i].second).c_str());
  }
}

void Hamiltonian::set_matrix(double** M,int ndets_)
{
  ndets = ndets_;

  matrix.clear();
  for(int mu = 0; mu < ndets; ++mu){
    std::vector<double> row(ndets,0.0);
    matrix.push_back(row);
  }

  for(int mu = 0; mu < ndets; ++mu){
    for(int nu = 0; nu < ndets; ++nu){
      matrix[mu][nu] = M[mu][nu];
    }
  }
}

void Hamiltonian::set_left_eigenvector(double* v,int ndets_)
{
  ndets = ndets_;
  left_eigenvector.assign(ndets,0.0);
  for(int mu = 0; mu < ndets; ++mu){
    left_eigenvector[mu] = v[mu];
  }
}

void Hamiltonian::set_right_eigenvector(double* v,int ndets_)
{
  ndets = ndets_;
  right_eigenvector.assign(ndets,0.0);
  for(int mu = 0; mu < ndets; ++mu){
    right_eigenvector[mu] = v[mu];
  }
}

void Hamiltonian::set_zeroth_order_eigenvector(double* v,int ndets_)
{
  ndets = ndets_;
  zeroth_order_eigenvector.assign(ndets,0.0);
  for(int mu = 0; mu < ndets; ++mu){
    zeroth_order_eigenvector[mu] = v[mu];
  }
}

double Hamiltonian::expectation_value()
{
  double value = 0.0;
  for(int mu = 0; mu < ndets; ++mu){
    for(int nu = 0; nu < ndets; ++nu){
      value += left_eigenvector[mu] * matrix[mu][nu] * right_eigenvector[nu];
    }
  }
  return(value);
}

double Hamiltonian::trace()
{
  double value = 0.0;
  for(int mu = 0; mu < ndets; ++mu){
      value += left_eigenvector[mu] * matrix[mu][mu] * right_eigenvector[mu];
  }
  return(value);
}


}} /* End Namespace */
