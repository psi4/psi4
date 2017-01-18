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

#include <cmath>
#include <cstdio>
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/memory_manager.h"

#include "algebra_interface.h"
#include "heff.h"

#include <algorithm>
#include <functional>
#include <utility>

namespace psi{

    namespace psimrcc{
    extern MemoryManager* memory_manager;


void sort_eigensystem(int ndets,double*& real,double*& imaginary,double**& left,double**& right);

double Hamiltonian::diagonalize(int root)
{
  double      energy;
  double*     real;
  double*     imaginary;
  double*     work;
  double**    left;
  double**    right;
  double**    H;

  int lwork = 6 * ndets * ndets;
  allocate1(double,work,lwork);
  allocate1(double,real,ndets);
  allocate1(double,imaginary,ndets);

  allocate2(double,H,ndets,ndets);
  allocate2(double,left,ndets,ndets);
  allocate2(double,right,ndets,ndets);

  for(int i=0;i<ndets;i++)
    for(int j=0;j<ndets;j++)
      H[j][i] = matrix[i][j];

  int info;

  F_DGEEV("V","V",&ndets, &(H[0][0]), &ndets, &(real[0]), &(imaginary[0]),
      &(left[0][0]), &ndets, &(right[0][0]), &ndets, &(work[0]), &lwork, &info);

  sort_eigensystem(ndets,real,imaginary,left,right);

//  if(initial){
//    if(ndets < 8){
//      outfile->Printf("\n\n  Heff Matrix\n");
//      for(int i=0;i<ndets;i++){
//        outfile->Printf("\n  ");
//        for(int j=0;j<ndets;j++)
//          outfile->Printf(" %22.12f",matrix[i][j]);
//      }
//
//      outfile->Printf("\n\n  Left Matrix\n");
//      for(int i=0;i<ndets;i++){
//        outfile->Printf("\n  ");
//        for(int j=0;j<ndets;j++)
//          outfile->Printf(" %22.12f",left[j][i]);
//      }
//
//      outfile->Printf("\n\n  Right Matrix\n");
//      for(int i=0;i<ndets;i++){
//        outfile->Printf("\n  ");
//        for(int j=0;j<ndets;j++)
//          outfile->Printf(" %22.12f",right[j][i]);
//      }
//
//      outfile->Printf("\n\n  Real                  Imaginary\n");
//      for(int i=0;i<ndets;i++)
//        outfile->Printf("\n  %22.12f   %22.12f",real[i],imaginary[i]);
//      outfile->Printf("\n");
//    }else{
//      outfile->Printf("\n\n  There are too many determinants to print the eigensystem");
//    }
//    outfile->Printf("\n\n  The eigenvalue for root %d is %.12f (%.12f)",root,real[root],imaginary[root]);
//  }

  bool initial = false;

  if(right_eigenvector.size() != static_cast<size_t>(ndets)){
    initial = true;
    right_eigenvector.assign(ndets,0.0);
    left_eigenvector.assign(ndets,0.0);
  }else{
    double norm = 0.0;
    for(int k = 0; k < ndets; ++k){
      norm += right_eigenvector[k] * right_eigenvector[k];
    }
    if(norm < 0.01){
      initial = true;
    }
  }

  // Select the eigenvector to follow
  if(initial){
    for(int k = 0; k < ndets; ++k){
      right_eigenvector[k] = right[root][k];
      left_eigenvector[k]  =  left[root][k];
    }
    energy = real[root];
  }
  else // find vector with maximum overlap
  {
    int    select_vect  = 0;
    double max_overlap  = 0.0;
    for(int i=0;i<ndets;i++){
      double overlap=0.0;
      for(int m=0;m<ndets;m++)
        overlap += right_eigenvector[m] * right[i][m];
      overlap = sqrt(overlap*overlap);
      if(overlap > max_overlap){
        select_vect = i;
        max_overlap = overlap;
      }
    }
    for(int m=0;m<ndets;m++){
      right_eigenvector[m] = right[select_vect][m];
      left_eigenvector[m]  =  left[select_vect][m];
    }
    energy = real[select_vect];
    root   = select_vect;
  }

  // Normalize the left-eigenvector to <L|R> = 1
  double lnorm = 0.0;
  for(int m = 0; m < ndets; m++){
    lnorm += right_eigenvector[m] * left_eigenvector[m];
  }

  for(int m = 0; m < ndets; m++){
    left_eigenvector[m] = left_eigenvector[m] / lnorm;
  }


  release1(work);
  release1(real);
  release1(imaginary);
  release2(H);
  release2(left);
  release2(right);
  return(energy);
}

void sort_eigensystem(int ndets,double*& real,double*& imaginary,double**& left,double**& right)
{
  std::vector<std::pair<double, int> > pairs;
  for(int i=0;i<ndets;i++)
      pairs.push_back(std::make_pair(real[i],i));
  sort(pairs.begin(),pairs.end());

  double*  tempv;
  double** tempm;
  allocate1(double,tempv,ndets);
  allocate2(double,tempm,ndets,ndets);

  for(int i=0;i<ndets;i++) tempv[i] = real[pairs[i].second];
  for(int i=0;i<ndets;i++) real[i]  = tempv[i];

  for(int i=0;i<ndets;i++) tempv[i]     = imaginary[pairs[i].second];
  for(int i=0;i<ndets;i++) imaginary[i] = tempv[i];

  for(int i=0;i<ndets;i++)
    for(int j=0;j<ndets;j++)
      tempm[i][j] = left[pairs[i].second][j];
  for(int i=0;i<ndets;i++)
    for(int j=0;j<ndets;j++)
      left[i][j] = tempm[i][j];

  for(int i=0;i<ndets;i++)
    for(int j=0;j<ndets;j++)
      tempm[i][j] = right[pairs[i].second][j];
  for(int i=0;i<ndets;i++)
    for(int j=0;j<ndets;j++)
      right[i][j] = tempm[i][j];

  release1(tempv);
  release2(tempm);
}

}} /* End Namespaces */
