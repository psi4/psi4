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

#ifndef _psi_src_bin_psimrcc_ccmanybody_h
#define _psi_src_bin_psimrcc_ccmanybody_h
/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

/*********************************************************
  CCManyBody Class
  1) Purpose
    This class is used to do the basic operations that any
    manybody code requires: compute the Fock matrix, denominators.
    However, this implementation is specifically designed to
    handle multireference cases
  2) Use
  3) Details
  4) Uses
    MOInfo class
    STL <vector>

*********************************************************/

#include <cmath>
#include <vector>
#include <string>

#include "psi4/liboptions/liboptions.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"

namespace psi{ namespace psimrcc{


enum SpinCase             {aaSpin,abSpin,bbSpin,aaaSpin,aabSpin,abbSpin,bbbSpin};
enum TriplesType          {pt2,ccsd,ccsd_t,ccsdt_1a,ccsdt_1b,ccsdt_2,ccsdt_3,ccsdt};
enum TriplesCouplingType  {nocoupling,linear,quadratic,cubic};


/**
 *	@author Francesco Evangelista <frank@ccc.uga.edu>
*/
class CCManyBody{
public:
  CCManyBody(SharedWavefunction ref_wfn, Options &options);
  virtual ~CCManyBody();
  void        generate_integrals();
  void        generate_denominators();
  void        compute_reference_energy();
  void        make_fock_matrix();
  void        make_denominators();
  void        print_method(const char* text);
//  void        zero_internal_amps();
//  void        zero_t1_internal_amps();
//  void        zero_internal_delta_amps();
protected:
  Options &options_;
  SharedWavefunction ref_wfn_;
  // Effective Hamiltonian and the correpsonding eigenvectors
  void        print_eigensystem(int ndets, double** Heff,double*& eigenvector);
  double      diagonalize_Heff(int root,int ndets, double** Heff,double*& right_eigenvector,double*& left_eigenvector, bool initial);
  void        sort_eigensystem(int ndets,double*& real,double*& imaginary,double**& left,double**& right);
  double      c_H_c(int ndets, double** H,double*& c);

  double*     zeroth_order_eigenvector;
  double*     right_eigenvector;
  double*     left_eigenvector;
  double**    Heff;
  double**    Heff_mrpt2;

  // Effective Hamiltonian and the correpsonding eigenvectors
  double      current_energy;
  double      delta_energy;
  double      cas_energy;
  double      old_energy;

  double      huge;

  double      total_time;

  double      norm_amps;
  double      delta_t1_amps;
  double      delta_t2_amps;

  bool        pert_cbs;
  bool        pert_cbs_coupling;
  TriplesType         triples_type;
  TriplesCouplingType triples_coupling_type;

  void generate_triples_denominators();
  void generate_d3_ijk(double***& d3,bool alpha_i,bool alpha_j,bool alpha_k);
  void generate_d3_abc(double***& d3,bool alpha_a,bool alpha_b,bool alpha_c);
  void deallocate_triples_denominators();

  double***  d3_ooo;
  double***  d3_ooO;
  double***  d3_oOO;
  double***  d3_OOO;
  double***  d3_vvv;
  double***  d3_vvV;
  double***  d3_vVV;
  double***  d3_VVV;
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_ccmanybody_h
