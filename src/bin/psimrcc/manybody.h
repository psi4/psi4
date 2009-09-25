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

#include <libciomr/libciomr.h>
#include <libqt/qt.h>

namespace psi{ namespace psimrcc{


enum SpinCase             {aaSpin,abSpin,bbSpin,aaaSpin,aabSpin,abbSpin,bbbSpin};
enum TriplesType          {pt2,ccsd,ccsd_t,ccsdt_1a,ccsdt_1b,ccsdt_2,ccsdt_3,ccsdt};
enum TriplesCouplingType  {nocoupling,linear,quadratic,cubic};


/**
 *	@author Francesco Evangelista <frank@ccc.uga.edu>
*/
class CCManyBody{
public:
  CCManyBody();
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
