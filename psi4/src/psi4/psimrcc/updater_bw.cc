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
//#include "mrcc.h"
//#include "matrix.h"
//#include "debugging.h"
//#include "psi4/libpsi4util/libpsi4util.h"

#include "blas.h"
#include "heff.h"
#include "updater.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

BWUpdater::BWUpdater(Options &options) :
        Updater(options)
{
}

BWUpdater::~BWUpdater()
{
}

void BWUpdater::update(int cycle,Hamiltonian* heff)
{
  blas->solve("d'2[oo][vv]{u}  = d2[oo][vv]{u}");
  blas->solve("d'2[oO][vV]{u}  = d2[oO][vV]{u}");
  blas->solve("d'2[OO][VV]{u}  = d2[OO][VV]{u}");

  for(int mu = 0; mu < moinfo->get_nunique(); ++mu){
    int mu_unique = moinfo->get_ref_number(mu,UniqueRefs);
    double denominator_shift = heff->get_eigenvalue() - heff->get_matrix(mu_unique,mu_unique);
    std::string shift = to_string(denominator_shift);
    std::string mu_str = to_string(mu_unique);
    blas->solve("d'2[oo][vv]{" + mu_str + "} += " + shift);
    blas->solve("d'2[oO][vV]{" + mu_str + "} += " + shift);
    blas->solve("d'2[OO][VV]{" + mu_str + "} += " + shift);
  }

  // (a) Compute eq. (20) of J. Chem. Phys. 110, 10275 (1999)
  // Comment : Look at eq. (21) of J. Chem. Phys. 110, 10275 (1999)
  blas->solve("t1_eqns[o][v]{u} += - d1[o][v]{u} * t1[o][v]{u}");
  blas->solve("t1_eqns[O][V]{u} += - d1[O][V]{u} * t1[O][V]{u}");

  // aaaa case
  // (b) Add PijPab (term from a) to the T2 equations
  blas->solve("t2_eqns[oo][vv]{u} += #1324#   t1[o][v]{u} X t1_eqns[o][v]{u}");
  blas->solve("t2_eqns[oo][vv]{u} += #2314# - t1[o][v]{u} X t1_eqns[o][v]{u}");
  blas->solve("t2_eqns[oo][vv]{u} += #1423# - t1[o][v]{u} X t1_eqns[o][v]{u}");
  blas->solve("t2_eqns[oo][vv]{u} += #2413#   t1[o][v]{u} X t1_eqns[o][v]{u}");
  // (c) Subtract (term from c) from the T2 equations
  for(int mu = 0; mu < moinfo->get_nunique(); ++mu){
    int mu_unique = moinfo->get_ref_number(mu,UniqueRefs);
    double denominator_shift = heff->get_eigenvalue() - heff->get_matrix(mu_unique,mu_unique);
    std::string neg_shift = to_string(-denominator_shift);
    std::string shift     = to_string(denominator_shift);
    std::string mu_str    = to_string(mu_unique);
    blas->solve("t2_eqns[oo][vv]{" + mu_str + "} += #1324# " + neg_shift + "  t1[o][v]{" + mu_str + "} X t1[o][v]{" + mu_str + "}");
    blas->solve("t2_eqns[oo][vv]{" + mu_str + "} += #2314# " + shift  + "  t1[o][v]{" + mu_str + "} X t1[o][v]{" + mu_str + "}");
  }

  // abab case
  // (b) Add PijPab (term from a) to the T2 equations
  blas->solve("t2_eqns[oO][vV]{u} += #1324# t1[o][v]{u} X t1_eqns[O][V]{u}");
  blas->solve("t2_eqns[oO][vV]{u} += #2413# t1[O][V]{u} X t1_eqns[o][v]{u}");
  // (c) Subtract (term from c) from the T2 equations
  for(int mu = 0; mu < moinfo->get_nunique(); ++mu){
      int mu_unique = moinfo->get_ref_number(mu,UniqueRefs);
      std::string mu_str    = to_string(mu_unique);
      double denominator_shift = heff->get_eigenvalue() - heff->get_matrix(mu_unique,mu_unique);
      std::string neg_shift = to_string(-denominator_shift);
      blas->solve("t2_eqns[oO][vV]{" + mu_str + "} += #1324# " + neg_shift + "  t1[o][v]{" + mu_str + "} X t1[O][V]{" + mu_str + "}");
    }

  // bbbb case
  // (b) Add PijPab (term from a) to the T2 equations
  blas->solve("t2_eqns[OO][VV]{u} += #1324#   t1[O][V]{u} X t1_eqns[O][V]{u}");
  blas->solve("t2_eqns[OO][VV]{u} += #2314# - t1[O][V]{u} X t1_eqns[O][V]{u}");
  blas->solve("t2_eqns[OO][VV]{u} += #1423# - t1[O][V]{u} X t1_eqns[O][V]{u}");
  blas->solve("t2_eqns[OO][VV]{u} += #2413#   t1[O][V]{u} X t1_eqns[O][V]{u}");
  // (c) Subtract (term from c) from the T2 equations
  for(int mu = 0; mu < moinfo->get_nunique(); ++mu){
    int mu_unique = moinfo->get_ref_number(mu,UniqueRefs);
    double denominator_shift = heff->get_eigenvalue() - heff->get_matrix(mu_unique,mu_unique);
    std::string neg_shift = to_string(-denominator_shift);
    std::string shift     = to_string(denominator_shift);
    std::string mu_str    = to_string(mu_unique);
    blas->solve("t2_eqns[OO][VV]{" + mu_str + "} += #1324# " + neg_shift + "  t1[O][V]{" + mu_str + "} X t1[O][V]{" + mu_str + "}");
    blas->solve("t2_eqns[OO][VV]{" + mu_str + "} += #2314# " + shift     + "  t1[O][V]{" + mu_str + "} X t1[O][V]{" + mu_str + "}");
  }

  blas->solve("t2_delta[oo][vv]{u} = t2_eqns[oo][vv]{u} / d'2[oo][vv]{u} - t2[oo][vv]{u}");
  blas->solve("t2_delta[oO][vV]{u} = t2_eqns[oO][vV]{u} / d'2[oO][vV]{u} - t2[oO][vV]{u}");
  blas->solve("t2_delta[OO][VV]{u} = t2_eqns[OO][VV]{u} / d'2[OO][VV]{u} - t2[OO][VV]{u}");

  blas->solve("t2[oo][vv]{u} = t2_eqns[oo][vv]{u} / d'2[oo][vv]{u}");
  blas->solve("t2[oO][vV]{u} = t2_eqns[oO][vV]{u} / d'2[oO][vV]{u}");
  blas->solve("t2[OO][VV]{u} = t2_eqns[OO][VV]{u} / d'2[OO][VV]{u}");

//  DEBUGGING(3,
//    blas->print("t2_eqns[oo][vv]{u}");
//    blas->print("t2[oo][vv]{u}");
//    blas->print("t2_eqns[oO][vV]{u}");
//    blas->print("t2[oO][vV]{u}");
//    blas->print("t2_eqns[OO][VV]{u}");
//    blas->print("t2[OO][VV]{u}");
//  );

  blas->solve("d'1[o][v]{u}  = d1[o][v]{u}");
  blas->solve("d'1[O][V]{u}  = d1[O][V]{u}");

  for(int mu = 0; mu < moinfo->get_nunique(); ++mu){
    int mu_unique = moinfo->get_ref_number(mu,UniqueRefs);
    double denominator_shift = heff->get_eigenvalue() - heff->get_matrix(mu_unique,mu_unique);
    std::string shift     = to_string(denominator_shift);
    std::string mu_str    = to_string(mu_unique);
    blas->solve("d'1[o][v]{" + mu_str + "} += " + shift);
    blas->solve("d'1[O][V]{" + mu_str + "} += " + shift);
  }

  blas->solve("t1_delta[o][v]{u}  =   t1_eqns[o][v]{u} / d'1[o][v]{u} - t1[o][v]{u}");
  blas->solve("t1_delta[O][V]{u}  =   t1_eqns[O][V]{u} / d'1[O][V]{u} - t1[O][V]{u}");

  blas->solve("t1[o][v]{u} = t1_eqns[o][v]{u} / d'1[o][v]{u}");
  blas->solve("t1[O][V]{u} = t1_eqns[O][V]{u} / d'1[O][V]{u}");

  blas->solve("t1_norm{u}  = t1[o][v]{u} . t1[o][v]{u}");
  blas->solve("t1_norm{u} += t1[O][V]{u} . t1[O][V]{u}");

  zero_t1_internal_amps();

  zero_internal_delta_amps();

  blas->solve("||Delta_t1||{u}  = t1_delta[o][v]{u} . t1_delta[o][v]{u}");
  blas->solve("||Delta_t1||{u} += t1_delta[O][V]{u} . t1_delta[O][V]{u}");

  blas->solve("||Delta_t2||{u}  = t2_delta[oo][vv]{u} . t2_delta[oo][vv]{u}");
  blas->solve("||Delta_t2||{u} += t2_delta[oO][vV]{u} . t2_delta[oO][vV]{u}");
  blas->solve("||Delta_t2||{u} += t2_delta[OO][VV]{u} . t2_delta[OO][VV]{u}");


//    DEBUGGING(3,
//      blas->print("t1[o][v]{u}");
//      blas->print("t1[O][V]{u}");
//    );

}

}} /* End Namespaces */
