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

/*!
   \file frag_natural.cc
   \ingroup optking
   \brief function in frag class to generate natural combinations.  Not yet complete.
*/

#include "frag.h"

//#include "mem.h"
#include "v3d.h"
#include "psi4/psi4-dec.h"
#include "print.h"
#define EXTERN
#include "globals.h"
#include <vector>
#include <algorithm>

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif

#if defined(OPTKING_PACKAGE_PSI)
#include "psi4/libparallel/ParallelPrinter.h"
#endif

namespace opt {

bool int_compare(int i, int j) {return i<j;}

using namespace v3d;

// Determine Pulay simple coordinate combinations.
int FRAG::form_natural_coord_combinations(void) {
  coords.clear_combos();

  int Ns = coords.Nsimples();

  oprintf_out(" Forming natural internal coordinate combinations (Pulay 1992).\n");
  oprintf_out(" Starting with %d simple coordinates.\n", Ns);

  // Make vector list of connected atoms.
  vector< vector<int> > bond;
  for (int i=0; i<natom; ++i) {
    vector<int> one;
    for (int j=0; j<natom; ++j)
      if (connectivity[i][j])
        one.push_back(j);
    bond.push_back(one);
  }

  // Make list of coordination numbers.
  int *CN = init_int_array(natom);
  for (int i=0; i<natom; ++i)
    CN[i] = (int) bond[i].size();

  if (Opt_params.print_lvl > 1) {
    oprintf_out("  -Atom-  -Coord. number-\n");
      for (int i=0; i<natom; ++i)
        oprintf_out("%8d%17d\n", i+1, CN[i]);
  }

  // *** Identify rings ***
  // First, remove terminal atoms from list
  vector< vector<int> > skeleton;
  for (int i=0; i<natom; ++i) {
    if (CN[i] > 1) { // not-terminal
      vector<int> one;
      for (int j=0; j<natom; ++j)
        if (CN[j] > 1) // not-terminal
          one.push_back(j);
      skeleton.push_back(one);
    }
  }

  // Identify rings. TODO: Fix for arbitrary ring sizes
  vector< vector<int> > rings_full;
  //int length;
  for (std::size_t i=0; i<skeleton.size(); ++i) {

    for (std::size_t j=0; j<skeleton[i].size(); ++j) {
      if (j == i) continue;

      for (std::size_t k=0; k<skeleton[j].size(); ++k) {
        if (k == j) continue;
        if (k == i ) {
          vector<int> one;
          one.push_back(i);
          one.push_back(j);
          one.push_back(k);
          rings_full.push_back(one);
        }
        for (std::size_t l=0; l<skeleton[k].size(); ++l) {
          if (l == k || l == j) continue;
          if (l == i ) {
            vector<int> one;
            one.push_back(i);
            one.push_back(j);
            one.push_back(k);
            one.push_back(l);
            rings_full.push_back(one);
          }
          for (std::size_t m=0; m<skeleton[l].size(); ++m) {
            if (m == l || m == k || m == j) continue;
            if (m == i ) {
              vector<int> one;
              one.push_back(i);
              one.push_back(j);
              one.push_back(k);
              one.push_back(l);
              one.push_back(m);
              rings_full.push_back(one);
            }
            for (std::size_t n=0; n<skeleton[m].size(); ++n) {
              if (n == m || n == l || n == k || n == j) continue;
              if (n == i ) {
                vector<int> one;
                one.push_back(i);
                one.push_back(j);
                one.push_back(k);
                one.push_back(l);
                one.push_back(m);
                one.push_back(n);
                rings_full.push_back(one);
              }
            }
          }
        }
      }
    }
  }

  // Put rings in a canonical order; minimize first atom, then second
  for (std::size_t i=0; i<rings_full.size(); ++i) {
    // Find iterator to atom value in ring.
    std::vector<int>::iterator it;
    it = std::min_element(rings_full[i].begin(), rings_full[i].end(), int_compare);

    // Rotate elements until minimum is first value.  Cool!
    while (rings_full[i][0] != *it)
      std::rotate(rings_full[i].begin(), rings_full[i].begin()+1, rings_full[i].end());

    // Now reverse to minimize the second atom
    if (rings_full[i][1] > *rings_full[i].end())
      std::reverse(rings_full[i].begin(), rings_full[i].end());
  }

  // Remove redundant rings.  Make a copy. Check for duplication
  vector< vector<int> > rings;
  for (std::size_t i=0; i<rings.size(); ++i) {
    bool match = false;
    for (std::size_t j=0; j<i; ++j) {
      if (rings_full[i] == rings[j]) { // does this correctly check vector contents?
        match = true;
        break;
      }
    }
    if (!match)
      rings.push_back(rings_full[i]);
  }

  oprintf_out(" %d rings were detected\n", rings.size());
  for (std::size_t i=0; i<rings.size(); ++i) {
    oprintf_out(" Ring %d : \n", i+1);
    for (std::size_t j=0; j<rings.size(); ++j)
      oprintf_out(" %d ", rings[i][j]);
  }

  // Array of ring / non-ring atoms.
  bool *inR = init_bool_array(natom);
  for (int i=0; i<natom; ++i)
    inR[i] = false;

  for (std::size_t i=0; i<rings.size(); ++i)
    for (std::size_t j=0; j<rings[i].size(); ++j)
      inR[ rings[i][j]] = true;

  // Make list of all atoms and the terminal atoms to which they are bonded
  vector< vector<int> > bond_to_T;
  for (int i=0; i<natom; ++i) {
    vector<int> one;
    for (int j=0; j<natom; ++j) {
      if (CN[j] == 1) // terminal
        one.push_back(j);
    }
    bond_to_T.push_back(one);
  }

  // Make list of all atoms and the non-terminal atoms to which they are bonded
  vector< vector<int> > bond_no_T;
  for (int i=0; i<natom; ++i) {
    vector<int> one;
    for (int j=0; j<natom; ++j) {
      if (CN[j] != 1) // terminal
        one.push_back(j);
    }
    bond_no_T.push_back(one);
  }

/*
  - name -  row dimension   - contents -
  bond      all             all
  CN        all             number of bonds

  bond_to_T  all            terminal atoms
  bond_no_T  all            non-terminal atoms

  skeleton  [non-terminal]  non-terminal atoms
  rings     # rings         ring atoms
  inR       all             is atom in ring?
*/
  oprintf_out(" %d simple stretches retained.\n");
  for (std::size_t i=0; i<coords.simples.size(); ++i)
    if (coords.simples[i]->g_type() == stre_type)
      add_trivial_coord_combination(i);

  // Combine angles and dihedrals

  vector<int> cc_index;
  vector<double> cc_coeff;
  bool ok;
  double val;

  // Loop over central atoms
  for (int i=0; i<natom; ++i) {
    if (CN[i] == 1) continue; // terminal atom

    // Number of terminal and non-terminal atoms
    std::size_t num_T = bond_to_T[i].size();
    int num_nonT = CN[i] - num_T;
    oprintf_out("Atom %d has %d non-terminal and %d terminal connections\n", i+1, num_T, num_nonT);

    // This is a ring atom; handle separately
    if (inR[i]) {
      oprintf_out("Ring atoms to be customized later.\n");
    }

    // *** Coordination == 2 cases ***

    // (2, 0) like H2O or CO2; Also (1, 1) like -O-H; Also (0, 2) like C-O-C
    if ( CN[i] == 2 ) {
      // This code mimics FRAG::add_bend_by_connectivity to check for need for linear bend complement
      if (v3d_angle(geom[bond_to_T[i][0]], geom[i], geom[bond_to_T[i][1]], val)) { // computable
        BEND *pA = new BEND( bond_to_T[i][0], i, bond_to_T[i][1]);
        if (val > Opt_params.linear_bend_threshold) // ~175 degrees
          pA->make_lb_normal();
        int A = find(pA);
        cc_index.push_back(A); cc_coeff.push_back(1.0);
        coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
        cc_index.clear(); cc_coeff.clear();

        if (val > Opt_params.linear_bend_threshold) { // ~175 degrees
          BEND *pB = new BEND(bond_to_T[i][0], i, bond_to_T[i][1]);
          pB->make_lb_complement();
          A = find(pB);
          cc_index.push_back(A); cc_coeff.push_back(1.0);
          coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
          cc_index.clear(); cc_coeff.clear();
        }
      }
    }
    // *** Coordination == 3 cases ***
    else if (num_T == 3 && num_nonT == 0) { // like NH3 or BF3
      // There are 3 simple bends for T-C-T.
      BEND *pX1 = new BEND( bond_to_T[i][0], i, bond_to_T[i][1]); // H-C-H, external ones
      int X1 = find(pX1);
      BEND *pX2 = new BEND( bond_to_T[i][0], i, bond_to_T[i][2]);
      int X2 = find(pX2);
      BEND *pX3 = new BEND( bond_to_T[i][1], i, bond_to_T[i][2]);
      int X3 = find(pX3);

      // if planar, add oofp angle
      ok = v3d_oofp(geom[bond_to_T[i][0]], geom[i], geom[bond_to_T[i][1]], geom[bond_to_T[i][2]], val);
      if (ok && fabs(val) < _pi/9) {  // less than 20 degrees?
        OOFP *pA = new OOFP(bond_to_T[i][0], i, bond_to_T[i][1], bond_to_T[i][2]);
        int A = find(pA);
        cc_index.push_back(A); cc_coeff.push_back(1.0);
        coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
        cc_index.clear(); cc_coeff.clear();
      }
      else { // otherwise add symmetric bend
        cc_index.push_back(X1); cc_coeff.push_back( 1.0);
        cc_index.push_back(X2); cc_coeff.push_back( 1.0);
        cc_index.push_back(X3); cc_coeff.push_back( 1.0);
        coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
        cc_index.clear(); cc_coeff.clear();
      }
      // add rock and scissor
      cc_index.push_back(X1); cc_coeff.push_back( 2.0);
      cc_index.push_back(X2); cc_coeff.push_back(-1.0);
      cc_index.push_back(X3); cc_coeff.push_back(-1.0);
      coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
      cc_index.clear(); cc_coeff.clear();

      cc_index.push_back(X2); cc_coeff.push_back( 1.0);
      cc_index.push_back(X3); cc_coeff.push_back(-1.0);
      coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
      cc_index.clear(); cc_coeff.clear();
    }
    else if (num_T == 2 && num_nonT == 1) { // =CH2; very similar to previous case but with a non-terminal atom
      // There is one external angle; and 2 internal angles
      // There are 3 simple bends for T-C-T.
      BEND *pX1 = new BEND( bond_no_T[i][0], i, bond_to_T[i][1]); // external angle H-C-H
      int X1 = find(pX1);
      BEND *pX2 = new BEND( bond_to_T[i][0], i, bond_no_T[i][0]);
      int X2 = find(pX2);
      BEND *pX3 = new BEND( bond_to_T[i][1], i, bond_no_T[i][1]);
      int X3 = find(pX3);

      // if planar, add oofp angle for one of the terminal atoms
      ok = v3d_oofp(geom[bond_to_T[i][0]], geom[i], geom[bond_to_T[i][1]], geom[bond_no_T[i][0]], val);
      if (ok && fabs(val) < _pi/9) {  // less than 20 degrees?
        OOFP *pA = new OOFP(bond_to_T[i][0], i, bond_to_T[i][1], bond_no_T[i][2]);
        int A = find(pA);
        cc_index.push_back(A); cc_coeff.push_back(1.0);
        coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
        cc_index.clear(); cc_coeff.clear();
      }
      else { // otherwise add symmetric bend
        cc_index.push_back(X1); cc_coeff.push_back( 1.0);
        cc_index.push_back(X2); cc_coeff.push_back( 1.0);
        cc_index.push_back(X3); cc_coeff.push_back( 1.0);
        coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
        cc_index.clear(); cc_coeff.clear();
      }
      // add rock and scissor
      cc_index.push_back(X1); cc_coeff.push_back( 2.0);
      cc_index.push_back(X2); cc_coeff.push_back(-1.0);
      cc_index.push_back(X3); cc_coeff.push_back(-1.0);
      coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
      cc_index.clear(); cc_coeff.clear();

      cc_index.push_back(X2); cc_coeff.push_back( 1.0);
      cc_index.push_back(X3); cc_coeff.push_back(-1.0);
      coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
      cc_index.clear(); cc_coeff.clear();
    }
    else if (num_T == 1 && num_nonT == 2) { // Rare?  very similar to previous case but with one terminal atom
      // There are two external angles; and 1 internal angles
      BEND *pX1 = new BEND( bond_to_T[i][0], i, bond_no_T[i][0]); // external
      int X1 = find(pX1);
      BEND *pX2 = new BEND( bond_to_T[i][0], i, bond_no_T[i][1]); // external
      int X2 = find(pX2);
      BEND *pX3 = new BEND( bond_no_T[i][0], i, bond_no_T[i][1]); // internal
      int X3 = find(pX3);

      // if planar, add oofp angle for one of the terminal atoms
      ok = v3d_oofp(geom[bond_to_T[i][0]], geom[i], geom[bond_to_T[i][1]], geom[bond_no_T[i][0]], val);
      if (ok && fabs(val) < _pi/9) {  // less than 20 degrees?
        OOFP *pA = new OOFP(bond_to_T[i][0], i, bond_to_T[i][1], bond_no_T[i][2]);
        int A = find(pA);
        cc_index.push_back(A); cc_coeff.push_back(1.0);
        coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
        cc_index.clear(); cc_coeff.clear();
      }
      else { // otherwise add symmetric bend
        cc_index.push_back(X1); cc_coeff.push_back( 1.0);
        cc_index.push_back(X2); cc_coeff.push_back( 1.0);
        cc_index.push_back(X3); cc_coeff.push_back( 1.0);
        coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
        cc_index.clear(); cc_coeff.clear();
      }
      // add rock and scissor
      cc_index.push_back(X1); cc_coeff.push_back( 2.0);
      cc_index.push_back(X2); cc_coeff.push_back(-1.0);
      cc_index.push_back(X3); cc_coeff.push_back(-1.0);
      coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
      cc_index.clear(); cc_coeff.clear();

      cc_index.push_back(X2); cc_coeff.push_back( 1.0);
      cc_index.push_back(X3); cc_coeff.push_back(-1.0);
      coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
      cc_index.clear(); cc_coeff.clear();
    }


    else if (num_T == 3 && num_nonT > 0) { // like - CH3
      // There are 6 simple bends; 3 external T-C-T and 3 internal.
      BEND *pX1 = new BEND( bond_to_T[i][0], i, bond_to_T[i][1]); // H-C-H, external ones
      int X1 = find(pX1);
      BEND *pX2 = new BEND( bond_to_T[i][0], i, bond_to_T[i][2]);
      int X2 = find(pX2);
      BEND *pX3 = new BEND( bond_to_T[i][1], i, bond_to_T[i][2]);
      int X3 = find(pX3);

      BEND *pI1 = new BEND( bond_to_T[i][0], i, bond_to_T[i][1]); // C-C-H, internal ones
      int I1 = find(pI1);
      BEND *pI2 = new BEND( bond_to_T[i][0], i, bond_to_T[i][2]);
      int I2 = find(pI2);
      BEND *pI3 = new BEND( bond_to_T[i][1], i, bond_to_T[i][2]);
      int I3 = find(pI3);

      // stupid line to stop compiler warnings until this code is done.
      if (I1 == I2 && I2 == I3) return 0;

      // if planar, add oofp angle
      ok = v3d_oofp(geom[bond_to_T[i][0]], geom[i], geom[bond_to_T[i][1]], geom[bond_to_T[i][2]], val);
      if (ok && fabs(val) < _pi/9) {  // less than 20 degrees?
        OOFP *pA = new OOFP(bond_to_T[i][0], i, bond_to_T[i][1], bond_to_T[i][2]);
        int A = find(pA);
        cc_index.push_back(A); cc_coeff.push_back(1.0);
        coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
        cc_index.clear(); cc_coeff.clear();
      }
      else { // otherwise add symmetric bend
        cc_index.push_back(X1); cc_coeff.push_back( 1.0);
        cc_index.push_back(X2); cc_coeff.push_back( 1.0);
        cc_index.push_back(X3); cc_coeff.push_back( 1.0);
        coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
        cc_index.clear(); cc_coeff.clear();
      }
      // add rock and scissor
      cc_index.push_back(X1); cc_coeff.push_back( 2.0);
      cc_index.push_back(X2); cc_coeff.push_back(-1.0);
      cc_index.push_back(X3); cc_coeff.push_back(-1.0);
      coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
      cc_index.clear(); cc_coeff.clear();

      cc_index.push_back(X2); cc_coeff.push_back( 1.0);
      cc_index.push_back(X3); cc_coeff.push_back(-1.0);
      coords.index.push_back(cc_index); coords.coeff.push_back(cc_coeff);
      cc_index.clear(); cc_coeff.clear();



      // choose one other atom to 'anchor' -CH3 group
      if (num_nonT > 0) {



      }
    }




    // X(t0)(t1)(c0)(c1)
    if (num_T == 2 && num_nonT == 2) {                         // CH2

      BEND *pA = new BEND( bond_to_T[i][0], i, bond_to_T[i][1]); // H-C-H
      int A = find(pA);
      BEND *pX1 = new BEND( bond_to_T[i][0], i, bond_no_T[i][0]);  // H-C-C
      int X1 = find(pX1);
      BEND *pX2 = new BEND( bond_to_T[i][0], i, bond_no_T[i][1]);  // H-C-C
      int X2 = find(pX2);
      BEND *pX3 = new BEND( bond_to_T[i][1], i, bond_no_T[i][0]);  // H-C-C
      int X3 = find(pX3);
      BEND *pX4 = new BEND( bond_to_T[i][1], i, bond_no_T[i][1]);  // H-C-C
      int X4 = find(pX4);

      // There are 6 simple bends; making 4 combinations; scissor
      cc_index.push_back(A);  cc_coeff.push_back( 4.0);
      cc_index.push_back(X1); cc_coeff.push_back(-1.0);
      cc_index.push_back(X2); cc_coeff.push_back(-1.0);
      cc_index.push_back(X3); cc_coeff.push_back(-1.0);
      cc_index.push_back(X4); cc_coeff.push_back(-1.0);
      coords.index.push_back(cc_index);
      coords.coeff.push_back(cc_coeff);
      cc_index.clear();
      cc_coeff.clear();

      // now rock, wag and twist
      cc_index.push_back(X1); cc_coeff.push_back( 1.0);
      cc_index.push_back(X2); cc_coeff.push_back( 1.0);
      cc_index.push_back(X3); cc_coeff.push_back(-1.0);
      cc_index.push_back(X4); cc_coeff.push_back(-1.0);
      coords.index.push_back(cc_index);
      coords.coeff.push_back(cc_coeff);
      cc_index.clear();
      cc_coeff.clear();

      cc_index.push_back(X1); cc_coeff.push_back( 1.0);
      cc_index.push_back(X2); cc_coeff.push_back(-1.0);
      cc_index.push_back(X3); cc_coeff.push_back( 1.0);
      cc_index.push_back(X4); cc_coeff.push_back(-1.0);
      coords.index.push_back(cc_index);
      coords.coeff.push_back(cc_coeff);
      cc_index.clear();
      cc_coeff.clear();

      cc_index.push_back(X1); cc_coeff.push_back( 1.0);
      cc_index.push_back(X2); cc_coeff.push_back(-1.0);
      cc_index.push_back(X3); cc_coeff.push_back(-1.0);
      cc_index.push_back(X4); cc_coeff.push_back( 1.0);
      coords.index.push_back(cc_index);
      coords.coeff.push_back(cc_coeff);
      cc_index.clear();
      cc_coeff.clear();

      // I think you omit this one?  bond_no_T[i][0], i, bond_no_T[1])  C-C-C
    }

  }
  return 0;
}

//void add(int a, int b, int c) {
  //combos.coords.
//}

}
