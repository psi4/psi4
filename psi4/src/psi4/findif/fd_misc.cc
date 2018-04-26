/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file util.cc
    \ingroup OPTKING
    \brief miscellaneous
*/

#include "psi4/findif/findif.h"

#include "psi4/physconst.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/cdsalclist.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/writer_file_prefix.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libmints/wavefunction.h"

#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
using SharedMatrix=std::shared_ptr<Matrix>;
namespace findif {

// displaces from a reference geometry: geom += salclist[salc_i] * disp_i * disp_size
// disp_size is in mass-weighted coordinates; cartesian displacement is DX/sqrt(mass)
void displace_cart(std::shared_ptr<Molecule> mol, SharedMatrix geom, const CdSalcList & salclist,
  int salc_i, int disp_factor, double disp_size) {

  geom->set_name("Coord: " + to_string(salc_i) + ", Disp: " + to_string(disp_factor));

  int nc = salclist[salc_i].ncomponent();

  for (int c=0; c<nc; ++c) {
    int a          = salclist[salc_i].component(c).atom;
    int xyz        = salclist[salc_i].component(c).xyz;
    double coef    = salclist[salc_i].component(c).coef;

    geom->add(0, a, xyz, disp_factor * disp_size * coef / sqrt(mol->mass(a)));
  }

  return;
}

// displaces from a reference geometry.
// geom += salclist[salc_i] * disp_i * disp_size + salclist[salc_j] * disp_j * disp_size
// disp_size is in mass-weighted coordinates; cartesian displacement is DX/sqrt(mass)
void displace_cart(std::shared_ptr<Molecule> mol, SharedMatrix geom, const CdSalcList & salclist,
  int salc_i, int salc_j, int disp_factor_i, int disp_factor_j, double disp_size) {

  geom->set_name("Coord: " + to_string(salc_i) + ", Disp: " + to_string(disp_factor_i)
    + "Coord: " + to_string(salc_j) + ", Disp: " + to_string(disp_factor_j));

  int a, xyz;
  double coef;

  for (int c=0; c<salclist[salc_i].ncomponent(); ++c) {
    a    = salclist[salc_i].component(c).atom;
    xyz  = salclist[salc_i].component(c).xyz;
    coef = salclist[salc_i].component(c).coef;

    geom->add(0, a, xyz, disp_factor_i * disp_size * coef / sqrt(mol->mass(a)));
  }

  for (int c=0; c<salclist[salc_j].ncomponent(); ++c) {
    a    = salclist[salc_j].component(c).atom;
    xyz  = salclist[salc_j].component(c).xyz;
    coef = salclist[salc_j].component(c).coef;

    geom->add(0, a, xyz, disp_factor_j * disp_size * coef / sqrt(mol->mass(a)));
  }

  return;
}

// it's assumed columns are cartesian dimensions
void mass_weight_columns_plus_one_half(std::shared_ptr<Molecule> mol, SharedMatrix B) {
  double u;

  for (int col=0; col<B->ncol(); ++col) {
    u = sqrt(mol->mass(col/3));
    for (int row=0; row<B->nrow(); ++row)
      B->set(row, col, B->get(row,col) * u);
  }
}

void displace_atom(SharedMatrix geom, const int atom, const int coord, const int sign, const double disp_size) {

  geom->add(0, atom, coord, sign * disp_size);

  return;
}

std::vector< SharedMatrix > atomic_displacements(std::shared_ptr<Molecule> mol, Options &options) {

  // This is the size in bohr because geometry is in bohr at this point
  double disp_size = options.get_double("DISP_SIZE");
  int pts = options.get_int("POINTS");

  int natom = mol->natom();

  Matrix ref_geom_temp = mol->geometry();
  SharedMatrix ref_geom(ref_geom_temp.clone());

  std::vector< SharedMatrix > disp_geoms;

  // Generate displacements
  if (pts == 3) {
    for(int atom=0; atom < natom; ++atom) {
      for(int coord=0; coord < 3; ++coord) {
        // minus displacement
        SharedMatrix m_geom(ref_geom->clone());
        displace_atom(m_geom, atom, coord, -1, disp_size);
        disp_geoms.push_back(m_geom);
        // plus displacement
        SharedMatrix p_geom(ref_geom->clone());
        displace_atom(p_geom, atom, coord, +1, disp_size);
        disp_geoms.push_back(p_geom);
      }
    }
  }
  else if (pts == 5) {
    for(int atom=0; atom < natom; ++atom) {
      for(int coord=0; coord < 3; ++coord) {
        // minus displacements
        SharedMatrix m2_geom(ref_geom->clone());
        displace_atom(m2_geom, atom, coord, -2, disp_size);
        disp_geoms.push_back(m2_geom);
        SharedMatrix m1_geom(ref_geom->clone());
        displace_atom(m1_geom, atom, coord, -1, disp_size);
        disp_geoms.push_back(m1_geom);
        // plus displacements
        SharedMatrix p1_geom(ref_geom->clone());
        displace_atom(p1_geom, atom, coord, +1, disp_size);
        disp_geoms.push_back(p1_geom);
        SharedMatrix p2_geom(ref_geom->clone());
        displace_atom(p2_geom, atom, coord, +2, disp_size);
        disp_geoms.push_back(p2_geom);
      }
    }
  }
  else {
    throw PsiException("FINDIF: Number of POINTS not supported", __FILE__, __LINE__);
  }
  return disp_geoms;
}

}}
