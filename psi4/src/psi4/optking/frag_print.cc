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
   \file frag_print.cc
   \ingroup optking
   \brief FRAG member functions for printing and string definitions.
*/

#include "frag.h"

#include "mem.h"
#include "v3d.h"
#include "atom_data.h"
#include "cov_radii.h"
#include "opt_data.h"
#include "psi4/optking/physconst.h"
#include "linear_algebra.h"
#include "psi4/psi4-dec.h"
#include "print.h"
#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif
#include "psi4/libparallel/ParallelPrinter.h"
namespace opt {

void FRAG::print_geom(std::string psi_fp, FILE *qc_fp, const int id, bool print_masses) {
  int i;
  oprintf(psi_fp, qc_fp, "\t---Fragment %d Geometry---\n", id+1);
  if (print_masses) {
    for (i=0; i<natom; ++i)
      oprintf(psi_fp, qc_fp, "\t %-4s%20.10lf%20.10lf%20.10lf%20.10lf\n",
        Z_to_symbol[(int) Z[i]], mass[i], geom[i][0], geom[i][1], geom[i][2]);
  }
  else {
    for (i=0; i<natom; ++i)
      oprintf(psi_fp, qc_fp, "\t %-4s%20.10lf%20.10lf%20.10lf\n",
        Z_to_symbol[(int) Z[i]], geom[i][0], geom[i][1], geom[i][2]);
  }
  oprintf(psi_fp, qc_fp,  "\n");
}

void FRAG::print_geom_grad(std::string psi_fp, FILE *qc_fp, const int id, bool print_masses) {
  int i;
  oprintf(psi_fp, qc_fp, "\t---Fragment %d Geometry and Gradient---\n", id+1);
  if (print_masses) {
    for (i=0; i<natom; ++i)
      oprintf(psi_fp, qc_fp, "\t %-4s%20.10lf%20.10lf%20.10lf%20.10lf\n",
        Z_to_symbol[(int) Z[i]], mass[i], geom[i][0], geom[i][1], geom[i][2]);
  }
  else {
    for (i=0; i<natom; ++i)
      oprintf(psi_fp, qc_fp, "\t %-4s%20.10lf%20.10lf%20.10lf\n",
        Z_to_symbol[(int) Z[i]], geom[i][0], geom[i][1], geom[i][2]);
  }
  for (i=0; i<natom; ++i)
    oprintf(psi_fp, qc_fp, "\t %24.10lf%20.10lf%20.10lf\n", grad[i][0], grad[i][1], grad[i][2]);
  oprintf(psi_fp, qc_fp,  "\n");
}

#if defined(OPTKING_PACKAGE_QCHEM)
void FRAG::print_geom(std::string psi_fp, FILE *qc_fp) {
  for (int i=0; i<natom; ++i)
    oprintf(psi_fp, qc_fp, "\t  %3s  %15.10lf%15.10lf%15.10lf\n",
      Z_to_symbol[(int) Z[i]], geom[i][0], geom[i][1], geom[i][2]);
  fflush(fp_geom);
}

void FRAG::print_geom_irc(std::string psi_fp, FILE *qc_fp) {
  for (int i=0; i<natom; ++i)
    oprintf(psi_fp, qc_fp, "@IRC     %3s  %15.10lf%15.10lf%15.10lf\n",
      Z_to_symbol[(int) Z[i]], geom[i][0], geom[i][1], geom[i][2]);
  fflush(fp_geom);
}
#elif defined(OPTKING_PACKAGE_PSI)
void FRAG::print_geom(std::string psi_fp, FILE *qc_fp) {
   for (int i=0; i<natom; ++i)
    oprintf(psi_fp, qc_fp, "\t  %3s  %15.10lf%15.10lf%15.10lf\n",
      Z_to_symbol[(int) Z[i]], geom[i][0] * _bohr2angstroms,
      geom[i][1] * _bohr2angstroms, geom[i][2] * _bohr2angstroms);
}

void FRAG::print_geom_irc(std::string psi_fp, FILE *qc_fp) {
   for (int i=0; i<natom; ++i)
    oprintf(psi_fp, qc_fp, "@IRC     %3s  %15.10lf%15.10lf%15.10lf\n",
      Z_to_symbol[(int) Z[i]], geom[i][0] * _bohr2angstroms,
      geom[i][1] * _bohr2angstroms, geom[i][2] * _bohr2angstroms);
}
#endif

void FRAG::print_simples(std::string psi_fp, FILE *qc_fp, int atom_offset) const {
  oprintf(psi_fp, qc_fp, "\t - Coordinate -           - BOHR/RAD -       - ANG/DEG -\n");
  for (std::size_t i=0; i<coords.simples.size(); ++i)
    coords.simples.at(i)->print(psi_fp,qc_fp,geom,atom_offset);
  oprintf(psi_fp, qc_fp,  "\n");
}

// obsolete?  not good for delocalized
void FRAG::print_coords(std::string psi_fp, FILE *qc_fp, int atom_offset) const {
  oprintf(psi_fp, qc_fp, "\t - Coordinate -           - BOHR/RAD -       - ANG/DEG -\n");
  for (int i=0; i<Ncoord(); ++i)
    coords.print(psi_fp, qc_fp, i, geom, atom_offset);
  oprintf(psi_fp, qc_fp,  "\n");
}

void FRAG::print_combinations(std::string psi_fp, FILE *qc_fp) const {
  oprintf(psi_fp, qc_fp, "\t-- Internal Coordinate Combinations\n");
  for (int cc=0; cc<Ncoord(); ++cc) {
    oprintf(psi_fp, qc_fp, " Coord %d:\n", cc+1);
    int cnt = 0;
    for (std::size_t s=0; s<coords.index[cc].size(); ++s) {
      oprintf(psi_fp, qc_fp, "%5d:%12.6f", coords.index[cc][s]+1, coords.coeff[cc][s]);
      ++cnt;
      if (cnt == 4) {
        oprintf(psi_fp, qc_fp,"\n");
        cnt = 0;
      }
    }
    if (cnt != 0) oprintf(psi_fp, qc_fp,"\n");
  }
}

// fetch string definition of intco; often for error message reporting
std::string FRAG::get_coord_definition(int coord_index, int atom_offset) {
  //oprintf_out("\tCoordinate index %d; Atom offset %d\n", coord_index, atom_offset);
  //oprintf_out("\t%15s%15s\n", "Coordinate", "Coefficient");
  //for (std::size_t s=0; s<coords.index[coord_index].size(); ++s)
    //oprintf_out("\t%15d%15.5lf", coords.index[coord_index][s]+1, coords.coeff[coord_index][s]);
  return coords.get_coord_definition(coord_index, atom_offset);
}

// fetch string definition of simple
std::string FRAG::get_simple_definition(int simple_index, int atom_offset) {
  oprintf_out("simple_index: %d; atom_offset: %d\n", simple_index, atom_offset);
  return coords.simples.at(simple_index)->get_definition_string(atom_offset);
}

void FRAG::print_intco_dat(std::string psi_fp, FILE *qc_fp, int atom_offset) const {
  for (std::size_t i=0; i<coords.simples.size(); ++i)
    coords.simples.at(i)->print_intco_dat(psi_fp,qc_fp,atom_offset);
  for (std::size_t cc=0; cc<coords.index.size(); ++cc) {
    oprintf(psi_fp, qc_fp, "C %6d\n", coords.index[cc].size());
    for (std::size_t s=0; s<coords.index[cc].size(); ++s)
      oprintf(psi_fp, qc_fp, "  %6d%12.6f\n", coords.index[cc].at(s)+1, coords.coeff[cc].at(s));
  }
}

void FRAG::print_connectivity(std::string psi_fp, FILE *qc_fp, const int id, const int offset) const {
  oprintf(psi_fp, qc_fp, "\t---Fragment %d Bond Connectivity---\n", id+1);
  int i,j;
  for (i=0; i<natom; ++i) {
    oprintf(psi_fp, qc_fp, "\t %d :", i+1+offset);
      for (j=0; j<natom; ++j)
        if (connectivity[i][j]) oprintf(psi_fp, qc_fp, " %d", j+1+offset);
    oprintf(psi_fp, qc_fp, "\n");
  }
  oprintf(psi_fp, qc_fp, "\n");
}

// computes and print B matrix
void FRAG::print_B(std::string psi_fp, FILE *qc_fp) const {
  double **B = compute_B();
  oprintf(psi_fp, qc_fp, "\t---B matrix---\n");
  oprint_matrix(psi_fp, qc_fp, B, coords.simples.size(), 3*natom);
  oprintf(psi_fp, qc_fp, "\n");
  free_matrix(B);
}

}
