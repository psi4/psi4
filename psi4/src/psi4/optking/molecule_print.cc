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

/*! \file molecule_print.cc
    \ingroup optking
    \brief molecule class functions for printing and string descriptions
*/

#include "molecule.h"

#include <iostream>
#include <sstream>

#include "linear_algebra.h"
#include "atom_data.h"
#include "psi4/optking/physconst.h"
#include "psi4/libparallel/ParallelPrinter.h"
#include "print.h"
#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
 #include "psi4/libmints/molecule.h"
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
 #include "EFP.h"
#endif

namespace opt {

using namespace std;

void MOLECULE::print_geom_out(void) {
#if defined(OPTKING_PACKAGE_QCHEM)
  oprintf_out("\tCartesian Geometry (au)\n");
#elif defined(OPTKING_PACKAGE_PSI)
  oprintf_out("\tCartesian Geometry (in Angstrom)\n");
#endif

  for (std::size_t i=0; i<fragments.size(); ++i)
    fragments[i]->print_geom(psi_outfile, qc_outfile);
}

void MOLECULE::print_geom_out_irc(void) {
#if defined(OPTKING_PACKAGE_QCHEM)
  oprintf_out("@IRC    Cartesian Geometry (au)\n");
#elif defined(OPTKING_PACKAGE_PSI)
  oprintf_out("@IRC    Cartesian Geometry (in Angstrom)\n");
#endif

  for (std::size_t i=0; i<fragments.size(); ++i)
    fragments[i]->print_geom_irc(psi_outfile, qc_outfile);
}

// This function is only used for an optional trajectory file.
// The awkward itershift is to decrement in the initial geometry to "iteration 0"
void MOLECULE::print_xyz(int iter_shift) {
  FILE *qc_fp = NULL;

#if defined(OPTKING_PACKAGE_QCHEM)
  qc_fp = fopen("geoms.xyz", "a");
#endif
  oprintf("geoms.xyz", qc_fp, "%d\n", g_natom());
  oprintf("geoms.xyz", qc_fp, "Geometry for iteration %d\n", p_Opt_data->g_iteration()+iter_shift);
  for (std::size_t i=0; i<fragments.size(); ++i)
    fragments[i]->print_geom("geoms.xyz", qc_fp);
#if defined(OPTKING_PACKAGE_QCHEM)
  fclose(qc_fp);
#endif
  return;
}

void MOLECULE::print_xyz_irc(int point, bool forward) {
  FILE *qc_fp = NULL;

  if(forward) {
#if defined(OPTKING_PACKAGE_QCHEM)
    qc_fp = fopen("irc_forward.xyz", "a");
#endif
    oprintf("irc_forward.xyz", qc_fp, "%d\n", g_natom());
    oprintf("irc_forward.xyz", qc_fp, "IRC point %d\n", point);
    for (std::size_t i=0; i<fragments.size(); ++i)
      fragments[i]->print_geom("irc_forward.xyz", qc_fp);
#if defined(OPTKING_PACKAGE_QCHEM)
    fclose(qc_fp);
#endif
  }
  else {
#if defined(OPTKING_PACKAGE_QCHEM)
    qc_fp = fopen("irc_backward.xyz", "a");
#endif
    oprintf("irc_backward.xyz", qc_fp, "%d\n", g_natom());
    oprintf("irc_backward.xyz", qc_fp, "IRC point %d\n", point);
    for (std::size_t i=0; i<fragments.size(); ++i)
      fragments[i]->print_geom("irc_backward.xyz", qc_fp);
#if defined(OPTKING_PACKAGE_QCHEM)
    fclose(qc_fp);
#endif
  }
  return;
}

// print internal coordinates to text file
void MOLECULE::print_intco_dat(std::string psi_fp, FILE *qc_fp) {
  for (std::size_t i=0; i<fragments.size(); ++i) {
    int first = g_atom_offset(i);
    if (fragments[i]->is_frozen())
        oprintf(psi_fp, qc_fp, "F* %d %d\n", first+1, first + fragments[i]->g_natom());
    else
        oprintf(psi_fp, qc_fp, "F %d %d\n", first+1, first + fragments[i]->g_natom());
    fragments[i]->print_intco_dat(psi_fp, qc_fp, g_atom_offset(i));
  }

  for (std::size_t I=0; I<interfragments.size(); ++I) {
    int frag_a = interfragments[I]->g_A_index();
    int frag_b = interfragments[I]->g_B_index();
    oprintf(psi_fp, qc_fp, "I %d %d\n", frag_a+1, frag_b+1);

    for (int i=0; i<6; ++i)
      oprintf(psi_fp, qc_fp, " %d", (int) interfragments[I]->coordinate_on(i));
    oprintf(psi_fp, qc_fp, "\n");

    interfragments[I]->print_intco_dat(psi_fp, qc_fp, g_atom_offset(frag_a), g_atom_offset(frag_b));
  }
}

// Fetches the string definition of an internal coordinate from global index
std::string MOLECULE::get_coord_definition_from_global_index(int index) const{
  int Nintra = Ncoord_intrafragment();
  int Ninter = Ncoord_interfragment();
  int Nfb = 0;
#if defined (OPTKING_PACKAGE_QCHEM)
  Nfb = Ncoord_fb_fragment();
#endif

  if ( index < 0 || index >= (Nintra + Ninter + Nfb) ) {
    oprintf_out( "get_coord_definition(): index %d out of range", index);
    throw(INTCO_EXCEPT("get_coord_definition(): index out of range"));
  }

  // coordinate is an intrafragment coordinate
  if (index < Nintra) {
    std::size_t f;

    // go to last fragment or first that isn't past the desired index
    for (f=0; f<fragments.size(); ++f)
      if (index < g_coord_offset(f))
        break;
    --f;

    return fragments[f]->get_coord_definition(index - g_coord_offset(f), g_atom_offset(f));
  }

  // coordinate is an interfragment coordinate
  if (index < Nintra + Ninter) {
    std::size_t f;

    for (f=0; f<interfragments.size(); ++f)
      if (index < g_interfragment_coord_offset(f))
        break;
    --f;

    return interfragments[f]->get_coord_definition(index - g_interfragment_coord_offset(f));
  }

#if defined (OPTKING_PACKAGE_QCHEM)
  if (index < Nintra + Ninter + Nfb) {
    for (std::size_t f=0; f<fb_fragments.size(); ++f)
      if (index < g_fb_fragment_coord_offset(f))
        break;
    --f;

    return fb_fragments[f]->get_coord_definition(index - Nintra - Nextra);
  }
#endif
  oprintf_out("Warning: impossible case in MOLECULE::get_coord_definition_from_global_index\n");
  return 0;
}

void MOLECULE::print_coords(std::string psi_fp, FILE *qc_fp) const {
  for (std::size_t i=0; i<fragments.size(); ++i) {
    oprintf(psi_fp, qc_fp, "\t---Fragment %d Intrafragment Coordinates---\n", i+1);
    offlush_out();
    fragments[i]->print_simples(psi_fp, qc_fp, g_atom_offset(i));

    if (Opt_params.coordinates == OPT_PARAMS::DELOCALIZED) {
      oprintf_out("\tThere are %d delocalized coordinates formed from these simples.\n\n", Ncoord());
      if ( (Opt_params.print_lvl > 1 && p_Opt_data->g_iteration()==1) ||
            Opt_params.print_lvl > 3 )
        fragments[i]->print_combinations(psi_fp, qc_fp);
    }
    else if ( Opt_params.coordinates == OPT_PARAMS::NATURAL)  {
      oprintf_out("\tThere are %d natural coordinates formed from these simples.\n");
    }
  }
  for (std::size_t i=0; i<interfragments.size(); ++i) {
    int a = interfragments[i]->g_A_index();
    int b = interfragments[i]->g_B_index();
    interfragments[i]->print_coords(psi_fp, qc_fp, g_atom_offset(a), g_atom_offset(b));
  }
  for (std::size_t i=0; i<fb_fragments.size(); ++i) {
    oprintf(psi_fp, qc_fp,"\t---Fragment %d FB fragment Coordinates---\n", i+1);
    fb_fragments[i]->print_coords(psi_fp, qc_fp);
  }
}

void MOLECULE::print_simples(std::string psi_fp, FILE *qc_fp) const {
  for (std::size_t i=0; i<fragments.size(); ++i) {
    oprintf(psi_fp, qc_fp, "\t---Fragment %d Intrafragment Coordinates---\n", i+1);
    fragments[i]->print_simples(psi_fp, qc_fp, g_atom_offset(i));
  }
  for (std::size_t i=0; i<interfragments.size(); ++i) {
    int a = interfragments[i]->g_A_index();
    int b = interfragments[i]->g_B_index();
    interfragments[i]->print_coords(psi_fp, qc_fp, g_atom_offset(a), g_atom_offset(b));
    // print just the simples here intead ?
    // interfragments[i]->inter_frag->print_simples(psi_fp, qc_fp, g_atom_offset(a), g_atom_offset(b));
  }

  for (std::size_t i=0; i<fb_fragments.size(); ++i) {
    oprintf(psi_fp, qc_fp,"\t---Fragment %d FB fragment Coordinates---\n", i+1);
    fb_fragments[i]->print_simples(psi_fp, qc_fp);
  }
}

}
