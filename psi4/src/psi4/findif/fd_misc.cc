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

bool ascending(const std::shared_ptr<VIBRATION> vib1, const std::shared_ptr<VIBRATION> vib2) {
  if (vib1->km < vib2->km)
    return true;
  else
    return false;
}

// function to print out (frequencies and normal modes) vector of vibrations
void print_vibrations(std::shared_ptr<Molecule> mol, std::vector<std::shared_ptr<VIBRATION>> modes) {

    outfile->Printf("Printing vibrations\n");
  char **irrep_lbls = mol->irrep_labels();
  int Natom = mol->natom();

  // compute harmonic frequencies, +/- in wavenumbers
  /* Convert evals from H/(kg bohr^2) to J/(kg m^2) = 1/s^2 */
  /* v = 1/(2 pi c) sqrt( eval ) */
  const double k_convert = pc_hartree2J/(pc_bohr2m * pc_bohr2m * pc_amu2kg);
  const double cm_convert = 1.0/(2.0 * pc_pi * pc_c * 100.0);

  for (int i=0; i<modes.size(); ++i) {
    outfile->Printf("modes %d is %f\n" , i , modes[i]->km);
      if(modes[i]->km < 0.0){
      outfile->Printf("<print if statement >modes %d is %f\n" , i , modes[i]->km);
      modes[i]->cm = -1*cm_convert * sqrt(-k_convert * modes[i]->km);
      }
              else {
      outfile->Printf("<print if statement >modes %d is %f\n" , i , modes[i]->km);
      modes[i]->cm =    cm_convert * sqrt( k_convert * modes[i]->km);
              }
  }

  for (int i=0; i<modes.size(); ++i)
    outfile->Printf("modes %d is %f\n" , i , modes[i]->km);
  
    
  // Sort modes by increasing eigenvalues.
  sort(modes.begin(), modes.end(), ascending);

  // Print out frequencies and irreps to output file.
  outfile->Printf( "\n\t  Irrep      Harmonic Frequency   \n");
  outfile->Printf(   "\t                  (cm-1)          \n");
  outfile->Printf(   "\t-----------------------------------------------\n");

  for(int i=0; i<modes.size(); ++i) {
    if(modes[i]->cm < 0.0)
      outfile->Printf( "\t  %5s   %15.4fi \n", irrep_lbls[modes[i]->irrep], -modes[i]->cm);
    else
      outfile->Printf( "\t  %5s   %15.4f  \n", irrep_lbls[modes[i]->irrep], modes[i]->cm);
  }

  outfile->Printf(   "\t-----------------------------------------------\n");


  // Return list of frequencies to wavefunction object.
  std::shared_ptr<Vector> freq_vector(new Vector(modes.size()));
  for (int i=0; i<modes.size(); ++i)
    freq_vector->set(i, modes[i]->cm);

  // Reture list of normal modes to wavefunction object.
  std::shared_ptr<Vector> nm_vector(new Vector(3*Natom*modes.size()));
  int count = 0;
  for (int i=0; i<modes.size(); ++i) {
    freq_vector->set(i, modes[i]->cm);
    for (int a=0; a<Natom; ++a) {
        for (int xyz=0; xyz<3; xyz++) {
            nm_vector->set(count, modes[i]->lx->get(3*a+xyz));
            count++;
        }
    }
  }

  Process::environment.set_frequencies(freq_vector);

  double sum = 0.0;
  for (int a=0; a<Natom; ++a)
     sum += mol->mass(a);

  // print out normal modes in format that WebMO likes
  outfile->Printf( "\n\tNormal Modes (non-mass-weighted).\n");
  outfile->Printf( "\tMolecular mass is %10.5f amu.\n", sum);
  outfile->Printf( "\tFrequencies in cm^-1; force constants in au.\n");

  for(int i=0; i<modes.size(); ++i) { // print descending order
    if (std::fabs(cm_convert * sqrt(k_convert * std::fabs(modes[i]->km))) < 5.0) continue;
    outfile->Printf("\n");
    if (modes[i]->km < 0.0)
      outfile->Printf( "   Frequency:      %8.2fi\n", cm_convert * sqrt(-k_convert * modes[i]->km));
    else
      outfile->Printf( "   Frequency:      %8.2f\n", cm_convert * sqrt(k_convert * modes[i]->km));

    outfile->Printf(   "   Force constant: %8.4f\n", modes[i]->km);

    //outfile->Printf(   "   IR Intensity: %8.2f\n", irint[i]*ir_prefactor);

    outfile->Printf( "\t     X       Y       Z           mass\t\n");
    for (int a=0; a<Natom; a++) {
      outfile->Printf( "  %s \t", mol->symbol(a).c_str() );

      for (int xyz=0; xyz<3; ++xyz)
        outfile->Printf( "%8.3f", modes[i]->lx->get(3*a+xyz));

      outfile->Printf("%15.6f", mol->mass(a));

      outfile->Printf( "\n");
    }
  }

  // awkward, but need nirrep to free labels
  int Nirrep = mol->point_group()->char_table().nirrep();

  for (int i=0; i<Nirrep; ++i)
    free(irrep_lbls[i]);
  free(irrep_lbls);
}

void save_normal_modes(std::shared_ptr<Molecule> mol,std::shared_ptr<Wavefunction> wfn, std::vector<std::shared_ptr<VIBRATION>> modes)
{
    wfn->set_normalmodes(modes);
    std::string normal_modes_fname = get_writer_file_prefix(mol->name()) + ".molden_normal_modes";
    std::shared_ptr <PsiOutStream> printer(new PsiOutStream(normal_modes_fname, std::ostream::trunc));

    printer->Printf("[Molden Format]\n[FREQ]\n");
    for(int i = 0; i < modes.size(); ++i) { // print descending order
        double frequency_cm = modes[i]->get_cm();
        printer->Printf("%.2f\n", frequency_cm);
    }
    printer->Printf("\n[FR-COORD]\n");
    int Natom = mol->natom();
    for (int a = 0; a < Natom; a++) {
        double x = mol->x(a);
        double y = mol->y(a);
        double z = mol->z(a);
        printer->Printf("%-3s %.6f %.6f %.6f\n",mol->symbol(a).c_str(),x,y,z);
    }
    printer->Printf("\n[FR-NORM-COORD]\n");
    for(int i = 0; i < modes.size(); ++i) {
        printer->Printf("vibration %d\n",i + 1);
        int Natom = mol->natom();
        double norm2 = 0.0;
        for (int a = 0; a < Natom; a++) {
            for (int xyz = 0; xyz < 3; ++xyz){
                norm2 += std::pow(modes[i]->lx->get(3 * a + xyz),2.0);
            }
        }
        double scaling_factor = 1.0 / std::sqrt(norm2);
        for (int a = 0; a < Natom; a++) {
            for (int xyz = 0; xyz < 3; ++xyz){
                double scaled_mode = scaling_factor * modes[i]->lx->get(3 * a + xyz);
                printer->Printf(" %.6f",scaled_mode);
            }
            printer->Printf("\n");
        }
    }
    printer->Printf("\n[INT]\n");
    for(int i = 0; i < modes.size(); ++i) {
        printer->Printf("1.0\n");
    }
}


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
