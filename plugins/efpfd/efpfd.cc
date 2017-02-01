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

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libefp_solver/efp_solver.h>

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace efpfd {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "EFPFD"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C" 
PsiReturnType efpfd(Options& options)
{
    int print = options.get_int("PRINT");
    outfile->Printf("Starting efpfd()\n");
    options.print();

    // 3-pt or 5-pt formula?
    int pts = options.get_int("POINTS");
    double disp_size = 0.1;

    // Create efp object
    psi::efp::EFP efp_(options);

    // Uses coordinates from molecule object to determing the location of the first
    // three points of each fragment; then calls efp_set_coordinates to initialize
    // the COM and rotation matrix for the fragments.
    efp_.SetGeometry();

    int nfrag = efp_.get_frag_count();
    outfile->Printf("\tNumber of fragments = %d\n", nfrag);

// The code below shows how to access the information for the full (>3) atoms
// in a fragment is the solvent is e.g., methanol.
/*
    // Number of atoms in each fragment
    int *atoms_per_frag = new int[nfrag];
    for (int i=0; i<nfrag; ++i)
      atoms_per_frag[i]= efp_.get_frag_atom_count(i);

    outfile->Printf("\n\tAtoms in fragment:\n");
    for (int i=0; i<nfrag; ++i)
      outfile->Printf("\t %d : %d\n", i+1, atoms_per_frag[i]);

    // Masses of all fragment atoms
    outfile->Printf("\tMasses in fragment:");
    for (int i=0; i<nfrag; ++i) {
      double *tmp_mass = efp_.get_frag_atom_mass(i);
      outfile->Printf("\n\t %d : ", i+1);
      for (int a=0; a<atoms_per_frag[i]; ++a)
        outfile->Printf("%8.5f ", tmp_mass[a]);
      outfile->Printf("\n");
      delete [] tmp_mass;
    }

    // Atomic numbers of all fragment atoms
    outfile->Printf("\n\tAtomic Numbers in Fragment:");
    for (int i=0; i<nfrag; ++i) {
      double *tmp_Zs = efp_.get_frag_atom_Z(i);
      outfile->Printf("\n\t %d :", i+1);
      for (int a=0; a<atoms_per_frag[i]; ++a)
        outfile->Printf(" %3.1f", tmp_Zs[a]);
      delete [] tmp_Zs;
    }
    delete [] atoms_per_frag;

    // How to access COM of full fragments
    outfile->Printf("\n\tCenters of mass of fragment:");
    for (int i=0; i<nfrag; ++i) {
      double *tmp_com = efp_.get_com(i);
      outfile->Printf("\n\t %d :", i+1);
      outfile->Printf("%15.10lf %15.10lf %15.10lf", tmp_com[0], tmp_com[1], tmp_com[2]);
      delete [] tmp_com;
    }
    outfile->Printf("\n");
*/

    // Compute and store analytic gradient
    SharedMatrix grad_ = psi::Process::environment.gradient();

    // Get reference geometry
    boost::shared_ptr<Molecule> mol = Process::environment.molecule();
    Matrix ref_geom_temp = mol->geometry();
    SharedMatrix ref_geom(ref_geom_temp.clone());
    ref_geom->set_name("Reference Geometry");

    // Get rotation axes ready.
    std::vector<Vector3> xyz_axes;
    xyz_axes.push_back( Vector3(1,0,0));
    xyz_axes.push_back( Vector3(0,1,0));
    xyz_axes.push_back( Vector3(0,0,1));
    char xyz_char[] = "xyz";

    // Build list of all the displacement coefficients
    std::vector<int> disp_coeff;
    if (pts == 3) {
      disp_coeff.push_back(-1);
      disp_coeff.push_back(+1);
    }
    else if (pts == 5) {
      disp_coeff.push_back(-2);
      disp_coeff.push_back(-1);
      disp_coeff.push_back(+1);
      disp_coeff.push_back(+2);
    }

    // Construct list of displaced geometries
    std::vector< SharedMatrix > disp_geoms;
    for(int f=0; f<nfrag; ++f) {

      // get first atom of this fragment in molecule object
      int first_atom = mol->fragment_atom_pair(f).first;

      // Build geometry of current fragment (assume we only need 3 atoms)
      Matrix frag_ref_geom_mat(3,3);
      for (int a=0; a<3; ++a)
        for (int xyz=0; xyz<3; ++xyz)
          frag_ref_geom_mat.set(a, xyz, ref_geom->get(first_atom+a, xyz));

      SharedMatrix frag_ref_geom(frag_ref_geom_mat.clone());
      frag_ref_geom->set_name("Reference geometry of frag " + to_string(f+1));
    frag_ref_geom->print();

      // First make list of displaced geometries for this fragment
      std::vector< SharedMatrix > frag_disp_geoms;

      // displace along COM coordinates
      for (int xyz=0; xyz<3; ++xyz) {

        for (int disp=0; disp<disp_coeff.size(); ++disp) { // (-1,+1) or (-2,-1,+1,+2)
          SharedMatrix geom1(frag_ref_geom->clone());
          geom1->set_name("Disp of frag " + to_string(f+1) + ": " +
            to_string(disp_coeff[disp]) + " COM " + xyz_char[xyz]);

          for (int a=0; a<3; ++a)
            geom1->add(a, xyz, disp_coeff[disp]*disp_size);

          frag_disp_geoms.push_back(geom1);
        }

      }

      // Find COM of this TOTAL efp fragment and move fragment to put its COM at origin
      double *frag_com = efp_.get_com(f);

      for (int a=0; a<3; ++a)
        for (int xyz=0; xyz<3; ++xyz)
          frag_ref_geom->add(a, xyz, -1.0*frag_com[xyz]);
    frag_ref_geom->set_name("Reference geometry moved to COM");
    frag_ref_geom->print();

      // displace along Rx, Ry, and Rz
      for (int xyz=0; xyz<3; ++xyz) {

        for (int disp=0; disp<disp_coeff.size(); ++disp) { // (-1,+1) or (-2,-1,+1,+2)
          SharedMatrix geom_rot_mat(frag_ref_geom->clone());
          SharedMatrix geom_rot = geom_rot_mat->matrix_3d_rotation(xyz_axes[xyz],
            disp_coeff[disp]*disp_size, 0);
          geom_rot->set_name("Disp of frag " + to_string(f+1) + ": " +
            to_string(disp_coeff[disp]) + " R" + xyz_char[xyz]);

          // Move geometry back to original COM
          for (int a=0; a<3; ++a)
            for (int xyz2=0; xyz2<3; ++xyz2)
              geom_rot->add(a, xyz2, frag_com[xyz2]);

          frag_disp_geoms.push_back(geom_rot);
        }

      } // end this Rx, Ry or Rz

      delete [] frag_com;

      // Insert geometries for fragment into full geometries and add to list
      for (int f_disp=0; f_disp<frag_disp_geoms.size(); ++f_disp) {

        SharedMatrix full_geom(ref_geom->clone());
        full_geom->set_name( frag_disp_geoms[f_disp]->name() );

        for (int a=0; a<3; ++a)
          for (int xyz=0; xyz<3; ++xyz)
            full_geom->set(first_atom+a, xyz, frag_disp_geoms[f_disp]->get(a, xyz));

        disp_geoms.push_back(full_geom);
      }

      frag_disp_geoms.clear();
    } // end frag

    disp_geoms.push_back(ref_geom); // add reference geometry to end

    for(int i=0; i<disp_geoms.size(); ++i)
      disp_geoms[i]->print();

    int Ndisp = disp_geoms.size();
    // check if expected Ndisp == disp_geoms.size()

    // Compute all energies.
    double *E = new double [Ndisp];
    int cnt = 0;
    for (int i=0; i<Ndisp; ++i) {
      Matrix mat_geom;
      mat_geom.copy(disp_geoms[i]);
      mol->set_geometry(mat_geom);
      efp_.SetGeometry(); // updates EFP coordinates
      efp_.Compute();
      E[cnt++] = psi::Process::environment.globals["CURRENT ENERGY"];
    }

    outfile->Printf("Displaced energies\n");
    for (int i=0; i<Ndisp; ++i)
      outfile->Printf("%15.10lf\n", E[i]);

    int Ncoord;
    if (pts == 3)
      Ncoord = disp_geoms.size()/2;
    else if (pts == 5)
      Ncoord = disp_geoms.size()/4;

    double *fd_grad = new double [Ncoord];

    if (pts == 3) {
      for (int i=0; i<Ncoord; ++i)
        fd_grad[i] = (E[2*i+1] - E[2*i]) / (2.0 * disp_size);
    }
    else if (pts == 5) {
      for (int i=0; i<Ncoord; ++i) {
        int I = 4*i;
        fd_grad[i] = (E[I] - 8.0*E[I+1] + 8.0*E[I+2] - E[I+3]) / (12.0 * disp_size);
      }
    }

    outfile->Printf("Finite-difference gradient\n");
    for (int f=0; f<nfrag; ++f) {
      outfile->Printf("\t %d: ", f+1);
      for (int i=0; i<6; ++i)
        outfile->Printf("%15.10lf", fd_grad[6*f+i]);
      outfile->Printf("\n");
    }

    outfile->Printf("Analytic gradient");
    grad_->print();

    return Success;
}

}} // End namespaces

