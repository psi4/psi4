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

/*! \file    geom_gradients_io.cc
    \ingroup optking
    \brief   functions to read and write the geometry, the gradient and the Hessian
*/

#include <cstring>
#include <sstream>
#include "molecule.h"
#include "psi4/psi4-dec.h"
#include "print.h"
#define EXTERN
#include "globals.h"

#include "io.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include "psi4/psi4-dec.h"
 #include "psi4/libmints/molecule.h"
 #include "psi4/libmints/matrix.h"
 #include "psi4/libmints/wavefunction.h"
 #include "psi4/libparallel/parallel.h"
 #include "psi4/libmints/writer_file_prefix.h"
#elif defined(OPTKING_PACKAGE_QCHEM)
 #include <qchem.h> // typedefs INTEGER
 #include "EFP.h"
 void get_carts(double** Carts, double** A, INTEGER** AtNo, INTEGER* NAtoms, bool noghosts);
#endif

namespace opt {

class BROKEN_SYMMETRY_EXCEPT;

using namespace std;

bool is_integer(const char *check) {
  for (ULI i=0; i<strlen(check); ++i) {
    if (!isdigit(check[i]))
      return false;
  }
  return true;
}

// read the number of atoms
int read_natoms(void) {
  int natom=0;

#if defined(OPTKING_PACKAGE_PSI)
  natom = psi::Process::environment.legacy_molecule()->natom();

#elif defined(OPTKING_PACKAGE_QCHEM)

  // read number of atoms from QChem
  natom = rem_read(REM_NATOMS);

  // now substract out EFP fragment atoms ?
  if (Opt_params.fb_fragments) {
    int n = ::EFP::GetInstance()->GetNumEFPatoms();
    natom -= n;
    oprintf_out( "\nNumber of atoms besides EFP fragments %d\n", natom);
  }

#endif
  return natom;
}


// Read geometry and gradient into a molecule object.
// The fragments must already exist with the number of atoms in each set
// and with allocated memory.
void MOLECULE::read_geom_grad(void) {
  int nfrag = fragments.size();
  //int nallatom = g_natom();

#if defined(OPTKING_PACKAGE_PSI)

  using namespace psi;

  SharedMatrix pgradient;
  pgradient = psi::Process::environment.gradient();

  Matrix& gradient = *pgradient.get();

  std::shared_ptr<Molecule> mol = psi::Process::environment.legacy_molecule();
  Matrix geometry = mol->geometry();

  energy = psi::Process::environment.globals["CURRENT ENERGY"];

  int atom =0;
  for (int f=0; f<nfrag; ++f) {
      double *Z = fragments[f]->g_Z_pointer();
      double **geom = fragments[f]->g_geom_pointer();
      double **grad = fragments[f]->g_grad_pointer();

      for (int i=0; i<fragments[f]->g_natom(); ++i) {
          Z[i] = mol->Z(atom);

          geom[i][0] = geometry(atom, 0);
          geom[i][1] = geometry(atom, 1);
          geom[i][2] = geometry(atom, 2);

          grad[i][0] = gradient(atom, 0);
          grad[i][1] = gradient(atom, 1);
          grad[i][2] = gradient(atom, 2);

          atom++;
      }
  }


#if 0
  // for now read from file 11
  std::ifstream fin;

  string lbl;
  char cline[256];
  int junk;

  try {
    // CDS: File11 isn't named this way anymore
    fin.open("psi.file11.dat", ios_base::in);
    fin.clear();
    fin.exceptions(ios_base::failbit);
    fin.seekg(0);

    // rewind, determine number of geometries
    int nlines = 0;
    try { // ignore exceptions only here
      while (fin.peek() != EOF) {
        fin.getline(cline, 256);
        ++nlines;
      }
    }
    catch (std::ios_base::failure & bf) { }

    int ngeom = nlines / (2*nallatom+2);

    // rewind through old geometries
    fin.clear();
    fin.seekg(0);
    for (int i=0; i< (ngeom-1)*(2*nallatom+2); ++i)
      fin.getline(cline, 256);

    // header line
    fin.getline(cline, 256);

    fin >> junk;
    if (junk != nallatom) throw(INTCO_EXCEPT("wrong number of atoms in text file"));
    fin >> energy;

    for (int f=0; f<nfrag; ++f) {
      double *Z = fragments[f]->g_Z_pointer();
      double **geom = fragments[f]->g_geom_pointer();

      for (int i=0; i<fragments[f]->g_natom(); ++i) {
        //// check to see if atomic number or symbol
        //if ( isalpha((lbl.c_str())[0]) )
          //Z[i] = Element_to_Z(lbl);
        //else
          fin >> Z[i];

        fin >> geom[i][0];
        fin >> geom[i][1];
        fin >> geom[i][2];
      }
    }

    for (int f=0; f<nfrag; ++f) {
      double **grad = fragments[f]->g_grad_pointer();
      for (int i=0; i<fragments[f]->g_natom(); ++i) {
        fin >> grad[i][0];
        fin >> grad[i][1];
        fin >> grad[i][2];
      }
    }
    fin.close();
  } // end try reading geometries
  catch (std::ios_base::failure & bf) {
    printf("Error reading molecular geometry and gradient\n");
    oprintf_out("Error reading molecular geometry and gradient\n");
    throw(INTCO_EXCEPT("Error reading molecular geometry and gradient"));
  }
#endif

#elif defined(OPTKING_PACKAGE_QCHEM)

  int FBfrag = 0;
  int FBatom = 0;
  if (Opt_params.fb_fragments) {
    FBfrag = EFP::GetInstance()->NFragments();
    std::vector<int> AtomsinFB = EFP::GetInstance()->NEFPAtoms();
    for (int i = 0; i < FBfrag; ++i)
      FBatom += AtomsinFB[i];
    oprintf_out("\t %d FB fragments containing %d atoms.\n", FBfrag, FBatom);
  }

  double *QX;
  INTEGER *QZ, QNATOMS;
  bool Qnoghosts = true;

  ::get_carts(NULL, &QX, &QZ, &QNATOMS, Qnoghosts);

  int QNATOMS_real = g_natom();
  if (QNATOMS_real != (QNATOMS-FBatom))
    QCrash("Number of computed real atoms is inconsistent.");

  oprintf_out( "\tNATOMS (total)=%d (minus EFP)=%d\n", QNATOMS, QNATOMS_real);

  if (Opt_params.print_lvl >= 3) {
    oprintf_out("\tCartesian coordinates from ::get_carts().\n");
    oprint_array_out(QX, 3*QNATOMS);
  }

  int Numgrad = QNATOMS_real*3 + 6*FBfrag;
  double* QGrad = init_array(Numgrad);

  ::FileMan_Open_Read(FILE_NUCLEAR_GRADIENT);
  ::FileMan(FM_READ, FILE_NUCLEAR_GRADIENT, FM_DP, Numgrad, 0, FM_BEG, QGrad);
  ::FileMan_Close(FILE_NUCLEAR_GRADIENT);

  // store geometry and gradient of real atoms in fragment objects
  int cnt=0;
  for (int f=0; f<nfrag; ++f) {

    double *Z = fragments[f]->g_Z_pointer();
    double **geom = fragments[f]->g_geom_pointer();
    double **grad = fragments[f]->g_grad_pointer();

    for (int i=0; i<fragments[f]->g_natom(); ++i) {
      Z[i] = QZ[cnt];
      for (int xyz=0; xyz<3; ++xyz) {
        geom[i][xyz] = QX[3*cnt+xyz];
        grad[i][xyz] = QGrad[3*cnt+xyz];
      }
      ++cnt;
    }

  }
  // Now read in gradients for EFP fragments
  for (int f=0; f<FBfrag; ++f) {

    double *efp_f = init_array(6);
    for (int i=0; i<6; ++i)
      efp_f[i] = -1 * QGrad[3*QNATOMS_real + 6*f + i];

    efp_fragments[f]->set_forces(efp_f);
    free_array(efp_f);
  }
  free_array(QGrad);

  if (Opt_params.efp_fragments) {
    energy = EFP::GetInstance()->GetEnergy();
  }
  else {
    double E_tmp;
    ::FileMan_Open_Read(FILE_ENERGY);
    ::FileMan(FM_READ, FILE_ENERGY, FM_DP, 1, FILE_POS_CRNT_TOTAL_ENERGY, FM_BEG, &E_tmp);
    ::FileMan_Close(FILE_ENERGY);
    energy = E_tmp;
  }

#endif

  // update interfragment reference points if they exist
  for (ULI i=0; i<interfragments.size(); ++i)
    interfragments[i]->update_reference_points();

  return;
}

void MOLECULE::symmetrize_geom(bool flexible) {

#if defined(OPTKING_PACKAGE_PSI)
  double symm_tol = Opt_params.symm_tol;

  // put matrix into environment.legacy_molecule and it will symmetrize it
  double **geom_2D = g_geom_2D();
  int max_iter = (flexible ? 10 : 1);
  bool symmetrized = false;
  int iter = 0;

  while (iter < max_iter && !symmetrized) {
    try {
      ++iter;
      psi::Process::environment.legacy_molecule()->set_geometry(geom_2D);
      psi::Process::environment.legacy_molecule()->symmetrize(symm_tol, true);
      oprintf_out("\tSuccessfully symmetrized geometry.\n");
      symmetrized = true;
      free_matrix(geom_2D);
    }
    catch (psi::PsiException exc) {
      oprintf_out("\tUnable to symmetrize geometry.\n");
      if (flexible && iter < 10) {
        symm_tol *= 1.5;
        oprintf_out("\tIncreasing symmetry tolerance to %10.4f\n", symm_tol);
      }
      else {
        free_matrix(geom_2D);
        throw(BROKEN_SYMMETRY_EXCEPT("Broken symmetry in OPT::MOLECULE::SYMMETRIZE_GEOM\n"));
      }
    }
  }

  psi::Matrix geom = psi::Process::environment.legacy_molecule()->geometry();
  geom_2D = geom.pointer(); // don't free; it's shared
  set_geom_array(geom_2D[0]);

#elif defined(OPTKING_PACKAGE_QCHEM)

  // not implemented yet

#endif
}

void MOLECULE::write_geom(void) {

#if defined(OPTKING_PACKAGE_PSI)

  double **geom_2D = g_geom_2D();
  psi::Process::environment.legacy_molecule()->set_geometry(geom_2D);
  psi::Process::environment.legacy_molecule()->update_geometry();
  free_matrix(geom_2D);

#elif defined(OPTKING_PACKAGE_QCHEM)

  // rotation matrix is already updated by displace()
  double *geom = g_geom_array();
  if (!Opt_params.efp_fragments_only) {
    int NoReor = rem_read(REM_NO_REORIENT);  // save state
    rem_write(1, REM_NO_REORIENT); // write cartesians ; do NOT reorient
    ::set_carts(geom, true);
    rem_write(NoReor, REM_NO_REORIENT); // replace to original state
  }

#endif
}

double ** OPT_DATA::read_cartesian_H(void) const {

  double **H_cart = init_matrix(Ncart, Ncart);

#if defined(OPTKING_PACKAGE_PSI)
  std::ifstream if_Hcart;
  // Need to enable exceptions in ifstream.
  if_Hcart.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    std::string hess_fname = psi::get_writer_file_prefix(psi::Process::environment.legacy_molecule()->name()) + ".hess";
    if_Hcart.open(hess_fname.c_str(), ios_base::in);
    int n;
    if_Hcart >> n; // read natom
    if_Hcart >> n; // read natom*6 (?)

    for (int i=0; i<Ncart; ++i)
      for (int j=0; j<Ncart; ++j)
        if_Hcart >> H_cart[i][j];

    if_Hcart.close();
  }
  catch (std::ios_base::failure & bf) {
    printf("Error reading cartesian Hessian matrix\n");
    oprintf_out("Error reading cartesian Hessian matrix\n");
    throw(INTCO_EXCEPT("Error reading cartesian Hessian matrix"));
  }

#elif defined(OPTKING_PACKAGE_QCHEM)

  double *Hess = init_array(Ncart*Ncart);

  FileMan_Open_Read(FILE_NUCLEAR_HESSIAN);
  FileMan(FM_READ, FILE_NUCLEAR_HESSIAN, FM_DP, Ncart*Ncart, 0, FM_BEG, Hess);
  FileMan_Close(FILE_NUCLEAR_HESSIAN);
  //oprint_matrix_out(&Hess, 1, Ncart * Ncart);

// instead read from $QCSCRATCH/scr_dir_name/HESS ?

  int cnt = 0;
  for (int i=0; i<Ncart; ++i)
    for (int j=0; j<Ncart; ++j)
      H_cart[i][j] = Hess[cnt++];

  free_array(Hess); // streamline this read later

#endif

  oprintf_out("\tCartesian Hessian matrix read in from external file: \n");
  oprint_matrix_out(H_cart, Ncart, Ncart);

  return H_cart;
}

}
