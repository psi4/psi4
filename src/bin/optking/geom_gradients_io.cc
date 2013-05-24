/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file    geom_gradients_io.cc
    \ingroup optking
    \brief   functions to read and write the geometry, the gradient and the Hessian
*/

#include <cstring>
#include <sstream>
#include "molecule.h"

#define EXTERN
#include "globals.h"

#include "io.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <psi4-dec.h>
 #include <libmints/molecule.h>
 #include <libmints/matrix.h>
 #include <libmints/wavefunction.h>
 #include <libparallel/parallel.h>
 #include <libmints/writer_file_prefix.h>
//****AVC****//
#include <libefp_solver/efp_solver.h>
#include <libmints/vector3.h>
//****AVC****//
#elif defined(OPTKING_PACKAGE_QCHEM)
 #include <qchem.h> // typedefs INTEGER
 #include "EFP.h"
 void get_carts(double** Carts, double** A, INTEGER** AtNo, INTEGER* NAtoms, bool noghosts);
#endif

namespace opt {

using namespace std;

bool is_integer(const char *check) {
  for (int i=0; i<strlen(check); ++i) {
    if (!isdigit(check[i]))
      return false;
  }
  return true;
}

// read the number of atoms
int read_natoms(void) {
  int natom=0;

#if defined(OPTKING_PACKAGE_PSI)
  natom = psi::Process::environment.molecule()->natom();

#elif defined(OPTKING_PACKAGE_QCHEM)

  // read number of atoms from QChem
  natom = rem_read(REM_NATOMS);

  // now substract out EFP fragment atoms ?
  if (Opt_params.efp_fragments) {
    int n = ::EFP::GetInstance()->GetNumEFPatoms();
    natom -= n;
    fprintf(outfile, "\nNumber of atoms besides EFP fragments %d\n", natom);
  }

#endif
  return natom;
}


// Read geometry and gradient into a molecule object.
// The fragments must already exist with the number of atoms in each set
// and with allocated memory.
void MOLECULE::read_geom_grad(void) {
  int nfrag = fragments.size();
  int nallatom = g_natom();
//****AVC****//
  int EFPfrag = efp_fragments.size();
//****AVC****//

fprintf(outfile,"in read_geom_grad\n");
fprintf(outfile,"nfrag %d\n", nfrag);
fprintf(outfile,"EFPfrag %d\n", EFPfrag);
fflush(outfile);

#if defined(OPTKING_PACKAGE_PSI)

  using namespace psi;

  SharedMatrix pgradient;
  if (psi::Process::environment.wavefunction()) {
    pgradient = psi::Process::environment.wavefunction()->gradient();
  } else {
    pgradient = psi::Process::environment.gradient();
  }

  Matrix& gradient = *pgradient.get();

gradient.print_out();
fflush(outfile);

/*
//  AVC - NONE OF THIS SHOULD END UP IN THE WORKING CODE -- JUST A HACK TO GET A REASONABLE GRADIENT
  for(int i=0; i<2; i++)
    for(int j=3; j<6; j++)
      gradient.set(i,j,0.00);
  gradient.set(0,5,0.000514);
  gradient.set(1,5,0.000818);

if(nfrag > 0) {
  SharedMatrix efp_grad = gradient.clone();
  SharedMatrix qm_grad(new Matrix(nallatom - EFPfrag*3, 6));
  vector<SharedMatrix> dummy_grad;
  dummy_grad.push_back(efp_grad);
  dummy_grad.push_back(qm_grad);

  SharedMatrix mixed_grad = mixed_grad->vertcat(dummy_grad);

  gradient = *mixed_grad.get();
  gradient.print();
}
//****AVC****//

  boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();
  Matrix geometry = mol->geometry();

fprintf(outfile,"geometry:\n");
geometry.print_out();
fflush(outfile);

  energy = psi::Process::environment.globals["CURRENT ENERGY"];

fprintf(outfile,"energy %15.10lf\n", energy);
fflush(outfile);

  int atom =0;
//****AVC****//
  for (int f=0; f<EFPfrag; ++f) {
    double **geom = efp_fragments[f]->get_xyz_pointer();
    double *force  = efp_fragments[f]->get_forces_pointer();
    double *com   = efp_fragments[f]->get_com_pointer();

fprintf(outfile, "index for fragment %d : %d\n", f, efp_fragments[f]->get_libmints_grad_index());

    for(int i=0; i<3; i++) {
      force[i]   = - gradient(efp_fragments[f]->get_libmints_grad_index(), i);  //forces
      force[i+3] = - gradient(efp_fragments[f]->get_libmints_grad_index(), i+3);//torques

      com[i] = p_efp->get_com(f)[i];  //center of mass

      for(int j=0; j<3; j++)          //geometry
        geom[i][j] = geometry(efp_fragments[f]->get_libmints_geom_index()+i, j);

      atom++;
    }

fprintf(outfile, "\nEFP Fragment %d:\n", f);
fprintf(outfile, "geom\n");
print_matrix(outfile, geom, 3, 3); fprintf(outfile, "\n");
fprintf(outfile, "force\n");
print_array(outfile, force, 6); fprintf(outfile, "\n");
fprintf(outfile, "com\n");
print_array(outfile, com, 3); fprintf(outfile, "\n");
fflush(outfile);

  }
//****AVC****//
  for (int f=0; f<nfrag; ++f) {
      double *Z = fragments[f]->g_Z_pointer();
      double **geom = fragments[f]->g_geom_pointer();
      double **grad = fragments[f]->g_grad_pointer();

      for (int i=0; i<fragments[f]->g_natom(); ++i) {
          Z[i] = mol->Z(fragments[f]->get_libmints_geom_index()+i);

          geom[i][0] = geometry(fragments[f]->get_libmints_geom_index()+i, 0);
          geom[i][1] = geometry(fragments[f]->get_libmints_geom_index()+i, 1);
          geom[i][2] = geometry(fragments[f]->get_libmints_geom_index()+i, 2);

          grad[i][0] = gradient(fragments[f]->get_libmints_grad_index()+i, 0);
          grad[i][1] = gradient(fragments[f]->get_libmints_grad_index()+i, 1);
          grad[i][2] = gradient(fragments[f]->get_libmints_grad_index()+i, 2);

          atom++;
      }
fprintf(outfile, "\nQM Fragment:\n"); fflush(outfile);
fprintf(outfile, "Z\n");
print_array(outfile, Z, fragments[f]->g_natom());
fprintf(outfile, "geom\n");
print_matrix(outfile, geom, fragments[f]->g_natom(), 3); fprintf(outfile, "\n");
fprintf(outfile, "grad\n");
print_matrix(outfile, grad, fragments[f]->g_natom(), 3); fprintf(outfile, "\n");

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
    psi::WorldComm->sync();
    fin.close();
  } // end try reading geometries
  catch (std::ios_base::failure & bf) {
    printf("Error reading molecular geometry and gradient\n");
    fprintf(outfile,"Error reading molecular geometry and gradient\n");
    throw(INTCO_EXCEPT("Error reading molecular geometry and gradient"));
  }
#endif

#elif defined(OPTKING_PACKAGE_QCHEM)

  int EFPfrag = 0;
  int EFPatom = 0;
  if (Opt_params.efp_fragments) {
    EFPfrag = EFP::GetInstance()->NFragments();
    std::vector<int> AtomsinEFP = EFP::GetInstance()->NEFPAtoms();
    for (int i = 0; i < EFPfrag; ++i)
      EFPatom += AtomsinEFP[i];
    fprintf(outfile,"\t %d EFP fragments containing %d atoms.\n", EFPfrag, EFPatom);
  }

  double *QX;
  INTEGER *QZ, QNATOMS;
  bool Qnoghosts = true;

  ::get_carts(NULL, &QX, &QZ, &QNATOMS, Qnoghosts);

  int QNATOMS_real = g_natom();
  if (QNATOMS_real != (QNATOMS-EFPatom))
    QCrash("Number of computed real atoms is inconsistent.");

  fprintf(outfile, "\tNATOMS (total)=%d (minus EFP)=%d\n", QNATOMS, QNATOMS_real);

  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"\tCartesian coordinates from ::get_carts().\n");
    print_array(outfile, QX, 3*QNATOMS);
  }

  int Numgrad = QNATOMS_real*3 + 6*EFPfrag;
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
  for (int f=0; f<EFPfrag; ++f) {

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
  for (int i=0; i<interfragments.size(); ++i)
    interfragments[i]->update_reference_points();

  return;
}



void MOLECULE::symmetrize_geom(void) {

#if defined(OPTKING_PACKAGE_PSI)

  // put matrix into environment molecule and it will symmetrize it
  double **geom_2D = g_geom_2D();
  psi::Process::environment.molecule()->set_geometry(geom_2D);
  psi::Process::environment.molecule()->symmetrize();
  free_matrix(geom_2D);

  psi::Matrix geom = psi::Process::environment.molecule()->geometry();
  geom_2D = geom.pointer(); // don't free; it's shared
  set_geom_array(geom_2D[0]);

#elif defined(OPTKING_PACKAGE_QCHEM)

  // not implemented yet

#endif
}

void MOLECULE::write_geom(void) {

#if defined(OPTKING_PACKAGE_PSI)

  double **geom_2D = g_geom_2D();
fprintf(outfile,"molecule::write_geom outputting to environment:\n");
print_matrix(outfile, geom_2D, efp_fragments.size()*3, 3); fflush(outfile);
  psi::Process::environment.molecule()->set_geometry(geom_2D);
  psi::Process::environment.molecule()->update_geometry();
  free_matrix(geom_2D);
fflush(outfile);

/*
  double *g = g_geom_array();
  p_efp->set_coordinates(1,g); // type is POINTS
  free_array(g);
*/

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
  try {
    std::string hess_fname = psi::get_writer_file_prefix() + ".hess";
    if_Hcart.open(hess_fname.c_str(), ios_base::in);
    int n;
    if_Hcart >> n; // read natom
    if_Hcart >> n; // read natom*6 (?)

    for (int i=0; i<Ncart; ++i)
      for (int j=0; j<Ncart; ++j)
        if_Hcart >> H_cart[i][j];

    psi::WorldComm->sync();
    if_Hcart.close();
  }
  catch (std::ios_base::failure & bf) {
    printf("Error reading cartesian Hessian matrix\n");
    fprintf(outfile,"Error reading cartesian Hessian matrix\n");
    throw(INTCO_EXCEPT("Error reading cartesian Hessian matrix"));
  }

#elif defined(OPTKING_PACKAGE_QCHEM)

  double *Hess = init_array(Ncart*Ncart);

  FileMan_Open_Read(FILE_NUCLEAR_HESSIAN);
  FileMan(FM_READ, FILE_NUCLEAR_HESSIAN, FM_DP, Ncart*Ncart, 0, FM_BEG, Hess);
  FileMan_Close(FILE_NUCLEAR_HESSIAN);
  //print_matrix(outfile, &Hess, 1, Ncart * Ncart);

// instead read from $QCSCRATCH/scr_dir_name/HESS ?

  int cnt = 0;
  for (int i=0; i<Ncart; ++i)
    for (int j=0; j<Ncart; ++j)
      H_cart[i][j] = Hess[cnt++];

  free_array(Hess); // streamline this read later

#endif

  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"\tCartesian Hessian matrix read: \n");
    print_matrix(outfile, H_cart, Ncart, Ncart);
    fflush(outfile);
  }
  else {
      fprintf(outfile, "\tCartesian Hessian matrix read in from external file.\n");
      fflush(outfile);
  }

  return H_cart;
}

}

