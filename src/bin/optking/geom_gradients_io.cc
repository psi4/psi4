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

// Reads in the number of QM atoms.
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
  }
#endif
  fprintf(outfile, "\nNumber of atoms besides EFP fragments %d\n", natom);
  return natom;
}


// Read geometry and gradient into a molecule object.
// The fragments must already exist with the number of atoms in each set
// and with allocated memory.
void MOLECULE::read_geom_grad(void) {
  int nfrag = fragments.size();
  int EFPfrag = efp_fragments.size();
  int nallatom = g_natom() + 3 * EFPfrag;

#if defined(OPTKING_PACKAGE_PSI)

  using namespace psi;

  SharedMatrix pgradient;

  // QM
  if (nfrag > 0) {
    if (psi::Process::environment.wavefunction())
      pgradient = psi::Process::environment.wavefunction()->gradient();
    else
      pgradient = psi::Process::environment.gradient();
  }

  // EFP
  if (EFPfrag > 0) {
    if (psi::Process::environment.get_efp()->get_frag_count() > 0)
      pgradient = psi::Process::environment.get_efp()->torque();
    else
      pgradient = psi::Process::environment.efp_torque();
  }

  SharedMatrix grad = pgradient->clone();
  grad->print_out();
fflush(outfile);

  boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();
  Matrix geometry = mol->geometry();
fprintf(outfile,"QM geometry:\n");
geometry.print_out();
fflush(outfile);

  energy = psi::Process::environment.globals["CURRENT ENERGY"];

fprintf(outfile,"energy %15.10lf\n", energy);
fflush(outfile);

  // Put center of mass in each efp fragment
  for (int f=0; f<EFPfrag; ++f) {
    double *f_com  = efp_fragments[f]->get_com_pointer();
    double *tmp_com = p_efp->get_com(f);
    for (int xyz=0; xyz<3; ++xyz)
      f_com[xyz] = tmp_com[xyz];
    delete [] tmp_com;
  }

  // Put forces into each efp fragment
  for (int f=0; f<EFPfrag; ++f) {
    double *f_force = efp_fragments[f]->get_forces_pointer();

    for(int i=0; i<6; i++)
      f_force[i] = - grad->get(f,i);
  }

  // Put geometry into each efp fragment
  int efp_atom_cnt = 0;
  for (int f=0; f<EFPfrag; ++f) {
    double **f_geom = efp_fragments[f]->get_xyz_pointer();
    double *tmp_geom  = p_efp->get_frag_atom_coord(f);

    for(int i=0; i<3; i++) // first 3 atoms wanted
      for(int xyz=0; xyz<3; ++xyz)
        f_geom[i][xyz] = tmp_geom[3*i+xyz];

    ++efp_atom_cnt;
    delete [] tmp_geom;
  }

/*
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
*/


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

fprintf(outfile,"Preparing to write geometry\n");

  double **geom_2D = g_geom_2D();
fprintf(outfile,"New coordinates for all atoms:\n");
print_matrix(outfile, geom_2D, efp_fragments.size()*3, 3); fflush(outfile);

  // Make compact QM geometry matrix
  int qm_natom = g_qm_natom();
  if (qm_natom > 0) {
    double **qm_geom = init_matrix(qm_natom,3);

    for (int a=0; a<qm_natom; ++a)
      for (int xyz=0; xyz<3; ++xyz)
        qm_geom[a][xyz] = geom_2D[a][xyz];

    psi::Process::environment.molecule()->set_geometry(qm_geom);
    psi::Process::environment.molecule()->update_geometry();

    free_matrix(qm_geom);
  }

  // Make compact EFP geometry array
  // Put new EFP geometry into environment efp pointer.
  if (efp_fragments.size()) {
    double *g = init_array(9*efp_fragments.size()); 

    int cnt = -1;
    for (int a=0; a<3*efp_fragments.size(); ++a)
      for (int xyz=0; xyz<3; ++xyz)
        g[++cnt] = geom_2D[qm_natom+a][xyz];

    fprintf(outfile,"\tNew coordinates for EFP fragments:\n");
    for (int f=0; f<efp_fragments.size(); ++f)
      for (int a=0; a<3; ++a)
        fprintf(outfile,"%15.10lf %15.10lf %15.10lf\n", g[9*f+3*a+0], g[9*f+3*a+1], g[9*f+3*a+2]);

    p_efp->set_coordinates(1,g); // type is POINTS
    free_array(g);
  }
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

