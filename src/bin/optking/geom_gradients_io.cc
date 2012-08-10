/*! \file    geom_gradients_io.cc
    \ingroup optking
    \brief   functions to read and write the geometry, the gradient and the Hessian
*/

#include <cstring>
#include <sstream>
#include "molecule.h"

#define EXTERN
#include "globals.h"
#include <libparallel/parallel.h>


#include "io.h"

#if defined(OPTKING_PACKAGE_PSI)
#include <psi4-dec.h>
#include <libmints/molecule.h>
#include <libmints/matrix.h>
#include <libmints/wavefunction.h>
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

#if defined(OPTKING_PACKAGE_PSI)

  using namespace psi;
    
  SharedMatrix pgradient;
  if (psi::Process::environment.wavefunction()) {
    pgradient = psi::Process::environment.wavefunction()->gradient();
  } else {
    pgradient = psi::Process::environment.gradient();
  }    

  Matrix& gradient = *pgradient.get();

  boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();
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
    fin.open(FILENAME_GEOM_GRAD_IN, ios_base::in);
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
    psi::Communicator::world->sync();
    fin.close();
  } // end try reading geometries
  catch (std::ios_base::failure & bf) {
    printf("Error reading molecular geometry and gradient\n");
    fprintf(outfile,"Error reading molecular geometry and gradient\n");
    throw(INTCO_EXCEPT("Error reading molecular geometry and gradient"));
  }
#endif

#elif defined(OPTKING_PACKAGE_QCHEM)

  double *QX;
  INTEGER *QZ, QNATOMS;
  bool Qnoghosts = true;

  ::get_carts(NULL, &QX, &QZ, &QNATOMS, Qnoghosts);

  if (QNATOMS != g_natom()) {
    fprintf(outfile,"read_geom_grad() QNATOMS=%d\n", QNATOMS);
    //QCrash("Number of atoms read inconsistent with REM variable.");
  }
  //fprintf(outfile, "QX read with get_carts()\n");
  //print_array(outfile, QX, 3*QNATOMS);

  double* QGrad = init_array(3*QNATOMS);
  ::FileMan_Open_Read(FILE_NUCLEAR_GRADIENT);
  ::FileMan(FM_READ, FILE_NUCLEAR_GRADIENT, FM_DP, 3*QNATOMS, 0, FM_BEG, QGrad);
  ::FileMan_Close(FILE_NUCLEAR_GRADIENT);
  //fprintf(outfile, "QGrad read with get_carts()\n");
  //print_array(outfile, QGrad, 3*QNATOMS);


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
  free_array(QGrad);

  double QEnergy;
  ::FileMan_Open_Read(FILE_ENERGY);
  ::FileMan(FM_READ, FILE_ENERGY, FM_DP, 1, FILE_POS_CRNT_TOTAL_ENERGY, FM_BEG, &QEnergy);
  ::FileMan_Close(FILE_ENERGY);
  energy = QEnergy;

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
  psi::Process::environment.molecule()->set_geometry(geom_2D);
  psi::Process::environment.molecule()->update_geometry();
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
    if_Hcart.open(FILENAME_CARTESIAN_H, ios_base::in);
    int n;
    if_Hcart >> n; // read natom
    if_Hcart >> n; // read natom*6 (?)

    for (int i=0; i<Ncart; ++i)
      for (int j=0; j<Ncart; ++j)
        if_Hcart >> H_cart[i][j];

    psi::Communicator::world->sync();
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

  return H_cart;
}

}

