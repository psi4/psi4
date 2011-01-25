/*! \file    geom_gradients_io.cc
    \ingroup OPTKING
    \brief   functions to read and write the geometry, the gradient and the Hessian
*/

#include <cstring>
#include <sstream>
#include "molecule.h"

#define EXTERN
#include "globals.h"

#include "io.h"

#if defined(OPTKING_PACKAGE_PSI)
#include <libchkpt/chkpt.h>
#include <psi4-dec.h>
#include <libmints/molecule.h>
#include <libmints/wavefunction.h>
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
  // for now read from file11
  std::ifstream fin;

  stringstream sline (stringstream::in | stringstream::out);
  char cline[256];

  try {
    fin.open(FILENAME_GEOM_GRAD_IN, ios_base::in);
    fin.clear();
    fin.getline(cline, 256);   // header row
    fin.getline(cline, 256);   // natom_i and energy row
    sline << cline;
    sline >> natom ;
    fin.close();
  }
  catch (std::ios_base::failure & bf) {
    printf("Error reading number of atoms.\n");
    fprintf(outfile,"Error reading number of atoms.\n");
    throw "Error reading number of atoms.";
  }

#elif defined(OPTKING_PACKAGE_QCHEM)

 // get number of atoms

#endif
  return natom;
}


// Read geometry and gradient into a molecule object.
// The fragments must already exist with the number of atoms in each set
// and with allocated memory.
void MOLECULE::read_geom_grad(void) {

#if defined(OPTKING_PACKAGE_PSI)
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

    int nfrag = fragments.size();
    int nallatom = g_natom();
    int ngeom = nlines / (2*nallatom+2);

    // rewind through old geometries
    fin.clear();
    fin.seekg(0);
    for (int i=0; i< (ngeom-1)*(2*nallatom+2); ++i)
      fin.getline(cline, 256);

    // header line
    fin.getline(cline, 256);

    fin >> junk;
    if (junk != nallatom) throw("wrong number of atoms in text file");
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
    fprintf(outfile,"Error reading molecular geometry and gradient\n");
    throw "Error reading molecular geometry and gradient\n";
  }

  // update interfragment reference points if they exist
  for (int i=0; i<interfragments.size(); ++i)
    interfragments[i]->update_reference_points();

#elif defined(OPTKING_PACKAGE_QCHEM)

  printf("Need to add code to read qchem geometry gradient and energy here\n");
  throw "Need to add code to read qchem geometry gradient and energy here\n";

#endif

  return;
}

void MOLECULE::write_geom(void) {
  double **geom_2D = g_geom_2D();

#if defined(OPTKING_PACKAGE_PSI)
  psi::Process::environment.molecule()->set_geometry(geom_2D);
#elif defined(OPTKING_PACKAGE_QCHEM)

  printf("Need to add code to write new geometry to qchem here\n");
  throw "Need to add code to write new geometry to qchem here\n";

#endif

  free_matrix(geom_2D);
/* double *x = g_geom_array();
    double **geom_2D = init_matrix(g_natom(), 3);
    int cnt = 0;
    for (int i=0; i<g_natom(); ++i)
      for (int xyz=0; xyz < 3; ++xyz)
        geom_2D[i][xyz] = x[cnt++];
    Process::environment.molecule()->set_geometry(geom_2D); */
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

    if_Hcart.close();
  }
  catch (std::ios_base::failure & bf) {
    printf("Error reading cartesian Hessian matrix\n");
    fprintf(outfile,"Error reading cartesian Hessian matrix\n");
    throw "Error reading cartesian Hessian matrix\n";
  }

  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"\tCartesian Hessian matrix read: \n");
    print_matrix(outfile, H_cart, Ncart, Ncart);
    fflush(outfile);
  }
#elif defined(OPTKING_PACKAGE_QCHEM)

#endif

  return H_cart;
}


}

