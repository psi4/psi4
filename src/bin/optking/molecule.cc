/*! \file molecule.cc
    \ingroup OPT10
    \brief molecule class (really, molecular system class)
*/

#include "molecule.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include "linear_algebra.h"
#include "print.h"
#include "atom_data.h"
#include "physconst.h"

#define EXTERN
#include "globals.h"

namespace opt {

bool is_integer(const char *check);

using namespace std;

// constructor reads geometries of fragments and creates fragment objects
// Currently, the assumed format is :
// header row
// natom_i natom_j ... [# in each fragment]  total energy
// atomic_symbol/atomic_number  x y z
// ...
// dE/dx  dE/dy   dE/dz
//
MOLECULE::MOLECULE(ifstream & fin) {
  string lbl;
  stringstream sline (stringstream::in | stringstream::out);
  stringstream sline2 (stringstream::in | stringstream::out);
  char cline[256];

  try {
    fin.clear();
    fin.exceptions(ios_base::failbit);

    fin.getline(cline, 256);   // header row
    fin.getline(cline, 256);   // natom_i and energy row
    sline << cline;

    // cnt number of fragments (integer entries)
    sline >> lbl;
    int nfrag = 0;
    while (is_integer(lbl.c_str())) {
      nfrag++;
      sline >> lbl; 
    }

    // record natoms for each fragment
    sline2 << cline;
    int * natom_i = init_int_array(nfrag);
    int i, frag, nallatom=0;
    for (frag=0; frag<nfrag; ++frag) {
      sline2 >> natom_i[frag];
      nallatom += natom_i[frag];
    }

    // rewind, determine number of geometries
    fin.seekg(0);
    int nlines = 0;
    try {
      //while ( !fin.eof() ) { // ignore exception
      while (fin.peek() != EOF) {
        fin.getline(cline, 256);
        ++nlines;
      }
    }
    catch (std::ios_base::failure & bf) { }
    int ngeom = nlines / (2*nallatom+2);

    fin.clear();
    fin.seekg(0);
    for (i=0; i< (ngeom-1)*(2*nallatom+2); ++i)
      fin.getline(cline, 256);

    fin.getline(cline, 256); // header line

    for (frag=0; frag<nfrag; ++frag)
      fin >> i;
    fin >> energy;

    FRAG * one_frag;

    for (frag=0; frag<nfrag; ++frag) {
      double *Z = init_array(natom_i[frag]);
      double **geom = init_matrix(natom_i[frag], 3);

      for (i=0; i<natom_i[frag]; ++i) {
        // check to see if atomic number or symbol
        if ( isalpha((lbl.c_str())[0]) )
          Z[i] = Element_to_Z(lbl);
        else
          fin >> Z[i];

        fin >> geom[i][0];
        fin >> geom[i][1];
        fin >> geom[i][2];
      }

      one_frag = new FRAG(natom_i[frag], Z, geom);
      fragments.push_back(one_frag);
    } // end loop nfrag

    for (frag=0; frag<nfrag; ++frag) {
      double **grad = init_matrix(natom_i[frag], 3);
      for (i=0; i<natom_i[frag]; ++i) {
        fin >> grad[i][0];
        fin >> grad[i][1];
        fin >> grad[i][2];
      }
      fragments[frag]->set_grad(grad);
      free_matrix(grad);
    }
    free_int_array(natom_i);

  } // end try reading geometries
  catch (std::ios_base::failure & bf) {
    printf("Error reading molecular geometry and gradient\n");
    fprintf(outfile,"Error reading molecular geometry and gradient\n");
    throw "Error reading molecular geometry and gradient\n";
  }

  return;
}


bool is_integer(const char *check) {
  for (int i=0; i<strlen(check); ++i) {
    if (!isdigit(check[i]))
      return false;
  }
  return true;
}


// compute forces in internal coordinates
// forces in internal coordinates, f_q = G_inv B u f_x
// if u is unit matrix, f_q = (BB^T)^(-1) * B f_x
void MOLECULE::forces (void) {
  double *f_x, *temp_arr, **B, **G, **G_inv;
  int Ncart = 3*g_natom();
  int Nintco = g_nintco();

  B = compute_B();
  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile, "B matrix\n");
    print_matrix(outfile, B, Nintco, Ncart);
  }

  f_x = g_grad_array(); // Hartree / bohr
  array_scm(f_x, -1, Ncart); // switch gradient -> forces

  temp_arr = init_array(Nintco);
  opt_matrix_mult(B, 0, &f_x, 1, &temp_arr, 1, Nintco, Ncart, 1, 0);
  free_array(f_x);

  G = init_matrix(Nintco, Nintco);
  opt_matrix_mult(B, 0, B, 1, G, 0, Nintco, Ncart, Nintco, 0);
  free_matrix(B);

  G_inv = symm_matrix_inv(G, Nintco, 1);
  free_matrix(G);

  double * f_q = p_Opt_data->g_forces_pointer();
  opt_matrix_mult(G_inv, 0, &temp_arr, 1, &f_q, 1, Nintco, Nintco, 1, 0);
  free_matrix(G_inv);
  free_array(temp_arr);

  if (Opt_params.print_lvl > 3) {
    fprintf(outfile,"Internal forces in au\n");
    print_matrix(outfile, &f_q, 1, Nintco);

    // test by transforming f_q back to cartesian forces and compare
    B = compute_B();
    temp_arr = init_array(Ncart);
    opt_matrix_mult(B, 1, &f_q, 1, &temp_arr, 1, Ncart, Nintco, 1, 0);
    fprintf(outfile,"Recomputed forces in cartesian coordinates\n");
    print_matrix(outfile, &temp_arr, 1, Ncart);
    free_array(temp_arr);
    free_matrix(B);
  }

  return ;
}

// project redundancies (and constraints) out of forces and Hessian matrix
// add constraints here later
void MOLECULE::project_f_and_H(void) {
  int Nintco = g_nintco();
  int Ncart = 3*g_natom();

  double **G = init_matrix(Nintco, Nintco);
  double **B = compute_B();
  opt_matrix_mult(B, 0, B, 1, G, 0, Nintco, Ncart, Nintco, 0);
  free_matrix(B);

  double **G_inv = symm_matrix_inv(G, Nintco, 1);
  double **P = init_matrix(Nintco, Nintco);
  opt_matrix_mult(G, 0, G_inv, 0, P, 0, Nintco, Nintco, Nintco, 0);
  free_matrix(G);
  free_matrix(G_inv);

  double *f_q = p_Opt_data->g_forces_pointer();
  double **H = p_Opt_data->g_H_pointer();

  // f_q~ = P f_q
  double * temp_arr = init_array(Nintco);
  opt_matrix_mult(P, 0, &f_q, 1, &temp_arr, 1, Nintco, Nintco, 1, 0);
  array_copy(temp_arr, f_q, Nintco);
  free_array(temp_arr);

  // H~ = PHP + 1000(1-P)
  double ** temp_mat = init_matrix(Nintco, Nintco);
  opt_matrix_mult(H, 0, P, 0, temp_mat, 0, Nintco, Nintco, Nintco, 0);
  opt_matrix_mult(P, 0, temp_mat, 0, H, 0, Nintco, Nintco, Nintco, 0);
  free_matrix(temp_mat);

  for (int i=0; i<Nintco; ++i)
    P[i][i] = (P[i][i] - 1) * 2000;

  for (int i=0; i<Nintco; ++i)
    for (int j=0; j<Nintco; ++j)
      H[i][j] -= P[i][j];

  free_matrix(P);
}

void MOLECULE::write_geom(void) {
  FILE *fp_geom = fopen(FILENAME_GEOM_OUT, "w");

  fprintf(fp_geom, "PSI: (\n");
  fprintf(fp_geom, "GEOMETRY = ( \n");
  for (int i=0; i<fragments.size(); ++i)
    fragments[i]->write_geom(fp_geom);
  fprintf(fp_geom, "  )\n");
  fprintf(fp_geom, ")\n");
  fclose(fp_geom);
}

}

#ifdef PSI4
#include <libmints/molecule.h>
#include <libchkpt/chkpt.h>
namespace opt {

void MOLECULE::write_geom_chkpt(void) {
  using namespace psi;
  double *x = g_geom_array();

  double **geom_2D = init_matrix(g_natom(), 3);
  int cnt = 0;
  for (int i=0; i<g_natom(); ++i)
    for (int xyz=0; xyz < 3; ++xyz)
      geom_2D[i][xyz] = x[cnt++];

fprintf(outfile,"Print new geom\n");
print_matrix(outfile, geom_2D, g_natom(), 3);

  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_geom(geom_2D);
  chkpt_close();
}

void MOLECULE::write_geom_to_active_molecule() {
    using namespace psi;
    double *x = g_geom_array();

    double **geom_2D = init_matrix(g_natom(), 3);
    int cnt = 0;
    for (int i=0; i<g_natom(); ++i)
      for (int xyz=0; xyz < 3; ++xyz)
        geom_2D[i][xyz] = x[cnt++];

    // Modify "active molecule" with new geometry
    Process::environment.molecule()->set_geometry(geom_2D);
}

}
#endif

namespace opt {

void MOLECULE::apply_intrafragment_step_limit(double * & dq) {
  int i, f;
  double scale = 1.0;
  double limit = Opt_params.intrafragment_step_limit;

  for (f=0; f<fragments.size(); ++f)
    for (i=0; i<fragments[f]->g_nintco(); ++i)
      if (scale * fabs(dq[g_intco_offset(f)+i]) > limit)
        scale = limit / fabs(dq[g_intco_offset(f)+i]);

  if (scale != 1.0) {
    fprintf(outfile,"\tChange in coordinate exceeds step limit of %10.5lf.\n", limit);
    fprintf(outfile,"\tScaling displacements by %10.5lf\n", scale);

    for (f=0; f<fragments.size(); ++f)
      for (i=0; i<fragments[f]->g_nintco(); ++i)
        dq[g_intco_offset(f)+i] *= scale;
  }
}

// don't let any angles get smaller than 0.0
void MOLECULE::check_intrafragment_zero_angles(double const * const dq) {
  for (int f=0; f<fragments.size(); ++f)
    fragments[f]->check_zero_angles(&(dq[g_intco_offset(f)]));
}

void MOLECULE::H_guess(void) const {
  double **H, **H_frag;
  int f, i, j;

  H = p_Opt_data->g_H_pointer();
  for (f=0; f<fragments.size(); ++f) {
    H_frag = fragments[f]->H_guess();
    for (i=0; i<fragments[f]->g_nintco(); ++i)
      for (j=0; j<fragments[f]->g_nintco(); ++j)
        H[g_intco_offset(f) + i][g_intco_offset(f) + j] = H_frag[i][j];
    free_matrix(H_frag);
  }

  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile,"Initial Hessian guess\n");
    print_matrix(outfile, H, g_nintco(), g_nintco());
  }
  return;
}

}
