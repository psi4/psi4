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
    try { // ignore exceptions only here
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


// compute forces in internal coordinates in au
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
  for (int i=0; i<fragments.size(); ++i)
    fragments[i]->write_geom(outfile);
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

  //fprintf(outfile,"Print new geom\n");
  //print_matrix(outfile, geom_2D, g_natom(), 3);

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
  double **H, **H_frag, **H_interfrag;

  H = p_Opt_data->g_H_pointer();

  for (int f=0; f<fragments.size(); ++f) {
    H_frag = fragments[f]->H_guess();

    for (int i=0; i<fragments[f]->g_nintco(); ++i)
      for (int j=0; j<fragments[f]->g_nintco(); ++j)
        H[g_intco_offset(f) + i][g_intco_offset(f) + j] = H_frag[i][j];

    free_matrix(H_frag);
  }

  for (int I=0; I<interfragments.size(); ++I) {
    H_interfrag = interfragments[I]->H_guess();

    for (int i=0; i<interfragments[I]->g_nintco(); ++i)
      for (int j=0; j<interfragments[I]->g_nintco(); ++j)
        H[g_interfragment_intco_offset(I) + i][g_interfragment_intco_offset(I) + j] =
          H_interfrag[i][j];

    free_matrix(H_interfrag);
  }

  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile,"Initial Hessian guess\n");
    print_matrix(outfile, H, g_nintco(), g_nintco());
  }
  return;
}

double ** MOLECULE::cartesian_H_to_internals(void) const {
  int Nintco = g_nintco();
  int Ncart = 3*g_natom();

  // compute A = u B^t (B u B^t)^-1 where u=unit matrix and -1 is generalized inverse
  double **B = compute_B();
  double **G = init_matrix(Nintco, Nintco);
  opt_matrix_mult(B, 0, B, 1, G, 0, Nintco, Ncart, Nintco, 0);

  double **G_inv = symm_matrix_inv(G, Nintco, true);
  free_matrix(G);

  double **A = init_matrix(Ncart, Nintco);
  opt_matrix_mult(B, 1, G_inv, 0, A, 0, Ncart, Nintco, Nintco, 0);
  free_matrix(G_inv);
  free_matrix(B);

  // compute gradient in internal coordinates, A^t g_x = g_q
  double *grad_x = g_grad_array();
  double *grad_q = init_array(Nintco);
  opt_matrix_mult(A, 1, &grad_x, 1, &grad_q, 1, Nintco, Ncart, 1, 0);
  free_array(grad_x);

  // read in cartesian H
  double **H_cart = p_Opt_data->read_cartesian_H();

  // transform cartesian H to internals; A^t (H_x - K) A
  // K_ij = sum_q ( grad_q[q] d^2(q)/(dxi dxj) )
  double **dq2dx2;

  for (int q=0; q<Nintco; ++q) {
    dq2dx2 = compute_derivative_B(q); // d^2(q)/ dx_i dx_j

    for (int i=0; i<Ncart; ++i)
      for (int j=0; j<Ncart; ++j)
        H_cart[i][j] -= grad_q[q] * dq2dx2[i][j];

    free_matrix(dq2dx2);
  }
  free_array(grad_q);

  double **temp_mat = init_matrix(Ncart, Nintco);
  opt_matrix_mult(H_cart, 0, A, 0, temp_mat, 0, Ncart, Ncart, Nintco, 0);
  free_matrix(H_cart);

  double **H_int = init_matrix(Nintco, Nintco);
  opt_matrix_mult(A, 1, temp_mat, 0, H_int, 0, Nintco, Ncart, Nintco, 0);
  free_matrix(temp_mat);

  free_matrix(A);

  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile, "Hessian transformed to internal coordinates:\n");
    print_matrix(outfile, H_int, Nintco, Nintco); fflush(outfile);
  }

  // Check by transforming internal coordinate Hessian back into cartesian coordinates:
/*
  B = compute_B();
  temp_mat = init_matrix(Ncart, Nintco);
  opt_matrix_mult(B, 1, H_int, 0, temp_mat, 0, Ncart, Nintco, Nintco, 0);
  H_cart = init_matrix(Ncart, Ncart);
  opt_matrix_mult(temp_mat, 0, B, 0, H_cart, 0, Ncart, Nintco, Ncart, 0);
  free_matrix(temp_mat);

  for (int q=0; q<Nintco; ++q) {
    dq2dx2 = compute_derivative_B(q); // d^2(q)/ dx_i dx_j
    for (int i=0; i<Ncart; ++i)
      for (int j=0; j<Ncart; ++j)
        H_cart[i][j] += grad_q[q] * dq2dx2[i][j];
    free_matrix(dq2dx2);
  }
  free_array(grad_q);
  fprintf(outfile, "Hessian transformed back into Cartesian coordinates\n");
  print_matrix(outfile, H_cart, Ncart, Ncart); fflush(outfile);
  free_matrix(B);
*/
  return H_int;
}

double ** MOLECULE::compute_B(void) const {
  double **B, **B_frag, **B_inter;
  int i, j;

  B = init_matrix(g_nintco(), 3*g_natom());

  for (int f=0; f<fragments.size(); ++f) {
    B_frag = fragments[f]->compute_B();

    for (i=0; i<fragments[f]->g_nintco(); ++i)
      for (j=0; j<3*fragments[f]->g_natom(); ++j)
        B[g_intco_offset(f)+i][3*g_atom_offset(f)+j] = B_frag[i][j];

    free_matrix(B_frag);
  }

  for (int I=0; I<interfragments.size(); ++I) {
    B_inter = interfragments[I]->compute_B(); // ->g_nintco() X (3*atom A)+3(natom_B)

    int iA = interfragments[I]->g_A_index();
    int iB = interfragments[I]->g_B_index();
    int nA = interfragments[I]->g_natom_A();
    int nB = interfragments[I]->g_natom_B();

    for (i=0; i<interfragments[I]->g_nintco(); ++i) { // for each of up to 6 coordinates

      for (j=0; j<3*nA; ++j)
        B[g_interfragment_intco_offset(I)+i][3*g_atom_offset(iA)+j] = B_inter[i][j];

      for (j=0; j<3*nB; ++j)
        B[g_interfragment_intco_offset(I)+i][3*g_atom_offset(iB)+j] = B_inter[i][3*nA+j];

    }
    free_matrix(B_inter);
  }
  return B;
}

double ** MOLECULE::compute_derivative_B(int intco_index) const {
  int cnt_intcos = 0;
  int fragment_index = -1;
  int coordinate_index = 0;
  bool is_interfragment = true;
  int f;

  for (f=0; f<fragments.size(); ++f) {
    for (int i=0; i<fragments[f]->g_nintco(); ++i) {
      if (cnt_intcos++ == intco_index) {
        fragment_index = f;
        coordinate_index = i;
        is_interfragment = false;
        break;
      }
    }
  }

  if (is_interfragment) {  // inntco_index not yet found
    for (f=0; f<interfragments.size(); ++f) {
      for (int i=0; i<interfragments[f]->g_nintco(); ++i) {
        if (cnt_intcos++ == intco_index) {
          fragment_index = f;
          coordinate_index = i;
          break;
        }
      }
    }
  }

  if (fragment_index == -1)
    throw("MOLECULE::compute_derivative_B() could not find intco_index");

  double **dq2dx2_frag;

  double **dq2dx2 = init_matrix(3*g_natom(), 3*g_natom());

  if (!is_interfragment) {
    dq2dx2_frag = fragments[fragment_index]->compute_derivative_B(coordinate_index);
    // dimension is 2x2, 3x3 or 4x4 and given by natom_intco
    int natom_intco = fragments[fragment_index]->g_intco_natom(coordinate_index);
    int atom_a, atom_b;

    for (int a=0; a<natom_intco; ++a) {
      atom_a = g_atom_offset(fragment_index) + fragments[fragment_index]->g_intco_atom(coordinate_index, a);
      for (int b=0; b<natom_intco; ++b) {
        atom_b = g_atom_offset(fragment_index) + fragments[fragment_index]->g_intco_atom(coordinate_index, b);

      for (int xyz_a=0; xyz_a<3; ++xyz_a)
        for (int xyz_b=0; xyz_b<3; ++xyz_b)
          dq2dx2[3*atom_a + xyz_a][3*atom_b + xyz_b] = dq2dx2_frag[3*a+xyz_a][3*b+xyz_b];

      }
    }
    free_matrix(dq2dx2_frag);
  }
  else { // interfragment cordinate
    dq2dx2_frag = interfragments[fragment_index]->compute_derivative_B(coordinate_index);

    // dimension is (3*natomA + 3*natomB) x (3*natomA + 3*natomB) -> 3*natom X 3*natom
    int nA = interfragments[fragment_index]->g_natom_A();
    int nB = interfragments[fragment_index]->g_natom_B();
    int iA = 3*g_atom_offset(interfragments[fragment_index]->g_A_index());
    int iB = 3*g_atom_offset(interfragments[fragment_index]->g_B_index());

    //print_matrix(outfile,dq2dx2_frag, 3*(nA+nB), 3*(nA+nB));
    for (int a=0; a<3*nA; ++a) // A-A block
      for (int aa=0; aa<3*nA; ++aa)
        dq2dx2[iA + a][iA + aa] = dq2dx2_frag[a][aa];

    for (int a=0; a<3*nA; ++a) // A-B block
      for (int bb=0; bb<3*nB; ++bb)
        dq2dx2[iA + a][iB + bb] = dq2dx2_frag[a][3*nA + bb];

    for (int b=0; b<3*nB; ++b) // B-A block
      for (int aa=0; aa<3*nA; ++aa)
        dq2dx2[iB + b][iA + aa] = dq2dx2_frag[3*nA + b][aa];

    for (int b=0; b<3*nB; ++b) // B-B block
      for (int bb=0; bb<3*nB; ++bb)
        dq2dx2[iB + b][iB + bb] = dq2dx2_frag[3*nA + b][3*nA + bb];

    free_matrix(dq2dx2_frag);
  }
  return dq2dx2;
}

void MOLECULE::print_intco_dat(FILE *fout) {

  for (int i=0; i<fragments.size(); ++i) {
    int first = g_atom_offset(i);
    fprintf(fout,"%d-%d\n", first+1, first + fragments[i]->g_natom());
    fragments[i]->print_intco_dat(fout, g_atom_offset(i));
  }

  for (int I=0; I<interfragments.size(); ++I) {
    int frag_a = interfragments[I]->g_A_index();
    int frag_b = interfragments[I]->g_B_index();
    fprintf(fout,"F %d %d\n", frag_a+1, frag_b+1);

    for (int i=0; i<6; ++i) 
      fprintf(fout," %d", (int) interfragments[I]->coordinate_on(i));
    fprintf(fout,"\n");

    interfragments[I]->print_intco_dat(fout, g_atom_offset(frag_a), g_atom_offset(frag_b));
  }
}

}


