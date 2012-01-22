/*! \file
    \ingroup OPTKING
    \brief fd_freq_0(): compute frequencies from energies
*/

#include "findif.h"

#include <boost/python.hpp>
#include <boost/python/list.hpp>
using namespace boost::python;

#include <physconst.h>

extern "C" {

#define F_DGEMM dgemm_

extern void F_DGEMM(char *transa, char *transb, int *m, int *n, int *k,
  double *alpha, double *A, int *lda, double *B, int *ldb,
  double *beta, double *C, int *ldc);

void opt_matrix_mult(double **A, bool tA, double **B, bool tB, double **C, bool tC,
    int nr, int nl, int nc, bool add);
}

namespace psi { namespace findif {

int iE0(std::vector<int> & Ndisp_pi, std::vector< std::vector<int> > & salcs_pi, int pts,
 int irrep, int ii, int jj, int disp_i, int disp_j);

PsiReturnType fd_freq_0(Options &options, const boost::python::list& python_energies, int freq_irrep_only)
{
  int pts = options.get_int("POINTS");
  double disp_size = options.get_double("DISP_SIZE");
  int print_lvl = 3;

  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();
  int Natom = mol->natom();
  boost::shared_ptr<MatrixFactory> fact;
  CdSalcList salc_list(mol, fact, 0xF, true, true);
  int Nirrep = salc_list.nirrep();
  int Nsalc_all = salc_list.ncd();

  // build vectors that list indices of salcs for each irrep
  std::vector< std::vector<int> > salcs_pi;
  for (int h=0; h<Nirrep; ++h)
    salcs_pi.push_back( std::vector<int>() );
  for (int i=0; i<Nsalc_all; ++i)
    salcs_pi[salc_list[i].irrep()].push_back(i);

  // Now remove irreps that are not requested
  if (freq_irrep_only != -1) {
    for (int h=0; h<Nirrep; ++h)
      if (h != freq_irrep_only)
        salcs_pi[h].clear();
  }

  // count displacements
  std::vector<int> Ndisp_pi (Nirrep);
  // diagonal for symmetric coordinates
  if (pts == 3)
    Ndisp_pi[0] = 2 * salcs_pi[0].size();
  else if (pts == 5)
    Ndisp_pi[0] = 4 * salcs_pi[0].size();
  // diagonal for asymmetric coordinates
  for (int h=1; h<Nirrep; ++h) {
    if (pts == 3)
      Ndisp_pi[h] = salcs_pi[h].size();
    else if (pts == 5)
      Ndisp_pi[h] = 2* salcs_pi[h].size();
  }
  // off-diagonal
  for (int h=0; h<Nirrep; ++h) {
    if (pts == 3)
      Ndisp_pi[h] += 2 * salcs_pi[h].size() * (salcs_pi[h].size() - 1) / 2;
    else if (pts == 5)
      Ndisp_pi[h] += 8 * salcs_pi[h].size() * (salcs_pi[h].size() - 1) / 2;
  }
  int Ndisp_all = 0;
  for (int h=0; h<Nirrep; ++h)
    Ndisp_all += Ndisp_pi[h];

  // Check number of energies and displacements
  std::vector<double> E;
  for (int i=0; i<len(python_energies); ++i)
    E.push_back( (double)extract<double>(python_energies[i]) );

  fprintf(outfile,"\n-------------------------------------------------------------\n\n");

  fprintf(outfile, "  Computing second-derivative from energies using projected, \n");
  fprintf(outfile, "  symmetry-adapted, cartesian coordinates (fd_freq_0).\n");

  fprintf(outfile, "\t%d energies passed in, including the reference energy.\n", (int) E.size());
  if (E.size() != Ndisp_all+1) { // last energy is the reference non-displaced energy
    throw PsiException("FINDIF: Incorrect number of energies passed in!",__FILE__,__LINE__);  }

  double energy_ref = E[Ndisp_all];
  fprintf(outfile, "\tUsing %d-point formula.\n", pts);
  fprintf(outfile, "\tEnergy without displacement: %15.10lf\n", energy_ref);
  fprintf(outfile, "\tCheck energies below for precision!\n");
  for (int i=0; i<Ndisp_all+1; ++i)
    fprintf(outfile,"\t%5d : %20.10lf\n", i+1, E[i]);
  fprintf(outfile,"\n");

  // Determine the number of translation and rotational coordinates projected out
  // and obtain them.  might be needed for cartesian hessian.
  // SharedMatrix rot_trans_out = salc_list.matrix_projected_out();
  //int num_projected_out = rot_trans_out->nrow();

  char **irrep_lbls = mol->irrep_labels();
  double **H = block_matrix(Nsalc_all, Nsalc_all);

  std::vector<VIBRATION *> modes;

  for (int h=0; h<Nirrep; ++h) {

    if (salcs_pi[h].size() == 0) continue;

    double **H_irr = block_matrix(salcs_pi[h].size(),salcs_pi[h].size());

    // do diagonal displacements
    for (int i=0; i<salcs_pi[h].size(); ++i) { // loop over salcs of this irrep

      if (h == 0) { // symmetric
        if (pts == 3) {
          H_irr[i][i] = ( + E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0, +1, 0)]
                          + E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0, -1, 0)]
                                - 2.0 * energy_ref) / (disp_size*disp_size);
        }
        else if (pts == 5) {
          H_irr[i][i] = (
            -  1.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0, -2, 0)]
            + 16.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0, -1, 0)]
            + 16.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0,  1, 0)]
            -  1.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0,  2, 0)]
            - 30.0 * energy_ref ) / (12.0*disp_size*disp_size);
        }
      }
      else {  // asymmetric
        if (pts == 3)
          H_irr[i][i] = 2.0 * (E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0, -1, 0)] - energy_ref) /
            (disp_size * disp_size);
        else if (pts == 5)
          H_irr[i][i] = (
            -  2.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0, -2, 0)]
            + 32.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0, -1, 0)]
            - 30.0 * energy_ref ) / (12.0 * disp_size * disp_size);
      }
    } // i, salc_i

    // off-diagonal displacements
    for (int i=0; i<salcs_pi[h].size(); ++i) { // loop over salcs of this irrep

      for (int j=0; j<i; ++j) {        // loop over salcs of this irrep

        if (pts == 3) {
          H_irr[i][j] = H_irr[j][i] = (
            + E[iE0(Ndisp_pi, salcs_pi, pts, h, i, j, +1, +1)]
            + E[iE0(Ndisp_pi, salcs_pi, pts, h, i, j, -1, -1)]
            + 2.0 * energy_ref
            - E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0, +1, 0)]
            - E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0, -1, 0)]
            - E[iE0(Ndisp_pi, salcs_pi, pts, h, j, 0, +1, 0)]
            - E[iE0(Ndisp_pi, salcs_pi, pts, h, j, 0, -1, 0)]
            ) / (2.0*disp_size*disp_size) ;
        }
        else if (pts == 5) {
          H_irr[i][j] = H_irr[j][i] = (
            - 1.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, j, -1, -2)]
            - 1.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, j, -2, -1)]
            + 9.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, j, -1, -1)]
            - 1.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, j, +1, -1)]
            - 1.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, j, -1,  1)]
            + 9.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, j, +1, +1)]
            - 1.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, j, +2, +1)]
            - 1.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, j, +1, +2)]
            + 1.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0, -2,  0)]
            - 7.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0, -1,  0)]
            - 7.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0, +1,  0)]
            + 1.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, i, 0, +2,  0)]
            + 1.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, j, 0, -2,  0)]
            - 7.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, j, 0, -1,  0)]
            - 7.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, j, 0, +1,  0)]
            + 1.0 * E[iE0(Ndisp_pi, salcs_pi, pts, h, j, 0, +2,  0)]
            + 12.0 * energy_ref) / (12.0 * disp_size * disp_size);
        }
      } // j, salc_j
    } // i, salc_i

    if (print_lvl >= 3) {
      fprintf(outfile, "\n\tForce Constants for irrep %s in mass-weighted, ", irrep_lbls[h]);
      fprintf(outfile, "symmetry-adapted cartesian coordinates.\n");
      mat_print(H_irr, salcs_pi[h].size(), salcs_pi[h].size(), outfile);
    }

    // Build complete hessian for transformation later to cartesians
    int start_irr = 0;
    for (int i=0; i<h; ++i)
      start_irr += salcs_pi[i].size();

    for (int i=0; i<salcs_pi[h].size(); ++i)
      for (int j=0; j<=i; ++j)
        H[start_irr + i][start_irr + j] = H[start_irr + j][start_irr + i] = H_irr[i][j];

    // diagonalize force constant matrix
    int dim = salcs_pi[h].size(); 
    double *evals= init_array(dim);
    double **evects = block_matrix(dim, dim);

    sq_rsp(dim, dim, H_irr, evals, 3, evects, 1e-14);

    // Build Bu^1/2 matrix for this irrep
    SharedMatrix B_irr_shared = salc_list.matrix_irrep(h);
    double **B_irr = B_irr_shared->pointer();

    for (int i=0; i<dim; ++i) 
      for (int a=0; a<Natom; ++a)
        for (int xyz=0; xyz<3; ++xyz)
          B_irr[i][3*a+xyz] /= sqrt(mol->mass(a));

    double **normal_irr = block_matrix(3*Natom, dim);
    C_DGEMM('t', 'n', 3*Natom, dim, dim, 1.0, B_irr[0], 3*Natom, evects[0],
      dim, 0, normal_irr[0], dim);

    if (print_lvl >= 3) {
      fprintf(outfile,"\n\tNormal coordinates (mass-weighted) for irrep %s:\n", irrep_lbls[h]);
      eivout(normal_irr, evals, 3*Natom, dim, outfile);
    }

    for (int i=0; i<salcs_pi[h].size(); ++i) {
      double *v = init_array(3*Natom);
      for (int x=0; x<3*Natom; ++x)
        v[x] = normal_irr[x][i];
      VIBRATION *vib = new VIBRATION(h, evals[i], v);
      modes.push_back(vib);
    }

    free(evals);
    free_block(evects);
    free_block(H_irr);
    free_block(normal_irr);
  }

  // This print function also saves frequencies in wavefunction.
  print_vibrations(modes);

  for (int i=0; i<modes.size(); ++i)
    delete modes[i];
  modes.clear();

  // Transform Hessian into cartesian coordinates

  if (print_lvl >= 3) {
    fprintf(outfile, "\n\tFull force constant matrix in mass-weighted SALCS.\n");
    mat_print(H, Nsalc_all, Nsalc_all, outfile);
  }

  int dim = Nsalc_all;
  // Build Bu^-1/2 matrix for the whole Hessian
  SharedMatrix B_shared = salc_list.matrix();
  double **B = B_shared->pointer();

  for (int i=0; i<dim; ++i) 
    for (int a=0; a<Natom; ++a)
      for (int xyz=0; xyz<3; ++xyz)
        B[i][3*a+xyz] *= sqrt(mol->mass(a));

  double **Hx = block_matrix(3*Natom, 3*Natom);

  // Hx = Bt H B
  for (int i=0; i<dim; ++i)
    for (int j=0; j<dim; ++j)
      for (int x1=0; x1<3*Natom; ++x1)
        for (int x2=0; x2<3*Natom; ++x2)
          Hx[x1][x2] += B[i][x1] * H[i][j] * B[j][x2];

  //for (int x1=0; x1<3*Natom; ++x1)
    //for (int x2=0; x2 < x1; ++x2)
      //Hx[x2][x1] = Hx[x1][x2];

  free_block(H);

  fprintf(outfile, "\n\tForce Constants in mass-weighted cartesian coordinates assuming zero gradient.\n");
  mat_print(Hx, 3*Natom, 3*Natom, outfile);

  // Compute the gradient. It's needed for the (probably) smaller derivative B term 
  // that vanishes at stationary points.

  // Gradient in mass-weighted SALCS
  double *g_q = init_array(salcs_pi[0].size());
  if (pts == 3) {
    for (int i=0; i<salcs_pi[0].size(); ++i)
      g_q[i] = ( E[iE0(Ndisp_pi, salcs_pi, pts, 0, i, 0, +1,  0)]
               - E[iE0(Ndisp_pi, salcs_pi, pts, 0, i, 0, -1,  0)]) / (2.0*disp_size);
  }
  else if (pts == 5) {
    for (int i=0; i<salcs_pi[0].size(); ++i)
      g_q[i] = (+1.0*E[iE0(Ndisp_pi, salcs_pi, pts, 0, i, 0, -2,  0)]
                -8.0*E[iE0(Ndisp_pi, salcs_pi, pts, 0, i, 0, -1,  0)]
                +8.0*E[iE0(Ndisp_pi, salcs_pi, pts, 0, i, 0, +1,  0)]
                -1.0*E[iE0(Ndisp_pi, salcs_pi, pts, 0, i, 0, +2,  0)]) / (12.0*disp_size);
  }

  fprintf(outfile,"Gradient in mass-weighted SALC coordinates:\n");
  for (int i=0; i<salcs_pi[0].size(); ++i)
    fprintf(outfile," %d  %20.10lf\n", i, g_q[i]);

  // linear, gradient term is H_kj += dE/dq . d2q/(dx_k dx_j); 
  // g_q is in mass-weighted SALCS
  // displacements done by displace_cart are also in mass-weighted coordinates
  double **ref_geom = mol->geometry().to_block_matrix();

  double dq2dx2;
  const double DX = 0.01;

  for (int atom_a=0; atom_a<Natom; ++atom_a) {
    for (int xyz_a=0; xyz_a<3; ++xyz_a) {

      ref_geom[atom_a][xyz_a] += DX / sqrt(mol->mass(atom_a));
      mol->set_geometry(ref_geom);
      CdSalcList B_p_list(mol, fact, 0x1, true, true);
      SharedMatrix B_p = B_p_list.matrix();
      //mass_weight_columns_plus_one_half(B_p);

      ref_geom[atom_a][xyz_a] += DX / sqrt(mol->mass(atom_a));
      mol->set_geometry(ref_geom);
      CdSalcList B_p2_list(mol, fact, 0x1, true, true);
      SharedMatrix B_p2 = B_p2_list.matrix();
      //mass_weight_columns_plus_one_half(B_p2);

      ref_geom[atom_a][xyz_a] -= 3.0*DX / sqrt(mol->mass(atom_a));
      mol->set_geometry(ref_geom);
      CdSalcList B_m_list(mol, fact, 0x1, true, true);
      SharedMatrix B_m = B_m_list.matrix();
      //mass_weight_columns_plus_one_half(B_m);

      ref_geom[atom_a][xyz_a] -= 1.0*DX / sqrt(mol->mass(atom_a));
      mol->set_geometry(ref_geom);
      CdSalcList B_m2_list(mol, fact, 0x1, true, true);
      SharedMatrix B_m2 = B_m2_list.matrix();
      //mass_weight_columns_plus_one_half(B_m2);

      ref_geom[atom_a][xyz_a] += 2.0*DX / sqrt(mol->mass(atom_a)); // restore to orig

      for (int atom_b=0; atom_b<Natom; ++atom_b)
        for (int xyz_b=0; xyz_b<3; ++xyz_b) {

          for (int i=0; i<salcs_pi[0].size(); ++i) {
              dq2dx2 = (B_m2->get(i,3*atom_b+xyz_b) - 8.0*B_m->get(i,3*atom_b+xyz_b)
                    +8.0*B_p->get(i,3*atom_b+xyz_b) -    B_p2->get(i,3*atom_b+xyz_b)) / 
                   (12.0*DX); //* sqrt(mol->mass(atom_b));
                   //fprintf(outfile,"%20.10lf\n", dq2dx2);

      //        Hx[3*atom_a+xyz_a][3*atom_b+xyz_b] += g_q[i] * dq2dx2;
          }
        }
      }
    }
  free_block(ref_geom);

  fprintf(outfile, "\n\tForce Constants in cartesian coordinates including derivative term b/f un-mw.\n");
  mat_print(Hx, 3*Natom, 3*Natom, outfile);

  // Un-mass-weight Hessian
/*
  for (int atom_a=0; atom_a<Natom; ++atom_a)
    for (int xyz_a=0; xyz_a<3; ++xyz_a)
      for (int atom_b=0; atom_b<Natom; ++atom_b)
        for (int xyz_b=0; xyz_b<3; ++xyz_b)
      Hx[3*atom_a+xyz_a][3*atom_b+xyz_b] *= sqrt(mol->mass(atom_a)) * sqrt(mol->mass(atom_b));
*/

  //for (int x1=0; x1<3*Natom; ++x1)
    //for (int x2=0; x2<3*Natom; ++x2)

  for (int atom_a=0; atom_a<Natom; ++atom_a)
    for (int xyz_a=0; xyz_a<3; ++xyz_a)
      for (int atom_b=0; atom_b<=atom_a; ++atom_b)
        for (int xyz_b=0; xyz_b<3; ++xyz_b)
          if (fabs(Hx[3*atom_a+xyz_a][3*atom_b+xyz_b] - Hx[3*atom_b+xyz_b][3*atom_a+xyz_a]) > 1e-10)
            printf("assymetry (%d,%d;%d,%d) vs transpose: %15.10f\n", atom_a,xyz_a,atom_b,xyz_b,
              Hx[3*atom_a+xyz_a][3*atom_b+xyz_b] - Hx[3*atom_b+xyz_b][3*atom_a+xyz_a]);
         

  fprintf(outfile, "\n\tForce Constants in cartesian coordinates including derivative term.\n");
  mat_print(Hx, 3*Natom, 3*Natom, outfile);

  free(g_q);

/*
  SharedMatrix B_sym_shared = salc_list.matrix_irrep(0);
  double **B_sym = B_sym_shared->pointer();

  for (int i=0; salcs_pi[0].size(); ++i) 
    for (int a=0; a<Natom; ++a)
      for (int xyz=0; xyz<3; ++xyz)
        B_sym[i][3*a+xyz] *= sqrt(mol->mass(a));

  // Gradient in ordinary cartesians g_q B -> g_cart
  double *g_cart = init_array(3*Natom);
  C_DGEMM('n', 'n', 1, 3*Natom, salcs_pi[0].size(), 1.0, g_q, salcs_pi[0].size(),
    B_sym[0], 3*Natom, 0, g_cart, 3*Natom);

  fprintf(outfile,"Forces in cartesian coordinates\n");
  for (int i=0; i<3*Natom; ++i)
    fprintf(outfile,"\t%15.10lf\n",g_cart[i]);
  free(g_cart);
*/


/*
  double **tmat = block_matrix(3*Natom, dim);
  //C_DGEMM('t', 'n', 3*Natom, dim, dim, 1.0, B[0], 3*Natom, H[0], dim, 0, tmat[0], dim);

  opt_matrix_mult(B, 1, H, 0, tmat, 0, 3*Natom, dim, dim, 0);

  //C_DGEMM('n', 'n', 3*Natom, 3*Natom, dim, 1.0, tmat[0], dim, B[0], 3*Natom, 0, Hx[0], 3*Natom);

  opt_matrix_mult(tmat, 0, B, 0, Hx, 0, 3*Natom, dim, 3*Natom, 0);

  // un mass weight the thing
  for (int i=0; i<Natom; ++i) {
    double mass_i = mol->mass(i);
    for (int j=0; j<Natom; ++j) {
      double mass_j = mol->mass(j);
      for (int i_xyz=0; i_xyz<3; ++i_xyz)
        for (int j_xyz=0; j_xyz<3; ++j_xyz)
          Hx[3*i+i_xyz][3*j+j_xyz] /= sqrt(mass_i) * sqrt(mass_j);
    }
  }
*/
/*
  // Do I need to project this thing?
  double **Px = block_matrix(3*Natom, 3*Natom);
  for (int i=0; i<3*Natom; ++i)
    Px[i][i] = 1.0;

  double **R = block_matrix(3*Natom, 3*Natom);

  fprintf(outfile, "\n\tProjected Force Constants in cartesian coordinates.\n");
  mat_print(Hx, 3*Natom, 3*Natom, outfile);
*/

  FILE *of_Hx = fopen("psi.file15.dat","w");
  fprintf(of_Hx,"%5d", Natom);
  fprintf(of_Hx,"%5d\n", 6*Natom);

  int cnt = -1;
  for (int i=0; i<3*Natom; ++i) {
    for (int j=0; j<3*Natom; ++j) {
      fprintf(of_Hx, "%20.10lf", Hx[i][j]);
      if (++cnt == 2) {
        fprintf(of_Hx,"\n");
        cnt = -1;
      }
    }
  }

  fclose(of_Hx);
  free_block(Hx);

  fprintf(outfile,"\n-------------------------------------------------------------\n");

  return Success;
}

/* iE0() returns index for the energy of a displacement, according to the order
generated in fd_geoms_freq_0()
ii and jj are coordinates, displaced by quantized steps disp_i and disp_j
disp_i,disp_j are {-1,0,+1} for a three-point formula
disp_j,disp_j are {-2,-1,0,+1,+2} for a five-point formula
It is assumed that ii >= jj .
For diagonal displacements disp_j=0 and jj is arbitrary/meaningless.
*/

int iE0(std::vector<int> & Ndisp_pi, std::vector< std::vector<int> > & salcs_pi, int pts,
  int irrep, int ii, int jj, int disp_i, int disp_j) {

  int ndiag_this_irrep;
  int rval=-1;

  // compute starting location of displacements for this irrep
  int start_irr = 0;
  for (int h=0; h<irrep; ++h)
    start_irr += Ndisp_pi[h];

  if (pts == 3) {
    if (disp_j == 0) {  // diagonal; all diagonals at beginning of irrep 
      if (irrep == 0) {
        if (disp_i == -1)
          rval = 2*ii;   // f(-1, 0)
        else if (disp_i == +1)
          rval = 2*ii+1; // f(+1, 0)
      }
      else if (disp_i == -1 || disp_i == +1)
        rval = start_irr + ii;     // f(+1,0) = f(-1, 0)
    }
    else {    // off_diagonal
      if (irrep == 0)
        ndiag_this_irrep = 2 * salcs_pi[0].size();
      else
        ndiag_this_irrep = salcs_pi[irrep].size();

      int ij_pair = 2 * ((ii*(ii-1))/2 + jj);

      if      (disp_i == +1 && disp_j == +1) rval = start_irr + ndiag_this_irrep + ij_pair;
      else if (disp_i == -1 && disp_j == -1) rval = start_irr + ndiag_this_irrep + ij_pair + 1;
    }
  }
  else if (pts == 5) {
    if (disp_j == 0) {   // diagonal 
      if (irrep == 0) {
        if (disp_i == -2)      rval = start_irr + 4*ii;     // f(-2, 0)
        else if (disp_i == -1) rval = start_irr + 4*ii+1;   // f(-1, 0)
        else if (disp_i ==  1) rval = start_irr + 4*ii+2;   // f(+1, 0)
        else if (disp_i ==  2) rval = start_irr + 4*ii+3;   // f(+2, 0)
      }
      else { // irrep != 0
        if      (disp_i == -2 || disp_i == 2) rval = start_irr + 2*ii;     // f(-2, 0)
        else if (disp_i == -1 || disp_i == 1) rval = start_irr + 2*ii+1;   // f(-1, 0)
      }
    }
    else {   //off-diagonal
      if (irrep == 0)
        ndiag_this_irrep = 4 * salcs_pi[0].size() ;
      else
        ndiag_this_irrep = 2 * salcs_pi[irrep].size() ;

      int ij_pair = 8 * ((ii*(ii-1))/2 + jj);

      if      (disp_i == -1 && disp_j == -2) rval = start_irr + ndiag_this_irrep + ij_pair;
      else if (disp_i == -2 && disp_j == -1) rval = start_irr + ndiag_this_irrep + ij_pair+1;
      else if (disp_i == -1 && disp_j == -1) rval = start_irr + ndiag_this_irrep + ij_pair+2;
      else if (disp_i == +1 && disp_j == -1) rval = start_irr + ndiag_this_irrep + ij_pair+3;
      else if (disp_i == -1 && disp_j == +1) rval = start_irr + ndiag_this_irrep + ij_pair+4;
      else if (disp_i == +1 && disp_j == +1) rval = start_irr + ndiag_this_irrep + ij_pair+5;
      else if (disp_i == +2 && disp_j == +1) rval = start_irr + ndiag_this_irrep + ij_pair+6;
      else if (disp_i == +1 && disp_j == +2) rval = start_irr + ndiag_this_irrep + ij_pair+7;
    }
  }

  if (rval < 0) {
    fprintf(outfile,"Problem finding displaced energy.\n");
    throw PsiException("FINDIF: Problem finding displaced energy.",__FILE__,__LINE__);
  }
  //fprintf(outfile,"irrep: %d, ii: %d, jj: %d, disp_i: %d, disp_j: %d, rval: %d \n",
  //  irrep, ii, jj, disp_i, disp_j, rval);
  return rval;
}
}}

