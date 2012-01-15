/*! \file
    \ingroup OPTKING
    \brief fd_hessian_0(): compute hessian from energies
*/

#include "findif.h"

#include <boost/python.hpp>
#include <boost/python/list.hpp>
using namespace boost::python;

#include <physconst.h>

namespace psi { namespace findif {

int iE0(std::vector<int> & symm_salcs, int pts, int ii, int jj, int disp_i, int disp_j);

PsiReturnType fd_hessian_0(Options &options, const boost::python::list& E_list)
{
  int pts = options.get_int("POINTS");
  double disp_size = options.get_double("DISP_SIZE");
  int print_lvl = options.get_int("PRINT");

  fprintf(outfile,"\n-------------------------------------------------------------\n\n");

  fprintf(outfile, "  Computing totally symmetric second-derivative force constants\n");
  fprintf(outfile, "  from symmetry-adapted, cartesian coordinates (fd_hessian_0).\n");

  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();
  boost::shared_ptr<MatrixFactory> fact;
  // do all for now at least
  CdSalcList salc_list(mol, fact, 0xF, true, true);

  int Natom = mol->natom();
  fprintf(outfile,"\tNumber of atoms is %d.\n", Natom);

  int Nirrep = salc_list.nirrep();
  fprintf(outfile,"\tNumber of irreps is %d.\n", Nirrep);
  
  int Nsalc_all = salc_list.ncd();
  fprintf(outfile,"\tNumber of SALCS is %d.\n", Nsalc_all);

  // build vectors that list indices of salcs for each irrep
  std::vector<int> symm_salcs;
  
  for (int i=0; i<Nsalc_all; ++i)
    if (salc_list[i].irrep() == 0)
      symm_salcs.push_back(i);

  fprintf(outfile,"\tNumber of symmetric SALC's is %d.\n", (int) symm_salcs.size());
  
  int Ndisp;

  // diagonal displacements for symmetric coordinates
  if (pts == 3)
    Ndisp = 2 * symm_salcs.size();
  else if (pts == 5)
    Ndisp = 4 * symm_salcs.size();

  // off-diagonal displacements
  if (pts == 3)
    Ndisp += 2 * symm_salcs.size() * (symm_salcs.size() - 1) / 2;
  else if (pts == 5)
    Ndisp += 8 * symm_salcs.size() * (symm_salcs.size() - 1) / 2;
  
  fprintf(outfile,"\tNumber of symmetric displacements (including reference) is %d.\n", Ndisp);

  if (options.get_int("PRINT") > 1)
    for (int i=0; i<salc_list.ncd(); ++i)
      salc_list[i].print();

  // Check number of energies and displacements
  std::vector<double> E;
  for (int i=0; i<len(E_list); ++i)
    E.push_back( (double)extract<double>(E_list[i]) );

  fprintf(outfile, "\t%d energies passed in, including the reference energy.\n", (int) E.size());
  if (E.size() != Ndisp+1) { // last energy is the reference non-displaced energy
    throw PsiException("FINDIF: Incorrect number of energies passed in!",__FILE__,__LINE__);  }

  double energy_ref = E[Ndisp];
  fprintf(outfile, "\tUsing %d-point formula.\n", pts);
  fprintf(outfile, "\tEnergy without displacement: %15.10lf\n", energy_ref);
  fprintf(outfile, "\tCheck energies below for precision!\n");
  for (int i=0; i<Ndisp+1; ++i)
    fprintf(outfile,"\t%5d : %20.10lf\n", i+1, E[i]);
  fprintf(outfile,"\n");

  //char **irrep_lbls = mol->irrep_labels();
  //std::vector<VIBRATION *> modes;

  double **H = block_matrix(symm_salcs.size(),symm_salcs.size());

  // do diagonal displacements
  for (int i=0; i<symm_salcs.size(); ++i) { // loop over salcs of this irrep

        if (pts == 3) {
          H[i][i] = ( + E[iE0(symm_salcs, pts, i, 0, +1, 0)]
                          + E[iE0(symm_salcs, pts, i, 0, -1, 0)]
                          - 2.0 * energy_ref) / (disp_size*disp_size);
        }
        else if (pts == 5) {
          H[i][i] = (
            -  1.0 * E[iE0(symm_salcs, pts, i, 0, -2, 0)]
            + 16.0 * E[iE0(symm_salcs, pts, i, 0, -1, 0)]
            + 16.0 * E[iE0(symm_salcs, pts, i, 0,  1, 0)]
            -  1.0 * E[iE0(symm_salcs, pts, i, 0,  2, 0)]
            - 30.0 * energy_ref ) / (12.0*disp_size*disp_size);
        }
      }

    // off-diagonal displacements
    for (int i=0; i<symm_salcs.size(); ++i) { // loop over salcs of this irrep

      for (int j=0; j<i; ++j) {        // loop over salcs of this irrep

        if (pts == 3) {
          H[i][j] = H[j][i] = (
            + E[iE0(symm_salcs, pts, i, j, +1, +1)]
            + E[iE0(symm_salcs, pts, i, j, -1, -1)]
            + 2.0 * energy_ref
            - E[iE0(symm_salcs, pts, i, 0, +1, 0)]
            - E[iE0(symm_salcs, pts, i, 0, -1, 0)]
            - E[iE0(symm_salcs, pts, j, 0, +1, 0)]
            - E[iE0(symm_salcs, pts, j, 0, -1, 0)]
            ) / (2.0*disp_size*disp_size) ;
        }
        else if (pts == 5) {
          H[i][j] = H[j][i] = (
            - 1.0 * E[iE0(symm_salcs, pts, i, j, -1, -2)]
            - 1.0 * E[iE0(symm_salcs, pts, i, j, -2, -1)]
            + 9.0 * E[iE0(symm_salcs, pts, i, j, -1, -1)]
            - 1.0 * E[iE0(symm_salcs, pts, i, j, +1, -1)]
            - 1.0 * E[iE0(symm_salcs, pts, i, j, -1,  1)]
            + 9.0 * E[iE0(symm_salcs, pts, i, j, +1, +1)]
            - 1.0 * E[iE0(symm_salcs, pts, i, j, +2, +1)]
            - 1.0 * E[iE0(symm_salcs, pts, i, j, +1, +2)]
            + 1.0 * E[iE0(symm_salcs, pts, i, 0, -2,  0)]
            - 7.0 * E[iE0(symm_salcs, pts, i, 0, -1,  0)]
            - 7.0 * E[iE0(symm_salcs, pts, i, 0, +1,  0)]
            + 1.0 * E[iE0(symm_salcs, pts, i, 0, +2,  0)]
            + 1.0 * E[iE0(symm_salcs, pts, j, 0, -2,  0)]
            - 7.0 * E[iE0(symm_salcs, pts, j, 0, -1,  0)]
            - 7.0 * E[iE0(symm_salcs, pts, j, 0, +1,  0)]
            + 1.0 * E[iE0(symm_salcs, pts, j, 0, +2,  0)]
            + 12.0 * energy_ref) / (12.0 * disp_size * disp_size);
        }
      } // j, salc_j
    } // i, salc_i

    //if (print_lvl >= 3) {
      fprintf(outfile, "\n\tSymmetric Force Constants in mass-weighted, cartesian coordinates.\n");
      mat_print(H, symm_salcs.size(), symm_salcs.size(), outfile);
    //}

    int dim = symm_salcs.size(); 

    // Build Bu^1/2 matrix for this irrep
    double **B = block_matrix(dim, 3*Natom);

    for (int i=0; i<dim; ++i) {
      int salc_i = symm_salcs[i];
      for (int c=0; c<salc_list[salc_i].ncomponent(); ++c) {
        int a          = salc_list[salc_i].component(c).atom;
        int xyz        = salc_list[salc_i].component(c).xyz;
        double coef    = salc_list[salc_i].component(c).coef;
        B[i][3*a+xyz] = coef / sqrt(mol->mass(a));
      }
    }

fprintf(outfile, "\n\tB u ^1/2:\n");
mat_print(B, dim, 3*Natom, outfile);

    double **tmat = block_matrix(3*Natom, dim);

    C_DGEMM('t', 'n', 3*Natom, dim, dim, 1.0, B[0], 3*Natom, H[0], dim, 0, tmat[0], dim);

fprintf(outfile, "\n\tTmat.\n");
mat_print(tmat, 3*Natom, dim, outfile);

    double **Hx = block_matrix(3*Natom, 3*Natom);

    C_DGEMM('n', 'n', 3*Natom, 3*Natom, dim, 1.0, tmat[0], dim, B[0], 3*Natom, 0, Hx[0], 3*Natom);

    fprintf(outfile, "\n\tSymmetric Force Constants in cartesian coordinates.\n");
    mat_print(Hx, 3*Natom, 3*Natom, outfile);

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
    free_block(tmat);

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

int iE0(std::vector<int> & symm_salcs, int pts, int ii, int jj, int disp_i, int disp_j) {

  int ndiag_this_irrep;
  int rval=-1;

  if (pts == 3) {
    if (disp_j == 0) {  // diagonal; all diagonals at beginning of irrep 
        if (disp_i == -1)
          rval = 2*ii;   // f(-1, 0)
        else if (disp_i == +1)
          rval = 2*ii+1; // f(+1, 0)
    }
    else {    // off_diagonal
      ndiag_this_irrep = 2 * symm_salcs.size();

      int ij_pair = 2 * ((ii*(ii-1))/2 + jj);

      if      (disp_i == +1 && disp_j == +1) rval = ndiag_this_irrep + ij_pair;
      else if (disp_i == -1 && disp_j == -1) rval = ndiag_this_irrep + ij_pair + 1;
    }
  }
  else if (pts == 5) {
    if (disp_j == 0) {   // diagonal 
        if (disp_i == -2)      rval = 4*ii;     // f(-2, 0)
        else if (disp_i == -1) rval = 4*ii+1;   // f(-1, 0)
        else if (disp_i ==  1) rval = 4*ii+2;   // f(+1, 0)
        else if (disp_i ==  2) rval = 4*ii+3;   // f(+2, 0)
    }
    else {   //off-diagonal
      ndiag_this_irrep = 4 * symm_salcs.size() ;

      int ij_pair = 8 * ((ii*(ii-1))/2 + jj);

      if      (disp_i == -1 && disp_j == -2) rval = ndiag_this_irrep + ij_pair;
      else if (disp_i == -2 && disp_j == -1) rval = ndiag_this_irrep + ij_pair+1;
      else if (disp_i == -1 && disp_j == -1) rval = ndiag_this_irrep + ij_pair+2;
      else if (disp_i == +1 && disp_j == -1) rval = ndiag_this_irrep + ij_pair+3;
      else if (disp_i == -1 && disp_j == +1) rval = ndiag_this_irrep + ij_pair+4;
      else if (disp_i == +1 && disp_j == +1) rval = ndiag_this_irrep + ij_pair+5;
      else if (disp_i == +2 && disp_j == +1) rval = ndiag_this_irrep + ij_pair+6;
      else if (disp_i == +1 && disp_j == +2) rval = ndiag_this_irrep + ij_pair+7;
    }
  }

  if (rval < 0) {
    fprintf(outfile,"Problem finding displaced energy.\n");
    throw PsiException("FINDIF: Problem finding displaced energy.",__FILE__,__LINE__);
  }
  return rval;
}

}}

