/*! \file
    \ingroup OPTKING
    \brief fd_grad_1_0(): compute gradient using energies and finite-differences
*/

#include <libmints/mints.h>
#include <libqt/qt.h>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>

#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include <fstream> // to write out file15

using namespace boost::python;

namespace psi { namespace findif {

int iE(int *ndisp, int *nsalc, int pts, int irr, int ii, int jj, int disp_i, int disp_j);

PsiReturnType fd_2_0(Options &options, const boost::python::list& E_list)
{
  int pts = options.get_int("POINTS");
  double Dx = options.get_double("DISP_SIZE");

  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();
  int Natom = mol->natom();
  boost::shared_ptr<MatrixFactory> fact;
  int nirreps = mol->point_group()->char_table().nirrep();

  // Compute number of displacements - check with number of energies passed in
  // Determine number of geometries (1 + # of displacements)
  int Ndisp_all;
  // diagonal
  if (pts == 3)
    Ndisp_all = 2 * (3*Natom);
  else if (pts == 5)
    Ndisp_all = 4 * (3*Natom);

  // off-diagonal 
  if (pts == 3)
    Ndisp_all += 2 * (3*Natom) * (3*Natom - 1) / 2;
  else if (pts == 5)
    Ndisp_all += 8 * (3*Natom) * (3*Natom - 1) / 2;
  else
    throw PSIEXCEPTION("fd_2_0: Cannot handle POINTS != (3 or 5).");

  // assume no symmetry for now
  int Ndisp[nirreps];
  for (int h=0; h<nirreps; ++h)
    Ndisp[h] = 0; 
  Ndisp[0] = Ndisp_all;

  int Nsalc[nirreps];
  for (int h=0; h<nirreps; ++h)
    Nsalc[h] = 0; 
  Nsalc[0] = 3*Natom;

  std::vector<double> E;
  for (int i=0; i<len(E_list); ++i)
    E.push_back( (double)extract<double>(E_list[i]) );

  fprintf(outfile, "\t%d energies passed in, including non-displaced energy.\n", (int) E.size());
  fprintf(outfile, "\t%d displaced energies expected.\n", Ndisp_all);

  if (E.size() != Ndisp_all+1) { // last energy is the reference non-displaced energy
    throw PsiException("FINDIF: Incorrect number of energies passed in!",__FILE__,__LINE__);
  }

  for (int i=0; i<Ndisp_all+1; ++i)
    fprintf(outfile,"%d %20.10lf\n", i+1, E[i]);

  double energy_ref = E[Ndisp_all];
  fprintf(outfile, "\tFinite difference computation of second-derivative of energy with respect\n");
  fprintf(outfile, "\tto cartesian coordinates using %d-point formula.  Check for precision!\n", pts);
  fprintf(outfile, "\tEnergy without displacement: %15.10lf\n", energy_ref);

  double **H = block_matrix(3*Natom, 3*Natom);

  int h=0; //irrep

  if (pts == 3) {

    for (int i=0; i<3*Natom; ++i) { // diagonal
      H[i][i] = ( + E[iE(Ndisp, Nsalc, pts, h, i, 0, +1, 0)]
                  + E[iE(Ndisp, Nsalc, pts, h, i, 0, -1, 0)]
               - 2.0 * energy_ref) / (Dx*Dx);
    }

    for (int i=0; i<3*Natom; ++i) { // off-diagonal
      for (int j=0; j<i; ++j) {
        H[i][j] = H[j][i] = (
          + E[iE(Ndisp, Nsalc, pts, h, i, j, +1, +1)]
          + E[iE(Ndisp, Nsalc, pts, h, i, j, -1, -1)]
          + 2.0 * energy_ref
          - E[iE(Ndisp, Nsalc, pts, h, i, 0, +1, 0)]
          - E[iE(Ndisp, Nsalc, pts, h, i, 0, -1, 0)]
          - E[iE(Ndisp, Nsalc, pts, h, j, 0, +1, 0)]
          - E[iE(Ndisp, Nsalc, pts, h, j, 0, -1, 0)]
          ) / (2.0*Dx*Dx) ;
      }
    }

  }
  else if (pts == 5) {

    for (int i=0; i<3*Natom; ++i) { // diagonal
      H[i][i] = (
        -  1.0 * E[iE(Ndisp, Nsalc, pts, h, i, 0, -2, 0)]
        + 16.0 * E[iE(Ndisp, Nsalc, pts, h, i, 0, -1, 0)]
        + 16.0 * E[iE(Ndisp, Nsalc, pts, h, i, 0,  1, 0)]
        -  1.0 * E[iE(Ndisp, Nsalc, pts, h, i, 0,  2, 0)]
        - 30.0 * energy_ref ) / (12.0*Dx*Dx);
    }

    for (int i=0; i<3*Natom; ++i) { // off-diagonal
      for (int j=0; j<i; ++j) {
        H[i][j] = H[j][i] = (
          - 1.0 * E[iE(Ndisp, Nsalc, pts, h, i, j, -1, -2)]
          - 1.0 * E[iE(Ndisp, Nsalc, pts, h, i, j, -2, -1)]
          + 9.0 * E[iE(Ndisp, Nsalc, pts, h, i, j, -1, -1)]
          - 1.0 * E[iE(Ndisp, Nsalc, pts, h, i, j, +1, -1)]
          - 1.0 * E[iE(Ndisp, Nsalc, pts, h, i, j, -1,  1)]
          + 9.0 * E[iE(Ndisp, Nsalc, pts, h, i, j, +1, +1)]
          - 1.0 * E[iE(Ndisp, Nsalc, pts, h, i, j, +2, +1)]
          - 1.0 * E[iE(Ndisp, Nsalc, pts, h, i, j, +1, +2)]
          + 1.0 * E[iE(Ndisp, Nsalc, pts, h, i, 0, -2,  0)]
          - 7.0 * E[iE(Ndisp, Nsalc, pts, h, i, 0, -1,  0)]
          - 7.0 * E[iE(Ndisp, Nsalc, pts, h, i, 0, +1,  0)]
          + 1.0 * E[iE(Ndisp, Nsalc, pts, h, i, 0, +2,  0)]
          + 1.0 * E[iE(Ndisp, Nsalc, pts, h, j, 0, -2,  0)]
          - 7.0 * E[iE(Ndisp, Nsalc, pts, h, j, 0, -1,  0)]
          - 7.0 * E[iE(Ndisp, Nsalc, pts, h, j, 0, +1,  0)]
          + 1.0 * E[iE(Ndisp, Nsalc, pts, h, j, 0, +2,  0)]
          + 12.0 * energy_ref) / (12.0*Dx*Dx);
      }
    }

  } //pts == 5

  FILE *fp_file15 = fopen("psi.file15.dat", "w");

  fprintf(fp_file15, "%5d%5d\n", Natom, 6*Natom);

  int cnt=0;
  for (int i=0; i<3*Natom; ++i) {
    for (int j=0; j<3*Natom; ++j, ++cnt) {
      fprintf(fp_file15,"%20.10lf", H[i][j]);
      if (cnt == 2) { 
        fprintf(fp_file15, "\n");
        cnt = -1;
      }
    }
  }

  fclose(fp_file15);

  free_block(H);

  return Success;
}

/* iE() returns index for the energy of a displacement, according to the order
generated in fd_geoms_2_0()
ii and jj are coordinates, displaced by quantized steps disp_i and disp_j
disp_i,disp_j are {-1,0,+1} for a three-point formula
disp_j,disp_j are {-2,-1,0,+1,+2} for a five-point formula
It is assumed that ii >= jj .
For diagonal displacements disp_j=0 and jj is arbitrary/meaningless.
*/

int iE(int *ndisp, int *nsalc, int pts, int irr, int ii, int jj, int disp_i, int disp_j) {

  if (pts == 3) {
    if (disp_j == 0) {          // diagonal 
      if (disp_i == -1)
        return (2*ii);          // f(-1, 0)
      else if (disp_i == 1)
        return (2*ii+1);        // f(+1, 0)
    }
    else {                      // off_diagonal
      int ndiag = 2 * (*nsalc);

      int ij_pair = ndiag + 2 * ((ii*(ii-1))/2 + jj);

      if ((disp_i == +1) && (disp_j == +1)) // f(+1,+1)
        return(ij_pair + 0);
      else if ((disp_i == -1) && (disp_j == -1)) // f(-1,-1)
        return(ij_pair + 1);
      else
        throw PsiException("FINDIF: Displacements invalid.",__FILE__,__LINE__);
    }

  }
  else if (pts == 5) {

    if (disp_j == 0) {   // diagonal 
      if (disp_i == -2)      return (4*ii);     // f(-2, 0)
      else if (disp_i == -1) return (4*ii+1);   // f(-1, 0)
      else if (disp_i ==  1) return (4*ii+2);   // f(+1, 0)
      else if (disp_i ==  2) return (4*ii+3);   // f(+2, 0)
    }
    else {              //off-diagonal
      int ndiag = 4 * (*nsalc);

      int ij_pair = ndiag + 8 * ((ii*(ii-1))/2 + jj);

      if      ((disp_i == -1) && (disp_j == -2)) return(ij_pair);
      else if ((disp_i == -2) && (disp_j == -1)) return(ij_pair+1);
      else if ((disp_i == -1) && (disp_j == -1)) return(ij_pair+2);
      else if ((disp_i == +1) && (disp_j == -1)) return(ij_pair+3);
      else if ((disp_i == -1) && (disp_j == +1)) return(ij_pair+4);
      else if ((disp_i == +1) && (disp_j == +1)) return(ij_pair+5);
      else if ((disp_i == +2) && (disp_j == +1)) return(ij_pair+6);
      else if ((disp_i == +1) && (disp_j == +2)) return(ij_pair+7);
      else
        throw PsiException("FINDIF: Displacements invalid.",__FILE__,__LINE__);
    }

  } //pts == 5

  fprintf(outfile,"Problem finding displaced energy.\n");
  throw PsiException("FINDIF: Problem finding displaced energy.",__FILE__,__LINE__);
}

}}

