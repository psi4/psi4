/*! \file
    \ingroup OPTKING
    \brief fd_freq_0(): compute frequencies from energies
*/

#include "findif.h"

#include <boost/python.hpp>
#include <boost/python/list.hpp>
using namespace boost::python;

#include <physconst.h>

namespace psi { namespace findif {

int iE0(std::vector<int> & Ndisp_pi, std::vector< std::vector<int> > & salcs_pi, int pts,
 int irrep, int ii, int jj, int disp_i, int disp_j);

PsiReturnType fd_freq_0(Options &options, const boost::python::list& E_list, int freq_irrep_only)
{
  int pts = options.get_int("POINTS");
  double disp_size = options.get_double("DISP_SIZE");
  int print_lvl = options.get_int("PRINT");

  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();
  int Natom = mol->natom();
  boost::shared_ptr<MatrixFactory> fact;
  // symmetric modes with rotations and translations projected out
  CdSalcList salc_list(mol, fact);
  int Nirrep = salc_list.nirrep();
  int Nsalc_all = salc_list.ncd();

  // build vectors that list indices of salcs for each irrep
  std::vector< std::vector<int> > salcs_pi; // salcs per irrep
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
  for (int i=0; i<len(E_list); ++i)
    E.push_back( (double)extract<double>(E_list[i]) );

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

  char **irrep_lbls = mol->irrep_labels();

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

    // diagonalize force constant matrix
    int dim = salcs_pi[h].size(); 
    double *evals= init_array(dim);
    double **evects = block_matrix(dim, dim);

    sq_rsp(dim, dim, H_irr, evals, 3, evects, 1e-14);

    // Build Bu^1/2 matrix for this irrep
    double **B_irr = block_matrix(salcs_pi[h].size(), 3*Natom);

    for (int i=0; i<salcs_pi[h].size(); ++i) {
      int salc_i = salcs_pi[h][i];
      for (int c=0; c<salc_list[salc_i].ncomponent(); ++c) {
        int a          = salc_list[salc_i].component(c).atom;
        int xyz        = salc_list[salc_i].component(c).xyz;
        double coef    = salc_list[salc_i].component(c).coef;
        B_irr[i][3*a+xyz] = coef / sqrt(mol->mass(a));
      }
    }

    double **normal_irr = block_matrix(3*Natom, dim);
    C_DGEMM('t', 'n', 3*Natom, dim, dim, 1.0, B_irr[0], 3*Natom, evects[0],
      dim, 0, normal_irr[0], dim);

    if (print_lvl >= 3) {
      fprintf(outfile,"\n\tNormal coordinates (mass-weighted) for irrep %s:\n", irrep_lbls[h]);
      eivout(normal_irr, evals, 3*Natom, dim, outfile);
    }

    for (int i=0; i<salcs_pi[h].size(); ++i) {
      // copy eigenvector into row
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

