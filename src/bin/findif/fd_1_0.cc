/*! \file
    \ingroup OPTKING
    \brief fd_1_0(): compute gradient using energies and finite-differences
*/

#include "findif.h"

#include <boost/python.hpp>
#include <boost/python/list.hpp>
using namespace boost::python;

namespace psi { namespace findif {

PsiReturnType
fd_1_0(Options &options, const boost::python::list& E)
{
  int pts = options.get_int("POINTS");
  double Dx = options.get_double("DISP_SIZE");

  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();
  int Natom = mol->natom();
  boost::shared_ptr<MatrixFactory> fact;

  CdSalcList cdsalc(mol, fact, 0x1, true, true);
  int Nsalc = cdsalc.ncd();

  // Compute number of displacements - check with number of energies passed in
  // Determine number of geometries (1 + # of displacements)
  int Ndisp = 1;
  if (pts == 3)
    Ndisp += 2 * Nsalc;
  else if (pts == 5)
    Ndisp += 4 * Nsalc;
  else
      throw PSIEXCEPTION("fd_1_0: Unable to handle requested point formula. 3 or 5-point formula are supported.");

  if (len(E) != Ndisp)
    throw PsiException("FINDIF: Incorrect number of energies passed in!",__FILE__,__LINE__);

  // Compute gradient in mass-weighted symmetry-adapted cartesians in ATOMIC units
  double *g_q = init_array(Nsalc);
  if (pts == 3) {
    for (int i=0; i<Nsalc; ++i)
      g_q[i] = (extract<double>(E[2*i+1]) - extract<double>(E[2*i])) / (2.0 * Dx);
  }
  else if (pts == 5) {
    for (int i=0; i<Nsalc; ++i)
      g_q[i] = (extract<double>(E[4*i])
                - 8.0*extract<double>(E[4*i+1])
                + 8.0*extract<double>(E[4*i+2])
                - extract<double>(E[4*i+3])) / (12.0 * Dx);
  }

  fprintf(outfile,"\n-------------------------------------------------------------\n\n");
  fprintf(outfile,"  Computing gradient from energies (fd_1_0).\n");

  // Print out energies and gradients
  double energy_ref = extract<double>(E[Ndisp-1]);
  fprintf(outfile, "\tUsing %d-point formula.\n", pts);
  fprintf(outfile, "\tEnergy without displacement: %15.10lf\n", energy_ref);
  fprintf(outfile, "\tCheck energies below for precision!\n");
  fprintf(outfile, "\tForces are for mass-weighted, symmetry-adapted cartesians (in au).\n");

  int cnt;
  if (pts == 3) {
    cnt = -2;
    fprintf(outfile,"\n\t Coord      Energy(-)        Energy(+)        Force\n");
    for (int i=0; i<Nsalc; ++i) {
      cnt += 2;
      fprintf(outfile,"\t%5d %17.10lf%17.10lf%17.10lf\n",
              i,
              (double)extract<double>(E[cnt]),
              (double)extract<double>(E[cnt+1]),
              g_q[i]);
    }
    fprintf(outfile,"\n");
  }
  else if (pts == 5) {
    cnt = -4;
    fprintf(outfile,
      "\n\t Coord      Energy(-2)        Energy(-1)        Energy(+1)        Energy(+2)            Force\n");
    for (int i=0; i<Nsalc; ++i) {
      cnt += 4;
      fprintf(outfile,"\t%5d %17.10lf %17.10lf %17.10lf %17.10lf %17.10lf\n",
              i,
              (double)extract<double>(E[cnt]),
              (double)extract<double>(E[cnt+1]),
              (double)extract<double>(E[cnt+2]),
              (double)extract<double>(E[cnt+3]),
              g_q[i]);
    }
    fprintf(outfile,"\n");
  }

  // Build B matrix of salc coefficients
  double **B = block_matrix(Nsalc, 3*Natom);

  for (int i=0; i<Nsalc; ++i) {
    int nc = cdsalc[i].ncomponent();
    for (int c=0; c<nc; ++c) {
      int a          = cdsalc[i].component(c).atom;
      int xyz        = cdsalc[i].component(c).xyz;
      double coef    = cdsalc[i].component(c).coef;
      B[i][3*a+xyz] = coef;
    }
  }

  // compute gradient in mass-weighted (non-SALC) cartesians
  double *g_cart = init_array(3*Natom);

  // B^t g_q^t = g_x^t -> g_q B = g_x
  C_DGEMM('n', 'n', 1, 3*Natom, Nsalc, 1.0, g_q, Nsalc, B[0], 3*Natom, 0, g_cart, 3*Natom);

  free(g_q);
  free_block(B);

  // The undisplaced geometry should be in the global molecule, and the undisplaced
  // energy in globals["CURRENT ENERGY"], since we did that one last.  Clever, huh.

  // Write out the geometry and gradient to file 11
  Matrix gradient_matrix("F-D gradient", Natom, 3);

  // Convert gradient to un-massweighted cartesians
  for (int a=0; a<Natom; ++a)
    for (int xyz=0; xyz<3; ++xyz)
      gradient_matrix.set(a, xyz, g_cart[3*a+xyz] * sqrt(mol->mass(a)));

  GradientWriter grad(mol, gradient_matrix);
  grad.write("psi.file11.dat");
  fprintf(outfile,"\tGradient written to file11.\n");

  SharedMatrix sgradient(gradient_matrix.clone());
  Process::environment.reference_wavefunction()->set_gradient(sgradient);
  fprintf(outfile,"\tGradent saved to wavefunction.\n");

  fprintf(outfile,"\n-------------------------------------------------------------\n");

  free(g_cart);

  return Success;
}

}}

