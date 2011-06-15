/*! \file
    \ingroup OPTKING
    \brief fd_grad_1_0(): compute gradient using energies and finite-differences
*/

#include <libmints/mints.h>
#include <libqt/qt.h>
#include <psi4-dec.h>

namespace psi { namespace findif {

PsiReturnType fd_grad_1_0(Options &options, boost::shared_ptr<Vector> E) {

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

  if (E->dim() != Ndisp)
    throw PsiException("FINDIF: Incorrect number of energies passed in!",__FILE__,__LINE__);

  // Compute gradient in mass-weighted symmetry-adapted cartesians in ATOMIC units
  double *g_q = init_array(Nsalc);
  if (pts == 3) {
    for (int i=0; i<Nsalc; ++i)
      g_q[i] = (E->get(0,2*i+1) - E->get(0,2*i)) / (2.0 * Dx);
  }
  else if (pts == 5) {
    for (int i=0; i<Nsalc; ++i)
      g_q[i] = (E->get(0,4*i) - 8.0*E->get(0,4*i+1) + 8.0*E->get(0,4*i+2)
              - E->get(0,4*i+3)) / (12.0 * Dx);
  }

  // Print out energies and gradients
  double energy_ref = E->get(0, Ndisp-1);
  fprintf(outfile, "\tFinite difference computation of gradient using %d-point formula\n", pts);
  fprintf(outfile, "\tCheck for precision!\n");
  fprintf(outfile, "\tEnergy without displacment: %15.10lf\n", energy_ref);
  fprintf(outfile, "\tGradients are in mass-weighted, symmetry-adapted cartesians (in au).\n");

  int cnt;
  if (pts == 3) {
    cnt = -2;
    fprintf(outfile," Coord      Energy(-)        Energy(+)        Force\n");
    for (int i=0; i<Nsalc; ++i) {
      cnt += 2;
      fprintf(outfile,"%5d %17.10lf%17.10lf%17.10lf\n", i, E->get(0,cnt), E->get(0,cnt+1), g_q[i]);
    }
    fprintf(outfile,"\n");
  }
  else if (pts == 5) {
    cnt = -4;
    fprintf(outfile,
      " Coord      Energy(-2)        Energy(-1)        Energy(+1)        Energy(+2)            Force\n");
    for (int i=0; i<Nsalc; ++i) {
      cnt += 4;
      fprintf(outfile,"%5d %17.10lf %17.10lf %17.10lf %17.10lf %17.10lf\n", i, E->get(0,cnt),
        E->get(0,cnt+1), E->get(0,cnt+2), E->get(0,cnt+3), g_q[i]);
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

  // compute gradient in ordinary cartesians
  double *g_cart = init_array(3*Natom);

  // B^t g_q^t = g_x^t -> g_q B = g_x
  C_DGEMM('n', 'n', 1, 3*Natom, Nsalc, 1.0, g_q, Nsalc, B[0], 3*Natom, 0, g_cart, 3*Natom);

  free(g_q);
  free_block(B);

  // The undisplaced geometry should be in the global molecule, and the undisplaced
  // energy in globals["CURRENT ENERGY"], since we did that one last.  Clever, huh.

  // Write out the geometry and gradient to file 11
  SimpleMatrix gradient_matrix("F-D gradient", Natom, 3);

  for (int a=0; a<Natom; ++a)
    for (int xyz=0; xyz<3; ++xyz)
      gradient_matrix.set(a, xyz, g_cart[3*a+xyz]);

  GradientWriter grad(mol, gradient_matrix);
  grad.write("psi.file11.dat");

  free(g_cart);

  return Success;
}

}}

