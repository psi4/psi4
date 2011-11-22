/*! \file
    \ingroup OPTKING
    \brief fd_geoms_1_0(): returns geometries necessary for finite-difference
     computation of gradients from energies; puts undisplaced geometry last in list
*/

#include "findif.h"

namespace psi { namespace findif {

std::vector< SharedMatrix > fd_geoms_1_0(Options &options) {

  // Print what we are doing
  fprintf(outfile,"\tUsing finite-differences of energies to determine gradients.\n");

  int pts = options.get_int("POINTS");
  fprintf(outfile,"\tGenerating geometries for use with %d-point formula.\n",pts);
  if (pts != 3 && pts != 5)
    throw PsiException("FINDIF: Invalid number of points!",__FILE__,__LINE__);

  double Dx = options.get_double("DISP_SIZE");
  fprintf(outfile,"\tDisplacement size will be %6.2e.\n", Dx);

  // read in molecular data: Natom, reference geometry, and SALC coordinates
  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();

  int Natom = mol->natom();
  fprintf(outfile,"\tNumber of atoms is %d.\n", Natom);

  // Get SALCS from libmints
  boost::shared_ptr<MatrixFactory> fact;
  CdSalcList cdsalc(mol, fact, 0x1, true, true);

  int Nsalc = cdsalc.ncd();
  fprintf(outfile,"\tNumber of SALC's is %d.\n", Nsalc);

  for (int i=0; i<cdsalc.ncd(); ++i)
    cdsalc[i].print();

  // Determine number of geometries (1 + # of displacements)
  int Ndisp = 1;
  if (pts == 3)
    Ndisp += 2 * Nsalc;
  else if (pts == 5)
    Ndisp += 4 * Nsalc;

  fprintf(outfile, "\tNumber of displacements (including reference) is %d.\n", Ndisp);

  // Get reference geometry
  Matrix ref_geom_temp = mol->geometry();
  SharedMatrix ref_geom(ref_geom_temp.clone());

  ref_geom->set_name("Reference geometry");

  // to be returned and converted into "matrix_vector" list in python
  std::vector< SharedMatrix > disp_geoms;

  if (pts == 3) {
    for (int i=0; i<Nsalc; ++i) {
      // - displacement
      SharedMatrix geom_m(ref_geom->clone());
      geom_m->set_name("Displacement - SALC #" + to_string(i+1));

      // + displacement
      SharedMatrix geom_p(ref_geom->clone());
      geom_p->set_name("Displacement + SALC #" + to_string(i+1));

      int nc = cdsalc[i].ncomponent();
      //fprintf(outfile, "Ncomponent of salc %d is %d\n", i, nc);

      for (int c=0; c<nc; ++c) {
        int a          = cdsalc[i].component(c).atom;
        int xyz        = cdsalc[i].component(c).xyz;
        double coef    = cdsalc[i].component(c).coef;

        geom_m->add(0, a, xyz, - Dx * coef / sqrt(mol->mass(a)));
        geom_p->add(0, a, xyz, + Dx * coef / sqrt(mol->mass(a)));
      }

      disp_geoms.push_back(geom_m);
      disp_geoms.push_back(geom_p);
    }
  } // pts 3
  else if (pts == 5) {
    for (int i=0; i<Nsalc; ++i) {
      SharedMatrix geom_m2(ref_geom->clone());
      geom_m2->set_name("Displacement - SALC #" + to_string(i+1) + " * 2");

      SharedMatrix geom_m1(ref_geom->clone());
      geom_m1->set_name("Displacement - SALC #" + to_string(i+1));

      SharedMatrix geom_p1(ref_geom->clone());
      geom_p1->set_name("Displacement + SALC #" + to_string(i+1));

      SharedMatrix geom_p2(ref_geom->clone());
      geom_p2->set_name("Displacement + SALC #" + to_string(i+1) + " * 2");

      int nc = cdsalc[i].ncomponent();

      for (int c=0; c<nc; ++c) {
        int a          = cdsalc[i].component(c).atom;
        int xyz        = cdsalc[i].component(c).xyz;
        double coef    = cdsalc[i].component(c).coef;

        geom_m2->add(0, a, xyz, - 2.0 * Dx * coef / sqrt(mol->mass(a)));
        geom_m1->add(0, a, xyz, - 1.0 * Dx * coef / sqrt(mol->mass(a)));
        geom_p1->add(0, a, xyz, + 1.0 * Dx * coef / sqrt(mol->mass(a)));
        geom_p2->add(0, a, xyz, + 2.0 * Dx * coef / sqrt(mol->mass(a)));
      }

      disp_geoms.push_back(geom_m2);
      disp_geoms.push_back(geom_m1);
      disp_geoms.push_back(geom_p1);
      disp_geoms.push_back(geom_p2);
    }
  } // pts 3

  // put reference geometry list in list
  disp_geoms.push_back(ref_geom);

  return disp_geoms;
}

}}

