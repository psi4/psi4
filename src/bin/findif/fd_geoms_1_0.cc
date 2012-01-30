/*! \file
    \ingroup OPTKING
    \brief fd_geoms_1_0(): returns geometries necessary for finite-difference
     computation of gradients from energies; puts undisplaced geometry last in list
*/

#include "findif.h"

namespace psi { namespace findif {

std::vector< SharedMatrix > fd_geoms_1_0(Options &options) {

  fprintf(outfile,"\n-------------------------------------------------------------\n\n");

  fprintf(outfile,"  Using finite-differences of energies to determine gradients (fd_geoms_1_0).\n");

  int pts = options.get_int("POINTS");
  fprintf(outfile,"\tGenerating geometries for use with %d-point formula.\n",pts);
  if (pts != 3 && pts != 5)
    throw PsiException("FINDIF: Invalid number of points!",__FILE__,__LINE__);

  double disp_size = options.get_double("DISP_SIZE");
  fprintf(outfile,"\tDisplacement size will be %6.2e.\n", disp_size);

  // read in molecular data: Natom, reference geometry, and SALC coordinates
  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();

  int Natom = mol->natom();
  fprintf(outfile,"\tNumber of atoms is %d.\n", Natom);

  // Get SALCS from libmints
  boost::shared_ptr<MatrixFactory> fact;
  CdSalcList cdsalc(mol, fact, 0x1, true, true);

  int Nsalc = cdsalc.ncd();
  fprintf(outfile,"\tNumber of symmetric SALC's is %d.\n", Nsalc);

  // Determine number of geometries (1 + # of displacements)
  int Ndisp = 1;
  if (pts == 3)
    Ndisp += 2 * Nsalc;
  else if (pts == 5)
    Ndisp += 4 * Nsalc;

  fprintf(outfile, "\tNumber of displacements (including reference) is %d.\n", Ndisp);

  if (options.get_int("PRINT") > 1)
    for (int i=0; i<cdsalc.ncd(); ++i)
      cdsalc[i].print();

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
      displace_cart(geom_m, cdsalc, i, -1, disp_size);
      disp_geoms.push_back(geom_m);

      // + displacement
      SharedMatrix geom_p(ref_geom->clone());
      geom_p->set_name("Displacement + SALC #" + to_string(i+1));
      displace_cart(geom_p, cdsalc, i, +1, disp_size);
      disp_geoms.push_back(geom_p);

    }
  } // pts 3
  else if (pts == 5) {
    for (int i=0; i<Nsalc; ++i) {

      SharedMatrix geom_m2(ref_geom->clone());
      geom_m2->set_name("Displacement - SALC #" + to_string(i+1) + " * 2");
      displace_cart(geom_m2, cdsalc, i, -2, disp_size);
      disp_geoms.push_back(geom_m2);

      SharedMatrix geom_m1(ref_geom->clone());
      geom_m1->set_name("Displacement - SALC #" + to_string(i+1));
      displace_cart(geom_m1, cdsalc, i, -1, disp_size);
      disp_geoms.push_back(geom_m1);

      SharedMatrix geom_p1(ref_geom->clone());
      geom_p1->set_name("Displacement + SALC #" + to_string(i+1));
      displace_cart(geom_p1, cdsalc, i, +1, disp_size);
      disp_geoms.push_back(geom_p1);

      SharedMatrix geom_p2(ref_geom->clone());
      geom_p2->set_name("Displacement + SALC #" + to_string(i+1) + " * 2");
      displace_cart(geom_p2, cdsalc, i, +2, disp_size);
      disp_geoms.push_back(geom_p2);

    }
  } // pts 3

  // put reference geometry list in list
  disp_geoms.push_back(ref_geom);

  fprintf(outfile,"\n-------------------------------------------------------------\n");

  return disp_geoms;
}

}}

