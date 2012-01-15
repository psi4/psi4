/*! \file
    \ingroup OPTKING
    \brief fd_geoms_hessian_0(): returns geometries necessary for finite-difference
     computation of hessian; formulas in fd_geoms_freq_0 

     For starters will work only in C1 and without translations/rotations projected out
*/

#include "findif.h"

namespace psi { namespace findif {

std::vector< SharedMatrix > fd_geoms_hessian_0(Options &options) {

  fprintf(outfile,"\n-------------------------------------------------------------\n\n");

  fprintf(outfile,"  Using finite-differences of energies to determine hessian.\n");

  int pts = options.get_int("POINTS");
  fprintf(outfile,"\tGenerating geometries for use with %d-point formula.\n", pts);
  if (pts != 3 && pts != 5)
    throw PsiException("FINDIF: Invalid number of points!",__FILE__,__LINE__);

  double disp_size = options.get_double("DISP_SIZE");
  fprintf(outfile,"\tDisplacement size will be %6.2e.\n", disp_size);

  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();

  // make all salcs for now, just in the case the symmetric ones don't come out identically
  // we'll try to restrict later
  boost::shared_ptr<MatrixFactory> fact;
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

  fprintf(outfile,"\tNumber of symmetric SALC's is %d\n", symm_salcs.size());

  int Ndisp=1; // for reference geometry

  // diagonal displacements for symmetric coordinates
  if (pts == 3)
    Ndisp += 2 * symm_salcs.size();
  else if (pts == 5)
    Ndisp += 4 * symm_salcs.size();

  // off-diagonal displacements
  if (pts == 3)
    Ndisp += 2 * symm_salcs.size() * (symm_salcs.size() - 1) / 2;
  else if (pts == 5)
    Ndisp += 8 * symm_salcs.size() * (symm_salcs.size() - 1) / 2;

  fprintf(outfile,"\tNumber of symmetric displacements (including reference) is %d.\n", Ndisp);

  if (options.get_int("PRINT") > 1)
    for (int i=0; i<salc_list.ncd(); ++i)
      salc_list[i].print();

  // Get reference geometry
  Matrix ref_geom_temp = mol->geometry();
  SharedMatrix ref_geom(ref_geom_temp.clone());
  ref_geom->set_name("Reference geometry");

  // to be returned and converted into "matrix_vector" list in python
  std::vector< SharedMatrix > disp_geoms;

  for (int i=0; i<symm_salcs.size(); ++i) { // loop over salcs of this irrep
    int salc_i = symm_salcs[i];   // index in cdsalc of this salc

    if (pts == 3) {
      SharedMatrix geom1(ref_geom->clone());
      displace_cart(geom1, salc_list, salc_i, -1, disp_size);
      disp_geoms.push_back(geom1);

      SharedMatrix geom2(ref_geom->clone());
      displace_cart(geom2, salc_list, salc_i, +1, disp_size);
      disp_geoms.push_back(geom2);
    }
    else if (pts == 5) {
      SharedMatrix geom1(ref_geom->clone());
      displace_cart(geom1, salc_list, salc_i, -2, disp_size);
      disp_geoms.push_back(geom1);

      SharedMatrix geom2(ref_geom->clone());
      displace_cart(geom2, salc_list, salc_i, -1, disp_size);
      disp_geoms.push_back(geom2);

      SharedMatrix geom3(ref_geom->clone());
      displace_cart(geom3, salc_list, salc_i, +1, disp_size);
      disp_geoms.push_back(geom3);

      SharedMatrix geom4(ref_geom->clone());
      displace_cart(geom4, salc_list, salc_i, +2, disp_size);
      disp_geoms.push_back(geom4);
    }
  } // i, salcs of this irrep

  // off-diagonal displacements
  for (int i=0; i<symm_salcs.size(); ++i) { // loop over salcs of this irrep
    int salc_i = symm_salcs[i];     // index in cdsalc of this salc

    for (int j=0; j<i; ++j) {        // loop over salcs of this irrep
    int salc_j = symm_salcs[j];   // index in cdsalc of this salc

      if (pts == 3) {
        SharedMatrix geom1(ref_geom->clone());
        displace_cart(geom1, salc_list, salc_i, salc_j, +1, +1, disp_size);
        disp_geoms.push_back(geom1);

        SharedMatrix geom2(ref_geom->clone());
        displace_cart(geom2, salc_list, salc_i, salc_j, -1, -1, disp_size);
        disp_geoms.push_back(geom2);
      }
      else if (pts == 5) {
        SharedMatrix geom1(ref_geom->clone());
        displace_cart(geom1, salc_list, salc_i, salc_j, -1, -2, disp_size);
        disp_geoms.push_back(geom1);

        SharedMatrix geom2(ref_geom->clone());
        displace_cart(geom2, salc_list, salc_i, salc_j, -2, -1, disp_size);
        disp_geoms.push_back(geom2);

        SharedMatrix geom3(ref_geom->clone());
        displace_cart(geom3, salc_list, salc_i, salc_j, -1, -1, disp_size);
        disp_geoms.push_back(geom3);

        SharedMatrix geom4(ref_geom->clone());
        displace_cart(geom4, salc_list, salc_i, salc_j, +1, -1, disp_size);
        disp_geoms.push_back(geom4);

        SharedMatrix geom5(ref_geom->clone());
        displace_cart(geom5, salc_list, salc_i, salc_j, -1, +1, disp_size);
        disp_geoms.push_back(geom5);

        SharedMatrix geom6(ref_geom->clone());
        displace_cart(geom6, salc_list, salc_i, salc_j, +1, +1, disp_size);
        disp_geoms.push_back(geom6);

        SharedMatrix geom7(ref_geom->clone());
        displace_cart(geom7, salc_list, salc_i, salc_j, +2, +1, disp_size);
        disp_geoms.push_back(geom7);

        SharedMatrix geom8(ref_geom->clone());
        displace_cart(geom8, salc_list, salc_i, salc_j, +1, +2, disp_size);
        disp_geoms.push_back(geom8);
      } // pts == 5
    } // m, salc_j
  } // i, salc_i 

  // put reference geometry list in list
  disp_geoms.push_back(ref_geom);

  fprintf(outfile,"\n-------------------------------------------------------------\n");

  return disp_geoms;
}

}}

