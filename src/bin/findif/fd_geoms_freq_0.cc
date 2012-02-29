/*! \file
    \ingroup OPTKING
    \brief fd_geoms_freq_0(): returns geometries necessary for finite-difference
     computation of second derivatives from energies; puts undisplaced geometry last in list

  **Finite-difference formulas
  3-point - diagonal
    O(1/h^2): [f(1,0) + f(-1,0) - 2f(0,0)]/(h^2)
     [which is same as : [f(2,0) + f(-2,0) - 2f(0,0)]/(4h^2)
  3-point - off-diagonal
    O(1/h^2): new-way: [f(1,1)+f(-1,-1)+2f(0,0) -f(1,0) -f(-1,0) -f(0,1) -f(0,-1)]/(2h^2)

  5-point formula - diagonal
    O(1/h^4): [-f(-2,0) + 16f(-1,0) + 16f(1,0) - f(2,0) - 30f(0,0)] / (12h^2)
  5-point formula - off-diagonal
    O(1/h^4): [-1f(-1,-2) - 1f(-2,-1) + 9f(-1,-1) - 1f(+1,-1)
      - 1f(-1,1) + 9f(+1,+1) - 1f(+2,+1) - 1f(1,2)
      + 1f(-2,0) - 7f(-1,0)  - 7f(+1,0) + 1f(+2,0)
      + 1f(0,-2) - 7f(0,-1)  - 7f(0,+1) + 1f(0,+2) + 12f(0,0)]/(12h^2)
*/

#include "findif.h"

namespace psi { namespace findif {

// freq_irrep_only designates single irrep for which to compute frequencies

std::vector< SharedMatrix > fd_geoms_freq_0(Options &options, int freq_irrep_only) {

  fprintf(outfile,"\n-------------------------------------------------------------\n\n");

  fprintf(outfile,"  Using finite-differences of energies to determine vibrational frequencies and \n");
  fprintf(outfile,"  normal modes.  Resulting frequencies are only valid at stationary points.\n");

  int pts = options.get_int("POINTS");
  fprintf(outfile,"\tGenerating geometries for use with %d-point formula.\n", pts);
  if (pts != 3 && pts != 5)
    throw PsiException("FINDIF: Invalid number of points!",__FILE__,__LINE__);

  double disp_size = options.get_double("DISP_SIZE");
  fprintf(outfile,"\tDisplacement size will be %6.2e.\n", disp_size);

  const boost::shared_ptr<Molecule> mol = psi::Process::environment.molecule();

  // Get SALCS from libmints: all modes with rotations and translations projected out
  boost::shared_ptr<MatrixFactory> fact;
  CdSalcList salc_list(mol, fact, 0xF, true, true);

  int Natom = mol->natom();
  fprintf(outfile,"\tNumber of atoms is %d.\n", Natom);

  int Nirrep = salc_list.nirrep();
  fprintf(outfile,"\tNumber of irreps is %d.\n", Nirrep);

  int Nsalc_all = salc_list.ncd();
  fprintf(outfile,"\tNumber of SALCS is %d.\n", Nsalc_all);

  if (freq_irrep_only >= Nirrep || freq_irrep_only < -1)
    throw PsiException("FINDIF: Irrep value not in valid range.",__FILE__,__LINE__);

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

  // From now on in code, salcs_pi establishes the canonical order of SALCs and displacements
  fprintf(outfile,"\t Number of SALC's per irrep:\n");
  for (int h=0; h<Nirrep; ++h)
    fprintf(outfile,"\t  Irrep %d: %d\n", h+1, (int) salcs_pi[h].size());

  if (options.get_int("PRINT") > 1) {
    fprintf(outfile,"\tIndex of active salcs per irrep:\n");
    for (int h=0; h<Nirrep; ++h) {
      if (salcs_pi[h].size()) {
        fprintf(outfile, "\t  Irrep %d : ", h+1);
        for (int i=0; i<salcs_pi[h].size(); ++i)
          fprintf(outfile," %d ", salcs_pi[h][i]);
        fprintf(outfile, "\n");
      }
    }
  }

  // Determine number of displacements
  std::vector<int> Ndisp_pi (Nirrep);

  // diagonal displacements for symmetric coordinates
  if (pts == 3)
    Ndisp_pi[0] = 2 * salcs_pi[0].size();
  else if (pts == 5)
    Ndisp_pi[0] = 4 * salcs_pi[0].size();

  // diagonal displacements for asymmetric coordinates
  for (int h=1; h<Nirrep; ++h) {
    if (pts == 3)
      Ndisp_pi[h] = salcs_pi[h].size();
    else if (pts == 5)
      Ndisp_pi[h] = 2* salcs_pi[h].size();
  }

  // off-diagonal displacements
  for (int h=0; h<Nirrep; ++h) {
    if (pts == 3)
      Ndisp_pi[h] += 2 * salcs_pi[h].size() * (salcs_pi[h].size() - 1) / 2;
    else if (pts == 5)
      Ndisp_pi[h] += 8 * salcs_pi[h].size() * (salcs_pi[h].size() - 1) / 2;
  }

  int Ndisp_all = 0;
  for (int h=0; h<Nirrep; ++h)
    Ndisp_all += Ndisp_pi[h];

  fprintf(outfile,"\tNumber of geometries (including reference) is %d.\n", Ndisp_all);
  fprintf(outfile,"\tNumber of displacements per irrep:\n");
  for (int h=0; h<Nirrep; ++h)
    fprintf(outfile,"\t  Irrep %d: %d\n", h+1, Ndisp_pi[h]);

  if (options.get_int("PRINT") > 2)
    for (int i=0; i<salc_list.ncd(); ++i)
      salc_list[i].print();

  // Get reference geometry
  Matrix ref_geom_temp = mol->geometry();
  SharedMatrix ref_geom(ref_geom_temp.clone());
  ref_geom->set_name("Reference geometry");

  // disp_geoms will be returned and converted into "matrix_vector" list in python
  std::vector< SharedMatrix > disp_geoms;

  for (int h=0; h<Nirrep; ++h) { // loop over irreps

    for (int i=0; i<salcs_pi[h].size(); ++i) { // loop over salcs of this irrep
      int salc_i = salcs_pi[h][i];   // index in cdsalc of this salc

      if (h == 0) { // symmetric diagonal displacements (first in list)
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
      }
      else { // h != 0; assymmetric diagonal displacements
        if (pts == 3) {
          SharedMatrix geom1(ref_geom->clone());
          displace_cart(geom1, salc_list, salc_i, -1, disp_size);
          disp_geoms.push_back(geom1);
        }
        else if (pts == 5) {
          SharedMatrix geom1(ref_geom->clone());
          displace_cart(geom1, salc_list, salc_i, -2, disp_size);
          disp_geoms.push_back(geom1);

          SharedMatrix geom2(ref_geom->clone());
          displace_cart(geom2, salc_list, salc_i, -1, disp_size);
          disp_geoms.push_back(geom2);
        }
      }
    } // i, salcs of this irrep

    // off-diagonal displacements
    for (int i=0; i<salcs_pi[h].size(); ++i) { // loop over salcs of this irrep
      int salc_i = salcs_pi[h][i];     // index in cdsalc of this salc

      for (int j=0; j<i; ++j) {        // loop over salcs of this irrep
        int salc_j = salcs_pi[h][j];   // index in cdsalc of this salc

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
  } // h, irreps

  // put reference geometry list in list
  disp_geoms.push_back(ref_geom);

  fprintf(outfile,"\n-------------------------------------------------------------\n");

  return disp_geoms;
}

}}

