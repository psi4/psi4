/*! \file
    \ingroup OPTKING
    \brief fd_geoms_2_0(): returns geometries necessary for finite-difference
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

void displace_simple_cart(boost::shared_ptr<Matrix> geom, int coord, int a, double disp_size);
void displace_simple_cart(boost::shared_ptr<Matrix> geom, int coord, int a, 
  int coord_2, int b, double disp_size);

std::vector< boost::shared_ptr<Matrix> > fd_geoms_2_0(Options &options) {

  fprintf(outfile,"\tUsing finite-differences of energies to determine second derivatives.\n");

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

  // count number of displacements
  int Ndisp;

  // diagonal
  if (pts == 3)
    Ndisp = 2 * (3*Natom);
  else if (pts == 5)
    Ndisp = 4 * (3*Natom);

  // off-diagonal 
  if (pts == 3)
    Ndisp += 2 * (3*Natom) * (3*Natom - 1) / 2;
  else if (pts == 5)
    Ndisp += 8 * (3*Natom) * (3*Natom - 1) / 2;

  fprintf(outfile,"\tNumber of displacements is: %d (plus reference).\n", Ndisp);

  // Get reference geometry
  Matrix ref_geom_temp = mol->geometry();
  SharedMatrix ref_geom(ref_geom_temp.clone());
  ref_geom->set_name("Reference geometry");

  // to be returned and converted into "matrix_vector" list in python
  std::vector< boost::shared_ptr<Matrix> > disp_geoms;

  if (pts == 3) { // 3-point formula

    // do all diagonal displacements first
    for (int i=0; i<3*Natom; ++i) {
      boost::shared_ptr<Matrix> geom1(ref_geom->clone());
      displace_simple_cart(geom1, i, -1, disp_size);
      disp_geoms.push_back(geom1);

      boost::shared_ptr<Matrix> geom2(ref_geom->clone());
      displace_simple_cart(geom2, i, +1, disp_size);
      disp_geoms.push_back(geom2);
    }

    // do all off-diagonal displacements
    for (int i=0; i<3*Natom; ++i) {
      for (int j=0; j<i; ++j) {
        boost::shared_ptr<Matrix> geom3(ref_geom->clone());
        displace_simple_cart(geom3, i, j, +1, +1, disp_size);
        disp_geoms.push_back(geom3);

        boost::shared_ptr<Matrix> geom4(ref_geom->clone());
        displace_simple_cart(geom4, i, j, -1, -1, disp_size);
        disp_geoms.push_back(geom4);
      }
    }

  }
  else if (pts == 5) {

    // do all diagonal displacements first
    for (int i=0; i<3*Natom; ++i) {
      boost::shared_ptr<Matrix> geom1(ref_geom->clone());
      displace_simple_cart(geom1, i, -2, disp_size);
      disp_geoms.push_back(geom1);

      boost::shared_ptr<Matrix> geom2(ref_geom->clone());
      displace_simple_cart(geom2, i, -1, disp_size);
      disp_geoms.push_back(geom2);

      boost::shared_ptr<Matrix> geom3(ref_geom->clone());
      displace_simple_cart(geom3, i, +1, disp_size);
      disp_geoms.push_back(geom3);

      boost::shared_ptr<Matrix> geom4(ref_geom->clone());
      displace_simple_cart(geom4, i, +2, disp_size);
      disp_geoms.push_back(geom4);
    }

    // do all off-diagonal displacements
    for (int i=0; i<3*Natom; ++i) {
      for (int j=0; j<i; ++j) {
        boost::shared_ptr<Matrix> geom5(ref_geom->clone());
        displace_simple_cart( geom5, i, j, -1, -2, disp_size);
        disp_geoms.push_back(geom5);

        boost::shared_ptr<Matrix> geom6(ref_geom->clone());
        displace_simple_cart( geom6, i, j, -2, -1, disp_size);
        disp_geoms.push_back(geom6);

        boost::shared_ptr<Matrix> geom7(ref_geom->clone());
        displace_simple_cart( geom7, i, j, -1, -1, disp_size);
        disp_geoms.push_back(geom7);

        boost::shared_ptr<Matrix> geom8(ref_geom->clone());
        displace_simple_cart( geom8, i, j, +1, -1, disp_size);
        disp_geoms.push_back(geom8);

        boost::shared_ptr<Matrix> geom9(ref_geom->clone());
        displace_simple_cart( geom9, i, j, -1, +1, disp_size);
        disp_geoms.push_back(geom9);

        boost::shared_ptr<Matrix> geom10(ref_geom->clone());
        displace_simple_cart(geom10, i, j, +1, +1, disp_size);
        disp_geoms.push_back(geom10);

        boost::shared_ptr<Matrix> geom11(ref_geom->clone());
        displace_simple_cart(geom11, i, j, +2, +1, disp_size);
        disp_geoms.push_back(geom11);

        boost::shared_ptr<Matrix> geom12(ref_geom->clone());
        displace_simple_cart(geom12, i, j, +1, +2, disp_size);
        disp_geoms.push_back(geom12);
      }
    }

  }

  // put reference geometry list in list
  disp_geoms.push_back(ref_geom);

  return disp_geoms;
}

// displaces from a reference geometry: geom[ii] += disp_i * disp_size
// fixes the name
void displace_simple_cart(boost::shared_ptr<Matrix> geom, int ii, int disp_i, double Dx) {
  geom->set_name("Coord: " + to_string(ii) + ", Disp: " + to_string(disp_i));

  int atom = ii / 3;
  int xyz = ii % 3;

  geom->add(0, atom, xyz, disp_i * Dx);

  return;
}

/* displaces from a reference geometry:
  geom[ii] += disp_i * disp_size
  geom[jj] += disp_j * disp_size
*/
void displace_simple_cart(boost::shared_ptr<Matrix> geom, int ii, int jj, int disp_i, int disp_j,
  double Dx) {

  geom->set_name("Coord: " + to_string(ii+1) + ", Disp: " + to_string(disp_i) + " Coord: "
  + to_string(jj+1) + ", Disp: " + to_string(disp_j));

  int atom = ii / 3;
  int xyz = ii % 3;
  geom->add(0, atom, xyz, disp_i * Dx);

  atom = jj / 3;
  xyz = jj % 3;
  geom->add(0, atom, xyz, disp_j * Dx);

  return;
}

}}

