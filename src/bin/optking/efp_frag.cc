/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*!
   \ingroup optking
*/

#include "efp_frag.h"
#include <libmints/matrix.h>
#include <libmints/vector3.h>

#define EXTERN
#include "globals.h"

#if defined (OPTKING_PACKAGE_QCHEM)
#include "EFP.h" // QChem's header for calling displace function
//****AVC****//
#endif
//****AVC****//


namespace opt {

void EFP_FRAG::set_values(double * values_in) {
  for (int i=0; i<6; ++i)
    values[i] = values_in[i];
}

void EFP_FRAG::set_forces(double * forces_in) {
  for (int i=0; i<6; ++i)
    forces[i] = forces_in[i];
}

void EFP_FRAG::set_xyz(double ** xyz_in) {
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      xyz_geom[i][j] = xyz_in[i][j];
}

void EFP_FRAG::set_com(double *com_in) {
  for(int i=0; i<3; i++)
    com[i] = com_in[i];
}

// we don't have a valid B matrix for these
// add 6 bogus stretches
void EFP_FRAG::add_dummy_intcos(int ndummy) {
  STRE *one_stre;
  for (int i=0; i<ndummy; ++i) {
    one_stre = new STRE(1, 2);
    intcos.push_back(one_stre);
  }
}

// We will assign a value of 0 to these on the first iteration; subsequently, the
// values will be calculated as total Delta(q) from the start of the optimization
void EFP_FRAG::print_intcos(FILE *fp) {
  double *v = get_values_pointer();
  fprintf(fp,"\t * Coordinate *           * BOHR/RAD *       * ANG/DEG *\n");
  fprintf(fp,"\t     COM X        %20.10lf%20.10lf \n", v[0], v[0] * _bohr2angstroms);
  fprintf(fp,"\t     COM Y        %20.10lf%20.10lf \n", v[1], v[1] * _bohr2angstroms);
  fprintf(fp,"\t     COM Z        %20.10lf%20.10lf \n", v[2], v[2] * _bohr2angstroms);
  fprintf(fp,"\t     alpha        %20.10lf%20.10lf \n", v[3], v[3] / _pi * 180.0);
  fprintf(fp,"\t     beta         %20.10lf%20.10lf \n", v[4], v[4] / _pi * 180.0);
  fprintf(fp,"\t     gamma        %20.10lf%20.10lf \n", v[5], v[5] / _pi * 180.0);
  fprintf(fp, "\n");
}

double **EFP_FRAG::H_guess(void) {
  double **H = init_matrix(6, 6);

  for (int xyz=0; xyz<3; xyz++)
  {
    H[xyz][xyz] = 0.01;
    H[xyz+3][xyz+3] = 0.03;
  }

  return H;
}

//****AVC****//
  // Tell QChem to update rotation matrix and com for EFP fragment
void EFP_FRAG::displace (int efp_frag_index, double *dq)
{
fprintf(outfile, "\nEFP_FRAG::displace()\n");
  using namespace psi;

  Vector3 T(dq[0],dq[1],dq[2]);
  Vector3 Phi(dq[3],dq[4],dq[5]);
  double phi = Phi.norm();

  SharedMatrix XYZ(new Matrix("XYZ", 3, 3));
  SharedMatrix R_XYZ(new Matrix("R_XYZ", 3, 3));
  SharedMatrix TR_XYZ(new Matrix("TR_XYZ", 3, 3));
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      XYZ->set(i,j,xyz_geom[i][j] - com[j]);

  XYZ->print_out();
  R_XYZ = XYZ->matrix_3d_rotation(Phi, phi, 0);
  R_XYZ->print_out();
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      TR_XYZ->set(i,j, R_XYZ->get(i,j) + T.get(j) + com[j]);
  TR_XYZ->print_out();

  double *xyz_array = init_array(3*3);

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      xyz_array[3*i+j] = TR_XYZ->get(i,j);

fprintf(outfile, "\nxyz_array\n");
  print_array(outfile, xyz_array, 3*3);
  p_efp->set_frag_coordinates(efp_frag_index, 1, xyz_array);
  free_array(xyz_array);

fprintf(outfile, "\n Phi            T\n");
for(int i=0; i<3; i++)
  fprintf(outfile, "%-15.8f%-15.8f\n", Phi.get(i), T.get(i));

}
//****AVC****//

}

//****AVC****//
//#endif
//****AVC****//

