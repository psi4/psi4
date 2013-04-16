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
  double **H = init_matrix(g_nintco(), g_nintco());

  for (int i=0; i<g_nintco(); ++i)
    H[i][i] = 0.01;

  return H;
}

//****AVC****//
#if defined (OPTKING_PACKAGE_QCHEM)
//****AVC****//
  // Tell QChem to update rotation matrix and com for EFP fragment
void EFP_FRAG::displace (int efp_frag_index, double *dq) {
  ::EFP::GetInstance()->displace(efp_frag_index, dq);
}
//****AVC****//
#endif
//****AVC****//

}

//****AVC****//
//#endif
//****AVC****//

