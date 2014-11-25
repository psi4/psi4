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

#include "fb_frag.h"

#include "print.h"
#define EXTERN
#include "globals.h"

#if defined (OPTKING_PACKAGE_QCHEM)
#include "EFP.h" // QChem's header for calling displace function
#endif

namespace opt {

void FB_FRAG::set_values(double * values_in) {
  for (int i=0; i<6; ++i)
    values[i] = values_in[i];
}

void FB_FRAG::set_forces(double * forces_in) {
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

double * EFP_FRAG::get_geom_array(void) {
    int cnt=0;
    double * geom_array = init_array(3*3);
    for (int i=0; i<3; ++i)
      for (int xyz=0; xyz<3; ++xyz)
        geom_array[cnt++] = xyz_geom[i][xyz];
    return geom_array;
}

void EFP_FRAG::set_geom_array( double * geom_array ) {
  int cnt = 0;
  for (int a=0; a<3; a++)
    for (int xyz=0; xyz<3; ++xyz)
      xyz_geom[a][xyz] = geom_array[cnt++];
}

// we don't have a valid B matrix for these
// add 6 bogus stretches
void FB_FRAG::add_dummy_coords(int ndummy) {
  STRE *one_stre;
  for (int i=0; i<ndummy; ++i) {
    one_stre = new STRE(1, 2);
    coords.simples.push_back(one_stre);
  }
}

// We will assign a value of 0 to these on the first iteration; subsequently, the
// values will be calculated as total Delta(q) from the start of the optimization
void FB_FRAG::print_intcos(std::string psi_fp, FILE *qc_fp) {
  double *v = get_values_pointer();
  oprintf(psi_fp, qc_fp,"\t * Coordinate *           * BOHR/RAD *       * ANG/DEG *\n");
  oprintf(psi_fp, qc_fp,"\t     COM X        %20.10lf%20.10lf \n", v[0], v[0] * _bohr2angstroms);
  oprintf(psi_fp, qc_fp,"\t     COM Y        %20.10lf%20.10lf \n", v[1], v[1] * _bohr2angstroms);
  oprintf(psi_fp, qc_fp,"\t     COM Z        %20.10lf%20.10lf \n", v[2], v[2] * _bohr2angstroms);
  oprintf(psi_fp, qc_fp,"\t     alpha        %20.10lf%20.10lf \n", v[3], v[3] / _pi * 180.0);
  oprintf(psi_fp, qc_fp,"\t     beta         %20.10lf%20.10lf \n", v[4], v[4] / _pi * 180.0);
  oprintf(psi_fp, qc_fp,"\t     gamma        %20.10lf%20.10lf \n", v[5], v[5] / _pi * 180.0);
  oprintf(psi_fp, qc_fp, "\n");
}

double **FB_FRAG::H_guess(void) {
  double **H = init_matrix(Ncoord(), Ncoord());

  for (int i=0; i<Ncoord(); ++i)
    H[i][i] = 0.01;

  return H;
}

  // Tell QChem to update rotation matrix and com for FB fragment
void FB_FRAG::displace (int fb_frag_index, double *dq) {
#if defined (OPTKING_PACKAGE_QCHEM)
  ::EFP::GetInstance()->displace(fb_frag_index, dq);
#endif
////****AVC****//
//  // Tell QChem to update rotation matrix and com for EFP fragment
//void EFP_FRAG::displace (int efp_frag_index, double *dq)
//{
//fprintf(outfile, "\nEFP_FRAG::displace()\n");
//  double *T = init_array(3);
//  double *Phi = init_array(3);
//  array_copy(dq, T, 3);
//  array_copy(dq+3, Phi, 3);
//  double phi = array_norm(Phi, 3);
//
//  double ** xyz_new = init_matrix(3,3);
//
//  for(int i=0; i<3; i++)
//    for(int j=0; j<3; j++)
//      xyz_new[i][j] = xyz_geom[i][j] - com[j];
//
//  rotate_vecs(Phi, phi, xyz_new, 3);
//
//  for(int i=0; i<3; i++)
//    for(int j=0; j<3; j++)
//      xyz_geom[i][j] = xyz_new[i][j] + com[j] + T[j];
//
//
//fprintf(outfile, "new geometry: \n");
//print_matrix(outfile, xyz_geom, 3, 3);
//
//fprintf(outfile, "\n         Phi            T\n");
//for(int i=0; i<3; i++)
//  fprintf(outfile, "%15.8f%15.8f\n", phi*Phi[i], T[i]);
}
//****AVC****//

}

