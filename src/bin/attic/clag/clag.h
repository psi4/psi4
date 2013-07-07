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

/*! \file
    \ingroup CLAG
    \brief Header file for CI Lagrangian program
*/

#ifndef _psi3_bin_clag_clag_h_
#define _psi3_bin_clag_clag_h_

namespace psi { namespace clag {

double **rdopdm(int nbf, int print_lvl, int opdm_file) ;  
double *rdtpdm(int nbf, int print_lvl, int tpdm_file) ; 
void init_io(int argc, char **argv) ; 
void close_io(void) ;
void trace_opdm(double **opdm, int nbf);
void trace_tpdm(double *tpdm, int nbf);
double lagcalc(double **OPDM, double *TPDM, double *h, double *TwoElec, 
  double **lag, int nmo, int npop, int print_lvl, int lag_file); 
void ci_energy(double **OPDM, double *TPDM, double *h, double *TwoElec,
  int nbf, double enuc, double eci_30, double lagtr); 
void compute_zorb (double *TwoElec,
  double **lag, double *epsilon, IndepPairs &IndPairs,
  int nmo, int nocc, int nvir, double *Zvec, double **Z_tilde);
void relax_pdms(double **onepdm, double *tpdm, double *tei,
  double **X_tilde, double *epsilon, IndepPairs &IndPairs,
  int nmo, int nocc, int npop, double *Zvec, int opdm_file, 
  int lag_file, int tpdm_file);
}} // end namespace psi::clag

#endif // header guard
