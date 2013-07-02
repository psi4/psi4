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

/*
** GENCI.H
** 
** Contains the function prototypes for the C routines in the GENCI
** module.
**
*/

#ifndef _psi_src_bin_detci_genci_h
#define _psi_src_bin_detci_genci_h

namespace psi { namespace detci {

int schmidt_addoc(double *buffer4, double *buffer5, int buf_size, 
                  int extra_buf, int num_buf, PSI_FPTR d_index, 
                  int N, int L, int b_file, int d_file);
void v_normalize(double *A, PSI_FPTR index, int buf_size, 
                 int extra_buf, int num_buf, int d_file);
double *v_schmidt(double *buffer4, double *buffer5, int buf_size, 
                  int extra_buf, int num_buf, int N, int L, int b_file);
void det2strings(BIGINT det, int *alp_code, int *alp_idx,
                 int *bet_code, int *bet_idx);
BIGINT strings2det(int alp_code, int alp_idx, int bet_code, int bet_idx);
void unit_guess(int alp_code, int alp_idx, int bet_code, int bet_idx,
                int switch_buf3, double *buffer, int buf_size,
                int num_buf, int extra_buf, PSI_FPTR b_file,
                PSI_FPTR b_writ, int M, int N);
void max_element(double *buffer, int num_elements, double *max, int *max_num);
void min_element(double *buffer, int num_elements, double *min, int *min_num);
void read_c(int switch_buf3, double *buffer, int buf_size, int num_buf,
            int extra_buf, int b_file, PSI_FPTR b_writ,
            int c_file, PSI_FPTR c_index);

}} // namespace psi::detci

#endif // header guard

