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

#ifndef _psi3_libchkpt_config_h_
#define _psi3_libchkpt_config_h_

#ifndef MAX_ELEMNAME
#define MAX_ELEMNAME 13
#endif

#ifndef CHKPT_PREFIX_LEN
#define CHKPT_PREFIX_LEN 32
#endif

#define CHKPT_ZMAT_LABEL_LEN 20
/*--- Z-matrix entry type ---*/
struct z_entry {
  int bond_atom;            /* first reference atom (bond) */
  int angle_atom;           /* second reference atom (angle) */
  int tors_atom;            /* third reference atom (torsion) */
  int bond_opt;             /* flags indicating to optimize values */
  int angle_opt; 
  int tors_opt; 
  double bond_val;          /* coordinate values */
  double angle_val; 
  double tors_val; 
  char bond_label[CHKPT_ZMAT_LABEL_LEN];      /* variable labels, if any */
  char angle_label[CHKPT_ZMAT_LABEL_LEN]; 
  char tors_label[CHKPT_ZMAT_LABEL_LEN]; 
};

/*--- Types of reference determinants ---*/
typedef enum {ref_rhf = 0, ref_uhf = 1, ref_rohf = 2, ref_tcscf = 3,
	      ref_rks = 4, ref_uks = 5} reftype;

#endif
