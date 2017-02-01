/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*!
   \file cov_radii.h : covalent radii
   \ingroup optking
*/

/*
* This header file contains the covalent radii of the atoms in Angstroms from:
* "Covalent radii revisited", Dalton Trans., 2008, 2832, by
* B. Cordero, V. Gómez, A. E. Platero-Prats, M. Revés, J. Echeverría,
* E. Cremades, F. Barragán and S. Alvarez
*
* The largest values have been chosen for C (sp3) and for
* Mn, Fe, and Co (high-spin).
*  - RAK, May 2008
*/

#ifndef _cov_radii_h_
#define _cov_radii_h_

#define LAST_COV_RADII_INDEX 96

const double cov_radii[] = { 2.0,  /* ghost ? */
0.31, //H
0.28, //He
1.28, //Li
0.96, //Be
0.84, //B
0.76, //C
0.71, //N
0.66, //O
0.57, //F
0.58, //Ne
1.66, //Na
1.41, //Mg
1.21, //Al
1.11, //Si
1.07, //P
1.05, //S
1.02, //Cl
1.06, //Ar
2.03, //K
1.76, //Ca
1.70, //Sc
1.60, //Ti
1.53, //V
1.39, //Cr
1.61, //Mn
1.52, //Fe
1.50, //Co
1.24, //Ni
1.32, //Cu
1.22, //Zn
1.22, //Ga
1.20, //Ge
1.19, //As
1.20, //Se
1.20, //Br
1.16, //Kr
2.20, //Rb
1.95, //Sr
1.90, //Y
1.75, //Zr
1.64, //Nb
1.54, //Mo
1.47, //Tc
1.46, //Ru
1.42, //Rh
1.39, //Pd
1.45, //Ag
1.44, //Cd
1.42, //In
1.39, //Sn
1.39, //Sb
1.38, //Te
1.39, //I
1.40, //Xe
2.44, //Cs
2.15, //Ba
2.07, //La
2.04, //Ce
2.03, //Pr
2.01, //Nd
1.99, //Pm
1.98, //Sm
1.98, //Eu
1.96, //Gd
1.94, //Tb
1.92, //Dy
1.92, //Ho
1.89, //Er
1.90, //Tm
1.87, //Yb
1.87, //Lu
1.75, //Hf
1.70, //Ta
1.62, //W
1.51, //Re
1.44, //Os
1.41, //Ir
1.36, //Pt
1.36, //Au
1.32, //Hg
1.45, //Tl
1.46, //Pb
1.48, //Bi
1.40, //Po
1.50, //At
1.50, //Rn
2.60, //Fr
2.21, //Ra
2.15, //Ac
2.06, //Th
2.00, //Pa
1.96, //U
1.90, //Np
1.87, //Pu
1.80, //Am
1.69 //Cm
};

#endif /* header guard */