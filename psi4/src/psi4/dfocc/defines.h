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

#ifndef DEFINES_
#define DEFINES_

#define ID(x) ints->DPD_ID(x)
#define index2(i,j) ((i>j) ? ((i*(i+1)/2)+j) : ((j*(j+1)/2)+i))
#define index4(i,j,k,l) index2(index2(i,j),index2(k,l))
#define index3(i,j,k) (i*(i+1)*(i+2)/6) + (j*(j+1)/2) + i
#define idx2(i,j,Nj) j + (i*Nj)
#define idx3(i,j,k,Nj,Nk) k + (j*Nk) + (i*Nj*Nk)
#define idx4(i,j,k,l,Nj,Nk,Nl) l + (k*Nl) + (j*Nk*Nl) + (i*Nj*Nk*Nl)
#define idx_asym(i,j) ((i>j) ? ((i*(i-1)/2)+j) : ((j*(j-1)/2)+i))

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

#define DIIS_MIN_DET 1.0E-16
#define DIVERGE 1.0E+3

#endif /* DEFINES_ */
