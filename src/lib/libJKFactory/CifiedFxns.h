/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#ifndef SRC_LIB_LIBJKFACTORY_CIFIEDFXNS_H_
#define SRC_LIB_LIBJKFACTORY_CIFIEDFXNS_H_


#ifdef _cplusplus
extern "C"{
#endif
   struct BasisSet;
   void timer_interface_on(char *);
   void timer_interface_off(char *);
   void SetUp();
   int ComputeShellQuartet(struct BasisSet*,int ThreadID,
         int M,int N,int P,int Q,
         double ** Ints);
#ifdef _cplusplus
   }
#endif



#endif /* SRC_LIB_LIBJKFACTORY_CIFIEDFXNS_H_ */
