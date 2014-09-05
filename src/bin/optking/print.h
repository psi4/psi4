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

/*! \file    print.h
    \ingroup optking
    \brief header for print functions
*/

#ifndef _opt_print_h_
#define _opt_print_h_

#include <cstdlib>
#include <cstdio>
#include <string>
namespace opt {

void print_matrix(const std::string OutFileRMR, double **A, const int x, const int y);

void print_array(const std::string OutFileRMR, double *A, const int x);

void print_geom_array(const std::string OutFileRMR, double *A, const int natom);

}

#endif

