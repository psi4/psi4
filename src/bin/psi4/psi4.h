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

#ifndef __psi_psi4_psi4_h_
#define __psi_psi4_psi4_h_

#include <stdio.h>
#include <string>

#ifdef MAIN
#define EXT
#else
#define EXT extern
#endif

namespace psi {

EXT FILE* infile;

/*! Directory location of the input file */
EXT std::string infile_directory;

/*! Verbosity */
EXT bool verbose;

/*! sanity check boolean */
EXT bool check_only;

/*! Leave psi temp file */
EXT bool messy;

/*! clean-up */
EXT bool clean_only;

}

#endif
