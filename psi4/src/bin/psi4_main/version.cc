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
#include <ctime>
#include "psi4/src/lib/libparallel/parallel.h"
#include <cstdio>
#include <cstring>
#include <sstream>

#include "psi4/include/psi4-dec.h"
#include "psi4/src/lib/libparallel/ParallelPrinter.h"
#include "gitversion.h"
namespace psi {
void print_version(std::string);
void print_version(std::string)
{
  boost::shared_ptr<PsiOutStream> printer=outfile;
  printer->Printf( "    -----------------------------------------------------------------------\n");
  printer->Printf( "          Psi4: An Open-Source Ab Initio Electronic Structure Package\n");
#ifdef PSI_VERSION
  printer->Printf( "                              Psi4 %s Driver\n", PSI_VERSION);
#endif

  // Are we using git? If so,what version string
#ifdef GIT_VERSION
  printer->Printf( "\n                          Git: Rev " GIT_VERSION "\n");
#endif

  printer->Printf( "\n");
  printer->Printf( "    J. M. Turney, A. C. Simmonett, R. M. Parrish, E. G. Hohenstein,\n");
  printer->Printf( "    F. A. Evangelista, J. T. Fermann, B. J. Mintz, L. A. Burns, J. J. Wilke,\n");
  printer->Printf( "    M. L. Abrams, N. J. Russ, M. L. Leininger, C. L. Janssen, E. T. Seidl,\n");
  printer->Printf( "    W. D. Allen, H. F. Schaefer, R. A. King, E. F. Valeev, C. D. Sherrill,\n");
  printer->Printf( "    and T. D. Crawford, WIREs Comput. Mol. Sci. 2, 556-565 (2012)\n");
  printer->Printf( "    (doi: 10.1002/wcms.93)\n");

  printer->Printf( "\n");
  printer->Printf( "                         Additional Contributions by\n");
  printer->Printf( "    A. E. DePrince, U. Bozkaya, A. Yu. Sokolov, D. G. A. Smith, R. Di Remigio,\n");
  printer->Printf( "    R. M. Richard, J. F. Gonthier, H. R. McAlexander, M. Saitow, and\n");
  printer->Printf( "    B. P. Pritchard\n");
  printer->Printf( "    -----------------------------------------------------------------------\n\n");
  printer->Printf("\n");

  std::time_t cur_time = time(0);
  printer->Printf( "    Psi4 started on: %s\n",ctime(&cur_time));

  pid_t pid = getpid();
  printer->Printf( "    Process ID: %6d\n",pid);
  printer->Printf( "    PSI4DATADIR: %s\n", Process::environment("PSIDATADIR").c_str());
}

}
