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

#ifndef psi_include_psi4_dec_h
#define psi_include_psi4_dec_h

#include <boost/shared_ptr.hpp>
#include <string>

#include "psi4/src/lib/libparallel/PsiOutStream.h"
#include "psi4/src/lib/libparallel/process.h"

namespace psi {
namespace LibParallel{
class ParallelEnvironment;
}
class PsiOutStream;
enum PsiReturnType {Success, Failure, Balk, EndLoop};
extern char *psi_file_prefix;
extern std::string outfile_name;
extern bool verbose;
extern std::string restart_id;

// Very useful regex for matching floating point numbers
#define NUMBER "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

extern boost::shared_ptr<LibParallel::ParallelEnvironment> WorldComm;
extern boost::shared_ptr<PsiOutStream> outfile;
void die_if_not_converged();
}//End namespace psi
#endif
