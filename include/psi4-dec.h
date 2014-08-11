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

#ifndef psi_include_psi4_dec_h
#define psi_include_psi4_dec_h
#include <boost/shared_ptr.hpp>
#include <boost/current_function.hpp>

#include <string>
#include <list>
#include <map>

#include <exception.h>
#include <libmints/typedefs.h>
#include "libparallel/local.h"
#include "libparallel/mpi_wrapper.h"
#include "process.h"
#include "libparallel/ParallelPrinter.h"

#include <cstdarg>
namespace psi {

enum PsiReturnType {Success, Failure, Balk, EndLoop};
/*extern int fprintf(FILE* __restrict __stream,
      const char * __restrict __format, ...);
extern int vfprintf(FILE* __restrict __stream,
      const char * __restrict __format,va_list& args);*/
//
//  extern PSIO *psio;
extern char *psi_file_prefix;
extern std::string outfile_name;
extern bool verbose;

// Very useful regex for matching floating point numbers
#define NUMBER "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

#if HAVE_MPI
   typedef MPICommunicator worldcomm;
#else
   typedef LocalCommWrapper worldcomm;
#endif

extern boost::shared_ptr<worldcomm> WorldComm;
extern boost::shared_ptr<PsiOutStream> outfile;
void die_if_not_converged();

}


#endif
