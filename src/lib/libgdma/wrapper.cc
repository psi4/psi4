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

#include <boost/python.hpp>
#include <libmints/mints.h>
#include <libparallel/parallel.h>
#include <libparallel/ParallelPrinter.h>
#include "psi4-dec.h"
#include <iostream>

using namespace boost;
using namespace std;

extern "C" {
  void run_gdma(const char* outfilename, const char*datfilename);
};

namespace psi { namespace gdma {

SharedWavefunction gdma(SharedWavefunction ref_wfn, Options & options, const std::string &datfilename)
{
    outfile->Flush();
    run_gdma(outfile_name.c_str(), datfilename.c_str());
    // Reopen the outfile
    if(outfile_name == "stdout"){
        outfile=boost::shared_ptr<PsiOutStream>(new PsiOutStream());
    }
    else{
       outfile=boost::shared_ptr<PsiOutStream>
          (new OutFile(outfile_name,(APPEND)));
    }
    return ref_wfn;
}

}}
