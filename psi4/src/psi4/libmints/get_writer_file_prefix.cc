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
 \file
 \ingroup MINTS
 */

#include "psi4/psi4-dec.h"
#include "psi4/libmints/molecule.h"

#include <string>

namespace psi {


/*!
** get_writer_file_prefix()
**
** For archival files such as the Hessian, Molden files, etc., we want
** to write a file in the current working directory with a name like
** benzene.molden, benzene.hess, etc.  This function returns the basename
** of such files.  The basename is determined as follows: If the option
** "WRITER_FILE_LABEL" exists, we will use that.  It could be user-specified,
** or created by a python driver doing a complex task (e.g.,
** S22-1-dimer-mp2-ccpvdz).  If the label option is blank (default),
** then we will just make up a prefix using the basename of the output
** file along with the active molecule name (e.g., test.h2o).
**
** C. David Sherrill and Lori Burns
**
** Returns: C++ standard string of the filename prefix
**
** \ingroup MINTS
*/

std::string get_writer_file_prefix(const std::string& molecule_name)
{
    const std::string pid="."+to_string(getpid());
    const std::string label = Process::environment.options.get_str("WRITER_FILE_LABEL");
    if (!label.empty()) {
        return label+pid;
    }

    // If no available options WRITER_FILE_LABEL, then we build a defult:
    // Get the basename of the output filename, append any active molecule name
    // to it, and return the resulting string
    std::string prefix=outfile_name.substr(0, outfile_name.find_last_of('.'));
    if (!molecule_name.empty()) {
        prefix += "." + molecule_name;
    }
    return prefix+pid;

}


}
