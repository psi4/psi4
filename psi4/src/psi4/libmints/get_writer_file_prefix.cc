/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
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
#include "psi4/libpsi4util/process.h"

#include <string>
#ifdef _MSC_VER
#include <process.h>
#define SYSTEM_GETPID ::_getpid
#else
#include <unistd.h>
#define SYSTEM_GETPID ::getpid
#endif

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

std::string PSI_API get_writer_file_prefix(const std::string& molecule_name) {
    const std::string pid = "." + std::to_string(SYSTEM_GETPID());
    const std::string label = Process::environment.options.get_str("WRITER_FILE_LABEL");
    if (!label.empty()) {
        const std::string label_with_pid = label + pid;

        // Keep user-specified paths as-is. For bare labels, mirror the directory
        // used by the output file so labeled and unlabeled writer files land together.
        if (label.find_last_of("/\\") != std::string::npos || outfile_name == "stdout") {
            return label_with_pid;
        }

        const auto out_sep = outfile_name.find_last_of("/\\");
        if (out_sep == std::string::npos) {
            return label_with_pid;
        }

        return outfile_name.substr(0, out_sep + 1) + label_with_pid;
    }

    // If no available options WRITER_FILE_LABEL, then we build a default:
    // Get the basename of the output filename, append any active molecule name
    // to it, and return the resulting string
    std::string prefix = outfile_name.substr(0, outfile_name.find_last_of('.'));
    if (!molecule_name.empty()) {
        prefix += "." + molecule_name;
    }
    return prefix + pid;
}

}  // namespace psi
