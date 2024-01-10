/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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
 ** \file
 ** \ingroup PSIO
 */

#include <iostream>
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsio/psio.h"

namespace psi {

/*!
 ** \ingroup PSIO
 **
 ** PSIO_ERROR(): Print out an error message for libpsio.
 **
 ** \param unit     = file number
 ** \param errval   = error code (defined symbolically, PSIO_ERROR_XXX)
 ** \param prev_msg = (optional) A previous error message that should be prepended to the text from this function.
 **
 */
void psio_error(size_t unit, size_t errval, std::string prev_msg /* = ""*/) {
    std::cerr << "PSIO_ERROR: unit = " << unit << ", errval = " << errval << std::endl;
    // Try to save the TOCs for all open units
    // psio_tocwrite() may end up indirectly calling psio_error() again if a write keeps failing, possibly leading to
    // infinite recursion. To avoid this, the TOCs are not saved if psio_error() is called with a write error.
    if (errval != PSIO_ERROR_WRITE){
        for (int i = 0; i < PSIO_MAXUNIT; i++) psio_tocwrite(i);
    }
    if ((prev_msg.length() > 0) && (prev_msg.back() != '\n')) {
        prev_msg += '\n';
    }
    switch (errval) {
        case PSIO_ERROR_INIT:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_INIT) + " (I/O inititalization failed)\n";
            break;
        case PSIO_ERROR_DONE:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_DONE) + " (I/O cleanup failed)\n";
            break;
        case PSIO_ERROR_MAXVOL:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_MAXVOL) + " (maximum number of volumes exceeded)\n";
            break;
        case PSIO_ERROR_NOVOLPATH:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_NOVOLPATH) + " (no volume path given)\n";
            break;
        case PSIO_ERROR_IDENTVOLPATH:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_IDENTVOLPATH) + " (two identical volume paths)\n";
            break;
        case PSIO_ERROR_OPEN:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_OPEN) + " (open call failed)\n\n"
                        " Check the location of your scratch directory which can be\n"
                        " specified via the $PSI_SCRATCH environment variable.\n"
                        " A fast, local (non-network) scratch disk is strongly preferred.\n\n"
                        " Please note that the scratch directory must exist, be\n"
                        " writable by Psi4, and have available space.\n";
            break;
        case PSIO_ERROR_REOPEN:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_REOPEN) + " (file is already open)\n";
            break;
        case PSIO_ERROR_CLOSE:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_CLOSE) + " (file close failed)\n";
            break;
        case PSIO_ERROR_RECLOSE:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_RECLOSE) + " (file is already closed)\n";
            break;
        case PSIO_ERROR_OSTAT:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_OSTAT) + " (invalid status flag for file open)\n";
            break;
        case PSIO_ERROR_LSEEK:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_LSEEK) + " (lseek failed)\n";
            break;
        case PSIO_ERROR_NOTOCENT:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_NOTOCENT) + " (no such TOC entry)\n";
            break;
        case PSIO_ERROR_TOCENTSZ:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_TOCENTSZ) + " (TOC entry size mismatch)\n";
            break;
        case PSIO_ERROR_KEYLEN:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_KEYLEN) + " (TOC key too long)\n";
            break;
        case PSIO_ERROR_BLKSIZ:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_BLKSIZ) + " (Requested blocksize invalid)\n";
            break;
        case PSIO_ERROR_BLKSTART:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_BLKSTART) + " (Incorrect block start address)\n";
            break;
        case PSIO_ERROR_BLKEND:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_BLKEND) + " (Incorrect block end address)\n";
            break;
        case PSIO_ERROR_WRITE:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_WRITE) + " (error writing to file)\n";
            break;
        case PSIO_ERROR_MAXUNIT:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_MAXUNIT) + " (Maximum unit number exceeded)\n";
            prev_msg += " Open failed because unit " + std::to_string(unit) + " exceeds PSIO_MAXUNIT = ";
            prev_msg += std::to_string(PSIO_MAXUNIT) + ".\n";
            break;
        case PSIO_ERROR_UNOPENED:
            prev_msg += "PSIO_ERROR: " + std::to_string(PSIO_ERROR_UNOPENED) + " (File not opened)\n";
            prev_msg += " You need to open file " + std::to_string(unit) +
                        " before you attempt this operation,\n"
                        " If you're a user, contact developers immediately. This is a bug.\n"
                        " If you're a developer, get yourself some coffee.\n";
            break;
    }
    std::cerr << prev_msg << std::endl;
    throw PSIEXCEPTION(prev_msg);
}

}  // namespace psi
