/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

#include <cstdio>
#include <cstdlib>
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/psi4-dec.h"

namespace psi {

/*!
 ** \ingroup PSIO
 **
 ** PSIO_ERROR(): Print out an error message for libpsio.
 **
 ** \param unit   = file number
 ** \param errval = error code (defined symbolically, PSIO_ERROR_XXX)
 **
 */
void psio_error(size_t unit, size_t errval) {
    int i;
    fprintf(stderr, "PSIO_ERROR: unit = %zu, errval = %zu\n", unit, errval);
    /* Try to save the TOCs for all open units */
    /* psio_tocwrite() does not call psio_error() so this is OK */
    for (i = 0; i < PSIO_MAXUNIT; i++) psio_tocwrite(i);
    auto msg = new char[800];
    switch (errval) {
        case PSIO_ERROR_INIT:
            sprintf(msg,"PSIO_ERROR: %d (I/O inititalization failed)\n", PSIO_ERROR_INIT);
            break;
        case PSIO_ERROR_DONE:
            sprintf(msg,"PSIO_ERROR: %d (I/O cleanup failed)\n", PSIO_ERROR_DONE);
            break;
        case PSIO_ERROR_MAXVOL:
            sprintf(msg,"PSIO_ERROR: %d (maximum number of volumes exceeded)\n", PSIO_ERROR_MAXVOL);
            break;
        case PSIO_ERROR_NOVOLPATH:
            sprintf(msg,"PSIO_ERROR: %d (no volume path given)\n", PSIO_ERROR_NOVOLPATH);
            break;
        case PSIO_ERROR_IDENTVOLPATH:
            sprintf(msg,"PSIO_ERROR: %d (two identical volume paths)\n", PSIO_ERROR_IDENTVOLPATH);
            break;
        case PSIO_ERROR_OPEN:
            sprintf(msg,"PSIO_ERROR: %d (open call failed)\n\n"
                        " Check the location of your scratch directory which can be\n"
                        " specified via the $PSI_SCRATCH environment variable.\n"
                        " A local (non-network) scratch disk is strongly preferred.\n\n"
                        " Please note that the scratch directory must exist, be\n"
                        " writable by Psi4, and have available space.\n", PSIO_ERROR_OPEN);
            break;
        case PSIO_ERROR_REOPEN:
            sprintf(msg,"PSIO_ERROR: %d (file is already open)\n", PSIO_ERROR_REOPEN);
            break;
        case PSIO_ERROR_CLOSE:
            sprintf(msg,"PSIO_ERROR: %d (file close failed)\n", PSIO_ERROR_CLOSE);
            break;
        case PSIO_ERROR_RECLOSE:
            sprintf(msg,"PSIO_ERROR: %d (file is already closed)\n", PSIO_ERROR_RECLOSE);
            break;
        case PSIO_ERROR_OSTAT:
            sprintf(msg,"PSIO_ERROR: %d (invalid status flag for file open)\n", PSIO_ERROR_OSTAT);
            break;
        case PSIO_ERROR_LSEEK:
            sprintf(msg,"PSIO_ERROR: %d (lseek failed)\n", PSIO_ERROR_LSEEK);
            break;
        case PSIO_ERROR_NOTOCENT:
            sprintf(msg,"PSIO_ERROR: %d (no such TOC entry)\n", PSIO_ERROR_NOTOCENT);
            // sprintf(msg,"PSIO_ERROR: %d (no such TOC entry)\n", PSIO_ERROR_NOTOCENT);
            break;
        case PSIO_ERROR_TOCENTSZ:
            sprintf(msg,"PSIO_ERROR: %d (TOC entry size mismatch)\n", PSIO_ERROR_TOCENTSZ);
            break;
        case PSIO_ERROR_KEYLEN:
            sprintf(msg,"PSIO_ERROR: %d (TOC key too long)\n", PSIO_ERROR_KEYLEN);
            break;
        case PSIO_ERROR_BLKSIZ:
            sprintf(msg,"PSIO_ERROR: %d (Requested blocksize invalid)\n", PSIO_ERROR_BLKSIZ);
            break;
        case PSIO_ERROR_BLKSTART:
            sprintf(msg,"PSIO_ERROR: %d (Incorrect block start address)\n", PSIO_ERROR_BLKSTART);
            break;
        case PSIO_ERROR_BLKEND:
            sprintf(msg,"PSIO_ERROR: %d (Incorrect block end address)\n", PSIO_ERROR_BLKEND);
            break;
        case PSIO_ERROR_WRITE:
            sprintf(msg,"PSIO_ERROR: %d (error writing to file)\n", PSIO_ERROR_WRITE);
            break;
        case PSIO_ERROR_MAXUNIT:
            sprintf(msg,"PSIO_ERROR: %d (Maximum unit number exceeded)\n"
                        " Open failed because unit %zu exceeds "
                        " PSIO_MAXUNIT = %d.\n",PSIO_ERROR_MAXUNIT, unit, PSIO_MAXUNIT);
            break;
        case PSIO_ERROR_UNOPENED:
            sprintf(msg,"PSIO_ERROR: %d (File not opened)\n" 
                        " You need to open file %zu before you attempt this operation,\n"
                        " If you're a user, contact developers immediately. This is a bug.\n"
                        " If you're a developer, get yourself some coffee.\n", PSIO_ERROR_UNOPENED,  unit);
            break;
    }
    throw PSIEXCEPTION(msg);
}

}  // namespace psi
