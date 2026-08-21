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
 \ingroup PSIO
 */

#include <cstdio>
#include <fcntl.h>
#include <cstring>
#include <cstdlib>
#include <string>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"
namespace psi {

void PSIO::open(size_t unit, int status) {
    char *name;
    psio_ud* this_unit;

    /* check for too large unit */
    if (unit > PSIO_MAXUNIT) psio_error(unit, PSIO_ERROR_MAXUNIT);

    this_unit = &(psio_unit[unit]);

    /* Check to see if this unit is already open */
    if (this_unit->vol.stream != -1) psio_error(unit, PSIO_ERROR_REOPEN);

    /* Get the file name prefix */
    get_filename(unit, &name);

    /* Build the file name and open the file */
    {
        std::string spath = PSIOManager::shared_object()->get_file_path(unit);
        const char* path = spath.c_str();

        char* fullpath = (char*)malloc((strlen(path) + strlen(name) + 80) * sizeof(char));
        sprintf(fullpath, "%s%s.%zu", path, name, unit);
        this_unit->vol.path = strdup(fullpath);
        free(fullpath);

        /* Register the file */
        PSIOManager::shared_object()->open_file(std::string(this_unit->vol.path), unit);

        /* Now open the file */
        if (status == PSIO_OPEN_OLD) {
            this_unit->vol.stream = SYSTEM_OPEN(this_unit->vol.path, PSIO_OPEN_OLD_FLAGS, PERMISSION_MODE);
        } else if (status == PSIO_OPEN_NEW) {
            this_unit->vol.stream = SYSTEM_OPEN(this_unit->vol.path, PSIO_OPEN_NEW_FLAGS, PERMISSION_MODE);
        } else
            psio_error(unit, PSIO_ERROR_OSTAT);

        if (this_unit->vol.stream == -1) psio_error(unit, PSIO_ERROR_OPEN);
    }

    if (status == PSIO_OPEN_OLD)
        tocread(unit);
    else if (status == PSIO_OPEN_NEW) {
        /* Init the TOC stats and write them to disk */
        this_unit->toclen = 0;
        this_unit->toc = nullptr;
        wt_toclen(unit, 0);
    } else
        psio_error(unit, PSIO_ERROR_OSTAT);

    free(name);
}

// Mirrors PSIO::open() but just check to see if the file is there
// status needs is assumed PSIO_OPEN_OLD if this is called
bool PSIO::exists(size_t unit) {
    char *name;
    psio_ud* this_unit;

    if (unit > PSIO_MAXUNIT) psio_error(unit, PSIO_ERROR_MAXUNIT);

    this_unit = &(psio_unit[unit]);

    /* If the unit is already open, the file exists */
    if (this_unit->vol.stream != -1) return (true);

    /* Get the file name prefix */
    get_filename(unit, &name);

    /* Build the file name and test whether it can be opened */
    std::string spath = PSIOManager::shared_object()->get_file_path(unit);
    const char* path = spath.c_str();

    char* fullpath = (char*)malloc((strlen(path) + strlen(name) + 80) * sizeof(char));
    sprintf(fullpath, "%s%s.%zu", path, name, unit);

    int stream = SYSTEM_OPEN(fullpath, O_RDWR);
    const bool file_exists = (stream != -1);
    /* and close it again, if opening worked */
    if (stream != -1) SYSTEM_CLOSE(stream);

    free(fullpath);
    free(name);
    return (file_exists);
}

void PSIO::rehash(size_t unit) {
    if (open_check(unit)) {
        close(unit, 1);
        open(unit, PSIO_OPEN_OLD);
    }
}

int psio_open(size_t unit, int status) {
    _default_psio_lib_->open(unit, status);
    return 1;
}

}  // namespace psi
