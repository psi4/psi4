/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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
#ifdef _MSC_VER
#include <io.h>
#define SYSTEM_OPEN ::_open
#define SYSTEM_CLOSE ::_close
#define PSIO_OPEN_OLD_FLAGS _O_BINARY | _O_CREAT | _O_RDWR
#define PSIO_OPEN_NEW_FLAGS _O_BINARY | _O_CREAT | _O_RDWR | _O_TRUNC
#define PERMISSION_MODE _S_IWRITE
#else
#include <unistd.h>
#define SYSTEM_OPEN ::open
#define SYSTEM_CLOSE ::close
#define PSIO_OPEN_OLD_FLAGS O_CREAT | O_RDWR
#define PSIO_OPEN_NEW_FLAGS O_CREAT | O_RDWR | O_TRUNC
#define PERMISSION_MODE S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH
#endif
#include <string>
#include <map>
#include <sstream>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"
namespace psi {

void PSIO::open(size_t unit, int status) {
    size_t i;
    char *name, *path;
    psio_ud* this_unit;

    /* check for too large unit */
    if (unit > PSIO_MAXUNIT) psio_error(unit, PSIO_ERROR_MAXUNIT);

    this_unit = &(psio_unit[unit]);

    /* Get number of volumes to stripe across */
    this_unit->numvols = get_numvols(unit);
    if (this_unit->numvols > PSIO_MAXVOL) psio_error(unit, PSIO_ERROR_MAXVOL);
    if (!(this_unit->numvols)) this_unit->numvols = 1;

    /* Check to see if this unit is already open */
    for (i = 0; i < this_unit->numvols; i++) {
        if (this_unit->vol[i].stream != -1) psio_error(unit, PSIO_ERROR_REOPEN);
    }

    /* Get the file name prefix */
    get_filename(unit, &name);

    // Check if any files will have the same name
    {
        using std::string;
        typedef std::map<string, int> Names;
        Names names;
        for (i = 0; i < this_unit->numvols; i++) {
            std::ostringstream oss;
            get_volpath(unit, i, &path);
            oss << path << name << "." << unit;
            const std::string fullpath = oss.str();
            typedef Names::const_iterator citer;
            citer n = names.find(fullpath);
            if (n != names.end()) psio_error(unit, PSIO_ERROR_IDENTVOLPATH);
            names[fullpath] = 1;
            free(path);
        }
    }

    /* Build the name for each volume and open the file */
    for (i = 0; i < this_unit->numvols; i++) {
        char* fullpath;
        get_volpath(unit, i, &path);

        // A bit of a hack in psio open at the moment, breaks volumes and some error checking
        std::string spath2 = PSIOManager::shared_object()->get_file_path(unit);
        const char* path2 = spath2.c_str();

        fullpath = (char*)malloc((strlen(path2) + strlen(name) + 80) * sizeof(char));
        sprintf(fullpath, "%s%s.%zu", path2, name, unit);
        this_unit->vol[i].path = strdup(fullpath);
        free(fullpath);

        /* Register the file */
        PSIOManager::shared_object()->open_file(std::string(this_unit->vol[i].path), unit);

        /* Now open the volume */
        if (status == PSIO_OPEN_OLD) {
            this_unit->vol[i].stream = SYSTEM_OPEN(this_unit->vol[i].path, PSIO_OPEN_OLD_FLAGS, PERMISSION_MODE);
        } else if (status == PSIO_OPEN_NEW) {
            this_unit->vol[i].stream = SYSTEM_OPEN(this_unit->vol[i].path, PSIO_OPEN_NEW_FLAGS, PERMISSION_MODE);
        } else
            psio_error(unit, PSIO_ERROR_OSTAT);

        if (this_unit->vol[i].stream == -1) psio_error(unit, PSIO_ERROR_OPEN);

        free(path);
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
    size_t i;
    char *name, *path;
    psio_ud* this_unit;

    if (unit > PSIO_MAXUNIT) psio_error(unit, PSIO_ERROR_MAXUNIT);

    this_unit = &(psio_unit[unit]);

    /* Get number of volumes to stripe across */
    this_unit->numvols = get_numvols(unit);
    if (this_unit->numvols > PSIO_MAXVOL) psio_error(unit, PSIO_ERROR_MAXVOL);
    if (!(this_unit->numvols)) this_unit->numvols = 1;

    /* Check to see if this unit is already open, if so, should be good.
       If every volume has a sream value other than -1, it's open */
    bool already_open = true;
    for (i = 0; i < this_unit->numvols; i++) {
        if (this_unit->vol[i].stream == -1) already_open = false;
    }
    if (already_open) return (true);

    /* Get the file name prefix */
    get_filename(unit, &name);
    // printf("%s\n",name);

    // Check if any files will have the same name
    {
        using std::string;
        typedef std::map<string, int> Names;
        Names names;
        for (i = 0; i < this_unit->numvols; i++) {
            std::ostringstream oss;
            get_volpath(unit, i, &path);
            oss << path << name << "." << unit;
            const std::string fullpath = oss.str();
            typedef Names::const_iterator citer;
            citer n = names.find(fullpath);
            if (n != names.end()) psio_error(unit, PSIO_ERROR_IDENTVOLPATH);
            names[fullpath] = 1;
            free(path);
        }
    }

    /* Build the name for each volume and open the file */
    bool file_exists = true;
    for (i = 0; i < this_unit->numvols; i++) {
        char* fullpath;
        int stream;
        get_volpath(unit, i, &path);

        // A bit of a hack in psio open at the moment, breaks volumes and some error checking
        std::string spath2 = PSIOManager::shared_object()->get_file_path(unit);
        const char* path2 = spath2.c_str();

        fullpath = (char*)malloc((strlen(path2) + strlen(name) + 80) * sizeof(char));
        sprintf(fullpath, "%s%s.%zu", path2, name, unit);

        /* Now open the volume */
        stream = SYSTEM_OPEN(fullpath, O_RDWR);
        /* and close it again, if opening worked */
        if (stream != -1) {
            SYSTEM_CLOSE(stream);
        }

        if (stream == -1) {
            file_exists = false;
        }

        free(path);
        free(fullpath);
    }

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
