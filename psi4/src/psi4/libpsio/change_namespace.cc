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
 \ingroup PSIO
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"
namespace psi {

void PSIO::change_file_namespace(unsigned int unit, const std::string & ns1, const std::string & ns2) {
    char *old_name, *new_name, *old_fullpath, *new_fullpath;
    _default_psio_lib_->get_filename(unit, &old_name, true);
    _default_psio_lib_->get_filename(unit, &new_name, true);
    //_default_psio_lib_->get_volpath(unit, 0, &path);
    std::string tpath = PSIOManager::shared_object()->get_file_path(unit);
    const char* path = tpath.c_str();

    old_fullpath = (char*) malloc( (strlen(path)+strlen(old_name)+80)*sizeof(char));
    new_fullpath = (char*) malloc( (strlen(path)+strlen(new_name)+80)*sizeof(char));

    if (ns1 == "") {
        sprintf(old_fullpath, "%s%s.%u", path, old_name, unit);
    } else {
        sprintf(old_fullpath, "%s%s.%s.%u", path, old_name, ns1.c_str(), unit);
    }
    if (ns2 == "") {
        sprintf(new_fullpath, "%s%s.%u", path, new_name, unit);
    } else {
        sprintf(new_fullpath, "%s%s.%s.%u", path, new_name, ns2.c_str(), unit);
    }

    //printf("%s\n",old_fullpath);
    //printf("%s\n",new_fullpath);

    PSIOManager::shared_object()->move_file(std::string(old_fullpath), std::string(new_fullpath));
        ::rename(old_fullpath,new_fullpath);

    free(old_fullpath);
    free(new_fullpath);
}

}
