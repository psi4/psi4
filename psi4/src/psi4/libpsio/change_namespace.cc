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
#include <cstdlib>
#include <cstring>
#include "psi4/pragma.h"
#include <memory>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"
namespace psi {

void PSIO::change_file_namespace(size_t unit, const std::string& ns1, const std::string& ns2) {
    // The namespace is spliced into the filename here (unlike get_unit_filename),
    // so compose the paths directly. Both endpoints share the same base name.
    const std::string name = _default_psio_lib_->get_filename(unit, true);
    const std::string path = PSIOManager::shared_object()->get_file_path(unit);
    const std::string suffix = "." + std::to_string(unit);

    const std::string old_fullpath = path + name + (ns1.empty() ? "" : "." + ns1) + suffix;
    const std::string new_fullpath = path + name + (ns2.empty() ? "" : "." + ns2) + suffix;

    PSIOManager::shared_object()->move_file(old_fullpath, new_fullpath);
    ::rename(old_fullpath.c_str(), new_fullpath.c_str());
}

}  // namespace psi
