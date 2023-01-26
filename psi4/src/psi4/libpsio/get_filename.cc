/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psi4-dec.h"

namespace psi {

void PSIO::get_filename(size_t unit, char **name, bool remove_namespace) {
    std::string kval;
    std::string dot(".");
    std::string ns = dot + pid_;
    ns += (default_namespace_ == "" || remove_namespace) ? "" : dot + default_namespace_;
    // std::string path = PSIOManager::shared_object()->get_file_path(unit);
    // printf("%s %s: %d\n", path.c_str(), ns.c_str(), unit);
    //*name = strdup((path + ns).c_str());
    // return;

    kval = filecfg_kwd("PSI", "NAME", unit);
    if (!kval.empty()) {
        kval = kval + ns;
        *name = strdup(kval.c_str());
        return;
    }
    kval = filecfg_kwd("PSI", "NAME", -1);
    if (!kval.empty()) {
        kval = kval + ns;
        *name = strdup(kval.c_str());
        return;
    }
    kval = filecfg_kwd("DEFAULT", "NAME", unit);
    if (!kval.empty()) {
        kval = kval + ns;
        *name = strdup(kval.c_str());
        return;
    }
    kval = filecfg_kwd("DEFAULT", "NAME", -1);
    if (!kval.empty()) {
        kval = kval + ns;
        *name = strdup(kval.c_str());
        return;
    }

    // assume that the default has been provided already
    abort();
}
}  // namespace psi
