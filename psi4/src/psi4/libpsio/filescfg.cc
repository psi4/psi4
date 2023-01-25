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

#include "psi4/pragma.h"
PRAGMA_WARNING_PUSH
PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
#include <memory>
PRAGMA_WARNING_POP
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsi4util/libpsi4util.h"

#include <string>
#include <map>
#include <sstream>
#include <algorithm>
#include <cctype>

#include <iostream>
#include <ostream>

namespace psi {

std::string fullkwd(const char* kwdgrp, const char* kwd, int unit) {
    std::string unitname;
    if (unit < 0)
        unitname = "DEFAULT";
    else {
        std::ostringstream oss;
        oss << "FILE" << unit;
        unitname = oss.str();
    }
    const std::string sep(":");

    std::string fkwd = sep + kwdgrp + sep + "FILES" + sep + unitname + sep + kwd;
    // convert to upper case
    to_upper(fkwd);
    return fkwd;
}

void PSIO::filecfg_kwd(const char* kwdgrp, const char* kwd, int unit, const char* kwdval) {
    std::string fkwd = fullkwd(kwdgrp, kwd, unit);
    files_keywords_[fkwd] = kwdval;
}

const std::string& PSIO::filecfg_kwd(const char* kwdgrp, const char* kwd, int unit) {
    static std::string nullstr;

    const std::string fkwd = fullkwd(kwdgrp, kwd, unit);
    KWDMap::const_iterator kwd_loc = files_keywords_.find(fkwd);
    if (kwd_loc != files_keywords_.end())
        return kwd_loc->second;
    else
        return nullstr;
}
}  // namespace psi
