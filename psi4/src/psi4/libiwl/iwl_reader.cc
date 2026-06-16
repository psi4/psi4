/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
  \ingroup IWL
*/
#include "iwl_reader.h"
#include "psi4/libpsio/psio.hpp"

namespace psi {

IWLReader::IWLReader(std::shared_ptr<PSIO> psio, int unit)
    : buf_(psio.get(), unit, /*cutoff*/ 0.0, /*oldfile*/ 1, /*readflag*/ 1) {}

bool IWLReader::next(Entry& e) {
    // Skip past any exhausted (or empty) buckets, fetching the next one until
    // we either find an entry or run out of buckets. The empty-bucket loop
    // guards the corner case the 1995 README warned about.
    while (cur_ >= buf_.buffer_count()) {
        if (buf_.last_buffer()) return false;
        buf_.fetch();
        cur_ = 0;
    }

    const Label* labels = buf_.labels();
    const int j = 4 * cur_;
    e.p = static_cast<int>(labels[j + 0]);
    e.q = static_cast<int>(labels[j + 1]);
    e.r = static_cast<int>(labels[j + 2]);
    e.s = static_cast<int>(labels[j + 3]);
    e.value = static_cast<double>(buf_.values()[cur_]);
    ++cur_;
    return true;
}

}  // namespace psi
