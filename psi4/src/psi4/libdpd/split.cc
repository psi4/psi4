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

/*! \file
    \ingroup LIBDPD
    \brief String parsing code for orbital indices, taken from Justin Turney's AMBIT code.
*/

#include <cctype>
#include <algorithm>
#include <iterator>
#include <string>
#include <sstream>
#include <vector>
#include <functional>

namespace psi {

// trim from start
static inline std::string &dpd_ltrim(std::string &s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(), [](int c) {return !std::isspace(c);}));
    return s;
}

// trim from end
static inline std::string &dpd_rtrim(std::string &s) {
    s.erase(find_if(s.rbegin(), s.rend(), [](int c) {return !std::isspace(c);}).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &dpd_trim(std::string &s) { return dpd_ltrim(dpd_rtrim(s)); }

/** Takes a string of indices and splits them into a vector of strings.
 *
 * If a comma is found in indices then they are split on the comma.
 * If no comma is found it assumes the indices are one character in length.
 *
 */
std::vector<std::string> dpd_split(const std::string &indices) {
    std::istringstream f(indices);
    std::string s;
    std::vector<std::string> v;

    if (indices.find(",") != std::string::npos) {
        while (std::getline(f, s, ',')) {
            std::string trimmed = dpd_trim(s);
            v.push_back(trimmed);
        }
    } else {
        // simply split the string up
        for (size_t i = 0; i < indices.size(); ++i) v.push_back(std::string(1, indices[i]));
    }

    return v;
}

}  // namespace psi
