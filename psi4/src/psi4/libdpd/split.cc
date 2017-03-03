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

/*! \file
    \ingroup LIBDPD
    \brief String parsing code for orbital indices, taken from Justin Turney's AMBIT code.
*/

#include <algorithm>
#include <iterator>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

namespace psi {

// trim from start
static inline string &dpd_ltrim(string &s)
{
    s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
    return s;
}

// trim from end
static inline string &dpd_rtrim(string &s)
{
    s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline string &dpd_trim(string &s)
{
    return dpd_ltrim(dpd_rtrim(s));
}

/** Takes a string of indices and splits them into a vector of strings.
 *
 * If a comma is found in indices then they are split on the comma.
 * If no comma is found it assumes the indices are one character in length.
 *
 */
vector<string> dpd_split(const string &indices)
{
    istringstream f(indices);
    string s;
    vector<string> v;

    if (indices.find(",") != string::npos) {
        while (getline(f, s, ',')) {
            string trimmed = dpd_trim(s);
            v.push_back(trimmed);
        }
    }
    else {
        // simply split the string up
        for (size_t i = 0; i < indices.size(); ++i)
            v.push_back(string(1, indices[i]));
    }

    return v;
}

} // namespace psi
