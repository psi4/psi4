/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#ifndef _psi_src_lib_libpsi4util_libpsi4util_h_
#define _psi_src_lib_libpsi4util_libpsi4util_h_

#include "psi4/pragma.h"

#include <cctype>
#include <algorithm>
#include <chrono>
#include <string>
#include <vector>

namespace psi {

typedef std::vector<std::string> strvec;

std::string file_to_string(std::string const &name);

bool space(char c);

bool not_space(char c);

std::vector<std::string> split(const std::string &str);

std::vector<std::string> split_indices(const std::string &str);

void to_lower(std::string &str);
std::string to_lower_copy(const std::string &str);

void to_upper(std::string &str);
std::string to_upper_copy(const std::string &str);

std::string to_string(const int val);

std::string to_string(const double val);

double to_double(const std::string str);

int to_integer(const std::string inString);

void append_reference(std::string &str, int reference);

std::string add_reference(std::string &str, int reference);

void append_reference(std::string &str, int reference);

std::string find_and_replace(std::string &source, const std::string &target, const std::string &replace);

void trim_spaces(std::string &str);

template <typename Range1T, typename Range2T>
bool iequals(const Range1T &Input, const Range2T &Test) {
    if (std::distance(std::begin(Input), std::end(Input)) != std::distance(std::begin(Test), std::end(Test)))
        return false;

    return std::equal(std::begin(Input), std::end(Input), std::begin(Test),
                      [](unsigned char a, unsigned char b) { return std::tolower(a) == std::tolower(b); });
}

std::vector<std::string> split(const std::string &input, const std::string &regex);

/**
 * @brief Compute the Levenshtein distance between two strings
 * @param s1 string to compute against
 * @param s2 string to compute against
 * @return the distance as an size_teger
 */
size_t edit_distance(const std::string &s1, const std::string &s2);

class PSI_API Timer {
   public:
    Timer();
    double get();

   private:
    std::chrono::high_resolution_clock::time_point start;
};

void generate_combinations(int n, int k, std::vector<std::vector<int>> &combinations);

/// @brief Converts an input type that is supported by std::to_string, to an std::string that is padded with spaces to a
/// specified minimum width. This mimics some of the functionality of C-style formatting, eg. %3d.
/// @tparam T : type of the variable to be converted
/// @param input : variable to be converted
/// @param width : pad the string with spaces until it is at least this wide, in a right-justified fashion
/// @return std::string representation of the input variable, padded if required
template <typename T>
std::string to_str_width(const T &input, const size_t width) {
    std::string str = std::to_string(input);
    while (str.length() < width) {
        str.insert(str.begin(), ' ');
    }
    return str;
}
}  // namespace psi

#endif  // _psi_src_lib_libpsi4util_libpsi4util_h_
