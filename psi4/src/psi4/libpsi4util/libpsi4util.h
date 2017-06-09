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

#include <string>
#include <vector>
//#include <algorithm>
#include <sys/time.h>
//#include <sys/types.h>
//#include <sys/stat.h>
//#include <unistd.h>

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

std::string find_and_replace(std::string &source, const std::string &target,
                             const std::string &replace);

void trim_spaces(std::string &str);

template <typename Range1T, typename Range2T>
bool iequals(const Range1T &Input, const Range2T &Test);

std::vector<std::string> split(const std::string &input, const std::string &regex);

/**
 * @brief Compute the Levenshtein distance between two strings
 * @param s1 string to compute against
 * @param s2 string to compute against
 * @return the distance as an unsigned integer
 */
unsigned int edit_distance(const std::string &s1, const std::string &s2);

class Timer {
   public:
    Timer();
    double get();

   private:
    struct timeval ___start, ___end;
    struct timezone ___dummy;
    double delta_time_seconds;
    double delta_time_hours;
    double delta_time_days;
};

void generate_combinations(int n, int k, std::vector<std::vector<int>> &combinations);
}

#endif // _psi_src_lib_libpsi4util_libpsi4util_h_
