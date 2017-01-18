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

#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <cstdio>
#include <sstream>
#include <libgen.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>

#include "libpsi4util.h"

namespace psi {

/*********************
 String manipulation
 ********************/

std::vector <std::string> split(const std::string &str)
{
    // Split a string
    typedef std::string::const_iterator iter;
    std::vector <std::string> splitted_string;
    iter i = str.begin();
    while (i != str.end()) {
        // Ignore leading blanks
        i = find_if(i, str.end(), not_space);
        // Find the end of next word
        iter j = find_if(i, str.end(), space);
        // Copy the characters in [i,j)
        if (i != str.end())
            splitted_string.push_back(std::string(i, j));
        i = j;
    }
    return (splitted_string);
}

bool opening_square_bracket(char c);

bool closing_square_bracket(char c);

std::vector <std::string> split_indices(const std::string &str)
{
    // Split a string
    typedef std::string::const_iterator iter;
    strvec splitted_string;
    iter i = str.begin();
    while (i != str.end()) {
        // Ignore leading blanks
        i = std::find_if(i, str.end(), opening_square_bracket);
        // Find the end of next word
        iter j = std::find_if(i, str.end(), closing_square_bracket);
        // Copy the characters in [i,j]
        if (i != str.end())
            splitted_string.push_back(std::string(i, j + 1));
        i = j;
    }
    return (splitted_string);
}

bool opening_square_bracket(char c)
{
    return (c == '[');
}

bool closing_square_bracket(char c)
{
    return (c == ']');
}

bool space(char c)
{
    return isspace(c);
}

bool not_space(char c)
{
    return !isspace(c);
}


std::string find_and_replace(std::string &source, const std::string &target, const std::string &replace)
{
    std::string str = source;
    std::string::size_type pos = 0;   // where we are now
    std::string::size_type found;     // where the found data is

    if (target.size() > 0)   // searching for nothing will cause a loop
    {
        while ((found = str.find(target, pos)) != std::string::npos) {
            str.replace(found, target.size(), replace);
            pos = found + replace.size();
        }
    }
    return str;
}

void trim_spaces(std::string &str)
{
    // Trim Both leading and trailing spaces
    size_t startpos = str.find_first_not_of(" \t"); // Find the first character position after excluding leading blank spaces
    size_t endpos = str.find_last_not_of(" \t"); // Find the first character position from reverse af

    // if all spaces or empty return an empty string
    if ((std::string::npos == startpos) || (std::string::npos == endpos)) {
        str = "";
    } else
        str = str.substr(startpos, endpos - startpos + 1);
}

unsigned int edit_distance(const std::string &s1, const std::string &s2)
{
    const size_t len1 = s1.size(), len2 = s2.size();
    std::vector < std::vector < unsigned
    int > > d(len1 + 1, std::vector < unsigned
    int > (len2 + 1));

    d[0][0] = 0;
    for (unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
    for (unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

    for (unsigned int i = 1; i <= len1; ++i) {
        for (unsigned int j = 1; j <= len2; ++j) {
            d[i][j] = std::min(std::min(d[i - 1][j] + 1, d[i][j - 1] + 1),
                               d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1));
        }
    }
    return d[len1][len2];
}

/*****************
 String conversion
 *****************/

void to_lower(std::string &str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
}

std::string to_lower_copy(const std::string& str)
{
    std::string cp = str;
    to_lower(cp);
    return cp;
}

void to_upper(std::string &str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

std::string to_upper_copy(const std::string& str)
{
    std::string cp = str;
    to_upper(cp);
    return cp;
}

double to_double(const std::string str)
{
    return std::atof(str.c_str());
}

std::string to_string(const int val)
{
    std::stringstream strm;
    strm << val;
    return strm.str();
}

std::string to_string(const double val)
{
    std::stringstream strm;
    strm << std::setprecision(25) << std::setw(35) << val;
    return strm.str();
}

int to_integer(const std::string inString)
{
    int i = 0;
    char *end;
    i = static_cast<int>(std::strtod(inString.c_str(), &end));
    return i;
}

std::string add_reference(std::string &str, int reference)
{
    return (str + "{" + to_string(reference) + "}");
}

void append_reference(std::string &str, int reference)
{
    str += "{" + to_string(reference) + "}";
}

/*********************************************************
 Memory Allocation
 *********************************************************/

/**
 * Convert the size of a doubles array in Mb using the definition 1Mb = 1048576 bytes
 * @param n size of the array
 * @return
 */
double to_MB(size_t n)
{
    return (double(n * sizeof(double)) / 1048576.0);
    // Using this definition 1 Mb has ca. 5% more than 1000000 bytes
}

unsigned long int init_smatrix(short **&matrix, int size1, int size2)
{
    unsigned long int size, uli_size1, uli_size2;
    uli_size1 = static_cast<unsigned long int>(size1);
    uli_size2 = static_cast<unsigned long int>(size2);
    size = uli_size1 * uli_size2;
    if (!uli_size1 || !uli_size2) {
        matrix = NULL;
    } else {
        matrix = new short *[uli_size1];
        short *vector = new short[size];
        for (unsigned long int i = 0; i < size; i++) vector[i] = 0;
        for (unsigned long int i = 0; i < uli_size1; i++)
            matrix[i] = &(vector[i * uli_size2]);
    }
    return (size * sizeof(short));
}

unsigned long int free_smatrix(short **&matrix, int size1, int size2)
{
    unsigned long int size, uli_size1, uli_size2;
    uli_size1 = static_cast<unsigned long int>(size1);
    uli_size2 = static_cast<unsigned long int>(size2);
    size = uli_size1 * uli_size2;
    if (matrix == NULL) return (0);
    delete[] matrix[0];
    delete[] matrix;
    return (size * sizeof(short));
}

unsigned long int init_smatrix(short ***&matrix, int size1, int size2, int size3)
{
    unsigned long int size = static_cast<unsigned long int>(size1 * size2 * size3);
    matrix = new short **[size1];
    for (int i = 0; i < size1; i++) {
        matrix[i] = new short *[size2];
    }
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++) {
            matrix[i][j] = new short[size3];
        }
    }
    return (size * sizeof(short));
}

unsigned long int free_smatrix(short ***matrix, int size1, int size2, int size3)
{
    unsigned long int size = static_cast<unsigned long int>(size1 * size2 * size3);
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++) {
            delete[] matrix[i][j];
        }
    }
    for (int i = 0; i < size1; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
    return (size * sizeof(short));
}

}
