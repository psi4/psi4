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
#include "psi4/libmints/basisset_parser.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/psi4-dec.h"

#include <cstdio>
#include <fstream>
#include <algorithm>
#include <ctype.h>
#include <memory>
#include <regex>

namespace psi {
namespace {
std::regex basis_separator("^\\s*\\[\\s*(.*?)\\s*\\]\\s*$");
}

// the third parameter of from_string() should be
// one of std::hex, std::dec or std::oct
template<class T>
bool from_string(T &t,
                 const std::string &s,
                 std::ios_base &(*f)(std::ios_base &))
{
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
}

BasisSetFileNotFound::BasisSetFileNotFound(std::string message,
                                           const char *_file,
                                           int _line) throw()
        : PsiException(message, _file, _line)
{
    std::stringstream sstr;
    sstr << "sanity check failed! " << message;
    rewrite_msg(sstr.str());
}

BasisSetFileNotFound::~BasisSetFileNotFound() throw()
{
}

BasisSetNotFound::BasisSetNotFound(std::string message,
                                   const char *_file,
                                   int _line) throw()
        : PsiException(message, _file, _line)
{
    std::stringstream sstr;
    sstr << "sanity check failed! " << message;
    rewrite_msg(sstr.str());
}

BasisSetNotFound::~BasisSetNotFound() throw()
{
}

BasisSetParser::BasisSetParser()
{
    force_puream_or_cartesian_ = false;
    forced_is_puream_ = false;
}

BasisSetParser::BasisSetParser(bool forced_puream)
{
    force_puream_or_cartesian_ = true;
    forced_is_puream_ = forced_puream;
}

BasisSetParser::~BasisSetParser()
{
}

std::vector <std::string> BasisSetParser::load_file(const std::string &filename,
                                                    const std::string &basisname)
{
    filename_ = filename;

    // Loads an entire file.
    std::vector <std::string> lines;

    std::smatch what;

    // temp variable
    std::string text;

    // Stream to use
    std::ifstream infile(filename.c_str());

    if (!infile)
        throw BasisSetFileNotFound("BasisSetParser::parse: Unable to open basis set file: " + filename, __FILE__, __LINE__);

    bool given_basisname = basisname.empty() ? false : true;
    bool found_basisname = false;

    while (infile.good()) {
        std::getline(infile, text);

        // If no basisname was given always save the line.
        if (given_basisname == false)
            lines.push_back(text);

        if (found_basisname) {

            // If we find another [*] we're done.
            if (std::regex_match(text, what, basis_separator))
                break;

            lines.push_back(text);
            continue;
        }

        // If the user gave a basisname AND text matches the basisname we want to trigger to retain
        if (given_basisname && regex_match(text, what, basis_separator)) {
            if (iequals(what[1].str(), basisname))
                found_basisname = true;
        }
    }


    return lines;
}

std::vector <std::string> BasisSetParser::string_to_vector(const std::string &input)
{
    std::regex re("\\n");
    std::sregex_token_iterator
            first{input.begin(), input.end(), re, -1},
            last;
    return {first, last};
}

std::vector <ShellInfo>
Gaussian94BasisSetParser::parse(const std::string &symbol, const std::vector <std::string> &lines)
{
    // Regular expressions that we'll be checking for.
    std::regex cartesian("^\\s*cartesian\\s*", std::regex_constants::icase);
    std::regex spherical("^\\s*spherical\\s*", std::regex_constants::icase);
    std::regex comment("^\\s*\\!.*");                                         // line starts with !
    std::regex separator("^\\s*\\*\\*\\*\\*");                                // line starts with ****
    std::regex shell("^\\s*(\\w+)\\s*(\\d+)\\s*(-?\\d+\\.\\d+)");             // Match beginning of contraction
    std::regex atom_array("^\\s*(([A-Z]{1,3})(?:(_\\w+)|(\\d+))?)\\s+0\\s*$", std::regex_constants::icase);  // atomic symbol/label terminated by 0

    // NUMBER is in psi4-dec.h
    std::regex primitives1("^\\s*" NUMBER "\\s+" NUMBER ".*");                // Match s, p, d, f, g, ... functions
    std::regex primitives2("^\\s*" NUMBER "\\s+" NUMBER "\\s+" NUMBER ".*");  // match sp functions

    // s, p and s, p, d can be grouped together in Pople-style basis sets
    const std::string sp("SP"), spd("SPD");

    char mo = (char) (-1);
    //                     a  b  c  d  e  f  g  h  i  j  k  l  m  n  o  p  q  r  s  t  u  v  w  x  y  z
    char shell_to_am[] = {mo, mo, mo, 2, mo, 3, 4, 5, 6, mo, 7, 8, 9, 10, 11, 1, 12, 13, 0, 14, 15, 16, 17, 18, 19, 20};

    // Hold the result of a regex_match
    std::smatch what;

    // Basis type.
    GaussianType gaussian_type = Pure;

    if (force_puream_or_cartesian_) {
        if (forced_is_puream_ == false) gaussian_type = Cartesian;
    }

    // Need a dummy center for the shell.
    Vector3 center;

    std::vector <ShellInfo> shell_list;

    size_t lineno = 0;
    bool found = false;

    while (lineno < lines.size()) {
        std::string line = lines[lineno++];

        // Ignore blank lines
        if (line.empty())
            continue;

        // Look for Cartesian or Spherical
        if (!force_puream_or_cartesian_) {
            if (regex_match(line, what, cartesian)) {
                gaussian_type = Cartesian;
                if (Process::environment.options.get_global("PUREAM").has_changed()) {
                    gaussian_type = ((Process::environment.options.get_global("PUREAM").to_integer()) ? Pure : Cartesian);
                }
                continue;
            } else if (regex_match(line, what, spherical)) {
                gaussian_type = Pure;
                if (Process::environment.options.get_global("PUREAM").has_changed()) {
                    gaussian_type = ((Process::environment.options.get_global("PUREAM").to_integer()) ? Pure : Cartesian);
                }
                continue;
            }
        } // end case where puream setting wasn't forced by caller

        // Do some matches
        if (regex_match(line, what, comment)) {
            continue;
        }
        if (regex_match(line, what, separator)) {
            continue;
        }

        // Match: H    0
        // or:    Al_99     0
        // 11 May 2014 LAB, dropped "H O 0" atom_array capability since code didn't follow through; permitted py-side
        if (regex_match(line, what, atom_array)) {
            // Check the captures and see if this basis set is for the atom we need.
            found = false;
            if (iequals(symbol, what[1].str())) {
                found = true;

                // Read in the next line
                line = lines[lineno++];

                // Need to do the following until we match a "****" which is the end of the basis set
                while (!regex_match(line, what, separator)) {
                    // Match shell information
                    if (regex_match(line, what, shell)) {
                        std::string shell_type(what[1].first, what[1].second);
                        std::transform(shell_type.begin(), shell_type.end(), shell_type.begin(), ::toupper);
                        int nprimitive;
                        double scale;
                        double exponent, contraction;

                        if (!from_string<int>(nprimitive, what[2], std::dec))
                            throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to convert number of primitives:\n" + line);
                        if (!from_string<double>(scale, what[3], std::dec))
                            throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to convert scale factor:\n" + line);

                        if (shell_type.size() == 1) {
                            int am = (int) shell_to_am[shell_type[0] - 'A'];

                            std::vector<double> exponents(nprimitive);
                            std::vector<double> contractions(nprimitive);

                            for (int p = 0; p < nprimitive; ++p) {
                                line = lines[lineno++];

                                int idx;
                                while ((idx = line.find_first_of('D')) >= 0) {
                                    line.replace(idx, 1, "e");
                                }
                                while ((idx = line.find_first_of('d')) >= 0) {
                                    line.replace(idx, 1, "e");
                                }

                                // Must match primitives1; will work on the others later
                                if (!regex_match(line, what, primitives1))
                                    throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to match an exponent with one contraction:\n" + line);

                                if (!from_string<double>(exponent, what[1], std::dec))
                                    throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to convert exponent:\n" + line);
                                if (!from_string<double>(contraction, what[2], std::dec))
                                    throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to convert contraction:\n" + line);

                                // Scale the contraction
                                contraction *= scale;

                                // Save the information.
                                exponents[p] = exponent;
                                contractions[p] = contraction;
                            }

                            //                                printf("Adding new shell. nprimitive = %d\n", nprimitive);
                            // We have a full shell, push it to the basis set
                            shell_list.push_back(ShellInfo(am, contractions, exponents, gaussian_type, 0, center, 0, Unnormalized));
                        } else if (shell_type.size() == 2) {
                            // This is to handle instances of SP, PD, DF, FG, ...
                            int am1 = (int) shell_to_am[shell_type[0] - 'A'];
                            int am2 = (int) shell_to_am[shell_type[1] - 'A'];

                            std::vector<double> exponents(nprimitive);
                            std::vector<double> contractions1(nprimitive);
                            std::vector<double> contractions2(nprimitive);

                            for (int p = 0; p < nprimitive; ++p) {
                                line = lines[lineno++];

                                int idx;
                                while ((idx = line.find_first_of('D')) >= 0) {
                                    line.replace(idx, 1, "e");
                                }
                                while ((idx = line.find_first_of('d')) >= 0) {
                                    line.replace(idx, 1, "e");
                                }

                                // Must match primitivies2;
                                if (!regex_match(line, what, primitives2))
                                    throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to match an exponent with two contractions:\n" + line);

                                if (!from_string<double>(exponent, what[1], std::dec))
                                    throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to convert exponent:\n" + line);
                                if (!from_string<double>(contraction, what[2], std::dec))
                                    throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to convert first contraction:\n" + line);

                                // Scale the contraction
                                contraction *= scale;

                                // Save the information
                                exponents[p] = exponent;
                                contractions1[p] = contraction;

                                // Do the other contraction
                                if (!from_string<double>(contraction, what[3], std::dec))
                                    throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to convert second contraction:\n" + line);

                                // Scale the contraction
                                contraction *= scale;

                                // Save the information
                                contractions2[p] = contraction;
                            }

                            //                                printf("Adding 2 new shells. nprimitive = %d\n", nprimitive);
                            shell_list.push_back(ShellInfo(am1, contractions1, exponents, gaussian_type, 0, center, 0, Unnormalized));
                            shell_list.push_back(ShellInfo(am2, contractions2, exponents, gaussian_type, 0, center, 0, Unnormalized));
                        } else {
                            throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Unable to parse basis sets with spd, or higher grouping\n");
                        }
                    } else {
                        throw PSIEXCEPTION("Gaussian94BasisSetParser::parse: Expected shell information, but got:\n" + line);
                    }
                    line = lines[lineno++];
                }
                break;
            }
        }
    }
    if (found == false)
        throw BasisSetNotFound("Gaussian94BasisSetParser::parser: Unable to find the basis set for " + symbol + " in " + filename_, __FILE__, __LINE__);

    // The constructor, or the caller, should refresh the basis set.
    return shell_list;
}

} // namespace psi
