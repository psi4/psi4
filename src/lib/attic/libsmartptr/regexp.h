/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef smartptr_regexp_h
#define smartptr_regexp_h

#include <string>
#include <vector>

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace smartptr {

enum RegexpFlags { IgnoreCase=1, UpperCase=2, FindAll=4, LowerCase=8, IncludeLineBreaks=16, StripWhitespace=32 };

/**
 * Extracts a double precision number from text.
 * @param regexp The regular expression containing the (captured) number.
 * @param text The text to apply the RegEx to.
 * @param flags An integer comprising bitwise options (see the PyFlags enum).
 * @return The value as a double.
 */
double get_regexp_double(const std::string &regexp, const std::string &text, int flags = 0);

/**
 * Extracts a string from text.
 * @param regexp The regular expression containing the (captured) string.
 * @param text The text to apply the RegEx to.
 * @param flags An integer comprising bitwise options (see the PyFlags enum).
 * @return The string.
 */
std::string get_regexp_string(const std::string &regexp, const std::string &text, int flags = 0);

/**
 * Extracts an array of doubles from text, according to a given RegEx.
 * Memory is allocated by the function and ownership is passed to calling routine.
 * @param regexp The regular expression containing the (captured) values.
 * @param text The text to apply the RegEx to.
 * @param flags An integer comprising bitwise options (see the PyFlags enum).
 * @return an array of all captured values.
 */
double* get_regexp_double_array(const std::string& regexp, const std::string &text, size_t& length, int flags = 0);

/**
 * Extracts strings from text, according to a given RegEx.
 * @param matches The vector of strings to hold the results.
 * @param regexp The regular expression containing the (captured) strings.
 * @param text The text to apply the RegEx to.
 * @param flags An integer comprising bitwise options (see the PyFlags enum).
 */
void
findmatch(std::vector<std::string>& matches, const std::string &regexp, const std::string &text, int flags = 0);

/**
 * Tests a string to see if it contains a certain RegEx.
 * @param regexp The regular expression to test.
 * @param text The text to apply the RegEx to.
 * @param flags An integer comprising bitwise options (see the PyFlags enum).
 * @return whether the text contains a match.
 */
bool has_regexp_match(const std::string& regexp, const std::string& text, int flags = 0);

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif
