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

#ifndef smartptr_printstream_h_
#define smartptr_printstream_h_

#include <sstream>
#include <stdio.h>
#include <iostream>
#include <map>

namespace std {

/**
    @param fmt A printf format string
    @param va_args The set of arguments
    @return A string with formatted values
*/
std::string 
stream_printf(const char* fmt, ...);

template <
  class Type
> void
print_keys(
    std::map<std::string, Type>& keymap,
    std::stringstream& sstr
)
{
    typename std::map<std::string, Type>::iterator it;
    for (it = keymap.begin(); it != keymap.end(); it++)
        sstr << it->first << std::endl;
}

/**
 @param keymap A string map
 @param os The stream to print the map to
*/
void
print_map(
    std::map<std::string, std::string>& keymap,
    std::ostream& os
);

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif
