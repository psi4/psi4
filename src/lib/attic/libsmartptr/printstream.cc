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

#include <stdarg.h>

#include "printstream.h"

using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

string
std::stream_printf(const char* fmt, ...)
{
    char str[1024];

    va_list args;
    va_start(args, fmt);
    vsprintf(str, fmt, args);
    va_end(args);

    string strobj = str;
    
    return strobj;
}

void
std::print_map(
    map<string, string>& keymap,
    ostream& os
)
{
    map<string, string>::iterator it;
    for (it = keymap.begin(); it != keymap.end(); ++it)
        os << stream_printf("%20s : %s", it->first.data(), it->second.data()) << endl;
}
