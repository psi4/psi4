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

#include "class.h"

using namespace yeti;
using namespace std;

const char* TypeInfo<double>::printf_str = "%12.8f";
const char* TypeInfo<float>::printf_str = "%12.8f";
const char* TypeInfo<quad>::printf_str = "%20.14Lf";
const char* TypeInfo<int>::printf_str = "%6d";

float TestEquals<float>::cutoff = 1e-12;
double TestEquals<double>::cutoff = 1e-12;
quad TestEquals<quad>::cutoff = 1e-12;

