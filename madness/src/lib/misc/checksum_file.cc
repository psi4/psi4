/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id: checksum_file.cc 1602 2009-12-27 19:53:06Z rjharrison $
*/


#include <cstdio>
using std::fopen;
using std::fgetc;
using std::fclose;

/// \file checksum_file.cc
/// \brief Miscellaneous useful stuffx

namespace madness {

    /// Simple checksum for ASCII characters in file
    unsigned  long checksum_file(const char* filename) {
        FILE *file = fopen(filename,"r");
        if (!file) return 0;

        unsigned long sum = 0;
        int c;
        while ((c = fgetc(file)) != EOF) {
            sum = (sum*31u + ((unsigned) c)) & 0xffffff;
        }
        fclose(file);
        return sum;
    }
}
