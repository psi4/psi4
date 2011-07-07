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

  $Id$
*/
#include <misc/misc.h>
#include <cstring>

namespace madness {
    std::istream& position_stream(std::istream& f, const std::string& tag) {
        std::string s;
        while (std::getline(f,s)) {
            std::string::size_type loc = s.find(tag, 0);
            if(loc != std::string::npos) return f;
        }
        std::string errmsg = std::string("position_stream: failed to locate ") + tag;
        MADNESS_EXCEPTION(errmsg.c_str(),0);
    }

    std::string lowercase(const std::string& s) {
        std::string r(s);
        for (unsigned int i=0; i<r.size(); ++i) r[i] = tolower(r[i]);
        return r;
    }

}

