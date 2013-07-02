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


  $Id: worldmpi.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#ifndef MADNESS_WORLD_WORLDMPI_H__INCLUDED
#define MADNESS_WORLD_WORLDMPI_H__INCLUDED

/// \file worldmpi.h
/// \brief Implements WorldMpiInterface

/*
// If include mpi.h BEFORE stdio/iostream should not need undefs
#ifdef SEEK_CUR
#undef SEEK_CUR
#endif
#ifdef SEEK_SET
#undef SEEK_SET
#endif
#ifdef SEEK_END
#undef SEEK_END
#endif
*/

#include <world/safempi.h>
#include <world/worldtypes.h>

namespace madness {

    static const Tag DYNAMIC_TAG_BASE = 1024;

    class WorldAmInterface;
    class WorldGopInterface;

    /// This class wraps/extends the MPI interface for World
    class WorldMpiInterface : public SafeMPI::Intracomm {
    public:
        WorldMpiInterface(MPI::Intracomm& comm) : SafeMPI::Intracomm(comm) {}

        /// Returns the associated SafeMPI communicator
        SafeMPI::Intracomm& comm() {
            return *static_cast<SafeMPI::Intracomm*>(this);
        }
    };

}

#endif // MADNESS_WORLD_WORLDMPI_H__INCLUDED
