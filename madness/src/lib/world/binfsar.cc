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


  $Id: $
*/

#include <world/binfsar.h>
#include <world/worldexc.h>
#include <cstring>

namespace madness {
    namespace archive {

        BinaryFstreamOutputArchive::BinaryFstreamOutputArchive(const char* filename, std::ios_base::openmode mode)
                : iobuf()
        {
            if (filename) open(filename, mode);
        }

        void BinaryFstreamOutputArchive::open(const char* filename, std::ios_base::openmode mode) {
            iobuf.reset(new char[IOBUFSIZE], &detail::checked_array_delete<char>);
            os.open(filename, mode);
#ifndef ON_A_MAC
            os.rdbuf()->pubsetbuf(iobuf.get(), IOBUFSIZE);
#endif

            store(ARCHIVE_COOKIE, strlen(ARCHIVE_COOKIE)+1);
        }

        void BinaryFstreamOutputArchive::close() {
            if (iobuf) {
                os.close();
                iobuf.reset();
            }
        };

        void BinaryFstreamOutputArchive::flush() {
            os.flush();
        }

        BinaryFstreamInputArchive::BinaryFstreamInputArchive(const char* filename, std::ios_base::openmode mode)
                : iobuf() {
            if (filename) open(filename, mode);
        }


        void BinaryFstreamInputArchive::open(const char* filename,  std::ios_base::openmode mode) {
            iobuf.reset(new char[IOBUFSIZE], &detail::checked_array_delete<char>);
            is.open(filename, mode);
            if (!is) MADNESS_EXCEPTION("BinaryFstreamInputArchive: open: failed", 1);
            is.rdbuf()->pubsetbuf(iobuf.get(), IOBUFSIZE);
            char cookie[255];
            int n = strlen(ARCHIVE_COOKIE)+1;
            load(cookie, n);
            if (strncmp(cookie,ARCHIVE_COOKIE,n) != 0)
                MADNESS_EXCEPTION("BinaryFstreamInputArchive: open: not an archive?", 1);
        }

        void BinaryFstreamInputArchive::close() {
            if (iobuf) {
                is.close();
                iobuf.reset();
            }
        }

    } // namespace archive
} // namespace madness
