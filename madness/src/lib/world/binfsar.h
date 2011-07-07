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


  $Id: binfsar.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#ifndef MADNESS_WORLD_BINFSAR_H__INCLUDED
#define MADNESS_WORLD_BINFSAR_H__INCLUDED

/// \file binfsar.h
/// \brief Implements archive wrapping a binary filestream

#include <fstream>
#include <world/archive.h>
#include <world/sharedptr.h>


namespace madness {
    namespace archive {

        /// Wraps an archive around a binary file stream for output
        class BinaryFstreamOutputArchive : public BaseOutputArchive {
            static const std::size_t IOBUFSIZE = 4*1024*1024;
            std::shared_ptr<char> iobuf;
            mutable std::ofstream os;
        public:
            BinaryFstreamOutputArchive(const char* filename = 0,
                                       std::ios_base::openmode mode = std::ios_base::binary | \
                                                                      std::ios_base::out | std::ios_base::trunc);

            template <class T>
            inline
            typename madness::enable_if< madness::is_serializable<T>, void >::type
            store(const T* t, long n) const {
                os.write((const char *) t, n*sizeof(T));
            }

            void open(const char* filename,
                      std::ios_base::openmode mode = std::ios_base::binary | \
                                                     std::ios_base::out |  std::ios_base::trunc);

            void close();

            void flush();
        };


        /// Wraps an archive around a binary file stream for input
        class BinaryFstreamInputArchive : public BaseInputArchive {
            static const std::size_t IOBUFSIZE = 4*1024*1024;
            std::shared_ptr<char> iobuf;
            mutable std::ifstream is;
        public:
            BinaryFstreamInputArchive(const char* filename = 0, std::ios_base::openmode mode = std::ios_base::binary | std::ios_base::in);

            template <class T>
            inline
            typename madness::enable_if< madness::is_serializable<T>, void >::type
            load(T* t, long n) const {
                is.read((char *) t, n*sizeof(T));
            }

            void open(const char* filename,  std::ios_base::openmode mode = std::ios_base::binary | std::ios_base::in);

            void close();
        };
    }
}
#endif // MADNESS_WORLD_BINFSAR_H__INCLUDED
