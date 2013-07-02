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


*/

#ifndef MADNESS_WORLD_POSIXMEM_H__INCLUDED
#define MADNESS_WORLD_POSIXMEM_H__INCLUDED

/// \file world/posixmem.h
/// \brief Implement dummy posix_memalign if it is missing on the system.

#include <madness_config.h>

#if !HAVE_POSIX_MEMALIGN
#include <sys/errno.h>
static inline int posix_memalign(void **memptr, std::size_t alignment, std::size_t size) {
    *memptr=malloc(size);
    if (*memptr) return 0;
    else return ENOMEM;
}
#elif MISSING_POSIX_MEMALIGN_PROTO
extern "C"  int posix_memalign(void **memptr, std::size_t alignment, std::size_t size);
#endif

#endif // MADNESS_WORLD_POSIXMEM_H__INCLUDED
