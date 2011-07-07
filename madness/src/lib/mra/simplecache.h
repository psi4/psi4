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

  $Id: simplecache.h 2244 2011-03-30 15:25:31Z justus.c79@gmail.com $
*/
#ifndef MADNESS_MRA_SIMPLECACHE_H__INCLUDED
#define MADNESS_MRA_SIMPLECACHE_H__INCLUDED

#include <mra/mra.h>

namespace madness {
    /// Simplified interface around hash_map to cache stuff for 1D

    /// This is a write once cache --- subsequent writes of elements
    /// have no effect (so that pointers/references to cached data
    /// cannot be invalidated)
    template <typename Q, std::size_t NDIM>
    class SimpleCache {
    private:
        typedef ConcurrentHashMap< Key<NDIM>, Q > mapT;
        typedef std::pair<Key<NDIM>, Q> pairT;
        mapT cache;

    public:
        SimpleCache() : cache() {};

        SimpleCache(const SimpleCache& c) : cache(c.cache) {};

        SimpleCache& operator=(const SimpleCache& c) {
            if (this != &c) {
                cache.clear();
                cache = c.cache;
            }
            return *this;
        }

        /// If key is present return pointer to cached value, otherwise return NULL
        inline const Q* getptr(const Key<NDIM>& key) const {
            typename mapT::const_iterator test = cache.find(key);
            if (test == cache.end()) return 0;
            return &(test->second);
        }


        /// If key=(n,l) is present return pointer to cached value, otherwise return NULL

        /// This for the convenience (backward compatibility) of 1D routines
        inline const Q* getptr(Level n, Translation l) const {
            Key<NDIM> key(n,Vector<Translation,NDIM>(l));
            return getptr(key);
        }


        /// If key=(n,l) is present return pointer to cached value, otherwise return NULL

        /// This for the convenience (backward compatibility) of 1D routines
        inline const Q* getptr(Level n, const Key<NDIM>& disp) const {
            Key<NDIM> key(n,disp.translation());
            return getptr(key);
        }


        /// Set value associated with key ... gives ownership of a new copy to the container
        inline void set(const Key<NDIM>& key, const Q& val) {
            cache.insert(pairT(key,val));
        }

        inline void set(Level n, Translation l, const Q& val) {
            Key<NDIM> key(n,Vector<Translation,NDIM>(l));
            set(key, val);
        }

        inline void set(Level n, const Key<NDIM>& disp, const Q& val) {
            Key<NDIM> key(n,disp.translation());
            set(key, val);
        }
    };
}
#endif // MADNESS_MRA_SIMPLECACHE_H__INCLUDED
