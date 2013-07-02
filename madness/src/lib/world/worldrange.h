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

  $Id: worldrange.h 2297 2011-05-04 12:45:23Z justus.c79@gmail.com $
*/
#ifndef MADNESS_WORLD_WORLDRANGE_H__INCLUDED
#define MADNESS_WORLD_WORLDRANGE_H__INCLUDED

#include <iterator>
#include <world/typestuff.h>

/// \file worldrange.h
/// \brief Implement Range class for parallel iteration

namespace madness {

    /// Dummy class a la Intel TBB used to distinguish splitting constructor
    class Split {};

    /// Range vaguely a la Intel TBB encapsulates random-access STL-like start and end iterators with chunksize
    template <typename iteratorT>
    class Range {
        long n;
        iteratorT start;
        iteratorT finish;
        int chunksize;
    public:
        typedef iteratorT iterator;

        /// Makes the range [start,finish)

        /// The motivated reader should look at the Intel TBB range,
        /// partitioner, split, concepts, etc..
        Range(const iterator& start, const iterator& finish, int chunk=1)
            : n(distance(start,finish))
            , start(start)
            , finish(finish)
            , chunksize(chunk)
        {
            if (chunksize < 1) chunksize = 1;
        }

        /// Copy constructor ... cost is O(1)
        Range(const Range& r)
                : n(r.n)
                , start(r.start)
                , finish(r.finish)
                , chunksize(r.chunksize)
        {}

        /// Splits range between new and old (r) objects ... cost is O(1)
        Range(Range& left, const Split& split)
                : n(0)
                , start(left.finish)
                , finish(left.finish)
                , chunksize(left.chunksize)
        {
            if (left.n > chunksize) {
                int nleft = (left.n+1)/2;

                start = left.start;
                advance(start,nleft);
                finish = left.finish;
                n = left.n - nleft;

                left.finish = start;
                left.n = nleft;
            }
        }

        /// Returns number of items in the range (cost is O(1))
        size_t size() const {
            return n;
        }

        /// Returns true if size=0
        bool empty() const {
            return n==0;
        }

        const iterator& begin() const {
            return start;
        }

        const iterator& end() const {
            return finish;
        }

        unsigned int get_chunksize() const {
            return chunksize;
        }

    private:
        template<typename integralT, typename distanceT>
        inline static typename enable_if<std::is_integral<integralT>, void>::type
        advance(integralT& i, distanceT n) { i += n; }

        template<typename iterT, typename distanceT>
        inline static typename disable_if<std::is_integral<iterT>, void>::type
        advance(iterT& it, distanceT n) { std::advance(it, n); }

        template<class integralT>
        inline static typename enable_if<std::is_integral<integralT>, integralT>::type
        distance(integralT first, integralT last) { return last - first; }

        template <class iterT>
        struct diff_type {
            typedef typename std::iterator_traits<iterT>::difference_type type;
        };

        template<class iterT>
        inline static typename lazy_disable_if<std::is_integral<iterT>, diff_type<iterT> >::type
        distance(iterT first, iterT last) { return std::distance(first, last); }
    };

} // namespace madness

#endif // MADNESS_WORLD_WORLDRANGE_H__INCLUDED
