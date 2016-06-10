/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#ifndef COMPILER_H
#define COMPILER_H

#include <psiconfig.h>

/*
 * Example of likely and unlikely
 *
 * Listing #5 from http://www.ibm.com/developerworks/linux/library/l-gcc-hacks/index.html
 *
 *  unsigned int __skb_checksum_complete(struct sk_buff *skb)
 *  {
 *      unsigned int sum;
 *
 *      sum = (u16)csum_fold(skb_checksum(skb, 0, skb->len, skb->csum));
 *      if (likely(!sum)) {
 *              if (unlikely(skb->ip_summed == CHECKSUM_HW))
 *                      netdev_rx_csum_fault(skb->dev);
 *              skb->ip_summed = CHECKSUM_UNNECESSARY;
 *      }
 *      return sum;
 *  }
 */

#if defined(HAVE_BUILTIN_EXPECT) && defined(HAVE_BUILTIN_CONSTANT_P)
#   define likely(x)   (__builtin_constant_p(x) ? !!(x) : __builtin_expect(!!(x), 1))
#   define unlikely(x) (__builtin_constant_p(x) ? !!(x) : __builtin_expect(!!(x), 0))
#elif defined(HAVE_BUILTIN_EXPECT)
#   define likely(x)   (__builtin_expect(!!(x), 1))
#   define unlikely(x) (__builtin_expect(!!(x), 0))
#else
#   define likely(x)   (x)
#   define unlikely(x) (x)
#endif

/*
 * Example of prefetching:
 *
 * Listing #6 from http://www.ibm.com/developerworks/linux/library/l-gcc-hacks/index.html
 *
 * static inline void prefetch_range(void *addr, size_t len)
 * {
 *    char *cp;
 *    char *end = addr + len;
 *
 *    for (cp = addr; cp < end; cp += PREFETCH_STRIDE)
 *        prefetch(cp);
 * }
 */

#ifdef HAVE_BUILTIN_PREFETCH
    /// Default prefetch; high degree of temporal locality with read access
#   define prefetch_r(x) __builtin_prefetch((x))
    /// High degree of temporal locality with read/write access
#   define prefetch_rw(x) __builtin_prefetch((x), 1)
    /// No temporal locality with (need not be kept in cache after access
#   define prefetch_r_nokeep(x) __builtin_prefetch((x), 0, 0)
#   define prefetch_rw_nokeep(x) __buitin_prefetch((x), 1, 0)
#else
#   define prefetch_r(x)
#   define prefetch_rw(x)
#   define prefetch_r_nokeep(x)
#   define prefetch_rw_nokeep(x)
#endif

#endif // COMPILER_H

