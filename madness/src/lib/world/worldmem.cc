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


  $Id: world.cc 454 2008-01-26 03:08:15Z rjharrison $
*/

#include <world/worldmem.h>
#include <cstdlib>
//#include <cstdio>
#include <limits.h>
#include <iostream>
#include <iomanip>

/*

  Memory allocation structure

  ---------
  base             <--- actual-pointer
  ...
  user-size        <--- p - sizeof(3*long)
  actual-pointer   <--- p - sizeof(2*long)
  pre-checksum     <--- p - sizeof(long)
  ---------
  user-size bytes  <--- p
  16-byte aligned
  ---------
  post-checksum
  ---------

 */


static madness::WorldMemInfo stats = {0, 0, 0, 0, 0, 0, ULONG_MAX, false};

namespace madness {
    WorldMemInfo* world_mem_info() {
        return &stats;
    }

    void WorldMemInfo::do_new(void *p, std::size_t size) {
        ++num_new_calls;
        ++cur_num_frags;
        if (cur_num_frags > max_num_frags) max_num_frags = cur_num_frags;
        cur_num_bytes += size;
        if (cur_num_bytes > max_num_bytes) max_num_bytes = cur_num_bytes;

        if (trace)
            std::cout << "WorldMemInfo: allocating " << p << " " << size << "\n";
    }

    void WorldMemInfo::do_del(void *p, std::size_t size) {
        ++num_del_calls;
        --cur_num_frags;
        cur_num_bytes -= size;

        if (trace)
            std::cout << "WorldMemInfo: deleting " << p << " " << size << "\n";
    }

    void WorldMemInfo::print() const {
        std::cout.flush();
        std::cout << "\n    MADNESS memory statistics\n";
        std::cout << "    -------------------------\n";
        std::cout << "      overhead bytes per frag " << std::setw(12)
            << overhead << "\n";
        std::cout << "         calls to new and del " << std::setw(12)
            << num_new_calls << " " << std::setw(12) << num_del_calls << "\n";
        std::cout << "  cur and max frags allocated " << std::setw(12)
            << cur_num_frags << " " << std::setw(12) << max_num_frags << "\n";
        std::cout << "  cur and max bytes allocated " << std::setw(12)
            << cur_num_bytes << " " << std::setw(12) << max_num_bytes << "\n";
    }

    void WorldMemInfo::reset() {
        num_new_calls = 0;
        num_del_calls = 0;
        cur_num_frags = 0;
        max_num_frags = 0;
        cur_num_bytes = 0;
        max_num_bytes = 0;
    }

}

#ifdef WORLD_GATHER_MEM_STATS

#include <exception>

static const unsigned long pre_checksum = 0xdeadbeef;
static const unsigned char post_checksum = 0xab;

namespace madness {
  // This is called from the delete operator when a buffer underflow is detected.
  // Any error code for this condition should be added here.
  void world_mem_buffer_underflow() {
      std::cerr << "WorldMemInfo: Buffer underflow detected.\n" <<
          "Set a breakpoint at madness::world_mem_buffer_underflow() to debug.\n";
  }

  // This is called from the delete operator when a buffer overflow is detected.
  // Any error code for this condition should be added here.
  void world_mem_buffer_overflow() {
      std::cerr << "WorldMemInfo: Buffer overflow detected.\n" <<
          "Set a breakpoint at madness::world_mem_buffer_overflow() to debug.\n";
  }
} // namespace madness

void* operator new(size_t size) throw (std::bad_alloc) {
    // user-size + actual_pointer + pre_checksum + post_checksum + padding
    std::size_t actual_size = size + madness::WorldMemInfo::overhead;

    if (size+stats.cur_num_bytes > stats.max_mem_limit)
        throw std::bad_alloc(); //MADNESS_EXCEPTION("WorldMemInfo: execeeded max_mem_limit", stats.max_mem_limit);

    unsigned long* actual_pointer= (unsigned long*) malloc(actual_size);
    if (actual_pointer==0) throw std::bad_alloc(); // ANSI/ISO compliant behavior
    unsigned char* p = (unsigned char*)(actual_pointer + 3);

    p += 16 - (((unsigned long) p)&0x0f);
    p[size] = post_checksum;

    unsigned long* lp = (unsigned long*) p;  // Non-portable alignment assumption here????

    lp[-1] = pre_checksum;
    lp[-2] = (unsigned long) actual_pointer;
    lp[-3] = size;

    stats.do_new(p, size);

    return (void *) p;
}

void operator delete(void *p) throw() {
    unsigned long* lp = (unsigned long*) p;
    if (lp[-1] != pre_checksum) madness::world_mem_buffer_underflow();
    unsigned long* actual_pointer = (unsigned long*) lp[-2];
    std::size_t size = lp[-3];
    unsigned char* cp = (unsigned char*) p;
    if (cp[size] != post_checksum) madness::world_mem_buffer_overflow();

    stats.do_del(p, size);

    free(actual_pointer);
}

#endif
