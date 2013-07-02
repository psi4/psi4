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
#include <world/array.h>
#include <madness_config.h>

#ifdef MADNESS_HAS_GOOGLE_TEST

#include <gtest/gtest.h>

namespace {
    TEST(VectorTest, SizeWorks) {
        madness::Vector<double,3> v;
        EXPECT_EQ(3u, v.size());
    }

    TEST(VectorTest, InitializerWorks) {
        madness::Vector<double,33> v(1.0);
        for (std::size_t i=0; i<v.size(); ++i) {
            EXPECT_EQ(1.0, v[i]);
        }
    }

    TEST(VectorTest, Hash) {
        madness::Vector<double,33> v(1.0);
        EXPECT_NE(0u, v.hash());
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#else

#include <iostream>
int main() {
    std::cout << "U need to build with Google test to enable the array test code\n";
    return 0;
}

#endif
