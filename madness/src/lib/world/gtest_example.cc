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

#include <madness_config.h>

#ifdef MADNESS_HAS_GOOGLE_TEST

#include <world/world.h>
#include <gtest/gtest.h>

madness::World* pworld;

namespace {
    class WorldTest : public ::testing::Test {
    public:
        WorldTest() {
            // You can do set-up work for each test here.
            madness::print("Inside WorldTest");
        }

        virtual ~WorldTest() {
            // You can do clean-up work that doesn't throw exceptions here.
            madness::print("Inside ~WorldTest");
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:
        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
            madness::print("Inside WorldSetUp");
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
            madness::print("Inside WorldTestDown");
        }

        // Objects declared here can be used by all tests in the test case for World.
        static madness::World* pWorld;
    };

    TEST_F(WorldTest, Something) {
        EXPECT_EQ(1,1);
    }

    TEST_F(WorldTest, Else) {
        EXPECT_EQ(1,2);
    }
}

int main(int argc, char **argv) {
    madness::initialize(argc,argv);
    madness::World world(MPI::COMM_WORLD);
    pworld = &world;

    if (world.rank()) madness::redirectio(world);
    world.args(argc,argv);
    world.gop.fence();

    ::testing::InitGoogleTest(&argc, argv);
    int status = RUN_ALL_TESTS();

    world.gop.fence();

    madness::finalize();
    return status;
}


#else

#include <iostream>
int main() {
    std::cout << "U need to build with Google test to enable the world test code\n";
    return 0;
}

#endif
