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

//#define MADNESS_DISPLAY_EXCEPTION_BREAK_MESSAGE 0
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/worldref.h>
#include <world/world.h>
#include <world/worldobj.h>
#include <gtest/gtest.h>

madness::World* pworld;

using madness::RemoteReference;
using madness::detail::RemoteCounter;

namespace {
    class WorldRefTest : public ::testing::Test {
    public:
        WorldRefTest() : r(*pworld, i), r0() {
        }

        virtual ~WorldRefTest() {
        }

        virtual void SetUp() {
        }

        virtual void TearDown() {
        }

        static std::shared_ptr<int> i;

        RemoteReference<int> r;
        RemoteReference<int> r0;

        template <typename T>
        class XferRef : public madness::WorldObject<XferRef<T> > {
            madness::Void set_ptr(const RemoteReference<T>& r, bool away) {
                if(away) {
                    //std::cout << pworld->rank() << ": Set remote_ref = " << r << "\n";
                    remote_ref.set(r);
                } else {
                    //std::cout << pworld->rank() << ": Set return_ref = " << r << "\n";
                    return_ref.set(r);
                }
                return madness::None;
            }
        public:
            madness::Future<RemoteReference<T> > remote_ref;
            madness::Future<RemoteReference<T> > return_ref;

            explicit XferRef(madness::World& w) :
                    madness::WorldObject<XferRef<T> >(w)
            {
                madness::WorldObject<XferRef<T> >::process_pending();
            }

            void xfer(const ProcessID dest, const RemoteReference<T>& r, bool away) const {
                //std::cout << pworld->rank() << ": Sending " << r << " to " << dest << "\n";
                madness::WorldObject<XferRef<T> >::send(dest, & XferRef<T>::set_ptr,
                        r, away);
            }
        };
    };

    std::shared_ptr<int> WorldRefTest::i(new int(1));

    TEST_F(WorldRefTest, DefaultConstructor) {
        // Check that default constructed reference gives you nothing
        EXPECT_FALSE(r0);
        EXPECT_EQ(-1, r0.owner());
        EXPECT_FALSE(r0.is_local());
        EXPECT_EQ(0, r0.use_count());
        EXPECT_FALSE(r0.unique());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(r0.get(), madness::MadnessException);
        EXPECT_THROW(r0.get_world(), madness::MadnessException);
#endif
    }

    TEST_F(WorldRefTest, PointerConstructor) {
        // Check that the reference points to the correct data
        EXPECT_TRUE(r);
        EXPECT_EQ(pworld->rank(), r.owner());
        EXPECT_TRUE(r.is_local());
        EXPECT_EQ(1, r.use_count());
        EXPECT_TRUE(r.unique());
        EXPECT_EQ(i.get(), r.get());
        EXPECT_EQ(pworld, & (r.get_world()));

        RemoteReference<int> r2(*pworld, i);

        // Check that the reference points the same counter as r
        EXPECT_TRUE(r2);
        EXPECT_EQ(pworld->rank(), r2.owner());
        EXPECT_TRUE(r2.is_local());
        EXPECT_EQ(2, r2.use_count());
        EXPECT_FALSE(r2.unique());
        EXPECT_EQ(i.get(), r2.get());
        EXPECT_EQ(pworld, & (r2.get_world()));

        // Check that the use counter for r was also increased
        EXPECT_EQ(2, r.use_count());
        EXPECT_FALSE(r.unique());
    }

    TEST_F(WorldRefTest, CopyConstructor) {
        // check copy a non-null pointer
        RemoteReference<int> c(r);

        // Check that the reference points to the correct data
        EXPECT_EQ(r,                c);
        EXPECT_EQ(r.owner(),        c.owner());
        EXPECT_EQ(r.is_local(),     c.is_local());
        EXPECT_EQ(r.use_count(),    c.use_count());
        EXPECT_EQ(r.unique(),       c.unique());
        EXPECT_EQ(r.get(),          c.get());
        EXPECT_EQ(&(r.get_world()), &(c.get_world()));

        // Check copy a null reference
        RemoteReference<int> c0(r0);

        EXPECT_EQ(r0,               c0);
        EXPECT_EQ(r0.owner(),       c0.owner());
        EXPECT_EQ(r0.is_local(),    c0.is_local());
        EXPECT_EQ(r0.use_count(),   c0.use_count());
        EXPECT_EQ(r0.unique(),      c0.unique());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(c0.get(),          madness::MadnessException);
        EXPECT_THROW(c0.get_world(),    madness::MadnessException);
        EXPECT_THROW(r0.get(),          madness::MadnessException);
        EXPECT_THROW(r0.get_world(),    madness::MadnessException);
#endif
    }

    TEST_F(WorldRefTest, CopyConversionConstructor) {
        RemoteReference<const int> c(r);

        // Check that the reference points to the correct data
        EXPECT_EQ(r,                c);
        EXPECT_EQ(r.owner(),        c.owner());
        EXPECT_EQ(r.is_local(),     c.is_local());
        EXPECT_EQ(r.use_count(),    c.use_count());
        EXPECT_EQ(r.unique(),       c.unique());
        EXPECT_EQ(r.get(),          c.get());
        EXPECT_EQ(&(r.get_world()), &(c.get_world()));

        // Check copy a null reference
        RemoteReference<const int> c0(r0);

        EXPECT_EQ(r0,               c0);
        EXPECT_EQ(r0.owner(),       c0.owner());
        EXPECT_EQ(r0.is_local(),    c0.is_local());
        EXPECT_EQ(r0.use_count(),   c0.use_count());
        EXPECT_EQ(r0.unique(),      c0.unique());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(c0.get(),          madness::MadnessException);
        EXPECT_THROW(c0.get_world(),    madness::MadnessException);
        EXPECT_THROW(r0.get(),          madness::MadnessException);
        EXPECT_THROW(r0.get_world(),    madness::MadnessException);
#endif
    }

    TEST_F(WorldRefTest, CopyAssignment) {
        RemoteReference<int> c;

        // Check assignment operator
        c = r;

        EXPECT_EQ(r,                c);
        EXPECT_EQ(r.owner(),        c.owner());
        EXPECT_EQ(r.is_local(),     c.is_local());
        EXPECT_EQ(r.use_count(),    c.use_count());
        EXPECT_EQ(r.unique(),       c.unique());
        EXPECT_EQ(r.get(),          c.get());
        EXPECT_EQ(&(r.get_world()), &(c.get_world()));

        // Check assign a null reference
        c = r0;
        EXPECT_EQ(r0,               c);
        EXPECT_EQ(r0.owner(),       c.owner());
        EXPECT_EQ(r0.is_local(),    c.is_local());
        EXPECT_EQ(r0.use_count(),   c.use_count());
        EXPECT_EQ(r0.unique(),      c.unique());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(c.get(),           madness::MadnessException);
        EXPECT_THROW(c.get_world(),     madness::MadnessException);
        EXPECT_THROW(r0.get(),          madness::MadnessException);
        EXPECT_THROW(r0.get_world(),    madness::MadnessException);
#endif
    }

    TEST_F(WorldRefTest, CopyConversionAssignment) {
        RemoteReference<const int> c;

        // Check conversion assignment operator
        c = r;

        EXPECT_EQ(r,                c);
        EXPECT_EQ(r.owner(),        c.owner());
        EXPECT_EQ(r.is_local(),     c.is_local());
        EXPECT_EQ(r.use_count(),    c.use_count());
        EXPECT_EQ(r.unique(),       c.unique());
        EXPECT_EQ(r.get(),          c.get());
        EXPECT_EQ(&(r.get_world()), &(c.get_world()));

        // Check conversion assign a null reference
        c = r0;

        EXPECT_EQ(r0,               c);
        EXPECT_EQ(r0.owner(),       c.owner());
        EXPECT_EQ(r0.is_local(),    c.is_local());
        EXPECT_EQ(r0.use_count(),   c.use_count());
        EXPECT_EQ(r0.unique(),      c.unique());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(c.get(),           madness::MadnessException);
        EXPECT_THROW(c.get_world(),     madness::MadnessException);
        EXPECT_THROW(r0.get(),          madness::MadnessException);
        EXPECT_THROW(r0.get_world(),    madness::MadnessException);
#endif
    }

    TEST_F(WorldRefTest, Counter) {

        EXPECT_EQ(1, r.use_count());
        EXPECT_TRUE(r.unique());

        {
            // Check that the use count goes up when a copy is made
            RemoteReference<int> c(r);
            EXPECT_EQ(2, r.use_count());
            EXPECT_FALSE(r.unique());
        }

        // Check that the use counter goes down when the copy is destroyed
        EXPECT_EQ(1, r.use_count());
        EXPECT_TRUE(r.unique());

        {
            // Check that the use counter goes up when another pointer is created
            // with to the same pointer.
            RemoteReference<int> c(*pworld, i);
            EXPECT_EQ(2, r.use_count());
            EXPECT_FALSE(r.unique());
        }

        // Check that the use counter goes down when the other reference is destroyed
        EXPECT_EQ(1, r.use_count());
        EXPECT_TRUE(r.unique());

        r.reset();

        // Check that the use counter goes goes to zero when nothing is referenced
        EXPECT_EQ(0, r.use_count());
        EXPECT_FALSE(r.unique());
    }

    TEST_F(WorldRefTest, Pointer) {
        // Check pointer access
        EXPECT_EQ(i.get(), r.get());
        EXPECT_EQ(*i, *r);

        // Check that a null pointer throws correctly.
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(r0.get(), madness::MadnessException);
        EXPECT_THROW(*r0, madness::MadnessException);
#endif
    }

    TEST_F(WorldRefTest, Ownership) {
        // Check that r has an owner and that it is local
        EXPECT_EQ(pworld->rank(), r.owner());
        EXPECT_TRUE(r.is_local());

        // Check that r0 has no owner and is not local
        EXPECT_EQ(-1, r0.owner());
        EXPECT_FALSE(r0.is_local());
    }

    TEST_F(WorldRefTest, World) {
        // Check world accessor
        EXPECT_EQ(pworld, & (r.get_world()));

        // Check world accessor for null reference
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(r0.get_world(), madness::MadnessException);
#endif
    }

    TEST_F(WorldRefTest, Reset) {
        // check reset to a null pointer
        r.reset();

        EXPECT_FALSE(r);
        EXPECT_EQ(-1, r.owner());
        EXPECT_FALSE(r.is_local());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(r.get(), madness::MadnessException);
        EXPECT_THROW(r.get_world(), madness::MadnessException);
#endif
    }

    TEST_F(WorldRefTest, Boolean) {
        // Check  boolean conversion operator
        EXPECT_FALSE(r0);
        EXPECT_TRUE(r);
    }

    TEST_F(WorldRefTest, Swap) {
        RemoteReference<int> r1(r);

        // Verify starting conditions for remote reference
        EXPECT_TRUE(r);
        EXPECT_EQ(pworld->rank(), r.owner());
        EXPECT_TRUE(r.is_local());
        EXPECT_EQ(2, r.use_count());
        EXPECT_FALSE(r.unique());
        EXPECT_EQ(i.get(), r.get());
        EXPECT_EQ(pworld, & (r.get_world()));


        std::shared_ptr<int> ii(new int(2));
        RemoteReference<int> r2(*pworld, ii);

        // Swap r and r2 with member function
        r.swap(r2);

        // Check that r2 has r values
        EXPECT_TRUE(r2);
        EXPECT_EQ(pworld->rank(), r2.owner());
        EXPECT_TRUE(r2.is_local());
        EXPECT_EQ(2, r2.use_count());
        EXPECT_FALSE(r2.unique());
        EXPECT_EQ(i.get(), r2.get());
        EXPECT_EQ(pworld, & (r2.get_world()));

        // Check that r has r2 values
        EXPECT_TRUE(r);
        EXPECT_EQ(pworld->rank(), r.owner());
        EXPECT_TRUE(r.is_local());
        EXPECT_EQ(1, r.use_count());
        EXPECT_TRUE(r.unique());
        EXPECT_EQ(ii.get(), r.get());
        EXPECT_EQ(pworld, & (r.get_world()));

        // Swap back with free function
        madness::swap(r, r2);

        // Check that r is back to original settings
        EXPECT_TRUE(r);
        EXPECT_EQ(pworld->rank(), r.owner());
        EXPECT_TRUE(r.is_local());
        EXPECT_EQ(2, r.use_count());
        EXPECT_FALSE(r.unique());
        EXPECT_EQ(i.get(), r.get());
        EXPECT_EQ(pworld, & (r.get_world()));

        // Check that r2 is back to original settings
        EXPECT_TRUE(r2);
        EXPECT_EQ(pworld->rank(), r2.owner());
        EXPECT_TRUE(r2.is_local());
        EXPECT_EQ(1, r2.use_count());
        EXPECT_TRUE(r2.unique());
        EXPECT_EQ(ii.get(), r2.get());
        EXPECT_EQ(pworld, & (r2.get_world()));

        // Swap r with a null reference
        madness::swap(r, r0);

        // Check that r0 is set to r values
        EXPECT_TRUE(r0);
        EXPECT_EQ(pworld->rank(), r0.owner());
        EXPECT_TRUE(r0.is_local());
        EXPECT_EQ(2, r0.use_count());
        EXPECT_FALSE(r0.unique());
        EXPECT_EQ(i.get(), r0.get());
        EXPECT_EQ(pworld, & (r0.get_world()));

        // Check that r is set to r0 valuse
        EXPECT_FALSE(r);
        EXPECT_EQ(-1, r.owner());
        EXPECT_FALSE(r.is_local());
        EXPECT_EQ(0, r.use_count());
        EXPECT_FALSE(r.unique());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(r.get(), madness::MadnessException);
        EXPECT_THROW(r.get_world(), madness::MadnessException);
#endif
    }

    TEST_F(WorldRefTest, Serialize) {
        // Serialize 2 pointers to a buffer
        unsigned char buf[10*sizeof(RemoteReference<int>)];
        madness::archive::BufferOutputArchive oar(buf,sizeof(buf));
        oar & r & r0;
        std::size_t nbyte = oar.size();
        oar.close();

        // Check that the use count goes up after serialization
        EXPECT_EQ(2, r.use_count());

        // Release the last local references
        r.reset();

        // Deserialize 2 pointers from a buffer
        std::shared_ptr<int> ii(new int(2));
        RemoteReference<int> a;
        RemoteReference<int> a0(*pworld, ii);
        madness::archive::BufferInputArchive iar(buf,nbyte);
        iar & a & a0;
        iar.close();

        // check that serialized and deserialized pointers match the original.
        EXPECT_TRUE(a);
        EXPECT_EQ(pworld->rank(), a.owner());
        EXPECT_TRUE(a.is_local());
        EXPECT_EQ(1, a.use_count());
        EXPECT_TRUE(a.unique());
        EXPECT_EQ(i.get(), a.get());
        EXPECT_EQ(pworld, & (a.get_world()));

        EXPECT_FALSE(a0);
        EXPECT_EQ(-1, a0.owner());
        EXPECT_FALSE(a0.is_local());
        EXPECT_EQ(0, a0.use_count());
        EXPECT_FALSE(a0.unique());
#ifdef MADNESS_ASSERTIONS_THROW
        EXPECT_THROW(a0.get(), madness::MadnessException);
        EXPECT_THROW(a0.get_world(), madness::MadnessException);
#endif
    }

    TEST_F(WorldRefTest, RemoteRef) {

        if(pworld->size() > 1) {
            ProcessID send = (pworld->rank() != (pworld->size() - 1) ? pworld->rank() + 1 : 0);
            ProcessID recv = (pworld->rank() != 0 ? pworld->rank() - 1 : (pworld->size() - 1) );

            EXPECT_EQ(1, r.use_count());

            // Send the pointer to the next process
            XferRef<int> xfer_wobj(*pworld);
            xfer_wobj.xfer(send, r, true);

            pworld->gop.barrier();

            // Check that there is one remote reference in addition to the
            // single local reference
            EXPECT_EQ(2, r.use_count());

            pworld->gop.barrier();

            // wait for the remote reference
            RemoteReference<int>& remote = xfer_wobj.remote_ref.get();

            // Check for locality
            EXPECT_TRUE(remote);
            EXPECT_FALSE(remote.is_local());
            EXPECT_EQ(recv, remote.owner());
            EXPECT_EQ(pworld, &(remote.get_world()));
            EXPECT_EQ(0, remote.use_count());
            EXPECT_FALSE(remote.unique());

            // Check that you cannot remotely access reference
#ifdef MADNESS_ASSERTIONS_THROW
            EXPECT_THROW(remote.get(), madness::MadnessException);
            EXPECT_THROW(*remote, madness::MadnessException);
            // Should check arrow operator too but the assertions are the same as
            // operator* so we should get the same throw conditions
#endif

            xfer_wobj.xfer(recv, remote, false);

            pworld->gop.barrier();
            // Check that the remote reference is NULL after sending it back
            EXPECT_FALSE(remote);


            EXPECT_EQ(2, r.use_count());

            // wait for the original to return.
            RemoteReference<int>& back = xfer_wobj.return_ref.get();

            EXPECT_EQ(2, r.use_count());

            // Check that the reference was returned with correct info
            EXPECT_TRUE(back);
            EXPECT_TRUE(back.is_local());
            EXPECT_EQ(pworld->rank(), back.owner());
            EXPECT_EQ(2, back.use_count());
            EXPECT_FALSE(back.unique());
            EXPECT_EQ(pworld, &(back.get_world()));
            EXPECT_EQ(i.get(), back.get());
            EXPECT_EQ(*i, *back);
        }
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
