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
#include <world/world.h>

using namespace madness;


void dave(int i) {
    PROFILE_FUNC;
    if (i&1) {
        PROFILE_BLOCK(mary);
    }
    if (i&3) {
        PROFILE_BLOCK(calvin);
    }
}

void fred(int i) {
    PROFILE_FUNC;
    dave(i);
}

void realmain(int argc, char** argv)
{
    World world(MPI::COMM_WORLD);
    for (int i=0; i<1000; ++i)
        fred(i);

    WorldProfile::print(world);
}

int main(int argc, char** argv) {
    MPI::Init(argc, argv);

    //This programming style guarantees that
    //the world instance is destroyed before
    //MPI::Finalize().
    //It is necessary if the destructor of the
    //class World also takes care of MPI status.
    realmain(argc, argv);

    MPI::Finalize();
    return 0;
}
