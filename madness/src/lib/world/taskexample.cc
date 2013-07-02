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
#include <string>

using namespace madness;
using namespace std;

class Fred {
public:
    string a(const string& input) const {
        return input + string("a");
    }

    string b(const string& input) const {
        return input + string("b");
    }
};


int main(int argc, char** argv) {
    madness::initialize(argc,argv);
    madness::World world(MPI::COMM_WORLD);

    TaskAttributes attr;
    attr.set_stealable(true);

    Fred fred;
    
    Future<string> r = world.taskq.add(fred, &Fred::a, string("Hello"),attr);
    Future<string> s = world.taskq.add(fred, &Fred::b, r, attr);

    Future<string> fff;
    Future<string> ggg = world.taskq.add(fred, &Fred::a, fff, attr);

    world.taskq.steal(100);

    fff.set("done");

    print(s.get(), ggg.get());

    madness::finalize();
    return 0;
}
