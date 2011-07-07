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
#include <world/safempi.h>
//#include <world/worldthread.h>
#include <world/worldexc.h>

namespace SafeMPI {

#ifdef SERIALIZE_MPI
    madness::SCALABLE_MUTEX_TYPE charon;
#endif
    
    void Intracomm::binary_tree_info(int root, int& parent, int& child0, int& child1) {
        int np = nproc();
        int me = (rank()+np-root)%np;   // Renumber processes so root has me=0
        parent = (((me-1)>>1)+root)%np;        // Parent in binary tree
        child0 = (me<<1)+1+root;        // Left child
        child1 = (me<<1)+2+root;        // Right child
        if (child0 >= np && child0<(np+root)) child0 -= np;
        if (child1 >= np && child1<(np+root)) child1 -= np;
        
        if (me == 0) parent = -1;
        if (child0 >= np) child0 = -1;
        if (child1 >= np) child1 = -1;
    }


} //namespace SafeMPI
