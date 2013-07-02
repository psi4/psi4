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


  $Id: hang.cc 1602 2009-12-27 19:53:06Z rjharrison $
*/


#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    int np = MPI::COMM_WORLD.Get_size();
    if (np != 2) MADNESS_EXCEPTION("2 only", np);

    int me = MPI::COMM_WORLD.Get_rank();
    int other = me? 0 : 1;

    int a=0, b=-1;
    MPI::Request rsend = MPI::COMM_WORLD.Isend(&a, sizeof(a), MPI::BYTE, other, 1);
    MPI::Request rrecv = MPI::COMM_WORLD.Irecv(&b, sizeof(b), MPI::BYTE, other, 1);

    MPI::Status status;

    while (!rsend.Get_status(status));
    while (!rrecv.Get_status(status));
    rsend.Test(status);
    rrecv.Test(status);

    //while (!rsend.Test(status)) ;
    //while (!rrecv.Test(status)) ;

    cout << me << " got " << b << endl;

    MPI::Finalize();
    return 0;
}
