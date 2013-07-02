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


  $Id: listbench.cc 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#include <list>
#include <iostream>
#include <time.h>
#include <sys/time.h>

using namespace std;

double cpu_time() {
    return double(clock())/CLOCKS_PER_SEC;
}


int main() {
    double used;
    int** p = new int*[1000000];

    used = cpu_time();
    for (int i=0; i<1000000; ++i) p[i] = new int(i);
    used = (cpu_time() - used);
    cout << "time to new 1M integers " << used << endl;

    used = cpu_time();
    for (int i=0; i<1000000; ++i) delete p[i];
    used = (cpu_time() - used);
    cout << "time to del 1M integers " << used << endl;


    list<int> a;

    used = cpu_time();
    for (int i=0; i<1000000; ++i) {
        a.push_back(i);
    }
    used = cpu_time() - used;
    cout << "time to push 1M integers in list " << used << endl;

    used = cpu_time();
    for (int i=0; i<1000000; ++i) {
        a.pop_front();
    }
    used = cpu_time() - used;
    cout << "time to pop  1M integers in list " << used << endl;


    list<int*> b;

    used = cpu_time();
    for (int i=0; i<1000000; ++i) {
        b.push_back(new int(i));
    }
    used = cpu_time() - used;
    cout << "time to new+push 1M integers in list " << used << endl;

    used = cpu_time();
    for (int i=0; i<1000000; ++i) {
        delete b.front();
        b.pop_front();
    }
    used = cpu_time() - used;
    cout << "time to del+pop  1M integers in list " << used << endl;

    return 0;
}


