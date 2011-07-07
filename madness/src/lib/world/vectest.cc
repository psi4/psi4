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


  $Id: vectest.cc 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#include <iostream>
#include <vector>
using namespace std;

class A {
    int val;
public:
    A() : val(0) {
        cout << "Default constructor" << endl;
    };

    A(int val) : val(val) {
        cout << "Value constructor " << val << endl;
    };

    A(const A& a) : val(a.val) {
        cout << "Copy constructor " << val << endl;
    };

    A& operator=(const A& a) {
        cout << "Assignment " << val << " " << a.val << endl;
        val = a.val;
        return *this;
    };

    void set(int a) {
        val = a;
    };

    ~A() {
        cout << "Destructor" << endl;
    }
};

int main() {
    cout << "Making vector(3)" << endl;
    vector<A> v(3);
    cout << "Finished making vector" << endl;

    cout << "Assigning vector values" << endl;
    for (int i=0; i<(int)v.size(); ++i) v[i].set(i+1);
    cout << "Finished assigning vector values" << endl;

    cout << "Making empty vector" << endl;
    vector<A> u;
    cout << "Finished making empty vector" << endl;

    cout << "Assigning to empty vector from vector(3)" << endl;
    u = v;
    cout << "Finished assigning to empty vector" << endl;

    cout << "Reassigning vector values" << endl;
    for (int i=0; i<(int)v.size(); ++i) v[i].set(i+5);
    cout << "Finished reassigning vector values" << endl;

    cout << "Assigning to existing vector from vector(3)" << endl;
    u = v;
    cout << "Finished assigning to existing vector" << endl;

    cout << "Vector copy constructor from vector(3)" << endl;
    vector<A> p(v);
    cout << "Finished vector copy constructor from vector(3)" << endl;

    return 0;
}






