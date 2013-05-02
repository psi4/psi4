/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <iostream>

#include "ref.h"

using namespace std;
using namespace smartptr;

#define concheck cout << "Constructor: " << __FILE__ << " " << __LINE__ << endl

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

void
boost::intrusive_ptr_add_ref(Countable* c)
{
    if (c == 0) return;

    unsigned int count = c->incref();
}

void
boost::intrusive_ptr_add_ref(const Countable* c)
{
    Countable* nonconst = const_cast<Countable*>(c);
    intrusive_ptr_add_ref(nonconst);
}

void
boost::intrusive_ptr_release(Countable* c)
{
    if (c == 0) return;
        

    unsigned int count = c->decref();

    if (count == 0)
    {
        delete c;
    }
    else if (count < 0)
    {
        cerr << "Pointer reference to " << c << " already zero.  Cannot release." << endl;
        abort();
    }
    else
    {
        //cerr << "Decrementing " << c << endl;
    }
}

void
boost::intrusive_ptr_release(const Countable* c)
{
    Countable* nonconst = const_cast<Countable*>(c);
    intrusive_ptr_release(nonconst);
}

Countable::Countable()
    : refcount_(0)
{
    //cerr << "Allocating " << this << endl;
}

unsigned int
Countable::incref()
{
    ++refcount_;
    return refcount_;
}

unsigned int
Countable::decref()
{
    --refcount_;
    return refcount_;
}

unsigned int
Countable::nref() const 
{
    return refcount_;
}

