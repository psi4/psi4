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


  $Id: print.h 2207 2011-03-04 16:32:32Z justus.c79@gmail.com $
*/


#ifndef MADNESS_WORLD_PRINT_H__INCLUDED
#define MADNESS_WORLD_PRINT_H__INCLUDED

/// \file print.h
/// \brief Defines simple templates for printing to std::cout "a la Python"


#include <iostream>
#include <complex>
#include <list>
#include <vector>
#include <world/typestuff.h>
#include <world/enable_if.h>

#ifdef BRAINDEAD
// Cray XT nonsense
#define ENDL "\n"

#else

#define ENDL std::endl

#endif


namespace madness {

    /// Easy printing of complex numbers
    template <typename T>
    std::ostream& operator<<(std::ostream& s, const std::complex<T>& c) {
        s << c.real() << "+" << c.imag() << "j";
        return s;
    }

    /// Easy printing of pairs
    template <typename T, typename U>
    std::ostream& operator<<(std::ostream& s, const std::pair<T,U>& p) {
        s << "(" << p.first << "," << p.second << ")";
        return s;
    }

    /// Easy printing of lists
    template <typename T>
    std::ostream& operator<<(std::ostream& s, const std::list<T>& c) {
        s << "[";
        typename std::list<T>::const_iterator it = c.begin();
        while (it != c.end()) {
            s << *it;
            ++it;
            if (it != c.end()) s << ", ";
        };
        s << "]";
        return s;
    }

    /// Easy printing of vectors
    template <typename T>
    std::ostream& operator<<(std::ostream& s, const std::vector<T>& c) {
        s << "[";
        typename std::vector<T>::const_iterator it = c.begin();
        while (it != c.end()) {
            s << *it;
            ++it;
            if (it != c.end()) s << ", ";
        };
        s << "]";
        return s;
    }

    /// Easy printing of fixed dimension arrays

    /// STL I/O already does char.
    template <typename T, std::size_t N>
    typename disable_if<std::is_same<T,char>, std::ostream&>::type
    operator<<(std::ostream& s, const T(&v)[N]) {
        s << "[";
        for (std::size_t i=0; i<N; ++i) {
            s << v[i];
            if (i != (N-1)) s << ",";
        }
        s << "]";
        return s;
    }

    /// Print a single item to std::cout terminating with new line
    template <class A>
    void print(const A& a) {
        std::cout << a << ENDL;
    }

    /// Print two items separated by spaces to std::cout terminating with new line
    template <class A, class B>
    void print(const A& a, const B& b) {
        std::cout << a << " " << b << ENDL;
    }

    /// Print three items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C>
    void print(const A& a, const B& b, const C& c) {
        std::cout << a << " " << b << " " << c << ENDL;
    }

    /// Print four items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D>
    void print(const A& a, const B& b, const C& c, const D& d) {
        std::cout << a << " " << b << " " << c << " " << d << ENDL;
    }

    /// Print five items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << ENDL;
    }

    /// Print six items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E, class F>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << ENDL;
    }

    /// Print seven items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << ENDL;
    }

    /// Print eight items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G, class H>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << ENDL;
    }

    /// Print nine items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G, class H, class I>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h, const I& i) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << " " << i << ENDL;
    }

    /// Print ten items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G, class H, class I, class J>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h, const I& i, const J& j) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << " " << i << " " << j << ENDL;
    }

    /// Print eleven items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G, class H, class I, class J, class K>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h, const I& i, const J& j, const K& k) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << " " << i << " " << j << " " << k << ENDL;
    }

    /// Print twelve items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G, class H, class I, class J, class K, class L>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h, const I& i, const J& j, const K& k, const L& l) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << " " << i << " " << j << " " << k << " " << l << ENDL;
    }

    /// Print thirteen items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G, class H, class I, class J, class K, class L, class M>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h, const I& i, const J& j, const K& k, const L& l, const M& m) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << " " << i << " " << j << " " << k << " " << l << " " << m << ENDL;
    }

    /// Print fourteen items separated by spaces to std::cout terminating with new line
    template <class A, class B, class C, class D, class E, class F, class G, class H, class I, class J, class K, class L, class M, class N>
    void print(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h, const I& i, const J& j, const K& k, const L& l, const M& m, const N& n) {
        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << " " << i << " " << j << " " << k << " " << l << " " << m << " " << n << ENDL;
    }

    /// Print a string justified on the left to start at the given column with optional underlining
    void print_justified(const char* s, int column=0, bool underline=true);

    /// Print a string centered at the given column with optional underlining
    void print_centered(const char* s, int column=40, bool underline=true);

}
#endif // MADNESS_WORLD_PRINT_H__INCLUDED
