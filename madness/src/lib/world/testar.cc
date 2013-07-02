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


  $Id: testar.cc 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#include <iostream>
using std::cout;
using std::endl;

#include <cmath>

/// \file test.cc
/// \brief Tests serialization by some of the archives

#define ARCHIVE_REGISTER_TYPE_INSTANTIATE_HERE

#include <world/textfsar.h>
using madness::archive::TextFstreamInputArchive;
using madness::archive::TextFstreamOutputArchive;

#include <world/binfsar.h>
using madness::archive::BinaryFstreamInputArchive;
using madness::archive::BinaryFstreamOutputArchive;

#include <world/vecar.h>
using madness::archive::VectorInputArchive;
using madness::archive::VectorOutputArchive;

#include <world/bufar.h>
using madness::archive::BufferInputArchive;
using madness::archive::BufferOutputArchive;


// A is a class that provides a symmetric serialize method
class A {
public:
    float a;
    A(float a = 0.0) : a(a) {};

    template <class Archive>
    inline void serialize(const Archive& ar) {
        ar & a;
    }
};


// B is a class without a serialize method but with symmetric serialization.
class B {
public:
    bool b;
    B(bool b = false) : b(b) {};
};

namespace madness {
    namespace archive {

        template <class Archive>
        struct ArchiveSerializeImpl<Archive,B> {
            static inline void serialize(const Archive& ar, B& b) {
                ar & b.b;
            };
        };
    }
}

// C is a class with asymmetric load/store.
class C {
public:
    long c;
    C(long c = 0) : c(c) {};
};

namespace madness {
    namespace archive {
        template <class Archive>
        struct ArchiveLoadImpl<Archive,C> {
            static inline void load(const Archive& ar, C& c) {
                ar >> c.c;
            };
        };

        template <class Archive>
        struct ArchiveStoreImpl<Archive,C> {
            static inline void store(const Archive& ar, const C& c) {
                ar << c.c;
            };
        };
    }
}


// A better example of a class with asym load/store
class linked_list {
    int value;
    linked_list *next;
public:
    linked_list(int value = 0) : value(value), next(0) {};

    void append(int value) {
        if (next) next->append(value);
        else next = new linked_list(value);
    };

    void set_value(int val) {
        value = val;
    };

    int get_value() const {
        return value;
    };

    linked_list* get_next() const {
        return next;
    };

    void print() {
        if (next) {
            cout << value << " --> ";
            next->print();
        }
        else cout << value << endl;
    };
};

namespace madness {
    namespace archive {
        template <class Archive>
        struct ArchiveStoreImpl<Archive,linked_list> {
            static void store(const Archive& ar, const linked_list& c) {
                ar & c.get_value() & bool(c.get_next());
                if (c.get_next()) ar & *c.get_next();
            };
        };

        template <class Archive>
        struct ArchiveLoadImpl<Archive,linked_list> {
            static void load(const Archive& ar, linked_list& c) {
                int value;
                bool flag;
                ar & value & flag;
                c.set_value(value);
                if (flag) {
                    c.append(0);
                    ar & *c.get_next();
                }
            };
        };
    }
}

namespace madness {
    namespace archive {
        typedef std::map< short,std::complex<double> > map_short_complex_double;
        typedef std::pair< short,std::complex<double> > pair_short_complex_double;
        typedef std::pair<int,double> pair_int_double;
        ARCHIVE_REGISTER_TYPE_AND_PTR(A,128);
        ARCHIVE_REGISTER_TYPE_AND_PTR(B,129);
        ARCHIVE_REGISTER_TYPE_AND_PTR(C,130);
        ARCHIVE_REGISTER_TYPE_AND_PTR(linked_list,131);
        ARCHIVE_REGISTER_TYPE_AND_PTR(pair_int_double,132);
        ARCHIVE_REGISTER_TYPE_AND_PTR(map_short_complex_double,133);
        ARCHIVE_REGISTER_TYPE_AND_PTR(pair_short_complex_double, 134);
    }
}

using namespace std;

using madness::archive::wrap;

typedef std::complex<double> double_complex;

template <class OutputArchive>
void test_out(const OutputArchive& oar) {
    const int n = 3;
    A a, an[n];
    B b, bn[n];
    C c, cn[n];
    int i, in[n];
    double *p = new double[n];
    A *q = new A[n];
    vector<int> v(3);
    pair<int,double> pp(33,99.0);
    map<short,double_complex> m;
    const char* teststr = "hello \n dude !";
    string str(teststr);
    linked_list list(0);
    double pi = atan(1.0)*4.0;
    double e = exp(1.0);

    // Initialize data
    a.a = b.b = c.c = i = 1;
    for (int k=0; k<n; ++k) {
        p[k] = q[k].a = an[k].a = v[k] = cn[k].c = in[k] = k;
        bn[k].b = k&1;
        m[k] = double_complex(k,k);
        list.append(k+1);
    }

    // test example list code
    list.print();
    oar & list;

    oar & pi & e;

    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " constant double value" << std::endl);
    oar & 1.0;
    oar << 1.0;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " int" << std::endl);
    oar & i;
    oar << i;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " A" << std::endl);
    oar & a;
    oar << a;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " B" << std::endl);
    oar & b;
    oar << b;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " C" << std::endl);
    oar & c;
    oar << c;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " int[]" << std::endl);
    oar & in;
    oar << in;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " A[]" << std::endl);
    oar & an;
    oar << an;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " B[]" << std::endl);
    oar << bn;
    oar & bn;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " C[]" << std::endl);
    oar << cn;
    oar & cn;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " double *p wrapped" << std::endl);
    oar << wrap(p,n);
    oar & wrap(p,n);
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " A *q wrapped" << std::endl);
    oar << wrap(q,n);
    oar & wrap(q,n);
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " vector<int>" << std::endl);
    oar << v;
    oar & v;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " pair<int,double>" << std::endl);
    oar << pp;
    oar & pp;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " map<short,double<complex>>" << std::endl);
    oar << m;
    oar & m;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " string" << std::endl);
    oar << str;
    oar & str;

    oar & 1.0 & i & a & b & c & in & an & bn & cn & wrap(p,n) & wrap(q,n) & pp & m & str;
}

template <class InputArchive>
void test_in(const InputArchive& iar) {
    const int n = 3;
    A a, an[n];
    B b, bn[n];
    C c, cn[n];
    int i, in[n];
    double *p = new double[n];
    A *q = new A[n];
    vector<int> v(3);
    pair<int,double> pp(33,99.0);
    map<short,double_complex> m;
    const char* teststr = "hello \n dude !";
    string str(teststr);
    linked_list list;
    double pi, e;

    // Destroy in-core data
    a.a = b.b = c.c = i = 0;
    for (int k=0; k<n; ++k) {
        p[k] = q[k].a = an[k].a = v[k] = cn[k].c = in[k] = -1;
        bn[k].b = (k+1)&1;
        m[k] = double_complex(0,0);
    }
    pp = pair<int,double>(0,0.0);
    str = string("");

    iar & list;
    list.print();

    iar & pi & e;
    cout.setf(std::ios::scientific);
    cout << "error in pi " << (pi - 4.0*atan(1.0)) << endl;
    cout << "error in e " << (e - exp(1.0)) << endl;

    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " constant double value" << std::endl);
    double val;
    iar & val;
    iar >> val;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " int" << std::endl);
    iar & i;
    iar >> i;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " A" << std::endl);
    iar & a;
    iar >> a;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " B" << std::endl);
    iar & b;
    iar >> b;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " C" << std::endl);
    iar & c;
    iar >> c;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " int[]" << std::endl);
    iar & in;
    iar >> in;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " A[]" << std::endl);
    iar & an;
    iar >> an;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " B[]" << std::endl);
    iar & bn;
    iar >> bn;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " C[]" << std::endl);
    iar & cn;
    iar >> cn;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " double *p wrapped" << std::endl);
    iar & wrap(p,n);
    iar >> wrap(p,n);
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " A *q wrapped" << std::endl);
    iar & wrap(q,n);
    iar >> wrap(q,n);
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " vector<int>" << std::endl);
    iar & v;
    iar >> v;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " pair<int,double>" << std::endl);
    iar & pp;
    iar >> pp;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " map<short,double<complex>>" << std::endl);
    iar & m;
    iar >> m;
    MAD_ARCHIVE_DEBUG(std::cout << std::endl << " string" << std::endl);
    iar & str;
    iar >> str;

    iar & 1.0 & i & a & b & c & in & an & bn & cn & wrap(p,n) & wrap(q,n) & pp & m & str;
    // Test data
    bool status = true;

#define TEST(cond) status &= cond;  \
                   if (!(cond)) std::cout << #cond << " failed" << std::endl
    TEST(a.a == 1);
    TEST(b.b == 1);
    TEST(c.c == 1);
    TEST(i == 1);
    for (int k=0; k<n; ++k) {
        TEST(an[k].a == k);
        TEST(bn[k].b == (k&1));
        TEST(cn[k].c == k);
        TEST(in[k] == k);
        TEST(p[k] == k);
        TEST(q[k].a == k);
        TEST(v[k] == k);
        TEST(m[k] == double_complex(k,k));
    }
    TEST(pp.first==33 && pp.second==99.0);
    TEST(str == string(teststr));

#undef TEST

    if (status)
        std::cout << "Serialization appears to work." << std::endl;
    else
        std::cout << "Sorry, back to the drawing board.\n";
}

int main() {
    madness::archive::archive_initialize_type_names();
    ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(A);
    ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(B);
    ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(C);
    ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(linked_list);
    ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(madness::archive::pair_int_double);
    ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(madness::archive::pair_short_complex_double);
    ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(madness::archive::map_short_complex_double);

    {
        const char* f = "test.dat";
        cout << endl << "testing binary fstream archive" << endl;
        BinaryFstreamOutputArchive oar(f);
        test_out(oar);
        oar.close();

        BinaryFstreamInputArchive iar(f);
        test_in(iar);
        iar.close();
    }

    {
        cout << endl << "testing vector archive" << endl;
        std::vector<unsigned char> f;
        VectorOutputArchive oar(f);
        test_out(oar);
        oar.close();

        VectorInputArchive iar(f);
        test_in(iar);
        iar.close();
    }

    {
        const char* f = "test.dat";
        cout << endl << "testing text fstream archive" << endl;
        TextFstreamOutputArchive oar(f);
        test_out(oar);
        oar.close();

        TextFstreamInputArchive iar(f);
        test_in(iar);
        iar.close();
    }

    {
        cout << endl << "testing buffer archive" << endl;
        unsigned char buf[32768];
        BufferOutputArchive oar(buf,sizeof(buf));
        test_out(oar);
        std::size_t nbyte = oar.size();
        oar.close();

        BufferInputArchive iar(buf,nbyte);
        test_in(iar);
        iar.close();
    }

    return 0;
}
