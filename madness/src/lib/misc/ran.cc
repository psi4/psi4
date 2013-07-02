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
#include <misc/ran.h>
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <limits.h>
using std::fopen;
using std::fwrite;


namespace madness {


    void Random::generate() {
        // Assume we have the lock when we come in here

        double * restrict ur = const_cast<double*>(u);
        double * restrict us = const_cast<double*>(u)+r-s;
        for (int i=0; i<s; ++i) {
            double t = ur[i] + beta*us[i];
            ur[i] = t - int(t);
        }

        ur = const_cast<double*>(u)+s;
        us = const_cast<double*>(u);
        int rs = r-s;
        for (int i=0; i<rs; ++i) {
            double t = ur[i] + beta*us[i];
            ur[i] = t - int(t);
        }

        cur = 0;
    }

    unsigned int Random::simple() {
        return simple_state = 1103515245U*simple_state + 12345U;
    }


    Random::Random(unsigned int seed) : r(1279), s(861), beta(7.0), cur(0), u(new double [r]) {
        // If you switch r don't forget to change size in RandomState

        // a) not clear if beta != 1 is an improvement.
        // b) must ensure s >= r/2.
        // c) r=19937, s=10095 or r=1279, s=861 seem worse than 4423/3004 ?
        // Random(unsigned int seed = 5461) : r(4423), s(3004), beta(7.0), cur(0) {

        setstate(seed);
    }

    void Random::setstate(unsigned int seed) {
        ScopedMutex<Mutex> safe(this);
        // Initialize startup generator
        if ((seed&1) == 0) seed += 1;
        simple_state = seed;
        for (int i=0; i<10*r; ++i) simple();

        // Initialize stream with 48 bit values by first generating
        // roughly 52 bit values and then truncating to exactly 48 bits
        double two21 = 2097152.0;
        double two52 = 4503599627370496.0;
        double two24 = 16777216.0;
        double rtwo24 = 5.9604644775390625e-08;
        for (int i=0; i<r; ++i) u[i] = double(simple());
        for (int i=0; i<r; ++i) u[i] += double(simple())*two21;
        for (int i=0; i<r; ++i) u[i] /= two52;
        // Next line breaks on Cray X1 CC 5.3 ... sigh
        //for (int i=0; i<r; ++i) u[i] -= int(u[i]);
        for (int i=0; i<r; ++i) {int tmp=int(u[i]); u[i] -= double(tmp);}
        for (int i=0; i<r; ++i) {
            int high = int(two24*u[i]);
            int lo = int(two24*(two24*u[i]-high));
            u[i] = rtwo24*(high + rtwo24*lo);
        }

        // Verify that we have only set 48 bits and that at least
        // some of the 48th bits are set.
        double n48 = 0;
        for (int i=0; i<r; ++i) {
            double rem = u[i]*two24;
            rem -= int(rem);
            double rem48 = rem*two24;
            rem48 -= int(rem48);
            if (rem48 != 0) throw "Random: bad bits?";

            double rem47 = rem*two24*0.5;
            rem47 -= int(rem47);
            if (rem47 != 0) ++n48;
        }
        if (n48 == 0) throw "Random: bad 48'th bit?";

        // Warm up
        for (int i=0; i<2000; ++i) generate();
    }

    Random::~Random() {delete [] u;}


    //     // Might be faster than the one-byte-per-word version
    //     void getbytes(int n, unsigned char * restrict v) {
    //         const double two24 = 16777216.0;
    //         int n6 = n%6;
    //         n -= n6;
    //         while (n>0) {
    //             if (cur >= r) generate();
    //             int ndo = std::min(n/6,r-cur);
    //             for (int i=0; i<ndo; ++i) {
    //                 unsigned int high = int(two24*u[i+cur]);
    //                 unsigned int lo = int(two24*(two24*u[i+cur]-high));
    //                 *v++ = (high>>16)&0xff;
    //                 *v++ = (high>>8 )&0xff;
    //                 *v++ = (high    )&0xff;
    //                 *v++ = (lo>>16 )&0xff;
    //                 *v++ = (lo>>8  )&0xff;
    //                 *v++ = (lo     )&0xff;
    //             };
    //             n -= ndo*6;
    //             cur += ndo;
    //         }
    //         for (int i=0; i<n6; ++i) *v++ = (unsigned char) (256*get());
    //     };

    void Random::getbytes(int n, unsigned char * restrict v) {
        ScopedMutex<Mutex> safe(this);
        while (n) {
            if (cur >= r) generate();
            int ndo = std::min(n,r-cur);
            const double* ucur = const_cast<const double*>(u) + cur;
            for (int i=0; i<ndo; ++i) v[i] = (unsigned char) (256.0*ucur[i]);
            n -= ndo;
            v += ndo;
            cur += ndo;
        }
    }

    RandomState Random::getstate() const {
        ScopedMutex<Mutex> safe(this);
        RandomState s;
        s.cur = cur;
        for (int i=0; i<r; ++i) s.u[i] = u[i];
        return s;
    }

    void Random::setstate(const RandomState &s) {
        ScopedMutex<Mutex> safe(this);
        cur = s.cur;
        for (int i=0; i<r; ++i) u[i] = s.u[i];
    }


    void Random::test() {
        // Crude tests just to find gross errors ... the file written
        // at the end is used as input to diehard.
        Random r;
        double sum1 = 0;
        double sum2 = 0;
        double x12 = 0;
        double two24 = 16777216.0;

        double n = 100000000.0;
        double nn = n;
        const int nbuf = 64;
        nn /= nbuf;
        n = nn*nbuf;

        double buf[nbuf];

        while (nn--) {
            r.getv(nbuf,buf);
            for (int i=0; i<nbuf; ++i) {
                double d1 = buf[i];
                sum1 += d1;
                double d2 = d1 * two24;
                d2 -= int(d2);
                sum2 += d2;
                x12 += (d1-0.5)*(d2-0.5);
            }
        }

        std::cout << "high   " << sum1/n << std::endl;
        std::cout << "lo     " << sum2/n << std::endl;
        std::cout << "hi-lo  " << x12/n << std::endl;

        const int nb = 100000000;
        unsigned char *b = new unsigned char[nb+1];

        b[nb-1] = 0;
        b[nb] = 99;
        r.getbytes(nb,b);

        std::cout << int(b[nb-1]) << std::endl;
        std::cout << int(b[nb]) << std::endl;

        FILE *f = fopen("stream","w+");
        if (!f) {std::cout << "fopen?\n"; std::exit(1);}
        fwrite(b, 1, nb, f);
        fclose(f);
    }

    Random default_random_generator;


    template <> double RandomValue<double> () {
        return default_random_generator.get();
    }

    template <> float RandomValue<float> () {
        return float(default_random_generator.get());
    }

    template <> double_complex RandomValue<double_complex> () {
        return double_complex(RandomValue<double>(),RandomValue<double>());
    }

    template <> float_complex RandomValue<float_complex> () {
        return float_complex(RandomValue<float>(),RandomValue<float>());
    }

    template <> int RandomValue<int> () {
        return int(2147483648.0 * RandomValue<double>());
    }

    template <> long RandomValue<long> () {
        return long(2147483648.0 * RandomValue<double>());
    }

    template <> void RandomVector<double>(int n, double* t) {
        default_random_generator.getv(n, t);
    }

    template <> void RandomVector<float>(int n, float* t) {
        default_random_generator.getv(n, t);
    }

    template <> void RandomVector<double_complex>(int n, double_complex* t) {
        default_random_generator.getv(2*n, (double*)(t));
    }

    template <> void RandomVector<float_complex>(int n, float_complex* t) {
        default_random_generator.getv(2*n, (float*)(t));
    }
}

// using namespace madness;

// int main() {
//     Random::test();
//     return 0;
// }
