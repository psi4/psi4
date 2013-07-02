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

  $Id: ran.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/
#ifndef MADNESS_MISC_RAN_H__INCLUDED
#define MADNESS_MISC_RAN_H__INCLUDED

#include <madness_config.h>
#include <world/worldthread.h>

#include <complex>
typedef std::complex<float> float_complex;
typedef std::complex<double> double_complex;


namespace madness {

    struct RandomState {
        int cur;
        double u[1279];
    };

    /// A random number generator (portable, vectorized, and thread-safe)

    /// Following Brent 1992, we use a 48-bit generalized Fibonacci generator
    /// \code
    ///     u[n] = alpha*u[n-r] + beta*u[n-s] mod m
    /// \endcode
    /// with alpha=1, beta=7, r=1279, s=861, m=2^48.  Double precision
    /// numbers are used to perform exact integer arithmetic.  48-bit
    /// because we have 52 bits of mantissa, alpha+1 is 3 bits and 1 bit spare.
    ///
    /// The period is nominally 2^m (2^r - 1) / 2 but if p is the period,
    /// X[n] and X[n+p/2k] differ by at most k bits (0 < k < 48) so usage
    /// should be limited to the first 2^r-1 entries (about 10^385 values).
    ///
    /// Each instance provides a separate stream, but it is up to the
    /// user to partition the sequence by selecting distinct seeds or
    /// other means.
    ///
    /// The streams are thread safe.
    ///
    /// A default stream is provided as madness::default_random_generator.
    class Random : private Mutex {
    private:
        const int r;
        const int s;
        const double beta;
        volatile int cur;
        volatile double* const u;
        unsigned int simple_state;

        void generate();

        unsigned int simple();

    public:
        Random(unsigned int seed = 5461);

        virtual ~Random();

        double get() {
            ScopedMutex<Mutex> safe(this);
            if (cur >= r) generate();
            return u[cur++];
        }

        /// Returns a vector of uniform doubles in [0,1)
        template <typename T>
        void getv(int n, T * restrict v) {
            ScopedMutex<Mutex> safe(this);
            while (n) {
                if (cur >= r) generate();
                int ndo = std::min(n,r-cur);
                const double* ucur = const_cast<const double*>(u) + cur;
                for (int i=0; i<ndo; ++i) v[i] = (T)(ucur[i]);
                n -= ndo;
                v += ndo;
                cur += ndo;
            }
        }

        /// Returns vector of random bytes in [0,256)
        void getbytes(int n, unsigned char * restrict v);

        /// Returns full state of the generator
        RandomState getstate() const;

        /// Restores state of the generator
        void setstate(const RandomState &s);

        /// Sets state of the generator from integer
        void setstate(unsigned int seed);

        /// Test the generator
        static void test();
    };


    /// The default random number stream
    extern Random default_random_generator;

    /// Random value that wraps the default Fibonacci generator
    template <class T> T RandomValue();

    /// Random double
    template <> double RandomValue<double> ();

    /// Random float
    template <> float RandomValue<float> ();

    /// Random int
    template <> int RandomValue<int> ();

    /// Random long
    template <> long RandomValue<long> ();

    /// Random double_complex
    template <> double_complex RandomValue<double_complex> ();

    /// Random float_complex
    template <> float_complex RandomValue<float_complex> ();

    template <class T> void RandomVector(int n, T* t) {
        for (int i=0; i<n; ++i) t[i] = RandomValue<T>();
    }

    template <> void RandomVector<double>(int n, double* t);

    template <> void RandomVector<float>(int n, float* t);

    template <> void RandomVector<double_complex>(int n, double_complex* t);

    template <> void RandomVector<float_complex>(int n, float_complex* t);
}

#endif // MADNESS_MISC_RAN_H__INCLUDED
