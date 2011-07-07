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


  $Id: twoscale.cc 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#include <world/world.h>
#include <iostream>
using std::cout;
using std::endl;
using std::ios;

#include <cstdio>
#include <cstdlib>

#include <cmath>
using std::abs;

#include <mra/twoscale.h>
#include <tensor/tensor.h>
#include <misc/misc.h>

/// \file twoscale.cc
/// \brief Routines to provide twoscale & correlation coeffs for Legendre basis

namespace madness {

    static const int kmax = 60;
    static const char *twoscale_filename = "coeffs";  // Will be overridden by load_coeffs
    static const char *autocorr_filename = "autocorr";  // Will be overriden by load_coeff


    static class twoscale_cache_class {
        /// This caches the two-scale coefficients
    public:
        Tensor<double> h0;
        Tensor<double> h1;
        Tensor<double> g0;
        Tensor<double> g1;
    }
    cache[kmax+1];

    static bool loaded = 0;


    static Tensor<double> readmat(int k, FILE* file) {
        Tensor<double> a(k,k);
        for (int i=0; i<k; ++i) {
            for (int j=0; j<k; ++j) {
                double c;
                if (fscanf(file,"%lf",&c) != 1) {
                    cout << "readmat: twoscale missing coeff?\n";
                    throw "readmat";
                }
                a(i,j) = c;
            }
        }
        return a;
    }

    static inline double phase(long i) {
        return (i&1) ? -1.0 : 1.0;
    }

    static bool read_twoscale(int kmax) {
        unsigned long correct = 6931979l;
        unsigned long computed = checksum_file(twoscale_filename);
        MADNESS_ASSERT(correct == computed);
        FILE* file = fopen(twoscale_filename,"r");
        if (!file) {
            cout << "twoscale: failed opening file with twoscale coefficients\n";
            return false;
        }
        for (int k=1; k<kmax+1; ++k) {
            Tensor <double> h0, g0;
            try {
                h0 = readmat(k,file);
                g0 = readmat(k,file);
            }
            catch (char *e) {
                fclose(file);
                return false;
            }

            Tensor<double> h1(k,k);
            Tensor<double> g1(k,k);

            for (int i=0; i<k; ++i) {
                for (int j=0; j<k; ++j) {
                    h1(i,j) = h0(i,j)*phase(i+j);
                    g1(i,j) = g0(i,j)*phase(i+j+k);
                }
            }

            cache[k].h0 = h0;
            cache[k].h1 = h1;
            cache[k].g0 = g0;
            cache[k].g1 = g1;
        }
        fclose(file);

        loaded = true;
        return true;
    }


    /// Return the two scale coefficients in the Legendre basis

    /// Returns true on success, false on failure (and prints error message)
    bool two_scale_coefficients(int k,
                                Tensor<double>* h0, Tensor<double>* h1,
                                Tensor<double>* g0, Tensor<double>* g1) {
        if (!loaded) {
            if (!read_twoscale(kmax)) return false;
        }

        if (k < 1 || k > kmax) return false;

        *h0 = copy(cache[k].h0);
        *h1 = copy(cache[k].h1);
        *g0 = copy(cache[k].g0);
        *g1 = copy(cache[k].g1);

        return true;
    }

    bool two_scale_hg(int k, Tensor<double>* hg) {
        Tensor<double> h0(k,k), h1(k,k), g0(k,k), g1(k,k);

        if (!two_scale_coefficients(k, &h0, &h1, &g0, &g1)) return false;

        *hg = Tensor<double>(2*k,2*k);

        Slice sk(0,k-1), sk2(k,-1);
        (*hg) = Tensor<double>(2*k,2*k);
        (*hg)(sk,sk)   = h0;
        (*hg)(sk,sk2)  = h1;
        (*hg)(sk2,sk)  = g0;
        (*hg)(sk2,sk2) = g1;

        return true;
    }

    bool test_two_scale_coefficients() {
        /// Test the two scale coefficients for orthogonality
        for (int k=1; k<kmax; ++k) {
            Tensor<double> hg;
            if (!two_scale_hg(k,&hg)) return false;
            Tensor<double> ident(2*k,2*k);
            for (int i=0; i<2*k; ++i) ident(i,i) = 1.0;

            double err0 = (inner(hg,hg,0,0)-ident).absmax();
            if (err0 > 9e-16) {
                std::cout << "twoscale failed 0: " << k << " " << err0 << std::endl;
                std::cout << (inner(hg,hg,0,0)-ident);
                return false;
            }

            double err1 = (inner(hg,hg,1,1)-ident).absmax();
            if (err1 > 9e-16) {
                std::cout << "twoscale failed 1: " << k << " " << err1 << std::endl;
                std::cout << (inner(hg,hg,1,1)-ident);
                return false;
            }
        }
        return true;
    }


    // BELOW HERE THE AUTOCORRELATION ROUTINES

    static const int kmax_autoc = 30;
    static int kread = -1;  // value of k for data read from disk into _cread
    static Tensor<double> _cread;

    static int kcur = -1; // current value of k in for data in _c
    static Tensor<double> _c;

    static bool read_data(int k);

    /// Return the autocorrelation coefficients for scaling functions of given order

    /// Returned is a view of the cached data ... do not modify.
    ///
    /// The autocorrelation functions are defined as
    /// \code
    ///    Phi_ij(z) = int(0,z+1) phi_i(x) phi_j(x-z)  z<=0
    ///    Phi_ij(z) = int(z,1) phi_i(x) phi_j(x-z)    z>=0
    /// \endcode
    /// and are expanded in the double order Legendre scaling functions on either
    /// side of the origin
    /// \code
    ///    Phi_ij(z) = sum(p) [phi_p(z+1)*cminus_ijp + phi_p(z)*cplus_ijp]
    ///
    ///    cplus_ijp  = int(-1,0) Phi_ij(z) phi_p(z+1)
    ///    cminus_ijp = int(0,1)  Phi_ij(z) phi_p(z)
    /// \endcode
    ///
    /// The expansion coefficients \c cminus and \c cplus have been precomputed
    /// with Maple to a precision of about 1e-30.  Only \c cplus is stored and we use
    /// \code
    ///    cplus(i,j,p) = (-1)**(i+j+p) * cminus(i,j,p)
    ///    cplus(i,j,p) = (-1)**(i+j)   * cplus(j,i,p)
    /// \endcode
    ///
    /// The returned tensor concatenates cminus and cplus.
    ///
    /// Return true on success, false on failure (which also prints error
    /// message to stdout).
    bool autoc(int k, Tensor<double>* c) {
        if (k < 1 || k > kread) { // was kmax_autoc
            cout << "autoc: invalid k " << k << endl;
            return false;
        }

        // The data is cached in _c ... reread (now from _cread) for each new value of k
        if (kcur != k) {
            _c = Tensor<double>(k,k,4*k);
            _c(Slice(0,k-1),Slice(0,k-1),Slice(0,2*k-1)) = _cread(Slice(0,k-1),Slice(0,k-1),Slice(0,2*k-1));
            _c(Slice(0,k-1),Slice(0,k-1),Slice(2*k,4*k-1)) = _cread(Slice(0,k-1),Slice(0,k-1),Slice(2*kread,2*kread+2*k-1));
            kcur = k;
        }

        *c = _c;
        return true;
    }

    /// Return +1 if i is even, -1 if i is odd ... perverse no if-test form
    static inline int parity(int i) {
        return 1 - ((i&1)<<1);
    }

    bool test_autoc() {
        unsigned long correct = 9056188; // 0x638a9b;
        unsigned long computed = checksum_file(autocorr_filename);
        if (correct != computed)
            cout << "test_autoc: file checksum invalid: correct="
                 << correct << " computed=" << computed << endl;

        return (correct == computed);
    }

    static bool read_data(int k) {
        if (!test_autoc()) return false;
        kread = -1;
        FILE *file = fopen(autocorr_filename,"r");
        if (!file) {
            cout << "autoc: failed opening file with autocorrelation coefficients" << endl;
            return false;
        }

        _cread = Tensor<double>(k,k,4*k);

        long twok = 2*k;
        while (1) {
            long i, j, p;
            double val;
            if (fscanf(file,"%ld %ld %ld %lf",&i,&j,&p,&val) != 4) {
                cout<<"autoc: failed reading file " << endl;
                fclose(file);
                return false;
            }
            if (i >= k) break;
            double ij = parity(i+j);
            double ijp= parity(i+j+p);

            _cread(i,j,p)      = val*ijp;	// c-
            _cread(j,i,p)      = val*ij*ijp;
            _cread(i,j,p+twok) = val;	// c+
            _cread(j,i,p+twok) = val*ij;
        }

        fclose(file);
        kread = k;
        return true;
    }

    /// Collective routine to load and cache twoscale & autorrelation coefficients

    /// Only process rank 0 will access the files.
    void load_coeffs(World& world, const char* dir) {
        if (!loaded) {
            int ktop = kmax_autoc;   // Plausible maximum value
            if (world.rank() == 0) {
                char buf[32768];
                buf[0] = 0;
                strcat(buf,dir);
                strcat(buf,"/");
                strcat(buf,twoscale_filename);
                twoscale_filename = strdup(buf);

                buf[0] = 0;
                strcat(buf,dir);
                strcat(buf,"/");
                strcat(buf,autocorr_filename);
                autocorr_filename = strdup(buf);

                if (!read_twoscale(kmax))
                    throw "load_coeffs: failed reading twoscale coeffs";

                if (!read_data(ktop))
                    throw "load_coeffs: failed reading coeffs";
            }
            else {
                for (int k=1; k<=kmax; ++k) {
                    cache[k].h0 = Tensor<double>(k,k);
                    cache[k].h1 = Tensor<double>(k,k);
                    cache[k].g0 = Tensor<double>(k,k);
                    cache[k].g1 = Tensor<double>(k,k);
                }
                _cread = Tensor<double>(ktop,ktop,4*ktop);
                kread = ktop;
            }

            for (int k=1; k<=kmax; ++k) {
                world.gop.broadcast(cache[k].h0.ptr(), k*k, 0);
                world.gop.broadcast(cache[k].h1.ptr(), k*k, 0);
                world.gop.broadcast(cache[k].g0.ptr(), k*k, 0);
                world.gop.broadcast(cache[k].g1.ptr(), k*k, 0);
            }

            world.gop.broadcast(_cread.ptr(), ktop*ktop*4*ktop, 0);

            loaded = true;
        }
    }

}



