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
#include <tensor/tensor.h>
#include <tensor/tensor.h>
#include <tensor/mtxmq.h>
#include <world/print.h>
#include <world/worldtime.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace madness;
using namespace std;

// void aligned_add(long n, double* restrict a, const double* restrict b) {
//     long n4 = (n>>2)<<2;
//     long rem = n-n4;
//     for (long i=0; i<n4; i+=4,a+=4,b+=4) {
//         a[0] += b[0];
//         a[1] += b[1];
//         a[2] += b[2];
//         a[3] += b[3];
//     }
//     for (long i=0; i<rem; ++i) *a++ += *b++;
// }

// void aligned_sub(long n, double* restrict a, const double* restrict b) {
//     long n4 = (n>>2)<<2;
//     long rem = n-n4;
//     for (long i=0; i<n4; i+=4,a+=4,b+=4) {
//         a[0] -= b[0];
//         a[1] -= b[1];
//         a[2] -= b[2];
//         a[3] -= b[3];
//     }
//     for (long i=0; i<rem; ++i) *a++ -= *b++;
// }

int main() {

    const long k=10; // The polynomial rank
    const long twok = 2*k;
    const long NBOX=100; // No. of boxes
    const long RANK=30;   // No. of terms in separated expansion

    Tensor<double> t1(twok,twok,twok); // Workspace
    Tensor<double> t2(twok,twok,twok);
    double* tmp1=t1.ptr();
    double* tmp2=t2.ptr();

    vector< Tensor<double> > Xlist, Ylist, Zlist;
    vector< Tensor<double> > XlistS, YlistS, ZlistS;
    for (long m=0; m<RANK*27; ++m) {
        Xlist.push_back(Tensor<double>(twok,twok));
        Ylist.push_back(Tensor<double>(twok,twok));
        Zlist.push_back(Tensor<double>(twok,twok));

        XlistS.push_back(Tensor<double>(k,k));
        YlistS.push_back(Tensor<double>(k,k));
        ZlistS.push_back(Tensor<double>(k,k));
    }

    double start = cpu_time();

    // Significant boxes and their coefficients are normally obtained
    // by iterating thru entries in a sparse array (hash table)
    for (long box=0; box<NBOX; ++box) {
        const Tensor<double> c(twok,twok,twok); // The coefficients for box

        if ((box%10) == 0) print("doing box", box);

        // Each of the neighbors of a cube in 3D including self
        for (long neighbor=0; neighbor<27; ++neighbor) {
            Tensor<double> r(twok,twok,twok); // This will hold the result
            Tensor<double> s(k,k,k);          // This will hold the result

            // Loop thru terms in the expansion of this operator
            for (long m=0; m<RANK; ++m) {
                // In practice low-rank approximations are used

                const Tensor<double>& X = Xlist[neighbor*27+m];
                const Tensor<double>& Y = Ylist[neighbor*27+m];
                const Tensor<double>& Z = Zlist[neighbor*27+m];

                mTxmq(twok*twok, twok, twok, tmp1, c.ptr(), X.ptr());
                mTxmq(twok*twok, twok, twok, tmp2, tmp1, Y.ptr());
                mTxmq(twok*twok, twok, twok, tmp1, tmp2, Z.ptr());

                aligned_add(twok*twok*twok, r.ptr(), tmp1);

                const Tensor<double>& XS = XlistS[neighbor*27+m];
                const Tensor<double>& YS = YlistS[neighbor*27+m];
                const Tensor<double>& ZS = ZlistS[neighbor*27+m];

                mTxmq(k*k, k, k, tmp1, c.ptr(), XS.ptr());
                mTxmq(k*k, k, k, tmp2, tmp1, YS.ptr());
                mTxmq(k*k, k, k, tmp1, tmp2, ZS.ptr());

                aligned_add(k*k*k, s.ptr(), tmp1);
            }
            aligned_sub(k*k*k, r.ptr(), s.ptr());

            // store r into result sparse hash table
        }
    }
    double used = cpu_time() - start;
    double nflops = NBOX*RANK*3.0*2.0*27.0*(twok*twok*twok*twok + k*k*k*k);
    double rate = nflops/used;
    print("used ", used,"s nflop", nflops, "rate", rate);

    return 0;
}
