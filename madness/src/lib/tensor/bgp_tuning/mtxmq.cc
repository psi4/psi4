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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>

using namespace std;

typedef complex<double> double_complex;

void mTxmq(long dimi, long dimj, long dimk, double* c, const double* a, const double* b) {
    for (long i=0; i<dimi; ++i,c+=dimj,++a) {
        for (long j=0; j<dimj; ++j) c[j] = 0.0;
        const double *ai = a;
        for (long k=0; k<dimk; ++k,ai+=dimi) {
            const double aki = *ai;
            for (long j=0; j<dimj; ++j) {
                c[j] += aki*b[k*dimj+j];
            }
        }
    }
}

void mTxmq(long dimi, long dimj, long dimk, double_complex* c, const double_complex* a, const double* b) {
    for (long i=0; i<dimi; ++i,c+=dimj,++a) {
        for (long j=0; j<dimj; ++j) c[j] = 0.0;
        const double_complex *ai = a;
        for (long k=0; k<dimk; ++k,ai+=dimi) {
            const double_complex aki = *ai;
            for (long j=0; j<dimj; ++j) {
                c[j] += aki*b[k*dimj+j];
            }
        }
    }
}

void mTxmq(long dimi, long dimj, long dimk, double_complex* c, const double_complex* a, const double_complex* b) {
    for (long i=0; i<dimi; ++i,c+=dimj,++a) {
        for (long j=0; j<dimj; ++j) c[j] = 0.0;
        const double_complex *ai = a;
        for (long k=0; k<dimk; ++k,ai+=dimi) {
            const double_complex aki = *ai;
            for (long j=0; j<dimj; ++j) {
                c[j] += aki*b[k*dimj+j];
            }
        }
    }
}

