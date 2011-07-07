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
#ifdef X86_64

#include <tensor/tensor.h>
#include <locale>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <xmmintrin.h>
//#include <pmmintrin.h>


#include <tensor/mtxmq.h>

using namespace madness;

void mtxmGooberSaysHey(long dimi, long dimj, long dimk,
                       std::complex<double>* restrict c,
                       const std::complex<double>* a, const double* b)
{
  const long itile = 14;
  for (long j = 0; j < dimj; ++j)
  {
    for (long ilo = 0; ilo < dimi; ilo += itile)
    {
      long ni = dimi - ilo;
      ni = (ni >= itile) ? itile : ni;
      if (ni == 14)
      {
        std::complex<double>* cij0  = c + (ilo     )*dimj + j;
        std::complex<double>* cij1  = c + (ilo +  1)*dimj + j;
        std::complex<double>* cij2  = c + (ilo +  2)*dimj + j;
        std::complex<double>* cij3  = c + (ilo +  3)*dimj + j;
        std::complex<double>* cij4  = c + (ilo +  4)*dimj + j;
        std::complex<double>* cij5  = c + (ilo +  5)*dimj + j;
        std::complex<double>* cij6  = c + (ilo +  6)*dimj + j;
        std::complex<double>* cij7  = c + (ilo +  7)*dimj + j;
        std::complex<double>* cij8  = c + (ilo +  8)*dimj + j;
        std::complex<double>* cij9  = c + (ilo +  9)*dimj + j;
        std::complex<double>* cij10 = c + (ilo + 10)*dimj + j;
        std::complex<double>* cij11 = c + (ilo + 11)*dimj + j;
        std::complex<double>* cij12 = c + (ilo + 12)*dimj + j;
        std::complex<double>* cij13 = c + (ilo + 13)*dimj + j;
//        *(cij0) = std::complex<double>(0.0,0.0);
//        *(cij1) = std::complex<double>(0.0,0.0);
//        *(cij2) = std::complex<double>(0.0,0.0);
//        *(cij3) = std::complex<double>(0.0,0.0);
//        *(cij4) = std::complex<double>(0.0,0.0);
//        *(cij5) = std::complex<double>(0.0,0.0);
//        *(cij6) = std::complex<double>(0.0,0.0);
//        *(cij7) = std::complex<double>(0.0,0.0);
//        *(cij8) = std::complex<double>(0.0,0.0);
//        *(cij9) = std::complex<double>(0.0,0.0);
//        *(cij10) = std::complex<double>(0.0,0.0);
//        *(cij11) = std::complex<double>(0.0,0.0);
//        *(cij12) = std::complex<double>(0.0,0.0);
//        *(cij13) = std::complex<double>(0.0,0.0);
        __asm__ __volatile__ (
            "pxor %xmm2,%xmm2;"
            "pxor %xmm3,%xmm3;"
            "pxor %xmm4,%xmm4;"
            "pxor %xmm5,%xmm5;"
            "pxor %xmm6,%xmm6;"
            "pxor %xmm7,%xmm7;"
            "pxor %xmm8,%xmm8;"
            "pxor %xmm9,%xmm9;"
            "pxor %xmm10,%xmm10;"
            "pxor %xmm11,%xmm11;"
            "pxor %xmm12,%xmm12;"
            "pxor %xmm13,%xmm13;"
            "pxor %xmm14,%xmm14;"
            "pxor %xmm15,%xmm15;"
        );
        for (long k = 0; k < dimk; ++k)
        {
//          const std::complex<double>* aki0  = a + k*dimi + ilo;
//          const std::complex<double>* aki1  = a + k*dimi + ilo + 1;
//          const std::complex<double>* aki2  = a + k*dimi + ilo + 2;
//          const std::complex<double>* aki3  = a + k*dimi + ilo + 3;
//          const std::complex<double>* aki4  = a + k*dimi + ilo + 4;
//          const std::complex<double>* aki5  = a + k*dimi + ilo + 5;
//          const std::complex<double>* aki6  = a + k*dimi + ilo + 6;
//          const std::complex<double>* aki7  = a + k*dimi + ilo + 7;
//          const std::complex<double>* aki8  = a + k*dimi + ilo + 8;
//          const std::complex<double>* aki9  = a + k*dimi + ilo + 9;
//          const std::complex<double>* aki10 = a + k*dimi + ilo + 10;
//          const std::complex<double>* aki11 = a + k*dimi + ilo + 11;
//          const std::complex<double>* aki12 = a + k*dimi + ilo + 12;
//          const std::complex<double>* aki13 = a + k*dimi + ilo + 13;
          const std::complex<double>* aki  = a + k*dimi + ilo;
          const double* bkj = b + k*dimj + j;
//          *(cij0) += *(aki0) * *(bkj);
//          *(cij1) += *(aki1) * *(bkj);
//          *(cij2) += *(aki2) * *(bkj);
//          *(cij3) += *(aki3) * *(bkj);
//          *(cij4) += *(aki4) * *(bkj);
//          *(cij5) += *(aki5) * *(bkj);
//          *(cij6) += *(aki6) * *(bkj);
//          *(cij7) += *(aki7) * *(bkj);
//          *(cij8) += *(aki8) * *(bkj);
//          *(cij9) += *(aki9) * *(bkj);
//          *(cij10) += *(aki10) * *(bkj);
//          *(cij11) += *(aki11) * *(bkj);
//          *(cij12) += *(aki12) * *(bkj);
//          *(cij13) += *(aki13) * *(bkj);
          __asm__ volatile(
            "movddup (%1), %%xmm0;"

            "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;"
            "movapd  16(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;"
            "movapd  32(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4;"
            "movapd  48(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm5;"
            "movapd  64(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm6;"
            "movapd  80(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm7;"
            "movapd  96(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm8;"
            "movapd 112(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm9;"
            "movapd 128(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm10;"
            "movapd 144(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm11;"
            "movapd 160(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm12;"
            "movapd 176(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm13;"
            "movapd 192(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm14;"
            "movapd 208(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm15;"

          :
          : "r"(aki), "r"(bkj)
          :
          );
        }
        __asm__ __volatile__ (
            "movapd   %%xmm2,  (%0);"
            "movapd   %%xmm3,  (%1);"
            "movapd   %%xmm4,  (%2);"
            "movapd   %%xmm5,  (%3);"
            "movapd   %%xmm6,  (%4);"
            "movapd   %%xmm7,  (%5);"
            "movapd   %%xmm8,  (%6);"
            "movapd   %%xmm9,  (%7);"
            "movapd   %%xmm10, (%8);"
            "movapd   %%xmm11, (%9);"
            "movapd   %%xmm12, (%10);"
            "movapd   %%xmm13, (%11);"
            "movapd   %%xmm14, (%12);"
            "movapd   %%xmm15, (%13);"
        :
        : "r"(cij0), "r"(cij1), "r"(cij2), "r"(cij3), "r"(cij4), "r"(cij5), "r"(cij6), "r"(cij7), "r"(cij8), "r"(cij9), "r"(cij10), "r"(cij11), "r"(cij12), "r"(cij13)
        :
        );
      }
      else
      {
        for (int i = ilo; i < dimi; ++i)
        {
          std::complex<double>* cij = c + i*dimj + j;
          *(cij) = std::complex<double>(0.0,0.0);
          for (long k = 0; k < dimk; ++k)
          {
            const std::complex<double>* aki = a + k*dimi + i;
            const double* bkj = b + k*dimj + j;
            *(cij) += *(aki) * *(bkj);
          }
        }
      }
    }
  }
}

void mTxmSCOTT(long dimi, long dimj, long dimk,
           std::complex<double>* restrict c, const std::complex<double>* a, const double* b) {

    const long jtile = 14;
    for (long i=0; i<dimi; ++i) {
    	std::complex<double>* cij = c + i*dimj;
    	for (long jlo=0; jlo<dimj; jlo+=jtile, cij+=jtile) {
    		int nj = dimj-jlo;
    		if (nj > jtile) nj = jtile;
    		switch (nj) {
    		case 14:
				{__asm__ __volatile__ (
						"pxor %xmm2,%xmm2;"
						"pxor %xmm3,%xmm3;"
						"pxor %xmm4,%xmm4;"
						"pxor %xmm5,%xmm5;"
						"pxor %xmm6,%xmm6;"
						"pxor %xmm7,%xmm7;"
						"pxor %xmm8,%xmm8;"
						"pxor %xmm9,%xmm9;"
						"pxor %xmm10,%xmm10;"
						"pxor %xmm11,%xmm11;"
						"pxor %xmm12,%xmm12;"
						"pxor %xmm13,%xmm13;"
						"pxor %xmm14,%xmm14;"
						"pxor %xmm15,%xmm15;"
				);
				const std::complex<double>* aki = a + i;
				const double* bkj = b + jlo;
				for (long k=0; k<dimk; ++k,aki+=dimi,bkj+=dimj) {
					__asm__ volatile(
						"movapd (%0), %%xmm0;"

						"movddup    (%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;"
						"movddup   8(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;"
						"movddup  16(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4;"
						"movddup  24(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm5;"
						"movddup  32(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm6;"
						"movddup  40(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm7;"
						"movddup  48(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm8;"
						"movddup  56(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm9;"
						"movddup  64(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm10;"
						"movddup  72(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm11;"
						"movddup  80(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm12;"
						"movddup  88(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm13;"
						"movddup  96(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm14;"
						"movddup 104(%1), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm15;"

					:
					: "r"(aki), "r"(bkj)
					:
					);
				}
				__asm__ __volatile__ (
						"movapd   %%xmm2,    (%0);"
						"movapd   %%xmm3,  16(%0);"
						"movapd   %%xmm4,  32(%0);"
						"movapd   %%xmm5,  48(%0);"
						"movapd   %%xmm6,  64(%0);"
						"movapd   %%xmm7,  80(%0);"
						"movapd   %%xmm8,  96(%0);"
						"movapd   %%xmm9, 112(%0);"
						"movapd  %%xmm10, 128(%0);"
						"movapd  %%xmm11, 144(%0);"
						"movapd  %%xmm12, 160(%0);"
						"movapd  %%xmm13, 176(%0);"
						"movapd  %%xmm14, 192(%0);"
						"movapd  %%xmm15, 208(%0);"
				:
				: "r"(cij)
				:
				);}
				break;

    		default:
    		    for (long j=jlo; j<dimj; ++j) {
					std::complex<double> cij(0.0,0.0);
					for (long k=0; k<dimk; ++k) {
						cij += a[k*dimi+i]*b[k*dimj+j];
					}
					c[i*dimj + j] = cij;
    		    }
    		    break;
    		}
        }
    }



}

int main(int argc, char** argv)
{
	int k = 28;

	Tensor< std::complex<double> > dc(k,k*k);
	Tensor< std::complex<double> > res(k,k*k), res2(k,k*k);
	Tensor< double > d(k,k);

	dc.fillrandom();
	d.fillrandom();

	mTxmq(k*k, k, k, res.ptr(), dc.ptr(), d.ptr());
	mtxmGooberSaysHey(k*k, k, k, res2.ptr(), dc.ptr(), d.ptr());
	//mTxmSCOTT(k*k, k, k, res2.ptr(), dc.ptr(), d.ptr());
	print("ERROR IS", (res-res2).normf());

//	double start = cpu_time();
//	for (int i = 0; i < 10000; ++i)
//	{
//		mTxmSCOTT(k*k, k, k, res.ptr(), dc.ptr(), d.ptr());
//	}
//	double tused = cpu_time() - start;
//	print("Cpu time used: ", tused, "flops/s: ", 2.0*k*k*k*k*1e4 / tused);

	return 0;
}

#else
#include <iostream>
int main() {
    std::cout << "only on x8664" << std::endl;
    return 0;
}
#endif
