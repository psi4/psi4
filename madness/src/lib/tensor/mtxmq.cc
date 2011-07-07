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
#include <madness_config.h>
#include <tensor/tensor.h>
#include <tensor/mtxmq.h>
#include <world/worldprofile.h>

// For x86-32/64 have assembly versions for double precision
// For x86-64 have assembly versions for complex double precision

#if defined(X86_32) || defined(X86_64)

#ifdef X86_64
extern "C" void mTxm26(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm24(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm22(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm20(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm18(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm16(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm14(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm12(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;
#endif // X86_64
extern "C" void mTxm10(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void mTxm8(long dimi, long dimj, long dimk,
                      double* c, const double* a, const double* b) ;

extern "C" void mTxm6(long dimi, long dimj, long dimk,
                      double* c, const double* a, const double* b) ;

extern "C" void mTxm4(long dimi, long dimj, long dimk,
                      double* c, const double* a, const double* b) ;

extern "C" void mTxm2(long dimi, long dimj, long dimk,
                      double* c, const double* a, const double* b) ;


#ifdef X86_64
extern "C" void TmTxm26(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm24(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm22(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm20(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm18(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm16(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm14(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm12(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;
#endif // X86_64
extern "C" void TmTxm10(long dimi, long dimj, long dimk,
                        double* c, const double* a, const double* b) ;

extern "C" void TmTxm8(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void TmTxm6(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void TmTxm4(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

extern "C" void TmTxm2(long dimi, long dimj, long dimk,
                       double* c, const double* a, const double* b) ;

namespace madness {

    template<>
    void mTxmq(const long dimi, const long dimj, const long dimk,
               double* restrict c, const double* a, const double* b) {
        PROFILE_BLOCK(mTxmq_double_asm);
        //std::cout << "IN DOUBLE ASM VERSION " << dimi << " " << dimj << " " << dimk << "\n";


        if (IS_ODD(dimi) || IS_ODD(dimj) || IS_ODD(dimk) ||
            IS_UNALIGNED(a) || IS_UNALIGNED(b) || IS_UNALIGNED(c)) {
            //std::cout << "slow\n";
            // CALL SLOW CODE
            for (long i=0; i<dimi; ++i,c+=dimj,++a) {
                for (long j=0; j<dimj; ++j) c[j] = 0.0;
                const double *ai = a;
                for (long k=0; k<dimk; ++k,ai+=dimi) {
                    double aki = *ai;
                    for (long j=0; j<dimj; ++j) {
                        c[j] += aki*b[k*dimj+j];
                    }
                }
            }
            return;
        }

        /*
           Choice is to unroll i or j
        */

#if   defined(AMD_QUADCORE_TUNE)
        bool test = dimj>=14 && dimj<=26;
#elif defined(OPTERON_TUNE)
        bool test = dimi <= dimj; /* Based on times from X86_64 Opteron ... an old one */
#elif defined(CORE_DUO_TUNE)
        bool test = true; /* Based on times from X86_32 Core Duo ... my old laptop */
#elif (defined(CORE2_TUNE) && defined(X86_32))
        bool test = false; /* Based on times from Core2 running in 32-bit mode ... a sad thing */
#elif (defined(CORE2_TUNE) && defined(X86_64))
        bool test = dimj > 12 || dimi <= dimj; /* Based on times from X86_64 Core2 */
#else
        bool test = dimj > 12 || dimi <= dimj; /* Based on times from X86_64 Core2 */
#endif
        if (test) {
            long nj = dimj;
            do {
#ifdef X86_64
                long numj = (nj>26) ? 26 : nj;
#else
                long numj = (nj>10) ? 10 : nj;
#endif

                switch (numj) {
#ifdef X86_64
                case 26:
                    TmTxm26(dimj, dimi, dimk, c, b, a) ;
                    break;

                case 24:
                    TmTxm24(dimj, dimi, dimk, c, b, a) ;
                    break;

                case 22:
                    TmTxm22(dimj, dimi, dimk, c, b, a) ;
                    break;

                case 20:
                    TmTxm20(dimj, dimi, dimk, c, b, a) ;
                    break;

                case 18:
                    TmTxm18(dimj, dimi, dimk, c, b, a) ;
                    break;

                case 16:
                    TmTxm16(dimj, dimi, dimk, c, b, a) ;
                    break;

                case 14:
                    TmTxm14(dimj, dimi, dimk, c, b, a) ;
                    break;

                case 12:
                    TmTxm12(dimj, dimi, dimk, c, b, a) ;
                    break;
#endif // X86_64

                case 10:
                    TmTxm10(dimj, dimi, dimk, c, b, a) ;
                    break;

                case 8:
                    TmTxm8(dimj, dimi, dimk, c, b, a) ;
                    break;

                case 6:
                    TmTxm6(dimj, dimi, dimk, c, b, a) ;
                    break;

                case 4:
                    TmTxm4(dimj, dimi, dimk, c, b, a) ;
                    break;

                case 2:
                    TmTxm2(dimj, dimi, dimk, c, b, a) ;
                    break;

                default:
                    throw "mtxmq_byj: should not be here";

                }
                nj -= numj;
                c += numj;
                b += numj;
            } while (nj);
        }
        else {
            long ni = dimi;
            do {
#ifdef X86_64
                long numi = (ni>26) ? 26 : ni;
#else
                long numi = (ni>10) ? 10 : ni;
#endif

                switch (numi) {
#ifdef X86_64
                case 26:
                    mTxm26(dimi, dimj, dimk, c, a, b) ;
                    break;

                case 24:
                    mTxm24(dimi, dimj, dimk, c, a, b) ;
                    break;

                case 22:
                    mTxm22(dimi, dimj, dimk, c, a, b) ;
                    break;

                case 20:
                    mTxm20(dimi, dimj, dimk, c, a, b) ;
                    break;

                case 18:
                    mTxm18(dimi, dimj, dimk, c, a, b) ;
                    break;

                case 16:
                    mTxm16(dimi, dimj, dimk, c, a, b) ;
                    break;

                case 14:
                    mTxm14(dimi, dimj, dimk, c, a, b) ;
                    break;

                case 12:
                    mTxm12(dimi, dimj, dimk, c, a, b) ;
                    break;
#endif // X86_64

                case 10:
                    mTxm10(dimi, dimj, dimk, c, a, b) ;
                    break;

                case 8:
                    mTxm8(dimi, dimj, dimk, c, a, b) ;
                    break;

                case 6:
                    mTxm6(dimi, dimj, dimk, c, a, b) ;
                    break;

                case 4:
                    mTxm4(dimi, dimj, dimk, c, a, b) ;
                    break;

                case 2:
                    mTxm2(dimi, dimj, dimk, c, a, b) ;
                    break;

                default:
                    throw "mtxmq: should not be here!";
                }
                ni -= numi;
                c += numi*dimj;
                a += numi;
            } while (ni);

        }
    }
}

#endif // defined(X86_32) || defined(X86_64)


/* all preprocessor ifdef-endif pairs are closed */


#if defined(X86_64)  && !defined(DISABLE_SSE3)
namespace madness {
    template <>
    void mTxmq(const long dimi, const long dimj, const long dimk,
               double_complex* restrict c, const double_complex* a, const double_complex* b) {

        PROFILE_BLOCK(mTxmq_complex_asm);
        const long dimi16 = dimi<<4;
        const long dimj16 = dimj<<4;

#define ZERO(c) "pxor " #c "," #c ";\n"

#ifdef AMD_QUADCORE_TUNE
#  define ENTRY(loop) "mov %q0,%%r9;  prefetcht0 (%%r9); mov %q1, %%r10; mov %q4,%%r11;.align 32;"#loop": "
#  define LOADA   "movaps (%%r9),%%xmm0; movaps %%xmm0, %%xmm1; shufpd $1,%%xmm1,%%xmm1;  mov %%r10,%%r8; add %q2,%%r9; add %q3,%%r10; prefetcht0 (%%r9);\n"
#  define DOIT(c) "movddup (%%r8),%%xmm2; mulpd %%xmm0,%%xmm2; addpd %%xmm2,"#c"; movddup 8(%%r8),%%xmm2; mulpd %%xmm1,%%xmm2; addsubpd %%xmm2,"#c"; \n"
#else
#  ifndef ON_A_MAC
#    define ENTRY(loop) "mov %q0,%%r9; mov %q1, %%r10; mov %q4,%%r11;.align 32;"#loop": "
#  else
#    define ENTRY(loop) "mov %q0,%%r9; mov %q1, %%r10; mov %q4,%%r11;"#loop": "
#  endif
#  define LOADA   "movddup  (%%r9), %%xmm0; mov %%r10,%%r8; movddup 8(%%r9), %%xmm1; add %q2,%%r9; add %q3,%%r10; prefetcht0 (%%r9);\n"
#  define DOIT(c) "movaps (%%r8),%%xmm2; movaps %%xmm2,%%xmm3; mulpd %%xmm0,%%xmm2; addpd %%xmm2,"#c"; shufpd $1,%%xmm3,%%xmm3; mulpd %%xmm1,%%xmm3; addsubpd %%xmm3, "#c"; \n"
#endif

// see comment in mtxmq_asm.S about movaps vs. movntps
#  define STORE(c) "movaps " #c ", (%%r8); add $16,%%r8;\n"

#define NEXT(loop) "sub $1,%%r11; jnz "#loop";"
#define INCB    "add $16,%%r8;\n"

        const long jtile = 12;
        const double_complex* asave = a;
        for (long jlo=0; jlo<dimj; jlo+=jtile,c+=jtile,b+=jtile) {
            long nj = std::min(dimj-jlo,jtile);
            double_complex* restrict ci = c;
            a = asave;
            for (long i=0; i<dimi; ++i,ci+=dimj,++a) {
                const double_complex *ai = a;
                const double_complex *bk = b;
                switch (nj) {
                case 1:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4)

                                          ENTRY(.KLOOP1)
                                          LOADA
                                          DOIT(%%xmm4)
                                          NEXT(.KLOOP1)

                                          "mov %q5, %%r8;\n"
                                          STORE(%%xmm4)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 2:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4)
                                          ZERO(%%xmm5)

                                          ENTRY(.KLOOP2)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)
                                          NEXT(.KLOOP2)

                                          "mov %q5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 3:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4)
                                          ZERO(%%xmm5)
                                          ZERO(%%xmm6)

                                          ENTRY(.KLOOP3)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)
                                          NEXT(.KLOOP3)

                                          "mov %q5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 4:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4)
                                          ZERO(%%xmm5)
                                          ZERO(%%xmm6)
                                          ZERO(%%xmm7)

                                          ENTRY(.KLOOP4)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)
                                          NEXT(.KLOOP4)

                                          "mov %q5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 5:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4)
                                          ZERO(%%xmm5)
                                          ZERO(%%xmm6)
                                          ZERO(%%xmm7)
                                          ZERO(%%xmm8)

                                          ENTRY(.KLOOP5)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)
                                          NEXT(.KLOOP5)

                                          "mov %q5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 6:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4)
                                          ZERO(%%xmm5)
                                          ZERO(%%xmm6)
                                          ZERO(%%xmm7)
                                          ZERO(%%xmm8)
                                          ZERO(%%xmm9)

                                          ENTRY(.KLOOP6)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)  INCB
                                          DOIT(%%xmm9)
                                          NEXT(.KLOOP6)

                                          "mov %q5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          STORE(%%xmm9)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );

                    break;

                case 7:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4)
                                          ZERO(%%xmm5)
                                          ZERO(%%xmm6)
                                          ZERO(%%xmm7)
                                          ZERO(%%xmm8)
                                          ZERO(%%xmm9)
                                          ZERO(%%xmm10)

                                          ENTRY(.KLOOP7)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)  INCB
                                          DOIT(%%xmm9)  INCB
                                          DOIT(%%xmm10)
                                          NEXT(.KLOOP7)

                                          "mov %q5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          STORE(%%xmm9)
                                          STORE(%%xmm10)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 8:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4)
                                          ZERO(%%xmm5)
                                          ZERO(%%xmm6)
                                          ZERO(%%xmm7)
                                          ZERO(%%xmm8)
                                          ZERO(%%xmm9)
                                          ZERO(%%xmm10)
                                          ZERO(%%xmm11)

                                          ENTRY(.KLOOP8)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)  INCB
                                          DOIT(%%xmm9)  INCB
                                          DOIT(%%xmm10) INCB
                                          DOIT(%%xmm11)
                                          NEXT(.KLOOP8)

                                          "mov %q5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          STORE(%%xmm9)
                                          STORE(%%xmm10)
                                          STORE(%%xmm11)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 9:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4)
                                          ZERO(%%xmm5)
                                          ZERO(%%xmm6)
                                          ZERO(%%xmm7)
                                          ZERO(%%xmm8)
                                          ZERO(%%xmm9)
                                          ZERO(%%xmm10)
                                          ZERO(%%xmm11)
                                          ZERO(%%xmm12)

                                          ENTRY(.KLOOP9)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)  INCB
                                          DOIT(%%xmm9)  INCB
                                          DOIT(%%xmm10) INCB
                                          DOIT(%%xmm11) INCB
                                          DOIT(%%xmm12)
                                          NEXT(.KLOOP9)

                                          "mov %q5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          STORE(%%xmm9)
                                          STORE(%%xmm10)
                                          STORE(%%xmm11)
                                          STORE(%%xmm12)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 10:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4)
                                          ZERO(%%xmm5)
                                          ZERO(%%xmm6)
                                          ZERO(%%xmm7)
                                          ZERO(%%xmm8)
                                          ZERO(%%xmm9)
                                          ZERO(%%xmm10)
                                          ZERO(%%xmm11)
                                          ZERO(%%xmm12)
                                          ZERO(%%xmm13)

                                          ENTRY(.KLOOP10)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)  INCB
                                          DOIT(%%xmm9)  INCB
                                          DOIT(%%xmm10) INCB
                                          DOIT(%%xmm11) INCB
                                          DOIT(%%xmm12) INCB
                                          DOIT(%%xmm13)
                                          NEXT(.KLOOP10)

                                          "mov %q5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          STORE(%%xmm9)
                                          STORE(%%xmm10)
                                          STORE(%%xmm11)
                                          STORE(%%xmm12)
                                          STORE(%%xmm13)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 11:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4)
                                          ZERO(%%xmm5)
                                          ZERO(%%xmm6)
                                          ZERO(%%xmm7)
                                          ZERO(%%xmm8)
                                          ZERO(%%xmm9)
                                          ZERO(%%xmm10)
                                          ZERO(%%xmm11)
                                          ZERO(%%xmm12)
                                          ZERO(%%xmm13)
                                          ZERO(%%xmm14)

                                          ENTRY(.KLOOP11)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)  INCB
                                          DOIT(%%xmm9)  INCB
                                          DOIT(%%xmm10) INCB
                                          DOIT(%%xmm11) INCB
                                          DOIT(%%xmm12) INCB
                                          DOIT(%%xmm13) INCB
                                          DOIT(%%xmm14)
                                          NEXT(.KLOOP11)

                                          "mov %q5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          STORE(%%xmm9)
                                          STORE(%%xmm10)
                                          STORE(%%xmm11)
                                          STORE(%%xmm12)
                                          STORE(%%xmm13)
                                          STORE(%%xmm14)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;

                case 12:
                    __asm__ __volatile__ (
                                          ZERO(%%xmm4)
                                          ZERO(%%xmm5)
                                          ZERO(%%xmm6)
                                          ZERO(%%xmm7)
                                          ZERO(%%xmm8)
                                          ZERO(%%xmm9)
                                          ZERO(%%xmm10)
                                          ZERO(%%xmm11)
                                          ZERO(%%xmm12)
                                          ZERO(%%xmm13)
                                          ZERO(%%xmm14)
                                          ZERO(%%xmm15)

                                          ENTRY(.KLOOP12)
                                          LOADA
                                          DOIT(%%xmm4)  INCB
                                          DOIT(%%xmm5)  INCB
                                          DOIT(%%xmm6)  INCB
                                          DOIT(%%xmm7)  INCB
                                          DOIT(%%xmm8)  INCB
                                          DOIT(%%xmm9)  INCB
                                          DOIT(%%xmm10) INCB
                                          DOIT(%%xmm11) INCB
                                          DOIT(%%xmm12) INCB
                                          DOIT(%%xmm13) INCB
                                          DOIT(%%xmm14) INCB
                                          DOIT(%%xmm15)
                                          NEXT(.KLOOP12)

                                          "mov %q5, %%r8;\n"
                                          STORE(%%xmm4)
                                          STORE(%%xmm5)
                                          STORE(%%xmm6)
                                          STORE(%%xmm7)
                                          STORE(%%xmm8)
                                          STORE(%%xmm9)
                                          STORE(%%xmm10)
                                          STORE(%%xmm11)
                                          STORE(%%xmm12)
                                          STORE(%%xmm13)
                                          STORE(%%xmm14)
                                          STORE(%%xmm15)
                                          :
                                          : "r"(ai),"r"(bk),"r"(dimi16),"r"(dimj16),"r"(dimk), "r"(ci)
                                          : "r8", "r9", "r10", "r11", "memory"
                                          );
                    break;
                }
            }
        }
    }

#ifndef __INTEL_COMPILER
    template <>
    void mTxmq(const long dimi, const long dimj, const long dimk,
               double_complex* restrict c, const double_complex* a, const double* b)
    {
      const long itile = 14;
      for (long ilo = 0; ilo < dimi; ilo += itile, a+=itile, c+=itile*dimj)
      {
        long ni = dimi - ilo;
        ni = (ni >= itile) ? itile : ni;
        if (ni == 1)
        {
          for (long j = 0; j < dimj; ++j)
          {
            __asm__ volatile
            (
                // save registers to be 'clobbered'
                "push %0; push %1; push %4; push %5;\n "
                // zero out mmx registers
                "pxor %%xmm2,%%xmm2;\n"

                "0:\n "
                // load the 'b' part into %1 and update pointer
                "movddup   (%1), %%xmm0; add %3,%1;\n"
                // begin i-tile
                "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2; add %2,%0;\n"
                "sub $1,%4; jnz 0b;\n"

                "movapd   %%xmm2, (%5);\n"

                "pop %5; pop %4; pop %1; pop %0;\n"

                :
                : "r"(a), "r"(b + j), "r"(dimi<<4), "r"(dimj<<3), "r"(dimk), "r"(c + j), "r"(dimj<<4)
                :
            );
          }
        }
        else if (ni == 2)
        {
          for (long j = 0; j < dimj; ++j)
          {
            __asm__ volatile
            (
                // save registers to be 'clobbered'
                "push %0; push %1; push %4; push %5;\n "
                // zero out mmx registers
                "pxor %%xmm2,%%xmm2;\n"
                "pxor %%xmm3,%%xmm3;\n"

                "0:\n "
                // load the 'b' part into %1 and update pointer
                "movddup   (%1), %%xmm0; add %3,%1;\n"
                // begin i-tile
                "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;\n"
                "movapd  16(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3; add %2,%0;\n"
                "sub $1,%4; jnz 0b;\n"

                "movapd   %%xmm2, (%5); add %6,%5;\n"
                "movapd   %%xmm3, (%5);\n"

                "pop %5; pop %4; pop %1; pop %0;\n"

                :
                : "r"(a), "r"(b + j), "r"(dimi<<4), "r"(dimj<<3), "r"(dimk), "r"(c + j), "r"(dimj<<4)
                :
            );
          }
        }
        else if (ni == 3)
        {
          for (long j = 0; j < dimj; ++j)
          {
            __asm__ volatile
            (
                // save registers to be 'clobbered'
                "push %0; push %1; push %4; push %5;\n "
                // zero out mmx registers
                "pxor %%xmm2,%%xmm2;\n"
                "pxor %%xmm3,%%xmm3;\n"
                "pxor %%xmm4,%%xmm4;\n"

                "0:\n "
                // load the 'b' part into %1 and update pointer
                "movddup   (%1), %%xmm0; add %3,%1;\n"
                // begin i-tile
                "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;\n"
                "movapd  16(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;\n"
                "movapd  32(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4; add %2,%0;\n"
                "sub $1,%4; jnz 0b;\n"

                "movapd   %%xmm2, (%5); add %6,%5;\n"
                "movapd   %%xmm3, (%5); add %6,%5;\n"
                "movapd   %%xmm4, (%5);\n"

                "pop %5; pop %4; pop %1; pop %0;\n"

                :
                : "r"(a), "r"(b + j), "r"(dimi<<4), "r"(dimj<<3), "r"(dimk), "r"(c + j), "r"(dimj<<4)
                :
            );
          }
        }
        else if (ni == 4)
        {
          for (long j = 0; j < dimj; ++j)
          {
            __asm__ volatile
            (
                // save registers to be 'clobbered'
                "push %0; push %1; push %4; push %5;\n "
                // zero out mmx registers
                "pxor %%xmm2,%%xmm2;\n"
                "pxor %%xmm3,%%xmm3;\n"
                "pxor %%xmm4,%%xmm4;\n"
                "pxor %%xmm5,%%xmm5;\n"

                "0:\n "
                // load the 'b' part into %1 and update pointer
                "movddup   (%1), %%xmm0; add %3,%1;\n"
                // begin i-tile
                "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;\n"
                "movapd  16(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;\n"
                "movapd  32(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4;\n"
                "movapd  48(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm5; add %2,%0;\n"
                "sub $1,%4; jnz 0b;\n"

                "movapd   %%xmm2, (%5); add %6,%5;\n"
                "movapd   %%xmm3, (%5); add %6,%5;\n"
                "movapd   %%xmm4, (%5); add %6,%5;\n"
                "movapd   %%xmm5, (%5);\n"

                "pop %5; pop %4; pop %1; pop %0;\n"

                :
                : "r"(a), "r"(b + j), "r"(dimi<<4), "r"(dimj<<3), "r"(dimk), "r"(c + j), "r"(dimj<<4)
                :
            );
          }
        }
        if (ni == 5)
        {
          for (long j = 0; j < dimj; ++j)
          {
            __asm__ volatile
            (
                // save registers to be 'clobbered'
                "push %0; push %1; push %4; push %5;\n "
                // zero out mmx registers
                "pxor %%xmm2,%%xmm2;\n"
                "pxor %%xmm3,%%xmm3;\n"
                "pxor %%xmm4,%%xmm4;\n"
                "pxor %%xmm5,%%xmm5;\n"
                "pxor %%xmm6,%%xmm6;\n"

                "0:\n "
                // load the 'b' part into %1 and update pointer
                "movddup   (%1), %%xmm0; add %3,%1;\n"
                // begin i-tile
                "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;\n"
                "movapd  16(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;\n"
                "movapd  32(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4;\n"
                "movapd  48(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm5;\n"
                "movapd  64(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm6; add %2,%0; \n"
                "sub $1,%4; jnz 0b;\n"

                "movapd   %%xmm2, (%5); add %6,%5;\n"
                "movapd   %%xmm3, (%5); add %6,%5;\n"
                "movapd   %%xmm4, (%5); add %6,%5;\n"
                "movapd   %%xmm5, (%5); add %6,%5;\n"
                "movapd   %%xmm6, (%5);\n"

                "pop %5; pop %4; pop %1; pop %0;\n"

                :
                : "r"(a), "r"(b + j), "r"(dimi<<4), "r"(dimj<<3), "r"(dimk), "r"(c + j), "r"(dimj<<4)
                :
            );
          }
        }
        else if (ni == 6)
        {
          for (long j = 0; j < dimj; ++j)
          {
            __asm__ volatile
            (
                // save registers to be 'clobbered'
                "push %0; push %1; push %4; push %5;\n "
                // zero out mmx registers
                "pxor %%xmm2,%%xmm2;\n"
                "pxor %%xmm3,%%xmm3;\n"
                "pxor %%xmm4,%%xmm4;\n"
                "pxor %%xmm5,%%xmm5;\n"
                "pxor %%xmm6,%%xmm6;\n"
                "pxor %%xmm7,%%xmm7;\n"

                "0:\n "
                // load the 'b' part into %1 and update pointer
                "movddup   (%1), %%xmm0; add %3,%1;\n"
                // begin i-tile
                "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;\n"
                "movapd  16(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;\n"
                "movapd  32(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4;\n"
                "movapd  48(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm5;\n"
                "movapd  64(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm6;\n"
                "movapd  80(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm7; add %2,%0; \n"
                "sub $1,%4; jnz 0b;\n"

                "movapd   %%xmm2, (%5); add %6,%5;\n"
                "movapd   %%xmm3, (%5); add %6,%5;\n"
                "movapd   %%xmm4, (%5); add %6,%5;\n"
                "movapd   %%xmm5, (%5); add %6,%5;\n"
                "movapd   %%xmm6, (%5); add %6,%5;\n"
                "movapd   %%xmm7, (%5); \n"

                "pop %5; pop %4; pop %1; pop %0;\n"

                :
                : "r"(a), "r"(b + j), "r"(dimi<<4), "r"(dimj<<3), "r"(dimk), "r"(c + j), "r"(dimj<<4)
                :
            );
          }
        }
        else if (ni == 7)
         {
           for (long j = 0; j < dimj; ++j)
           {
             __asm__ volatile
             (
               // save registers to be 'clobbered'
               "push %0; push %1; push %4; push %5;\n "
               // zero out mmx registers
               "pxor %%xmm2,%%xmm2;\n"
               "pxor %%xmm3,%%xmm3;\n"
               "pxor %%xmm4,%%xmm4;\n"
               "pxor %%xmm5,%%xmm5;\n"
               "pxor %%xmm6,%%xmm6;\n"
               "pxor %%xmm7,%%xmm7;\n"
               "pxor %%xmm8,%%xmm8;\n"

              "0:\n "
              // load the 'b' part into %1 and update pointer
              "movddup   (%1), %%xmm0; add %3,%1;\n"
              // begin i-tile
              "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;\n"
              "movapd  16(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;\n"
              "movapd  32(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4;\n"
              "movapd  48(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm5;\n"
              "movapd  64(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm6;\n"
              "movapd  80(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm7;\n"
              "movapd  96(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm8; add %2,%0; \n"
              "sub $1,%4; jnz 0b;\n"

              "movapd   %%xmm2, (%5); add %6,%5;\n"
              "movapd   %%xmm3, (%5); add %6,%5;\n"
              "movapd   %%xmm4, (%5); add %6,%5;\n"
              "movapd   %%xmm5, (%5); add %6,%5;\n"
              "movapd   %%xmm6, (%5); add %6,%5;\n"
              "movapd   %%xmm7, (%5); add %6,%5;\n"
              "movapd   %%xmm8, (%5);\n"

              "pop %5; pop %4; pop %1; pop %0;\n"

              :
              : "r"(a), "r"(b + j), "r"(dimi<<4), "r"(dimj<<3), "r"(dimk), "r"(c + j), "r"(dimj<<4)
              :
             );
           }
         }
       else if (ni == 8)
        {
          for (long j = 0; j < dimj; ++j)
          {
            __asm__ volatile
            (
              // save registers to be 'clobbered'
              "push %0; push %1; push %4; push %5;\n "
              // zero out mmx registers
              "pxor %%xmm2,%%xmm2;\n"
              "pxor %%xmm3,%%xmm3;\n"
              "pxor %%xmm4,%%xmm4;\n"
              "pxor %%xmm5,%%xmm5;\n"
              "pxor %%xmm6,%%xmm6;\n"
              "pxor %%xmm7,%%xmm7;\n"
              "pxor %%xmm8,%%xmm8;\n"
              "pxor %%xmm9,%%xmm9;\n"

             "0:\n "
             // load the 'b' part into %1 and update pointer
             "movddup   (%1), %%xmm0; add %3,%1;\n"
             // begin i-tile
             "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;\n"
             "movapd  16(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;\n"
             "movapd  32(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4;\n"
             "movapd  48(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm5;\n"
             "movapd  64(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm6;\n"
             "movapd  80(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm7;\n"
             "movapd  96(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm8;\n"
             "movapd 112(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm9; add %2,%0;\n"
             "sub $1,%4; jnz 0b;\n"

             "movapd   %%xmm2, (%5); add %6,%5;\n"
             "movapd   %%xmm3, (%5); add %6,%5;\n"
             "movapd   %%xmm4, (%5); add %6,%5;\n"
             "movapd   %%xmm5, (%5); add %6,%5;\n"
             "movapd   %%xmm6, (%5); add %6,%5;\n"
             "movapd   %%xmm7, (%5); add %6,%5;\n"
             "movapd   %%xmm8, (%5); add %6,%5;\n"
             "movapd   %%xmm9, (%5);\n"

             "pop %5; pop %4; pop %1; pop %0;\n"

             :
             : "r"(a), "r"(b + j), "r"(dimi<<4), "r"(dimj<<3), "r"(dimk), "r"(c + j), "r"(dimj<<4)
             :
            );
          }
        }
        else if (ni == 9)
        {
          for (long j = 0; j < dimj; ++j)
          {
            __asm__ volatile
            (
              // save registers to be 'clobbered'
              "push %0; push %1; push %4; push %5;\n "
              // zero out mmx registers
              "pxor %%xmm2,%%xmm2;\n"
              "pxor %%xmm3,%%xmm3;\n"
              "pxor %%xmm4,%%xmm4;\n"
              "pxor %%xmm5,%%xmm5;\n"
              "pxor %%xmm6,%%xmm6;\n"
              "pxor %%xmm7,%%xmm7;\n"
              "pxor %%xmm8,%%xmm8;\n"
              "pxor %%xmm9,%%xmm9;\n"
              "pxor %%xmm10,%%xmm10;\n"

             "0:\n "
             // load the 'b' part into %1 and update pointer
             "movddup   (%1), %%xmm0; add %3,%1;\n"
             // begin i-tile
             "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;\n"
             "movapd  16(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;\n"
             "movapd  32(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4;\n"
             "movapd  48(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm5;\n"
             "movapd  64(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm6;\n"
             "movapd  80(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm7;\n"
             "movapd  96(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm8;\n"
             "movapd 112(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm9;\n"
             "movapd 128(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm10; add %2,%0;\n"
             "sub $1,%4; jnz 0b;\n"

             "movapd   %%xmm2, (%5); add %6,%5;\n"
             "movapd   %%xmm3, (%5); add %6,%5;\n"
             "movapd   %%xmm4, (%5); add %6,%5;\n"
             "movapd   %%xmm5, (%5); add %6,%5;\n"
             "movapd   %%xmm6, (%5); add %6,%5;\n"
             "movapd   %%xmm7, (%5); add %6,%5;\n"
             "movapd   %%xmm8, (%5); add %6,%5;\n"
             "movapd   %%xmm9, (%5); add %6,%5;\n"
             "movapd  %%xmm10, (%5);\n"

             "pop %5; pop %4; pop %1; pop %0;\n"

             :
             : "r"(a), "r"(b + j), "r"(dimi<<4), "r"(dimj<<3), "r"(dimk), "r"(c + j), "r"(dimj<<4)
             :
            );
          }
        }
        if (ni == 10)
        {
          for (long j = 0; j < dimj; ++j)
          {
            __asm__ volatile
            (
              // save registers to be 'clobbered'
              "push %0; push %1; push %4; push %5;\n "
              // zero out mmx registers
              "pxor %%xmm2,%%xmm2;\n"
              "pxor %%xmm3,%%xmm3;\n"
              "pxor %%xmm4,%%xmm4;\n"
              "pxor %%xmm5,%%xmm5;\n"
              "pxor %%xmm6,%%xmm6;\n"
              "pxor %%xmm7,%%xmm7;\n"
              "pxor %%xmm8,%%xmm8;\n"
              "pxor %%xmm9,%%xmm9;\n"
              "pxor %%xmm10,%%xmm10;\n"
              "pxor %%xmm11,%%xmm11;\n"

             "0:\n "
             // load the 'b' part into %1 and update pointer
             "movddup   (%1), %%xmm0; add %3,%1;\n"
             // begin i-tile
             "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;\n"
             "movapd  16(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;\n"
             "movapd  32(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4;\n"
             "movapd  48(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm5;\n"
             "movapd  64(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm6;\n"
             "movapd  80(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm7;\n"
             "movapd  96(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm8;\n"
             "movapd 112(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm9;\n"
             "movapd 128(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm10;\n"
             "movapd 144(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm11; add %2,%0;\n"
             "sub $1,%4; jnz 0b;\n"

             "movapd   %%xmm2, (%5); add %6,%5;\n"
             "movapd   %%xmm3, (%5); add %6,%5;\n"
             "movapd   %%xmm4, (%5); add %6,%5;\n"
             "movapd   %%xmm5, (%5); add %6,%5;\n"
             "movapd   %%xmm6, (%5); add %6,%5;\n"
             "movapd   %%xmm7, (%5); add %6,%5;\n"
             "movapd   %%xmm8, (%5); add %6,%5;\n"
             "movapd   %%xmm9, (%5); add %6,%5;\n"
             "movapd  %%xmm10, (%5); add %6,%5;\n"
             "movapd  %%xmm11, (%5);\n"

             "pop %5; pop %4; pop %1; pop %0;\n"

             :
             : "r"(a), "r"(b + j), "r"(dimi<<4), "r"(dimj<<3), "r"(dimk), "r"(c + j), "r"(dimj<<4)
             :
            );
          }
        }
        else if (ni == 11)
        {
          for (long j = 0; j < dimj; ++j)
          {
            __asm__ volatile
            (
              // save registers to be 'clobbered'
              "push %0; push %1; push %4; push %5;\n "
              // zero out mmx registers
              "pxor %%xmm2,%%xmm2;\n"
              "pxor %%xmm3,%%xmm3;\n"
              "pxor %%xmm4,%%xmm4;\n"
              "pxor %%xmm5,%%xmm5;\n"
              "pxor %%xmm6,%%xmm6;\n"
              "pxor %%xmm7,%%xmm7;\n"
              "pxor %%xmm8,%%xmm8;\n"
              "pxor %%xmm9,%%xmm9;\n"
              "pxor %%xmm10,%%xmm10;\n"
              "pxor %%xmm11,%%xmm11;\n"
              "pxor %%xmm12,%%xmm12;\n"

             "0:\n "
             // load the 'b' part into %1 and update pointer
             "movddup   (%1), %%xmm0; add %3,%1;\n"
             // begin i-tile
             "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;\n"
             "movapd  16(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;\n"
             "movapd  32(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4;\n"
             "movapd  48(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm5;\n"
             "movapd  64(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm6;\n"
             "movapd  80(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm7;\n"
             "movapd  96(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm8;\n"
             "movapd 112(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm9;\n"
             "movapd 128(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm10;\n"
             "movapd 144(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm11;\n"
             "movapd 160(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm12; add %2,%0;\n"
             "sub $1,%4; jnz 0b;\n"

             "movapd   %%xmm2, (%5); add %6,%5;\n"
             "movapd   %%xmm3, (%5); add %6,%5;\n"
             "movapd   %%xmm4, (%5); add %6,%5;\n"
             "movapd   %%xmm5, (%5); add %6,%5;\n"
             "movapd   %%xmm6, (%5); add %6,%5;\n"
             "movapd   %%xmm7, (%5); add %6,%5;\n"
             "movapd   %%xmm8, (%5); add %6,%5;\n"
             "movapd   %%xmm9, (%5); add %6,%5;\n"
             "movapd  %%xmm10, (%5); add %6,%5;\n"
             "movapd  %%xmm11, (%5); add %6,%5;\n"
             "movapd  %%xmm12, (%5);\n"

             "pop %5; pop %4; pop %1; pop %0;\n"

             :
             : "r"(a), "r"(b + j), "r"(dimi<<4), "r"(dimj<<3), "r"(dimk), "r"(c + j), "r"(dimj<<4)
             :
            );
          }
        }
        else if (ni == 12)
        {
          for (long j = 0; j < dimj; ++j)
          {
            __asm__ volatile
            (
              // save registers to be 'clobbered'
              "push %0; push %1; push %4; push %5;\n "
              // zero out mmx registers
              "pxor %%xmm2,%%xmm2;\n"
              "pxor %%xmm3,%%xmm3;\n"
              "pxor %%xmm4,%%xmm4;\n"
              "pxor %%xmm5,%%xmm5;\n"
              "pxor %%xmm6,%%xmm6;\n"
              "pxor %%xmm7,%%xmm7;\n"
              "pxor %%xmm8,%%xmm8;\n"
              "pxor %%xmm9,%%xmm9;\n"
              "pxor %%xmm10,%%xmm10;\n"
              "pxor %%xmm11,%%xmm11;\n"
              "pxor %%xmm12,%%xmm12;\n"
              "pxor %%xmm13,%%xmm13;\n"

             "0:\n "
             // load the 'b' part into %1 and update pointer
             "movddup   (%1), %%xmm0; add %3,%1;\n"
             // begin i-tile
             "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;\n"
             "movapd  16(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;\n"
             "movapd  32(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4;\n"
             "movapd  48(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm5;\n"
             "movapd  64(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm6;\n"
             "movapd  80(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm7;\n"
             "movapd  96(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm8;\n"
             "movapd 112(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm9;\n"
             "movapd 128(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm10;\n"
             "movapd 144(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm11;\n"
             "movapd 160(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm12;\n"
             "movapd 176(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm13; add %2,%0;\n"
             "sub $1,%4; jnz 0b;\n"

             "movapd   %%xmm2, (%5); add %6,%5;\n"
             "movapd   %%xmm3, (%5); add %6,%5;\n"
             "movapd   %%xmm4, (%5); add %6,%5;\n"
             "movapd   %%xmm5, (%5); add %6,%5;\n"
             "movapd   %%xmm6, (%5); add %6,%5;\n"
             "movapd   %%xmm7, (%5); add %6,%5;\n"
             "movapd   %%xmm8, (%5); add %6,%5;\n"
             "movapd   %%xmm9, (%5); add %6,%5;\n"
             "movapd  %%xmm10, (%5); add %6,%5;\n"
             "movapd  %%xmm11, (%5); add %6,%5;\n"
             "movapd  %%xmm12, (%5); add %6,%5;\n"
             "movapd  %%xmm13, (%5);\n"

             "pop %5; pop %4; pop %1; pop %0;\n"

             :
             : "r"(a), "r"(b + j), "r"(dimi<<4), "r"(dimj<<3), "r"(dimk), "r"(c + j), "r"(dimj<<4)
             :
            );
          }
        }
        else if (ni == 13)
        {
          for (long j = 0; j < dimj; ++j)
          {
            __asm__ volatile
            (
              // save registers to be 'clobbered'
              "push %0; push %1; push %4; push %5;\n "
              // zero out mmx registers
              "pxor %%xmm2,%%xmm2;\n"
              "pxor %%xmm3,%%xmm3;\n"
              "pxor %%xmm4,%%xmm4;\n"
              "pxor %%xmm5,%%xmm5;\n"
              "pxor %%xmm6,%%xmm6;\n"
              "pxor %%xmm7,%%xmm7;\n"
              "pxor %%xmm8,%%xmm8;\n"
              "pxor %%xmm9,%%xmm9;\n"
              "pxor %%xmm10,%%xmm10;\n"
              "pxor %%xmm11,%%xmm11;\n"
              "pxor %%xmm12,%%xmm12;\n"
              "pxor %%xmm13,%%xmm13;\n"
              "pxor %%xmm14,%%xmm14;\n"

             "0:\n "
             // load the 'b' part into %1 and update pointer
             "movddup   (%1), %%xmm0; add %3,%1;\n"
             // begin i-tile
             "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;\n"
             "movapd  16(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;\n"
             "movapd  32(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4;\n"
             "movapd  48(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm5;\n"
             "movapd  64(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm6;\n"
             "movapd  80(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm7;\n"
             "movapd  96(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm8;\n"
             "movapd 112(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm9;\n"
             "movapd 128(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm10;\n"
             "movapd 144(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm11;\n"
             "movapd 160(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm12;\n"
             "movapd 176(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm13;\n"
             "movapd 192(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm14; add %2,%0;\n"
             "sub $1,%4; jnz 0b;\n"

             "movapd   %%xmm2, (%5); add %6,%5;\n"
             "movapd   %%xmm3, (%5); add %6,%5;\n"
             "movapd   %%xmm4, (%5); add %6,%5;\n"
             "movapd   %%xmm5, (%5); add %6,%5;\n"
             "movapd   %%xmm6, (%5); add %6,%5;\n"
             "movapd   %%xmm7, (%5); add %6,%5;\n"
             "movapd   %%xmm8, (%5); add %6,%5;\n"
             "movapd   %%xmm9, (%5); add %6,%5;\n"
             "movapd  %%xmm10, (%5); add %6,%5;\n"
             "movapd  %%xmm11, (%5); add %6,%5;\n"
             "movapd  %%xmm12, (%5); add %6,%5;\n"
             "movapd  %%xmm13, (%5); add %6,%5;\n"
             "movapd  %%xmm14, (%5);"

             "pop %5; pop %4; pop %1; pop %0;\n"

             :
             : "r"(a), "r"(b + j), "r"(dimi<<4), "r"(dimj<<3), "r"(dimk), "r"(c + j), "r"(dimj<<4)
             :
            );
          }
        }
        else if (ni == 14)
        {
          for (long j = 0; j < dimj; ++j)
          {
            __asm__ volatile
            (
              // save registers to be 'clobbered'
              "push %0; push %1; push %4; push %5;\n "
              // zero out mmx registers
              "pxor %%xmm2,%%xmm2;\n"
              "pxor %%xmm3,%%xmm3;\n"
              "pxor %%xmm4,%%xmm4;\n"
              "pxor %%xmm5,%%xmm5;\n"
              "pxor %%xmm6,%%xmm6;\n"
              "pxor %%xmm7,%%xmm7;\n"
              "pxor %%xmm8,%%xmm8;\n"
              "pxor %%xmm9,%%xmm9;\n"
              "pxor %%xmm10,%%xmm10;\n"
              "pxor %%xmm11,%%xmm11;\n"
              "pxor %%xmm12,%%xmm12;\n"
              "pxor %%xmm13,%%xmm13;\n"
              "pxor %%xmm14,%%xmm14;\n"
              "pxor %%xmm15,%%xmm15;\n"

             "0:\n "
             // load the 'b' part into %1
             "movddup   (%1), %%xmm0; add %3,%1;\n"
             // begin i-loop
             "movapd    (%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm2;\n"
             "movapd  16(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm3;\n"
             "movapd  32(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm4;\n"
             "movapd  48(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm5;\n"
             "movapd  64(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm6;\n"
             "movapd  80(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm7;\n"
             "movapd  96(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm8;\n"
             "movapd 112(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm9;\n"
             "movapd 128(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm10;\n"
             "movapd 144(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm11;\n"
             "movapd 160(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm12;\n"
             "movapd 176(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm13;\n"
             "movapd 192(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm14;\n"
             "movapd 208(%0), %%xmm1; mulpd %%xmm0, %%xmm1; addpd %%xmm1,%%xmm15; add %2,%0;\n"
             "sub $1,%4; jnz 0b;\n"

             "movapd   %%xmm2, (%5); add %6,%5;\n"
             "movapd   %%xmm3, (%5); add %6,%5;\n"
             "movapd   %%xmm4, (%5); add %6,%5;\n"
             "movapd   %%xmm5, (%5); add %6,%5;\n"
             "movapd   %%xmm6, (%5); add %6,%5;\n"
             "movapd   %%xmm7, (%5); add %6,%5;\n"
             "movapd   %%xmm8, (%5); add %6,%5;\n"
             "movapd   %%xmm9, (%5); add %6,%5;\n"
             "movapd  %%xmm10, (%5); add %6,%5;\n"
             "movapd  %%xmm11, (%5); add %6,%5;\n"
             "movapd  %%xmm12, (%5); add %6,%5;\n"
             "movapd  %%xmm13, (%5); add %6,%5;\n"
             "movapd  %%xmm14, (%5); add %6,%5;\n"
             "movapd  %%xmm15, (%5);\n"

             "pop %5; pop %4; pop %1; pop %0;\n"

             :
             : "r"(a), "r"(b + j), "r"(dimi<<4), "r"(dimj<<3), "r"(dimk), "r"(c + j), "r"(dimj<<4)
             :
            );
          }
        }
      }
    }
#endif // __INTEL_COMPILER
}
#endif // defined(X86_64)  && !defined(DISABLE_SSE3)





