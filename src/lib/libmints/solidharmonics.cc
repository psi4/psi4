#include <cstdlib>
#include <cstring>
#include <cmath>
#include <stdint.h>

#include "matrix.h"
#include "integral.h"

#include <psiconfig.h>

using namespace psi;

uint64_t binomial(int n, int c1)
{
    uint64_t num = 1;
    uint64_t den = 1;
    int c2 = n - c1;
    int i;
    for (i=c2+1; i<=n; i++) {
        num *= i;
    }
    for (i=2; i<=c1; i++) {
        den *= i;
    }
    return num/den;
}

uint64_t fact(int n)
{
    uint64_t r = 1;
    for (int i=2; i<=n; i++) {
        r *= i;
    }
    return r;
}

uint64_t factoverfact(int nnum,int nden)
{
    uint64_t r = 1;
    for (int i=nden+1; i<=nnum; i++) {
        r *= i;
    }
    return r;
}

uint64_t factfact(int n)
{
    uint64_t result;
    int i;

    result = 1;
    if (n&1) {
        for (i=3; i<=n; i+=2) {
            result *= i;
        }
    }
    else {
        for (i=2; i<=n; i+=2) {
            result *= i;
        }
    }
    return result;
}

void reduce(uint64_t &num, uint64_t &den)
{
    if (num > den) {
        for (uint64_t i=2; i<=den;) {
            if (num%i == 0UL && den%i == 0UL) {
                num /= i;
                den /= i;
            }
            else i++;
        }
    }
    else {
        for (uint64_t i=2; i<=num;) {
            if (num%i == 0UL && den%i == 0UL) {
                num /= i;
                den /= i;
            }
            else i++;
        }
    }
}

uint64_t powll(uint64_t n, unsigned long p)
{
    uint64_t result = 1;
    for (unsigned long i=0; i<p; i++) result *= n;
    return result;
}

// there ordering here is arbitrary and doesn't have to match the
// basis set ordering
static inline int ncart(int l) { return (l>=0)?((((l)+2)*((l)+1))>>1):0; }
static inline int npure(int l) { return 2*l+1; }
static inline int icart(int a, int b, int c)
{
  return (((((a+b+c+1)<<1)-a)*(a+1))>>1)-b-1;
}
static inline int ipure(int l, int m) { return m<0?2*-m:(m==0?0:2*m-1); }

void solidharmcontrib(int sign,
                      const uint64_t &bin,const uint64_t &den,
                      uint64_t norm2num,uint64_t norm2den,
                      int r2,int x,int y,int z,
                      SimpleMatrix &coefmat, int pureindex)
{
    if (r2>0) {
        solidharmcontrib(sign,bin,den,norm2num,norm2den,r2-1,x+2,y,z,
                         coefmat,pureindex);
        solidharmcontrib(sign,bin,den,norm2num,norm2den,r2-1,x,y+2,z,
                         coefmat,pureindex);
        solidharmcontrib(sign,bin,den,norm2num,norm2den,r2-1,x,y,z+2,
                         coefmat,pureindex);
    }
    else {
        double coef = sign*double(bin)/double(den);
        double norm = sqrt(double(norm2num)/double(norm2den));
        coefmat.add(icart(x,y,z), pureindex, coef*norm);
    }
}

// l is the total angular momentum
// m is the z component
// r2 is the number of factors of r^2 that are included
void solidharm(unsigned int l, int m, unsigned int r2, SimpleMatrix& coefmat)
{
//    printf("in solidharm(unsigned int l, int m, unsigned int r2, RefSCMatrix coefmat\n");
//    printf("l = %d, m = %d, r2 = %d\n", l, m, r2);

    int pureindex = ipure(l,m);
//    printf("pureindex = %d\n", pureindex);
    for (unsigned int i=1; i<=r2; i++) pureindex += npure(l+2*i);
//    printf("pureindex = %d\n", pureindex);

    unsigned int absm = abs(m);

    // this overflows 32bits for l=9
    uint64_t norm2num = factoverfact(l+absm,l);
    uint64_t norm2den = factoverfact(l,l-absm);
    reduce(norm2num,norm2den);
    norm2num *= fact(absm);
    norm2den *= factfact(2*absm);
    reduce(norm2num,norm2den);
    norm2num *= fact(absm);
    norm2den *= factfact(2*absm);
    if (m != 0) norm2num *= 2;
    reduce(norm2num,norm2den);

    for (unsigned int t=0; t <= (l - absm)/2; t++) {
        for (unsigned int u=0; u<=t; u++) {
            int v2m;
            if (m >= 0) v2m = 0;
            else v2m = 1;
            for (unsigned int v2 = v2m; v2 <= absm; v2+=2) {
                int x = 2*t + absm - 2*u - v2;
                int y = 2*u + v2;
                int z = l - x - y;
                uint64_t bin = binomial(l,t)
                        *binomial(l-t,absm+t)
                        *binomial(t,u)
                        *binomial(absm,v2);
                uint64_t den = powll(4,t);
                int sign;
                if ((t + (v2-v2m)/2)%2) sign = -1;
                else sign = 1;
                reduce(bin,den);
                solidharmcontrib(sign,bin,den,norm2num,norm2den,
                                 r2,x,y,z,coefmat,pureindex);
            }
        }
    }
}

void solidharmonic(int l, SimpleMatrix &coefmat)
{
    solidharm(l,0,0,coefmat);
    for (int m=1; m<=l; m++) {
        solidharm(l, m,0,coefmat);
        solidharm(l,-m,0,coefmat);
    }
    for (int r=2; r<=l; r+=2) {
        solidharm(l-r,0,r/2,coefmat);
        for (int m=1; m<=l-r; m++) {
            solidharm(l-r, m,r/2,coefmat);
            solidharm(l-r,-m,r/2,coefmat);
        }
    }
}

