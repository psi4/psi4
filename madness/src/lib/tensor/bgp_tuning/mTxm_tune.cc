#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>

using namespace std;

typedef complex<double> double_complex;

void mTxm_tune(long dimi, long dimj, long dimk, double* c, const double* a, const double* b) {
    int i, j, k;
    for (k=0; k<dimk; ++k) {
        for (j=0; j<dimj; ++j) {
            for (i=0; i<dimi; ++i) {
                c[i*dimj+j] += a[k*dimi+i]*b[k*dimj+j];
            }
        }
    }
}

void mTxm_tune(long dimi, long dimj, long dimk, double_complex* c, const double_complex* a, const double_complex* b) {
    int i, j, k;
    for (k=0; k<dimk; ++k) {
        for (j=0; j<dimj; ++j) {
            for (i=0; i<dimi; ++i) {
                c[i*dimj+j] += a[k*dimi+i]*b[k*dimj+j];
            }
        }
    }
}
