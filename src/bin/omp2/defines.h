#ifndef DEFINES_
#define DEFINES_

#define PSIF_OMP2_DPD 100
#define PSIF_OMP2_DENSITY 101

#define ID(x) ints->DPD_ID(x)

#define index2(i,j) ((i>j) ? ((i*(i+1)/2)+j) : ((j*(j+1)/2)+i))
#define index4(i,j,k,l) index2(index2(i,j),index2(k,l))

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

#define DIIS_MIN_DET 1.0E-16

#endif /* DEFINES_ */

