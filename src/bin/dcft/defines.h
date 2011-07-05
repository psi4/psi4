#ifndef _DEFINES_H
#define	_DEFINES_H

#define GIT_ID "d12f233900069eb274854278e3aa1c733c34c9e6"

#define PSIF_DCFT_DPD 100
#define PSIF_DCFT_DENSITY 101
#define PRINT_ENERGY_COMPONENTS 0
#define ZERO 1.0E-16

#define REFACTORED 1

#define ID(x) _ints->DPD_ID(x)

#ifndef INDEX
#define INDEX(i,j) (i > j ? i * (i + 1) / 2 + j : j + (j + 1) / 2 + i)
#endif

#endif	/* _DEFINES_H */

