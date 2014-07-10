/** @file global.h
 *
 * This is a private header file which defines all Fortran functions.
 * The names of these functions should mirror those found in global.h.
 */

#ifndef GLOBAL_H
#define GLOBAL_H

#include "typesf2c.h"

extern DoubleComplex   *DCPL_MB;
extern SingleComplex   *SCPL_MB;
extern DoublePrecision *DBL_MB;
extern float           *FLT_MB;
extern Integer         *INT_MB;

#endif /* GLOBAL_H */
