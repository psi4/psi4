/*! \file
**  \ingroup DETCI
**  \brief Contains C code for vector operations
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia
** 
*/

#include <cstdio>

namespace psi { namespace detci {

/*
** xey
** 
** Perform the operation X[] = Y[] for vectors 'x' and 'y'
** of length 'size'
**
*/
void xey(double *x, double *y, int size) {
   int i;

   for (i=0; i<size; i++) {
      x[i] = y[i];
      }
}

/*
** xeay
**
** Perform the operation X[] = a * Y[] for vectors 'x' and 'y' 
**   (of length 'size') and constant 'a'.
**
** David Sherrill, November 1995
*/
void xeay(double *x, double a, double *y, int size)
{
   int i;

   for (i=0; i<size; i++) {
      x[i] = a * y[i];
      }
}



/*
** xpeay
**
** Perform the operation X[] += A * Y[] for vectors 'x' and 'y' 
**   (of length 'size') and constant 'a'.
**
** David Sherrill, November 1995
*/
void xpeay(double *x, double a, double *y, int size)
{
   int i;

   for (i=0; i<size; i++) {
      x[i] += a * y[i];
      }
}

/*
** xeax
**
** Perform the operation X[] = A * X[] for vector 'X' and constant 'A'
**
** David Sherrill, February 1996
**
*/
void xeax(double *x, double a, int size)
{
   int i;

   for (i=0; i<size; i++) {
      x[i] *= a;
      }
}


/*
** xeaxmy
**
** Perform X[] = A * X[] - Y[]
**
** David Sherrill, February 1996
**
*/
void xeaxmy(double *x, double *y, double a, int size)
{
   int i;

   for (i=0; i<size; i++) {
      x[i] = x[i] * a - y[i];
      }
}

/*
** xeaxpby
**
** Perform X[] = A * X[] - B * Y[]
** 
** David Sherrill, March 1996
**
*/
void xeaxpby(double *x, double *y, double a, double b, int size)
{
   int i;

   for (i=0; i<size; i++) {
      x[i] = a * x[i] + b * y[i];
      }
}
/*
** xexy
**
** Perform X[] = X[] * Y[]
**
** Matt Leininger, September 1998
**
*/
void xexy(double *x, double *y, int size)
{
  int i;
  
  for (i=0; i<size; i++) {
     x[i] *= y[i];
     }
}

/*
** xexmy
**
** Perform X[] = X[] - Y[]
**
** Matt Leininger and Nick Petraco, February 1999
**
*/
void xexmy(double *x, double *y, int size)
{
   int i;

   for (i=0; i<size; i++) {
      x[i] -= y[i];
      }
}
 
/*
** xpey
**
** Perform X[] = X[] + Y[]
**
** Matt Leininger February 1999
**
*/
void xpey(double *x, double *y, int size)
{
   int i;

   for (i=0; i<size; i++) {
      x[i] += y[i];
      }
}

}} // namespace psi::detci

