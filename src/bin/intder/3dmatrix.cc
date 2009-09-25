/*! 
** \file
** \ingroup INTDER
** \brief Enter brief description of file here 
*/

/******************************************************************

3dmatrix.cc

A small batch of functions in order to manipulate 3d/4d/5d matrices

*******************************************************************/

#include <cstdlib>
#include "3dmatrix.h"

#include <libciomr/libciomr.h>

using namespace psi::intder;

C3DMatrix::C3DMatrix()
{
  array = NULL;
}

C3DMatrix::C3DMatrix(int x, int y, int z)
{
  Create(x, y, z);
}

C3DMatrix::~C3DMatrix()
{
  if (array)
    free(array);
  array = NULL;
}

void C3DMatrix::Create(int x, int y, int z)
{
  Size(x,y,z);
  array = (double*) malloc(size * sizeof(double));
  nx = x;
  ny = y;
  nz = z;
}

C4DMatrix::C4DMatrix()
{
  array = NULL;
}
                                                                                                
C4DMatrix::C4DMatrix(int w, int x, int y, int z)
{
  Create(w, x, y, z);
}
                                                                                                
C4DMatrix::~C4DMatrix()
{
  if (array)
    free(array);
  array = NULL;
}
                                                                                                
void C4DMatrix::Create(int w, int x, int y, int z)
{
  Size(w, x, y, z);
  array = (double*) malloc(size * sizeof(double));
	nw = w;
  nx = x;
  ny = y;
  nz = z;
}

C5DMatrix::C5DMatrix()
{
  array = NULL;
}
                                                                                                
C5DMatrix::C5DMatrix(int v, int w, int x, int y, int z)
{
  Create(v, w, x, y, z);
}
                                                                                                
C5DMatrix::~C5DMatrix()
{
  if (array)
    free(array);
  array = NULL;
}
                                                                                                
void C5DMatrix::Create(int v, int w, int x, int y, int z)
{
  Size(v, w, x, y, z);
  array = (double*) malloc(size * sizeof(double));
  nv = v;
  nw = w;
  nx = x;
  ny = y;
  nz = z;
}

