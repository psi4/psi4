/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/

#ifndef _psi_bin_intder_3dmatrix_h_
#define _psi_bin_intder_3dmatrix_h_

#include <libciomr/libciomr.h>

namespace psi { namespace intder {

class C3DMatrix
{
  int nx, ny, nz;
  int size;

public:
  double *array;

  C3DMatrix();
  C3DMatrix(int, int, int);
  ~C3DMatrix();

  void Create(int, int, int);

  inline double Get(int x, int y, int z)
    { return array[Position(x, y, z)]; }
  inline int Size(int x, int y, int z) 
    { return size = Position(x, y, z); }
  inline void Set(int x, int y, int z, double value)
    { array[Position(x, y, z)] = value; }
	inline void PlusEq(int x, int y, int z, double value)
		{ array[Position(x, y, z)] += value; }
  inline int Position(int x, int y, int z)
    { return (nz*ny*x + nz*y + z); }
  inline void Init(int x, int y, int z)
  { array = init_array(size); }
};

class C4DMatrix
{
  int nw, nx, ny, nz;
  int size;

public:
  double *array;
                                                                                            
  C4DMatrix();
  C4DMatrix(int, int, int, int);
  ~C4DMatrix();
                                                                                                
  void Create(int, int, int, int);

  inline int Size(int w,int x, int y, int z) 
    { return size = Position(w, x, y, z); }
  inline double Get(int w, int x, int y, int z)
    { return array[Position(w, x, y, z)]; }
  inline void Set(int w, int x, int y, int z, double value)
    { array[Position(w, x, y, z)] = value; }
  inline int Position(int w, int x, int y, int z)
    { return (nz*ny*nx*w + nz*ny*x + nz*y + z); }
  inline void Init(int w,int x, int y, int z)
  { array = init_array(size); }
};

class C5DMatrix
{
  int nv, nw, nx, ny, nz;
  int size;

public:
  double *array;
                                                                                             
  C5DMatrix();
  C5DMatrix(int, int, int, int, int);
  ~C5DMatrix();
                                                                                                
  void Create(int, int, int, int, int);

  inline int Size(int v, int w,int x, int y, int z) 
    { return size = Position(v, w, x, y, z); }
  inline double Get(int v, int w, int x, int y, int z)
    { return array[Position(v, w, x, y, z)]; }
  inline void Set(int v, int w, int x, int y, int z, double value)
    { array[Position(v, w, x, y, z)] = value; }
  inline int Position(int v, int w, int x, int y, int z)
    { return (nz*ny*nx*nw*v + nz*ny*nx*w + nz*ny*x + nz*y + z); }
  inline void Init(int v, int w, int x, int y, int z)
  { array = init_array(size); }
};

}} // namespace psi::intder

#endif // heaer guard
