#ifndef _psi_src_bin_omp2_arrays_h_
#define _psi_src_bin_omp2_arrays_h_

/** Standard library includes */
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string> 
#include <iomanip>
#include <vector> 

/** Required PSI4 includes */
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>


/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>


using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace omp2wave{
  
class Array1d;
class Array2d;
class Array3d;
class Array1i;
class Array3i;


  
class Array1d
{

  private:
  double *A1d_;
  int dim1_;
  string name_;      // Name of the array
    
  public:
  Array1d(int d1);  
  Array1d(int d1, string name);
  Array1d();			   //default constructer
  ~Array1d(); 		   	   //destructer
  
  Array1d* generate(int d1);
  Array1d* generate(int d1, string name);
  void init(int d1, string name);
  void init(int d1);
  void memalloc();
  void zero();
  void print();
  void print(FILE *out);
  void release();
  void set(int i, double value);
  void set(double *vec);
  void set(const Array1d  *vec);
  void add(const Array1d* Adum);
  void add(int i, double value);// add value to ith element of the vector
  void subtract(const Array1d* Adum);
  void subtract(int i, double value);
  double get(int i); 
  double rms(); 
  double rms(const Array1d* Atemp);//  rms of (A1d_ - Atemp)
  double dot(const Array1d *y); // return result of A1d_' * y
  void gemv(bool transa, double alpha, const Array2d* a, const Array1d* b, double beta);
  void gbmv(bool transa, double alpha, const Array2d* a, const Array1d* b, double beta);//This function may NOT working correctly!!!!
  double xay(const Array2d *a, const Array1d *y); // return result of A1d_' * A * y
  void scale(double a);
  void copy(double *x);
  void copy(const Array1d *x);
  void row_vector(Array2d *A, int n); // set A1d to nth row of A, dim1_ = A->dim2
  void column_vector(Array2d *A, int n); // set A1d to nth column of A, dim1_ = A->dim1
  int dim1() const { return dim1_; }
  
  friend class Array2d;
  friend class Array3d;
};  


class Array2d
{

  private:
  double **A2d_;
  int dim1_,dim2_;
  string name_;      // Name of the array
    
  public:
  Array2d(int d1,int d2);  
  Array2d(int d1,int d2, string name);
  Array2d();			   //default constructer
  ~Array2d(); 		   	   //destructer
  
  Array2d* generate(int d1,int d2);
  Array2d* generate(int d1,int d2, string name);
  void init(int d1,int d2, string name);
  void init(int d1,int d2);
  void memalloc();
  void zero();
  void zero_diagonal();
  void print();
  void print(FILE *out);
  void release();
  void set(int i, int j, double value);
  double get(int i, int j);
  void add(const Array2d* Adum);
  void add(int i, int j, double value);
  void subtract(const Array2d* Adum);
  void subtract(int i, int j, double value);
  Array2d* transpose();
  void copy(const Array2d* Adum);
  void copy(double **a);
  void diagonalize(Array2d* eigvectors, Array1d* eigvalues, double cutoff);
  void cdsyev(char jobz, char uplo, Array2d* eigvectors, Array1d* eigvalues); // diagonalize via acml
  void davidson(int n_eigval, Array2d* eigvectors, Array1d* eigvalues, double cutoff, int print); 
  void cdgesv(Array1d* Xvec); // solve lineq via acml
  void cdgesv(double* Xvec); // solve lineq via acml
  void lineq_flin(Array1d* Xvec, double *det); // solve lineq via flin
  void lineq_flin(double* Xvec, double *det); // solve lineq via flin
  void lineq_pople(Array1d* Xvec, int num_vecs, double cutoff); // solve lineq via pople  
  void lineq_pople(double* Xvec, int num_vecs, double cutoff); // solve lineq via pople  
  void gemm(bool transa, bool transb, double alpha, const Array2d* a, const Array2d* b, double beta);
  void level_shift(double value);
  void outer_product(const Array1d *x, const Array1d *y); // A = x * y' 
  void scale(double a);
  void scale_row(int m, double a);// scales mth row with a
  void scale_column(int n, double a);// scales nth column with a
  void identity();
  double trace();
  void transform(const Array2d* a, const Array2d* transformer);// A = L' * B * L
  void back_transform(const Array2d* a, const Array2d* transformer);// A = L * B * L'
  void pseudo_transform(const Array2d* a, const Array2d* transformer);// A = L * B * L
  void triple_gemm(const Array2d* a, const Array2d* b, const Array2d* c);// A2d_ = a * b * c
  double vector_dot(Array2d *rhs);
  double vector_dot(double **rhs);
  /*
  void write(psi::PSIO* psio, unsigned int fileno);
  void write(shared_ptr<psi::PSIO> psio, unsigned int fileno);
  void write(psi::PSIO& psio, unsigned int fileno);
  void read(psi::PSIO* psio, unsigned int fileno);
  void read(psi::PSIO& psio, unsigned int fileno);
  bool read(PSIO* psio, int itap, const char *label, int dim);
  bool read(shared_ptr<psi::PSIO> psio, int itap, const char *label, int dim);
  */
  double **to_block_matrix(); 
  double *to_lower_triangle();
  void mgs();
  void gs();
  int dim1() const { return dim1_; }
  int dim2() const { return dim2_; }
  double *row_vector(int n);// return nth row as a vector
  double *column_vector(int n);// return nth column as a vector
  
  
  friend class Array1d;
  friend class Array3d;
};

class Array3d
{

  private:
  double ***A3d_;
  int dim1_,dim2_,dim3_;
  string name_;      // Name of the array
    
  public:
  Array3d(int d1,int d2, int d3);  
  Array3d(int d1,int d2, int d3, string name);
  Array3d();			   //default constructer
  ~Array3d(); 		   	   //destructer
  
  Array3d* generate(int d1,int d2, int d3);
  Array3d* generate(int d1,int d2, int d3, string name);
  void init(int d1,int d2, int d3, string name);
  void init(int d1,int d2, int d3);
  void memalloc();
  void zero();
  void print();
  void release();
  void set(int h, int i, int j, double value);
  double get(int h, int i, int j); 
  
  friend class Array1d;
  friend class Array2d;
};

class Array1i
{

  private:
  int *A1i_;
  int dim1_;
  string name_;      // Name of the array
    
  public:
  Array1i(int d1);  
  Array1i(int d1, string name);
  Array1i();			   //default constructer
  ~Array1i(); 		   	   //destructer
  
  Array1i* generate(int d1);
  Array1i* generate(int d1, string name);
  void init(int d1, string name);
  void init(int d1);
  void memalloc();
  void zero();
  void print();
  void release();
  void set(int i, int value);
  int get(int i); 
  void add(const Array1i* Adum);
  void add(int i, int value);
  void subtract(const Array1i* Adum);
  void subtract(int i, int value);

}; 

class Array3i
{

  private:
  int ***A3i_;
  int dim1_,dim2_,dim3_;
  string name_;      // Name of the array
    
  public:
  Array3i(int d1,int d2, int d3);  
  Array3i(int d1,int d2, int d3, string name);
  Array3i();			   //default constructer
  ~Array3i(); 		           //destructer
  
  Array3i* generate(int d1,int d2, int d3);
  Array3i* generate(int d1,int d2, int d3, string name);
  void init(int d1,int d2, int d3, string name);
  void init(int d1,int d2, int d3);
  void memalloc();
  void zero();
  void print();
  void release();
  void set(int h, int i, int j, int value);
  int get(int h, int i, int j); 
};
}} // End Namespaces
#endif // _psi_src_bin_omp2_arrays_h_

