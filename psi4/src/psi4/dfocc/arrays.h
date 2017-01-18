/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _dfocc_arrays_h_
#define _dfocc_arrays_h_

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

class Array1d;
class Array2d;
class Array3d;
class Array1i;
class Array2i;
class Array3i;

class Array1d
{

  private:
  double *A1d_;
  int dim1_;
  string name_;      // Name of the array

  public:
  Array1d(int d1);
  Array1d(string name, int d1);
  Array1d();			   //default constructer
  ~Array1d(); 		   	   //destructer

  Array1d* generate(int d1);
  Array1d* generate(string name, int d1);
  void init(string name, int d1);
  void init(int d1);
  void memalloc();
  void zero();
  void print();
  void print(std::string OutFileRMR);
  void release();
  void set(int i, double value);
  void set(double *vec);
  void set(const Array1d  *vec);
  void add(const Array1d* Adum);
  void add(int i, double value);// add value to ith element of the vector
  void subtract(const Array1d* Adum);
  void subtract(int i, double value);
  double get(int i);
  // rms:  rms of A1d_
  double rms();
  // rms:  rms of (A1d_ - Atemp)
  double rms(const Array1d* Atemp);
  // dot: return result of A1d_' * y
  double dot(const Array1d *y);
  // gemv: A1d_ = alpha * A * b + beta, where A is a general matrix
  void gemv(bool transa, double alpha, const Array2d* a, const Array1d* b, double beta);
  // gbmv: This function may NOT working correctly!!!!
  void gbmv(bool transa, double alpha, const Array2d* a, const Array1d* b, double beta);
  // xay: return result of A1d_' * A * y
  double xay(const Array2d *a, const Array1d *y);
  void scale(double a);
  void copy(double *x);
  void copy(const Array1d *x);
  // row_vector: set A1d to nth row of A, dim1_ = A->dim2
  void row_vector(Array2d *A, int n);
  // column_vector: set A1d to nth column of A, dim1_ = A->dim1
  void column_vector(Array2d *A, int n);
  int dim1() const { return dim1_; }
  // dirprd: A1d_[i] = a[i] * b[i]
  void dirprd(Array1d *a, Array1d *b);

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
  Array2d(string name, int d1,int d2);
  Array2d(psi::PSIO* psio, unsigned int fileno, string name, int d1,int d2);
  Array2d(std::shared_ptr<psi::PSIO> psio, unsigned int fileno, string name, int d1,int d2);
  Array2d(psi::PSIO& psio, unsigned int fileno, string name, int d1,int d2);
  Array2d();			   //default constructer
  ~Array2d(); 		   	   //destructer

  Array2d* generate(int d1,int d2);
  Array2d* generate(string name, int d1,int d2);
  void init(string name, int d1,int d2);
  void init(int d1,int d2);
  void memalloc();
  void zero();
  void zero_diagonal();
  void print();
  void print(std::string OutFileRMR);
  void release();
  void set(int i, int j, double value);
  void set(double **A);
  void set(Array2d *A);
  void set(SharedMatrix A);
  double get(int i, int j);
  void add(const Array2d* Adum);
  // A2d = alpha * Adum
  void add(double alpha, const Array2d* Adum);
  void add(int i, int j, double value);
  void subtract(const Array2d* Adum);
  void subtract(int i, int j, double value);
  Array2d* transpose();
  void copy(const Array2d* Adum);
  void copy(double **a);
  // diagonalize: diagonalize via rsp
  void diagonalize(Array2d* eigvectors, Array1d* eigvalues, double cutoff);
  // cdsyev: diagonalize via lapack
  void cdsyev(char jobz, char uplo, Array2d* eigvectors, Array1d* eigvalues);
  // davidson: diagonalize via davidson algorithm
  void davidson(int n_eigval, Array2d* eigvectors, Array1d* eigvalues, double cutoff, int print);
  // cdgesv: solve a linear equation via lapack
  void cdgesv(Array1d* Xvec);
  void cdgesv(double* Xvec);
  void cdgesv(Array1d* Xvec, int errcod);
  void cdgesv(double* Xvec, int errcod);
  // lineq_flin: solve a linear equation via FLIN
  void lineq_flin(Array1d* Xvec, double *det);
  void lineq_flin(double* Xvec, double *det);
  // lineq_pople: solve a linear equation via Pople's algorithm
  void lineq_pople(Array1d* Xvec, int num_vecs, double cutoff);
  void lineq_pople(double* Xvec, int num_vecs, double cutoff);
  // gemm: matrix multiplication C = A * B
  void gemm(bool transa, bool transb, const Array2d* a, const Array2d* b, double alpha, double beta);
  // contract: general contraction C(m,n) = \sum_{k} A(m,k) * B(k,n)
  void contract(bool transa, bool transb, int m, int n, int k, const Array2d* a, const Array2d* b, double alpha, double beta);
  // contract323: C[Q](m,n) = \sum_{k} A[Q](m,k) * B(k,n)
  void contract323(bool transa, bool transb, int m, int n, const Array2d* a, const Array2d* b, double alpha, double beta);
  // contract233: C[Q](m,n) = \sum_{k} A(m,k) * B[Q](k,n)
  void contract233(bool transa, bool transb, int m, int n, const Array2d* a, const Array2d* b, double alpha, double beta);
  // level_shift: A[i][i] = A[i][i] - value
  void level_shift(double value);
  // outer_product: A = x * y'
  void outer_product(const Array1d *x, const Array1d *y);
  void scale(double a);
  // scale_row: scales mth row with a
  void scale_row(int m, double a);
  // scale_column: scales nth column with a
  void scale_column(int n, double a);
  // identity: A = I
  void identity();
  double trace();
  // transform: A = L' * B * L
  void transform(const Array2d* a, const Array2d* transformer);
  // back_transform: A = L * B * L'
  void back_transform(const Array2d* a, const Array2d* transformer);
  // pseudo_transform: A = L * B * L
  void pseudo_transform(const Array2d* a, const Array2d* transformer);
  // triple_gemm: A2d_ = a * b * c
  void triple_gemm(const Array2d* a, const Array2d* b, const Array2d* c);
  // vector_dot: value = Tr(A' * B)
  double vector_dot(Array2d *rhs);
  double vector_dot(double **rhs);
  double **to_block_matrix();
  double *to_lower_triangle();
  void to_shared_matrix(SharedMatrix A);
  // mgs: orthogonalize with a Modified Gram-Schmid algorithm
  void mgs();
  // gs: orthogonalize with a Classical Gram-Schmid algorithm
  void gs();
  // row_vector: return nth row as a vector
  double *row_vector(int n);
  // column_vector: return nth column as a vector
  double *column_vector(int n);
  int dim1() const { return dim1_; }
  int dim2() const { return dim2_; }

  void write(std::shared_ptr<psi::PSIO> psio, unsigned int fileno);
  void write(psi::PSIO* const psio, unsigned int fileno);
  void write(psi::PSIO& psio, unsigned int fileno);
  void read(psi::PSIO* psio, unsigned int fileno);
  void read(std::shared_ptr<psi::PSIO> psio, unsigned int fileno);
  void read(psi::PSIO& psio, unsigned int fileno);
  bool read(PSIO* psio, int itap, const char *label, int dim);
  bool read(std::shared_ptr<psi::PSIO> psio, int itap, const char *label, int dim);
  void save(std::shared_ptr<psi::PSIO> psio, unsigned int fileno);
  void save(psi::PSIO* const psio, unsigned int fileno);
  void save(psi::PSIO& psio, unsigned int fileno);
  void load(std::shared_ptr<psi::PSIO> psio, unsigned int fileno, string name, int d1,int d2);
  void load(psi::PSIO* const psio, unsigned int fileno, string name, int d1,int d2);
  void load(psi::PSIO& psio, unsigned int fileno, string name, int d1,int d2);

  // sort1432: A2d_(ps,rq) = A(pq,rs): d1 = num p, d2 = num q, d3 = num r, d4 = num s
  // col/row_pair_idx are belong to A, while col/row_pairidx2 are belong to A2d_
  void sort1432(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx, Array2i *row_pair_idx2, Array2i *col_pair_idx2);
  // sort2134: A2d_(qp,rs) = A(pq,rs): d1 = num p, d2 = num q, d3 = num r, d4 = num s
  void sort2134(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx, Array2i *row_pair_idx2);
  // sort1243: A2d_(pq,sr) = A(pq,sr): d1 = num p, d2 = num q, d3 = num r, d4 = num s
  void sort1243(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx, Array2i *col_pair_idx2);
  // sort2413: A2d_(qs,pr) = A(pq,sr): d1 = num p, d2 = num q, d3 = num r, d4 = num s
  void sort2413(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx, Array2i *row_pair_idx2, Array2i *col_pair_idx2);
  // sort2143: A2d_(qp,sr) = A(pq,sr): d1 = num p, d2 = num q, d3 = num r, d4 = num s
  void sort2143(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx, Array2i *row_pair_idx2, Array2i *col_pair_idx2);
  // sort4231: A2d_(sq,rp) = A(pq,sr): d1 = num p, d2 = num q, d3 = num r, d4 = num s
  void sort4231(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx, Array2i *row_pair_idx2, Array2i *col_pair_idx2);
  // sort3142: A2d_(rp,sq) = A(pq,sr): d1 = num p, d2 = num q, d3 = num r, d4 = num s
  void sort3142(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx, Array2i *row_pair_idx2, Array2i *col_pair_idx2);

  // sortp1432: A2d_(qp,rs) += A(pq,rs): d1 = num p, d2 = num q, d3 = num r, d4 = num s
  void sortp1432(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx, Array2i *row_pair_idx2, Array2i *col_pair_idx2);
  // sortp2134: A2d_(qp,rs) += A(pq,rs): d1 = num p, d2 = num q, d3 = num r, d4 = num s
  void sortp2134(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx, Array2i *row_pair_idx2);
  // sortp1243: A2d_(pq,sr) += A(pq,sr): d1 = num p, d2 = num q, d3 = num r, d4 = num s
  void sortp1243(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx, Array2i *col_pair_idx2);
  // sortp2413: A2d_(qs,pr) += A(pq,sr): d1 = num p, d2 = num q, d3 = num r, d4 = num s
  void sortp2413(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx, Array2i *row_pair_idx2, Array2i *col_pair_idx2);
  // sortp2143: A2d_(qp,sr) += A(pq,sr): d1 = num p, d2 = num q, d3 = num r, d4 = num s
  void sortp2143(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx, Array2i *row_pair_idx2, Array2i *col_pair_idx2);
  // sortp4231: A2d_(sq,rp) += A(pq,sr): d1 = num p, d2 = num q, d3 = num r, d4 = num s
  void sortp4231(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx, Array2i *row_pair_idx2, Array2i *col_pair_idx2);
  // sortp3142: A2d_(rp,sq) += A(pq,sr): d1 = num p, d2 = num q, d3 = num r, d4 = num s
  void sortp3142(int d1, int d2, int d3, int d4, Array2d *A, Array2i *row_pair_idx, Array2i *col_pair_idx, Array2i *row_pair_idx2, Array2i *col_pair_idx2);


  friend class Array1d;
  friend class Array3d;
  friend class Array1i;
  friend class Array2i;
};

class Array3d
{

  private:
  double ***A3d_;
  int dim1_,dim2_,dim3_;
  string name_;      // Name of the array

  public:
  Array3d(int d1,int d2, int d3);
  Array3d(string name, int d1,int d2, int d3);
  Array3d();			   //default constructer
  ~Array3d(); 		   	   //destructer

  Array3d* generate(int d1,int d2, int d3);
  Array3d* generate(string name, int d1,int d2, int d3);
  void init(string name, int d1,int d2, int d3);
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
  Array1i(string name, int d1);
  Array1i();			   //default constructer
  ~Array1i(); 		   	   //destructer

  Array1i* generate(int d1);
  Array1i* generate(string name, int d1);
  void init(string name, int d1);
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

class Array2i
{

  private:
  int **A2i_;
  int dim1_,dim2_;
  string name_;      // Name of the array

  public:
  Array2i(int d1,int d2);
  Array2i(string name, int d1,int d2);
  Array2i();			   //default constructer
  ~Array2i(); 		   	   //destructer

  Array2i* generate(int d1,int d2);
  Array2i* generate(string name, int d1,int d2);
  void init(string name, int d1,int d2);
  void init(int d1,int d2);
  void memalloc();
  void zero();
  void zero_diagonal();
  void print();
  void print(std::string OutFileRMR);
  void release();
  void set(int i, int j, int value);
  void set(int **A);
  double get(int i, int j);
  void add(const Array2i* Adum);
  void add(int i, int j, int value);
  void subtract(const Array2i* Adum);
  void subtract(int i, int j, int value);
  Array2i* transpose();
  void copy(const Array2i* Adum);
  void copy(int **a);
  void identity();
  int trace();
  int **to_int_matrix();
  int dim1() const { return dim1_; }
  int dim2() const { return dim2_; }

  friend class Array1i;
  friend class Array3i;
  friend class Array1d;
  friend class Array2d;
};


class Array3i
{

  private:
  int ***A3i_;
  int dim1_,dim2_,dim3_;
  string name_;      // Name of the array

  public:
  Array3i(int d1,int d2, int d3);
  Array3i(string name, int d1,int d2, int d3);
  Array3i();			   //default constructer
  ~Array3i(); 		           //destructer

  Array3i* generate(int d1,int d2, int d3);
  Array3i* generate(string name, int d1,int d2, int d3);
  void init(string name, int d1,int d2, int d3);
  void init(int d1,int d2, int d3);
  void memalloc();
  void zero();
  void print();
  void release();
  void set(int h, int i, int j, int value);
  int get(int h, int i, int j);
};
}} // End Namespaces
#endif // _dfocc_arrays_h_
