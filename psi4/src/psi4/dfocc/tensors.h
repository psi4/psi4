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

#ifndef _dfocc_tensors_h_
#define _dfocc_tensors_h_

#include "psi4/libpsio/psio.hpp"

#define index2(i,j) ((i>j) ? ((i*(i+1)/2)+j) : ((j*(j+1)/2)+i))
#define index4(i,j,k,l) index2(index2(i,j),index2(k,l))
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define idx_asym(i,j) ((i>j) ? ((i*(i-1)/2)+j) : ((j*(j-1)/2)+i))


using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

class Tensor1d;
class Tensor2d;
class Tensor3d;
class Tensor1i;
class Tensor2i;
class Tensor3i;

typedef std::shared_ptr<Tensor1d> SharedTensor1d;
typedef std::shared_ptr<Tensor2d> SharedTensor2d;
typedef std::shared_ptr<Tensor3d> SharedTensor3d;
typedef std::shared_ptr<Tensor1i> SharedTensor1i;
typedef std::shared_ptr<Tensor2i> SharedTensor2i;
typedef std::shared_ptr<Tensor3i> SharedTensor3i;

class Tensor1d
{

  private:
  double *A1d_;
  int dim1_;
  string name_;      // Name of the array

  public:
  Tensor1d(int d1);
  Tensor1d(string name, int d1);
  Tensor1d();			   //default constructer
  ~Tensor1d(); 		   	   //destructer

  void init(string name, int d1);
  void init(int d1);
  void memalloc();
  void zero();
  void print();
  void print(std::string OutFileRMR);
  void release();
  void set(int i, double value);
  void set(double *vec);
  void set(const SharedTensor1d  &vec);
  void add(const SharedTensor1d& Adum);
  void add(int i, double value);// add value to ith element of the vector
  void subtract(const SharedTensor1d& Adum);
  void subtract(int i, double value);
  double get(int i);
  void to_shared_vector(SharedVector A);
  // rms:  rms of A1d_
  double rms();
  // rms:  rms of (A1d_ - Atemp)
  double rms(const SharedTensor1d& Atemp);
  // dot: return result of A1d_' * y
  double dot(const SharedTensor1d &y);
  // gemv: A1d_ = alpha * A * b + beta, where A is a general matrix
  // gemv: C(m) = \sum_{n} A(m,n) b(n)
  void gemv(bool transa, const SharedTensor2d& a, const SharedTensor1d& b, double alpha, double beta);
  // gemv: C(m) = \sum_{n} A(m,n) b(n)
  void gemv(bool transa, int m, int n, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta);
  void gemv(bool transa, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta);
  void gemv(bool transa, int m, int n, const SharedTensor2d& a, const SharedTensor2d& b, int start_a, int start_b, double alpha, double beta);
  void gemv(bool transa, int m, int n, const SharedTensor2d& a, const SharedTensor2d& b, int start_a, int start_b, int start_c, double alpha, double beta);
  // gbmv: This function may NOT working correctly!!!!
  void gbmv(bool transa, const SharedTensor2d& a, const SharedTensor1d& b, double alpha, double beta);
  // xay: return result of A1d_' * A * y
  double xay(const SharedTensor2d &a, const SharedTensor1d &y);
  // axpy: Y <-- a * X + Y
  void axpy(const SharedTensor1d &a, double alpha);
  void scale(double a);
  void copy(double *x);
  void copy(const SharedTensor1d &x);
  // row_vector: set A1d to nth row of A, dim1_ = A->dim2
  void row_vector(SharedTensor2d &A, int n);
  // column_vector: set A1d to nth column of A, dim1_ = A->dim1
  void column_vector(SharedTensor2d &A, int n);
  int dim1() const { return dim1_; }
  // dirprd: A1d_[i] = a[i] * b[i]
  void dirprd(SharedTensor1d &a, SharedTensor1d &b);
  // symm_packed A(p>=q) = A(p,q) * (2 - \delta_{pq})
  void symm_packed(const SharedTensor2d &A);
  // ltm: lower triangular matrix: A(p>=q) = A(p,q)
  void ltm(const SharedTensor2d &A);

  friend class Tensor2d;
  friend class Tensor3d;
};


class Tensor2d
{

  private:
  double **A2d_;
  int dim1_, dim2_, d1_, d2_, d3_, d4_;
  int **row_idx_, **col_idx_;
  int *row2d1_, *row2d2_, *col2d1_, *col2d2_;
  string name_;      // Name of the array

  public:
  Tensor2d(int d1,int d2);
  Tensor2d(string name, int d1,int d2);
  Tensor2d(psi::PSIO* psio, unsigned int fileno, string name, int d1,int d2);
  Tensor2d(std::shared_ptr<psi::PSIO> psio, unsigned int fileno, string name, int d1,int d2);
  Tensor2d(psi::PSIO& psio, unsigned int fileno, string name, int d1,int d2);
  Tensor2d(string name, int d1, int d2, int d3, int d4);
  Tensor2d(string name, int d1, int d2, int d3);
  Tensor2d();			   //default constructer
  ~Tensor2d(); 		   	   //destructer

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
  void set(double *A);
  void set(SharedTensor2d &A);
  void set(SharedMatrix A);
  void set2(SharedMatrix A);
  void set(SharedTensor1d &A);
  // A2d_[n][ij] = A(i,j)
  void set_row(const SharedTensor2d &A, int n);
  // A2d_[ij][n] = A(i,j)
  void set_column(const SharedTensor2d &A, int n);
  double get(int i, int j);
  // A2d_[ij] = A(n, ij)
  void get_row(const SharedTensor2d &A, int n);
  // A2d_[ij] = A(ij, n)
  void get_column(const SharedTensor2d &A, int n);
  // A2d = alpha * Adum
  void add(const SharedTensor2d &a);
  void add(double **a);
  void add(double alpha, const SharedTensor2d &a);
  void add(int i, int j, double value);
  // A2d_[n][ij] += A(i,j)
  void add2row(const SharedTensor2d &A, int n);
  // A2d_[ij][n] += A(i,j)
  void add2col(const SharedTensor2d &A, int n);
  void subtract(const SharedTensor2d &a);
  void subtract(int i, int j, double value);
  // axpy: Y <-- a * X + Y
  void axpy(double **a, double alpha);
  void axpy(const SharedTensor2d &a, double alpha);
  void axpy(ULI length, int inc_a, const SharedTensor2d &a, int inc_2d, double alpha);
  void axpy(ULI length, int start_a, int inc_a, const SharedTensor2d &A, int start_2d, int inc_2d, double alpha);
  double **transpose2();
  SharedTensor2d transpose();
  void trans(const SharedTensor2d &A);
  void trans(double **A);
  void copy(double **a);
  void copy(const SharedTensor2d &Adum);
  void copy(ULI length, const SharedTensor2d &A, int inc_a, int inc_2d);
  void copy(const SharedTensor2d &A, int start);
  // partial copy
  void pcopy(const SharedTensor2d &A, int dim_copy, int dim_skip);
  void pcopy(const SharedTensor2d &A, int dim_copy, int dim_skip, int start);
  double get_max_element();
  // diagonalize: diagonalize via rsp
  void diagonalize(const SharedTensor2d &eigvectors, const SharedTensor1d &eigvalues, double cutoff);
  // cdsyev: diagonalize via lapack
  void cdsyev(char jobz, char uplo, const SharedTensor2d& eigvectors, const SharedTensor1d& eigvalues);
  // davidson: diagonalize via davidson algorithm
  void davidson(int n_eigval, const SharedTensor2d& eigvectors, const SharedTensor1d& eigvalues, double cutoff, int print);
  // cdgesv: solve a linear equation via lapack
  void cdgesv(const SharedTensor1d& Xvec);
  void cdgesv(double* Xvec);
  void cdgesv(const SharedTensor1d& Xvec, int errcod);
  void cdgesv(double* Xvec, int errcod);
  // lineq_flin: solve a linear equation via FLIN
  void lineq_flin(const SharedTensor1d& Xvec, double *det);
  void lineq_flin(double* Xvec, double *det);
  // pople: solve a linear equation via Pople's algorithm
  void lineq_pople(const SharedTensor1d& Xvec, int num_vecs, double cutoff);
  void lineq_pople(double* Xvec, int num_vecs, double cutoff);

  // gemm: matrix multiplication C = A * B
  void gemm(bool transa, bool transb, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta);
  // contract: general contraction C(m,n) = \sum_{k} A(m,k) * B(k,n)
  void contract(bool transa, bool transb, int m, int n, int k, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta);
  void contract(bool transa, bool transb, int m, int n, int k, const SharedTensor2d& a, const SharedTensor2d& b, int start_a, int start_b, double alpha, double beta);
  void contract(bool transa, bool transb, int m, int n, int k, const SharedTensor2d& a, const SharedTensor2d& b,
                int start_a, int start_b, int start_c, double alpha, double beta);
  // contract323: C[Q](m,n) = \sum_{k} A[Q](m,k) * B(k,n). Note: contract332 should be called with beta=1.0
  void contract323(bool transa, bool transb, int m, int n, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta);
  // contract233: C[Q](m,n) = \sum_{k} A(m,k) * B[Q](k,n)
  void contract233(bool transa, bool transb, int m, int n, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta);
  // contract332: C(m,n) = \sum_{k} A[Q](m,k) * B[Q](k,n)
  void contract332(bool transa, bool transb, int k, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta);
  // contract424: Z(pq,rs) = \sum_{o} X(pq,ro) * Y(o,s): where target_x = 4, target_y = 1, and Z = A2d_
  void contract424(int target_x, int target_y, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta);
  // contract442: C(p,q) \sum_{rst} A(pr,st) B(qr,st) , where row/col pair indices are related to B and sorted B.
  void contract442(int target_a, int target_b, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta);
  // gemv: C(mn) = \sum_{k} A(mn,k) * b(k)
  void gemv(bool transa, const SharedTensor2d& a, const SharedTensor1d& b, double alpha, double beta);
  // gemv: C(mn) = \sum_{kl} A(mn,kl) * b(kl)
  void gemv(bool transa, const SharedTensor2d& a, const SharedTensor2d& b, double alpha, double beta);

  // level_shift: A[i][i] = A[i][i] - value
  void level_shift(double value);
  // outer_product: A = x * y'
  void outer_product(const SharedTensor1d &x, const SharedTensor1d &y);
  void scale(double a);
  // scale_row: scales mth row with a
  void scale_row(int m, double a);
  // scale_column: scales nth column with a
  void scale_column(int n, double a);
  // identity: A = I
  void identity();
  double trace();
  double norm();
  // transform: A = L' * B * L
  void transform(const SharedTensor2d& a, const SharedTensor2d& transformer);
  // back_transform: A = L * B * L'
  void back_transform(const SharedTensor2d& a, const SharedTensor2d& transformer);
  // back_transform: A = alpha* L * B * L' + (beta*A)
  void back_transform(const SharedTensor2d& a, const SharedTensor2d& transformer, double alpha, double beta);
  // pseudo_transform: A = L * B * L
  void pseudo_transform(const SharedTensor2d& a, const SharedTensor2d& transformer);
  // triple_gemm: A2d_ = a * b * c
  void triple_gemm(const SharedTensor2d& a, const SharedTensor2d& b, const SharedTensor2d& c);
  // vector_dot: value = Tr(A' * B)
  double vector_dot(const SharedTensor2d &rhs);
  double vector_dot(double **rhs);
  double **to_block_matrix();
  double *to_lower_triangle();
  void to_shared_matrix(SharedMatrix A);
  void to_matrix(SharedMatrix A);
  void to_pointer(double *A);
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
  void write(std::shared_ptr<psi::PSIO> psio, unsigned int fileno, psio_address start, psio_address *end);
  void write(psi::PSIO* const psio, unsigned int fileno);
  void write(psi::PSIO* psio, unsigned int fileno, psio_address start, psio_address *end);
  void write(psi::PSIO& psio, unsigned int fileno);
  void write(psi::PSIO& psio, unsigned int fileno, psio_address start, psio_address *end);
  void write(std::shared_ptr<psi::PSIO> psio, const string& filename, unsigned int fileno);
  void write(std::shared_ptr<psi::PSIO> psio, unsigned int fileno, bool three_index, bool symm);
  void write(std::shared_ptr<psi::PSIO> psio, const string& filename, unsigned int fileno, bool three_index, bool symm);
  void write_symm(std::shared_ptr<psi::PSIO> psio, unsigned int fileno);
  void write_anti_symm(std::shared_ptr<psi::PSIO> psio, unsigned int fileno);

  void read(psi::PSIO* psio, unsigned int fileno);
  void read(psi::PSIO* psio, unsigned int fileno, psio_address start, psio_address *end);
  void read(std::shared_ptr<psi::PSIO> psio, unsigned int fileno);
  void read(std::shared_ptr<psi::PSIO> psio, unsigned int fileno, psio_address start, psio_address *end);
  void read(psi::PSIO& psio, unsigned int fileno);
  void read(psi::PSIO& psio, unsigned int fileno, psio_address start, psio_address *end);
  void read(std::shared_ptr<psi::PSIO> psio, unsigned int fileno, bool three_index, bool symm);
  void read_symm(std::shared_ptr<psi::PSIO> psio, unsigned int fileno);
  void read_anti_symm(std::shared_ptr<psi::PSIO> psio, unsigned int fileno);

  bool read(PSIO* psio, int itap, const char *label, int dim);
  bool read(std::shared_ptr<psi::PSIO> psio, int itap, const char *label, int dim);
  void save(std::shared_ptr<psi::PSIO> psio, unsigned int fileno);
  void save(psi::PSIO* const psio, unsigned int fileno);
  void save(psi::PSIO& psio, unsigned int fileno);
  void load(std::shared_ptr<psi::PSIO> psio, unsigned int fileno, string name, int d1,int d2);
  void load(psi::PSIO* const psio, unsigned int fileno, string name, int d1,int d2);
  void load(psi::PSIO& psio, unsigned int fileno, string name, int d1,int d2);

  void mywrite(const string& filename);
  void mywrite(int fileno);
  void mywrite(int fileno, bool append);

  void myread(const string& filename);
  void myread(int fileno);
  void myread(int fileno, bool append);
  void myread(int fileno, ULI start);

  // sort (for example 1432 sort): A2d_(ps,rq) = A(pq,rs)
  // A2d_ = alpha*A + beta*A2d_
  void sort(int sort_type, const SharedTensor2d &A, double alpha, double beta);
  // A2d_[p][qr] = sort(A[p][qr])
  void sort3a(int sort_type, int d1, int d2, int d3, const SharedTensor2d &A, double alpha, double beta);
  // A2d_[pq][r] = sort(A[pq][r])
  void sort3b(int sort_type, int d1, int d2, int d3, const SharedTensor2d &A, double alpha, double beta);
  // apply_denom: T(ij,ab) /= D(ij,ab)
  void apply_denom(int frzc, int occ, const SharedTensor2d &fock);
  // apply_denom_os: T(Ij,Ab) /= D(Ij,Ab)
  void apply_denom_os(int frzc, int occA, int occB, const SharedTensor2d &fockA, const SharedTensor2d &fockB);
  // apply_denom_chem: T(ia,jb) /= D(ij,ab)
  void apply_denom_chem(int frzc, int occ, const SharedTensor2d &fock);
  void reg_denom(int frzc, int occ, const SharedTensor2d &fock, double reg);
  void reg_denom_os(int frzc, int occA, int occB, const SharedTensor2d &fockA, const SharedTensor2d &fockB, double reg);
  void reg_denom_chem(int frzc, int occ, const SharedTensor2d &fock, double reg);
  // dirprd: A2d_[i][j] = a[i][j] * b[i][j]
  void dirprd(const SharedTensor2d &a, const SharedTensor2d &b);
  // dirprd123: A2d_[Q][ij] = a[Q] * b[i][j]
  void dirprd123(const SharedTensor1d &a, const SharedTensor2d &b, double alpha, double beta);
  void dirprd123(bool transb, const SharedTensor1d &a, const SharedTensor2d &b, double alpha, double beta);
  // dirprd112: A2d_[i][j] = a[i] * b[j]
  void dirprd112(const SharedTensor1d &a, const SharedTensor1d &b);
  // dirprd112: A2d_[i][j] = alpha *a[i] * b[j] + beta * A2d_[i][j]
  void dirprd112(const SharedTensor1d &a, const SharedTensor1d &b, double alpha, double beta);
  // dirprd224: A2d_[ij][kl] = a[i][j] * b[k][l]
  void dirprd224(const SharedTensor2d &a, const SharedTensor2d &b);
  void dirprd224(const SharedTensor2d &a, const SharedTensor2d &b, double alpha, double beta);
  double* to_vector(const SharedTensor2i &pair_idx);
  double* to_vector();
  double rms();
  double rms(const SharedTensor2d& a);

  void set_act_oo(int aocc, const SharedTensor2d &a);
  void set_act_oo(int frzc, int aocc, const SharedTensor2d &a);
  void set_act_vv(int occ, int avir, const SharedTensor2d &A);
  void set_act_vv(const SharedTensor2d &A);
  void set_oo(const SharedTensor2d &a);
  void set_ov(const SharedTensor2d &A);
  void set_vo(const SharedTensor2d &A);
  void set_vv(int occ, const SharedTensor2d &A);
  void set_act_ov(int frzc, const SharedTensor2d &A);
  void set_act_vo(int frzc, const SharedTensor2d &A);

  void add_oo(const SharedTensor2d &A, double alpha, double beta);
  void add_vv(int occ, const SharedTensor2d &A, double alpha, double beta);
  void add_ov(const SharedTensor2d &A, double alpha, double beta);
  void add_vo(const SharedTensor2d &A, double alpha, double beta);
  void add_aocc_fc(const SharedTensor2d &A, double alpha, double beta);
  void add_fc_aocc(const SharedTensor2d &A, double alpha, double beta);

  void add3_oo(const SharedTensor2d &A, double alpha, double beta);
  void set3_oo(const SharedTensor2d &A);
  void set3_ov(const SharedTensor2d &A);
  void set3_vo(const SharedTensor2d &A);
  void set3_vv(const SharedTensor2d &A, int occ);

  void set3_act_ov(int frzc, int aocc, int avir, int vir, const SharedTensor2d &a);
  void set3_act_oo(int frzc, const SharedTensor2d &A);
  void set3_act_vv(const SharedTensor2d &A);

  void swap_3index_col(const SharedTensor2d &A);

  void form_oo(const SharedTensor2d &A);
  void form_act_oo(int frzc, const SharedTensor2d &A);
  void form_vv(int occ, const SharedTensor2d &A);
  void form_act_vv(int occ, const SharedTensor2d &A);
  void form_vo(const SharedTensor2d &A);
  void form_vo(int occ, const SharedTensor2d &A);
  void form_act_vo(int frzc, const SharedTensor2d &A);
  void form_act_vo(int frzc, int occ, const SharedTensor2d &A);
  void form_ov(const SharedTensor2d &A);
  void form_ov(int occ, const SharedTensor2d &A);
  void form_act_ov(int frzc, const SharedTensor2d &A);
  void form_act_ov(int frzc, int occ, const SharedTensor2d &A);
  void form_ooAB(const SharedTensor2d &A);

  void form_b_ij(int frzc, const SharedTensor2d &A);
  void form_b_ia(int frzc, const SharedTensor2d &A);
  void form_b_ab(const SharedTensor2d &A);
  // form_b_kl: k is active occupied, and l is frozen core
  void form_b_kl(const SharedTensor2d &A);
  // form_b_ki: k is active occupied, and i is all occupied
  void form_b_ki(const SharedTensor2d &A);
  // form_b_li: l is frozen core, and i is all occupied
  void form_b_li(const SharedTensor2d &A);
  // form_b_il: l is frozen core, and i is all occupied
  void form_b_il(const SharedTensor2d &A);
  // form_b_ka: k is active occupied, and a is all virtual
  void form_b_ka(const SharedTensor2d &A);
  // form_b_la: l is frozen core, and a is all virtual
  void form_b_la(const SharedTensor2d &A);

  // B_pq = 1/2 (A_pq + A_qp)
  void symmetrize();
  void symmetrize(const SharedTensor2d &A);
  // B(Q,pq) = 1/2 [ A(Q,pq) + A(Q,qp) ]
  void symmetrize3(const SharedTensor2d &A);
  // A(Q, p>=q) = A(Q,pq) * (2 -\delta_{pq})
  void symm_packed(const SharedTensor2d &A);
  // A(Q, p>=q) = A(Q,pq)
  void ltm(const SharedTensor2d &A);
  // A(p,qr) = A(p,q>=r)
  void expand23(int d1, int d2, int d3, const SharedTensor2d &A);

  // (+)A(p>=q, r>=s) = 1/2 [A(pq,rs) + A(qp,rs)]
  void symm4(const SharedTensor2d &a);
  // (-)A(p>=q, r>=s) = 1/2 [A(pq,rs) - A(qp,rs)]
  void antisymm4(const SharedTensor2d &a);
  // (+)A(p>=q, r>=s) = 1/2 [A(pq,rs) + A(pq,sr)]
  void symm_col4(const SharedTensor2d &a);
  // (-)A(p>=q, r>=s) = 1/2 [A(pq,rs) - A(pq,sr)]
  void antisymm_col4(const SharedTensor2d &a);
  // (+)At(p>=q, r>=s) = 1/2 [A(pq,rs) + A(qp,rs)] * (2 -\delta_{pq})
  void symm_row_packed4(const SharedTensor2d &a);
  // (+)At(p>=q, r>=s) = 1/2 [A(pq,rs) + A(qp,rs)] * (2 -\delta_{rs})
  void symm_col_packed4(const SharedTensor2d &a);
  // (-)At(p>=q, r>=s) = 1/2 [A(pq,rs) - A(qp,rs)] * (2 -\delta_{pq})
  void antisymm_row_packed4(const SharedTensor2d &a);
  // (-)At(p>=q, r>=s) = 1/2 [A(pq,rs) - A(qp,rs)] * (2 -\delta_{rs})
  void antisymm_col_packed4(const SharedTensor2d &a);

  // A2d_(pq,rs) = 2 <pq|rs> - <pq|sr>
  void tei_cs1_anti_symm(const SharedTensor2d &J, const SharedTensor2d &K);
  // A2d_(pq,rs) = 2 <pq|rs> - <qp|rs>
  void tei_cs2_anti_symm(const SharedTensor2d &J, const SharedTensor2d &K);
  // A2d_(pq,rs) = 2 (pq|rs) - (ps|rq)
  void tei_cs3_anti_symm(const SharedTensor2d &J, const SharedTensor2d &K);
  // A2d_(pq,rs) = 2 (pq|rs) - (rq|ps)
  void tei_cs4_anti_symm(const SharedTensor2d &J, const SharedTensor2d &K);

  // A2d(ij,ab) += P_(ij) * P_(ab) A(ia,jb)
  void P_ijab(const SharedTensor2d &A);


  // General tensor contractions over
  // C(pq,rs) = \sum_{tu} A(pq,tu) B(tu,rs)
  // t_a1: t; t_a2: u; f_a1: p; f_a2: q
  void cont444(int t_a1, int t_a2, int f_a1, int f_a2, const SharedTensor2d& A, int t_b1, int t_b2, int f_b1, int f_b2, const SharedTensor2d& B, double alpha, double beta);
  void cont444(bool delete_a, int t_a1, int t_a2, int f_a1, int f_a2, SharedTensor2d& A,
	     bool delete_b, int t_b1, int t_b2, int f_b1, int f_b2, SharedTensor2d& B,
	     double alpha, double beta);
  void cont444(string idx_c, string idx_a, string idx_b, bool delete_a, bool delete_b, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta);
  // C(pq) = \sum_{rst} A(pr,st) B(rs,tq)
  void cont442(string idx_c, string idx_a, string idx_b, bool delete_a, bool delete_b, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta);
  // C(pq,rs) = \sum_{t} A(pq,rt) B(t,s)
  void cont424(string idx_c, string idx_a, string idx_b, bool delete_a, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta);
  // C(pq,rs) = \sum_{t} A(p,t) B(tq,rs)
  void cont244(string idx_c, string idx_a, string idx_b, bool delete_b, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta);
  // C(Q,pq) = \sum_{rs} A(Q,rs) B(rs,pq)
  // where dim(idx_c) & dim(idx_a)=2 but dim(idx_b)=4
  void cont343(string idx_c, string idx_a, string idx_b, bool delete_b, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta);
  // C(Q,pq) = \sum_{r} A(p,r) B(Q,rq)
  void cont233(string idx_c, string idx_a, string idx_b, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta);
  // C(Q,pq) = \sum_{r} A(Q,pr) B(r,q)
  void cont323(string idx_c, string idx_a, string idx_b, bool delete_a, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta);
  // C(pq) = \sum_{Qr} A(Q,rp) B(Q,rq)
  void cont332(string idx_c, string idx_a, string idx_b, bool delete_a, bool delete_b, SharedTensor2d& A, SharedTensor2d& B, double alpha, double beta);

  friend class Tensor1d;
  friend class Tensor3d;
  friend class Tensor1i;
  friend class Tensor2i;
};

class Tensor3d
{

  private:
  double ***A3d_;
  int dim1_,dim2_,dim3_;
  string name_;      // Name of the array

  public:
  Tensor3d(int d1,int d2, int d3);
  Tensor3d(string name, int d1,int d2, int d3);
  Tensor3d();			   //default constructer
  ~Tensor3d(); 		   	   //destructer

  void init(string name, int d1,int d2, int d3);
  void init(int d1,int d2, int d3);
  void memalloc();
  void zero();
  void print();
  void release();
  void set(int h, int i, int j, double value);
  double get(int h, int i, int j);

  friend class Tensor1d;
  friend class Tensor2d;
};

class Tensor1i
{

  private:
  int *A1i_;
  int dim1_;
  string name_;      // Name of the array

  public:
  Tensor1i(int d1);
  Tensor1i(string name, int d1);
  Tensor1i();			   //default constructer
  ~Tensor1i(); 		   	   //destructer

  void init(string name, int d1);
  void init(int d1);
  void memalloc();
  void zero();
  void print();
  void release();
  void set(int i, int value);
  int get(int i);
  void add(const SharedTensor1i& Adum);
  void add(int i, int value);
  void subtract(const SharedTensor1i& Adum);
  void subtract(int i, int value);

};

class Tensor2i
{

  private:
  int **A2i_;
  int dim1_,dim2_;
  string name_;      // Name of the array

  public:
  Tensor2i(int d1,int d2);
  Tensor2i(string name, int d1,int d2);
  Tensor2i();			   //default constructer
  ~Tensor2i(); 		   	   //destructer

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
  void add(const SharedTensor2i& Adum);
  void add(int i, int j, int value);
  void subtract(const SharedTensor2i& Adum);
  void subtract(int i, int j, int value);
  SharedTensor2i transpose();
  void copy(const SharedTensor2i& Adum);
  void copy(int **a);
  void identity();
  int trace();
  int **to_int_matrix();
  int dim1() const { return dim1_; }
  int dim2() const { return dim2_; }

  friend class Tensor1i;
  friend class Tensor3i;
  friend class Tensor1d;
  friend class Tensor2d;
};


class Tensor3i
{

  private:
  int ***A3i_;
  int dim1_,dim2_,dim3_;
  string name_;      // Name of the array

  public:
  Tensor3i(int d1,int d2, int d3);
  Tensor3i(string name, int d1,int d2, int d3);
  Tensor3i();			   //default constructer
  ~Tensor3i(); 		           //destructer

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
#endif // _dfocc_tensors_h_
