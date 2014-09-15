/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef _dfocc_tensors_h_
#define _dfocc_tensors_h_

#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

#define index2(i,j) ((i>j) ? ((i*(i+1)/2)+j) : ((j*(j+1)/2)+i))
#define index4(i,j,k,l) index2(index2(i,j),index2(k,l))
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define idx_asym(i,j) ((i>j) ? ((i*(i-1)/2)+j) : ((j*(j-1)/2)+i))

using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

class Tensor1d;
class Tensor2d;
class Tensor3d;
class Tensor1i;
class Tensor2i;
class Tensor3i;

typedef boost::shared_ptr<Tensor1d> SharedTensor1d;
typedef boost::shared_ptr<Tensor2d> SharedTensor2d;
typedef boost::shared_ptr<Tensor3d> SharedTensor3d;
typedef boost::shared_ptr<Tensor1i> SharedTensor1i;
typedef boost::shared_ptr<Tensor2i> SharedTensor2i;
typedef boost::shared_ptr<Tensor3i> SharedTensor3i;

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
  Tensor2d(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno, string name, int d1,int d2);
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
  void set(SharedTensor2d &A);
  void set(SharedMatrix A);
  void set(SharedTensor1d &A);
  double get(int i, int j);
  // A2d = alpha * Adum
  void add(const SharedTensor2d &a);
  void add(double **a);
  void add(double alpha, const SharedTensor2d &a);
  void add(int i, int j, double value);
  void subtract(const SharedTensor2d &a);
  void subtract(int i, int j, double value);
  // axpy: Y <-- a * X + Y
  void axpy(double **a, double alpha);
  void axpy(const SharedTensor2d &a, double alpha);
  double **transpose2();
  SharedTensor2d transpose();
  void copy(const SharedTensor2d &Adum);
  void copy(double **a);
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
  // contract323: C[Q](m,n) = \sum_{k} A[Q](m,k) * B(k,n)
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

  void write(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno);
  void write(psi::PSIO* const psio, unsigned int fileno);
  void write(psi::PSIO& psio, unsigned int fileno);
  void write(boost::shared_ptr<psi::PSIO> psio, const string& filename, unsigned int fileno);
  void write(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno, bool three_index, bool symm);
  void write(boost::shared_ptr<psi::PSIO> psio, const string& filename, unsigned int fileno, bool three_index, bool symm);
  void write_symm(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno);
  void write_anti_symm(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno);

  void read(psi::PSIO* psio, unsigned int fileno);
  void read(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno);
  void read(psi::PSIO& psio, unsigned int fileno);
  void read(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno, bool three_index, bool symm);
  void read_symm(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno);
  void read_anti_symm(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno);

  bool read(PSIO* psio, int itap, const char *label, int dim);
  bool read(boost::shared_ptr<psi::PSIO> psio, int itap, const char *label, int dim);
  void save(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno);
  void save(psi::PSIO* const psio, unsigned int fileno);
  void save(psi::PSIO& psio, unsigned int fileno);
  void load(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno, string name, int d1,int d2);
  void load(psi::PSIO* const psio, unsigned int fileno, string name, int d1,int d2);
  void load(psi::PSIO& psio, unsigned int fileno, string name, int d1,int d2);

  // sort (for example 1432 sort): A2d_(ps,rq) = A(pq,rs)
  // A2d_ = alpha*A + beta*A2d_
  void sort(int sort_type, const SharedTensor2d &A, double alpha, double beta);
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

  void set3_oo(const SharedTensor2d &A);
  void add3_oo(const SharedTensor2d &A, double alpha, double beta);
  void set3_act_ov(int frzc, int aocc, int avir, int vir, const SharedTensor2d &a);
  void set3_ov(const SharedTensor2d &A);
  void set3_vo(const SharedTensor2d &A);
  void set3_vv(const SharedTensor2d &A, int occ);
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
  // form_b_la: k is frozen core, and a is all virtual
  void form_b_la(const SharedTensor2d &A);

  // B_pq = 1/2 (A_pq + A_qp)
  void symmetrize();
  // C(Q,pq) = 1/2 [ A(Q,pq) + B(Q,qp) ]
  void symmetrize3(const SharedTensor2d &A);
  // A(Q, p>=q) = A(Q,pq) * (2 -\delta_{pq})
  void symm_packed(const SharedTensor2d &A);
  // A(Q, p>=q) = A(Q,pq)
  void ltm(const SharedTensor2d &A);

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

