/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

// #ifndef _dfocc_tensors_h_
// #define _dfocc_tensors_h_

#include "psi4/libpsio/psio.h"
#include "psi4/libmints/typedefs.h"
#include "tensors.h"

#define index2(i, j) ((i > j) ? ((i * (i + 1) / 2) + j) : ((j * (j + 1) / 2) + i))
#define index4(i, j, k, l) index2(index2(i, j), index2(k, l))
#define MIN0(a, b) (((a) < (b)) ? (a) : (b))
#define MAX0(a, b) (((a) > (b)) ? (a) : (b))
#define idx_asym(i, j) ((i > j) ? ((i * (i - 1) / 2) + j) : ((j * (j - 1) / 2) + i))

using namespace psi;
using namespace dfoccwave;

namespace psi {

class PSIO;

namespace dfoccwave {

class Tensor1f;
class Tensor2f;
class Tensor3f;

typedef std::shared_ptr<Tensor1f> SharedTensor1f;
typedef std::shared_ptr<Tensor2f> SharedTensor2f;
typedef std::shared_ptr<Tensor3f> SharedTensor3f;

class Tensor1f {
   private:
    float *A1d_;
    int dim1_;
    std::string name_;  // Name of the array

   public:
    Tensor1f(int d1);
    Tensor1f(std::string name, int d1);
    Tensor1f();   // default constructer
    ~Tensor1f();  // destructer

    void init(std::string name, int d1);
    void init(int d1);
    void memalloc();
    void zero();
    void print();
    void print(std::string out_fname);
    void print(const char *outfile);
    void print(FILE *out);
    void release();
    void set(int i, float value);
    void set(float *vec);
    void set(const SharedTensor1f &vec);
    void add(const SharedTensor1f &Adum);
    void add(int i, float value);  // add value to ith element of the vector
    void subtract(const SharedTensor1f &Adum);
    void subtract(int i, float value);
    float get(int i);
    void to_shared_vector(SharedVector A);
    // rms:  rms of A1d_
    float rms();
    // rms:  rms of (A1d_ - Atemp)
    float rms(const SharedTensor1f &Atemp);
    // dot: return result of A1d_' * y
    float dot(const SharedTensor1f &y);
    // gemv: A1d_ = alpha * A * b + beta, where A is a general matrix
    // gemv: C(m) = \sum_{n} A(m,n) b(n)
    void gemv(bool transa, const SharedTensor2f &a, const SharedTensor1f &b, float alpha, float beta);
    // gemv: C(m) = \sum_{n} A(m,n) b(n)
    void gemv(bool transa, int m, int n, const SharedTensor2f &a, const SharedTensor2f &b, float alpha, float beta);
    void gemv(bool transa, const SharedTensor2f &a, const SharedTensor2f &b, float alpha, float beta);
    void gemv(bool transa, int m, int n, const SharedTensor2f &a, const SharedTensor2f &b, int start_a, int start_b,
              float alpha, float beta);
    void gemv(bool transa, int m, int n, const SharedTensor2f &a, const SharedTensor2f &b, int start_a, int start_b,
              int start_c, float alpha, float beta);
    // gbmv: This function may NOT working correctly!!!!
    void gbmv(bool transa, const SharedTensor2f &a, const SharedTensor1f &b, float alpha, float beta);
    // xay: return result of A1d_' * A * y
    float xay(const SharedTensor2f &a, const SharedTensor1f &y);
    // axpy: Y <-- a * X + Y
    void axpy(const SharedTensor1f &a, float alpha);
    void scale(float a);
    void copy(float *x);
    void copy(const SharedTensor1f &x);
    // row_vector: set A1d to nth row of A, dim1_ = A->dim2
    void row_vector(SharedTensor2f &A, int n);
    // column_vector: set A1d to nth column of A, dim1_ = A->dim1
    void column_vector(SharedTensor2f &A, int n);
    int dim1() const { return dim1_; }
    // dirprd: A1d_[i] = a[i] * b[i]
    void dirprd(SharedTensor1f &a, SharedTensor1f &b);
    // symm_packed A(p>=q) = A(p,q) * (2 - \delta_{pq})
    void symm_packed(const SharedTensor2f &A);
    // ltm: lower triangular matrix: A(p>=q) = A(p,q)
    void ltm(const SharedTensor2f &A);

    friend class Tensor2f;
    friend class Tensor3f;
};

class Tensor2f {
   private:
    float **A2d_;
    int dim1_, dim2_, d1_, d2_, d3_, d4_;
    int **row_idx_, **col_idx_;
    int *row2d1_, *row2d2_, *col2d1_, *col2d2_;
    std::string name_;  // Name of the array

   public:
    Tensor2f(int d1, int d2);
    Tensor2f(std::string name, int d1, int d2);
    Tensor2f(psi::PSIO *psio, size_t fileno, std::string name, int d1, int d2);
    Tensor2f(std::shared_ptr<psi::PSIO> psio, size_t fileno, std::string name, int d1, int d2);
    Tensor2f(psi::PSIO &psio, size_t fileno, std::string name, int d1, int d2);
    Tensor2f(std::string name, int d1, int d2, int d3, int d4);
    Tensor2f(std::string name, int d1, int d2, int d3);
    Tensor2f();   // default constructer
    ~Tensor2f();  // destructer

    void init(std::string name, int d1, int d2);
    void init(int d1, int d2);
    void memalloc();
    void zero();
    void zero_diagonal();
    void zero_off_diagonal();
    void print();
    void print(std::string out_fname);
    void print(const char *outfile);
    void print(FILE *out);
    void release();
    void set(int i, int j, float value);
    void set(float **A);
    void set(float *A);
    void set(SharedTensor2f &A);
    void set(SharedMatrix A);
    void set2(SharedMatrix A);
    void set(SharedTensor1f &A);
    // A2d_[n][ij] = A(i,j)
    void set_row(const SharedTensor2f &A, int n);
    // A2d_[ij][n] = A(i,j)
    void set_column(const SharedTensor2f &A, int n);
    float get(int i, int j);
    // A2d_[ij] = A(n, ij)
    void get_row(const SharedTensor2f &A, int n);
    // A2d_[ij] = A(ij, n)
    void get_column(const SharedTensor2f &A, int n);
    // A2d = alpha * Adum
    void add(const SharedTensor2f &a);
    void add(float **a);
    void add(float alpha, const SharedTensor2f &a);
    void add(int i, int j, float value);
    // A2d_[n][ij] += A(i,j)
    void add2row(const SharedTensor2f &A, int n);
    // A2d_[ij][n] += A(i,j)
    void add2col(const SharedTensor2f &A, int n);
    void subtract(const SharedTensor2f &a);
    void subtract(int i, int j, float value);
    // axpy: Y <-- a * X + Y
    void axpy(float **a, float alpha);
    void axpy(const SharedTensor2f &a, float alpha);
    void axpy(size_t length, int inc_a, const SharedTensor2f &a, int inc_2d, float alpha);
    void axpy(size_t length, int start_a, int inc_a, const SharedTensor2f &A, int start_2d, int inc_2d, float alpha);
    float **transpose2();
    SharedTensor2f transpose();
    void trans(const SharedTensor2f &A);
    void trans(float **A);
    void copy(float **a);
    void copy(const SharedTensor2f &Adum);
    void copy(size_t length, const SharedTensor2f &A, int inc_a, int inc_2d);
    void copy(const SharedTensor2f &A, int start);
    // partial copy
    void pcopy(const SharedTensor2f &A, int dim_copy, int dim_skip);
    void pcopy(const SharedTensor2f &A, int dim_copy, int dim_skip, int start);
    float get_max_element();
    // diagonalize: diagonalize via rsp
    void diagonalize(const SharedTensor2f &eigvectors, const SharedTensor1f &eigvalues, float cutoff);
    // cdsyev: diagonalize via lapack
    void cdsyev(char jobz, char uplo, const SharedTensor2f &eigvectors, const SharedTensor1f &eigvalues);
    // davidson: diagonalize via davidson algorithm
    void davidson(int n_eigval, const SharedTensor2f &eigvectors, const SharedTensor1f &eigvalues, float cutoff,
                  int print);
    // cdgesv: solve a linear equation via lapack
    void cdgesv(const SharedTensor1f &Xvec);
    void cdgesv(float *Xvec);
    void cdgesv(const SharedTensor1f &Xvec, int errcod);
    void cdgesv(float *Xvec, int errcod);
    // lineq_flin: solve a linear equation via FLIN
    void lineq_flin(const SharedTensor1f &Xvec, float *det);
    void lineq_flin(float *Xvec, float *det);
    // pople: solve a linear equation via Pople's algorithm
    void lineq_pople(const SharedTensor1f &Xvec, int num_vecs, float cutoff);
    void lineq_pople(float *Xvec, int num_vecs, float cutoff);

    // gemm: matrix multiplication C = A * B
    void gemm(bool transa, bool transb, const SharedTensor2f &a, const SharedTensor2f &b, float alpha, float beta);
    // contract: general contraction C(m,n) = \sum_{k} A(m,k) * B(k,n)
    void contract(bool transa, bool transb, int m, int n, int k, const SharedTensor2f &a, const SharedTensor2f &b,
                  float alpha, float beta);
    void contract(bool transa, bool transb, int m, int n, int k, const SharedTensor2f &a, const SharedTensor2f &b,
                  int start_a, int start_b, float alpha, float beta);
    void contract(bool transa, bool transb, int m, int n, int k, const SharedTensor2f &a, const SharedTensor2f &b,
                  int start_a, int start_b, int start_c, float alpha, float beta);
    // contract323: C[Q](m,n) = \sum_{k} A[Q](m,k) * B(k,n). Note: contract332 should be called with beta=1.0
    void contract323(bool transa, bool transb, int m, int n, const SharedTensor2f &a, const SharedTensor2f &b,
                     float alpha, float beta);
    // contract233: C[Q](m,n) = \sum_{k} A(m,k) * B[Q](k,n)
    void contract233(bool transa, bool transb, int m, int n, const SharedTensor2f &a, const SharedTensor2f &b,
                     float alpha, float beta);
    // contract332: C(m,n) = \sum_{k} A[Q](m,k) * B[Q](k,n)
    void contract332(bool transa, bool transb, int k, const SharedTensor2f &a, const SharedTensor2f &b, float alpha,
                     float beta);
    // contract424: Z(pq,rs) = \sum_{o} X(pq,ro) * Y(o,s): where target_x = 4, target_y = 1, and Z = A2d_
    void contract424(int target_x, int target_y, const SharedTensor2f &a, const SharedTensor2f &b, float alpha,
                     float beta);
    // contract442: C(p,q) \sum_{rst} A(pr,st) B(qr,st) , where row/col pair indices are related to B and sorted B.
    void contract442(int target_a, int target_b, const SharedTensor2f &a, const SharedTensor2f &b, float alpha,
                     float beta);
    // gemv: C(mn) = \sum_{k} A(mn,k) * b(k)
    void gemv(bool transa, const SharedTensor2f &a, const SharedTensor1f &b, float alpha, float beta);
    // gemv: C(mn) = \sum_{kl} A(mn,kl) * b(kl)
    void gemv(bool transa, const SharedTensor2f &a, const SharedTensor2f &b, float alpha, float beta);

    // level_shift: A[i][i] = A[i][i] - value
    void level_shift(float value);
    // outer_product: A = x * y'
    void outer_product(const SharedTensor1f &x, const SharedTensor1f &y);
    void scale(float a);
    // scale_row: scales mth row with a
    void scale_row(int m, float a);
    // scale_column: scales nth column with a
    void scale_column(int n, float a);
    // identity: A = I
    void identity();
    float trace();
    float norm();
    // transform: A = L' * B * L
    void transform(const SharedTensor2f &a, const SharedTensor2f &transformer);
    // back_transform: A = L * B * L'
    void back_transform(const SharedTensor2f &a, const SharedTensor2f &transformer);
    // back_transform: A = alpha* L * B * L' + (beta*A)
    void back_transform(const SharedTensor2f &a, const SharedTensor2f &transformer, float alpha, float beta);
    // pseudo_transform: A = L * B * L
    void pseudo_transform(const SharedTensor2f &a, const SharedTensor2f &transformer);
    // triple_gemm: A2d_ = a * b * c
    void triple_gemm(const SharedTensor2f &a, const SharedTensor2f &b, const SharedTensor2f &c);
    // vector_dot: value = Tr(A' * B)
    float vector_dot(const SharedTensor2f &rhs);
    float vector_dot(float **rhs);
    float **to_block_matrix();
    float *to_lower_triangle();
    void to_shared_matrix(SharedMatrix A);
    void to_matrix(SharedMatrix A);
    void to_pointer(float *A);
    // mgs: orthogonalize with a Modified Gram-Schmid algorithm
    void mgs();
    // gs: orthogonalize with a Classical Gram-Schmid algorithm
    void gs();
    // row_vector: return nth row as a vector
    float *row_vector(int n);
    // column_vector: return nth column as a vector
    float *column_vector(int n);
    int dim1() const { return dim1_; }
    int dim2() const { return dim2_; }

    void write(std::shared_ptr<psi::PSIO> psio, size_t fileno);
    void write(std::shared_ptr<psi::PSIO> psio, size_t fileno, psio_address start, psio_address *end);
    void write(psi::PSIO *const psio, size_t fileno);
    void write(psi::PSIO *psio, size_t fileno, psio_address start, psio_address *end);
    void write(psi::PSIO &psio, size_t fileno);
    void write(psi::PSIO &psio, size_t fileno, psio_address start, psio_address *end);
    void write(std::shared_ptr<psi::PSIO> psio, const std::string &filename, size_t fileno);
    void write(std::shared_ptr<psi::PSIO> psio, size_t fileno, bool three_index, bool symm);
    void write(std::shared_ptr<psi::PSIO> psio, const std::string &filename, size_t fileno, bool three_index,
               bool symm);
    void write_symm(std::shared_ptr<psi::PSIO> psio, size_t fileno);
    void write_anti_symm(std::shared_ptr<psi::PSIO> psio, size_t fileno);

    void read(psi::PSIO *psio, size_t fileno);
    void read(psi::PSIO *psio, size_t fileno, psio_address start, psio_address *end);
    void read(std::shared_ptr<psi::PSIO> psio, size_t fileno);
    void read(std::shared_ptr<psi::PSIO> psio, size_t fileno, psio_address start, psio_address *end);
    void read(psi::PSIO &psio, size_t fileno);
    void read(psi::PSIO &psio, size_t fileno, psio_address start, psio_address *end);
    void read(std::shared_ptr<psi::PSIO> psio, size_t fileno, bool three_index, bool symm);
    void read_symm(std::shared_ptr<psi::PSIO> psio, size_t fileno);
    void read_anti_symm(std::shared_ptr<psi::PSIO> psio, size_t fileno);

    bool read(PSIO *psio, int itap, const char *label, int dim);
    bool read(std::shared_ptr<psi::PSIO> psio, int itap, const char *label, int dim);
    void save(std::shared_ptr<psi::PSIO> psio, size_t fileno);
    void save(psi::PSIO *const psio, size_t fileno);
    void save(psi::PSIO &psio, size_t fileno);
    void load(std::shared_ptr<psi::PSIO> psio, size_t fileno, std::string name, int d1, int d2);
    void load(psi::PSIO *const psio, size_t fileno, std::string name, int d1, int d2);
    void load(psi::PSIO &psio, size_t fileno, std::string name, int d1, int d2);

    void mywrite(const std::string &filename);
    void mywrite(int fileno);
    void mywrite(int fileno, bool append);

    void myread(const std::string &filename);
    void myread(int fileno);
    void myread(int fileno, bool append);
    void myread(int fileno, size_t start);

    // sort (for example 1432 sort): A2d_(ps,rq) = A(pq,rs)
    // A2d_ = alpha*A + beta*A2d_
    void sort(int sort_type, const SharedTensor2f &A, float alpha, float beta);
    // A2d_[p][qr] = sort(A[p][qr])
    void sort3a(int sort_type, int d1, int d2, int d3, const SharedTensor2f &A, float alpha, float beta);
    // A2d_[pq][r] = sort(A[pq][r])
    void sort3b(int sort_type, int d1, int d2, int d3, const SharedTensor2f &A, float alpha, float beta);
    // apply_denom: T(ij,ab) /= D(ij,ab)
    void apply_denom(int frzc, int occ, const SharedTensor2f &fock);
    // apply_denom_os: T(Ij,Ab) /= D(Ij,Ab)
    void apply_denom_os(int frzc, int occA, int occB, const SharedTensor2f &fockA, const SharedTensor2f &fockB);
    // apply_denom_chem: T(ia,jb) /= D(ij,ab)
    void apply_denom_chem(int frzc, int occ, const SharedTensor2f &fock);
    void reg_denom(int frzc, int occ, const SharedTensor2f &fock, float reg);
    void reg_denom_os(int frzc, int occA, int occB, const SharedTensor2f &fockA, const SharedTensor2f &fockB,
                      float reg);
    void reg_denom_chem(int frzc, int occ, const SharedTensor2f &fock, float reg);
    // dirprd: A2d_[i][j] = a[i][j] * b[i][j]
    void dirprd(const SharedTensor2f &a, const SharedTensor2f &b);
    // dirprd123: A2d_[Q][ij] = a[Q] * b[i][j]
    void dirprd123(const SharedTensor1f &a, const SharedTensor2f &b, float alpha, float beta);
    void dirprd123(bool transb, const SharedTensor1f &a, const SharedTensor2f &b, float alpha, float beta);
    // dirprd112: A2d_[i][j] = a[i] * b[j]
    void dirprd112(const SharedTensor1f &a, const SharedTensor1f &b);
    // dirprd112: A2d_[i][j] = alpha *a[i] * b[j] + beta * A2d_[i][j]
    void dirprd112(const SharedTensor1f &a, const SharedTensor1f &b, float alpha, float beta);
    // dirprd224: A2d_[ij][kl] = a[i][j] * b[k][l]
    void dirprd224(const SharedTensor2f &a, const SharedTensor2f &b);
    void dirprd224(const SharedTensor2f &a, const SharedTensor2f &b, float alpha, float beta);
    float *to_vector(const SharedTensor2i &pair_idx);
    float *to_vector();
    float rms();
    float rms(const SharedTensor2f &a);

    void set_act_oo(int aocc, const SharedTensor2f &a);
    void set_act_oo(int frzc, int aocc, const SharedTensor2f &a);
    void set_act_vv(int occ, int avir, const SharedTensor2f &A);
    void set_act_vv(const SharedTensor2f &A);
    void set_oo(const SharedTensor2f &a);
    void set_ov(const SharedTensor2f &A);
    void set_vo(const SharedTensor2f &A);
    void set_vv(int occ, const SharedTensor2f &A);
    void set_act_ov(int frzc, const SharedTensor2f &A);
    void set_act_vo(int frzc, const SharedTensor2f &A);

    void add_oo(const SharedTensor2f &A, float alpha, float beta);
    void add_vv(int occ, const SharedTensor2f &A, float alpha, float beta);
    void add_ov(const SharedTensor2f &A, float alpha, float beta);
    void add_vo(const SharedTensor2f &A, float alpha, float beta);
    void add_aocc_fc(const SharedTensor2f &A, float alpha, float beta);
    void add_fc_aocc(const SharedTensor2f &A, float alpha, float beta);

    void add3_oo(const SharedTensor2f &A, float alpha, float beta);
    void set3_oo(const SharedTensor2f &A);
    void set3_ov(const SharedTensor2f &A);
    void set3_vo(const SharedTensor2f &A);
    void set3_vv(const SharedTensor2f &A, int occ);

    void set3_act_ov(int frzc, int aocc, int avir, int vir, const SharedTensor2f &a);
    void set3_act_oo(int frzc, const SharedTensor2f &A);
    void set3_act_vv(const SharedTensor2f &A);

    void swap_3index_col(const SharedTensor2f &A);

    void form_oo(const SharedTensor2f &A);
    void form_act_oo(int frzc, const SharedTensor2f &A);
    void form_vv(int occ, const SharedTensor2f &A);
    void form_act_vv(int occ, const SharedTensor2f &A);
    void form_vo(const SharedTensor2f &A);
    void form_vo(int occ, const SharedTensor2f &A);
    void form_act_vo(int frzc, const SharedTensor2f &A);
    void form_act_vo(int frzc, int occ, const SharedTensor2f &A);
    void form_ov(const SharedTensor2f &A);
    void form_ov(int occ, const SharedTensor2f &A);
    void form_act_ov(int frzc, const SharedTensor2f &A);
    void form_act_ov(int frzc, int occ, const SharedTensor2f &A);
    void form_ooAB(const SharedTensor2f &A);

    void form_b_ij(int frzc, const SharedTensor2f &A);
    void form_b_ia(int frzc, const SharedTensor2f &A);
    void form_b_ab(const SharedTensor2f &A);
    // form_b_kl: k is active occupied, and l is frozen core
    void form_b_kl(const SharedTensor2f &A);
    // form_b_ki: k is active occupied, and i is all occupied
    void form_b_ki(const SharedTensor2f &A);
    // form_b_li: l is frozen core, and i is all occupied
    void form_b_li(const SharedTensor2f &A);
    // form_b_il: l is frozen core, and i is all occupied
    void form_b_il(const SharedTensor2f &A);
    // form_b_ka: k is active occupied, and a is all virtual
    void form_b_ka(const SharedTensor2f &A);
    // form_b_la: l is frozen core, and a is all virtual
    void form_b_la(const SharedTensor2f &A);

    // B_pq = 1/2 (A_pq + A_qp)
    void symmetrize();
    void symmetrize(const SharedTensor2f &A);
    // B(Q,pq) = 1/2 [ A(Q,pq) + A(Q,qp) ]
    void symmetrize3(const SharedTensor2f &A);
    // A(Q, p>=q) = A(Q,pq) * (2 -\delta_{pq})
    void symm_packed(const SharedTensor2f &A);
    // A(Q, p>=q) = A(Q,pq)
    void ltm(const SharedTensor2f &A);
    // A(p,qr) = A(p,q>=r)
    void expand23(int d1, int d2, int d3, const SharedTensor2f &A);

    // (+)A(p>=q, r>=s) = 1/2 [A(pq,rs) + A(qp,rs)]
    void symm4(const SharedTensor2f &a);
    // (-)A(p>=q, r>=s) = 1/2 [A(pq,rs) - A(qp,rs)]
    void antisymm4(const SharedTensor2f &a);
    // (+)A(p>=q, r>=s) = 1/2 [A(pq,rs) + A(pq,sr)]
    void symm_col4(const SharedTensor2f &a);
    // (-)A(p>=q, r>=s) = 1/2 [A(pq,rs) - A(pq,sr)]
    void antisymm_col4(const SharedTensor2f &a);
    // (+)At(p>=q, r>=s) = 1/2 [A(pq,rs) + A(qp,rs)] * (2 -\delta_{pq})
    void symm_row_packed4(const SharedTensor2f &a);
    // (+)At(p>=q, r>=s) = 1/2 [A(pq,rs) + A(qp,rs)] * (2 -\delta_{rs})
    void symm_col_packed4(const SharedTensor2f &a);
    // (-)At(p>=q, r>=s) = 1/2 [A(pq,rs) - A(qp,rs)] * (2 -\delta_{pq})
    void antisymm_row_packed4(const SharedTensor2f &a);
    // (-)At(p>=q, r>=s) = 1/2 [A(pq,rs) - A(qp,rs)] * (2 -\delta_{rs})
    void antisymm_col_packed4(const SharedTensor2f &a);

    // A2d_(pq,rs) = 2 <pq|rs> - <pq|sr>
    void tei_cs1_anti_symm(const SharedTensor2f &J, const SharedTensor2f &K);
    // A2d_(pq,rs) = 2 <pq|rs> - <qp|rs>
    void tei_cs2_anti_symm(const SharedTensor2f &J, const SharedTensor2f &K);
    // A2d_(pq,rs) = 2 (pq|rs) - (ps|rq)
    void tei_cs3_anti_symm(const SharedTensor2f &J, const SharedTensor2f &K);
    // A2d_(pq,rs) = 2 (pq|rs) - (rq|ps)
    void tei_cs4_anti_symm(const SharedTensor2f &J, const SharedTensor2f &K);

    // A2d(ij,ab) += P_(ij) * P_(ab) A(ia,jb)
    void P_ijab(const SharedTensor2f &A);

    // General tensor contractions over
    // C(pq,rs) = \sum_{tu} A(pq,tu) B(tu,rs)
    // t_a1: t; t_a2: u; f_a1: p; f_a2: q
    void cont444(int t_a1, int t_a2, int f_a1, int f_a2, const SharedTensor2f &A, int t_b1, int t_b2, int f_b1,
                 int f_b2, const SharedTensor2f &B, float alpha, float beta);
    void cont444(bool delete_a, int t_a1, int t_a2, int f_a1, int f_a2, SharedTensor2f &A, bool delete_b, int t_b1,
                 int t_b2, int f_b1, int f_b2, SharedTensor2f &B, float alpha, float beta);
    void cont444(std::string idx_c, std::string idx_a, std::string idx_b, bool delete_a, bool delete_b,
                 SharedTensor2f &A, SharedTensor2f &B, float alpha, float beta);
    // Tensors A and B will be Deleted!
    void cont444(std::string idx_c, std::string idx_a, std::string idx_b, SharedTensor2f &A, SharedTensor2f &B,
                 float alpha, float beta);
    // C(pq) = \sum_{rst} A(pr,st) B(rs,tq)
    void cont442(std::string idx_c, std::string idx_a, std::string idx_b, bool delete_a, bool delete_b,
                 SharedTensor2f &A, SharedTensor2f &B, float alpha, float beta);
    // C(pq,rs) = \sum_{t} A(pq,rt) B(t,s)
    void cont424(std::string idx_c, std::string idx_a, std::string idx_b, bool delete_a, SharedTensor2f &A,
                 SharedTensor2f &B, float alpha, float beta);
    // C(pq,rs) = \sum_{t} A(p,t) B(tq,rs)
    void cont244(std::string idx_c, std::string idx_a, std::string idx_b, bool delete_b, SharedTensor2f &A,
                 SharedTensor2f &B, float alpha, float beta);
    // C(Q,pq) = \sum_{rs} A(Q,rs) B(rs,pq)
    // where dim(idx_c) & dim(idx_a)=2 but dim(idx_b)=4
    void cont343(std::string idx_c, std::string idx_a, std::string idx_b, bool delete_b, SharedTensor2f &A,
                 SharedTensor2f &B, float alpha, float beta);
    // C(Q,pq) = \sum_{r} A(p,r) B(Q,rq)
    void cont233(std::string idx_c, std::string idx_a, std::string idx_b, SharedTensor2f &A, SharedTensor2f &B,
                 float alpha, float beta);
    // C(Q,pq) = \sum_{r} A(Q,pr) B(r,q)
    void cont323(std::string idx_c, std::string idx_a, std::string idx_b, bool delete_a, SharedTensor2f &A,
                 SharedTensor2f &B, float alpha, float beta);
    // C(pq) = \sum_{Qr} A(Q,rp) B(Q,rq)
    void cont332(std::string idx_c, std::string idx_a, std::string idx_b, bool delete_a, bool delete_b,
                 SharedTensor2f &A, SharedTensor2f &B, float alpha, float beta);

    friend class Tensor1f;
    friend class Tensor3f;
    friend class Tensor1i;
    friend class Tensor2i;
};

class Tensor3f {
   private:
    float ***A3d_;
    int dim1_, dim2_, dim3_;
    std::string name_;  // Name of the array

   public:
    Tensor3f(int d1, int d2, int d3);
    Tensor3f(std::string name, int d1, int d2, int d3);
    Tensor3f();   // default constructer
    ~Tensor3f();  // destructer

    void init(std::string name, int d1, int d2, int d3);
    void init(int d1, int d2, int d3);
    void memalloc();
    void zero();
    void print();
    void release();
    void set(int h, int i, int j, float value);
    float get(int h, int i, int j);

    friend class Tensor1f;
    friend class Tensor2f;
};

}  // namespace dfoccwave
}  // namespace psi
// #endif  // _dfocc_tensors_h_
