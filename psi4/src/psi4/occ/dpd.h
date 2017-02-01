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

#ifndef _psi_src_bin_occ_dpd_h_
#define _psi_src_bin_occ_dpd_h_


using namespace psi;
using namespace std;

namespace psi{ namespace occwave{

class SymBlockVector;

class SymBlockMatrix
{

  private:
  double ***matrix_; // Object
  int *rowspi_;      // Rows per irrep
  int *colspi_;      // Columns per irrep
  std::string name_;      // Name of the matrix
  int nirreps_;     // Number of irreps


  public:
  SymBlockMatrix();  //default constructer
  SymBlockMatrix(std::string name);
  SymBlockMatrix(int nirreps, int *ins_rowspi, int *ins_colspi);
  SymBlockMatrix(std::string name, int nirreps, int *ins_rowspi, int *ins_colspi);
  ~SymBlockMatrix(); //destructer

  SymBlockMatrix* generate(int nirreps, int *ins_rowspi, int *ins_colspi);
  SymBlockMatrix* generate(std::string name, int nirreps, int *ins_rowspi, int *ins_colspi);
  SymBlockMatrix* transpose();
  void init(std::string name, int nirreps, int *ins_rowspi, int *ins_colspi);
  void memalloc();
  void release();
  void zero();
  void zero_diagonal();
  double trace();
  SymBlockMatrix* transpose(int *ins_rowspi, int *ins_colspi);
  void copy(const SymBlockMatrix* Adum);
  void add(const SymBlockMatrix* Adum);
  void add(int h, int i, int j, double value);
  void subtract(const SymBlockMatrix* Adum);
  void subtract(int h, int i, int j, double value);
  void scale(double a);
  void scale_row(int h, int m, double a);
  void scale_column(int h, int n, double a);
  double sum_of_squares();
  double rms();
  double rms(SymBlockMatrix* Atemp);
  void set(double value);
  void set(int h, int i, int j, double value);
  void set(double **Asq);
  void set(dpdbuf4 G);
  double get(int h, int m, int n);
  double *to_lower_triangle();
  double **to_block_matrix();
  void print(std::string OutFileRMR);
  void print();
  void set_to_identity();
  void gemm(bool transa, bool transb, double alpha, const SymBlockMatrix* a, const SymBlockMatrix* b, double beta);
  int *rowspi();
  int *colspi();
  bool load(PSIO* psio, int itap, const char *label, int dim);
  bool load(std::shared_ptr<psi::PSIO> psio, int itap, const char *label, int dim);
  void write(PSIO* psio, int itap, bool saveSubBlocks);
  void write(std::shared_ptr<psi::PSIO> psio, int itap, bool saveSubBlocks);
  void read(std::shared_ptr<psi::PSIO> psio, int itap, bool readSubBlocks);
  void read(std::shared_ptr<psi::PSIO> psio, int itap, const char *label, bool readSubBlocks);
  void mgs();// Modified Gram-Schmidt
  void gs();// Gram-Schmidt
  void diagonalize(SymBlockMatrix* eigvectors, SymBlockVector* eigvalues);
  void cdsyev(char jobz, char uplo, SymBlockMatrix* eigvectors, SymBlockVector* eigvalues); // diagonalize via acml
  void davidson(int n_eigval, SymBlockMatrix* eigvectors, SymBlockVector* eigvalues, double cutoff, int print); // diagonalize via davidson alg.
  void cdgesv(SymBlockVector* Xvec); // solve lineq via acml
  void lineq_flin(SymBlockVector* Xvec, double *det); // solve lineq via flin
  void lineq_pople(SymBlockVector* Xvec, int num_vecs, double cutoff); // solve lineq via pople
  void read_oooo(std::shared_ptr<psi::PSIO> psio, int itap, int *mosym, int *qt2pitzer, int *occ_off, int *occpi, Array3i *oo_pairidx);
  void read_oovv(std::shared_ptr<psi::PSIO> psio, int itap, int nocc, int *mosym, int *qt2pitzer, int *occ_off, int *vir_off, int *occpi,
                 int *virpi, Array3i *oo_pairidx, Array3i *vv_pairidx);


  friend class SymBlockVector;
};


class SymBlockVector
{

  private:
  double **vector_; // Object
  int *dimvec_;      // dimensions per irrep
  string name_;      // Name of the vector
  int nirreps_;     // Number of irreps


  public:
  SymBlockVector();  //default constructer
  SymBlockVector(string name);
  SymBlockVector(int nirreps, int *ins_dimvec);
  SymBlockVector(string name, int nirreps, int *ins_dimvec);
  ~SymBlockVector(); //destructer

  int *dimvec();
  void memalloc();
  void release();
  void zero();
  double trace();
  void copy(const SymBlockVector* Adum);
  void add(const SymBlockVector* Adum);
  void add(int h, int i, double value);
  void subtract(const SymBlockVector* Adum);
  void subtract(int h, int i, double value);
  void scale(double a);
  double sum_of_squares();
  double rms();
  double rms(SymBlockVector* Atemp);
  double norm();
  void set(double value);
  void set(int h, int i, double value);
  void set(double *Avec);
  double get(int h, int m);
  double *to_vector();
  void print(std::string OutFileRMR);
  void print();
  void set_to_unit();
  void gemv(bool transa, double alpha, SymBlockMatrix* A, SymBlockVector* X, double beta);
  double dot(SymBlockVector* X);

  friend class SymBlockMatrix;
};

}} // End Namespaces
#endif // _psi_src_bin_occ_dpd_h_
