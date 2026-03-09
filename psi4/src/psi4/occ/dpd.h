/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#ifndef _psi_src_bin_occ_dpd_h_
#define _psi_src_bin_occ_dpd_h_

using namespace psi;

namespace psi {
namespace occwave {

class SymBlockVector;

class SymBlockMatrix {
   private:
    double ***matrix_;  // Object
    int *rowspi_;       // Rows per irrep
    int *colspi_;       // Columns per irrep
    std::string name_;  // Name of the matrix
    int nirreps_;       // Number of irreps

   public:
    SymBlockMatrix();  // default constructer
    SymBlockMatrix(std::string name, int nirreps, int *ins_rowspi, int *ins_colspi);
    ~SymBlockMatrix();  // destructer

    void memalloc();
    void release();
    void zero();
    void set(dpdbuf4 G);
    double get(int h, int m, int n);

    friend class SymBlockVector;
};

class SymBlockVector {
   private:
    double **vector_;   // Object
    int *dimvec_;       // dimensions per irrep
    std::string name_;  // Name of the vector
    int nirreps_;       // Number of irreps

   public:
    SymBlockVector();  // default constructer
    SymBlockVector(std::string name);
    SymBlockVector(int nirreps, int *ins_dimvec);
    SymBlockVector(std::string name, int nirreps, int *ins_dimvec);
    ~SymBlockVector();  // destructer

    int *dimvec();
    void memalloc();
    void release();
    void zero();
    double trace();
    void copy(const SymBlockVector *Adum);
    void add(const SymBlockVector *Adum);
    void add(int h, int i, double value);
    void subtract(const SymBlockVector *Adum);
    void subtract(int h, int i, double value);
    void scale(double a);
    double sum_of_squares();
    double rms();
    double rms(SymBlockVector *Atemp);
    double norm();
    void set(double value);
    void set(int h, int i, double value);
    void set(double *Avec);
    double get(int h, int m);
    double *to_vector();
    void print(std::string out_fname);
    void print();
    void set_to_unit();
    void gemv(bool transa, double alpha, SymBlockMatrix *A, SymBlockVector *X, double beta);
    double dot(SymBlockVector *X);

    friend class SymBlockMatrix;
};
}
}  // End Namespaces
#endif  // _psi_src_bin_occ_dpd_h_
