/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#ifndef _psi_src_bin_occ_arrays_h_
#define _psi_src_bin_occ_arrays_h_

#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.h"

using namespace psi;

namespace psi {
namespace occwave {

class Array1d;
class Array2d;
class Array3d;
class Array1i;
class Array3i;

class Array1d {
   private:
    double* A1d_;
    int dim1_;
    std::string name_;  // Name of the array

   public:
    Array1d(int d1);
    Array1d(std::string name, int d1);
    Array1d();   // default constructer
    ~Array1d();  // destructer

    Array1d* generate(int d1);
    Array1d* generate(std::string name, int d1);
    void init(std::string name, int d1);
    void init(int d1);
    void memalloc();
    void zero();
    void print();
    void print(std::string out_fname);
    void release();
    void set(int i, double value);
    void set(double* vec);
    void set(const Array1d* vec);
    void add(const Array1d* Adum);
    void add(int i, double value);  // add value to ith element of the vector
    void subtract(const Array1d* Adum);
    void subtract(int i, double value);
    double get(int i);
    // rms:  rms of A1d_
    double rms();
    // rms:  rms of (A1d_ - Atemp)
    double rms(const Array1d* Atemp);
    // dot: return result of A1d_' * y
    double dot(const Array1d* y);
    // gemv: A1d_ = alpha * A * b + beta, where A is a general matrix
    void gemv(bool transa, double alpha, const Array2d* a, const Array1d* b, double beta);
    // xay: return result of A1d_' * A * y
    double xay(const Array2d* a, const Array1d* y);
    void scale(double a);
    void copy(double* x);
    void copy(const Array1d* x);
    int dim1() const { return dim1_; }
    // dirprd: A1d_[i] = a[i] * b[i]
    void dirprd(Array1d* a, Array1d* b);
    std::string name() const { return name_; }
    double* nonconst_array() const { return A1d_; }
    const double* array() const { return A1d_; }

    friend class Array2d;
    friend class Array3d;
};

class Array2d {
   private:
    double** A2d_;
    int dim1_, dim2_;
    std::string name_;  // Name of the array

   public:
    Array2d(int d1, int d2);
    Array2d(std::string name, int d1, int d2);
    Array2d();   // default constructer
    ~Array2d();  // destructer

    Array2d* generate(int d1, int d2);
    Array2d* generate(std::string name, int d1, int d2);
    void init(std::string name, int d1, int d2);
    void init(int d1, int d2);
    void memalloc();
    void zero();
    void zero_diagonal();
    void print();
    void print(std::string out_fname);
    void release();
    void set(int i, int j, double value);
    void set(double** A);
    double get(int i, int j);
    void add(const Array2d* Adum);
    void add(int i, int j, double value);
    void subtract(const Array2d* Adum);
    void subtract(int i, int j, double value);
    Array2d* transpose();
    void copy(const Array2d* Adum);
    void copy(double** a);
    // cdgesv: solve a linear equation via lapack
    void cdgesv(Array1d* Xvec, int errcod);
    // gemm: matrix multiplication
    void gemm(bool transa, bool transb, double alpha, const Array2d* a, const Array2d* b, double beta);
    /*
    void write(psi::PSIO* psio, size_t fileno);
    void write(shared_ptr<psi::PSIO> psio, size_t fileno);
    void write(psi::PSIO& psio, size_t fileno);
    void read(psi::PSIO* psio, size_t fileno);
    void read(psi::PSIO& psio, size_t fileno);
    bool read(PSIO* psio, int itap, const char *label, int dim);
    bool read(shared_ptr<psi::PSIO> psio, int itap, const char *label, int dim);
    */

    friend class Array1d;
    friend class Array3d;
};

class Array3d {
   private:
    double*** A3d_;
    int dim1_, dim2_, dim3_;
    std::string name_;  // Name of the array

   public:
    Array3d(int d1, int d2, int d3);
    Array3d(std::string name, int d1, int d2, int d3);
    Array3d();   // default constructer
    ~Array3d();  // destructer

    Array3d* generate(int d1, int d2, int d3);
    Array3d* generate(std::string name, int d1, int d2, int d3);
    void init(std::string name, int d1, int d2, int d3);
    void init(int d1, int d2, int d3);
    void memalloc();
    void zero();
    void print();
    void release();
    void set(int h, int i, int j, double value);
    double get(int h, int i, int j);

    friend class Array1d;
    friend class Array2d;
};

class Array1i {
   private:
    int* A1i_;
    int dim1_;
    std::string name_;  // Name of the array

   public:
    Array1i(int d1);
    Array1i(std::string name, int d1);
    Array1i();   // default constructer
    ~Array1i();  // destructer

    Array1i* generate(int d1);
    Array1i* generate(std::string name, int d1);
    void init(std::string name, int d1);
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

class Array3i {
   private:
    int*** A3i_;
    int dim1_, dim2_, dim3_;
    std::string name_;  // Name of the array

   public:
    Array3i(int d1, int d2, int d3);
    Array3i(std::string name, int d1, int d2, int d3);
    Array3i();   // default constructer
    ~Array3i();  // destructer

    Array3i* generate(int d1, int d2, int d3);
    Array3i* generate(std::string name, int d1, int d2, int d3);
    void init(std::string name, int d1, int d2, int d3);
    void init(int d1, int d2, int d3);
    void memalloc();
    void zero();
    void print();
    void release();
    void set(int h, int i, int j, int value);
    int get(int h, int i, int j);
};
}
}  // End Namespaces
#endif  // _psi_src_bin_occ_arrays_h_
