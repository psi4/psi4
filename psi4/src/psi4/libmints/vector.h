/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#pragma once

#include <vector>
#include <memory>

#include "dimension.h"

namespace psi {

class Matrix;

class Vector;
using SharedVector = std::shared_ptr<Vector>;

/*! \ingroup MINTS */
class PSI_API Vector final {
   protected:
    /// Actual data, of size dimpi_.sum()
    std::vector<double> v_;
    /// Pointer offsets into v_, of size dimpi_.n()
    std::vector<double *> vector_;
    /// Number of irreps
    int nirrep_;
    /// Dimensions per irrep
    Dimension dimpi_;
    /// Name
    std::string name_;

    /// Allocates vector_
    void alloc();

    /// Releases vector_
    void release();

    /// Copies data to vector_
    void copy_from(const Vector &other);

    /// Assign pointer offsets in vector_ from v_.
    void assign_pointer_offsets();

    /// Numpy Shape
    std::vector<int> numpy_shape_;

   public:
    /// Default constructor, zeros everything out
    Vector();

    /// Copy constructor
    Vector(const Vector &copy);

    /// Constructor, allocates memory
    /// (this should be deprecated in favor of the Dimension-based version)
    explicit Vector(int nirrep, int *dimpi);

    /// Constructor, convenience for 1 irrep
    explicit Vector(int dim);

    /// Constructor, allocates memory
    /// (this should be deprecated in favor of the Dimension-based version)
    explicit Vector(const std::string &name, int nirrep, int *dimpi);

    /// Constructor, convenience for 1 irrep
    explicit Vector(const std::string &name, int dim);

    /// Constructor, takes Dimension object
    explicit Vector(const Dimension &dimpi);

    /// Constructor, takes Dimension object
    explicit Vector(const std::string &name, const Dimension &dimpi);

    /// Destructor, frees memory
    ~Vector();

    void init(int nirrep, int *dimpi);

    void init(int nirrep, const int *dimpi, const std::string &name = "");

    void init(const Dimension &v);

    Vector *clone();

    /// Returns a pointer to irrep h
    double *pointer(int h = 0) { return vector_[h]; }

    const double *pointer(int h = 0) const { return vector_[h]; }

    /// Returns a single element value
    double get(int m) const { return vector_[0][m]; }

    /// Sets a single element value
    void set(int m, double val) { vector_[0][m] = val; }

    /// Returns a single element value
    double get(int h, int m) const { return vector_[h][m]; }

    /// Sets a single element value
    void set(int h, int m, double val) { vector_[h][m] = val; }

    void add(int m, double val) { vector_[0][m] += val; }

    void add(int h, int m, double val) { vector_[h][m] += val; }

    /// Adds other vector to this
    void add(const SharedVector &other);
    void add(const Vector &other);

    /// Subtracts other vector from this
    void subtract(const SharedVector &other);
    void subtract(const Vector &other);

    void axpy(double scale, const SharedVector &other);
    void axpy(double scale, const Vector &other);

    /// Zeros the vector out
    void zero();

    /**
     * Get a vector block
     *
     * @param slice Vector slice
     * @return SharedVector object
     */
    SharedVector get_block(const Slice &slice);

    /**
     * Set a vector block
     *
     * @param slice Vector slice
     * @param block the SharedVector object block to set
     */
    void set_block(const Slice &slice, SharedVector block);

    double &operator()(int i) { return vector_[0][i]; }

    const double &operator()(int i) const { return vector_[0][i]; }

    double &operator[](int i) { return vector_[0][i]; }

    const double &operator[](int i) const { return vector_[0][i]; }

    /// Returns the dimension per irrep h
    int dim(int h = 0) const { return dimpi_[h]; }

    /// Returns the dimension array
    const Dimension &dimpi() const { return dimpi_; }

    /// Returns the number of irreps
    int nirrep() const { return nirrep_; }

    /**
     * Sets the name of the vector, used in print(...)
     *
     * @param name New name to use.
     */
    void set_name(const std::string &name) { name_ = name; }

    /**
     * Gets the name of the matrix.
     */
    std::string name() const { return name_; }

    /// Python compatible printer
    void print_out() { print("outfile"); }

    /**
     * Print the matrix using print_mat
     *
     * @param outfile File point to use, defaults to Psi4's outfile.
     * @param extra When printing the name of the 'extra' will be printing after the name.
     */
    void print(std::string outfile = "outfile", const char *extra = nullptr) const;

    /// Copies rhs to this
    void copy(const Vector *rhs);

    /// Copies rhs to this
    void copy(const Vector &rhs);

    /// General matrix vector multiplication
    void gemv(bool transa, double alpha, Matrix *A, Vector *X, double beta);

    /// Vector dot product
    double vector_dot(const SharedVector &other);
    double vector_dot(const Vector &other);
    double dot(Vector *X);

    /// Vector norm
    double norm();
    double sum_of_squares();
    double rms();

    /// Scale the elements of the vector
    void scale(double sc);

    /**
     * Adds accessability to the matrix shape for numpy
     */
    void set_numpy_shape(std::vector<int> shape) { numpy_shape_ = shape; }
    std::vector<int> numpy_shape() { return numpy_shape_; }

    PSI_DEPRECATED(
        "Using `Vector::create` instead of `auto my_vec = std::make_shared<Vector>(name, dim);` "
        "is deprecated, and in 1.4 it will "
        "stop working")
    static SharedVector create(const std::string &name, const Dimension &dim) {
        return std::make_shared<Vector>(name, dim);
    }
};

/*! \ingroup MINTS */
class PSI_API IntVector {
   protected:
    /// IntVector data
    int **vector_;
    /// Number of irreps
    int nirrep_;
    /// Dimensions per irrep
    int *dimpi_;
    /// Name of the IntVector
    std::string name_;

    /// Allocates vector_
    void alloc();

    /// Releases vector_
    void release();

    /// Copies data to vector_
    void copy_from(int **);

   public:
    /// Default constructor, zeros everything out
    IntVector();

    /// Copy constructor
    IntVector(const IntVector &copy);

    /// Constructor, allocates memory
    IntVector(int nirrep, int *dimpi);

    /// Constructor, convenience for 1 irrep
    IntVector(int dim);

    /// Constructor, allocates memory
    IntVector(const std::string &name, int nirrep, int *dimpi);

    /// Constructor, convenience for 1 irrep
    IntVector(const std::string &name, int dim);

    /// Destructor, frees memory
    virtual ~IntVector();

    void init(int nirrep, int *dimpi);

    /// Sets the vector_ to the data in vec
    void set(int *vec);

    /// Returns a pointer to irrep h
    int *pointer(int h = 0) { return vector_[h]; }

    /// Returns a single element value
    int get(int h, int m) { return vector_[h][m]; }

    /// Sets a single element value
    void set(int h, int m, int val) { vector_[h][m] = val; }

    /// Returns a copy of the vector_
    int *to_block_vector();

    /// Returns the dimension per irrep h
    int dim(int h = 0) const { return dimpi_[h]; }

    /// Returns the dimension array
    int *dimpi() const { return dimpi_; }

    /// Returns the number of irreps
    int nirrep() const { return nirrep_; }

    /**
     * Sets the name of the vector, used in print(...)
     *
     * @param name New name to use.
     */
    void set_name(const std::string &name) { name_ = name; }

    /**
     * Gets the name of the matrix.
     */
    std::string name() const { return name_; }

    /// Python compatible printer
    void print_out() { print("outfile"); }

    /**
     * Print the matrix using print_mat
     *
     * @param outfile File point to use, defaults to Psi4's outfile.
     * @param extra When printing the name of the 'extra' will be printing after the name.
     */
    void print(std::string outfile = "outfile", const char *extra = nullptr) const;

    /// Copies rhs to this
    void copy(const IntVector *rhs);

    /// Copies rhs to this
    void copy(const IntVector &rhs);
};

}  // namespace psi
