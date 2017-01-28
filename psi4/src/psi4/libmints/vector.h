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

#ifndef _psi_src_lib_libmints_vector_h
#define _psi_src_lib_libmints_vector_h

#include "psi4/libmints/dimension.h"

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iterator>
#include <memory>

#include "psi4/pybind11.h"

namespace psi {

class Matrix;

class VectorIterator;

/*! \ingroup MINTS */
class Vector
{
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
    Vector(int nirrep, int *dimpi);

    /// Constructor, convenience for 1 irrep
    Vector(int dim);

    /// Constructor, allocates memory
    /// (this should be deprecated in favor of the Dimension-based version)
    Vector(const std::string &name, int nirrep, int *dimpi);

    /// Constructor, convenience for 1 irrep
    Vector(const std::string &name, int dim);

    /// Constructor, takes Dimension object
    Vector(const Dimension &dimpi);

    /// Constructor, takes Dimension object
    Vector(const std::string &name, const Dimension &dimpi);

    /// Destructor, frees memory
    virtual ~Vector();

    /**
     * Convenient creation function return SharedMatrix
     */
    static std::shared_ptr<Vector> create(const std::string &name,
                                          const Dimension &dim);

    void init(int nirrep, int *dimpi);

    void init(int nirrep, const int *dimpi, const std::string &name = "");

    void init(const Dimension &v);

    Vector *clone();

    /// Sets the vector_ to the data in vec
    void set(double *vec);

    /// Returns a pointer to irrep h
    double *pointer(int h = 0) { return vector_[h]; }

    const double *pointer(int h = 0) const { return vector_[h]; }

    /// Returns a single element value
    double get(int h, int m) { return vector_[h][m]; }

    /// Sets a single element value
    void set(int h, int m, double val) { vector_[h][m] = val; }

    /// Returns a single element value
    double get(int m) { return vector_[0][m]; }

    /// Sets a single element value
    void set(int m, double val) { vector_[0][m] = val; }

    void add(int m, double val) { vector_[0][m] += val; }

    void add(int h, int m, double val) { vector_[h][m] += val; }

    void add(const std::vector<double> &rhs);

    /// Adds other vector to this
    void add(const std::shared_ptr<Vector> &other);
    void add(const Vector &other);

    /// Subtracts other vector from this
    void subtract(const std::shared_ptr<Vector> &other);
    void subtract(const Vector &other);

    void axpy(double scale, const std::shared_ptr<Vector> &other);
    void axpy(double scale, const Vector &other);

    /// Zeros the vector out
    void zero();


    double &operator()(int i) { return vector_[0][i]; }

    const double &operator()(int i) const { return vector_[0][i]; }

    double &operator[](int i) { return vector_[0][i]; }

    const double &operator[](int i) const { return vector_[0][i]; }

    double pyget(const py::tuple &key);

    void pyset(const py::tuple &key, double value);

    double pyget(int key);

    void pyset(int key, double value);

    /// Returns a copy of the vector_
    double *to_block_vector();

    /// Returns the dimension per irrep h
    int dim(int h = 0) const { return dimpi_[h]; }

    /// Returns the dimension array
    const Dimension& dimpi() const { return dimpi_; }

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
    void print(std::string outfile = "outfile", const char *extra = NULL) const;

    /// Copies rhs to this
    void copy(const Vector *rhs);

    /// Copies rhs to this
    void copy(const Vector &rhs);

    /// General matrix vector multiplication
    void gemv(bool transa, double alpha, Matrix *A, Vector *X, double beta);

    /// Vector dot product
    double vector_dot(const std::shared_ptr<Vector> &other);
    double vector_dot(const Vector &other);
    double dot(Vector *X);

    /// Vector norm
    double norm();
    double sum_of_squares();
    double rms();

    /// Scale the elements of the vector
    void scale(const double &sc);

    // Serializable pure virtual functions:
    void send();

    void recv();

    void bcast(int broadcaster);

    /**
     * Performs element-by-element sum of all data from all nodes.
     */
    void sum();

    typedef std::vector<double>::iterator iterator;
    typedef std::vector<double>::const_iterator const_iterator;

    /// @{
    /** Returns the starting iterator for the entire v_. */
    iterator begin() { return v_.begin(); }

    const_iterator begin() const { return v_.begin(); }
    /// @}

    /// @{
    /** Returns the ending iterator for the entire v_. */
    iterator end() { return v_.end(); }

    const_iterator end() const { return v_.end(); }
    /// @}

    /// @{
    /** Returns the starting iterator for irrep h. */
    iterator begin_irrep(int h) {
        iterator it = v_.begin();
        for (int g = 0; g < h; ++g) it += dimpi_[h];
        return it;
    }
    // The following won't compile with clang++ and c++11
    // const_iterator begin_irrep(int h) const
    //  { return const_iterator(vector_[h]); }
    /// @}

    /// @{
    /** Returns the starting iterator for irrep h. */
    iterator end_irrep(int h) {
        iterator it = v_.begin();
        for (int g = 0; g <= h; ++g) it += dimpi_[h];
        return it;
    }
    // The following won't compile with clang++ and c++11
    // const_iterator end_irrep(int h) const
    //    { return const_iterator(vector_[h]) + dimpi_[h]; }
    /// @}

    /**
    * Adds accessability to the matrix shape for numpy
    */
    void set_numpy_shape(std::vector<int> shape) { numpy_shape_ = shape; }
    std::vector<int> numpy_shape() { return numpy_shape_; }
    std::vector<py::buffer_info> array_interface();

    friend class Matrix;
};

/*! \ingroup MINTS */
class IntVector
{
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
    int *pointer(int h = 0)
    {
        return vector_[h];
    }

    /// Returns a single element value
    int get(int h, int m)
    {
        return vector_[h][m];
    }

    /// Sets a single element value
    void set(int h, int m, int val)
    {
        vector_[h][m] = val;
    }

    /// Returns a copy of the vector_
    int *to_block_vector();

    /// Returns the dimension per irrep h
    int dim(int h = 0) const
    {
        return dimpi_[h];
    }

    /// Returns the dimension array
    int *dimpi() const
    {
        return dimpi_;
    }

    /// Returns the number of irreps
    int nirrep() const
    {
        return nirrep_;
    }

    /**
     * Sets the name of the vector, used in print(...)
     *
     * @param name New name to use.
     */
    void set_name(const std::string &name)
    {
        name_ = name;
    }

    /**
     * Gets the name of the matrix.
     */
    std::string name() const
    {
        return name_;
    }

    /// Python compatible printer
    void print_out()
    { print("outfile"); }

    /**
     * Print the matrix using print_mat
     *
     * @param outfile File point to use, defaults to Psi4's outfile.
     * @param extra When printing the name of the 'extra' will be printing after the name.
     */
    void print(std::string outfile = "outfile", const char *extra = NULL) const;

    /// Copies rhs to this
    void copy(const IntVector *rhs);

    /// Copies rhs to this
    void copy(const IntVector &rhs);

    friend class VectorIterator;
};

//class VectorIterator : public std::iterator<std::forward_iterator_tag, double>
//{
//    pointer v_;

//public:
//    VectorIterator(pointer v) : v_(v) {}

//    VectorIterator(const VectorIterator& mit) : v_(mit.v_) {}

//    VectorIterator& operator++() {v_++; return *this;}  // prefix (++a)
//    VectorIterator operator++(int) { VectorIterator tmp(*this); operator++(); return tmp; } // suffix (a++)

//    bool operator==(const VectorIterator& rhs) { return v_ == rhs.v_; }
//    bool operator!=(const VectorIterator& rhs) { return v_ != rhs.v_; }

//    reference operator*() { return *v_; }
//};

typedef std::shared_ptr <Vector> SharedVector;
typedef std::shared_ptr <IntVector> SharedIntVector;

}

#endif
