/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include <algorithm>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"
#include "dimension.h"

namespace psi {

class PSIO;
class Matrix;

class Vector;
using SharedVector = std::shared_ptr<Vector>;
class IntVector;

namespace occwave {
class Array1d;
}

/**
 * Set a vector block
 *
 * @param slice Vector slice
 * @param block the Vector object block to set
 */
template <class T>
T get_block(const T& self, const Slice& slice) {
    // check if slice is within bounds
    for (int h = 0; h < self.nirrep(); h++) {
        if (slice.end()[h] > self.dimpi()[h]) {
            std::string msg =
                "Invalid call to Vector::get_block(): Slice is out of bounds. Irrep = " + std::to_string(h);
            throw PSIEXCEPTION(msg);
        }
    }
    const Dimension &slice_begin = slice.begin();
    auto slice_dim = slice.end() - slice.begin();
    auto block = T("Block", slice_dim);
    for (int h = 0; h < self.nirrep(); h++) {
        int max_p = slice_dim[h];
        // TODO: Copy contiguous memory instead.
        for (int p = 0; p < max_p; p++) {
            double value = self.get(h, p + slice_begin[h]);
            block.set(h, p, value);
        }
    }
    return block;
}

template <class T>
class IrreppedVector {
   protected:
    /// Actual data, of size dimpi_.sum()
    std::vector<T> v_;
    /// Pointer offsets into v_, of size dimpi_.n()
    std::vector<T*> vector_;
    /// Dimensions per irrep
    Dimension dimpi_;
    /// Name
    std::string name_;

    /// Allocates vector_
    void alloc() {
        if (vector_.size()) release();

        int total = dimpi_.sum();
        v_.resize(total);

        release();

        assign_pointer_offsets();
    };

    /// Releases vector_
    void release() {
        std::fill(vector_.begin(), vector_.end(), (T *) 0);
        std::fill(v_.begin(), v_.end(), 0);
    };

    /// Assign pointer offsets in vector_ from v_.
    void assign_pointer_offsets() {
        // Resize just to be sure it's the correct size
        vector_.resize(dimpi_.n(), 0);

        size_t offset = 0;
        for (int h = 0; h < nirrep(); ++h) {
            if (dimpi_[h])
                vector_[h] = &(v_[0]) + offset;
            else
                vector_[h] = nullptr;
            offset += dimpi_[h];
        }
    };

   public:
    explicit IrreppedVector(int dim) : dimpi_(1) {
        dimpi_[0] = dim;
        alloc();
    };

    explicit IrreppedVector(const std::string & name, int dim) : dimpi_(1) {
        dimpi_[0] = dim;
        alloc();
        name_ = name;
    };

    explicit IrreppedVector(const Dimension &dimpi) {
        dimpi_ = dimpi;
        alloc();
        name_ = dimpi.name();
    };

    explicit IrreppedVector(const std::string &name, const Dimension &dimpi) {
        dimpi_ = dimpi;
        alloc();
        name_ = name;
    };

    /// Copy constructor
    IrreppedVector(const IrreppedVector<T>& vector) {
        name_ = vector.name_;
        copy(vector);
    }

    ~IrreppedVector() { release(); }

    /// Re-initialize the vector. Should usually prefer creating new object.
    void init(const Dimension &v) {
        name_ = v.name();
        dimpi_ = v;
        alloc();
    };

    /// Copy IrreppedVector's data into this.
    void copy(const IrreppedVector<T>& vector) {
        dimpi_ = vector.dimpi_;
        v_ = vector.v_;
        assign_pointer_offsets();
    }

    IrreppedVector<T> clone() const { return IrreppedVector<T>(*this) ; }

    T &operator()(int i) { return vector_[0][i]; }
    const T &operator()(int i) const { return vector_[0][i]; }
    T &operator[](int i) { return vector_[0][i]; }
    const T &operator[](int i) const { return vector_[0][i]; }

    /// Returns a pointer to irrep h
    T *pointer(int h = 0) { return vector_[h]; }
    const T *pointer(int h = 0) const { return vector_[h]; }

    // Exploit contiguous memory storage to handle elements across all irreps.
    T get(int m) const {
        if (m >= dimpi_.sum()) {
            throw PSIEXCEPTION("Cannot get element " + std::to_string(m) +  " since there are only " + std::to_string(dimpi_.sum()) +  " elements.");
        }
        return v_[m];
    }
    T get(int h, int m) const {
        if (h >= nirrep()) {
            throw PSIEXCEPTION("Cannot get an element of irrep " + std::to_string(h) +  ", since there are only " + std::to_string(nirrep()) +  " irreps.");
        }
        if (m >= dim(h)) {
            throw PSIEXCEPTION("Cannot get element " + std::to_string(m) +  " of irrep " + std::to_string(h) +  " which only has " + std::to_string(dim(h)) +  " elements.");
        }
        return vector_[h][m];
    }

    // Exploit contiguous memory storage to handle elements across all irreps.
    void set(int m, T val) {
        if (m >= dimpi_.sum()) {
            throw PSIEXCEPTION("Cannot set element " + std::to_string(m) +  " since there are only " + std::to_string(dimpi_.sum()) +  " elements.");
        }
        v_[m] = val;
    }
    void set(int h, int m, T val) {
        if (h >= nirrep()) {
            throw PSIEXCEPTION("Cannot set an element of irrep " + std::to_string(h) +  ", since there are only " + std::to_string(nirrep()) +  " irreps.");
        }
        if (m >= dim(h)) {
            throw PSIEXCEPTION("Cannot set element " + std::to_string(m) +  " of irrep " + std::to_string(h) +  " which only has " + std::to_string(dim(h)) +  " elements.");
        }
        vector_[h][m] = val;
    }

    // Exploit contiguous memory storage to handle elements across all irreps.
    void add(int m, T val) {
        if (m >= dimpi_.sum()) {
            throw PSIEXCEPTION("Cannot add to element " + std::to_string(m) +  " since there are only " + std::to_string(dimpi_.sum()) +  " elements.");
        }
        v_[m] += val;
    }
    void add(int h, int m, T val) {
        if (h >= nirrep()) {
            throw PSIEXCEPTION("Cannot add to an element of irrep " + std::to_string(h) +  ", since there are only " + std::to_string(nirrep()) +  " irreps.");
        }
        if (m >= dim(h)) {
            throw PSIEXCEPTION("Cannot add to element " + std::to_string(m) +  " of irrep " + std::to_string(h) +  " which only has " + std::to_string(dim(h)) +  " elements.");
        }
        vector_[h][m] += val;
    }

    void zero() { std::fill(v_.begin(), v_.end(), 0); }

    int dim(int h = 0) const { return dimpi_[h]; }
    const Dimension &dimpi() const { return dimpi_; }
    int nirrep() const { return dimpi_.n(); }

    // Getter and setter for name, used in print functions
    std::string name() const { return name_; }
    void set_name(const std::string &name) { name_ = name; }

    void sort(std::function<bool(int, T, T)> func) {
        sort(func, Slice(Dimension(nirrep()), dimpi_));
    }
    void sort(std::function<bool(int, T, T)> func, Slice slice) {
        for (int h = 0; h < nirrep(); h++) {
            std::stable_sort(vector_[h] + slice.begin()[h], vector_[h] + slice.end()[h], std::bind(func, h, std::placeholders::_1, std::placeholders::_2));
        }
    }

    void sort(const IrreppedVector<int>& idxs) {
        auto orig = clone();
        if (dimpi_ != idxs.dimpi()) {
            throw PSIEXCEPTION("Indexing vector and vector to sort must have the same dimension.");
        }
        // WARNING! Function also requires each irrepped block to be a permutation of 0, 1, 2...
        for (int h = 0; h < nirrep(); h++) {
            for (int i = 0; i < dimpi_[h]; i++) {
                set(h, i, orig.get(h, idxs.get(h, i)));
            }
        }
    }

    void print(const std::string& out, const std::string& fmt) const {
        const auto printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
        printer->Printf("\n # %s #\n", name_.c_str());
        const std::string temp_str = "   %4d: " + fmt + "\n"; 
        for (int h = 0; h < nirrep(); ++h) {
            printer->Printf(" Irrep: %d\n", h + 1);
            for (int i = 0; i < dimpi_[h]; ++i) printer->Printf(temp_str.c_str(), i + 1, vector_[h][i]);
            printer->Printf("\n");
        }
    }

    IrreppedVector<T> get_block(const Slice &slice) const { return psi::get_block(*this, slice); };
    void set_block(const Slice &slice, const IrreppedVector<T>& block) {
        // check if slice is within bounds
        for (int h = 0; h < nirrep(); h++) {
            if (slice.end()[h] > dimpi_[h]) {
                std::string msg =
                    "Invalid call to get_block(): Slice is out of bounds. Irrep = " + std::to_string(h);
                throw PSIEXCEPTION(msg);
            }
        }
        const Dimension &slice_begin = slice.begin();
        Dimension slice_dim = slice.end() - slice.begin();
        for (int h = 0; h < nirrep(); h++) {
            int max_p = slice_dim[h];
            for (int p = 0; p < max_p; p++) {
                double value = block.get(h, p);
                set(h, p + slice_begin[h], value);
            }
        }
    }

};

/*! \ingroup MINTS */
class PSI_API Vector final : public IrreppedVector<double> {
   protected:

    /// Numpy Shape
    std::vector<int> numpy_shape_;

   public:
    explicit Vector(int dim) : IrreppedVector<double>(dim) {};
    explicit Vector(const std::string& name, int dim) : IrreppedVector<double>(name, dim) {};
    explicit Vector(const Dimension &dimpi) : IrreppedVector<double>(dimpi) {};
    explicit Vector(const std::string& name, const Dimension &dimpi) : IrreppedVector<double>(name, dimpi) {};
    Vector(const Vector& vector) : IrreppedVector<double>(vector) {};

    // Convert occ's Array1d into Vector.
    // Defined in occ/arrays.cc. Remove when no longer needed.
    explicit Vector(const Dimension &dimpi, const occwave::Array1d &array);

    Vector clone() const { return Vector(*this); }
    Vector get_block(const Slice &slice) const { return psi::get_block(*this, slice); };

    /// Adds other elt/vector to this
    void add(int m, double val) { IrreppedVector<double>::add(m, val); }
    void add(int h, int m, double val) { IrreppedVector<double>::add(h, m, val); }
    void add(const Vector &other);

    /// Subtracts other vector from this
    void subtract(const Vector &other);

    void axpy(double scale, const Vector &other);
    void axpby(double alpha, double beta, const Vector &other);

    void print(std::string outfile = "outfile") const { IrreppedVector<double>::print(outfile, "%20.15f"); };

    /**
     * General matrix vector multiplication into this, alpha * AX + beta Y -> Y
     *
     * @ transa Do transpose A?
     * @ alpha Scaling factor
     * @ A Matrix to multiply by.
     * @ X Vector to multiply by.
     * @ beta Scaling factor for current input.
     */
    void gemv(bool transa, double alpha, const Matrix& A, const Vector& X, double beta);

    /// Vector dot product
    double vector_dot(const Vector &other) const;

    /// Vector norm
    double norm() const;
    double sum_of_squares() const;
    double rms() const;

    /// Scale the elements of the vector
    void scale(double sc);

    /// Save the Vector to disk
    void save(psi::PSIO* const psio, size_t fileno) const;
    /// Load a Vector from disk
    void load(psi::PSIO* const psio, size_t fileno);

    /**
     * Adds accessability to the matrix shape for numpy
     */
    void set_numpy_shape(std::vector<int> shape) { numpy_shape_ = shape; }
    std::vector<int> numpy_shape() const { return numpy_shape_; }

};

/*! \ingroup MINTS */
class PSI_API IntVector : public IrreppedVector<int> {
   public:
    explicit IntVector(int dim) : IrreppedVector<int>(dim) {};
    explicit IntVector(const std::string& name, int dim) : IrreppedVector<int>(name, dim) {};
    explicit IntVector(const Dimension &dimpi) : IrreppedVector<int>(dimpi) {};
    explicit IntVector(const std::string& name, const Dimension &dimpi) : IrreppedVector<int>(name, dimpi) {};
    IntVector(const IntVector& vector) : IrreppedVector<int>(vector) {};

    IntVector clone() const { return IntVector(*this) ; }
    IntVector get_block(const Slice &slice) const { return psi::get_block(*this, slice); };

    void print(std::string outfile = "outfile") const { IrreppedVector<int>::print(outfile, "%10d"); };

    static IntVector iota(const Dimension &dim);
};

}  // namespace psi
