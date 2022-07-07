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

#pragma once

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
class IrrepedVector {
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
    explicit IrrepedVector(int dim) : dimpi_(1) {
        dimpi_[0] = dim;
        alloc();
    };

    explicit IrrepedVector(const std::string & name, int dim) : dimpi_(1) {
        dimpi_[0] = dim;
        alloc();
        name_ = name;
    };

    explicit IrrepedVector(const Dimension &dimpi) {
        dimpi_ = dimpi;
        alloc();
        name_ = dimpi.name();
    };

    explicit IrrepedVector(const std::string &name, const Dimension &dimpi) {
        dimpi_ = dimpi;
        alloc();
        name_ = name;
    };

    /// Copy constructor
    explicit IrrepedVector(const IrrepedVector<T>& vector) {
        name_ = vector.name_;
        copy(vector);
    }

    ~IrrepedVector() { release(); }

    /// Re-initialize the vector. Should usually prefer creating new object.
    void init(const Dimension &v) {
        name_ = v.name();
        dimpi_ = v;
        alloc();
    };

    /// Copy IrrepedVector's data into this.
    void copy(const IrrepedVector<T>& vector) {
        dimpi_ = vector.dimpi_;
        v_ = vector.v_;
        assign_pointer_offsets();
    }

    IrrepedVector<T> clone() { return IrrepedVector<T>(*this) ; }

    T &operator()(int i) { return vector_[0][i]; }
    const T &operator()(int i) const { return vector_[0][i]; }
    T &operator[](int i) { return vector_[0][i]; }
    const T &operator[](int i) const { return vector_[0][i]; }

    /// Returns a pointer to irrep h
    T *pointer(int h = 0) { return vector_[h]; }
    const T *pointer(int h = 0) const { return vector_[h]; }

    T get(int m) const { return vector_[0][m]; }
    T get(int h, int m) const { return vector_[h][m]; }

    void set(int m, T val) { vector_[0][m] = val; }
    void set(int h, int m, T val) { vector_[h][m] = val; }

    void add(int m, T val) { vector_[0][m] += val; }
    void add(int h, int m, T val) { vector_[h][m] += val; }

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

    IrrepedVector<T> get_block(const Slice &slice) const { return psi::get_block(*this, slice); };
    void set_block(const Slice &slice, const IrrepedVector<T>& block) {
        // check if slice is within bounds
        for (int h = 0; h < nirrep(); h++) {
            if (slice.end()[h] > dimpi_[h]) {
                std::string msg =
                    "Invalid call to Vector::set_block(): Slice is out of bounds. Irrep = " + std::to_string(h);
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
class PSI_API Vector final : public IrrepedVector<double> {
   protected:

    /// Numpy Shape
    std::vector<int> numpy_shape_;

   public:
    explicit Vector(int dim) : IrrepedVector<double>(dim) {}; 
    explicit Vector(const std::string& name, int dim) : IrrepedVector<double>(name, dim) {}; 
    explicit Vector(const Dimension &dimpi) : IrrepedVector<double>(dimpi) {}; 
    explicit Vector(const std::string& name, const Dimension &dimpi) : IrrepedVector<double>(name, dimpi) {}; 
    Vector(const Vector& vector) : IrrepedVector<double>(vector) {};

    // Convert occ's Array1d into Vector.
    // Defined in occ/arrays.cc. Remove when no longer needed.
    explicit Vector(const Dimension &dimpi, const occwave::Array1d &array);

    Vector clone() { return Vector(*this); }
    Vector get_block(const Slice &slice) const { return psi::get_block(*this, slice); };

    /// Adds other elt/vector to this
    void add(int m, double val) { IrrepedVector<double>::add(m, val); }
    void add(int h, int m, double val) { IrrepedVector<double>::add(h, m, val); }
    void add(const Vector &other);

    /// Subtracts other vector from this
    void subtract(const Vector &other);

    void axpy(double scale, const Vector &other);

    void print(std::string outfile = "outfile") const { IrrepedVector<double>::print(outfile, "%20.15f"); };

    /// Sort the vector according to the re-indexing vector.
    void sort(const IntVector& idxs);

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
class PSI_API IntVector : public IrrepedVector<int> {
   public:
    explicit IntVector(int dim) : IrrepedVector<int>(dim) {}; 
    explicit IntVector(const std::string& name, int dim) : IrrepedVector<int>(name, dim) {};
    explicit IntVector(const Dimension &dimpi) : IrrepedVector<int>(dimpi) {}; 
    explicit IntVector(const std::string& name, const Dimension &dimpi) : IrrepedVector<int>(name, dimpi) {};
    IntVector(const IntVector& vector) : IrrepedVector<int>(vector) {};

    IntVector clone() { return IntVector(*this) ; }

    void print(std::string outfile = "outfile") const { IrrepedVector<int>::print(outfile, "%10d"); };

    static IntVector iota(const Dimension &dim);
};

}  // namespace psi
