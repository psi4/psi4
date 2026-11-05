/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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
#include <concepts>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
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
        if (static_cast<size_t>(slice.end()[h]) > self.block_size(h)) {
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
    /// Contiguous storage for all irreps.
    std::vector<T> v_;
    /// One pointer per irrep (null for an empty block).
    std::vector<T*> vector_;
    /// Legacy int dimensions. Empty for a flat vector whose length only fits in size_t.
    Dimension dimpi_;
    /// Name
    std::string name_;

    template <std::integral I>
    static size_t checked_length(I length) {
        if constexpr (std::is_signed_v<I>) {
            if (length < 0) throw PSIEXCEPTION("Vector length cannot be negative.");
        }
        if constexpr (std::numeric_limits<I>::digits > std::numeric_limits<size_t>::digits) {
            if (length > static_cast<I>(std::numeric_limits<size_t>::max()))
                throw PSIEXCEPTION("Vector length does not fit in size_t.");
        }
        return static_cast<size_t>(length);
    }

    void alloc(const Dimension& dims) {
        if (dims.n() > std::numeric_limits<int>::max()) throw PSIEXCEPTION("Too many Vector irreps.");
        size_t total = 0;
        for (int length : dims.blocks()) {
            const auto n = checked_length(length);
            if (n > v_.max_size() - total) throw PSIEXCEPTION("Vector exceeds the maximum allocation size.");
            total += n;
        }
        v_.assign(total, T{});
        dimpi_ = dims;
        assign_pointer_offsets();
    }

    /// Assign pointer offsets in vector_ from v_.
    void assign_pointer_offsets() {
        // A nonempty flat vector with no representable Dimension still has exactly one irrep.
        if (!dimpi_.n() && !v_.empty()) {
            vector_.assign(1, v_.data());
            return;
        }
        vector_.resize(dimpi_.n());

        size_t offset = 0;
        for (int h = 0; h < nirrep(); ++h) {
            vector_[h] = dimpi_[h] ? v_.data() + offset : nullptr;
            offset += dimpi_[h];
        }
    };

    IrreppedVector(const std::string& name, size_t length)
        : v_(length), dimpi_(length <= std::numeric_limits<int>::max() ? 1 : 0), name_(name) {
        if (dimpi_.n()) dimpi_[0] = static_cast<int>(length);
        assign_pointer_offsets();
    }

   public:
    explicit IrreppedVector(int dim) : IrreppedVector("", checked_length(dim)) {}

    explicit IrreppedVector(const std::string& name, int dim) : IrreppedVector(name, checked_length(dim)) {}

    explicit IrreppedVector(const Dimension &dimpi) {
        alloc(dimpi);
        name_ = dimpi.name();
    };

    explicit IrreppedVector(const std::string &name, const Dimension &dimpi) {
        alloc(dimpi);
        name_ = name;
    };

    /// Copy constructor
    IrreppedVector(const IrreppedVector<T>& vector) {
        name_ = vector.name_;
        copy(vector);
    }

    IrreppedVector& operator=(const IrreppedVector& other) {
        if (this != &other) {
            copy(other);
            name_ = other.name_;
        }
        return *this;
    }

    /// Re-initialize the vector. Should usually prefer creating new object.
    void init(const Dimension &v) {
        alloc(v);
        name_ = v.name();
    };

    /// Copy IrreppedVector's data into this.
    void copy(const IrreppedVector<T>& vector) {
        dimpi_ = vector.dimpi_;
        v_ = vector.v_;
        assign_pointer_offsets();
    }

    IrreppedVector<T> clone() const { return IrreppedVector<T>(*this) ; }

    T &operator()(size_t i) { return vector_[0][i]; }
    const T &operator()(size_t i) const { return vector_[0][i]; }
    T &operator[](size_t i) { return vector_[0][i]; }
    const T &operator[](size_t i) const { return vector_[0][i]; }

    /// Returns a pointer to irrep h
    T *pointer(int h = 0) { return vector_[h]; }
    const T *pointer(int h = 0) const { return vector_[h]; }

    // Exploit contiguous memory storage to handle elements across all irreps.
    T get(size_t m) const {
        if (m >= size()) {
            throw PSIEXCEPTION("Cannot get element " + std::to_string(m) +  " since there are only " + std::to_string(size()) +  " elements.");
        }
        return v_[m];
    }
    T get(int h, size_t m) const {
        if (m >= block_size(h)) {
            throw PSIEXCEPTION("Cannot get element " + std::to_string(m) +  " of irrep " + std::to_string(h) +  " which only has " + std::to_string(block_size(h)) +  " elements.");
        }
        return vector_[h][m];
    }

    // Exploit contiguous memory storage to handle elements across all irreps.
    void set(size_t m, T val) {
        if (m >= size()) {
            throw PSIEXCEPTION("Cannot set element " + std::to_string(m) +  " since there are only " + std::to_string(size()) +  " elements.");
        }
        v_[m] = val;
    }
    void set(int h, size_t m, T val) {
        if (m >= block_size(h)) {
            throw PSIEXCEPTION("Cannot set element " + std::to_string(m) +  " of irrep " + std::to_string(h) +  " which only has " + std::to_string(block_size(h)) +  " elements.");
        }
        vector_[h][m] = val;
    }

    // Exploit contiguous memory storage to handle elements across all irreps.
    void add(size_t m, T val) {
        if (m >= size()) {
            throw PSIEXCEPTION("Cannot add to element " + std::to_string(m) +  " since there are only " + std::to_string(size()) +  " elements.");
        }
        v_[m] += val;
    }
    void add(int h, size_t m, T val) {
        if (m >= block_size(h)) {
            throw PSIEXCEPTION("Cannot add to element " + std::to_string(m) +  " of irrep " + std::to_string(h) +  " which only has " + std::to_string(block_size(h)) +  " elements.");
        }
        vector_[h][m] += val;
    }

    void zero() { std::fill(v_.begin(), v_.end(), 0); }

    /// Total number of elements, including all irreps.
    size_t size() const { return v_.size(); }
    size_t block_size(int h = 0) const {
        if (h < 0 || h >= nirrep()) throw PSIEXCEPTION("Vector irrep index is out of range.");
        return nirrep() == 1 ? size() : static_cast<size_t>(dimpi_[h]);
    }
    /// Legacy int-sized shape accessors; never return a truncated length.
    int dim(int h = 0) const {
        auto n = block_size(h);
        if (n > std::numeric_limits<int>::max()) throw PSIEXCEPTION("Vector::dim: use block_size() for lengths above INT_MAX.");
        return static_cast<int>(n);
    }
    const Dimension &dimpi() const {
        if (vector_.size() != dimpi_.n()) throw PSIEXCEPTION("Vector length exceeds the int range of Dimension; use block_size().");
        return dimpi_;
    }
    int nirrep() const { return static_cast<int>(vector_.size()); }

    // Getter and setter for name, used in print functions
    std::string name() const { return name_; }
    void set_name(const std::string &name) { name_ = name; }

    void sort(std::function<bool(int, T, T)> func) {
        for (int h = 0; h < nirrep(); ++h) {
            const size_t n = block_size(h);
            if (n > 1)
                std::stable_sort(pointer(h), pointer(h) + n, std::bind(func, h, std::placeholders::_1, std::placeholders::_2));
        }
    }
    void sort(std::function<bool(int, T, T)> func, Slice slice) {
        for (int h = 0; h < nirrep(); h++) {
            std::stable_sort(vector_[h] + slice.begin()[h], vector_[h] + slice.end()[h], std::bind(func, h, std::placeholders::_1, std::placeholders::_2));
        }
    }

    void sort(const IrreppedVector<int>& idxs) {
        if (dimpi() != idxs.dimpi()) {
            throw PSIEXCEPTION("Indexing vector and vector to sort must have the same dimension.");
        }
        auto orig = clone();
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
        const std::string temp_str = "   %4zu: " + fmt + "\n";
        for (int h = 0; h < nirrep(); ++h) {
            printer->Printf(" Irrep: %d\n", h + 1);
            for (size_t i = 0; i < block_size(h); ++i) printer->Printf(temp_str.c_str(), i + 1, vector_[h][i]);
            printer->Printf("\n");
        }
    }

    IrreppedVector<T> get_block(const Slice &slice) const { return psi::get_block(*this, slice); };
    void set_block(const Slice &slice, const IrreppedVector<T>& block) {
        // check if slice is within bounds
        for (int h = 0; h < nirrep(); h++) {
            if (static_cast<size_t>(slice.end()[h]) > block_size(h)) {
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
    /// Construct one contiguous irrep, including lengths above INT_MAX.
    explicit Vector(size_t dim) : Vector("", dim) {}
    explicit Vector(const std::string& name, size_t dim) : IrreppedVector<double>(name, dim) {}
    // Keep calls with other integral types unambiguous between int and size_t.
    template <std::integral I>
    explicit Vector(I dim) : Vector(checked_length(dim)) {}
    template <std::integral I>
    explicit Vector(const std::string& name, I dim) : Vector(name, checked_length(dim)) {}
    explicit Vector(const Dimension &dimpi) : IrreppedVector<double>(dimpi) {};
    explicit Vector(const std::string& name, const Dimension &dimpi) : IrreppedVector<double>(name, dimpi) {};
    Vector(const Vector& vector) : IrreppedVector<double>(vector) {};

    // Convert occ's Array1d into Vector.
    // Defined in occ/arrays.cc. Remove when no longer needed.
    explicit Vector(const Dimension &dimpi, const occwave::Array1d &array);

    Vector clone() const { return Vector(*this); }
    Vector get_block(const Slice &slice) const { return psi::get_block(*this, slice); };

    /// Adds other elt/vector to this
    void add(size_t m, double val) { IrreppedVector<double>::add(m, val); }
    void add(int h, size_t m, double val) { IrreppedVector<double>::add(h, m, val); }
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
