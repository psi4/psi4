/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

#ifndef _PSI_SRC_LIB_LIBDIIS_DIISENTRY_H_
#define _PSI_SRC_LIB_LIBDIIS_DIISENTRY_H_

#include <map>
#include <memory>
#include <string>

namespace psi {

class PSIO;

/**
 * @brief The DIISManager class is used to manage DIIS quantities and
 * their corresponding error vectors.
 */

class DIISEntry {
   public:
    /**
     * @brief The type of quantity to add to either an error vector or vector.
     *
     * DPDBuf4  - A pointer to a four-index DPD quantity, assumed to be initialized
     *            via dpd_buf4_init() already.
     * DPDFile2 - A pointer to a two-index DPD quality, assumed to be initialized
     *            via dpd_file2_init() already.
     * Matrix   - A pointer to a Psi Matrix object.
     * Vector   - A pointer to a Psi Vector object.
     * Pointer  - The quantity will be taken from contiguous memory pointed to by a
     *            double* pointer.  In this case, the next parameter passed to the
     *            set_size.. routines should be the number of elements to read.
     * Psio     - The PSIO object to use for I/O
     */
    enum class InputType { DPDBuf4, DPDFile2, Matrix, Vector, Pointer };
    DIISEntry(std::string label, int ID, int count, std::vector<double> vector,
              std::vector<double> errorVector, std::shared_ptr<PSIO> psio);
    /// Whether the dot product of this entry's and the nth entry's error vector is known
    bool dot_is_known_with(int n) { return _knownDotProducts[n]; }
    /// The dot product of this entry's and the nth entry's error vectors
    double dot_with(int n) { return _dotProducts[n]; }
    /// The RMS error of this entry
    double rmsError() { return _rmsError; }
    /// The absolute number of this entry
    int orderAdded() { return _orderAdded; }
    /// Sets the dot product with vector n to val
    void set_dot_with(int n, double val) {
        _knownDotProducts[n] = true;
        _dotProducts[n] = val;
    }
    /// Marks the dot product with vector n as invalid
    void invalidate_dot(int n) { _knownDotProducts[n] = false; }
    /// Set the vector
    void set_vector(std::vector<double> vec) { _vector = std::move(vec); }
    /// Put this vector entry on disk and free the memory
    void dump_vector_to_disk();
    /// Allocate vector memory and read from disk
    void read_vector_from_disk();
    /// Put this error vector entry on disk and free the memory
    void dump_error_vector_to_disk();
    /// Allocate error vector memory and read from disk
    void read_error_vector_from_disk();
    /// Free vector memory
    void free_vector_memory();
    /// Free error vector memory
    void free_error_vector_memory();
    /// Returns the error vector
    const double *errorVector() {
        read_error_vector_from_disk();
        return _errorVector.data();
    }
    /// Returns the vector
    const double *vector() {
        read_vector_from_disk();
        return _vector.data();
    }
    /// Open the psi file, if needed.
    void open_psi_file();
    /// Close the psi file, if needed.
    void close_psi_file();

   protected:
    /// The list of which dot products, with other DIISEntries, are known
    std::map<int, bool> _knownDotProducts;
    /// The list of known dot products with other DIISEntries
    std::map<int, double> _dotProducts;
    /// The length of the error vector
    int _errorVectorSize;
    /// The length of the vector
    int _vectorSize;
    /// The absolute number of this entry
    int _orderAdded;
    /// The number of this entry in the current subspace
    int _ID;
    /// The RMS error for this entry
    double _rmsError;
    /// The error vector
    std::vector<double> _errorVector;
    /// The error vector
    std::vector<double> _vector;
    /// The label used for disk storage
    std::string _label;
    /// PSIO object
    std::shared_ptr<PSIO> _psio;
};

}  // namespace psi

#endif  // Header guard
