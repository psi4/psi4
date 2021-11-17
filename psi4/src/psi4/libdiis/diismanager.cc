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

#include "diismanager.h"

#include <cmath>
#include <cstdarg>
#include <memory>

#include "psi4/psifiles.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"

using namespace psi;

namespace psi {
/**
 *
 * @param maxSubspaceSize Maximum number of vectors allowed in the subspace
 * @param label: the base part of the label used to store the vectors to disk
 * @param removalPolicy: How to decide which vectors to remove when the subspace is full
 * @param storagePolicy: How to store the DIIS vectors
 * @param psio: the PSIO object to use for I/O.  Do not specify if DPD is being used.
 */
DIISManager::DIISManager(int maxSubspaceSize, const std::string &label, RemovalPolicy removalPolicy,
                         StoragePolicy storagePolicy)
    : _maxSubspaceSize(maxSubspaceSize),
      _removalPolicy(removalPolicy),
      _storagePolicy(storagePolicy),
      _errorVectorSize(0),
      _vectorSize(0),
      _psio(_default_psio_lib_),
      _entryCount(0),
      _label(label) {}

int DIISManager::subspace_size() { return _subspace.size(); }

/**
 * Determines the size of the error vector from a list of input quantities.  This function should
 * not be called until set_error_vector_size() has been called.
 */
void DIISManager::set_vector_size(int numQuantities, ...) {
    if (_vectorSize)
        throw SanityCheckError("DIISManager: The size of the DIIS vector has already been set", __FILE__, __LINE__);
    if (_errorVectorSize == 0)
        throw SanityCheckError("DIISManager: The error vector size must be set before the vector size", __FILE__,
                               __LINE__);
    _numVectorComponents = numQuantities;
    dpdfile2 *file2;
    dpdbuf4 *buf4;
    Vector *vector;
    Matrix *matrix;
    va_list args;
    va_start(args, numQuantities);
    for (int i = 0; i < numQuantities; ++i) {
        DIISEntry::InputType type = static_cast<DIISEntry::InputType>(va_arg(args, int));
        _componentTypes.push_back(type);
        size_t size = 0;
        switch (type) {
            case DIISEntry::InputType::Pointer:
                size += va_arg(args, int);
                break;
            case DIISEntry::InputType::DPDBuf4:
                buf4 = va_arg(args, dpdbuf4 *);
                for (int h = 0; h < buf4->params->nirreps; ++h) {
                    size += static_cast<unsigned long> (buf4->params->rowtot[h]) * buf4->params->coltot[h];
                }
                break;
            case DIISEntry::InputType::DPDFile2:
                file2 = va_arg(args, dpdfile2 *);
                for (int h = 0; h < file2->params->nirreps; ++h) {
                    size += static_cast<unsigned long> (file2->params->rowtot[h]) * file2->params->coltot[h];
                }
                break;
            case DIISEntry::InputType::Matrix:
                matrix = va_arg(args, Matrix *);
                for (int h = 0; h < matrix->nirrep(); ++h) {
                    size += static_cast<unsigned long> (matrix->rowspi()[h]) * matrix->colspi()[h];
                }
                break;
            case DIISEntry::InputType::Vector:
                vector = va_arg(args, Vector *);
                size = vector->dimpi().sum();
                break;
            default:
                throw SanityCheckError("Unknown input type", __FILE__, __LINE__);
        }
        _componentSizes.push_back(size);
        _vectorSize += size;
    }
    va_end(args);
}

/**
 * Determines the size of the error vector from a list of input quantities.
 */
void DIISManager::set_error_vector_size(int numQuantities, ...) {
    if (_errorVectorSize)
        throw SanityCheckError("The size of the DIIS error vector has already been set", __FILE__, __LINE__);
    _numErrorVectorComponents = numQuantities;
    dpdfile2 *file2;
    dpdbuf4 *buf4;
    Vector *vector;
    Matrix *matrix;
    va_list args;
    va_start(args, numQuantities);
    for (int i = 0; i < numQuantities; ++i) {
        DIISEntry::InputType type = static_cast<DIISEntry::InputType>(va_arg(args, int));
        _componentTypes.push_back(type);
        size_t size = 0;
        switch (type) {
            case DIISEntry::InputType::Pointer:
                size += va_arg(args, int);
                break;
            case DIISEntry::InputType::DPDBuf4:
                buf4 = va_arg(args, dpdbuf4 *);
                for (int h = 0; h < buf4->params->nirreps; ++h) {
                    size += static_cast<unsigned long> (buf4->params->rowtot[h]) * buf4->params->coltot[h];
                }
                break;
            case DIISEntry::InputType::DPDFile2:
                file2 = va_arg(args, dpdfile2 *);
                for (int h = 0; h < file2->params->nirreps; ++h) {
                    size += static_cast<unsigned long> (file2->params->rowtot[h]) * file2->params->coltot[h];
                }
                break;
            case DIISEntry::InputType::Matrix:
                matrix = va_arg(args, Matrix *);
                for (int h = 0; h < matrix->nirrep(); ++h) {
                    size += static_cast<unsigned long> (matrix->rowspi()[h]) * matrix->colspi()[h];
                }
                break;
            case DIISEntry::InputType::Vector:
                vector = va_arg(args, Vector *);
                size = vector->dimpi().sum();
                break;
            default:
                throw SanityCheckError("Unknown input type", __FILE__, __LINE__);
        }
        _componentSizes.push_back(size);
        _errorVectorSize += size;
    }
    va_end(args);
}

/**
 * Adds a new vector and error vector to the DIIS subspace.
 * @param numQuantities - the number of components used to construct the
 *        vector and error vector.  Must match the sum of the numbers passed
 *        to the set_vector_size() and set_error_vector_size() functions.
 * The remaining parameters are then the components of the vector and then the error
 * vector, in that order, and in the same order they were passed to the
 * set_vector_size() and set_error_vector_size() functions.  N.B. Unlike the set_size
 * functions, the types of each component should not be specified here.  If the component
 * is an array, the pointer to the start of that array should be passed, in contrast to the
 * set_size functions, which takes only the size of that array.
 * @return Whether the subspace was updated
 */
bool DIISManager::add_entry(int numQuantities, ...) {
    if (!_maxSubspaceSize) return false;
    if (_componentSizes.size() != numQuantities)
        throw SanityCheckError(
            "The number of parameters passed to the set_size routines"
            " and add_entry are inconsistent",
            __FILE__, __LINE__);

    timer_on("DIISManager::add_entry");
    dpdfile2 *file2;
    dpdbuf4 *buf4;
    Vector *vector;
    Matrix *matrix;
    double *array;
    va_list args;
    va_start(args, numQuantities);
    auto errorVector = std::vector<double>(_errorVectorSize);
    auto paramVector = std::vector<double>(_vectorSize);
    double *arrayPtr = errorVector.data();

    for (int i = 0; i < numQuantities; ++i) {
        DIISEntry::InputType type = _componentTypes[i];
        // If we've filled the error vector, start filling the vector
        if (i == _numErrorVectorComponents) arrayPtr = paramVector.data();
        switch (type) {
            case DIISEntry::InputType::Pointer:
            {
                array = va_arg(args, double *);
                auto size = static_cast<int>(_componentSizes[i]);
                if (size) {
                    std::copy(array, array + size, arrayPtr);
                    arrayPtr += size;
                }
                break;
            }
            case DIISEntry::InputType::DPDBuf4:
            {
                buf4 = va_arg(args, dpdbuf4 *);
                for (int h = 0; h < buf4->params->nirreps; ++h) {
                    global_dpd_->buf4_mat_irrep_init(buf4, h);
                    global_dpd_->buf4_mat_irrep_rd(buf4, h);
                    auto size = buf4->params->rowtot[h] * buf4->params->coltot[h ^ buf4->file.my_irrep];
                    if (size) {
                        std::copy(buf4->matrix[h][0], buf4->matrix[h][0] + size, arrayPtr);
                        arrayPtr += size;
                    }
                    global_dpd_->buf4_mat_irrep_close(buf4, h);
                }
                break;
            }
            case DIISEntry::InputType::DPDFile2:
            {
                file2 = va_arg(args, dpdfile2 *);
                global_dpd_->file2_mat_init(file2);
                global_dpd_->file2_mat_rd(file2);
                for (int h = 0; h < file2->params->nirreps; ++h) {
                    auto size = file2->params->rowtot[h] * file2->params->coltot[file2->my_irrep ^ h];
                    if (size) {
                        std::copy(file2->matrix[h][0], file2->matrix[h][0] + size, arrayPtr);
                        arrayPtr += size;
                    }
                }
                global_dpd_->file2_mat_close(file2);
                break;
            }
            case DIISEntry::InputType::Matrix:
            {
                matrix = va_arg(args, Matrix *);
                for (int h = 0; h < matrix->nirrep(); ++h) {
                    auto size = matrix->rowdim(h) * matrix->coldim(h ^ matrix->symmetry());
                    if (size) {
                        std::copy(matrix->pointer(h)[0], matrix->pointer(h)[0] + size, arrayPtr);
                        arrayPtr += size;
                    }
                }
                break;
            }
            case DIISEntry::InputType::Vector:
            {
                vector = va_arg(args, Vector *);
                auto size = vector->dimpi().sum();
                if (size) {
                    std::copy(vector->pointer(), vector->pointer() + size, arrayPtr);
                    arrayPtr += size;
                }
                break;
            }
            default:
                throw SanityCheckError("Unknown input type", __FILE__, __LINE__);
        }
    }
    va_end(args);

    int entryID = get_next_entry_id();
    if (_subspace.size() < _maxSubspaceSize) {
        _subspace.emplace_back(_label, entryID, _entryCount++, std::move(errorVector), std::move(paramVector), _psio);
    } else {
        _subspace[entryID] = DIISEntry(_label, entryID, _entryCount++, std::move(errorVector), std::move(paramVector), _psio);
    }

    if (_storagePolicy == StoragePolicy::OnDisk) {
        _subspace[entryID].dump_vector_to_disk();
        _subspace[entryID].dump_error_vector_to_disk();
    }

    // Clear all inner products with this entry that may be cached
    for (int i = 0; i < _subspace.size(); ++i)
        if (i != entryID) _subspace[i].invalidate_dot(entryID);

    timer_off("DIISManager::add_entry");

    return true;
}

/**
 * Figures out the ID of the next entry to be added by determining whether an entry
 * must be removed in order to add a new one.
 */
int DIISManager::get_next_entry_id() {
    int entry = 0;
    if (_subspace.size() < _maxSubspaceSize) {
        entry = _subspace.size();
    } else {
        if (_removalPolicy == RemovalPolicy::OldestAdded) {
            int oldest = _subspace[0].orderAdded();
            for (int i = 1; i < _subspace.size(); ++i) {
                if (_subspace[i].orderAdded() < oldest) {
                    oldest = _subspace[i].orderAdded();
                    entry = i;
                }
            }
        } else if (_removalPolicy == RemovalPolicy::LargestError) {
            double largest = _subspace[0].rmsError();
            for (int i = 1; i < _subspace.size(); ++i) {
                if (_subspace[i].rmsError() > largest) {
                    largest = _subspace[i].rmsError();
                    entry = i;
                }
            }
        } else {
            throw SanityCheckError("Unknown RemovalPolicy", __FILE__, __LINE__);
        }
    }
    return entry;
}

/**
 * Performs the extapolation, based on the current subspace
 * @param numQuantitites - the number of quantities that the vector comprises.
 * Then pass these quantities in the order they were passed to the set_vector_size()
 * function.  The types of each component should not be specified here, unlike the
 * set_vector_size() function call.
 */
bool DIISManager::extrapolate(int numQuantities, ...) {
    if (!_subspace.size()) return false;

    timer_on("DIISManager::extrapolate");

    auto dimension = _subspace.size() + 1;
    auto B = Matrix("B (DIIS Connectivity Matrix", dimension, dimension);
    auto Bp = B.pointer();
    auto coefficients = std::vector<double>(dimension, 0.0);
    auto force = std::vector<double>(dimension, 0.0);

    timer_on("bMatrix setup");

    for (int i = 0; i < _subspace.size(); ++i) {
        Bp[i][_subspace.size()] = Bp[_subspace.size()][i] = 1.0;
        auto& entryI = _subspace[i];
        for (int j = 0; j < _subspace.size(); ++j) {
            auto& entryJ = _subspace[j];
            if (entryI.dot_is_known_with(j)) {
                Bp[i][j] = entryI.dot_with(j);
            } else {
                double dot = C_DDOT(_errorVectorSize, const_cast<double *>(entryI.errorVector()), 1,
                                    const_cast<double *>(entryJ.errorVector()), 1);
                Bp[i][j] = dot;
                entryI.set_dot_with(j, dot);
                entryJ.set_dot_with(i, dot);
                if (_storagePolicy == StoragePolicy::OnDisk) {
                    entryI.free_error_vector_memory();
                    entryJ.free_error_vector_memory();
                }
            }
        }
    }
    force[_subspace.size()] = 1.0;
    Bp[_subspace.size()][_subspace.size()] = 0.0;

    timer_off("bMatrix setup");
    timer_on("bMatrix pseudoinverse");

    // => Balance <= //

    auto S = std::vector<double>(dimension);

    // Trap an explicit zero
    bool is_zero = false;
    for (int i = 0; i < dimension - 1; i++) {
        if (Bp[i][i] <= 0.0) {
            is_zero = true;
        }
    }

    if (is_zero) {
        for (int i = 0; i < dimension; i++) {
            S[i] = 1.0;
        }
    } else {
        for (int i = 0; i < dimension - 1; i++) {
            S[i] = pow(Bp[i][i], -1.0 / 2.0);
        }
        S[dimension - 1] = 1.0;
    }

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            Bp[i][j] *= S[i] * S[j];
        }
    }

    // => S [S^-1 B S^-1] S \ f <= //

    B.power(-1.0, 1.0E-12);
    C_DGEMV('N', dimension, dimension, 1.0, Bp[0], dimension, force.data(), 1, 0.0, coefficients.data(), 1);
    for (int i = 0; i < dimension; i++) {
        coefficients[i] *= S[i];
    }

    timer_off("bMatrix pseudoinverse");

    timer_on("New vector");

    dpdfile2 *file2;
    dpdbuf4 *buf4;
    Vector *vector;
    Matrix *matrix;
    double *array;
    va_list args;
    int print = Process::environment.options.get_int("PRINT");
    if (print > 2) outfile->Printf("DIIS coefficients: ");
    for (int n = 0; n < _subspace.size(); ++n) {
        double coefficient = coefficients[n];
        if (print > 2) outfile->Printf(" %.3f ", coefficient);
        const double *arrayPtr = _subspace[n].vector();
        va_start(args, numQuantities);
        for (int i = 0; i < numQuantities; ++i) {
            // The indexing arrays contain the error vector, then the vector, so they
            // need to be offset by the number of components in the error vector
            int componentIndex = i + _numErrorVectorComponents;
            DIISEntry::InputType type = _componentTypes[componentIndex];
            switch (type) {
                case DIISEntry::InputType::Pointer:
                {
                    array = va_arg(args, double *);
                    auto size = _componentSizes[componentIndex];
                    if (!n) ::memset(array, 0, size * sizeof(double));
                    if (size) {
                        C_DAXPY(size, coefficient, arrayPtr, 1, array, 1);
                        arrayPtr += static_cast<int>(size);
                    }
                    break;
                }
                case DIISEntry::InputType::DPDBuf4:
                {
                    buf4 = va_arg(args, dpdbuf4 *);
                    if (!n) global_dpd_->buf4_scm(buf4, 0.0);
                    for (int h = 0; h < buf4->params->nirreps; ++h) {
                        global_dpd_->buf4_mat_irrep_init(buf4, h);
                        global_dpd_->buf4_mat_irrep_rd(buf4, h);
                        auto size = static_cast<size_t>(buf4->params->rowtot[h]) * buf4->params->coltot[h ^ buf4->file.my_irrep];
                        if (size) {
                            C_DAXPY(size, coefficient, arrayPtr, 1, buf4->matrix[h][0], 1);
                            arrayPtr += static_cast<int>(size);
                        }
                        global_dpd_->buf4_mat_irrep_wrt(buf4, h);
                        global_dpd_->buf4_mat_irrep_close(buf4, h);
                    }
                    break;
                }
                case DIISEntry::InputType::DPDFile2:
                {
                    file2 = va_arg(args, dpdfile2 *);
                    if (!n) global_dpd_->file2_scm(file2, 0.0);
                    global_dpd_->file2_mat_init(file2);
                    global_dpd_->file2_mat_rd(file2);
                    for (int h = 0; h < file2->params->nirreps; ++h) {
                        auto size = static_cast<size_t>(file2->params->rowtot[h]) * file2->params->coltot[file2->my_irrep ^ h];
                        if (size) {
                            C_DAXPY(size, coefficient, arrayPtr, 1, file2->matrix[h][0], 1);
                            arrayPtr += static_cast<int>(size);
                        }
                    }
                    global_dpd_->file2_mat_wrt(file2);
                    global_dpd_->file2_mat_close(file2);
                    break;
                }
                case DIISEntry::InputType::Matrix:
                {
                    matrix = va_arg(args, Matrix *);
                    if (!n) matrix->zero();
                    for (int h = 0; h < matrix->nirrep(); ++h) {
                        auto size = static_cast<size_t>(matrix->rowdim(h)) * matrix->colspi(h ^ matrix->symmetry());
                        if (size) {
                            C_DAXPY(size, coefficient, arrayPtr, 1, matrix->pointer(h)[0], 1);
                            arrayPtr += static_cast<int>(size);
                        }
                    }
                    break;
                }
                case DIISEntry::InputType::Vector:
                {
                    vector = va_arg(args, Vector *);
                    if (!n) vector->zero();
                    auto size = static_cast<size_t>(vector->dimpi().sum());
                    if (size) {
                        C_DAXPY(size, coefficient, arrayPtr, 1, vector->pointer(), 1);
                        arrayPtr += static_cast<int>(size);
                    }
                    break;
                }
                default:
                    throw SanityCheckError("Unknown input type", __FILE__, __LINE__);
            }
        }
        if (_storagePolicy == StoragePolicy::OnDisk) _subspace[n].free_vector_memory();
        va_end(args);
    }

    timer_off("New vector");

    if (print > 2) outfile->Printf("\n");
    timer_off("DIISManager::extrapolate");

    return true;
}

/**
 * Removes any vectors existing in the DIIS subspace.
 */
void DIISManager::reset_subspace() {
    _subspace.clear();
}

/**
 * Deletes the DIIS scratch file
 */
void DIISManager::delete_diis_file() {
    if (_psio->open_check(PSIF_LIBDIIS) == 0) {
        _psio->open(PSIF_LIBDIIS, PSIO_OPEN_OLD);
    }
    _psio->close(PSIF_LIBDIIS, 0);
}

DIISManager::~DIISManager() {
    if (_psio->open_check(PSIF_LIBDIIS)) _psio->close(PSIF_LIBDIIS, 1);
}

}  // namespace psi
