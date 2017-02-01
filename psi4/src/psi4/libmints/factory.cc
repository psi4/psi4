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

/*
 *  factory.cpp
 *  matrix
 *
 *  Created by Justin Turney on 4/8/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "psi4/libmints/factory.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/libciomr/libciomr.h"

;
using namespace psi;

MatrixFactory::MatrixFactory()
{
    nirrep_ = 0;
}

MatrixFactory::MatrixFactory(const MatrixFactory& copy)
{
    nirrep_ = copy.nirrep_;
    rowspi_ = copy.rowspi_;
    colspi_ = copy.colspi_;
}

MatrixFactory::~MatrixFactory()
{
}

bool MatrixFactory::init_with(int nirreps, int *rowspi, int *colspi)
{
    nirrep_ = nirreps;
    rowspi_ = Dimension(nirrep_);
    colspi_ = Dimension(nirrep_);

    nso_ = 0;
    for (int i=0; i<nirrep_; ++i) {
        rowspi_[i] = rowspi[i];
        colspi_[i] = colspi[i];
        nso_ += rowspi_[i];
    }

    return true;
}

bool MatrixFactory::init_with(const Dimension& rows, const Dimension& cols)
{
    nirrep_ = rows.n();

    if (rows.n() != cols.n())
        throw PSIEXCEPTION("MatrixFactory can only handle same symmetry for rows and cols.");

    rowspi_ = rows;
    colspi_ = cols;

    nso_ = 0;
    for (int i=0; i<nirrep_; ++i) {
        nso_ += rowspi_[i];
    }

    return true;
}

bool MatrixFactory::init_with(const std::shared_ptr<SOBasisSet>& sobasis)
{
    return init_with(sobasis->dimension(), sobasis->dimension());
}

/// Returns number of irreps
int MatrixFactory::nirrep() const {
    return nirrep_;
}

/// Returns the rows per irrep array
const Dimension& MatrixFactory::rowspi() const {
    return rowspi_;
}

/// Returns the number of rows in irrep h
int MatrixFactory::nrow(int h) const {
    return rowspi_[h];
}

/// Returns the columns per irrep array
const Dimension& MatrixFactory::colspi() const {
    return colspi_;
}

/// Returns the number of columns in irrep h
int MatrixFactory::ncol(int h) const {
    return colspi_[h];
}

/// Returns the number of orbitals
int MatrixFactory::norb() const {
    return nso_;
}

/// Returns a new Matrix object with default dimensions
Matrix * MatrixFactory::create_matrix(int symmetry)
{
    return new Matrix(nirrep_, rowspi_, colspi_, symmetry);
}

/// Returns a new Matrix object with default dimensions
SharedMatrix MatrixFactory::create_shared_matrix()
{
    return SharedMatrix(new Matrix(nirrep_, rowspi_, colspi_));
}

void MatrixFactory::create_matrix(Matrix& mat, int symmetry)
{
    mat.init(nirrep_, rowspi_, colspi_, "", symmetry);
}

/// Returns a new Matrix object named name with default dimensions
Matrix * MatrixFactory::create_matrix(std::string name, int symmetry)
{
    return new Matrix(name, nirrep_, rowspi_, colspi_, symmetry);
}

SharedMatrix MatrixFactory::create_shared_matrix(const std::string& name)
{
    return SharedMatrix(new Matrix(name, nirrep_, rowspi_, colspi_));
}

SharedMatrix MatrixFactory::create_shared_matrix(const std::string& name, int symmetry)
{
    return SharedMatrix(new Matrix(name, nirrep_, rowspi_, colspi_, symmetry));
}

SharedMatrix MatrixFactory::create_shared_matrix(const std::string& name, int rows, int cols)
{
    return SharedMatrix(new Matrix(name, rows, cols));
}

void MatrixFactory::create_matrix(Matrix& mat, std::string name, int symmetry)
{
    mat.init(nirrep_, rowspi_, colspi_, name, symmetry);
}

/// Returns a new Vector object with default dimensions
Vector * MatrixFactory::create_vector()
{
    return new Vector(rowspi_);
}

void MatrixFactory::create_vector(Vector& vec)
{
    vec.init(rowspi_);
}

SharedVector MatrixFactory::create_shared_vector(const std::string& name)
{
    return SharedVector(new Vector(name, rowspi_));
}
