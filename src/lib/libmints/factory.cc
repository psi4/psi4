/*
 *  factory.cpp
 *  matrix
 *
 *  Created by Justin Turney on 4/8/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "factory.h"
#include "sobasis.h"
#include <libparallel/parallel.h>
#include "libciomr/libciomr.h"

using namespace boost;
using namespace psi;

MatrixFactory::MatrixFactory()
{
    nirrep_ = 0;
    rowspi_ = NULL;
    colspi_ = NULL;
}

MatrixFactory::MatrixFactory(const MatrixFactory& copy)
{
    nirrep_ = copy.nirrep_;
    rowspi_ = new int[nirrep_];
    colspi_ = new int[nirrep_];

    memcpy(rowspi_, copy.rowspi_, sizeof(int) * nirrep_);
    memcpy(colspi_, copy.colspi_, sizeof(int) * nirrep_);
}

MatrixFactory::~MatrixFactory()
{
    if (rowspi_)
        Chkpt::free(rowspi_);
    if (colspi_)
        Chkpt::free(colspi_);
}

bool MatrixFactory::init_with_chkpt(boost::shared_ptr<PSIO> psio)
{
    boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio.get(), PSIO_OPEN_OLD));
    bool result = init_with_chkpt(chkpt);
    return result;
}

bool MatrixFactory::init_with_chkpt(boost::shared_ptr<Chkpt> chkpt)
{
    nirrep_ = chkpt->rd_nirreps();
    rowspi_  = chkpt->rd_sopi();
    colspi_  = chkpt->rd_sopi();
    nso_     = chkpt->rd_nso();

    return true;
}

bool MatrixFactory::init_with(int nirreps, int *rowspi, int *colspi)
{
    nirrep_ = nirreps;
    rowspi_ = new int[nirrep_];
    colspi_ = new int[nirrep_];

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

    rowspi_ = new int[rows.n()];
    colspi_ = new int[cols.n()];

    nso_ = 0;
    for (int i=0; i<nirrep_; ++i) {
        rowspi_[i] = rows[i];
        colspi_[i] = cols[i];
        nso_ += rowspi_[i];
    }

    return true;
}

bool MatrixFactory::init_with(const boost::shared_ptr<SOBasisSet>& sobasis)
{
    return init_with(sobasis->dimension(), sobasis->dimension());
}

/// Returns number of irreps
int MatrixFactory::nirrep() const {
    return nirrep_;
}

/// Returns the rows per irrep array
int *MatrixFactory::rowspi() const {
    return rowspi_;
}

/// Returns the number of rows in irrep h
int MatrixFactory::nrow(int h) const {
    return rowspi_[h];
}

/// Returns the columns per irrep array
int *MatrixFactory::colspi() const {
    return colspi_;
}

/// Returns the number of columns in irrep h
int MatrixFactory::ncol(int h) const {
    return colspi_[h];
}

/// Returns the number of orbitals
int MatrixFactory::nso() const {
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
    return new Vector(nirrep_, rowspi_);
}

void MatrixFactory::create_vector(Vector& vec)
{
    vec.init(nirrep_, rowspi_);
}

/// Returns a new SimpleMatrix object with default dimensions
SimpleMatrix * MatrixFactory::create_simple_matrix()
{
    return new SimpleMatrix(nso_, nso_);
}
void MatrixFactory::create_simple_matrix(SimpleMatrix& mat)
{
    mat.init(nso_, nso_);
}

/// Returns a new SimpleMatrix object named name with default dimensions
SimpleMatrix * MatrixFactory::create_simple_matrix(std::string name)
{
    return new SimpleMatrix(name, nso_, nso_);
}
void MatrixFactory::create_simple_matrix(SimpleMatrix& mat, std::string name)
{
    mat.init(nso_, nso_, name);
}

/// Returns a new SimpleMatrix object named name of size m x n
SimpleMatrix * MatrixFactory::create_simple_matrix(std::string name, int m, int n)
{
    return new SimpleMatrix(name, m, n);
}
void MatrixFactory::create_simple_matrix(SimpleMatrix& mat, std::string name, int m, int n)
{
    mat.init(m, n, name);
}

/// Returns a new SimpleMatrix object with size m x n
SimpleMatrix * MatrixFactory::create_simple_matrix(int m, int n)
{
    return new SimpleMatrix(m, n);
}
void MatrixFactory::create_simple_matrix(SimpleMatrix& mat, int m, int n)
{
    mat.init(m, n);
}

/// Returns a new SimpleVector object with default dimension
SimpleVector * MatrixFactory::create_simple_vector()
{
    return new SimpleVector(nso_);
}

/// Returns a new SimpleVector object with size m
SimpleVector * MatrixFactory::create_simple_vector(int m)
{
    return new SimpleVector(m);
}
