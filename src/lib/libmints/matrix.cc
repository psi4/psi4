/*
 *  matrix.cc
 *  matrix
 *
 *  Created by Justin Turney on 4/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <exception.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <libparallel/parallel.h>
#include <libmints/integral.h>
#include <libdpd/dpd.h>
#include "factory.h"
#include "wavefunction.h"
#include "dimension.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <ctype.h>
#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/xpressive/regex_actions.hpp>
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>

using namespace boost;
using namespace psi;
using namespace std;

// In molecule.cc
namespace psi {
extern int str_to_int(const std::string& s);
extern double str_to_double(const std::string& s);
}

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

namespace {
string to_string(const int val)
{
    stringstream strm;
    strm <<  val;
    return strm.str();
}
}

Matrix::Matrix()
{
    matrix_ = NULL;
    rowspi_ = NULL;
    colspi_ = NULL;
    nirrep_ = 0;
    symmetry_ = 0;
}

Matrix::Matrix(const string& name, int symmetry)
    : matrix_(0), rowspi_(0), colspi_(0), nirrep_(0),
      name_(name), symmetry_(symmetry)
{
}

Matrix::Matrix(const Matrix& c)
{
    matrix_ = NULL;
    nirrep_ = c.nirrep_;
    symmetry_ = c.symmetry_;
    rowspi_ = new int[nirrep_];
    colspi_ = new int[nirrep_];
    for (int i=0; i<nirrep_; ++i) {
        rowspi_[i] = c.rowspi_[i];
        colspi_[i] = c.colspi_[i];
    }
    alloc();
    copy_from(c.matrix_);
}

Matrix::Matrix(const boost::shared_ptr<Matrix>& c)
{
    matrix_ = NULL;
    nirrep_ = c->nirrep_;
    symmetry_ = c->symmetry_;
    rowspi_ = new int[nirrep_];
    colspi_ = new int[nirrep_];
    for (int i=0; i<nirrep_; ++i) {
        rowspi_[i] = c->rowspi_[i];
        colspi_[i] = c->colspi_[i];
    }
    alloc();
    copy_from(c->matrix_);
}

Matrix::Matrix(const Matrix* c)
{
    matrix_ = NULL;
    nirrep_ = c->nirrep_;
    symmetry_ = c->symmetry_;
    rowspi_ = new int[nirrep_];
    colspi_ = new int[nirrep_];
    for (int i=0; i<nirrep_; ++i) {
        rowspi_[i] = c->rowspi_[i];
        colspi_[i] = c->colspi_[i];
    }
    alloc();
    copy_from(c->matrix_);
}

Matrix::Matrix(int l_nirreps, const int *l_rowspi, const int *l_colspi, int symmetry)
{
    matrix_ = NULL;
    nirrep_ = l_nirreps;
    symmetry_ = symmetry;
    rowspi_ = new int[nirrep_];
    colspi_ = new int[nirrep_];
    for (int i=0; i<nirrep_; ++i) {
        rowspi_[i] = l_rowspi[i];
        colspi_[i] = l_colspi[i];
    }
    alloc();
}

Matrix::Matrix(const string& name, int l_nirreps, const int *l_rowspi, const int *l_colspi, int symmetry) : name_(name)
{
    matrix_ = NULL;
    nirrep_ = l_nirreps;
    symmetry_ = symmetry;
    rowspi_ = new int[nirrep_];
    colspi_ = new int[nirrep_];
    for (int i=0; i<nirrep_; ++i) {
        rowspi_[i] = l_rowspi[i];
        colspi_[i] = l_colspi[i];
    }
    alloc();
}

Matrix::Matrix(const string& name, int rows, int cols) : name_(name)
{
    matrix_ = NULL;
    nirrep_ = 1;
    symmetry_ = 0;
    rowspi_ = new int[nirrep_];
    colspi_ = new int[nirrep_];
    rowspi_[0] = rows;
    colspi_[0] = cols;
    alloc();
}

Matrix::Matrix(int rows, int cols)
{
    matrix_ = NULL;
    nirrep_ = 1;
    symmetry_ = 0;
    rowspi_ = new int[nirrep_];
    colspi_ = new int[nirrep_];
    rowspi_[0] = rows;
    colspi_[0] = cols;
    alloc();
}

Matrix::Matrix(int nirrep, int rows, const int *colspi)
{
    matrix_ = NULL;
    symmetry_ = 0;
    nirrep_ = nirrep;

    rowspi_ = new int[nirrep_];
    colspi_ = new int[nirrep_];

    for (int i=0; i<nirrep_; ++i) {
        rowspi_[i] = rows;
        colspi_[i] = colspi[i];
    }
    alloc();
}

Matrix::Matrix(int nirrep, const int *rowspi, int cols)
{
    matrix_ = NULL;
    symmetry_ = 0;
    nirrep_ = nirrep;

    rowspi_ = new int[nirrep_];
    colspi_ = new int[nirrep_];

    for (int i=0; i<nirrep_; ++i) {
        rowspi_[i] = rowspi[i];
        colspi_[i] = cols;
    }
    alloc();
}

Matrix::Matrix(const string& name, const Dimension& rows, const Dimension& cols, int symmetry)
{
    name_ = name;
    matrix_ = NULL;
    symmetry_ = symmetry;

    // This will happen in PetiteList::aotoso()
    if (rows.n() == 1) {
        nirrep_ = cols.n();
        rowspi_ = new int[nirrep_];
        colspi_ = new int[nirrep_];
        for (int i=0; i<nirrep_; ++i) {
            rowspi_[i] = rows[0];
            colspi_[i] = cols[i];
        }
    }
    else {
        nirrep_ = rows.n();
        rowspi_ = new int[nirrep_];
        colspi_ = new int[nirrep_];
        for (int i=0; i<nirrep_; ++i) {
            rowspi_[i] = rows[i];
            colspi_[i] = cols[i];
        }
    }

    alloc();
}

Matrix::Matrix(dpdfile2 *inFile) : name_(inFile->label)
{
    dpd_file2_mat_init(inFile);
    dpd_file2_mat_rd(inFile);
    matrix_ = NULL;
    symmetry_ = inFile->my_irrep;
    nirrep_ = inFile->params->nirreps;
    rowspi_ = new int[nirrep_];
    colspi_ = new int[nirrep_];
    for (int i=0; i<nirrep_; ++i) {
        rowspi_[i] = inFile->params->rowtot[i];
        colspi_[i] = inFile->params->coltot[i];
    }
    alloc();
    copy_from(inFile->matrix);
    dpd_file2_mat_close(inFile);
}

Matrix::~Matrix()
{
    release();
    if (rowspi_)
        delete[] rowspi_;
    if (colspi_)
        delete[] colspi_;
}

/// allocate a block matrix -- analogous to libciomr's block_matrix
double** Matrix::matrix(int nrow, int ncol)
{
    double** mat = (double**) malloc(sizeof(double*)*nrow);
    const size_t size = sizeof(double)*nrow*ncol;
    mat[0] = (double*) malloc(size);
    memset((void *)mat[0], 0, size);
    for(int r=1; r<nrow; ++r) mat[r] = mat[r-1] + ncol;
    return mat;
}
/// free a (block) matrix -- analogous to libciomr's free_block
void Matrix::free(double** Block)
{
    ::free(Block[0]);  ::free(Block);
}

void Matrix::init(int l_nirreps, const int *l_rowspi, const int *l_colspi, const string& name, int symmetry)
{
    if (rowspi_) delete[] rowspi_;
    if (colspi_) delete[] colspi_;
    name_ = name;
    symmetry_ = symmetry;
    nirrep_ = l_nirreps;
    rowspi_ = new int[nirrep_];
    colspi_ = new int[nirrep_];
    for (int i=0; i<nirrep_; ++i) {
        rowspi_[i] = l_rowspi[i];
        colspi_[i] = l_colspi[i];
    }
    alloc();
}

boost::shared_ptr<Matrix> Matrix::create(const std::string& name,
                                 int nirrep,
                                 int* rows,
                                 int *cols)
{
    return boost::shared_ptr<Matrix>(new Matrix(name, nirrep, rows, cols));
}

Matrix* Matrix::clone() const
{
    Matrix *temp = new Matrix(this);
    return temp;
}

void Matrix::copy(const Matrix* cp)
{
    // Make sure we are the same size as cp
    bool same = true;
    if (nirrep_ != cp->nirrep_ || symmetry_ != cp->symmetry_) {
        same = false;
    } else {
        for (int h=0; h<nirrep_; ++h)
            if (colspi_[h] != cp->colspi_[h] || rowspi_[h] != cp->rowspi_[h])
                same = false;
    }

    if (same == false) {
        release();
        if (rowspi_)
            delete[] rowspi_;
        if (colspi_)
            delete[] colspi_;
        nirrep_ = cp->nirrep_;
        symmetry_ = cp->symmetry_;
        rowspi_ = new int[nirrep_];
        colspi_ = new int[nirrep_];
        for (int i=0; i<nirrep_; ++i) {
            rowspi_[i] = cp->rowspi_[i];
            colspi_[i] = cp->colspi_[i];
        }
        alloc();
    }

    // When here we are the same size
    for (int h=0; h<nirrep_; ++h) {
        if (rowspi_[h] != 0 && colspi_[h^symmetry_] != 0)
            memcpy(&(matrix_[h][0][0]), &(cp->matrix_[h][0][0]), rowspi_[h] * colspi_[h^symmetry_] * sizeof(double));
    }
}

void Matrix::copy(const Matrix& cp)
{
    copy(&cp);
}

void Matrix::copy(const boost::shared_ptr<Matrix>& cp)
{
    copy(cp.get());
}

void Matrix::alloc()
{
    if (matrix_)
        release();

    matrix_ = (double***)malloc(sizeof(double***) * nirrep_);
    for (int i=0; i<nirrep_; ++i) {
        if (rowspi_[i] != 0 && colspi_[i^symmetry_] != 0)
            matrix_[i] = Matrix::matrix(rowspi_[i], colspi_[i^symmetry_]);
        else
            matrix_[i] = NULL;
    }
}

void Matrix::release()
{
    if (!matrix_)
        return;

    for (int h=0; h<nirrep_; ++h) {
        if (matrix_[h])
            Matrix::free(matrix_[h]);
    }
    ::free(matrix_);
    matrix_ = NULL;
}

void Matrix::copy_from(double ***c) {
    int size;

    for (int h=0; h<nirrep_; ++h) {
        size = rowspi_[h] * colspi_[h^symmetry_] * sizeof(double);
        if (size)
            memcpy(&(matrix_[h][0][0]), &(c[h][0][0]), size);
    }
}

// Sets all elements of matrix to val
void Matrix::set(double val)
{
    for (int h=0; h < nirrep_; ++h) {
        size_t size = rowspi_[h] * colspi_[h^symmetry_];

        for (size_t i=0; i<size; ++i) {
            matrix_[h][0][i] = val;
        }
    }
}

void Matrix::set(const double * const tri)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::set called on a non-totally symmetric matrix.");
    }
    int h, i, j, ii, jj;
    int offset;

    offset = 0;
    for (h=0; h<nirrep_; ++h) {
        for (i=0; i<rowspi_[h]; ++i) {
            ii = i + offset;
            for (j=0; j<=i; ++j) {
                jj = j + offset;
                matrix_[h][i][j] = matrix_[h][j][i] = tri[ii*(ii+1)/2 + jj];
            }
        }
        offset += rowspi_[h];
    }
}

void Matrix::set(const double * const * const sq)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::set called on a non-totally symmetric matrix.");
    }

    int h, i, j, ii, jj;
    int offset;

    if (sq == NULL) {
        throw PSIEXCEPTION("Matrix::set: Set call with a NULL double** matrix");
    }
    offset = 0;
    for (h=0; h<nirrep_; ++h) {
        for (i=0; i<rowspi_[h]; ++i) {
            ii = i + offset;
            for (j=0; j<=i; ++j) {
                jj = j + offset;
                matrix_[h][i][j] = sq[ii][jj];
                matrix_[h][j][i] = sq[jj][ii];
            }
        }
        offset += rowspi_[h];
    }
}

void Matrix::set(const SimpleMatrix * const sq)
{
    set(sq->const_pointer());
}

void Matrix::set(const boost::shared_ptr<SimpleMatrix>& sq)
{
    set(sq->const_pointer());
}

void Matrix::set_diagonal(const Vector * const vec)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::set_diagonal called on a non-totally symmetric matrix.");
    }

    int h, i, size;
    zero();
    for (h=0; h<nirrep_; ++h) {
        size = rowspi_[h];
        if (size) {
            for (i=0; i<size; ++i)
                matrix_[h][i][i] = vec->vector_[h][i];
        }
    }
}

void Matrix::set_diagonal(const Vector& vec)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::set_diagonal called on a non-totally symmetric matrix.");
    }

    int h, i, size;
    zero();
    for (h=0; h<nirrep_; ++h) {
        size = rowspi_[h];
        if (size) {
            for (i=0; i<size; ++i)
                matrix_[h][i][i] = vec.vector_[h][i];
        }
    }
}

void Matrix::set_diagonal(const boost::shared_ptr<Vector>& vec)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::set_diagonal called on a non-totally symmetric matrix.");
    }

    int h, i, size;
    zero();
    for (h=0; h<nirrep_; ++h) {
        size = rowspi_[h];
        if (size) {
            for (i=0; i<size; ++i)
                matrix_[h][i][i] = vec->vector_[h][i];
        }
    }
}

double *Matrix::to_lower_triangle() const
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::set called on a non-totally symmetric matrix.");
    }

    int sizer=0, sizec=0;
    for (int h=0; h<nirrep_; ++h) {
        sizer += rowspi_[h];
        sizec += colspi_[h];
    }
    if (sizer != sizec)
        return NULL;

    double *tri = new double[ioff[sizer]];
    double **temp = to_block_matrix();
    sq_to_tri(temp, tri, sizer);
    free_block(temp);
    return tri;
}

double **Matrix::to_block_matrix() const
{
    int sizer=0, sizec=0;
    for (int h=0; h<nirrep_; ++h) {
        sizer += rowspi_[h];
        sizec += colspi_[h^symmetry_];
    }

    int *col_offset = new int[nirrep_];
    col_offset[0] = 0;
    for (int h=1; h<nirrep_; ++h) {
        col_offset[h] = col_offset[h-1] + colspi_[h-1];
    }

    double **temp = block_matrix(sizer,sizec);
    int offsetr = 0, offsetc=0;
    for (int h=0; h <nirrep_; ++h) {
        offsetc = col_offset[h^symmetry_];
        for (int i=0; i<rowspi_[h]; ++i) {
            for (int j=0; j<colspi_[h^symmetry_]; ++j) {
                temp[i+offsetr][j+offsetc] = matrix_[h][i][j];
            }
        }
        offsetr += rowspi_[h];
//        offsetc += colspi_[h^symmetry_];
    }

    delete[] col_offset;
    return temp;
}

SimpleMatrix *Matrix::to_simple_matrix() const
{
    return new SimpleMatrix(this);
}

void Matrix::print_mat(const double *const *const a, int m, int n, FILE *out) const
{
    const int print_ncol = 5;
    int num_frames = int(n/print_ncol);
    int num_frames_rem = n%print_ncol; //adding one for changing 0->1 start
    int num_frame_counter = 0;
    //for each frame
    for(num_frame_counter=0;num_frame_counter<num_frames;num_frame_counter++){
        fprintf(out,"\n");
        for(int j=print_ncol*num_frame_counter+1;j<print_ncol*num_frame_counter+print_ncol+1;j++){
            if(j==print_ncol*num_frame_counter+1){ fprintf(out,"%18d",j); }
            else{ fprintf(out,"        %5d",j); }
        }
        fprintf(out,"\n\n");

        for(int k=1; k<=m; ++k){
            for(int j=print_ncol*num_frame_counter+1;j<print_ncol*num_frame_counter+print_ncol+2;j++){
                if(j==print_ncol*num_frame_counter+1){ fprintf(out,"%5d",k);}
                else{ fprintf(out," %12.7f",a[k-1][j-2]); }
            }
            fprintf(out,"\n");
        }
    }

    // ALREADY DID THE FULL FRAMES BY THIS POINT
    // NEED TO TAKE CARE OF THE REMAINDER
    if(num_frames_rem != 0){
        fprintf(out,"\n");
        for(int j=print_ncol*num_frame_counter+1;j<=n;j++){
            if(j==print_ncol*num_frame_counter+1){ fprintf(out,"%18d",j); }
            else{ fprintf(out,"        %5d",j); }
        }
        fprintf(out,"\n\n");

        for(int k=1; k<=m; ++k){
            for(int j=print_ncol*num_frame_counter+1;j<n+2;j++){
                if(j==print_ncol*num_frame_counter+1){ fprintf(out,"%5d",k); }
                else{ fprintf(out," %12.7f",a[k-1][j-2]); }
            }
            fprintf(out,"\n");
        }
    }
    fprintf(out,"\n\n");
    //R.I.P. goto statements - Aug 4th 2010 - MSM
}

void Matrix::print(FILE *out, const char *extra) const
{
    int h;

    if (name_.length()) {
        if (extra == NULL)
            fprintf(out, "  ## %s (Symmetry %d) ##\n", name_.c_str(), symmetry_);
        else
            fprintf(out, "  ## %s %s (Symmetry %d)##\n", name_.c_str(), extra, symmetry_);
    }

    for (h=0; h<nirrep_; ++h) {
        fprintf(out, "  Irrep: %d\n", h+1);
        if (rowspi_[h] == 0 || colspi_[h^symmetry_] == 0)
            fprintf(out, "\n\t(empty)\n");
        else
            print_mat(matrix_[h], rowspi_[h], colspi_[h^symmetry_], out);
        fprintf(out, "\n");
    }
    fflush(out);
}

void Matrix::print_atom_vector(FILE *out)
{
    int i;

    if (name_.length()) {
        fprintf(out,"\n  -%s:\n", name_.c_str());
    }
    fprintf(out,"     Atom            X                  Y                   Z\n");
    fprintf(out,"    ------   -----------------  -----------------  -----------------\n");

    for(i=0;i<nrow();i++) {
        fprintf(out,"    %4d   ",i+1);
        fprintf(out,"  %17.12lf  %17.12lf  %17.12lf", matrix_[0][i][0], matrix_[0][i][1], matrix_[0][i][2]);
        fprintf(out,"\n");
    }
    fprintf(out,"\n");
}


void Matrix::eivprint(const Vector * const values, FILE *out)
{
    if (symmetry_)
        throw PSIEXCEPTION("Matrix::eivprint: This print does not make sense for non-totally symmetric matrices.");

    int h;

    if (name_.length()) {
        fprintf(out, "  ## %s with eigenvalues ##\n", name_.c_str());
    }

    for (h=0; h<nirrep_; ++h) {
        fprintf(out, " Irrep: %d\n", h+1);
        eivout(matrix_[h], values->vector_[h], rowspi_[h], colspi_[h^symmetry_], out);
        fprintf(out, "\n");
    }
    fflush(out);
}

void Matrix::eivprint(const Vector& values, FILE *out)
{
    eivprint(&values, out);
}

void Matrix::eivprint(const boost::shared_ptr<Vector>& values, FILE *out)
{
    eivprint(values.get(), out);
}

void Matrix::identity()
{
    if (symmetry_)
        return;

    int h;
    size_t size;

    for (h=0; h<nirrep_; ++h) {
        size = rowspi_[h] * colspi_[h] * sizeof(double);

        if (size) {
            memset(&(matrix_[h][0][0]), 0, size);
            for (int i=0; i<MIN(rowspi_[h], colspi_[h]); ++i)
                matrix_[h][i][i] = 1.0;
        }
    }
}

void Matrix::zero()
{
    size_t size;
    int h;

    for (h=0; h<nirrep_; ++h) {
        size = rowspi_[h] * colspi_[h^symmetry_] * sizeof(double);

        if (size) {
            memset(&(matrix_[h][0][0]), 0, size);
        }
    }
}

void Matrix::zero_diagonal()
{
    if (symmetry_)
        return;

    int h, i;

    for (h=0; h<nirrep_; ++h) {
        for (i=0; i<MIN(rowspi_[h], colspi_[h]); ++i) {
            matrix_[h][i][i] = 0.0;
        }
    }
}

double Matrix::trace()
{
    if (symmetry_)
        return 0.0;

    int i, h;
    double val = (double)0.0;

    for (h=0; h<nirrep_; ++h) {
        for (i=0; i<MIN(rowspi_[h], colspi_[h]); ++i) {
            val += matrix_[h][i][i];
        }
    }

    return val;
}

Matrix* Matrix::transpose()
{
    Matrix *temp = new Matrix(name_, nirrep_, colspi_, rowspi_, symmetry_);

    int h, i, j;
    for (h=0; h<nirrep_; ++h) {
        for (i=0; i<rowspi_[h]; ++i) {
            for (j=0; j<colspi_[h^symmetry_]; ++j) {
                temp->matrix_[h][j][i] = matrix_[h][i][j];
            }
        }
    }
    return temp;
}

void Matrix::add(const Matrix * const plus)
{
    double *lhs, *rhs;
    for (int h=0; h<nirrep_; ++h) {
        size_t size = rowspi_[h] * colspi_[h^symmetry_];
        if (size) {
            lhs = matrix_[h][0];
            rhs = plus->matrix_[h][0];
            for (size_t ij=0; ij<size; ++ij) {
                *lhs += *rhs;
                lhs++; rhs++;
            }
        }
    }
}

void Matrix::add(const Matrix& plus)
{
    add(&plus);
}

void Matrix::add(const boost::shared_ptr<Matrix>& plus)
{
    add(plus.get());
}

void Matrix::subtract(const Matrix* const plus)
{
    double *lhs, *rhs;
    for (int h=0; h<nirrep_; ++h) {
        size_t size = rowspi_[h] * colspi_[h^symmetry_];
        if (size) {
            lhs = matrix_[h][0];
            rhs = plus->matrix_[h][0];
            for (size_t ij=0; ij<size; ++ij) {
                *lhs -= *rhs;
                lhs++; rhs++;
            }
        }
    }
}

void Matrix::subtract(const boost::shared_ptr<Matrix>& sub)
{
    subtract(sub.get());
}

void Matrix::apply_denominator(const Matrix * const plus)
{
    double *lhs, *rhs;
    for (int h=0; h<nirrep_; ++h) {
        size_t size = rowspi_[h] * colspi_[h^symmetry_];
        if (size) {
            lhs = matrix_[h][0];
            rhs = plus->matrix_[h][0];
            for (size_t ij=0; ij<size; ++ij) {
                *lhs /= *rhs;
                lhs++; rhs++;
            }
        }
    }
}

void Matrix::apply_denominator(const Matrix& plus)
{
    apply_denominator(&plus);
}

void Matrix::apply_denominator(const boost::shared_ptr<Matrix>& plus)
{
    apply_denominator(plus.get());
}

void Matrix::accumulate_product(const Matrix* const a, const Matrix* const b)
{
    gemm(false, false, 1.0, a, b, 1.0);
}

void Matrix::accumulate_product(const boost::shared_ptr<Matrix>& a,
                                const boost::shared_ptr<Matrix>& b)
{
    gemm(false, false, 1.0, a, b, 1.0);
}

void Matrix::scale(double a)
{
    int h;
    size_t size;
    for (h=0; h<nirrep_; ++h) {
        size = rowspi_[h] * colspi_[h^symmetry_];
        if (size)
            C_DSCAL(size, a, &(matrix_[h][0][0]), 1);
    }
}

void Matrix::scale_row(int h, int m, double a)
{
    C_DSCAL(colspi_[h^symmetry_], a, &(matrix_[h][m][0]), 1);
}

void Matrix::scale_column(int h, int n, double a)
{
    C_DSCAL(rowspi_[h], a, &(matrix_[h][0][n]), colspi_[h^symmetry_]);
}

double Matrix::sum_of_squares()
{
    double sum = (double)0.0;
    for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<rowspi_[h]; ++i) {
            for (int j=0; j<colspi_[h^symmetry_]; ++j) {
                sum += matrix_[h][i][j] * matrix_[h][i][j];
            }
        }
    }

    return sum;
}

double Matrix::rms()
{
    double sum = (double)0.0;
    long terms = 0;
    for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<rowspi_[h]; ++i) {
            for (int j=0; j<colspi_[h^symmetry_]; ++j) {
                sum += matrix_[h][i][j] * matrix_[h][i][j];
                terms++;
            }
        }
    }

    return sqrt(sum/terms);
}

void Matrix::transform(const Matrix* const a, const Matrix* const transformer)
{
    Matrix temp(a);

    temp.gemm(false, false, 1.0, a, transformer, 0.0);
    gemm(true, false, 1.0, transformer, &temp, 0.0);
}

void Matrix::transform(const boost::shared_ptr<Matrix>& a, const boost::shared_ptr<Matrix>& transformer)
{
    transform(a.get(), transformer.get());
}

void Matrix::transform(const Matrix* const transformer)
{
    bool square = true;
    int h = 0;

    while(h < nirrep_ && square){
        if(transformer->rowspi()[h] != transformer->colspi()[h]){
            square = false;
        }
        ++h;
    }

    if(square){
        Matrix temp(this);
        temp.gemm(false, false, 1.0, this, transformer, 0.0);
        gemm(true, false, 1.0, transformer, &temp, 0.0);
    }
    else{
        Matrix temp(nirrep_, rowspi_, transformer->colspi());
        Matrix result(nirrep_, transformer->colspi(), transformer->colspi());
        temp.gemm(false, false, 1.0, this, transformer, 0.0);
        result.gemm(true, false, 1.0, transformer, &temp, 0.0);
        copy(&result);
    }
}

void Matrix::transform(const boost::shared_ptr<Matrix>& transformer)
{
    transform(transformer.get());
}

void Matrix::transform(const boost::shared_ptr<Matrix>& L,
                       const boost::shared_ptr<Matrix>& F,
                       const boost::shared_ptr<Matrix>& R)
{
    Matrix temp(nirrep_, F->rowspi_, R->colspi_, F->symmetry_ ^ R->symmetry_);
    temp.gemm(false, false, 1.0, F, R, 0.0);
    gemm(true, false, 1.0, L, temp, 0.0);
}

void Matrix::back_transform(const Matrix* const a, const Matrix* const transformer)
{
    Matrix temp(a);

    temp.gemm(false, true, 1.0, a, transformer, 0.0);
    gemm(false, false, 1.0, transformer, &temp, 0.0);
}

void Matrix::back_transform(const boost::shared_ptr<Matrix>& a, const boost::shared_ptr<Matrix>& transformer)
{
    back_transform(a.get(), transformer.get());
}

void Matrix::back_transform(const Matrix* const transformer)
{
    bool square = true;
    int h = 0;

    while(h < nirrep_ && square){
        if(transformer->rowspi()[h] != transformer->colspi()[h]){
            square = false;
        }
        h++;
    }

    if(square){
        Matrix temp(this);
        temp.gemm(false, true, 1.0, this, transformer, 0.0);
        gemm(false, false, 1.0, transformer, &temp, 0.0);
    }
    else{
        Matrix temp(nirrep_, rowspi_, transformer->rowspi());
        Matrix result(nirrep_, transformer->rowspi(), transformer->rowspi());
        temp.gemm(false, true, 1.0, this, transformer, 0.0);
        result.gemm(false, false, 1.0, transformer, &temp, 0.0);
        copy(&result);
    }
}

void Matrix::back_transform(const boost::shared_ptr<Matrix>& transformer)
{
    back_transform(transformer.get());
}

void Matrix::gemm(bool transa, bool transb, double alpha, const Matrix* const a,
                  const Matrix* const b, double beta)
{
    // Check symmetry
    if (symmetry_ != (a->symmetry_ ^ b->symmetry_)) {
        fprintf(outfile, "Matrix::gemm error: Input symmetries will not result in target symmetry.\n");
        fprintf(outfile, " Asym %d ^ Bsym %d != Csym %d\n", a->symmetry(), b->symmetry(), symmetry());
        fprintf(outfile, "Result is %d\n", a->symmetry_ ^ b->symmetry_);
        throw PSIEXCEPTION("Matrix::gemm error: Input symmetries will not result in target symmetry.");
    }

    if (transa && a->symmetry_)
        throw PSIEXCEPTION("Matrix::gemm error: a is non totally symmetric and you're trying to transpose it");
    if (transb && b->symmetry_)
        throw PSIEXCEPTION("Matrix::gemm error: b is non totally symmetric and you're trying to transpose it");

    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int h, m, n, k, lda, ldb, ldc;

    for (h=0; h<nirrep_; ++h) {
        m = rowspi_[h];
        n = colspi_[h^symmetry_];
        k = transa ? a->rowspi_[h] : a->colspi_[h^a->symmetry_];
        lda = transa ? m : k;
        ldb = transb ? k : n;
        ldc = n;

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, &(a->matrix_[h][0][0]),
                    lda, &(b->matrix_[h^symmetry_^b->symmetry_][0][0]), ldb, beta, &(matrix_[h][0][0]),
                    ldc);
        }
    }
}

void Matrix::gemm(bool transa, bool transb, double alpha,
                  const boost::shared_ptr<Matrix>& a, const boost::shared_ptr<Matrix>& b,
                  double beta)
{
    gemm(transa, transb, alpha, a.get(), b.get(), beta);
}

void Matrix::gemm(bool transa, bool transb, double alpha,
                  const Matrix& a, const boost::shared_ptr<Matrix>& b,
                  double beta)
{
    gemm(transa, transb, alpha, &a, b.get(), beta);
}

void Matrix::gemm(bool transa, bool transb, double alpha,
                  const boost::shared_ptr<Matrix>& a, const Matrix& b,
                  double beta)
{
    gemm(transa, transb, alpha, a.get(), &b, beta);
}

void Matrix::gemm(bool transa, bool transb, double alpha,
                  const Matrix& a, const Matrix& b, double beta)
{
    gemm(transa, transb, alpha, &a, &b, beta);
}

bool Matrix::schmidt_add(int rows, Vector& v) throw()
{
    if (nirrep_ > 1 || v.nirrep() > 1)
        throw PSIEXCEPTION("Matrix::schmidt_add: This function needs to be adapted to handle symmetry blocks.");

    double dotval, normval;
    int i, I;

    for (i=0; i<rows; ++i) {
        dotval = C_DDOT(colspi_[0], matrix_[0][i], 1, v.pointer(), 1);
        for (I=0; I<colspi_[0]; ++I)
            v(I) -= dotval * matrix_[0][i][I];
    }

    normval = C_DDOT(colspi_[0], v.pointer(), 1, v.pointer(), 1);
    normval = sqrt(normval);

    if (normval > 1.0e-5) {
        for (I=0; I<colspi_[0]; ++I)
            matrix_[0][rows][I] = v(I) / normval;
        return true;
    }
    else
        return false;
}

bool Matrix::schmidt_add(int rows, double* v) throw()
{
    if (nirrep_ > 1)
        throw PSIEXCEPTION("Matrix::schmidt_add: This function needs to be adapted to handle symmetry blocks.");

    double dotval, normval;
    int i, I;

    for (i=0; i<rows; ++i) {
        dotval = C_DDOT(colspi_[0], matrix_[0][i], 1, v, 1);
        for (I=0; I<colspi_[0]; ++I)
            v[I] -= dotval * matrix_[0][i][I];
    }

    normval = C_DDOT(colspi_[0], v, 1, v, 1);
    normval = sqrt(normval);

    if (normval > 1.0e-5) {
        for (I=0; I<colspi_[0]; ++I)
            matrix_[0][rows][I] = v[I] / normval;
        return true;
    }
    else
        return false;
}

void Matrix::project_out(Matrix &constraints)
{
    // We're going to work through temp and add to this
    Matrix temp = *this;
    zero();

    double *v = new double[coldim()];
    for (int i=0; i<rowdim(); ++i) {
        memcpy(v, temp[0][i], sizeof(double)*coldim());
        for (int j=0; j<constraints.rowdim(); ++j) {
            double dotval = C_DDOT(coldim(), temp[0][i], 1, constraints[0][j], 1);
            for (int I=0; I<coldim(); ++I)
                v[I] -= dotval * constraints[0][j][I];
        }

        // At this point all constraints have been projected out of "v"
        // Normalize it add Schmidt orthogonalize it against this
        double normval = C_DDOT(coldim(), v, 1, v, 1);
        if (normval > 1.0E-10) {
            for (int j=0; j<coldim(); ++j)
                v[j] /= normval;

            schmidt_add(i, v);
        }
    }
    delete[] v;
}

double Matrix::vector_dot(const Matrix* const rhs)
{
    if (symmetry_ != rhs->symmetry_)
        return 0.0;

    double sum = 0.0;
    int h;
    size_t size;

    for (h=0; h<nirrep_; ++h) {
        size = rowspi_[h] * colspi_[h^symmetry_];
        // Check the size of the other
        if (size != rhs->rowdim(h) * rhs->coldim(h^symmetry_))
            throw PSIEXCEPTION("Matrix::vector_dot: Dimensions do not match!\n");

        if (size)
            sum += C_DDOT(size, (&matrix_[h][0][0]), 1, &(rhs->matrix_[h][0][0]), 1);
    }

    return sum;
}

double Matrix::vector_dot(const boost::shared_ptr<Matrix>& rhs)
{
    return vector_dot(rhs.get());
}

void Matrix::diagonalize(Matrix* eigvectors, Vector* eigvalues, DiagonalizeOrder nMatz)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::diagonalize: Matrix is non-totally symmetric.");
    }
    int h;
    for (h=0; h<nirrep_; ++h) {
        if (rowspi_[h]) {
            sq_rsp(rowspi_[h], colspi_[h], matrix_[h], eigvalues->vector_[h], static_cast<int>(nMatz), eigvectors->matrix_[h], 1.0e-14);
        }
    }
}

void Matrix::diagonalize(boost::shared_ptr<Matrix>& eigvectors, boost::shared_ptr<Vector>& eigvalues, DiagonalizeOrder nMatz)
{
    diagonalize(eigvectors.get(), eigvalues.get(), nMatz);
}

void Matrix::diagonalize(boost::shared_ptr<Matrix>& eigvectors, Vector& eigvalues, DiagonalizeOrder nMatz)
{
    diagonalize(eigvectors.get(), &eigvalues, nMatz);
}

void Matrix::diagonalize(boost::shared_ptr<Matrix>& metric, boost::shared_ptr<Matrix>& eigvectors, boost::shared_ptr<Vector>& eigvalues, DiagonalizeOrder nMatz)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::diagonalize: Matrix non-totally symmetric.");
    }

    // this and metric are destroyed in the process, so let's make a copy
    // that we work with.
    Matrix t(*this);
    Matrix m(metric);

    int lwork = 3*max_nrow();
    double *work = new double[lwork];

    for (int h=0; h<nirrep_; ++h) {
        if (!rowspi_[h] && !colspi_[h])
            continue;

        int err = C_DSYGV(1, 'V', 'U',
                          rowspi_[h], t.matrix_[h][0],
                          rowspi_[h], m.matrix_[h][0],
                          rowspi_[h], eigvalues->pointer(h),
                          work, lwork);

        if (err != 0) {
            if (err < 0) {
                fprintf(outfile, "Matrix::diagonalize with metric: C_DSYGV: argument %d has invalid parameter.\n", -err);
                fflush(outfile);
                abort();
            }
            if (err > 0) {
                fprintf(outfile, "Matrix::diagonalize with metric: C_DSYGV: error value: %d\n", err);
                fflush(outfile);
                abort();
            }
        }

        // TODO: Sort the data according to eigenvalues.
    }
    delete[] work;
}

void Matrix::swap_rows(int h, int i, int j)
{
    C_DSWAP(colspi_[h], &(matrix_[h][i][0]), 1, &(matrix_[h][j][0]), 1);
}

void Matrix::swap_columns(int h, int i, int j)
{
    C_DSWAP(rowspi_[h], &(matrix_[h][0][i]), colspi_[h], &(matrix_[h][0][j]), colspi_[h]);
}

void Matrix::cholesky_factorize()
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::cholesky_factorize: Matrix is non-totally symmetric.");
    }
    for (int h=0; h<nirrep_; ++h) {
        if (rowspi_[h]) {
            int err = C_DPOTRF('L', rowspi_[h], matrix_[h][0], rowspi_[h]);
            if (err != 0) {
                if (err < 0) {
                    fprintf(outfile, "cholesky_factorize: C_DPOTRF: argument %d has invalid paramter.\n", -err);
                    fflush(outfile);
                    abort();
                }
                if (err > 1) {
                    fprintf(outfile, "cholesky_factorize: C_DPOTRF: the leading minor of order %d is not "
                            "positive definite, and the factorization could not be "
                            "completed.", err);
                    fflush(outfile);
                    abort();
                }
            }
        }
    }
}

boost::shared_ptr<Matrix> Matrix::partial_cholesky_factorize(double delta, bool throw_if_negative)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::partial_cholesky_factorize: Matrix is non-totally symmetric.");
    }

    // Temporary cholesky factor (full memory)
    boost::shared_ptr<Matrix> K(new Matrix("L Temp", nirrep_, rowspi_, rowspi_));

    // Significant Cholesky columns per irrep
    int *sigpi = new int[nirrep_];
    ::memset(static_cast<void*>(sigpi), '\0', nirrep_*sizeof(int));

    for (int h=0; h<nirrep_; ++h) {

        if (!rowspi_[h]) continue;

        // Dimension
        int n = rowspi_[h];
        // Cholesky factor
        double** Kp = K->pointer(h);
        // Original matrix
        double** Ap = matrix_[h];

        // Diagonal (or later Schur complement diagonal)
        double* Dp = new double[n];
        ::memset(static_cast<void*>(Dp), '\0', nirrep_*sizeof(double));
        for (int i = 0; i < n; i++)
            Dp[i] = Ap[i][i];

        // Vector of completed columns (absolute)
        std::vector<int> order;

        int nQ = 0;
        int Q = -1;
        while (nQ < n) {

            // Find max error on diagonal
            // (Should always be in the Schur complement)
            int imax = 0;
            for (int i = 0; i < n; i++)
                if (fabs(Dp[i]) > fabs(Dp[imax]))
                    imax = i;

            double dmax = Dp[imax];
            if (fabs(dmax) <= delta)
                break;

            if (dmax <= 0.0) {
                if (throw_if_negative)
                    throw PSIEXCEPTION("Matrix::partial_cholesky_factorize: Pivot is numerically negative or zero");
                else
                    break;
            }

            // New vector!
            nQ++;
            Q++;

            // Find the diagonal
            double diag = sqrt(dmax);

            // Update the vector
            C_DCOPY(n,&Ap[0][imax],n,&Kp[0][Q],n);
            C_DGEMV('N',n,nQ-1,-1.0,Kp[0],n,Kp[imax],1,1.0,&Kp[0][Q],n);
            C_DSCAL(n,1.0 / diag, &Kp[0][Q], n);

            // Explicitly zero out elements of the vector
            // Which are psychologically upper triangular
            for (int i = 0; i < order.size(); i++)
                Kp[order[i]][Q] = 0.0;

            // Place the diagonal
            Kp[imax][Q] = diag;

            // Update the Schur complement
            for (int i = 0; i < n; i++)
                Dp[i] -= Kp[i][Q] * Kp[i][Q];

            // Explicitly zero out elements of the Schur complement
            // Which are already exact, and do not really belong
            // This prevents false selection due to roundoff
            Dp[imax] = 0.0;

            // Add the diagonal index to the list of completed indices
            order.push_back(imax);
        }
        sigpi[h] = nQ;
    }

    // Copy out to properly sized array
    boost::shared_ptr<Matrix> L(new Matrix("Partial Cholesky Factor", nirrep_, rowspi_, sigpi));
    for (int h = 0; h < nirrep_; h++) {
        if (!rowspi_[h]) continue;
        double** Kp = K->pointer(h);
        double** Lp = L->pointer(h);

        for (int i = 0; i < rowspi_[h]; i++) {
            ::memcpy(static_cast<void*>(Lp[i]),static_cast<void*>(Kp[i]),sizeof(double)*sigpi[h]);
        }
    }

    delete[] sigpi;
    return L;
}

void Matrix::invert()
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::invert: Matrix is non-totally symmetric.");
    }

    for (int h=0; h<nirrep_; ++h) {
        if (rowspi_[h]) {
            int err = C_DPOTRI('L', rowspi_[h], matrix_[h][0], rowspi_[h]);
            if (err != 0) {
                if (err < 0) {
                    fprintf(outfile, "invert: C_DPOTRI: argument %d has invalid paramter.\n", -err);
                    fflush(outfile);
                    abort();
                }
                if (err > 1) {
                    fprintf(outfile, "invert: C_DPOTRI: the (%d,%d) element of the factor U or L is "
                            "zero, and the inverse could not be computed.\n", err, err);
                    fflush(outfile);
                    abort();
                }
            }
        }
    }
    copy_lower_to_upper();
}

void Matrix::power(double alpha, double cutoff)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::power: Matrix is non-totally symmetric.");
    }

    for (int h=0; h<nirrep_; ++h) {
        if (rowspi_[h] == 0) continue;

        int n = rowspi_[h];
        double** A = matrix_[h];

        double** A1 = Matrix::matrix(n,n);
        double** A2 = Matrix::matrix(n,n);
        double* a  = new double[n];

        memcpy(static_cast<void*>(A1[0]), static_cast<void*>(A[0]), sizeof(double)*n*n);

        // Eigendecomposition
        double lwork;
        int stat = C_DSYEV('V','U',n,A1[0],n,a,&lwork,-1);
        double* work = new double[(int)lwork];
        stat = C_DSYEV('V','U',n,A1[0],n,a,work,(int)lwork);
        delete[] work;

        if (stat)
            throw PSIEXCEPTION("Matrix::power: C_DSYEV failed");

        memcpy(static_cast<void*>(A2[0]), static_cast<void*>(A1[0]), sizeof(double)*n*n);

        double max_a = (fabs(a[n-1]) > fabs(a[0]) ? fabs(a[n-1]) : fabs(a[0]));
        for (int i=0; i<n; i++) {

            if (alpha < 0.0 && fabs(a[i]) < cutoff * max_a)
                a[i] = 0.0;
            else
                a[i] = pow(a[i],alpha);

            C_DSCAL(n, a[i], A2[i], 1);
        }

        C_DGEMM('T','N',n,n,n,1.0,A2[0],n,A1[0],n,0.0,A[0],n);

        delete[] a;
        Matrix::free(A1);
        Matrix::free(A2);

    }
}

void Matrix::exp()
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::exp: Matrix is non-totally symmetric.");
    }

    for (int h=0; h<nirrep_; ++h) {
        if (rowspi_[h] == 0) continue;

        int n = rowspi_[h];
        double** A = matrix_[h];

        double** A1 = Matrix::matrix(n,n);
        double** A2 = Matrix::matrix(n,n);
        double* a  = new double[n];

        memcpy(static_cast<void*>(A1[0]), static_cast<void*>(A[0]), sizeof(double)*n*n);

        // Eigendecomposition
        double lwork;
        int stat = C_DSYEV('V','U',n,A1[0],n,a,&lwork,-1);
        double* work = new double[(int)lwork];
        stat = C_DSYEV('V','U',n,A1[0],n,a,work,(int)lwork);
        delete[] work;

        if (stat)
            throw PSIEXCEPTION("Matrix::exp: C_DSYEV failed");

        memcpy(static_cast<void*>(A2[0]), static_cast<void*>(A1[0]), sizeof(double)*n*n);

        for (int i=0; i<n; i++) {

            a[i] = ::exp(a[i]);

            C_DSCAL(n, a[i], A2[i], 1);
        }

        C_DGEMM('T','N',n,n,n,1.0,A2[0],n,A1[0],n,0.0,A[0],n);

        delete[] a;
        Matrix::free(A1);
        Matrix::free(A2);

    }
}

void Matrix::zero_lower()
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::zero_lower: Matrix is non-totally symmetric.");
    }

    for (int h=0; h<nirrep_; ++h) {
        for (int m=0; m<rowspi_[h]; ++m) {
            for (int n=0; n<m; ++n) {
                matrix_[h][m][n] = 0.0;
            }
        }
    }

}

void Matrix::zero_upper()
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::zero_upper: Matrix is non-totally symmetric.");
    }

    for (int h=0; h<nirrep_; ++h) {
        for (int m=0; m<rowspi_[h]; ++m) {
            for (int n=0; n<m; ++n) {
                matrix_[h][n][m] = 0.0;
            }
        }
    }

}

void Matrix::copy_lower_to_upper()
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::copy_lower_to_upper: Matrix is non-totally symmetric.");
    }

    for (int h=0; h<nirrep_; ++h) {
        for (int m=0; m<rowspi_[h]; ++m) {
            for (int n=0; n<m; ++n) {
                matrix_[h][n][m] = matrix_[h][m][n];
            }
        }
    }
}

void Matrix::copy_upper_to_lower()
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::copy_upper_to_lower: Matrix is non-totally symmetric.");
    }

    for (int h=0; h<nirrep_; ++h) {
        for (int m=0; m<rowspi_[h]; ++m) {
            for (int n=0; n<m; ++n) {
                matrix_[h][m][n] = matrix_[h][n][m];
            }
        }
    }
}

// Reference versions of the above functions:

void Matrix::transform(const Matrix& a, const Matrix& transformer)
{
    Matrix temp(a);

    temp.gemm(false, false, 1.0, a, transformer, 0.0);
    gemm(true, false, 1.0, transformer, temp, 0.0);
}

void Matrix::apply_symmetry(const boost::shared_ptr<Matrix>& a, const boost::shared_ptr<Matrix>& transformer)
{
    // Check dimensions of the two matrices and symmetry
    if(a->nirrep() > 1)
        {
            throw PSIEXCEPTION("Matrix::apply_symmetry: first matrix must have no symmetry.\n");
        }
    if (a->nrow() != transformer->rowdim(0)
            || a->ncol() != transformer->ncol()) {
        a->print();
        transformer->print();
        throw PSIEXCEPTION("Matrix::apply_symmetry: simple to regular. Sizes are not compatible.\n");
    }

    // Create temporary matrix of proper size.
    Matrix temp(nirrep(), a->nrow(), transformer->colspi());

    char ta = 'n';
    char tb = 'n';
    int h, m, n, k, nca, ncb, ncc;

    // Solve F = T^ M T

    // temp = M T
    for(int h=0; h<nirrep_; ++h) {
        m = temp.rowdim(h);
        n = temp.coldim(h);
        k = a->ncol();
        nca = k;
        ncb = n;
        ncc = n;

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, 1.0, &(a->matrix_[0][0][0]),
                    nca, &(transformer->matrix_[h][0][0]), ncb,
                    0.0, &(temp.matrix_[h][0][0]), ncc);
        }
    }

    // F = T^ temp
    ta = 't';
    for (int h=0; h<nirrep_; ++h) {
        m = rowdim(h);
        n = coldim(h);
        k = transformer->rowdim(h);
        nca = m;
        ncb = n;
        ncc = n;

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, 1.0, &(transformer->matrix_[h][0][0]),
                    nca, &(temp.matrix_[h][0][0]), ncb, 0.0, &(matrix_[h][0][0]),
                    ncc);
        }
    }
}

void Matrix::remove_symmetry(const boost::shared_ptr<Matrix>& a, const boost::shared_ptr<Matrix>& transformer)
{
    // Check dimensions of the two matrices and symmetry
    if(a->nirrep() !=  transformer->nirrep())
        {
            throw PSIEXCEPTION("Matrix::remove_symmetry: matrices must have same symmetry.\n");
        }
    if(nirrep() != 1)
        {
            throw PSIEXCEPTION("Matrix::remove_symmetry: result matrix must not have symmetry. \n");
        }
    if (ncol() != transformer->coldim(0)
            || a->nrow() != transformer->nrow()) {
        a->print();
        transformer->print();
        throw PSIEXCEPTION("Matrix::remove_symmetry: Sizes are not compatible.\n");
    }

    zero();

    // Create temporary matrix of proper size.
    Matrix temp(transformer->nirrep(), transformer->rowspi(), transformer->colspi());

    char ta = 'n';
    char tb = 'n';
    int h, m, n, k, nca, ncb, ncc;

    // Solve F = T^ M T

    // temp = M T
    for(int h=0; h<transformer->nirrep(); ++h) {
        m = temp.rowdim(h);
        n = temp.coldim(h);
        k = a->coldim(h);
        nca = k;
        ncb = n;
        ncc = n;

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, 1.0, &(a->matrix_[h][0][0]),
                    nca, &(transformer->matrix_[h][0][0]), ncb,
                    1.0, &(temp.matrix_[h][0][0]), ncc);
        }
    }

    // F = T^ temp
    ta = 't';
    for (int h=0; h<transformer->nirrep(); ++h) {
        m = nrow();
        n = ncol();
        k = temp.rowdim(h);
        nca = m; //k
        ncb = n; //k
        ncc = n; //k

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, 1.0, &(transformer->matrix_[h][0][0]),
                    nca, &(temp.matrix_[h][0][0]), ncb, 1.0, &(matrix_[0][0][0]),
                    ncc);
        }
    }
}

void Matrix::transform(const Matrix& transformer)
{
    Matrix temp(this);

    temp.gemm(false, false, 1.0, *this, transformer, 0.0);
    gemm(true, false, 1.0, transformer, temp, 0.0);
}

void Matrix::back_transform(const Matrix& a, const Matrix& transformer)
{
    Matrix temp(a);

    temp.gemm(false, true, 1.0, a, transformer, 0.0);
    gemm(false, false, 1.0, transformer, temp, 0.0);
}

void Matrix::back_transform(const Matrix& transformer)
{
    Matrix temp(this);

    temp.gemm(false, true, 1.0, *this, transformer, 0.0);
    gemm(false, false, 1.0, transformer, temp, 0.0);
}

double Matrix::vector_dot(const Matrix& rhs)
{
    return vector_dot(&rhs);
}

void Matrix::diagonalize(Matrix& eigvectors, Vector& eigvalues, int nMatz)
{
    int h;
    for (h=0; h<nirrep_; ++h) {
        if (rowspi_[h]) {
            sq_rsp(rowspi_[h], colspi_[h], matrix_[h], eigvalues.vector_[h], nMatz, eigvectors.matrix_[h], 1.0e-14);
        }
    }
}

void Matrix::write_to_dpdfile2(dpdfile2 *outFile)
{
    dpd_file2_mat_init(outFile);

    if(outFile->params->nirreps != nirrep_) {
        char *str = new char[100];
        sprintf(str, "Irrep count mismatch.  Matrix class has %d irreps, but dpdfile2 has %d.",
                nirrep_, outFile->params->nirreps);
        throw SanityCheckError(str, __FILE__, __LINE__);
    }

    if(outFile->my_irrep != 0) {
        throw FeatureNotImplemented("libmints Matrix class",
                "Matrices whose irrep is not the symmetric one", __FILE__, __LINE__);
    }

    for(int h = 0; h < nirrep_; ++h){
        if(outFile->params->rowtot[h] != rowspi_[h]){
            char *str = new char[100];
            sprintf(str, "Row count mismatch for irrep %d.  Matrix class has %d rows, but dpdfile2 has %d.",
                    h, rowspi_[h], outFile->params->rowtot[h]);
            throw SanityCheckError(str, __FILE__, __LINE__);
        }
        if(outFile->params->coltot[h] != colspi_[h]){
            char *str = new char[100];
            sprintf(str, "Column count mismatch for irrep %d.  Matrix class has %d columns, but dpdfile2 has %d.",
                    h, colspi_[h], outFile->params->coltot[h]);
            throw SanityCheckError(str, __FILE__, __LINE__);
        }

        // TODO: optimize this with memcopys
        for(int row = 0; row < rowspi_[h]; ++row){
            for(int col = 0; col < colspi_[h]; ++col){
                outFile->matrix[h][row][col] = matrix_[h][row][col];
            }
        }
    }

    dpd_file2_mat_wrt(outFile);
    dpd_file2_mat_close(outFile);
}

void Matrix::save(const string& filename, bool append, bool saveLowerTriangle, bool saveSubBlocks)
{
    static const char *str_block_format = "%3d %3d %3d %20.15f\n";
    static const char *str_full_format  = "%3d %3d %20.15f\n";

    // We can only save lower triangle if symmetry_ if 0
    if (symmetry_ && saveLowerTriangle)
        throw PSIEXCEPTION("Matrix::save: Unable to save lower triangle for non-totally symmetric matrix.");

    FILE *out = NULL;
    if (append == true) {
        out = fopen(filename.c_str(), "a");
    } else {
        out = fopen(filename.c_str(), "w");
    }

    fprintf(out, "%s\n", name_.c_str());
    fprintf(out, "symmetry %d\n", symmetry_);

    if (saveSubBlocks == false) {
        // Convert the matrix to a full matrix
        double **fullblock = to_block_matrix();

        // Need to know the size
        int sizer=0, sizec=0;
        for (int h=0; h<nirrep_; ++h) {
            sizer += rowspi_[h];
            sizec += colspi_[h^symmetry_];
        }

        if (saveLowerTriangle) {
            // Count the number of non-zero element
            int count=0;
            for (int i=0; i<sizer; ++i) {
                for (int j=0; j<=i; ++j) {
                    if (fabs(fullblock[i][j]) > 1.0e-14) {
                        count++;
                    }
                }
            }
            fprintf(out, "%5d\n", count);
            for (int i=0; i<sizer; ++i) {
                for (int j=0; j<=i; ++j) {
                    if (fabs(fullblock[i][j]) > 1.0e-14) {
                        fprintf(out, str_full_format, i, j, fullblock[i][j]);
                    }
                }
            }
        } else {
            // Count the number of non-zero element
            int count=0;
            for (int i=0; i<sizer; ++i) {
                for (int j=0; j<sizec; ++j) {
                    if (fabs(fullblock[i][j]) > 1.0e-14) {
                        count++;
                    }
                }
            }
            fprintf(out, "%5d\n", count);
            for (int i=0; i<sizer; ++i) {
                for (int j=0; j<sizec; ++j) {
                    if (fabs(fullblock[i][j]) > 1.0e-14) {
                        fprintf(out, str_full_format, i, j, fullblock[i][j]);
                    }
                }
            }
        }
        Matrix::free(fullblock);
    } else {
        if (saveLowerTriangle) {
            // Count the number of non-zero elements
            int count=0;
            for (int h=0; h<nirrep_; ++h) {
                for (int i=0; i<rowspi_[h]; ++i) {
                    for (int j=0; j<=i; ++j) {
                        if (fabs(matrix_[h][i][j]) > 1.0e-14) {
                            count++;
                        }
                    }
                }
            }
            fprintf(out, "%5d\n", count);
            for (int h=0; h<nirrep_; ++h) {
                for (int i=0; i<rowspi_[h]; ++i) {
                    for (int j=0; j<=i; ++j) {
                        if (fabs(matrix_[h][i][j]) > 1.0e-14) {
                            fprintf(out, str_block_format, h, i, j, matrix_[h][i][j]);
                        }
                    }
                }
            }
        } else {
            // Count the number of non-zero elements
            int count=0;
            for (int h=0; h<nirrep_; ++h) {
                for (int i=0; i<rowspi_[h]; ++i) {
                    for (int j=0; j<colspi_[h^symmetry_]; ++j) {
                        if (fabs(matrix_[h][i][j]) > 1.0e-14) {
                            count++;
                        }
                    }
                }
            }
            fprintf(out, "%5d\n", count);
            for (int h=0; h<nirrep_; ++h) {
                for (int i=0; i<rowspi_[h]; ++i) {
                    for (int j=0; j<colspi_[h^symmetry_]; ++j) {
                        if (fabs(matrix_[h][i][j]) > 1.0e-14) {
                            fprintf(out, str_block_format, h, i, j, matrix_[h][i][j]);
                        }
                    }
                }
            }
        }
    }

    fclose(out);
}

void Matrix::save(boost::shared_ptr<psi::PSIO>& psio,
                  unsigned int fileno,
                  SaveType savetype)
{
    save(psio.get(), fileno, savetype);
}

void Matrix::save(psi::PSIO* const psio, unsigned int fileno, SaveType st)
{
    if (symmetry_ && st == LowerTriangle)
        throw PSIEXCEPTION("Matrix::save: Unable to save triangle of non-totally symmetric matrix.");

    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) {
        already_open = true;
    } else {
        psio->open(fileno, PSIO_OPEN_OLD);
    }

    // Need to know the size
    int sizer=0, sizec=0;
    for (int h=0; h<nirrep_; ++h) {
        sizer += rowspi_[h];
        sizec += colspi_[h^symmetry_];
    }

    if (st == SubBlocks) {
        for (int h=0; h<nirrep_; ++h) {
            string str(name_);
            str += " Symmetry " + to_string(symmetry_) + " Irrep " + to_string(h);

            // Write the sub-blocks
            if (colspi_[h^symmetry_] > 0 && rowspi_[h] > 0)
                psio->write_entry(fileno, const_cast<char*>(str.c_str()), (char*)matrix_[h][0],
                                  sizeof(double) * colspi_[h^symmetry_] * rowspi_[h]);
        }
    }
    else if (st == Full) {
        double **fullblock = to_block_matrix();

        // Write the full block
        if (sizer > 0 && sizec > 0)
            psio->write_entry(fileno, const_cast<char*>(name_.c_str()), (char*)fullblock[0],
                              sizeof(double) * sizer * sizec);

        Matrix::free(fullblock);
    }
    else if (st == LowerTriangle) {
        double *lower = to_lower_triangle();

        if (sizer > 0)
            psio->write_entry(fileno, const_cast<char*>(name_.c_str()), (char*)lower, sizeof(double)*ioff[sizer]);
        delete[] lower;
    }
    else {
        throw PSIEXCEPTION("Matrix::save: Unknown SaveType\n");
    }
    if (!already_open)
        psio->close(fileno, 1);     // Close and keep
}

bool Matrix::load(boost::shared_ptr<psi::PSIO>& psio,
                  unsigned int fileno,
                  const std::string& tocentry,
                  int nso)
{
    return load(psio.get(), fileno, tocentry, nso);
}

bool Matrix::load(psi::PSIO* const psio, unsigned int fileno, const string& tocentry, int nso)
{
    if (symmetry_) {
        throw PSIEXCEPTION("Matrix::load: Matrix is non-totally symmetric.");
    }

    double *integrals = init_array(ioff[nso]);

    // If psi fails to read in the data this will abort out.
    if (!tocentry.empty())
        psi::IWL::read_one(psio, fileno, tocentry.c_str(), integrals, ioff[nso], 0, 0, outfile);
    else
        psi::IWL::read_one(psio, fileno, name_.c_str(), integrals, ioff[nso], 0, 0, outfile);

    set(integrals);

    ::free(integrals);

    return true;
}

void Matrix::load(psi::PSIO* const psio, unsigned int fileno, SaveType st)
{
    if (symmetry_ && st == LowerTriangle)
        throw PSIEXCEPTION("Matrix::load: Unable to read lower triangle of non-totally symmetric matrix.");

    // The matrix must be sized correctly first
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) {
        already_open = true;
    } else {
        psio->open(fileno, PSIO_OPEN_OLD);
    }

    // Need to know the size
    int sizer=0, sizec=0;
    for (int h=0; h<nirrep_; ++h) {
        sizer += rowspi_[h];
        sizec += colspi_[h^symmetry_];
    }

    if (st == SubBlocks) {
        for (int h=0; h<nirrep_; ++h) {
            string str(name_);
            str += " Symmetry " + to_string(symmetry_) + " Irrep " + to_string(h);

            // Read the sub-blocks
            if (colspi_[h] > 0 && rowspi_[h] > 0)
                psio->read_entry(fileno, str.c_str(), (char*)matrix_[h][0],
                                 sizeof(double) * colspi_[h^symmetry_] * rowspi_[h]);
        }
    }
    else if (st == Full) {
        double **fullblock = to_block_matrix();

        // Read the full block
        if (sizer > 0 && sizec > 0)
            psio->read_entry(fileno, name_.c_str(), (char*)fullblock[0], sizeof(double) * sizer * sizec);

        set(fullblock);
        Matrix::free(fullblock);
    }
    else if (st == LowerTriangle) {
        double *lower = to_lower_triangle();

        if (sizer > 0)
            psio->read_entry(fileno, name_.c_str(), (char*)lower, sizeof(double)*ioff[sizer]);
        set(lower);
        delete[] lower;
    }
    else {
        throw PSIEXCEPTION("Matrix::load: Unknown SaveType\n");
    }
    if (!already_open)
        psio->close(fileno, 1);     // Close and keep // Idempotent, win!
}

void Matrix::load(boost::shared_ptr<psi::PSIO>& psio,
                  unsigned int fileno,
                  SaveType savetype)
{
    load(psio.get(), fileno, savetype);
}

void Matrix::load(const std::string &filename)
{
    // Entire file.
    vector<string> lines;
    // temp
    string text;

    // Stream to use
    ifstream infile(filename.c_str());
    if (!infile)
        throw PSIEXCEPTION("Matrix::load: Unable to open file " + filename);

    // stupid way of reading in entire file
    while (infile.good()) {
        getline(infile, text);
        trim(text);
        if (!text.empty())
            lines.push_back(text);
    }

    // File MUST be at least 3 lines
    if (lines.size() < 3)
        throw PSIEXCEPTION("Matrix::load: " + filename + " insufficient length.");

    // First line is label (name of the matrix)
    set_name(lines[0]);

    // Second line is symmetry information
    smatch match;
    int infile_symm=-1;
    regex symmetry_line("^\\s*symmetry\\s*(\\d+)\\s*", regbase::icase);
    if (regex_match(lines[1], match, symmetry_line))
        infile_symm = str_to_int(match[1]);
    if (infile_symm != symmetry_) {
        release();
        symmetry_ = infile_symm;
        alloc();
    }

    // Third line is number of nonzero elements in the matrix.
    int infile_nonzero=0;
    regex nonzero_line("^\\s*(\\d+)\\s*");
    if (regex_match(lines[2], match, nonzero_line))
        infile_nonzero = str_to_int(match[1]);

    // Check the number of lines with the number of nonzero elements
    int nline = lines.size();
    if (nline-3 != infile_nonzero)
        throw PSIEXCEPTION("Matrix::load: Specified number of nonzero elements does not match number of elements given.");

    // Clear out this matrix
    zero();

    // Go through the file grabbing the data.
    regex element_line("^\\s*(\\d+)\\s*(\\d+)\\s*(\\d+)\\s*" NUMBER);
    int h, m, n;
    double val;
    for (int elem=0; elem<infile_nonzero; ++elem) {
        if (regex_match(lines[elem+3], match, element_line)) {
            h = str_to_int(match[1]);
            m = str_to_int(match[2]);
            n = str_to_int(match[3]);
            val = str_to_double(match[4]);
        }
        else
            throw PSIEXCEPTION("Matrix::load: Unable to match the following line:\n" + lines[elem+3]);

        // Check the info to make sure it is good
        if (h >= nirrep_)
            throw PSIEXCEPTION("Matrix::load: Irrep number is too large:\n" + lines[elem+3]);
        if (m >= rowspi_[h])
            throw PSIEXCEPTION("Matrix::load: Row number is too large:\n" + lines[elem+3]);
        if (n >= colspi_[h^symmetry_])
            throw PSIEXCEPTION("Matrix::load: Column number is too large:\n" + lines[elem+3]);

        // Set the data
        set(h, m, n, val);
    }
}

void Matrix::send()
{
}

void Matrix::recv()
{
}

void Matrix::bcast(int broadcaster)
{
    // Assume the user allocated the matrix to the correct size first.
    for (int h=0; h<nirrep_; ++h) {
        if (rowspi_[h] > 0 && colspi_[h] > 0)
            Communicator::world->bcast(matrix_[h][0], rowspi_[h] * colspi_[h^symmetry_], broadcaster);
    }
}

void Matrix::sum()
{
    for (int h=0; h<nirrep_; ++h)
        if (rowspi_[h] > 0 && colspi_[h] > 0)
            Communicator::world->sum(matrix_[h][0], rowspi_[h] * colspi_[h^symmetry_]);
}

bool Matrix::equal(const Matrix& rhs)
{
    return equal(&rhs);
}

bool Matrix::equal(const boost::shared_ptr<Matrix>& rhs)
{
    return equal(rhs.get());
}

bool Matrix::equal(const Matrix* rhs)
{
    // Check dimensions
    if (rhs->nirrep() != nirrep())
        return false;

    if (symmetry_ != rhs->symmetry_)
        return false;

    for (int h=0; h<nirrep(); ++h)
        if ((rowspi()[h] != rhs->rowspi()[h]) ||
                (colspi()[h] != rhs->colspi()[h]))
            return false;

    // Check element by element
    for (int h=0; h<nirrep(); ++h) {
        for (int m = 0; m < rowspi()[h]; ++m) {
            for (int n = 0; n < colspi()[h^symmetry_]; ++n) {
                if (get(h, m, n) != rhs->get(h, m, n))
                    return false;
            }
        }
    }

    return true;
}

double Matrix::pyget(const boost::python::tuple &key)
{
    return get(boost::python::extract<int>(key[0]),
               boost::python::extract<int>(key[1]),
               boost::python::extract<int>(key[2]));
}

void Matrix::pyset(const boost::python::tuple &key, double value)
{
    return set(boost::python::extract<int>(key[0]),
               boost::python::extract<int>(key[1]),
               boost::python::extract<int>(key[2]),
               value);
}

//
// SimpleMatrix
//
SimpleMatrix::SimpleMatrix() : matrix_(0), rows_(0), cols_(0)
{

}

SimpleMatrix::SimpleMatrix(string name) :
        matrix_(0), rows_(0), cols_(0), name_(name)
{

}

SimpleMatrix::SimpleMatrix(const SimpleMatrix& c) : matrix_(0)
{
    rows_ = c.rows_;
    cols_ = c.cols_;
    name_ = c.name_ + " Copy";
    alloc();
    copy_from(c.matrix_);
}

SimpleMatrix::SimpleMatrix(const SimpleMatrix* c) : matrix_(0)
{
    rows_ = c->rows_;
    cols_ = c->cols_;
    name_ = c->name_ + " Copy";
    alloc();
    copy_from(c->matrix_);
}

SimpleMatrix::SimpleMatrix(boost::shared_ptr<SimpleMatrix> c) : matrix_(0)
{
    rows_ = c->rows_;
    cols_ = c->cols_;
    name_ = c->name_ + " Copy";
    alloc();
    copy_from(c->matrix_);
}

SimpleMatrix::SimpleMatrix(int l_rows, int l_cols) : matrix_(0)
{
    rows_ = l_rows;
    cols_ = l_cols;
    alloc();
}

SimpleMatrix::SimpleMatrix(string name, int l_rows, int l_cols) : matrix_(0)
{
    rows_ = l_rows;
    cols_ = l_cols;
    name_ = name;
    alloc();
}

SimpleMatrix::SimpleMatrix(const Matrix* c) : matrix_(0), rows_(0), cols_(0)
{
    for (int i=0; i<c->nirrep(); ++i) {
        rows_ += c->rowspi()[i];
        cols_ += c->colspi()[i];
    }
    matrix_ = c->to_block_matrix();
}

SimpleMatrix::SimpleMatrix(const Matrix& c) : matrix_(0), rows_(0), cols_(0)
{
    for (int i=0; i<c.nirrep(); ++i) {
        rows_ += c.rowspi()[i];
        cols_ += c.colspi()[i];
    }
    matrix_ = c.to_block_matrix();
}

SimpleMatrix::SimpleMatrix(string& name, const Dimension& rows, const Dimension& cols)
{
    name_ = name;
    rows_ = rows[0];
    cols_ = cols[0];
    alloc();
}

SimpleMatrix::~SimpleMatrix()
{
    release();
}

double** SimpleMatrix::matrix(int nrow, int ncol)
{
    double** mat = (double**) malloc(sizeof(double*)*nrow);
    const size_t size = sizeof(double)*nrow*ncol;
    mat[0] = (double*) malloc(size);
    //bzero((void*)mat[0],size);
    memset((void *)mat[0], '\0', size);
    for(int r=1; r<nrow; ++r) mat[r] = mat[r-1] + ncol;
    return mat;
}

void SimpleMatrix::free(double** Block)
{
    ::free(Block[0]);  ::free(Block);
}

void SimpleMatrix::init(int rowspi, int colspi, string name)
{
    rows_ = rowspi;
    cols_ = colspi;
    name_ = name;
    alloc();
}

SimpleMatrix* SimpleMatrix::clone() const
{
    SimpleMatrix *temp = new SimpleMatrix(this);
    return temp;
}

void SimpleMatrix::copy(SimpleMatrix* cp)
{
    // Make sure we are the same size
    bool same = true;
    if (rows_ != cp->rows_ || cols_ != cp->cols_)
        same = false;

    if (same == false) {
        release();
        rows_ = cp->rows_;
        cols_ = cp->cols_;
        alloc();
    }

    memcpy(&(matrix_[0][0]), &(cp->matrix_[0][0]), rows_ * cols_ * sizeof(double));
}

void SimpleMatrix::copy(boost::shared_ptr<SimpleMatrix> cp)
{
    copy(cp.get());
}

void SimpleMatrix::alloc()
{
    if (matrix_)
        release();

    matrix_ = SimpleMatrix::matrix(rows_, cols_);
}

void SimpleMatrix::release()
{
    if (!matrix_)
        return;

    SimpleMatrix::free(matrix_);
    matrix_ = NULL;
}

void SimpleMatrix::copy_from(double **c)
{
    size_t size = rows_ * cols_ * sizeof(double);
    if (size)
        memcpy(&(matrix_[0][0]), &(c[0][0]), size);
}

void SimpleMatrix::set(double val)
{
    size_t size = rows_ * cols_;

    for (size_t i=0; i<size; ++i) {
        matrix_[0][i] = val;
    }
}

void SimpleMatrix::set(const double *tri)
{
    int ij = 0;
    for (int i=0; i<rows_; ++i) {
        for (int j=0; j<=cols_; ++j) {
            matrix_[i][j] = matrix_[j][i] = tri[ij++];
        }
    }
}

void SimpleMatrix::set(double **mat)
{
    copy_from(mat);
}

void SimpleMatrix::set(SimpleVector* vec)
{
    zero();
    for (int i=0; i<rows_; ++i)
        matrix_[i][i] = vec->vector_[i];
}

void SimpleMatrix::set(boost::shared_ptr<SimpleVector> vec)
{
    set(vec.get());
}

double ** SimpleMatrix::to_block_matrix() const
{
    double **temp = SimpleMatrix::matrix(rows_, cols_);
    memcpy(&(temp[0][0]), &(matrix_[0][0]), rows_*cols_*sizeof(double));
    return temp;
}

void SimpleMatrix::print(FILE *out)
{
    if (name_.length()) {
        fprintf(out, "  ## %s ##\n", name_.c_str());
    }

    print_mat(matrix_, rows_, cols_, out);
    fprintf(out, "\n");
}

void SimpleMatrix::print_atom_vector(FILE *out)
{
    int i;

    if (name_.length()) {
        fprintf(out,"\n  -%s:\n", name_.c_str());
    }
    fprintf(out,"     Atom            X                  Y                   Z\n");
    fprintf(out,"    ------   -----------------  -----------------  -----------------\n");

    for(i=0;i<nrow();i++) {
        fprintf(out,"    %4d   ",i+1);
        fprintf(out,"  %17.12lf  %17.12lf  %17.12lf", matrix_[i][0], matrix_[i][1], matrix_[i][2]);
        fprintf(out,"\n");
    }
    fprintf(out,"\n");
}

void SimpleMatrix::eivprint(SimpleVector *values, FILE *out)
{
    if (name_.length()) {
        fprintf(out, "  ## %s with eigenvalues ##\n", name_.c_str());
    }

    eivout(matrix_, values->vector_, rows_, cols_, out);
    fprintf(out, "\n");
}

void SimpleMatrix::eivprint(boost::shared_ptr<SimpleVector> values, FILE *out)
{
    eivprint(values.get(), out);
}

void SimpleMatrix::identity()
{
    zero();
    for (int i=0; i<MIN(rows_, cols_); ++i)
        matrix_[i][i] = 1.0;
}

void SimpleMatrix::zero()
{
    memset(&(matrix_[0][0]), 0, rows_*cols_*sizeof(double));
}

void SimpleMatrix::zero_diagonal()
{
    for (int i=0; i<MIN(rows_, cols_); ++i)
        matrix_[i][i] = 0.0;
}

double SimpleMatrix::trace() const
{
    double val = 0.0;
    for (int i=0; i<MIN(rows_, cols_); ++i) {
        val += matrix_[i][i];
    }
    return val;
}

void SimpleMatrix::lu_factorize()
{
    if (rows_) {
        int *ipiv = new int[rows_];
        int err = C_DGETRF(rows_, cols_, matrix_[0], rows_, ipiv);
        if (err != 0) {
            if (err < 0) {
                fprintf(outfile, "cholesky_factorize: C_DPOTRF: argument %d has invalid paramter.\n", -err);
                fflush(outfile);
                abort();
            }
            if (err > 1) {
                fprintf(outfile, "cholesky_factorize: C_DPOTRF: the leading minor of order %d is not "
                        "positive definite, and the factorization could not be "
                        "completed.", err);
                fflush(outfile);
                abort();
            }
        }
        delete[] ipiv;
    }
}

void SimpleMatrix::invert()
{
    if (rows_) {
        double **work = block_matrix(rows_, cols_);
        invert_matrix(matrix_, work, rows_, outfile);
        copy_from(work);
        free_block(work);
//        double *work = new double[rows_*cols_];
//        int *ipiv = new int[rows_];
//        int err = C_DGETRI(rows_, matrix_[0], rows_, ipiv, work, rows_*cols_);
//        if (err != 0) {
//            if (err < 0) {
//                fprintf(outfile, "invert: C_DPOTRI: argument %d has invalid paramter.\n", -err);
//                fflush(outfile);
//                abort();
//            }
//            if (err > 1) {
//                fprintf(outfile, "invert: C_DPOTRI: the (%d,%d) element of the factor U or L is "
//                        "zero, and the inverse could not be computed.\n", err, err);
//                fflush(outfile);
//                abort();
//            }
//        }
//        delete[] work;
//        delete[] ipiv;
    }
}

void SimpleMatrix::copy_lower_to_upper()
{
    for (int m=0; m<rows_; ++m) {
        for (int n=0; n<m; ++n) {
            matrix_[m][n] = matrix_[n][m];
        }
    }
}

void SimpleMatrix::copy_upper_to_lower()
{
    for (int m=0; m<rows_; ++m) {
        for (int n=0; n<m; ++n) {
            matrix_[n][m] = matrix_[m][n];
        }
    }
}

void SimpleMatrix::transpose_this()
{
    double temp;

    for (int i=0; i<rows_; ++i) {
        for (int j=0; j<i; ++j) {
            temp = matrix_[i][j];
            matrix_[i][j] = matrix_[j][i];
            matrix_[j][i] = temp;
        }
    }
}

SimpleMatrix* SimpleMatrix::transpose()
{
    SimpleMatrix* temp = new SimpleMatrix(this);

    for (int i=0; i<rows_; ++i) {
        for (int j=0; j<cols_; ++j) {
            temp->matrix_[i][j] = matrix_[j][i];
        }
    }

    return temp;
}

void SimpleMatrix::add(const SimpleMatrix* plus)
{
    double *lhs, *rhs;
    size_t size = rows_ * cols_;
    if (size) {
        lhs = matrix_[0];
        rhs = plus->matrix_[0];
        for (size_t ij=0; ij<size; ++ij) {
            *lhs += *rhs;
            lhs++; rhs++;
        }
    }
}

void SimpleMatrix::add(boost::shared_ptr<SimpleMatrix> plus)
{
    add(plus.get());
}

void SimpleMatrix::element_add_mirror()
{
    for (int i=0; i<rows_; ++i) {
        for (int j=0; j<i; ++j) {
            matrix_[i][j] = matrix_[j][i] = (matrix_[i][j] + matrix_[j][i]);
        }
    }
}

void SimpleMatrix::subtract(const SimpleMatrix* plus)
{
    double *lhs, *rhs;
    size_t size = rows_ * cols_;
    if (size) {
        lhs = matrix_[0];
        rhs = plus->matrix_[0];
        for (size_t ij=0; ij<size; ++ij) {
            *lhs -= *rhs;
            lhs++; rhs++;
        }
    }
}

void SimpleMatrix::subtract(boost::shared_ptr<SimpleMatrix> sub)
{
    subtract(sub.get());
}

void SimpleMatrix::accumulate_product(const SimpleMatrix* a, const SimpleMatrix* b)
{
    gemm(false, false, 1.0, a, b, 1.0);
}

void SimpleMatrix::accumulate_product(boost::shared_ptr<SimpleMatrix> a, boost::shared_ptr<SimpleMatrix> b)
{
   gemm(false, false, 1.0, a, b, 1.0);
}

void SimpleMatrix::scale(double a)
{
    size_t size;
    size = rows_ * cols_;
    if (size)
        C_DSCAL(size, a, &(matrix_[0][0]), 1);
}

void SimpleMatrix::scale_row(int m, double a)
{
    C_DSCAL(cols_, a, &(matrix_[m][0]), 1);
}

void SimpleMatrix::scale_column(int n, double a)
{
    C_DSCAL(rows_, a, &(matrix_[0][n]), cols_);
}

double SimpleMatrix::sum_of_squares()
{
    double sum = (double)0.0;
    for (int i=0; i<rows_; ++i) {
        for (int j=0; j<cols_; ++j) {
            sum += matrix_[i][j] * matrix_[i][j];
        }
    }

    return sum;
}

void SimpleMatrix::transform(SimpleMatrix* a, SimpleMatrix* transformer)
{
    SimpleMatrix temp(a);

    temp.gemm(false, false, 1.0, a, transformer, 0.0);
    gemm(true, false, 1.0, transformer, &temp, 0.0);
}

void SimpleMatrix::transform(boost::shared_ptr<SimpleMatrix> a, boost::shared_ptr<SimpleMatrix> transformer)
{
    transform(a.get(), transformer.get());
}

void SimpleMatrix::transform(SimpleMatrix* transformer)
{
    SimpleMatrix temp(this);

    temp.gemm(false, false, 1.0, this, transformer, 0.0);
    gemm(true, false, 1.0, transformer, &temp, 0.0);
}

void SimpleMatrix::transform(boost::shared_ptr<SimpleMatrix> transformer)
{
    transform(transformer.get());
}

void SimpleMatrix::back_transform(SimpleMatrix* a, SimpleMatrix* transformer)
{
    SimpleMatrix temp(a);

    temp.gemm(false, true, 1.0, a, transformer, 0.0);
    gemm(false, false, 1.0, transformer, &temp, 0.0);
}

void SimpleMatrix::back_transform(boost::shared_ptr<SimpleMatrix> a, boost::shared_ptr<SimpleMatrix> transformer)
{
     back_transform(a.get(), transformer.get());
}

void SimpleMatrix::back_transform(SimpleMatrix* transformer)
{
    SimpleMatrix temp(this);

    temp.gemm(false, true, 1.0, this, transformer, 0.0);
    gemm(false, false, 1.0, transformer, &temp, 0.0);
}

void SimpleMatrix::back_transform(boost::shared_ptr<SimpleMatrix> transformer)
{
    back_transform(transformer.get());
}

void SimpleMatrix::gemm(bool transa, bool transb, double alpha, const SimpleMatrix* a, const SimpleMatrix* b, double beta)
{
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int m, n, k, nca, ncb, ncc;

    m = rows_;
    n = cols_;
    k = transa ? a->rows_ : a->cols_;
    nca = transa ? m : k;
    ncb = transb ? k : n;
    ncc = n;

    if (m && n && k) {
        C_DGEMM(ta, tb, m, n, k, alpha, &(a->matrix_[0][0]),
                nca, &(b->matrix_[0][0]), ncb, beta, &(matrix_[0][0]),
                ncc);
    }
}

void SimpleMatrix::gemm(bool transa, bool transb, double alpha, boost::shared_ptr<SimpleMatrix> a, boost::shared_ptr<SimpleMatrix> b, double beta)
{
    gemm(transa, transb, alpha, a.get(), b.get(), beta);
}

double SimpleMatrix::vector_dot(SimpleMatrix* rhs)
{
    double sum = 0.0;
    size_t size;

    size = rows_ * cols_;
    if (size)
        sum += C_DDOT(size, (&matrix_[0][0]), 1, &(rhs->matrix_[0][0]), 1);

    return sum;
}

double SimpleMatrix::vector_dot(boost::shared_ptr<SimpleMatrix> rhs)
{
    return vector_dot(rhs.get());
}

void SimpleMatrix::diagonalize(SimpleMatrix* eigvectors, SimpleVector* eigvalues, int sort)
{
    if (rows_) {
        sq_rsp(rows_, cols_, matrix_, eigvalues->vector_, sort, eigvectors->matrix_, 1.0e-14);
    }
}

void SimpleMatrix::diagonalize(boost::shared_ptr<SimpleMatrix> eigvectors, boost::shared_ptr<SimpleVector> eigvalues, int sort)
{
    diagonalize(eigvectors.get(), eigvalues.get(), sort);
}

void SimpleMatrix::save(psi::PSIO* psio, unsigned int fileno)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) {
        already_open = true;
    } else {
        psio->open(fileno, PSIO_OPEN_OLD);
    }

    psio->write_entry(fileno, const_cast<char*>(name_.c_str()), (char*)matrix_[0], sizeof(double) * rows_ * cols_);

    if (!already_open)
        psio->close(fileno, 1);     // Close and keep
}

void SimpleMatrix::load(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno)
{
    // Check to see if the file is open
    bool already_open = false;
    if (psio->open_check(fileno)) {
        already_open = true;
    } else {
        psio->open(fileno, PSIO_OPEN_OLD);
    }

    psio->read_entry(fileno, const_cast<char*>(name_.c_str()), (char*)matrix_[0], sizeof(double) * rows_ * cols_);

    if (!already_open)
        psio->close(fileno, 1);     // Close and keep
}

void SimpleMatrix::save(psi::PSIO& psio, unsigned int fileno)
{
    save(&psio, fileno);
}

void SimpleMatrix::save(boost::shared_ptr<psi::PSIO> psio, unsigned int fileno)
{
    save(psio.get(), fileno);
}

void SimpleMatrix::save(const char *filename, bool append, bool saveLowerTriangle)
{
    static const char *str_full_format  = "%3d %3d %20.15f\n";

    FILE *out = NULL;
    if (append == true) {
        out = fopen(filename, "a");
    } else {
        out = fopen(filename, "w");
    }

    fprintf(out, "%s\n", name_.c_str());

    if (saveLowerTriangle) {
        // Count the number of non-zero element
        int count=0;
        for (int i=0; i<rows_; ++i) {
            for (int j=0; j<=i; ++j) {
                if (fabs(matrix_[i][j]) > 1.0e-14) {
                    count++;
                }
            }
        }
        fprintf(out, "%5d\n", count);
        for (int i=0; i<rows_; ++i) {
            for (int j=0; j<=i; ++j) {
                if (fabs(matrix_[i][j]) > 1.0e-14) {
                    fprintf(out, str_full_format, i, j, matrix_[i][j]);
                }
            }
        }
    } else {
        // Count the number of non-zero element
        int count=0;
        for (int i=0; i<rows_; ++i) {
            for (int j=0; j<cols_; ++j) {
                if (fabs(matrix_[i][j]) > 1.0e-14) {
                    count++;
                }
            }
        }
        fprintf(out, "%5d\n", count);
        for (int i=0; i<rows_; ++i) {
            for (int j=0; j<cols_; ++j) {
                if (fabs(matrix_[i][j]) > 1.0e-14) {
                    fprintf(out, str_full_format, i, j, matrix_[i][j]);
                }
            }
        }
    }

    fclose(out);
}

void SimpleMatrix::swap_rows(int i, int j)
{
    C_DSWAP(cols_, &(matrix_[i][0]), 1, &(matrix_[j][0]), 1);
}

void SimpleMatrix::swap_columns(int i, int j)
{
    C_DSWAP(rows_, &(matrix_[0][i]), cols_, &(matrix_[0][j]), cols_);
}

