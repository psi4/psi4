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
#include "factory.h"
#include "wavefunction.h"
#include "dimension.h"

#include <libdpd/dpd.h>

#include <cmath>
#include <sstream>

 using namespace psi;

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

std::string to_string(const int val)
{
    std::stringstream strm;
    strm <<  val;
    return strm.str();
}

Matrix::Matrix()
{
    matrix_ = NULL;
    rowspi_ = NULL;
    colspi_ = NULL;
    nirreps_ = 0;
}

Matrix::Matrix(std::string name)
{
    matrix_ = NULL;
    rowspi_ = NULL;
    colspi_ = NULL;
    nirreps_ = 0;
    name = name_;
}

Matrix::Matrix(const Matrix& c)
{
    matrix_ = NULL;
    nirreps_ = c.nirreps_;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
    for (int i=0; i<nirreps_; ++i) {
        rowspi_[i] = c.rowspi_[i];
        colspi_[i] = c.colspi_[i];
    }
    alloc();
    copy_from(c.matrix_);
}

Matrix::Matrix(shared_ptr<Matrix> c)
{
    matrix_ = NULL;
    nirreps_ = c->nirreps_;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
    for (int i=0; i<nirreps_; ++i) {
        rowspi_[i] = c->rowspi_[i];
        colspi_[i] = c->colspi_[i];
    }
    alloc();
    copy_from(c->matrix_);
}

Matrix::Matrix(Matrix& c)
{
    matrix_ = NULL;
    nirreps_ = c.nirreps_;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
    for (int i=0; i<nirreps_; ++i) {
        rowspi_[i] = c.rowspi_[i];
        colspi_[i] = c.colspi_[i];
    }
    alloc();
    copy_from(c.matrix_);
}

Matrix::Matrix(const Matrix* c)
{
    matrix_ = NULL;
    nirreps_ = c->nirreps_;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
    for (int i=0; i<nirreps_; ++i) {
        rowspi_[i] = c->rowspi_[i];
        colspi_[i] = c->colspi_[i];
    }
    alloc();
    copy_from(c->matrix_);
}

Matrix::Matrix(int l_nirreps, int *l_rowspi, int *l_colspi)
{
    matrix_ = NULL;
    nirreps_ = l_nirreps;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
    for (int i=0; i<nirreps_; ++i) {
        rowspi_[i] = l_rowspi[i];
        colspi_[i] = l_colspi[i];
    }
    alloc();
}

Matrix::Matrix(std::string name, int l_nirreps, int *l_rowspi, int *l_colspi) : name_(name)
{
    matrix_ = NULL;
    nirreps_ = l_nirreps;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
    for (int i=0; i<nirreps_; ++i) {
        rowspi_[i] = l_rowspi[i];
        colspi_[i] = l_colspi[i];
    }
    alloc();
}

Matrix::Matrix(const std::string& name, const Dimension& rows, const Dimension& cols)
{
    matrix_ = NULL;

    // This will happen in PetiteList::aotoso()
    if (rows.n() == 1) {
        nirreps_ = cols.n();
        rowspi_ = new int[nirreps_];
        colspi_ = new int[nirreps_];
        for (int i=0; i<nirreps_; ++i) {
            rowspi_[i] = rows[0];
            colspi_[i] = cols[i];
        }
    }
    else {
        nirreps_ = rows.n();
        rowspi_ = new int[nirreps_];
        colspi_ = new int[nirreps_];
        for (int i=0; i<nirreps_; ++i) {
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
    nirreps_ = inFile->params->nirreps;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
    for (int i=0; i<nirreps_; ++i) {
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

void Matrix::init(int l_nirreps, int *l_rowspi, int *l_colspi, std::string name)
{
    if (rowspi_) delete[] rowspi_;
    if (colspi_) delete[] colspi_;
    name_ = name;
    nirreps_ = l_nirreps;
    rowspi_ = new int[nirreps_];
    colspi_ = new int[nirreps_];
    for (int i=0; i<nirreps_; ++i) {
        rowspi_[i] = l_rowspi[i];
        colspi_[i] = l_colspi[i];
    }
    alloc();
}

Matrix* Matrix::clone() const
{
    Matrix *temp = new Matrix(this);
    return temp;
}

void Matrix::copy(Matrix* cp)
{
    // Make sure we are the same size as cp
    bool same = true;
    if (nirreps_ != cp->nirreps_) {
        same = false;
    } else {
        for (int h=0; h<nirreps_; ++h)
            if (colspi_[h] != cp->colspi_[h] || rowspi_[h] != cp->colspi_[h])
                same = false;
    }

    if (same == false) {
        release();
        if (rowspi_)
            delete[] rowspi_;
        if (colspi_)
            delete[] colspi_;
        nirreps_ = cp->nirreps_;
        rowspi_ = new int[nirreps_];
        colspi_ = new int[nirreps_];
        for (int i=0; i<nirreps_; ++i) {
            rowspi_[i] = cp->rowspi_[i];
            colspi_[i] = cp->colspi_[i];
        }
        alloc();
    }

    // When here we are the same size
    for (int h=0; h<nirreps_; ++h) {
        if (rowspi_[h] != 0 && colspi_[h] != 0)
            memcpy(&(matrix_[h][0][0]), &(cp->matrix_[h][0][0]), rowspi_[h] * colspi_[h] * sizeof(double));
    }
}

void Matrix::copy(Matrix& cp)
{
    copy(&cp);
}

void Matrix::copy(shared_ptr<Matrix> cp)
{
    copy(cp.get());
}

void Matrix::copy(const Matrix& cp)
{
    copy(const_cast<Matrix&>(cp));
}

void Matrix::copy(const Matrix* cp)
{
    copy(const_cast<Matrix*>(cp));
}

void Matrix::alloc()
{
    if (matrix_)
        release();

    matrix_ = (double***)malloc(sizeof(double***) * nirreps_);
    for (int i=0; i<nirreps_; ++i) {
        if (rowspi_[i] != 0 && colspi_[i] != 0)
            matrix_[i] = Matrix::matrix(rowspi_[i], colspi_[i]);
        else
            matrix_[i] = NULL;
    }
}

void Matrix::release()
{
    if (!matrix_)
        return;

    for (int h=0; h<nirreps_; ++h) {
        if (matrix_[h])
            Matrix::free(matrix_[h]);
    }
    ::free(matrix_);
    matrix_ = NULL;
}

void Matrix::copy_from(double ***c) {
    int size;

    for (int h=0; h<nirreps_; ++h) {
        size = rowspi_[h] * colspi_[h] * sizeof(double);
        if (size)
            memcpy(&(matrix_[h][0][0]), &(c[h][0][0]), size);
    }
}

// Sets all elements of matrix to val
void Matrix::set(double val)
{
    for (int h=0; h < nirreps_; ++h) {
        size_t size = rowspi_[h] * colspi_[h];

        for (size_t i=0; i<size; ++i) {
            matrix_[h][0][i] = val;
        }
    }
}

void Matrix::set(const double *tri)
{
    int h, i, j, ii, jj;
    int offset;

    offset = 0;
    for (h=0; h<nirreps_; ++h) {
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

void Matrix::set(const double **sq)
{
    int h, i, j, ii, jj;
    int offset;

    if (sq == NULL) {
        zero();
        // TODO: Need to throw an exception here.
        return;
    }
    offset = 0;
    for (h=0; h<nirreps_; ++h) {
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

void Matrix::set(double **sq)
{
    int h, i, j, ii, jj;
    int offset;

    if (sq == NULL) {
        zero();
        // TODO: Need to throw an exception here.
        return;
    }
    offset = 0;
    for (h=0; h<nirreps_; ++h) {
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

void Matrix::set(SimpleMatrix *sq)
{
    set(const_cast<const double**>(sq->matrix_));
}

void Matrix::set(shared_ptr<SimpleMatrix> sq)
{
    set(const_cast<const double**>(sq.get()->matrix_));
}

void Matrix::set(Vector* vec)
{
    int h, i, size;
    zero();
    for (h=0; h<nirreps_; ++h) {
        size = rowspi_[h];
        if (size) {
            for (i=0; i<size; ++i)
                matrix_[h][i][i] = vec->vector_[h][i];
        }
    }
}

void Matrix::set(Vector& vec)
{
    set(&vec);
}

void Matrix::set(shared_ptr<Vector> vec)
{
    set(vec.get());
}

double *Matrix::to_lower_triangle() const
{
    int sizer=0, sizec=0;
    for (int h=0; h<nirreps_; ++h) {
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
    for (int h=0; h<nirreps_; ++h) {
        sizer += rowspi_[h];
        sizec += colspi_[h];
    }

    double **temp = block_matrix(sizer,sizec);
    int offsetr = 0, offsetc=0;
    for (int h=0; h <nirreps_; ++h) {
        for (int i=0; i<rowspi_[h]; ++i) {
            for (int j=0; j<colspi_[h]; ++j) {
                temp[i+offsetr][j+offsetc] = matrix_[h][i][j];
            }
        }
        offsetr += rowspi_[h];
        offsetc += colspi_[h];
    }

    return temp;
}

SimpleMatrix *Matrix::to_simple_matrix()
{
    return new SimpleMatrix(this);
}

void Matrix::print_mat(double **a, int m, int n, FILE *out)
{
    #if 0
    int ii,jj,kk,nn,ll;
    int i,j,k;

    ii=0;jj=0;
L200:
    ii++;
    jj++;
    kk=10*jj;
    nn=n;
    if (nn > kk) nn=kk;
    ll = 2*(nn-ii+1)+1;
    fprintf (out,"\n");
    for (i=ii; i <= nn; i++) fprintf(out,"       %5d",i);
    fprintf (out,"\n");
    for (i=0; i < m; i++) {
       fprintf (out,"\n%5d",i+1);
       for (j=ii-1; j < nn; j++) {
          fprintf (out,"%20.15f",a[i][j]);
          }
       }
    fprintf (out,"\n");
    if (n <= kk) {
       fflush(out);
       return;
       }
    ii=kk; goto L200;
    #endif
    ::print_mat(a, m, n, out);
}

void Matrix::print(FILE *out, const char *extra)
{
    int h;

    if (name_.length()) {
        if (extra == NULL)
            fprintf(out, "  ## %s ##\n", name_.c_str());
        else
            fprintf(out, "  ## %s %s ##\n", name_.c_str(), extra);
    }

    for (h=0; h<nirreps_; ++h) {
        fprintf(out, "  Irrep: %d\n", h+1);
        print_mat(matrix_[h], rowspi_[h], colspi_[h], out);
        fprintf(out, "\n");
    }
    fflush(out);
}

void Matrix::eivprint(Vector *values, FILE *out)
{
    int h;

    if (name_.length()) {
        fprintf(out, "  ## %s with eigenvalues ##\n", name_.c_str());
    }

    for (h=0; h<nirreps_; ++h) {
        fprintf(out, " Irrep: %d\n", h+1);
        eivout(matrix_[h], values->vector_[h], rowspi_[h], colspi_[h], out);
        fprintf(out, "\n");
    }
    fflush(out);
}

void Matrix::eivprint(Vector& values, FILE *out)
{
    eivprint(&values, out);
}

void Matrix::eivprint(shared_ptr<Vector> values, FILE *out)
{
    eivprint(values.get(), out);
}

void Matrix::identity()
{
    int h;
    size_t size;

    for (h=0; h<nirreps_; ++h) {
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

    for (h=0; h<nirreps_; ++h) {
        size = rowspi_[h] * colspi_[h] * sizeof(double);

        if (size) {
            memset(&(matrix_[h][0][0]), 0, size);
        }
    }
}

void Matrix::zero_diagonal()
{
    int h, i;

    for (h=0; h<nirreps_; ++h) {
        for (i=0; i<MIN(rowspi_[h], colspi_[h]); ++i) {
            matrix_[h][i][i] = 0.0;
        }
    }
}

double Matrix::trace()
{
    int i, h;
    double val = (double)0.0;

    for (h=0; h<nirreps_; ++h) {
        for (i=0; i<MIN(rowspi_[h], colspi_[h]); ++i) {
            val += matrix_[h][i][i];
        }
    }

    return val;
}

Matrix* Matrix::transpose()
{
    Matrix *temp = new Matrix(this);

    int h, i, j;
    for (h=0; h<nirreps_; ++h) {
        for (i=0; i<rowspi_[h]; ++i) {
            for (j=0; j<colspi_[h]; ++j) {
                temp->matrix_[h][i][j] = matrix_[h][j][i];
            }
        }
    }
    return temp;
}

void Matrix::add(const Matrix* plus)
{
    double *lhs, *rhs;
    for (int h=0; h<nirreps_; ++h) {
        size_t size = rowspi_[h] * colspi_[h];
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

void Matrix::add(shared_ptr<Matrix> plus)
{
    add(plus.get());
}

void Matrix::subtract(const Matrix* plus)
{
    double *lhs, *rhs;
    for (int h=0; h<nirreps_; ++h) {
        size_t size = rowspi_[h] * colspi_[h];
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

void Matrix::subtract(shared_ptr<Matrix> sub)
{
    subtract(sub.get());
}

void Matrix::accumulate_product(const Matrix* a, const Matrix* b)
{
    gemm(false, false, 1.0, a, b, 1.0);
}

void Matrix::accumulate_product(shared_ptr<Matrix> a, shared_ptr<Matrix> b)
{
    gemm(false, false, 1.0, a, b, 1.0);
}

void Matrix::scale(double a)
{
    int h;
    size_t size;
    for (h=0; h<nirreps_; ++h) {
        size = rowspi_[h] * colspi_[h];
        if (size)
            C_DSCAL(size, a, &(matrix_[h][0][0]), 1);
    }
}

void Matrix::scale_row(int h, int m, double a)
{
    C_DSCAL(rowspi_[h], a, &(matrix_[h][m][0]), 1);
}

void Matrix::scale_column(int h, int n, double a)
{
    C_DSCAL(colspi_[h], a, &(matrix_[h][0][n]), rowspi_[h]);
}

double Matrix::sum_of_squares()
{
    double sum = (double)0.0;
    for (int h=0; h<nirreps_; ++h) {
        for (int i=0; i<rowspi_[h]; ++i) {
            for (int j=0; j<colspi_[h]; ++j) {
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
    for (int h=0; h<nirreps_; ++h) {
        for (int i=0; i<rowspi_[h]; ++i) {
            for (int j=0; j<colspi_[h]; ++j) {
                sum += matrix_[h][i][j] * matrix_[h][i][j];
                terms++;
            }
        }
    }

    return sqrt(sum/terms);
}

void Matrix::transform(Matrix* a, Matrix* transformer)
{
    Matrix temp(a);

    temp.gemm(false, false, 1.0, a, transformer, 0.0);
    gemm(true, false, 1.0, transformer, &temp, 0.0);
}

void Matrix::transform(shared_ptr<Matrix> a, shared_ptr<Matrix> transformer)
{
    transform(a.get(), transformer.get());
}

void Matrix::transform(Matrix* transformer)
{
    bool square = true;
    int h = 0;

    while(h < nirreps_ && square){
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
        Matrix temp(nirreps_, rowspi_, transformer->colspi());
        Matrix result(nirreps_, transformer->colspi(), transformer->colspi());
        temp.gemm(false, false, 1.0, this, transformer, 0.0);
        result.gemm(true, false, 1.0, transformer, &temp, 0.0);
        copy(&result);
    }
}

void Matrix::transform(shared_ptr<Matrix> transformer)
{
    transform(transformer.get());
}

void Matrix::back_transform(Matrix* a, Matrix* transformer)
{
    Matrix temp(a);

    temp.gemm(false, true, 1.0, a, transformer, 0.0);
    gemm(false, false, 1.0, transformer, &temp, 0.0);
}

void Matrix::back_transform(shared_ptr<Matrix> a, shared_ptr<Matrix> transformer)
{
    back_transform(a.get(), transformer.get());
}

void Matrix::back_transform(Matrix* transformer)
{
    bool square = true;
    int h = 0;

    while(h < nirreps_ && square){
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
        Matrix temp(nirreps_, rowspi_, transformer->rowspi());
        Matrix result(nirreps_, transformer->rowspi(), transformer->rowspi());
        temp.gemm(false, true, 1.0, this, transformer, 0.0);
        result.gemm(false, false, 1.0, transformer, &temp, 0.0);
        copy(&result);
    }
}

void Matrix::back_transform(shared_ptr<Matrix> transformer)
{
    back_transform(transformer.get());
}

void Matrix::gemm(bool transa, bool transb, double alpha, const Matrix* a, const Matrix* b, double beta)
{
    char ta = transa ? 't' : 'n';
    char tb = transb ? 't' : 'n';
    int h, m, n, k, nca, ncb, ncc;

    for (h=0; h<nirreps_; ++h) {
        m = rowspi_[h];
        n = colspi_[h];
        k = transa ? a->rowspi_[h] : a->colspi_[h];
        nca = transa ? m : k;
        ncb = transb ? k : n;
        ncc = n;

        if (m && n && k) {
            C_DGEMM(ta, tb, m, n, k, alpha, &(a->matrix_[h][0][0]),
                    nca, &(b->matrix_[h][0][0]), ncb, beta, &(matrix_[h][0][0]),
                    ncc);
        }
    }
}

void Matrix::gemm(bool transa, bool transb, double alpha, shared_ptr<Matrix> a, shared_ptr<Matrix> b, double beta)
{
    gemm(transa, transb, alpha, a.get(), b.get(), beta);
}

void Matrix::gemm(bool transa, bool transb, double alpha, Matrix& a, shared_ptr<Matrix> b, double beta)
{
    gemm(transa, transb, alpha, const_cast<const Matrix*>(&a), b.get(), beta);
}

void Matrix::gemm(bool transa, bool transb, double alpha, shared_ptr<Matrix> a, Matrix& b, double beta)
{
    gemm(transa, transb, alpha, a.get(), const_cast<const Matrix*>(&b), beta);
}

void Matrix::gemm(bool transa, bool transb, double alpha, const Matrix& a, const Matrix& b, double beta)
{
    gemm(transa, transb, alpha, &a, &b, beta);
}

double Matrix::vector_dot(Matrix* rhs)
{
    double sum = 0.0;
    int h;
    size_t size;

    for (h=0; h<nirreps_; ++h) {
        size = rowspi_[h] * colspi_[h];
        if (size)
            sum += C_DDOT(size, (&matrix_[h][0][0]), 1, &(rhs->matrix_[h][0][0]), 1);
    }

    return sum;
}

double Matrix::vector_dot(shared_ptr<Matrix> rhs)
{
    return vector_dot(rhs.get());
}

void Matrix::diagonalize(Matrix* eigvectors, Vector* eigvalues)
{
    int h;
    for (h=0; h<nirreps_; ++h) {
        if (rowspi_[h]) {
            sq_rsp(rowspi_[h], colspi_[h], matrix_[h], eigvalues->vector_[h], 1, eigvectors->matrix_[h], 1.0e-14);
        }
    }
}

void Matrix::diagonalize(shared_ptr<Matrix> eigvectors, shared_ptr<Vector> eigvalues)
{
    diagonalize(eigvectors.get(), eigvalues.get());
}

void Matrix::diagonalize(shared_ptr<Matrix> eigvectors, Vector& eigvalues)
{
    diagonalize(eigvectors.get(), &eigvalues);
}

void Matrix::cholesky_factorize()
{
    for (int h=0; h<nirreps_; ++h) {
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

void Matrix::invert()
{
    for (int h=0; h<nirreps_; ++h) {
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

void Matrix::copy_lower_to_upper()
{
    for (int h=0; h<nirreps_; ++h) {
        for (int m=0; m<rowspi_[h]; ++m) {
            for (int n=0; n<m; ++n) {
                matrix_[h][n][m] = matrix_[h][m][n];
            }
        }
    }
}

void Matrix::copy_upper_to_lower()
{
    for (int h=0; h<nirreps_; ++h) {
        for (int m=0; m<rowspi_[h]; ++m) {
            for (int n=0; n<m; ++n) {
                matrix_[h][m][n] = matrix_[h][n][m];
            }
        }
    }
}

// Reference versions of the above functions:

void Matrix::transform(Matrix& a, Matrix& transformer)
{
    Matrix temp(a);

    temp.gemm(false, false, 1.0, a, transformer, 0.0);
    gemm(true, false, 1.0, transformer, temp, 0.0);
}

void Matrix::transform(Matrix& transformer)
{
    Matrix temp(this);

    temp.gemm(false, false, 1.0, *this, transformer, 0.0);
    gemm(true, false, 1.0, transformer, temp, 0.0);
}

void Matrix::back_transform(Matrix& a, Matrix& transformer)
{
    Matrix temp(a);

    temp.gemm(false, true, 1.0, a, transformer, 0.0);
    gemm(false, false, 1.0, transformer, temp, 0.0);
}

void Matrix::back_transform(Matrix& transformer)
{
    Matrix temp(this);

    temp.gemm(false, true, 1.0, *this, transformer, 0.0);
    gemm(false, false, 1.0, transformer, temp, 0.0);
}

double Matrix::vector_dot(Matrix& rhs)
{
    return vector_dot(&rhs);
}

void Matrix::diagonalize(Matrix& eigvectors, Vector& eigvalues)
{
    int h;
    for (h=0; h<nirreps_; ++h) {
        if (rowspi_[h]) {
            sq_rsp(rowspi_[h], colspi_[h], matrix_[h], eigvalues.vector_[h], 1, eigvectors.matrix_[h], 1.0e-14);
        }
    }
}

void Matrix::write_to_dpdfile2(dpdfile2 *outFile)
{

    dpd_file2_mat_init(outFile);

    if(outFile->params->nirreps != nirreps_) {
        char *str = new char[100];
        sprintf(str, "Irrep count mismatch.  Matrix class has %d irreps, but dpdfile2 has %d.",
                nirreps_, outFile->params->nirreps);
        throw SanityCheckError(str, __FILE__, __LINE__);
    }

    if(outFile->my_irrep != 0) {
        throw FeatureNotImplemented("libmints Matrix class",
                "Matrices whose irrep is not the symmetric one", __FILE__, __LINE__);
    }

    for(int h = 0; h < nirreps_; ++h){
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


void Matrix::save(const char *filename, bool append, bool saveLowerTriangle, bool saveSubBlocks)
{
    static const char *str_block_format = "%3d %3d %3d %20.15f\n";
    static const char *str_full_format  = "%3d %3d %20.15f\n";

    FILE *out = NULL;
    if (append == true) {
        out = fopen(filename, "a");
    } else {
        out = fopen(filename, "w");
    }

    fprintf(out, "%s", name_.c_str());
    fprintf(out, "\n");

    if (saveSubBlocks == false) {
        // Convert the matrix to a full matrix
        double **fullblock = to_block_matrix();

        // Need to know the size
        int sizer=0, sizec=0;
        for (int h=0; h<nirreps_; ++h) {
            sizer += rowspi_[h];
            sizec += colspi_[h];
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
            for (int h=0; h<nirreps_; ++h) {
                for (int i=0; i<rowspi_[h]; ++i) {
                    for (int j=0; j<=i; ++j) {
                        if (fabs(matrix_[h][i][j]) > 1.0e-14) {
                            count++;
                        }
                    }
                }
            }
            fprintf(out, "%5d\n", count);
            for (int h=0; h<nirreps_; ++h) {
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
            for (int h=0; h<nirreps_; ++h) {
                for (int i=0; i<rowspi_[h]; ++i) {
                    for (int j=0; j<colspi_[h]; ++j) {
                        if (fabs(matrix_[h][i][j]) > 1.0e-14) {
                            count++;
                        }
                    }
                }
            }
            fprintf(out, "%5d\n", count);
            for (int h=0; h<nirreps_; ++h) {
                for (int i=0; i<rowspi_[h]; ++i) {
                    for (int j=0; j<colspi_[h]; ++j) {
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

bool Matrix::load(psi::PSIO* psio, unsigned int fileno, char *tocentry, int nso)
{
    double *integrals = init_array(ioff[nso]);

    // If psi fails to read in the data this will abort out.
    if (tocentry != NULL)
        psi::IWL::read_one(psio, fileno, tocentry, integrals, ioff[nso], 0, 0, outfile);
    else
        psi::IWL::read_one(psio, fileno, const_cast<char*>(name_.c_str()), integrals, ioff[nso], 0, 0, outfile);

    set(integrals);

    ::free(integrals);

    return true;
}

bool Matrix::load(shared_ptr<psi::PSIO> psio, unsigned int fileno, char *tocentry, int nso)
{
    double *integrals = init_array(ioff[nso]);

    // If psi fails to read in the data this will abort out.
    if (tocentry != NULL)
        psi::IWL::read_one(psio.get(), fileno, tocentry, integrals, ioff[nso], 0, 0, outfile);
    else
        psi::IWL::read_one(psio.get(), fileno, const_cast<char*>(name_.c_str()), integrals, ioff[nso], 0, 0, outfile);

    set(integrals);

    ::free(integrals);

    return true;
}

void Matrix::save(psi::PSIO* psio, unsigned int fileno, bool saveSubBlocks)
{
    if(Communicator::world->me() == 0) {
        // Check to see if the file is open
        bool already_open = false;
        if (psio->open_check(fileno)) {
            already_open = true;
        } else {
            psio->open(fileno, PSIO_OPEN_OLD);
        }

        if (saveSubBlocks) {
            for (int h=0; h<nirreps_; ++h) {
                std::string str(name_);
                str += " Irrep " + to_string(h);

                // Write the sub-blocks
                if (colspi_[h] > 0 && rowspi_[h] > 0)
                    psio->write_entry(fileno, const_cast<char*>(str.c_str()), (char*)matrix_[h][0], sizeof(double) * colspi_[h] * rowspi_[h]);
            }
        } else {
            double **fullblock = to_block_matrix();
            // Need to know the size
            int sizer=0, sizec=0;
            for (int h=0; h<nirreps_; ++h) {
                sizer += rowspi_[h];
                sizec += colspi_[h];
            }

            // Write the full block
            if (sizer > 0 && sizec > 0)
                psio->write_entry(fileno, const_cast<char*>(name_.c_str()), (char*)fullblock[0], sizeof(double) * sizer * sizec);
            Matrix::free(fullblock);
        }

        if (!already_open)
            psio->close(fileno, 1);     // Close and keep
    }
}

void Matrix::save(shared_ptr<psi::PSIO> psio, unsigned int fileno, bool saveSubBlocks)
{
    save(psio.get(), fileno, saveSubBlocks);
}

void Matrix::send(Communicator* comm)
{
}

void Matrix::recv(Communicator* comm)
{
}

void Matrix::bcast(Communicator* comm, int broadcaster)
{
    // Assume the user allocated the matrix to the correct size first.
    for (int h=0; h<nirreps_; ++h) {
        comm->bcast(matrix_[h][0], rowspi_[h] * colspi_[h], broadcaster);
    }
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
    if (rhs->nirreps() != nirreps())
        return false;

    for (int h=0; h<nirreps(); ++h)
        if ((rowspi()[h] != rhs->rowspi()[h]) ||
                (colspi()[h] != rhs->colspi()[h]))
            return false;

    // Check element by element
    for (int h=0; h<nirreps(); ++h) {
        for (int m = 0; m < rowspi()[h]; ++m) {
            for (int n = 0; n < colspi()[h]; ++n) {
                if (get(h, m, n) != rhs->get(h, m, n))
                    return false;
            }
        }
    }

    return true;
}

//
// SimpleMatrix
//
SimpleMatrix::SimpleMatrix() : matrix_(0), rows_(0), cols_(0)
{

}

SimpleMatrix::SimpleMatrix(std::string name) :
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

SimpleMatrix::SimpleMatrix(shared_ptr<SimpleMatrix> c) : matrix_(0)
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

SimpleMatrix::SimpleMatrix(std::string name, int l_rows, int l_cols) : matrix_(0)
{
    rows_ = l_rows;
    cols_ = l_cols;
    name_ = name;
    alloc();
}

SimpleMatrix::SimpleMatrix(const Matrix* c) : matrix_(0), rows_(0), cols_(0)
{
    for (int i=0; i<c->nirreps(); ++i) {
        rows_ += c->rowspi()[i];
        cols_ += c->colspi()[i];
    }
    matrix_ = c->to_block_matrix();
}

SimpleMatrix::SimpleMatrix(const Matrix& c) : matrix_(0), rows_(0), cols_(0)
{
    for (int i=0; i<c.nirreps(); ++i) {
        rows_ += c.rowspi()[i];
        cols_ += c.colspi()[i];
    }
    matrix_ = c.to_block_matrix();
}

SimpleMatrix::SimpleMatrix(std::string& name, const Dimension& rows, const Dimension& cols)
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

void SimpleMatrix::init(int rowspi, int colspi, std::string name)
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

void SimpleMatrix::copy(shared_ptr<SimpleMatrix> cp)
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

void SimpleMatrix::set(shared_ptr<SimpleVector> vec)
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

    for(i=0;i<rows();i++) {
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

void SimpleMatrix::eivprint(shared_ptr<SimpleVector> values, FILE *out)
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

void SimpleMatrix::add(shared_ptr<SimpleMatrix> plus)
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

void SimpleMatrix::subtract(shared_ptr<SimpleMatrix> sub)
{
    subtract(sub.get());
}

void SimpleMatrix::accumulate_product(const SimpleMatrix* a, const SimpleMatrix* b)
{
    gemm(false, false, 1.0, a, b, 1.0);
}

void SimpleMatrix::accumulate_product(shared_ptr<SimpleMatrix> a, shared_ptr<SimpleMatrix> b)
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
    C_DSCAL(rows_, a, &(matrix_[m][0]), 1);
}

void SimpleMatrix::scale_column(int n, double a)
{
    C_DSCAL(cols_, a, &(matrix_[0][n]), rows_);
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

void SimpleMatrix::transform(shared_ptr<SimpleMatrix> a, shared_ptr<SimpleMatrix> transformer)
{
    transform(a.get(), transformer.get());
}

void SimpleMatrix::transform(SimpleMatrix* transformer)
{
    SimpleMatrix temp(this);

    temp.gemm(false, false, 1.0, this, transformer, 0.0);
    gemm(true, false, 1.0, transformer, &temp, 0.0);
}

void SimpleMatrix::transform(shared_ptr<SimpleMatrix> transformer)
{
    transform(transformer.get());
}

void SimpleMatrix::back_transform(SimpleMatrix* a, SimpleMatrix* transformer)
{
    SimpleMatrix temp(a);

    temp.gemm(false, true, 1.0, a, transformer, 0.0);
    gemm(false, false, 1.0, transformer, &temp, 0.0);
}

void SimpleMatrix::back_transform(shared_ptr<SimpleMatrix> a, shared_ptr<SimpleMatrix> transformer)
{
     back_transform(a.get(), transformer.get());
}

void SimpleMatrix::back_transform(SimpleMatrix* transformer)
{
    SimpleMatrix temp(this);

    temp.gemm(false, true, 1.0, this, transformer, 0.0);
    gemm(false, false, 1.0, transformer, &temp, 0.0);
}

void SimpleMatrix::back_transform(shared_ptr<SimpleMatrix> transformer)
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

void SimpleMatrix::gemm(bool transa, bool transb, double alpha, shared_ptr<SimpleMatrix> a, shared_ptr<SimpleMatrix> b, double beta)
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

double SimpleMatrix::vector_dot(shared_ptr<SimpleMatrix> rhs)
{
    return vector_dot(rhs.get());
}

void SimpleMatrix::diagonalize(SimpleMatrix* eigvectors, SimpleVector* eigvalues, int sort)
{
    if (rows_) {
        sq_rsp(rows_, cols_, matrix_, eigvalues->vector_, sort, eigvectors->matrix_, 1.0e-14);
    }
}

void SimpleMatrix::diagonalize(shared_ptr<SimpleMatrix> eigvectors, shared_ptr<SimpleVector> eigvalues, int sort)
{
    diagonalize(eigvectors.get(), eigvalues.get(), sort);
}

void SimpleMatrix::save(psi::PSIO* psio, unsigned int fileno)
{
    if(Communicator::world->me() == 0) {
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
}

void SimpleMatrix::load(shared_ptr<psi::PSIO> psio, unsigned int fileno)
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

void SimpleMatrix::save(shared_ptr<psi::PSIO> psio, unsigned int fileno)
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

//////////////////////////////////////////////////////////////////////////////
//
//  Class: View
//
//////////////////////////////////////////////////////////////////////////////

View::~View()
{
    nirrep_ = 0;
    delete[] row_offset_per_irrep_; row_offset_per_irrep_ = 0;
    delete[] col_offset_per_irrep_; col_offset_per_irrep_ = 0;
    delete[] rows_per_irrep_;       rows_per_irrep_ = 0;
    delete[] cols_per_irrep_;       cols_per_irrep_ = 0;
}

View::View(int nirrep, int *rows, int *cols)
    : nirrep_(nirrep), row_offset_per_irrep_(0), col_offset_per_irrep_(0),
      rows_per_irrep_(0), cols_per_irrep_(0)
{
    if (nirrep_ <= 0)
        throw InputException("Number of irreps is less than or equal to zero.", "nirrep", nirrep, __FILE__, __LINE__);
    if (rows == 0)
        throw InputException("Array of row sizes is 0.", "rows", 0, __FILE__, __LINE__);
    if (cols == 0)
        throw InputException("Array of column sizes is 0.", "cols", 0, __FILE__, __LINE__);

    rows_per_irrep_ = new int[nirrep_];
    cols_per_irrep_ = new int[nirrep_];
    row_offset_per_irrep_ = new int[nirrep_];
    col_offset_per_irrep_ = new int[nirrep_];

    for (int h=0; h<nirrep_; ++h) {
        rows_per_irrep_[h] = rows[h];
        cols_per_irrep_[h] = cols[h];
        row_offset_per_irrep_[h] = 0;
        col_offset_per_irrep_[h] = 0;
    }
}

View::View(int nirrep, int *rows, int *cols, int *row_offsets, int *col_offsets)
    : nirrep_(nirrep), row_offset_per_irrep_(0), col_offset_per_irrep_(0),
      rows_per_irrep_(0), cols_per_irrep_(0)
{
    if (nirrep_ <= 0)
        throw InputException("Number of irreps is less than or equal to zero.", "nirrep", nirrep, __FILE__, __LINE__);
    if (rows == 0)
        throw InputException("Array of row sizes is 0.", "rows", 0, __FILE__, __LINE__);
    if (cols == 0)
        throw InputException("Array of column sizes is 0.", "cols", 0, __FILE__, __LINE__);
    if (row_offsets == 0)
        throw InputException("Array of row offsets is 0.", "row_offsets", 0, __FILE__, __LINE__);
    if (col_offsets == 0)
        throw InputException("Array of column offsets is 0.", "col_offsets", 0, __FILE__, __LINE__);

    rows_per_irrep_ = new int[nirrep_];
    cols_per_irrep_ = new int[nirrep_];
    row_offset_per_irrep_ = new int[nirrep_];
    col_offset_per_irrep_ = new int[nirrep_];

    for (int h=0; h<nirrep_; ++h) {
        rows_per_irrep_[h] = rows[h];
        cols_per_irrep_[h] = cols[h];
        row_offset_per_irrep_[h] = row_offsets[h];
        col_offset_per_irrep_[h] = col_offsets[h];
    }
}

View::View(boost::shared_ptr<Matrix> matrix, int *rows, int *cols)
    : matrix_(matrix), nirrep_(0),
      row_offset_per_irrep_(0), col_offset_per_irrep_(0),
      rows_per_irrep_(0), cols_per_irrep_(0)

{
    nirrep_ = matrix_->nirreps();

    if (rows == 0)
        throw InputException("Array of row sizes is 0.", "rows", 0, __FILE__, __LINE__);
    if (cols == 0)
        throw InputException("Array of column sizes is 0.", "cols", 0, __FILE__, __LINE__);

    rows_per_irrep_ = new int[nirrep_];
    cols_per_irrep_ = new int[nirrep_];
    row_offset_per_irrep_ = new int[nirrep_];
    col_offset_per_irrep_ = new int[nirrep_];

    for (int h=0; h<nirrep_; ++h) {
        rows_per_irrep_[h] = rows[h];
        cols_per_irrep_[h] = cols[h];
        row_offset_per_irrep_[h] = 0;
        col_offset_per_irrep_[h] = 0;
    }
}

View::View(boost::shared_ptr<Matrix> matrix, int *rows, int *cols, int *row_offsets, int *col_offsets)
    : matrix_(matrix), nirrep_(0),
      row_offset_per_irrep_(0), col_offset_per_irrep_(0),
      rows_per_irrep_(0), cols_per_irrep_(0)

{
    nirrep_ = matrix_->nirreps();

    if (rows == 0)
        throw InputException("Array of row sizes is 0.", "rows", 0, __FILE__, __LINE__);
    if (cols == 0)
        throw InputException("Array of column sizes is 0.", "cols", 0, __FILE__, __LINE__);
    if (row_offsets == 0)
        throw InputException("Array of row offsets is 0.", "row_offsets", 0, __FILE__, __LINE__);
    if (col_offsets == 0)
        throw InputException("Array of column offsets is 0.", "col_offsets", 0, __FILE__, __LINE__);

    rows_per_irrep_ = new int[nirrep_];
    cols_per_irrep_ = new int[nirrep_];
    row_offset_per_irrep_ = new int[nirrep_];
    col_offset_per_irrep_ = new int[nirrep_];

    for (int h=0; h<nirrep_; ++h) {
        rows_per_irrep_[h] = rows[h];
        cols_per_irrep_[h] = cols[h];
        row_offset_per_irrep_[h] = row_offsets[h];
        col_offset_per_irrep_[h] = col_offsets[h];
    }
}

shared_ptr<Matrix> View::operator ()()
{
    // Create a new matrix with needed size
    shared_ptr<Matrix> matrix = shared_ptr<Matrix>(new Matrix(nirrep_, rows_per_irrep_, cols_per_irrep_));

    // Copy over the data
    for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<rows_per_irrep_[h]; ++i) {
            int ioffset = i + row_offset_per_irrep_[h];
            for (int j=0; j<cols_per_irrep_[h]; ++j) {
                int joffset = j + col_offset_per_irrep_[h];

                matrix->set(h, i, j, matrix_->get(h, ioffset, joffset));
            }
        }
    }

    return matrix;
}

shared_ptr<Matrix> View::view(boost::shared_ptr<Matrix> matrix)
{
    shared_ptr<Matrix> old = matrix_;
    matrix_ = matrix;
    return old;
}

void SimpleMatrix::swap_rows(int i, int j)
{
    C_DSWAP(cols_, &(matrix_[i][0]), 1, &(matrix_[j][0]), 1);
}

void SimpleMatrix::swap_columns(int i, int j)
{
    C_DSWAP(rows_, &(matrix_[0][i]), cols_, &(matrix_[0][j]), cols_);
}

