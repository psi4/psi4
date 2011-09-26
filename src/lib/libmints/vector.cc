#include <cstdlib>
#include <string.h>
#include <libqt/qt.h>
#include "matrix.h"
#include "vector.h"
#include "dimension.h"

#include <libparallel/parallel.h>
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>

using namespace boost;
using namespace psi;

Vector::Vector()
{
    vector_ = NULL;
    dimpi_ = NULL;
    nirrep_ = 0;
    name_ = "";
}

Vector::Vector(const Vector& c)
{
    vector_ = NULL;
    nirrep_ = c.nirrep_;
    dimpi_ = new int[nirrep_];
    for (int h=0; h<nirrep_; ++h)
        dimpi_[h] = c.dimpi_[h];
    alloc();
    copy_from(c.vector_);
    name_ = c.name_;
}

Vector::Vector(int nirreps, int *dimpi)
{
    vector_ = NULL;
    nirrep_ = nirreps;
    dimpi_ = new int[nirrep_];
    for (int h=0; h<nirrep_; ++h)
        dimpi_[h] = dimpi[h];
    alloc();
}

Vector::Vector(int dim)
{
    vector_ = NULL;
    nirrep_ = 1;
    dimpi_ = new int[nirrep_];
    dimpi_[0] = dim;
    alloc();
}

Vector::Vector(const std::string& name, int nirreps, int *dimpi)
{
    vector_ = NULL;
    nirrep_ = nirreps;
    dimpi_ = new int[nirrep_];
    for (int h=0; h<nirrep_; ++h)
        dimpi_[h] = dimpi[h];
    alloc();
    name_ = name;
}

Vector::Vector(const std::string& name, int dim)
{
    vector_ = NULL;
    nirrep_ = 1;
    dimpi_ = new int[nirrep_];
    dimpi_[0] = dim;
    alloc();
    name_ = name;
}

Vector::Vector(const Dimension& v)
{
    nirrep_ = v.n();
    vector_ = NULL;
    dimpi_ = new int[nirrep_];
    for (int i=0; i<nirrep_; ++i)
        dimpi_[i] = v[i];
    alloc();
    name_ = v.name();
}

Vector::~Vector() {
    release();
    if (dimpi_)
        delete[] dimpi_;
}

void Vector::init(int nirreps, int *dimpi)
{
    if (dimpi_) delete[] dimpi_;
    nirrep_ = nirreps;
    dimpi_ = new int[nirrep_];
    for (int h=0; h<nirrep_; ++h)
        dimpi_[h] = dimpi[h];
    alloc();
}

void Vector::init(int nirreps, const int *dimpi, const std::string& name) {
    name_ = name;
    if (dimpi_) delete[] dimpi_;
    nirrep_ = nirreps;
    dimpi_ = new int[nirrep_];
    for (int h=0; h<nirrep_; ++h)
        dimpi_[h] = dimpi[h];
    alloc();
}

void Vector::init(const Dimension &v)
{
    name_ = v.name();
    if (dimpi_) delete[] dimpi_;
    nirrep_ = v.n();
    dimpi_ = new int[nirrep_];
    for (int h=0; h<nirrep_; ++h)
        dimpi_[h] = v[h];
    alloc();
}

Vector* Vector::clone()
{
    Vector *temp = new Vector(nirrep_, dimpi_);
    temp->copy(this);
    return temp;
}

void Vector::alloc()
{
    if (vector_)
        release();

    vector_ = (double**)malloc(sizeof(double*) * nirrep_);
    for (int h=0; h<nirrep_; ++h) {
        if (dimpi_[h]) {
            vector_[h] = new double[dimpi_[h]];
            ::memset(vector_[h], 0, sizeof(double)*dimpi_[h]);
        }
    }
}

void Vector::release()
{
    if (!vector_)
        return;

    for (int h=0; h<nirrep_; ++h) {
        if (dimpi_[h])
            delete[] (vector_[h]);
    }
    free(vector_);
    vector_ = NULL;
}

void Vector::copy_from(double **c)
{
    size_t size;
    for (int h=0; h<nirrep_; ++h) {
        size = dimpi_[h] * sizeof(double);
        if (size)
            memcpy(&(vector_[h][0]), &(c[h][0]), size);
    }
}

void Vector::copy(const Vector *rhs)
{
    if (nirrep_ != rhs->nirrep_) {
        release();
        if (dimpi_)
            delete[] dimpi_;
        nirrep_ = rhs->nirrep_;
        dimpi_ = new int[nirrep_];
        for (int h=0; h<nirrep_; ++h)
            dimpi_[h] = rhs->dimpi_[h];
        alloc();
    }
    copy_from(rhs->vector_);
}

void Vector::copy(const Vector &rhs)
{
    if (nirrep_ != rhs.nirrep_) {
        release();
        if (dimpi_)
            delete[] dimpi_;
        nirrep_ = rhs.nirrep_;
        dimpi_ = new int[nirrep_];
        for (int h=0; h<nirrep_; ++h)
            dimpi_[h] = rhs.dimpi_[h];
        alloc();
    }
    copy_from(rhs.vector_);
}

void Vector::set(double *vec)
{
    int h, i, ij;

    ij = 0;
    for (h=0; h<nirrep_; ++h) {
        for (i=0; i<dimpi_[h]; ++i) {
            vector_[h][i] = vec[ij++];
        }
    }
}

double Vector::pyget(const boost::python::tuple &key)
{
    int h = 0, elem = 0;
    h = python::extract<int>(key[0]);
    elem = python::extract<int>(key[1]);

    return get(h, elem);
}

void Vector::pyset(const boost::python::tuple &key, double value)
{
    int h = 0, elem = 0;
    h = python::extract<int>(key[0]);
    elem = python::extract<int>(key[1]);

    set(h, elem, value);
}

double Vector::pyget(int key)
{
    int h = 0, elem = key;
    return get(h, elem);
}

void Vector::pyset(int key, double value)
{
    int h = 0, elem = key;
    set(h, elem, value);
}

void Vector::print(FILE* out, const char* extra) const
{
    int h;
    if (extra == NULL) {
        fprintf(outfile, "\n # %s #\n", name_.c_str());
    } else {
        fprintf(outfile, "\n # %s %s #\n", name_.c_str(), extra);
    }
    for (h=0; h<nirrep_; ++h) {
        fprintf(outfile, " Irrep: %d\n", h+1);
        for (int i=0; i<dimpi_[h]; ++i)
            fprintf(outfile, "   %4d: %10.7f\n", i+1, vector_[h][i]);
        fprintf(outfile, "\n");
    }
}

double *Vector::to_block_vector() {
    size_t size=0;
    for (int h=0; h<nirrep_; ++h)
        size += dimpi_[h];

    double *temp = new double[size];
    size_t offset = 0;
    for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<dimpi_[h]; ++i) {
            temp[i+offset] = vector_[h][i];
        }
        offset += dimpi_[h];
    }

    return temp;
}

void Vector::gemv(bool transa, double alpha, Matrix* A, Vector* X, double beta)
{
    char trans = transa ? 't' : 'n';

    for (int h =0; h < nirrep_; ++h) {
        C_DGEMV(trans, A->rowspi_[h], A->colspi_[h], alpha, &(A->matrix_[h][0][0]),
                A->rowspi_[h], &(X->vector_[h][0]), 1, beta, &(vector_[h][0]), 1);
    }
}

double Vector::dot(Vector* X)
{
    double tmp = 0.0;

    for (int h=0; h<nirrep_; ++h) {
        if (dimpi_[h] != X->dimpi_[h]) {
            printf("Vector::dot: Vectors are not of the same size.\n");
            return 0.0;
        }
        for (int i=0; i<dimpi_[h]; ++i) {
            tmp += vector_[h][i] * X->vector_[h][i];
        }
    }

    return tmp;
}

void Vector::scale(double sc)
{
    for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<dimpi_[h]; ++i) {
            vector_[h][i] *= sc;
        }
    }
}

void Vector::send()
{
}

void Vector::recv()
{
}

void Vector::bcast(int broadcaster)
{
    // Assume the user allocated the matrix to the correct size first.
    for (int h=0; h<nirrep_; ++h) {
        if (dimpi_[h] > 0)
            Communicator::world->bcast(vector_[h], dimpi_[h], broadcaster);
    }
}

void Vector::sum()
{
    for (int h=0; h<nirrep_; ++h)
        if (dimpi_[h] > 0)
            Communicator::world->sum(vector_[h], dimpi_[h]);
}

//
// SimpleVector class
//
SimpleVector::SimpleVector() {
    vector_ = NULL;
    dim_ = 0;
}

SimpleVector::SimpleVector(const SimpleVector& c) {
    vector_ = NULL;
    dim_ = c.dim_;
    alloc();
    copy_from(c.vector_);
}

SimpleVector::SimpleVector(int dim) {
    vector_ = NULL;
    dim_ = dim;
    alloc();
}

SimpleVector::~SimpleVector() {
    release();
}

void SimpleVector::init(int dim)
{
    dim_ = dim;
    alloc();
}

void SimpleVector::alloc() {
    if (vector_)
        release();

    vector_ = new double[dim_];
    memset(vector_, 0, sizeof(double) * dim_);
}

void SimpleVector::release() {
    if (!vector_)
        return;

    delete[] (vector_);
    vector_ = NULL;
}

void SimpleVector::copy_from(double *c) {
    size_t size;
    size = dim_ * sizeof(double);
    if (size)
        memcpy(&(vector_[0]), &(c[0]), size);
}

void SimpleVector::copy(const SimpleVector *rhs) {
    release();
    dim_ = rhs->dim_;
    alloc();
    copy_from(rhs->vector_);
}

void SimpleVector::set(double *vec) {
    int i;

    for (i=0; i<dim_; ++i) {
        vector_[i] = vec[i];
    }
}

void SimpleVector::print() const {
    for (int i=0; i<dim_; ++i)
        fprintf(outfile, "   %4d: %10.7f\n", i+1, vector_[i]);
    fprintf(outfile, "\n");
}

double *SimpleVector::to_block_vector() {
    double *temp = new double[dim_];
    for (int i=0; i<dim_; ++i) {
        temp[i] = vector_[i];
    }

    return temp;
}

void SimpleVector::scale(double a) {
    for (int i=0; i<dim_; ++i) {
        vector_[i] *= a;
    }
}

double SimpleVector::magnitude() const
{
    double mag=0.0;
    for (int i=0; i<dim_; ++i)
        mag += (vector_[i] * vector_[i]);
    return sqrt(mag);
}
