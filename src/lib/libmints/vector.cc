#include <cstdlib>
#include <string.h>
#include <libqt/qt.h>
#include "matrix.h"
#include "vector.h"
#include "dimension.h"

#include <libparallel/parallel.h>
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>

#include <algorithm>
#include <numeric>

using namespace boost;
using namespace psi;

Vector::Vector()
{
    dimpi_ = NULL;
    nirrep_ = 0;
    name_ = "";
}

Vector::Vector(const Vector& c)
{
    nirrep_ = c.nirrep_;
    dimpi_ = c.dimpi_;
    alloc();
    copy_from(c);
    name_ = c.name_;
}

Vector::Vector(int nirreps, int *dimpi)
    : dimpi_(nirreps)
{
    nirrep_ = nirreps;
    dimpi_ = dimpi;
    alloc();
}

Vector::Vector(int dim)
    : dimpi_(1)
{
    nirrep_ = 1;
    dimpi_[0] = dim;
    alloc();
}

Vector::Vector(const std::string& name, int nirreps, int *dimpi)
    : dimpi_(nirreps)
{
    nirrep_ = nirreps;
    dimpi_ = new int[nirrep_];
    for (int h=0; h<nirrep_; ++h)
        dimpi_[h] = dimpi[h];
    alloc();
    name_ = name;
}

Vector::Vector(const std::string& name, int dim)
    : dimpi_(1)
{
    nirrep_ = 1;
    dimpi_[0] = dim;
    alloc();
    name_ = name;
}

Vector::Vector(const Dimension& v)
{
    nirrep_ = v.n();
    dimpi_ = v;
    alloc();
    name_ = v.name();
}

Vector::Vector(const std::string& name, const Dimension& v)
{
    nirrep_ = v.n();
    dimpi_ = v;
    alloc();
    name_ = name;
}

Vector::~Vector()
{
    release();
}

void Vector::init(int nirreps, int *dimpi)
{
    dimpi_.init("", nirreps);
    nirrep_ = nirreps;
    dimpi_ = dimpi;
    alloc();
}

void Vector::init(int nirreps, const int *dimpi, const std::string& name)
{
    name_ = name;
    dimpi_.init("", nirreps);
    dimpi_ = dimpi;
    alloc();
}

void Vector::init(const Dimension &v)
{
    name_ = v.name();
    nirrep_ = v.n();
    dimpi_ = v;
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
    if (vector_.size())
        release();

    v_.resize(dimpi_.sum());

    std::fill(vector_.begin(), vector_.end(), (double*)0);
    std::fill(v_.begin(), v_.end(), 0.0);

    assign_pointer_offsets();
}

void Vector::assign_pointer_offsets()
{
    // Resize just to be sure it's the correct size
    vector_.resize(dimpi_.n(), 0);

    size_t offset = 0;
    for (int h=0; h<nirrep_; ++h) {
        if (dimpi_[h])
            vector_[h] = v_.data() + offset;
        else
            vector_[h] = NULL;
        offset += dimpi_[h];
    }
}

void Vector::release()
{
    std::fill(vector_.begin(), vector_.end(), (double*)0);
    std::fill(v_.begin(), v_.end(), 0.0);
}

void Vector::copy_from(const Vector& other)
{
    nirrep_ = other.dimpi_.n();
    dimpi_ = other.dimpi_;
    v_ = other.v_;
    assign_pointer_offsets();
}

void Vector::copy(const Vector *rhs)
{
    copy_from(*rhs);
}

void Vector::copy(const Vector &rhs)
{
    copy_from(rhs);
}

void Vector::set(double *vec)
{
    std::copy(vec, vec + dimpi_.sum(), v_.begin());
}

void Vector::zero()
{
    std::fill(v_.begin(), v_.end(), 0.0);
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
        fprintf(out, "\n # %s #\n", name_.c_str());
    } else {
        fprintf(out, "\n # %s %s #\n", name_.c_str(), extra);
    }
    for (h=0; h<nirrep_; ++h) {
        fprintf(out, " Irrep: %d\n", h+1);
        for (int i=0; i<dimpi_[h]; ++i)
            fprintf(out, "   %4d: %10.7f\n", i+1, vector_[h][i]);
        fprintf(out, "\n");
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

template<class T>
struct scale_vector {
    typedef T data_type;
    typedef data_type result_type;

    const data_type scalar;
    scale_vector(const data_type& sc) : scalar(sc) {}

    result_type operator()(data_type value) const { value *= scalar; return value; }
};

void Vector::scale(const double& sc)
{
    std::transform(v_.begin(), v_.end(), v_.begin(), scale_vector<double>(sc));
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
