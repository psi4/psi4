#include <libyeti/gigmatrix.h>
#include <cstdarg>
#include <algorithm>
#include <cstring>

#include <libsmartptr/timer.h>
#include <libsmartptr/printstream.h>

#define VALGRIND 0

#include "blas.h"
extern "C" {

extern void F_DGEMM(const char*, const char*, const int*,
  const int*, const int*, const double*, const double*, const int*,
  const double*, const int*, const double*, double*, const int*);

extern void F_DGEEV(const char*, const char*, const int*, const double*,
  const int*, const double*, const double*, const double*, const int*,
  const double*, const int*, const double*, const int*, const int*);

extern void F_DSYEV(const char*, const char*, const int*, const double*,
  const int*, const double*, const double*, const int*, const int*);


#if NEED_MPQC_DGEEV_WRAPPER
void dgeev(const char* a, const char* b, const int* c, const double* d,
  const int* e, const double* f, const double* g, const double* h, const int* i,
  const double* j, const int* k, const double* l, const int* m, const int* n)
{
    F_DGEEV(
        a,b,c,d,e,f,g,h,i,j,k,l,m,n
    );
}
#endif

}//EndExternC

using namespace yeti;
using namespace std;

SerialDeclare(Matrix);
SerialDeclare(Vector);

#define X 0
#define Y 1
#define Z 2

void
yeti::error(
    const Matrix* l,
    const Matrix* r,
    const std::string& op
)
{
    cerr << stream_printf("Matrices of size %dx%d and %dx%d not aligned for %s",
                  l->nrow(), l->ncol(),
                  r->nrow(), r->ncol(),
                  op.c_str()) << endl;
    abort();
}

void
yeti::error(
    const Vector* l,
    const Vector* r,
    const std::string& op
)
{
    cerr << stream_printf("Vector of size %d and %d not aligned for %s",
                  l->n(),
                  r->n(),
                  op.c_str()) << endl;
    abort();
}

void
yeti::scale_array(
    double* v,
    double s,
    uli n
)
{
    double* vptr = v;
    for (uli i=0; i < n; ++i, ++vptr)
        (*vptr) *= s;
}

void
yeti::assign_array(
    double* v,
    double d,
    uli n
)
{
    double* vptr = v;
    for (uli i=0; i < n; ++i, ++vptr)
        (*vptr) = d;
}

double 
yeti::max_array(
    const double* v,
    uli n
)
{
    double max = 0;
    for (uli i=0; i < n; ++i)
    {
        double newval = fabs(v[i]);
        if (newval > max)
            max = newval;
    }
    return max;
}

void
yeti::transpose(
    double* vals,
    uli nrow,
    uli ncol
)
{
    double* scratch = new double[nrow*ncol];
    double* scratchptr = scratch;
    double* ptr;
    for (uli i=0; i < ncol; ++i)
    {
        ptr = vals + i;
        for (uli j=0; j < nrow; ++j, ptr += ncol, ++scratchptr)
        {
            *scratchptr = *ptr;
        }
    }
    memcpy(vals, scratch, nrow*ncol*sizeof(double));

    delete[] scratch;
}

void
yeti::eigenvalues(
    const Matrix* matrix,
    Vector* evals,
    Matrix* evecs
)
{

    uli n = matrix->nrow();
    uli nblock = n * n;

    //evecs end up here
    double* m = new double[nblock];
    memcpy(m, matrix->data(), nblock * sizeof(double));
    double* vals = new double[n];
    int worksize = 8*n;
    double* work = new double[worksize];


#if VALGRIND
    double* matptr = m;
    for (uli i=0; i < n; ++i)
    {
        vals[i] = 1;
        for (uli j=0; j < n; ++j, ++matptr)
        {
            if (i==j)
            {
                *matptr = 1;
            }
            else
            {
                *matptr = 0;
            }
        }
    }
#else
    //loss of precision... but what can you do
    int n_ = (int) n;
    int info = 0;
    F_DSYEV("V", "U", &n_, m, &n_, vals, work, &worksize, &info);
    if (info != 0)
    {
        cerr << "Eigenvalue routined failed" << endl;
        matrix->print("matrix", cerr);
        abort();
    }
#endif


    transpose(m, n, n);
    evecs->assign(m);
    evals->assign(vals);

    delete[] work;
    delete[] m;
    delete[] vals;
}

void
yeti::eigenvalues(
    const Matrix* matrix,
    Vector* rvals,
    Vector* ivals,
    Matrix* Levecs,
    Matrix* Revecs
)
{
    if (matrix->nrow() != matrix->ncol())
    {
        cerr << "Attempting to compute eigenvalues of non-square matrix" << endl;
        abort();
    }

    uli n = matrix->nrow();
    uli nblock = n * n;
    //must memcpy since dgeev overwrites array
    double* vals = new double[nblock]; memcpy(vals, matrix->data(), nblock * sizeof(double));
    double* levecs = new double[nblock];
    double* revecs = new double[nblock];

    double* revals = new double[n];
    double* ievals = new double[n];
    int worksize = 8*n;
    double* work = new double[worksize];


#if VALGRIND
    double* revecptr = revecs;
    double* levecptr = levecs;
    double* matptr = vals;
    for (uli i=0; i < n; ++i)
    {
        revals[i] = 1;
        ievals[i] = 0;
        for (uli j=0; j < n; ++j, ++revecptr, ++levecptr, ++matptr)
        {
            if (i==j)
            {
                *revecptr = 1;
                *levecptr = 1;
            }
            else
            {
                *revecptr = 0;
                *levecptr = 0;
            }
        }
    }
#else
    //loss of precision... but what can you do
    int n_ = (int) n;
    int info = 0;
    F_DGEEV("V", "V", &n_, vals, &n_, revals, ievals, levecs, &n_, revecs, &n_, work, &worksize, &info);

    if (info != 0)
    {
        cerr << "Eigenvalue routined failed" << endl;
        matrix->print("matrix", cerr);
        abort();
    }
#endif


    transpose(revecs, n, n);
    transpose(levecs, n, n);

    //these are transposed
    Levecs->assign(revecs);
    Revecs->assign(levecs);
    rvals->assign(revals);
    ivals->assign(ievals);

    delete[] revals;
    delete[] ievals;
    delete[] levecs;
    delete[] revecs;
    delete[] vals;
    delete[] work;
}

double*
yeti::subtract_arrays(
    const double* l,
    const double* r,
    uli n
)
{
    const double* lptr = l;
    const double* rptr = r;
    double* vals = new double[n];
    double* valptr = vals;
    
    for (uli i=0; i < n; ++i, ++lptr, ++rptr, ++valptr)
        (*valptr) = (*lptr) - (*rptr);

    return vals;
}

double*
yeti::add_arrays(
    const double* l,
    const double* r,
    uli n
)
{
    const double* lptr = l;
    const double* rptr = r;
    double* vals = new double[n];
    double* valptr = vals;
    
    for (uli i=0; i < n; ++i, ++lptr, ++rptr, ++valptr)
        (*valptr) = (*lptr) + (*rptr);

    return vals;
}

void
yeti::accumulate_array(
    double* target,
    const double* src,
    uli n
)
{
    double* tptr = target;
    const double* sptr = src;
    for (uli i=0; i < n; ++i, ++tptr, ++sptr)
        (*tptr) += (*sptr);
}

Matrix*
yeti::multiply(
    const Matrix* l,
    const Matrix* r,
    bool transpose_l,
    bool transpose_r
)
{
    const double* ldata = l->data();
    const double* rdata = r->data();

    const char* opl;
    const char* opr;
    uli nrow, ncol, nlink_l, nlink_r;

    //compute l * r
    if (transpose_l)
    {
        opl = "N";
        nrow = l->ncol();
        nlink_l = l->nrow();
    }
    else
    {
        opl = "T";
        nrow = l->nrow();
        nlink_l = l->ncol();
    }

    if (transpose_r)
    {
        opr = "N";
        ncol = r->nrow();
        nlink_r = r->ncol();
    }
    else
    {
        opr = "T";
        ncol = r->ncol();
        nlink_r = r->nrow();
    }

    if (nlink_l != nlink_r)
    {
        stringstream sstr;
        sstr << "multiplication(" << opl << "," << opr << ")";
        error(l, r, sstr.str());
    }

    uli nlink = nlink_r = nlink_l;
    uli ldl = l->ncol();
    uli ldr = r->ncol();

    uli nblock = nrow * ncol;

    double* prod_data = new double[nblock];
    memset(prod_data, 0, nblock * sizeof(double));
    double alpha = 1.0;
    double beta = 0.0;

#if VALGRIND
    memset(prod_data, 0, ncol * nrow * sizeof(double));
    double* prodptr = prod_data;
    for (uli i=0; i < nrow; ++i)
    {
        for (uli j=0; j < ncol; ++j, ++prodptr)
        {
            double Xij = 0;
            const double* jptr = rdata + j * ncol;
            const double* iptr = ldata + i * nrow;
            for (uli k=0; k < nlink; ++k, ++iptr, ++jptr)
            {
                Xij += (*iptr) * (*jptr);
            }
            *prodptr = Xij;
        }
    }
#else
    int nlink_ = (int) nlink;
    int nrow_ = (int) nrow;
    int ncol_ = (int) ncol;
    int ldl_ = (int) ldl;
    int ldr_ = (int) ldr;
    F_DGEMM(opl, opr, &nrow_, &ncol_, &nlink_, &alpha, ldata, &ldl_, rdata, &ldr_, &beta, prod_data, &nrow_);
#endif

    //this is all transposed
    transpose(prod_data, ncol, nrow);

    Matrix* prod = new Matrix(nrow, ncol);
    prod->assign(prod_data);

    delete[] prod_data;

    return prod;
}

Matrix::Matrix(uli nrow, uli ncol)
    : nrow_(nrow), ncol_(ncol)
{
    SetRuntime(Matrix);

    data_ = new double[nrow * ncol];
    memset(data_, 0, nrow * ncol * sizeof(double));
}

Matrix::Matrix(double* d, uli nrow, uli ncol)
    : nrow_(nrow), ncol_(ncol), data_(d)
{
    SetRuntime(Matrix);
}

Matrix::Matrix(Vector* v, uli nrow, uli ncol)
    : nrow_(nrow), ncol_(ncol)
{
    SetRuntime(Matrix);

    data_ = new double[nrow * ncol];
    assign(v->data());
}

Matrix::Matrix(const XMLArchivePtr& arch)
{
    SetRuntime(Matrix);

    uli size;
    arch->getBinary<double>(data_, size, "data");
    serial_load(nrow);
    serial_load(ncol);

    if (nrow_ * ncol_ != size)
    {
        cerr << "Matrix data not aligned with nrow x ncol" << endl;
        abort();
    }
}

void
Matrix::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);
    uli size = nrow_ * ncol_;
    arch->setBinary<double>(data_, size, "data");
    serial_save(nrow);
    serial_save(ncol);
}

Matrix::~Matrix()
{
    delete[] data_;
}

Vector*
Matrix::toVector() const
{
    uli nblock = nrow_ * ncol_;
    Vector* v = new Vector(nblock);
    v->assign(data());
    return v;
}

Matrix*
Matrix::add(const Matrix* r) const
{
    if (nrow() != r->nrow() || ncol() != r->ncol())
        error(this, r, "addition");

    const double* ldata = data();
    const double* rdata = r->data();

    double* vals = add_arrays(ldata, rdata, nrow() * ncol());
    Matrix* newptr = new Matrix(nrow(), ncol());
    newptr->assign(vals);

    delete[] vals;

    return newptr;
}


Matrix*
Matrix::subtract(const Matrix* r) const
{
    if (nrow() != r->nrow() || ncol() != r->ncol())
        error(this, r, "subtraction");

    const double* ldata = data();
    const double* rdata = r->data();

    double* vals = subtract_arrays(ldata, rdata, nrow_ * ncol_);
    Matrix* newptr = new Matrix(nrow(), ncol());
    newptr->assign(vals);

    delete[] vals;

    return newptr;
}


Matrix*
Matrix::mult(const Matrix* r) const
{
    Matrix* newptr = multiply(this, r);
    return newptr;
}

Matrix*
Matrix::mult(double d) const
{
    const double* vals = data();
    const double* ptr = vals;
    uli n = nrow() * ncol();

    double* newvals = new double[n];
    double* newptr = newvals;
    for (uli i=0; i < n; ++i, ++ptr, ++newptr)
        (*newptr) = d * (*ptr);
    Matrix* m = new Matrix(nrow(), ncol());
    m->assign(newvals);
    
    delete[] newvals;

    return m;
}

Vector*
Matrix::mult(const Vector* v) const
{
    Matrix* m = v->toMatrix(v->n(), 1);
    Matrix* t = mult(m);
    Vector* newv = t->toVector();
    delete m;
    delete t;
    return newv;
}

Matrix*
Matrix::t() const 
{
    uli nblock = nrow_ * ncol_;
    double* vals = new double[nblock];
    memcpy(vals, data_, nblock * sizeof(double));
    transpose(vals, nrow_, ncol_);
    Matrix* m = new Matrix(vals, ncol_, nrow_);
    return m;
}

const double*
Matrix::data() const
{
    return data_;
}

void
Matrix::accumulate(const Matrix* r)
{
    accumulate_array(data_, r->data(), nrow_ * ncol_);
}

void
Matrix::canonicalize()
{
    double* newdata = new double[nrow_ * ncol_];
    uli* cols_taken = new uli[ncol_];
    for (uli i=0; i < ncol_; ++i)
        cols_taken[i] = 0;

    for (uli col=0; col < ncol_; ++col)
    {
        uli maxrow = 0;
        double maxabs = -1;
        double maxval = 0;
        //figure out the maximum element in the column
        for (uli row=0; row < nrow_; ++row)
        {
            double val = get_element(row,col);
            double abs = fabs(val);
            if (abs > maxabs && !cols_taken[row])
            {
                maxabs = abs;
                maxval = val;
                maxrow = row;
            }
        }

        if (maxabs == -1)
        {
            cerr << "Eigenvector switch. Cannot proceed!" << endl;
            abort();
        }

        double scale = maxval < 0 ? -1 : 1;

        uli idx = 0;
        for (uli row=0; row < nrow_; ++row)
        {
            newdata[row * ncol_ + maxrow] = data_[row * ncol_ + col] * scale;
        }
        cols_taken[maxrow] = 1;
    }
    delete[] data_;
    delete[] cols_taken;
    data_ = newdata;
}

uli
Matrix::nrow() const
{
    return nrow_;
}

uli
Matrix::ncol() const
{
    return ncol_;
}

void
Matrix::scale(double s)
{
    scale_array(data_, s, nrow_ * ncol_);
}

void
Matrix::assign(const double* newvals)
{
    memcpy(data_, newvals, ncol_ * nrow_ * sizeof(double));
}

void
Matrix::assign(double d)
{
    assign_array(data_, d, ncol_ * nrow_);
}

void
Matrix::set_element(uli i, uli j, double val)
{
    if (i >= nrow_ || j >= ncol_)
    {
        cerr << stream_printf("Invalid indices (%d,%d)", i, j) << endl;
        abort();
    }
    uli index = i * ncol_ + j;
    data_[index] = val;
}

double
Matrix::get_element(uli i, uli j) const
{

    uli index = i * ncol_ + j;
    return data_[index];
}

void
Matrix::accumulate_element(uli i, uli j, double val)
{
    uli index = i * ncol_ + j;
    data_[index] += val;
}

void
Matrix::printValues(std::ostream& os) const
{
    for (uli colstart=0; colstart < ncol_; colstart += 5)
    {
        uli colstop = colstart + 5;
        if (colstop > ncol_) colstop = ncol_;

        os << "  ";
        for (uli col=colstart; col < colstop; ++col)
        {
            os << stream_printf("  %16d", col);
        }
        os << endl;

        for (uli row=0; row < nrow_; ++row)
        {
            os << "  ";
            for (uli col=colstart; col < colstop; ++col)
            {
                double val = get_element(row,col);
                if (fabs(val) > 100)
                    os << stream_printf("  %16.8e", val);
                else
                    os << stream_printf("  %16.10f", val);
            }
            os << endl;
        }
    }
}

void
Matrix::print(const std::string& title, std::ostream& os) const
{
    os << stream_printf("%d x %d Matrix %s", nrow_, ncol_, title.data()) << endl;
    printValues(os);
}

void
Matrix::print(std::ostream& os) const
{
    print("no title", os);
}

Vector::Vector(uli n) :
    n_(n)
{
    SetRuntime(Vector);

    data_ = new double[n];
    memset(data_, 0, n * sizeof(double));
}

Vector::Vector(double* d, uli n) :
    n_(n), data_(d)
{
    SetRuntime(Vector);
}

Vector::Vector(const XMLArchivePtr& arch)
{
    SetRuntime(Vector);

    uli size;
    arch->getBinary(data_, size, "data");
    serial_load(n);

    if (size != n_)
    {
        cerr << stream_printf("Data size %ld does not match dim %ld", size, n_) << endl;
        abort();
    }
}

void
Vector::serialize(const XMLArchivePtr& arch) const
{
    Serializable::serialize(arch);

    arch->setBinary(data_, n_, "data");
    serial_save(n);
}

Matrix*
Vector::toMatrix(uli nrow, uli ncol) const
{
    if (nrow * ncol != n())
    {
        cerr << stream_printf("Vector of size %d cannot be converted to matrix of size %d x %d",
                              n(), nrow, ncol) << endl;
        abort();
    }
        
    Matrix* m = new Matrix(nrow, ncol);
    m->assign(data());
    return m;
}

Vector::~Vector()
{
    delete[] data_;
}

Vector*
Vector::add(const Vector* v) const
{
    if (n() != v->n())
        error(this, v, "addition");

    double* vals = add_arrays(data_, v->data(), n());

    Vector* newptr = new Vector(n());
    newptr->assign(vals);

    delete[] vals;

    return newptr;
}

Vector*
Vector::subtract(const Vector* v) const
{
    if (n() != v->n())
        error(this, v, "subtraction");

    double* newvals = subtract_arrays(data(), v->data(), n()); 
    Vector* vec = new Vector(n());
    vec->assign(newvals);

    delete[] newvals;

    return vec;
}

Vector*
Vector::mult(double d) const
{
    double* vals = new double[n()];
    for (uli i=0; i < n(); ++i)
        vals[i] = data_[i] * d;

    Vector* v = new Vector(n());
    v->assign(vals);

    delete[] vals;

    return v;
}

const double*
Vector::data() const
{
    return data_;
}

void
Vector::accumulate(const Vector* v)
{
    accumulate_array(data_, v->data(), n());
}

void
Vector::scale(double s)
{
    scale_array(data_, s, n_);
}

void
Vector::assign(const double* vals)
{
    memcpy(data_, vals, n() * sizeof(double));
}

void
Vector::assign(double d)
{
    assign_array(data_, d, n_);
}

double
Vector::dot(const Vector* v) const
{
    if (n() != v->n())
        error(this, v, "dot product");

    const double* ldata = data();
    const double* rdata = v->data();
    const double* lptr = ldata;
    const double* rptr = rdata;
    
    double sum = 0;
    for (uli i=0; i < n(); ++i, ++lptr, ++rptr)
        sum += (*lptr) * (*rptr);
        
    return sum;
}

double
Vector::get_element(uli i) const
{
    return data_[i];
}

void
Vector::set_element(uli i, double val)
{
    data_[i] = val;
}

void
Vector::accumulate_element(uli i, double val)
{
    data_[i] += val;
}

void
Vector::sort()
{
    std::sort(data_, data_ + n_);
}

void
Vector::printValues(std::ostream& os) const
{
    for (uli i=0; i < n_; ++i)
    {
        os << stream_printf("  %16.10f", get_element(i));
        os << endl;
    }
}

void
Vector::print(std::ostream& os) const
{
    print("no title", os);
}

void
Vector::print(const std::string& title, std::ostream& os) const
{
    os << stream_printf("Vector %s of size %d", title.data(), n_) << endl;
    printValues(os);
}

uli
Vector::n() const
{
    return n_;
}

Vector1::Vector1(double a0)
    : VectorPtr(1)
{
    set_element(0, a0);
}

Vector2::Vector2(double a0, double a1)
    : VectorPtr(2)
{
    set_element(0, a0);
    set_element(1, a1);
}

Vector3::Vector3(double a0, double a1, double a2)
    : VectorPtr(3)
{
    set_element(0, a0);
    set_element(1, a1);
    set_element(2, a2);
}

Vector4::Vector4(double a0, double a1, double a2, double a3)
    : VectorPtr(4)
{
    set_element(0, a0);
    set_element(1, a1);
    set_element(2, a2);
    set_element(3, a3);
}

Vector8::Vector8(double a0, double a1, double a2, double a3,
                 double a4, double a5, double a6, double a7)
    : VectorPtr(8)
{
    set_element(0, a0);
    set_element(1, a1);
    set_element(2, a2);
    set_element(3, a3);
    set_element(4, a4);
    set_element(5, a5);
    set_element(6, a6);
    set_element(7, a7);
}


VectorPtr
yeti::cross(const VectorPtr& v, const VectorPtr& w)
{
    v.nullcheck();
    w.nullcheck();

    VectorPtr cp = v.clone();
    double newx = v[Y]*w[Z] - v[Z]*w[Y];
    double newy = v[Z]*w[X] - v[X]*w[Z];
    double newz = v[X]*w[Y] - v[Y]*w[X];
    cp.set_element(X, newx);
    cp.set_element(Y, newy);
    cp.set_element(Z, newz);
    return cp;
}

bool
yeti::equals(const RectMatrixPtr& l, const RectMatrixPtr& r, double tol)
{
    l.nullcheck();
    r.nullcheck();

    RectMatrixPtr diff = l - r;
    return diff.maxabs() < tol;
}

bool
yeti::equals(const VectorPtr& l, const VectorPtr& r, double tol)
{
    l.nullcheck();
    r.nullcheck();

    VectorPtr diff = l - r;
    return diff.maxabs() < tol;
}

RectMatrixPtr yeti::operator*(const RectMatrixPtr& l, const RectMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();

    Matrix* m = l.get()->mult(r.get());    

    return m;
}

RectMatrixPtr yeti::operator*(const RectMatrixPtr& l, const SymmMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();
    Matrix* m = l.get()->mult(r.get());    

    return m;
}

RectMatrixPtr yeti::operator*(const SymmMatrixPtr& l, const RectMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();
    Matrix* m = l.get()->mult(r.get());    

    return m;
}

RectMatrixPtr yeti::operator*(const SymmMatrixPtr& l, const SymmMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();
    Matrix* m = l.get()->mult(r.get());    

    return m;
}

RectMatrixPtr yeti::operator*(const RectMatrixPtr& m, double d)
{
    m.nullcheck();

    Matrix* p = m.get()->mult(d);

    return p;
}

RectMatrixPtr yeti::operator*(double d, const RectMatrixPtr& m)
{
    m.nullcheck();

    Matrix* p = m.get()->mult(d);

    return p;
}

RectMatrixPtr yeti::operator+(const RectMatrixPtr& l, const RectMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();

    Matrix* m = l.get()->subtract(r.get());

    return m;
}

RectMatrixPtr yeti::operator-(const RectMatrixPtr& l, const RectMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();

    Matrix* m = l.get()->subtract(r.get());

    return m;
}

SymmMatrixPtr yeti::operator+(const SymmMatrixPtr& l, const SymmMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();

    Matrix* m = l.get()->add(r.get());

    return m;
}

SymmMatrixPtr yeti::operator-(const SymmMatrixPtr& l, const SymmMatrixPtr& r)
{
    l.nullcheck();
    r.nullcheck();

    Matrix* m = l.get()->subtract(r.get());

    return m;
}

SymmMatrixPtr yeti::operator*(const SymmMatrixPtr& m, double d)
{
    m.nullcheck();

    Matrix* p = m.get()->mult(d);

    return p;
}

SymmMatrixPtr yeti::operator*(double d, const SymmMatrixPtr& m)
{
    m.nullcheck();

    Matrix* p = m.get()->mult(d);

    return p;
}

VectorPtr yeti::operator*(const RectMatrixPtr& m, const VectorPtr& v)
{
    m.nullcheck();
    v.nullcheck();
    Vector* p =  m.get()->mult(v.get());

    return p;
}

VectorPtr yeti::operator*(const SymmMatrixPtr& m, const VectorPtr& v)
{
    m.nullcheck();
    v.nullcheck();
    Vector* p =  m.get()->mult(v.get());

    return p;
}

VectorPtr yeti::operator*(double d, const VectorPtr& v)
{
    v.nullcheck();
    Vector* p = v.get()->mult(d);

    return p;
}

VectorPtr yeti::operator*(const VectorPtr& v, double d)
{
    v.nullcheck();

    Vector* p = v.get()->mult(d);

    return p;
}

VectorPtr yeti::operator+(const VectorPtr& l, const VectorPtr& r)
{
    l.nullcheck();
    r.nullcheck();

    Vector* v = l.get()->add(r.get());

    return v;
}

VectorPtr yeti::operator-(const VectorPtr& l, const VectorPtr& r)
{
    l.nullcheck();
    r.nullcheck();

    Vector* v = l.get()->subtract(r.get());

    return v;
}

SymmMatrixPtr::SymmMatrixPtr()
    : Parent()
{
}

SymmMatrixPtr::SymmMatrixPtr(Matrix* m)
    : Parent(m)
{
}

SymmMatrixPtr::SymmMatrixPtr(uli n)
    : Parent(n)
{
}


SymmMatrixPtr::SymmMatrixPtr(const SymmMatrixPtr& m)
    : Parent(m)
{
}

RectMatrixPtr::RectMatrixPtr()
    : Parent()
{
}

RectMatrixPtr::RectMatrixPtr(Matrix* m)
    : Parent(m)
{
}

RectMatrixPtr::RectMatrixPtr(uli nrow, uli ncol)
    : Parent(nrow, ncol)
{
}

RectMatrixPtr::RectMatrixPtr(const RectMatrixPtr& m)
    : Parent(m)
{
}

VectorPtr::VectorPtr()
    : Parent()
{
}

VectorPtr::VectorPtr(uli n)
    : Parent(n)
{
}

VectorPtr::VectorPtr(Vector* v)
    : Parent(v)
{
}

VectorPtr::VectorPtr(const VectorPtr& v)
    : Parent(v.get())
{
}

#if 0
void
yeti::instantiate()
{
    RectMatrixPtr m;
    VectorPtr v(3);
    SymmMatrixPtr s;
    const double* vals;

    RectMatrixTemplate<Matrix>* mptr = new RectMatrixTemplate<Matrix>;
    mptr->print("test");
    delete mptr;

    m.t();
    m.get_subblock(0,0,0,0);
    m.copy();
    m.clone();
    m.set_element(0,0,0);
    m.accumulate_element(0,0,0);
    m.assign_subblock(m,0,0);
    m.assign(vals);
    m.assign(m);
    m.assign_row(v,0);
    m.assign_column(v,0);
    m.accumulate(m);
    m.eigen(v,v,m,m);
    m.print("test");
    m.printValues();
    m.scale(0.0);
    m.toVector();
    m.nonnull();
    m.null();
    m.get_row(0);
    m.get_column(0);
    m.zero();
    m.nrow();
    m.ncol();
    vals = m.data();

    ConstRectMatrixPtr cm(m);
    cm.t();
    cm.get_subblock(0,0,0,0);
    cm.copy();
    cm.clone();
    cm.eigen(v,v,m,m);
    cm.print("test");
    cm.printValues();
    cm.toVector();
    cm.nonnull();
    cm.null();
    cm.get_row(0);
    cm.get_column(0);
    m.nrow();
    m.ncol();
    vals = m.data();


    s.i();
    s.sqrt_matrix();
    s.invsqrt_matrix();
    s.n();
    s.set_element(0,0,0);
    s.accumulate_element(0,0,0);
    s.assign(s);
    s.eigen(v,m);
    s.accumulate_symmetric_product(m);
    s.accumulate_transform(m,s);
    s.accumulate_transform(s,s);
    s.clone();
    s.copy();

    ConstSymmMatrixPtr cs(s);
    cs.i();
    cs.sqrt_matrix();
    cs.invsqrt_matrix();
    cs.n();
    cs.eigen(v,m);
    cs.clone();
    cs.copy();

    v[0];
    v.dot(v);
    v.get_element(0);
    v.set_element(0,0);
    v.accumulate_element(0,0);
    v.assign(vals);
    v.assign(v);
    v.printValues();
    v.print("test");
    v.null();
    v.nonnull();
    v.n();
    v.normalize();
    v.sort();
    v.accumulate(v);
    v.scale(0.5);
    v.assign(0.0);
    v.clone();
    v.copy();
    v.zero();
    v.maxabs();
    v.norm();
    v.symmetric_outer_product();
    v.data();
    v.nullcheck();
    v.sort();
    
    ConstVectorPtr cv(v);
    cv[0];
    cv.dot(v);
    cv.get_element(0);
    cv.printValues();
    cv.print("test");
    cv.null();
    cv.nonnull();
    cv.n();
    cv.clone();
    cv.copy();
    cv.maxabs();
    cv.norm();
    cv.symmetric_outer_product();
    cv.data();
    cv.nullcheck();

}
#endif


