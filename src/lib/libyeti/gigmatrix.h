#ifndef yeti_gigmatrix_h
#define yeti_gigmatrix_h

#include <libsmartptr/serialize.h>
#include <iostream>
#include <math.h>
#include "gigmatrix.hpp"
#include "class.h"

namespace yeti {


/**
    @class Matrix

    Class encapsulating a 2-D matrix
*/
class Matrix : 
    public smartptr::Serializable 
{

    protected:
        double* data_;
        uli nrow_;
        uli ncol_;
    
        void swap_column(uli i, uli j, double* tmp);


    public:
        /**
            @param nrow
            @param ncol
        */
        Matrix(uli nrow, uli ncol);

        /**
            Assign pointer constructor.  This does NOT copy values.
            This directly assigns the pointer for efficients reasons.

            @param vals The pointer to assign to the matrix
            @param nrow
            @param ncol
        */
        Matrix(double* vals, uli nrow, uli ncol);

        /**
            Take the values from a vector and treat them as a matrix.  This copies values.
            @param v
            @param nrow
            @param ncol
        */
        Matrix(Vector* v, uli nrow, uli ncol);

        Matrix(const XMLArchivePtr& arch);

        void serialize(const XMLArchivePtr& arch) const;

        ~Matrix();

        double get_element(uli i, uli j) const;

        void set_element(uli i, uli j, double val);

        void accumulate_element(uli i, uli j, double val);

        void canonicalize();

        /**
            Scale all elements of the matrix
            @param s The scale factor
        */
        void scale(double s);

        /**
            Assign values to the matrix.  This copies values
            @param vals The values to copy assign to the matrix
        */
        void assign(const double* vals);

        /**
            Assign a uniform value to all elements of the matrix
            @param s The value to assign to each matrix element
        */
        void assign(double s);

        /**
            @return The linear data array that stores the matrix
        */
        const double* data() const;

        void printValues(std::ostream& os = std::cout) const;

        /**
            @param title An informative title for the matrix
            @param os The output stream to print to
        */
        void print(const std::string& title, std::ostream& os = std::cout) const;

        void print(std::ostream& os = std::cout) const;

        /**
            Add the values from the parameter matrix into the current values
            and return a new matrix
            @param r The matrix to add (right matrix in summation)
            @return The sum matrix
        */
        Matrix* add(const Matrix* r) const;

        /**
            Subtract the values from the parameter matrix into the current values
            and return a new matrix
            @param r The matrix to subtract (right matrix in summation)
            @return The sum matrix
        */
        Matrix* subtract(const Matrix* r) const;

        /**
            Multiply all values of the current matrix and return a new matrix
            @param d The scale factor for the matrix
            @return The product matrix
        */
        Matrix* mult(double d) const;

        /**
            Multiply this matrix by the paramater matrix and return a new matrix
            @param m The right matrix in the multiplication
            @return The product matrix
        */
        Matrix* mult(const Matrix* m) const;

        /**
            Multiply this matrix by the paramater matrix and return a new vector
            @param m The right vector in the multiplication
            @return The product vector
        */
        Vector* mult(const Vector* m) const;

        /**
            @return The number of rows in the matrix
        */
        uli nrow() const;

        /**
            @return The number of cols in the matrix
        */
        uli ncol() const;

        /**
            Accumulate element by element the values of the parameter
            matrix into the current matrix
            @param m The matrix to accumulate
        */
        void accumulate(const Matrix* m);

        /**
            Transpose the current matrix and return a new matrix
            @return The transpose matrix
        */
        Matrix* t() const;

        /**
            Create a 1-D vector from the matrix. This copies values!
            @return A vector with the same values as the matrix
        */
        Vector* toVector() const;

};

/**
    @class Vector

    Encapsulates a vector
*/
class Vector : 
    public smartptr::Serializable 
{

    protected:
        double* data_;
        uli n_;

    public:
        /**
            Create a vector of size n and assign pointer.  This does NOT
            copy values. It directly assigns the pointer.
            @param d
            @param n
        */
        Vector(double* d, uli n);

        /**
            Create a zero vector of size n
            @param n
        */
        Vector(uli n);

        Vector(const XMLArchivePtr& arch);

        void serialize(const XMLArchivePtr& arch) const;

        ~Vector();

        /**
            Add the values from the parameter vector into the current vector
            and return new vector
            @param r The right vector in sum
            @return Sum vector
        */
        Vector* add(const Vector* r) const;

        /**
            Subtract the values from the parameter vector from the current vector
            and return new vector
            @param r The right vector in the subtraction
            @return Difference vector
        */
        Vector* subtract(const Vector* r) const;

        /**
            Multiply the current vector by the scale factor and return a new vector
            @param d The scale factor
            @return The scaled vector
        */
        Vector* mult(double d) const;

        /**
            @return The number of elements in the vector
        */
        uli n() const;

        /**
            @return The dot product of two vectors
        */
        double dot(const Vector* v) const;

        /**
            @param v The vector to accumulate into the current vector
        */
        void accumulate(const Vector* v);

        void sort();

        /**
            @param s The scale factor
        */
        void scale(double s);

        /**
            Copy values from array into current vector
            @param vals The values to copy into the current vector
        */
        void assign(const double* vals);

        /**
            Assign a uniform value to all elements in the vector
            @param s The value to assign to all vector elements
        */
        void assign(double s);

        void print(const std::string& title, std::ostream& os = std::cout) const;

        void print(std::ostream& os = std::cout) const;

        void printValues(std::ostream& os = std::cout) const;

        double get_element(uli i) const;

        void set_element(uli i, double val);

        void accumulate_element(uli i, double val);

        /**
            @return The data array that stores the vector values
        */
        const double* data() const;

        /**
            Reinterpret the vector as a matrix and return a matrix
            which holds a copy of the vector values
            @param nrow
            @param ncol
            @return A matrix interpretation of the vector
        */
        Matrix* toMatrix(uli nrow, uli ncol) const;

};

template <class T>
class MatrixTemplate : 
    public boost::intrusive_ptr<T> 
{
    
    protected:
        typedef boost::intrusive_ptr<T> Parent;

    public:
        using Parent::get;
    
    public:
        /**
            Create a matrix with a null pointer
        */
        MatrixTemplate();

        /**
            Create matrix smart pointer for the given matrix
            @param m The matrix pointer to store a refcount object for
        */
        MatrixTemplate(T* m);

        /**
            @return Whether the matrix pointer is null
        */
        bool null() const;

        /**
            @return Whether the matrix pointer is non-null
        */
        bool nonnull() const;

        void scale(double s);

        double get_element(uli i, uli j) const;

        void printValues(std::ostream& os = std::cout) const;

        void print(const std::string& title, std::ostream& os = std::cout) const;

        void zero();

        uli nrow() const;

        uli ncol() const;

        const double* data() const;

        VectorPtr get_column(uli col) const;

        VectorPtr get_row(uli row) const;

        double maxabs() const;

        VectorPtr toVector() const;

        /**
            Ensure that the matrix is non-null before doing matrix operations.
            This aborts on null.
        */
        void nullcheck() const;

};

/**
    @class RectMatrixTemplate

    Template class for rectangular matrices.  The template parameter is essentially
    for whether the underlying matrix is const or not.
*/
template <class T>
class RectMatrixTemplate : public MatrixTemplate<T> {
    
    private:
        typedef MatrixTemplate<T> Parent;

    public:
        using Parent::get;
        using Parent::nullcheck;

        RectMatrixTemplate();

        ~RectMatrixTemplate();

        RectMatrixTemplate(T* m);

        /**
            Create and allocate a zero matrix with the given dimensions
            @param nrow
            @param ncol
        */
        RectMatrixTemplate(uli nrow, uli ncol);

        /**
            Copy constructor
        */
        RectMatrixTemplate(const RectMatrixPtr& m);

        /**
            @return The transpose matrix
        */
        RectMatrixPtr t() const;

        /**
            Return a subblock of the current matrix
            @param rowstart
            @param rowstop Inclusive index
            @param colstart
            @param colstop Inclusive index
        */
        RectMatrixPtr get_subblock(uli rowstart, uli rowstop, uli colstart, uli colstop) const;

        /**
            @return An exact copy of the given matrix
        */
        RectMatrixPtr copy() const;

        /**
            @return A zero matrix with the same dimensions
        */
        RectMatrixPtr clone() const;

        void set_element(uli i, uli j, double val);

        void accumulate_element(uli i, uli j, double val);

        void assign_subblock(const RectMatrixPtr& block, uli rowstart, uli colstart);

        void assign(const double* vals);

        void assign(const RectMatrixPtr& m);

        void assign_row(const VectorPtr& v, uli row);

        void assign_column(const VectorPtr& v, uli col);

        void accumulate(const RectMatrixPtr&);

        /**
            @param revals Reference return of the real part of the eigenvalues. The vector
                           can be null and will be properly allocated.
            @param ievels Reference return of the imaginary part of the eigenvalues. The vector
                            can be null and will be properly allocated.
            @param Levecs Reference return of the left eigenvectors. The matrix can be null
                            and will be properly allocated.
            @param Revecs Reference return of the right eigenvectors. The matrix can be null
                            and will be properly allocated.
        */
        void eigen(VectorPtr& revals, VectorPtr& ievals, RectMatrixPtr& Levecs, RectMatrixPtr& Revecs) const;

};

/**
    @class SymmMatrixTemplate

    Template class for symmetric matrices.  The template parameter is essentially
    for whether the underlying matrix is const or not.
*/
template <class T>
class SymmMatrixTemplate : public MatrixTemplate<T> {
    
    private:
        typedef MatrixTemplate<T> Parent;
    
    public:
        using Parent::get;
        using Parent::nullcheck;

        SymmMatrixTemplate();

        SymmMatrixTemplate(uli n);

        SymmMatrixTemplate(T* m);

        SymmMatrixTemplate(const SymmMatrixPtr& m);

        SymmMatrixPtr i(double tol = 1e-8) const;

        SymmMatrixPtr sqrt_matrix(double tol = 1e-8) const;

        SymmMatrixPtr invsqrt_matrix(double tol = 1e-8) const;

        uli n() const;

        void set_element(uli i, uli j, double val);

        void accumulate_element(uli i, uli j, double val);

        void accumulate(const SymmMatrixPtr&);

        void assign(const SymmMatrixPtr& m);

        void eigen(VectorPtr& evals, RectMatrixPtr& evecs) const;

        void accumulate_symmetric_product(const RectMatrixPtr& m);

        void accumulate_transform(const RectMatrixPtr& t, const SymmMatrixPtr& m);

        void accumulate_transform(const SymmMatrixPtr& t, const SymmMatrixPtr& m);

        SymmMatrixPtr clone() const;

        SymmMatrixPtr copy() const;

};

/**
    @class VectorTemplate

    Template class for vectors.  The template parameter is essentially
    for whether the underlying vector is const or not.
*/
template <class T>
class VectorTemplate : 
    public boost::intrusive_ptr<T> 
{

    protected:
        typedef boost::intrusive_ptr<T> Parent;

    public:
        using Parent::get;

        VectorTemplate(); 

        VectorTemplate(uli n);

        VectorTemplate(T* v);

        VectorTemplate(const VectorPtr& v);

        double operator[](uli i) const;

        double dot(const VectorPtr& v) const;

        double get_element(uli i) const;

        void set_element(uli i, double val);

        void accumulate_element(uli i, double val);

        void assign(const double* vals);

        void assign(const VectorPtr& v);

        void print(const std::string& title, std::ostream& os = std::cout) const;

        void printValues(std::ostream& os = std::cout) const;

        bool null() const;

        bool nonnull() const;

        uli n() const;

        void normalize();

        void sort();

        void accumulate(const VectorPtr& v);

        void scale(double s);

        void assign(double d);

        VectorPtr clone() const;

        VectorPtr copy() const;

        void zero();

        double maxabs() const;

        double norm() const;

        SymmMatrixPtr symmetric_outer_product() const;

        const double* data() const;

        void nullcheck() const;

};

/**
    @class RectMatrixPtr

    Non-template class for easier use of rectangular matrices
*/
class RectMatrixPtr : 
    public RectMatrixTemplate<Matrix> 
{

    private:
        typedef RectMatrixTemplate<Matrix> Parent;

    public:
        using Parent::get;

        /**
            Create a null pointer rectangular matrix
        */
        RectMatrixPtr();

        RectMatrixPtr(Matrix* m);

        RectMatrixPtr(uli nrow, uli ncol);

        RectMatrixPtr(const RectMatrixPtr& m);

        RectMatrixPtr(const boost::intrusive_ptr<Matrix>& m);

};

/**
    @class SymmMatrixPtr

    Non-template class for easier use of symmetric matrices
*/
class SymmMatrixPtr : 
    public SymmMatrixTemplate<Matrix> 
{
    
    private:
        typedef SymmMatrixTemplate<Matrix> Parent;

    public:
        using Parent::get;
        using Parent::nullcheck;

        SymmMatrixPtr();

        SymmMatrixPtr(Matrix* m);

        SymmMatrixPtr(uli n);

        SymmMatrixPtr(const SymmMatrixPtr& m);

        SymmMatrixPtr(const boost::intrusive_ptr<Matrix>& m);

};

/**
    @class VectorPtr
    Non-template class for easier use of vectors
*/
class VectorPtr : public VectorTemplate<Vector> {

    private:
        typedef VectorTemplate<Vector> Parent;

    public:
        VectorPtr(); 

        VectorPtr(uli n);

        VectorPtr(Vector* v);

        VectorPtr(const VectorPtr& v);

        VectorPtr(const boost::intrusive_ptr<Vector>& v);

};

class Vector1 : public VectorPtr {

    public:
        Vector1(double a0);

};

class Vector2 : public VectorPtr {

    public:
        Vector2(double a0, double a1);

};

class Vector3 : public VectorPtr {

    public:
        Vector3(double a0, double a1, double a2);

};

class Vector4 : public VectorPtr {

    public:
        Vector4(double a0, double a1, double a2, double a3);

};

class Vector8 : public VectorPtr {

    public:
        Vector8(double a0, double a1, double a2, double a3, 
                double a4, double a5, double a6, double a7);

};


RectMatrixPtr operator*(const RectMatrixPtr& l, const RectMatrixPtr& r);

RectMatrixPtr operator*(const RectMatrixPtr& l, const SymmMatrixPtr& r);

RectMatrixPtr operator*(const SymmMatrixPtr& l, const RectMatrixPtr& r);

RectMatrixPtr operator*(const SymmMatrixPtr& l, const SymmMatrixPtr& r);

RectMatrixPtr operator*(const RectMatrixPtr&, double d);

RectMatrixPtr operator*(double d, const RectMatrixPtr&);

RectMatrixPtr operator+(const RectMatrixPtr& l, const RectMatrixPtr& r);

RectMatrixPtr operator-(const RectMatrixPtr& l, const RectMatrixPtr& r);

SymmMatrixPtr operator+(const SymmMatrixPtr& l, const SymmMatrixPtr& r);

SymmMatrixPtr operator-(const SymmMatrixPtr& l, const SymmMatrixPtr& r);

SymmMatrixPtr operator*(const SymmMatrixPtr&, double d);

SymmMatrixPtr operator*(double d, const SymmMatrixPtr&);

VectorPtr operator*(const RectMatrixPtr& m, const VectorPtr& v);

VectorPtr operator*(const SymmMatrixPtr& m, const VectorPtr& v);

VectorPtr operator*(double d, const VectorPtr& v);

VectorPtr operator*(const VectorPtr&, double d);

VectorPtr operator+(const VectorPtr& l, const VectorPtr& r);

VectorPtr operator-(const VectorPtr& l, const VectorPtr& r);

VectorPtr cross(const VectorPtr& l, const VectorPtr& r);

bool equals(const RectMatrixPtr& l, const RectMatrixPtr& r, double tol = 1e-8);

bool equals(const VectorPtr& l, const VectorPtr& r, double tol = 1e-8);

void instantiate();

void
error(
    const Matrix* l,
    const Matrix* r,
    const std::string& op
);

void
error(
    const Vector* l,
    const Vector* r,
    const std::string& op
);

void
scale_array(
    double* v,
    double s,
    uli n
);

void
assign_array(
    double* v,
    double d,
    uli n
);

double 
max_array(
    const double* v,
    uli n
);

void
transpose(
    double* vals,
    uli nrow,
    uli ncol
);

void
eigenvalues(
    const Matrix* matrix,
    Vector* evals,
    Matrix* evecs
);

void
eigenvalues(
    const Matrix* matrix,
    Vector* rvals,
    Vector* ivals,
    Matrix* Levecs,
    Matrix* Revecs
);

double*
subtract_arrays(
    const double* l,
    const double* r,
    uli n
);

double*
add_arrays(
    const double* l,
    const double* r,
    uli n
);

void
accumulate_array(
    double* target,
    const double* src,
    uli n
);

Matrix*
multiply(
    const Matrix* l,
    const Matrix* r,
    bool transpose_l = false,
    bool transpose_r = false
);

template <class T> 
MatrixTemplate<T>::MatrixTemplate()
    : Parent(0)
{
}

template <class T> 
MatrixTemplate<T>::MatrixTemplate(T* m)
    : Parent(m)
{
}

template <class T> 
bool
MatrixTemplate<T>::null() const
{
    return !get();
}

template <class T> 
bool
MatrixTemplate<T>::nonnull() const
{
    return get();
}

template <class T> 
void
MatrixTemplate<T>::scale(double s)
{
    get()->scale(s);
}

template <class T> 
void
MatrixTemplate<T>::zero()
{
    get()->assign(0.0);
}

template <class T> 
double
MatrixTemplate<T>::get_element(uli i, uli j) const
{
    return get()->get_element(i,j);
}

template <class T> 
VectorPtr
MatrixTemplate<T>::toVector() const
{
    Vector* v = get()->toVector();
    return v;
}

template <class T> 
const double*
MatrixTemplate<T>::data() const
{
    return get()->data();
}

template <class T> 
void
MatrixTemplate<T>::printValues(std::ostream& os) const
{
    nullcheck();
    get()->printValues(os);
}

template <class T> 
void
MatrixTemplate<T>::print(const std::string& title, std::ostream& os) const
{
    if (get() == 0)
        os << "null matrix" << std::endl;
    else
        get()->print(title, os);
}

template <class T> 
void
MatrixTemplate<T>::nullcheck() const
{
    if (null())
    {
        std::cerr << "Called method on null matrix" << std::endl;
        abort();
    }
}

template <class T> 
double
MatrixTemplate<T>::maxabs() const
{
    nullcheck();
    uli nrow = get()->nrow();
    uli ncol = get()->ncol();
    return max_array(data(), nrow * ncol);
}

template <class T> 
VectorPtr
MatrixTemplate<T>::get_row(uli row) const
{
    nullcheck();
    uli ncol = get()->ncol();
    Vector* v = new Vector(ncol);
    for (uli col=0; col < ncol; ++col)
        v->set_element(col, get_element(row, col));
    return v;
}

template <class T> 
VectorPtr
MatrixTemplate<T>::get_column(uli col) const
{
    nullcheck();
    uli nrow = get()->nrow();
    Vector* v = new Vector(nrow);
    for (uli row=0; row < nrow; ++row)
        v->set_element(row, get_element(row, col));
    return v;
}


template <class T> 
uli
MatrixTemplate<T>::nrow() const
{
    nullcheck();
    return get()->nrow();
}


template <class T> 
uli
MatrixTemplate<T>::ncol() const
{
    nullcheck();
    return get()->ncol();
}

template <class T> 
RectMatrixTemplate<T>::RectMatrixTemplate()
    : Parent(0)
{
}

template <class T> 
RectMatrixTemplate<T>::~RectMatrixTemplate()
{
}

template <class T> 
RectMatrixTemplate<T>::RectMatrixTemplate(T* m)
    : Parent(m)
{
}

template <class T> 
RectMatrixTemplate<T>::RectMatrixTemplate(uli nrow, uli ncol)
    : Parent(new Matrix(nrow, ncol))
{
}

template <class T> 
RectMatrixTemplate<T>::RectMatrixTemplate(const RectMatrixPtr& m)
    : Parent(m.get())
{
}

template <class T> 
void
RectMatrixTemplate<T>::accumulate_element(uli i, uli j, double val)
{
    nullcheck();
    get()->accumulate_element(i,j,val);
}

template <class T>
void
RectMatrixTemplate<T>::set_element(uli i, uli j, double val)
{
    nullcheck();
    get()->set_element(i,j,val);
}

template <class T> 
void
RectMatrixTemplate<T>::accumulate(const RectMatrixPtr& m)
{
    nullcheck();
    get()->accumulate(m.get());
}

template <class T> 
void
RectMatrixTemplate<T>::assign(const RectMatrixPtr& m)
{
    nullcheck();
    get()->assign(m.data());
}

template <class T> 
void
RectMatrixTemplate<T>::assign(const double* vals)
{
    nullcheck();
    get()->assign(vals);
}


template <class T> 
RectMatrixPtr
RectMatrixTemplate<T>::t() const
{
    nullcheck();
    return get()->t();
}

template <class T> 
void
RectMatrixTemplate<T>::assign_row(const VectorPtr& v, uli row)
{
    nullcheck();
    for (uli col=0; col < Parent::ncol(); ++col)
        set_element(row, col, v[col]);
}

template <class T> 
void
RectMatrixTemplate<T>::assign_column(const VectorPtr& v, uli col)
{
    nullcheck();
    for (uli row=0; row < Parent::nrow(); ++row)
        set_element(row, col, v[row]);
}

template <class T> 
RectMatrixPtr
RectMatrixTemplate<T>::clone() const
{
    nullcheck();
    return new Matrix(Parent::nrow(), Parent::ncol());
}

template <class T> 
RectMatrixPtr
RectMatrixTemplate<T>::copy() const
{
    nullcheck();
    Matrix* m = new Matrix(Parent::nrow(), Parent::ncol());
    m->assign(Parent::data());
    return m;
}

template <class T> 
RectMatrixPtr
RectMatrixTemplate<T>::get_subblock(uli rowstart, uli rowstop, uli colstart, uli colstop) const
{
    nullcheck();
    uli nr = rowstop - rowstart + 1;
    uli nc = colstop - colstart + 1;
    Matrix* ptr = new Matrix(nr, nc);
    uli r = 0;
    for (uli row=rowstart; row <= rowstop; ++row, ++r)
    {
        uli c = 0;
        for (uli col=colstart; col <= colstop; ++col, ++c)
        {
            ptr->set_element(r, c, Parent::get_element(row, col)); 
        }
    }
    return ptr;
}

template <class T> 
void
RectMatrixTemplate<T>::assign_subblock(const RectMatrixPtr& block, uli rowstart, uli colstart)
{
    nullcheck();
    uli nr = block.nrow();
    uli nc = block.ncol();
    uli row, r, col, c;
    for (r=0, row=rowstart; r < nr; ++r, ++row)
    {
        for (c=0, col=colstart; c < nc; ++c, ++col)
        {
            set_element(row, col, block.get_element(r,c));
        }
    }
}

template <class T> 
void
RectMatrixTemplate<T>::eigen(VectorPtr& revals, VectorPtr& ievals, RectMatrixPtr& Levecs, RectMatrixPtr& Revecs) const
{
    nullcheck();

    if (revals.null()) revals = new Vector(Parent::nrow());
    if (ievals.null()) ievals = new Vector(Parent::nrow());
    if (Levecs.null()) Levecs = new Matrix(Parent::nrow(), Parent::nrow());
    if (Revecs.null()) Revecs = new Matrix(Parent::nrow(), Parent::nrow());

    eigenvalues(get(), revals.get(), ievals.get(), Levecs.get(), Revecs.get());
}

template <class T>
VectorTemplate<T>::VectorTemplate(uli n)
    : Parent(new Vector(n))
{
}

template <class T>
VectorTemplate<T>::VectorTemplate()
    : Parent(0)
{
}

template <class T>
VectorTemplate<T>::VectorTemplate(T* v)
    : Parent(v)
{
}

template <class T>
VectorTemplate<T>::VectorTemplate(const VectorPtr& v)
    : Parent(v.get())
{
}

template <class T>
uli
VectorTemplate<T>::n() const
{
    nullcheck();
    return get()->n();
}

template <class T>
double
VectorTemplate<T>::dot(const VectorPtr& v) const
{
    nullcheck();
    return get()->dot(v.get());
}

template <class T>
void
VectorTemplate<T>::nullcheck() const
{
    if (null())
    {
        std::cerr << "Called method on null matrix" << std::endl;
        abort();
    }
}

template <class T>
void
VectorTemplate<T>::printValues(std::ostream& os) const
{
    nullcheck();
    get()->printValues(os);
}

template <class T>
void
VectorTemplate<T>::print(const std::string& title, std::ostream& os) const
{
    nullcheck();
    get()->print(title, os);
}

template <class T>
void
VectorTemplate<T>::assign(const VectorPtr& v)
{
    nullcheck();
    get()->assign(v.data());
}

template <class T>
void
VectorTemplate<T>::assign(const double* vals)
{
    nullcheck();
    get()->assign(vals);
}

template <class T>
const double*
VectorTemplate<T>::data() const
{
    nullcheck();
    return get()->data();
}

template <class T>
double
VectorTemplate<T>::operator[](uli i) const
{
    nullcheck();
    return get()->get_element(i);
}

template <class T>
double
VectorTemplate<T>::get_element(uli i) const
{
    nullcheck();
    return get()->get_element(i);
}

template <class T>
void
VectorTemplate<T>::set_element(uli i, double val)
{
    nullcheck();
    get()->set_element(i, val);
}

template <class T>
void
VectorTemplate<T>::accumulate_element(uli i, double val)
{
    nullcheck();
    get()->accumulate_element(i, val);
}

template <class T>
void
VectorTemplate<T>::sort()
{
    get()->sort();
}

template <class T>
bool
VectorTemplate<T>::null() const
{
    return !get();
}

template <class T>
bool
VectorTemplate<T>::nonnull() const
{
    return get();
}

template <class T>
void
VectorTemplate<T>::zero()
{
    nullcheck();
    get()->assign(0.0);
}

template <class T>
double
VectorTemplate<T>::norm() const
{
    nullcheck();
    double n2 = get()->dot(get());
    return sqrt(n2);
}

template <class T>
void
VectorTemplate<T>::accumulate(const VectorPtr& v)
{
    nullcheck();
    v.nullcheck();
    get()->accumulate(v.get());
}

template <class T>
void
VectorTemplate<T>::scale(double d)
{
    nullcheck();
    get()->scale(d);
}

template <class T>
void
VectorTemplate<T>::assign(double d)
{
    nullcheck();
    get()->assign(d);
}

template <class T>
void
VectorTemplate<T>::normalize()
{
    nullcheck();
    double normsq = get()->dot(get()); 
    double oonorm = 1.0/sqrt(normsq);
    scale(oonorm);
}

template <class T>
VectorPtr
VectorTemplate<T>::copy() const
{
    nullcheck();
    Vector* v = new Vector(n());
    v->assign(data());
    return v;
}

template <class T>
VectorPtr
VectorTemplate<T>::clone() const
{
    nullcheck();
    Vector* v = new Vector(n());
    return v;
}

template <class T>
double
VectorTemplate<T>::maxabs() const
{
    nullcheck();
    return max_array(data(), n());
}

template <class T>
SymmMatrixPtr
VectorTemplate<T>::symmetric_outer_product() const
{
    nullcheck();
    const double* v = data();
    Matrix* m = new Matrix(n(), n());
    for (uli i=0; i < n(); ++i)
    {
        for (uli j=0; j < n(); ++j)
        {
            m->set_element(i, j, v[i] * v[j]);
        }
    }
    return m;
}

template <class T> 
SymmMatrixTemplate<T>::SymmMatrixTemplate()
    : Parent(0)
{
}

template <class T> 
SymmMatrixTemplate<T>::SymmMatrixTemplate(uli n)
    : Parent(new Matrix(n,n))
{
}

template <class T> 
SymmMatrixTemplate<T>::SymmMatrixTemplate(T* m)
    : Parent(m)
{
}

template <class T> 
SymmMatrixTemplate<T>::SymmMatrixTemplate(const SymmMatrixPtr& m)
    : Parent(m.get())
{
    m.nullcheck();
}

template <class T> 
uli
SymmMatrixTemplate<T>::n() const
{
    nullcheck();
    return get()->nrow();
}

template <class T> 
void
SymmMatrixTemplate<T>::set_element(uli i, uli j, double val)
{
    nullcheck();
    get()->set_element(i,j,val);
    get()->set_element(j,i,val);
}

template <class T> 
void
SymmMatrixTemplate<T>::accumulate_element(uli i, uli j, double val)
{
    nullcheck();
    get()->accumulate_element(i, j, val);
}

template <class T> 
void
SymmMatrixTemplate<T>::accumulate(const SymmMatrixPtr& m)
{
    nullcheck();
    get()->accumulate(m.get());
}

template <class T> 
SymmMatrixPtr
SymmMatrixTemplate<T>::clone() const
{
    nullcheck();
    return new Matrix(n(), n());
}

template <class T> 
SymmMatrixPtr
SymmMatrixTemplate<T>::copy() const
{
    nullcheck();
    Matrix* m = new Matrix(n(), n());
    m->assign(Parent::data());
    return m;
}

template <class T> 
void
SymmMatrixTemplate<T>::assign(const SymmMatrixPtr& m)
{
    get()->assign(m->data());
}

template <class T> 
void
SymmMatrixTemplate<T>::eigen(VectorPtr& evals, RectMatrixPtr& evecs) const
{
    nullcheck();

    if (evals.null()) evals = new Vector(n());
    if (evecs.null()) evecs = new Matrix(n(), n());

    eigenvalues(get(), evals.get(), evecs.get());
}

template <class T> 
SymmMatrixPtr
SymmMatrixTemplate<T>::sqrt_matrix(double tol) const
{
    nullcheck();
    RectMatrixPtr evecs;
    VectorPtr evals;
    SymmMatrixPtr epsilon(n());

    eigen(evals, evecs);
    for (uli i=0; i < n(); ++i)
    {
        double val = evals[i];
        if ( fabs(val) > tol)
            epsilon.set_element(i, i, sqrt(val));
    }

    SymmMatrixPtr sqrtmat = clone();
    sqrtmat.accumulate_transform(evecs, epsilon);

    return sqrtmat;
}

template <class T> 
SymmMatrixPtr
SymmMatrixTemplate<T>::invsqrt_matrix(double tol) const
{
    nullcheck();
    RectMatrixPtr evecs;
    VectorPtr evals;
    SymmMatrixPtr epsilon(n());

    eigen(evals, evecs);
    for (uli i=0; i < n(); ++i)
    {
        double val = evals[i];
        if ( fabs(val) > tol)
            epsilon.set_element(i, i, 1.0/sqrt(val));
    }

    SymmMatrixPtr sqrtmat = clone();
    sqrtmat.accumulate_transform(evecs, epsilon);

    return sqrtmat;
}

template <class T> 
SymmMatrixPtr
SymmMatrixTemplate<T>::i(double tol) const
{
    nullcheck();
    RectMatrixPtr evecs;
    VectorPtr evals;
    SymmMatrixPtr epsilon(n());

    eigen(evals, evecs);
    for (uli i=0; i < n(); ++i)
    {
        double val = evals[i];
        if ( fabs(val) > tol )
            epsilon.set_element(i, i, 1.0/val);
    }

    SymmMatrixPtr inv = clone();
    inv.accumulate_transform(evecs, epsilon);

    return inv;
}

template <class T> 
void
SymmMatrixTemplate<T>::accumulate_transform(const RectMatrixPtr& r, const SymmMatrixPtr& s)
{
    nullcheck();
    r.nullcheck();
    s.nullcheck();

    Matrix* temp = multiply(r.get(), s.get());
    Matrix* final = multiply(temp, r.get(), false, true);
    get()->accumulate(final);
    delete temp;
    delete final;
}

template <class T> 
void
SymmMatrixTemplate<T>::accumulate_transform(const SymmMatrixPtr& r, const SymmMatrixPtr& s)
{
    nullcheck();
    r.nullcheck();
    s.nullcheck();
    Matrix* temp = multiply(r.get(), s.get());

    Matrix* final = multiply(temp, r.get());
    get()->accumulate(final);
    delete temp;
    delete final;
}

template <class T> 
void
SymmMatrixTemplate<T>::accumulate_symmetric_product(const RectMatrixPtr& m)
{
    nullcheck();
    m.nullcheck();

    Matrix* prod = multiply(m.get(), m.get(), false, true); //multiply by transpose
    get()->accumulate(prod);
    delete prod;
}


}

#if 0
namespace smartptr {
    SerialDecideSubptr(yeti::RectMatrixPtr);
    SerialDecideSubptr(yeti::SymmMatrixPtr);
    SerialDecideSubptr(yeti::VectorPtr);
}
#endif

#endif
