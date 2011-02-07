#ifndef _psi_src_lib_libmints_vector_h
#define _psi_src_lib_libmints_vector_h

#include <cstdlib>
#include <cstdio>
#include <psi4-dec.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class Matrix;

/*! \ingroup MINTS */
class Vector {
protected:
    /// Vector data
    double **vector_;
    /// Number of irreps
    int nirrep_;
    /// Dimensions per irrep
    int *dimpi_;

    /// Allocates vector_
    void alloc();
    /// Releases vector_
    void release();

    /// Copies data to vector_
    void copy_from(double **);

public:
    /// Default constructor, zeros everything out
    Vector();
    /// Copy constructor
    Vector(const Vector& copy);
    /// Constructor, allocates memory
    Vector(int nirrep, int *dimpi);
    /// Constructor, convenience for 1 irrep
    Vector(int dim);

    /// Destructor, frees memory
    ~Vector();

    void init(int nirrep, int *dimpi);

    /// Sets the vector_ to the data in vec
    void set(double *vec);

    /// Returns a pointer to irrep h
    double* pointer(int h = 0) {
        return vector_[h];
    }

    /// Returns a single element value
    double get(int h, int m) {
        return vector_[h][m];
    }
    /// Sets a single element value
    void set(int h, int m, double val) {
        vector_[h][m] = val;
    }

    /// Returns a copy of the vector_
    double *to_block_vector();

    /// Returns the dimension per irrep h
    int dim(int h = 0) const {
        return dimpi_[h];
    }

    /// Returns the dimension array
    int *dimpi() const {
        return dimpi_;
    }
    /// Returns the number of irreps
    int nirrep() const {
        return nirrep_;
    }

    /// Python compatible printer
    void print_out() { print(outfile); }

    /// Prints the vector
    void print(FILE *);
    /// Copies rhs to this
    void copy(const Vector* rhs);
    /// Copies rhs to this
    void copy(const Vector& rhs);

    /// General matrix vector multiplication
    void gemv(bool transa, double alpha, Matrix* A, Vector* X, double beta);

    /// Vector dot product
    double dot(Vector* X);

    friend class Matrix;
};

/*! \ingroup MINTS */
class SimpleVector
{
protected:
    /// Vector data
    double *vector_;
    /// Dimension of the vector
    int dim_;

    /// Allocate memory
    void alloc();
    /// Free memory
    void release();

    /// Copy data to this
    void copy_from(double *);

public:
    /// Default constructor, zeroes everything out
    SimpleVector();
    /// Copy constructor
    SimpleVector(const SimpleVector& copy);
    /// Constructor, creates the vector
    SimpleVector(int dim);

    /// Destructor, frees memory
    ~SimpleVector();

    /// Set vector_ to vec
    void set(double *vec);
    /// Returns a pointer to the vector's contents
    double* pointer() {
        return vector_;
    }
    /// Returns an element value
    double get(int m) {
        return vector_[m];
    }
    /// Sets an element value
    void set(int m, double val) {
        vector_[m] = val;
    }
    /// Returns a copy of vector_
    double *to_block_vector();

    /// Returns the dimension of the vector
    int dim() const {
        return dim_;
    }

    double& operator[](int i) { return vector_[i]; }

    void operator=(const SimpleVector& x) {
        for (int i=0; i<dim_; ++i)
            vector_[i] = x.vector_[i];
    }

    /// Python compatible printer
    void print_out() { print(outfile); }

    /// Prints the vector
    void print(FILE *);
    /// Copy rhs to this
    void copy(const SimpleVector* rhs);

    /// Scale the vector
    void scale(double a);

    friend class SimpleMatrix;
};
/*! \ingroup MINTS */
class IntVector {
protected:
    /// IntVector data
    int **vector_;
    /// Number of irreps
    int nirrep_;
    /// Dimensions per irrep
    int *dimpi_;

    /// Allocates vector_
    void alloc();
    /// Releases vector_
    void release();

    /// Copies data to vector_
    void copy_from(int **);

public:
    /// Default constructor, zeros everything out
    IntVector();
    /// Copy constructor
    IntVector(const IntVector& copy);
    /// Constructor, allocates memory
    IntVector(int nirrep, int *dimpi);
    /// Constructor, convenience for 1 irrep
    IntVector(int dim);

    /// Destructor, frees memory
    ~IntVector();

    void init(int nirrep, int *dimpi);

    /// Sets the vector_ to the data in vec
    void set(int *vec);

    /// Returns a pointer to irrep h
    int* pointer(int h = 0) {
        return vector_[h];
    }

    /// Returns a single element value
    int get(int h, int m) {
        return vector_[h][m];
    }
    /// Sets a single element value
    void set(int h, int m, int val) {
        vector_[h][m] = val;
    }
    /// Sets a single element value
    void set_python(int h, int m, int val) {
        vector_[h][m] = val;
    }

    /// Returns a copy of the vector_
    int *to_block_vector();

    /// Returns the dimension per irrep h
    int dim(int h = 0) const {
        return dimpi_[h];
    }

    /// Returns the dimension array
    int *dimpi() const {
        return dimpi_;
    }
    /// Returns the number of irreps
    int nirrep() const {
        return nirrep_;
    }

    /// Python compatible printer
    void print_out() { print(outfile); }

    /// Prints the vector
    void print(FILE *);
    /// Copies rhs to this
    void copy(const IntVector* rhs);
    /// Copies rhs to this
    void copy(const IntVector& rhs);

};

/*! \ingroup MINTS */
class SimpleIntVector
{
protected:
    /// IntVector data
    int *vector_;
    /// Dimension of the vector
    int dim_;

    /// Allocate memory
    void alloc();
    /// Free memory
    void release();

    /// Copy data to this
    void copy_from(int *);

public:
    /// Default constructor, zeroes everything out
    SimpleIntVector();
    /// Copy constructor
    SimpleIntVector(const SimpleIntVector& copy);
    /// Constructor, creates the vector
    SimpleIntVector(int dim);

    /// Destructor, frees memory
    ~SimpleIntVector();

    /// Set vector_ to vec
    void set(int *vec);
    /// Returns a pointer to the vector's contents
    int* pointer() {
        return vector_;
    }
    /// Returns an element value
    int get(int m) {
        return vector_[m];
    }
    /// Sets an element value
    void set(int m, int val) {
        vector_[m] = val;
    }
    /// Returns a copy of vector_
    int *to_block_vector();

    /// Returns the dimension of the vector
    int dim() const {
        return dim_;
    }

    int& operator[](int i) { return vector_[i]; }

    void operator=(const SimpleIntVector& x) {
        for (int i=0; i<dim_; ++i)
            vector_[i] = x.vector_[i];
    }

    /// Python compatible printer
    void print_out() { print(outfile); }

    /// Prints the vector
    void print(FILE *);
    /// Copy rhs to this
    void copy(const SimpleIntVector* rhs);

};

typedef boost::shared_ptr<Vector> SharedVector;
typedef boost::shared_ptr<SimpleVector> SharedSimpleVector;
typedef boost::shared_ptr<IntVector> SharedIntVector;
typedef boost::shared_ptr<SimpleIntVector> SharedSimpleIntVector;

}

#endif
