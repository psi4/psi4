#ifndef _psi_src_lib_libmints_vector_h
#define _psi_src_lib_libmints_vector_h

#include <cstdlib>
#include <cstdio>
#include <psi4-dec.h>

#include <libparallel/serialize.h>

namespace boost {
template<class T> class shared_ptr;

namespace python {
class tuple;
}}

namespace psi {

class Dimension;
class Matrix;

/*! \ingroup MINTS */
class Vector : public Serializable {
protected:
    /// Vector data
    double **vector_;
    /// Number of irreps
    int nirrep_;
    /// Dimensions per irrep
    int *dimpi_;
    /// Name
    std::string name_;

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
    /// Constructor, allocates memory
    Vector(const std::string& name, int nirrep, int *dimpi);
    /// Constructor, convenience for 1 irrep
    Vector(const std::string& name, int dim);
    /// Constructor, takes Dimension object
    Vector(const Dimension& dimpi);
    /// Constructor, takes Dimension object
    Vector(const std::string& name, const Dimension& dimpi);

    /// Destructor, frees memory
    ~Vector();

    void init(int nirrep, int *dimpi);
    void init(int nirrep, const int *dimpi, const std::string& name = "");
    void init(const Dimension& v);

    Vector* clone();

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

    /// Returns a single element value
    double get(int m) {
        return vector_[0][m];
    }
    /// Sets a single element value
    void set(int m, double val) {
        vector_[0][m] = val;
    }

    void add(int m, double val) {
        vector_[0][m] += val;
    }

    void add(int h, int m, double val) {
        vector_[h][m] += val;
    }

    /// Adds other vector to this
    void add(const boost::shared_ptr<Vector>& other) {
        for (int h=0; h<nirrep_; ++h) {
            for (int m=0; m<dimpi_[h]; ++m) {
                vector_[h][m] += other->vector_[h][m];
            }
        }
    }

    /// Adds other vector to this
    void add(const Vector& other) {
        for (int h=0; h<nirrep_; ++h) {
            for (int m=0; m<dimpi_[h]; ++m) {
                vector_[h][m] += other.vector_[h][m];
            }
        }
    }

    double& operator()(int i) { return vector_[0][i]; }
    const double& operator()(int i) const { return vector_[0][i]; }

    double pyget(const boost::python::tuple& key);
    void pyset(const boost::python::tuple& key, double value);
    double pyget(int key);
    void pyset(int key, double value);

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

    /**
     * Sets the name of the vector, used in print(...)
     *
     * @param name New name to use.
     */
    void set_name(const std::string& name) {
        name_ = name;
    }

    /**
     * Gets the name of the matrix.
     */
    std::string name() const {
        return name_;
    }

    /// Python compatible printer
    void print_out() { print(outfile); }
    /**
     * Print the matrix using print_mat
     *
     * @param outfile File point to use, defaults to Psi4's outfile.
     * @param extra When printing the name of the 'extra' will be printing after the name.
     */
    void print(FILE *out = outfile, const char *extra=NULL) const;

    /// Copies rhs to this
    void copy(const Vector* rhs);
    /// Copies rhs to this
    void copy(const Vector& rhs);

    /// General matrix vector multiplication
    void gemv(bool transa, double alpha, Matrix* A, Vector* X, double beta);

    /// Vector dot product
    double dot(Vector* X);

    /// Scale the elements of the vector
    void scale(double sc);

    // Serializable pure virtual functions:
    void send();
    void recv();
    void bcast(int broadcaster);
    /**
     * Performs element-by-element sum of all data from all nodes.
     */
    void sum();

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

    /// Other vector initializer
    void init(int dim);

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

    SimpleVector& operator+=(const SimpleVector& b) {
        for (int i=0; i<dim_; ++i)
            vector_[i] += b.vector_[i];
        return *this;
    }

    /// Computes the magnitude of the vector.
    double magnitude() const;

    /// Prints the vector
    void print() const;

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
    /// Name of the IntVector
    std::string name_;

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
    /// Constructor, allocates memory
    IntVector(const std::string& name, int nirrep, int *dimpi);
    /// Constructor, convenience for 1 irrep
    IntVector(const std::string& name, int dim);

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

    /**
     * Sets the name of the vector, used in print(...)
     *
     * @param name New name to use.
     */
    void set_name(const std::string& name) {
        name_ = name;
    }

    /**
     * Gets the name of the matrix.
     */
    std::string name() const {
        return name_;
    }
    /// Python compatible printer
    void print_out() { print(outfile); }

    /**
     * Print the matrix using print_mat
     *
     * @param outfile File point to use, defaults to Psi4's outfile.
     * @param extra When printing the name of the 'extra' will be printing after the name.
     */
    void print(FILE *out = outfile, const char *extra=NULL) const;
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

    /// Other vector initializer
    void init(int dim);

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
