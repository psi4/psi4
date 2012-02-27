#ifndef _psi_src_lib_libmints_vector_h
#define _psi_src_lib_libmints_vector_h

#include <cstdlib>
#include <cstdio>
#include <psi4-dec.h>

#include <libparallel/serialize.h>

#include <vector>
#include <iterator>

namespace boost {
template<class T> class shared_ptr;

namespace python {
class tuple;
}}

namespace psi {

class Dimension;
class Matrix;
class VectorIterator;

/*! \ingroup MINTS */
class Vector : public Serializable {
protected:
    /// Actual data, of size dimpi_.sum()
    std::vector<double> v_;
    /// Pointer offsets into v_, of size dimpi_.n()
    std::vector<double* > vector_;
    /// Number of irreps
    int nirrep_;
    /// Dimensions per irrep
    Dimension dimpi_;
    /// Name
    std::string name_;

    /// Allocates vector_
    void alloc();
    /// Releases vector_
    void release();

    /// Copies data to vector_
    void copy_from(const Vector& other);

    /// Assign pointer offsets in vector_ from v_.
    void assign_pointer_offsets();

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

    /// Zeros the vector out
    void zero();

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
    double& operator[](int i) { return vector_[0][i]; }
    const double& operator[](int i) const { return vector_[0][i]; }

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
    void scale(const double& sc);

    // Serializable pure virtual functions:
    void send();
    void recv();
    void bcast(int broadcaster);
    /**
     * Performs element-by-element sum of all data from all nodes.
     */
    void sum();

    typedef std::vector<double>::iterator       iterator;
    typedef std::vector<double>::const_iterator const_iterator;

    /// @{
    /** Returns the starting iterator for the entire v_. */
    iterator begin()
        { return v_.begin(); }
    const_iterator begin() const
        { return v_.begin(); }
    /// @}

    /// @{
    /** Returns the ending iterator for the entire v_. */
    iterator end()
        { return v_.end(); }
    const_iterator end() const
        { return v_.end(); }
    /// @}

    /// @{
    /** Returns the starting iterator for irrep h. */
    iterator begin_irrep(int h)
        { return iterator(vector_[h]); }
    const_iterator begin_irrep(int h) const
        { return const_iterator(vector_[h]); }
    /// @}

    /// @{
    /** Returns the starting iterator for irrep h. */
    iterator end_irrep(int h)
        { return iterator(vector_[h]) + dimpi_[h]; }
    const_iterator end_irrep(int h) const
        { return const_iterator(vector_[h]) + dimpi_[h]; }
    /// @}

    friend class Matrix;
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

    friend class VectorIterator;
};

//class VectorIterator : public std::iterator<std::forward_iterator_tag, double>
//{
//    pointer v_;

//public:
//    VectorIterator(pointer v) : v_(v) {}

//    VectorIterator(const VectorIterator& mit) : v_(mit.v_) {}

//    VectorIterator& operator++() {v_++; return *this;}  // prefix (++a)
//    VectorIterator operator++(int) { VectorIterator tmp(*this); operator++(); return tmp; } // suffix (a++)

//    bool operator==(const VectorIterator& rhs) { return v_ == rhs.v_; }
//    bool operator!=(const VectorIterator& rhs) { return v_ != rhs.v_; }

//    reference operator*() { return *v_; }
//};

typedef boost::shared_ptr<Vector> SharedVector;
typedef boost::shared_ptr<IntVector> SharedIntVector;

}

#endif
