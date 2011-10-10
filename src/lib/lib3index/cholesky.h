#ifndef THREE_INDEX_CHOLESKY
#define THREE_INDEX_CHOLESKY

namespace psi {

class Matrix;
class Vector;
class TwoBodyAOInt;

class Cholesky {

protected:
    /// Maximum Chebyshev error allowed in the decomposition
    double delta_;
    /// Maximum memory to use, in doubles
    unsigned long int memory_;
    /// Full L (Q x n), if choleskify() called()
    boost::shared_ptr<Matrix> L_;
    /// Number of columns required, if choleskify() called
    int Q_;

public:
    /*!
     * Constructor, does not build decomposition.
     * \param delta: maximum Chebyshev error allowed in the decomposition
     * \param memory: maximum memory allowed, in doubles
     **/
    Cholesky(double delta, unsigned long int memory);
    /// Destructor, resets L_
    ~Cholesky();

    /// Perform the cholesky decomposition (requires 2QN memory)
    virtual void choleskify();

    /// Shared pointer to decomposition (Q x N), if choleskify() called
    boost::shared_ptr<Matrix> L() const { return L_; }
    /// Number of columns required to reach accuracy delta, if choleskify() called
    int Q() const { return Q_; }
    /// Dimension of the original square tensor, provided by the subclass
    virtual int N() = 0;
    /// Maximum Chebyshev error allowed in the decomposition
    double delta() const { return delta_; }

    /// Diagonal of the original square tensor, provided by the subclass
    virtual void compute_diagonal(double* target) = 0;
    /// Row row of the original square tensor, provided by the subclass
    virtual void compute_row(int row, double* target) = 0;
};

class CholeskyMatrix : public Cholesky {

protected:
    boost::shared_ptr<Matrix> A_;
public:
    CholeskyMatrix(boost::shared_ptr<Matrix> A, double delta, unsigned long int memory);
    ~CholeskyMatrix();

    virtual int N();
    virtual void compute_diagonal(double* target);
    virtual void compute_row(int row, double* target);
};

class CholeskyERI : public Cholesky {

protected:
    double schwarz_;
    boost::shared_ptr<BasisSet> basisset_;
    boost::shared_ptr<TwoBodyAOInt> integral_;
public:
    CholeskyERI(boost::shared_ptr<TwoBodyAOInt> integral, double schwarz, double delta, unsigned long int memory);
    ~CholeskyERI();

    virtual int N();
    virtual void compute_diagonal(double* target);
    virtual void compute_row(int row, double* target);
};

class CholeskyMP2 : public Cholesky {

protected:
    bool symmetric_;
    boost::shared_ptr<Matrix> Qia_;
    boost::shared_ptr<Vector> eps_aocc_;
    boost::shared_ptr<Vector> eps_avir_;
public:
    CholeskyMP2(boost::shared_ptr<Matrix> Qia, boost::shared_ptr<Vector> eps_aocc,
        boost::shared_ptr<Vector> eps_avir, bool symmetric,
        double delta, unsigned long int memory);
    ~CholeskyMP2();

    virtual int N();
    virtual void compute_diagonal(double* target);
    virtual void compute_row(int row, double* target);
};

class CholeskyDelta : public Cholesky {

protected:
    boost::shared_ptr<Vector> eps_aocc_;
    boost::shared_ptr<Vector> eps_avir_;
public:
    CholeskyDelta(boost::shared_ptr<Vector> eps_aocc,
        boost::shared_ptr<Vector> eps_avir,
        double delta, unsigned long int memory);
    ~CholeskyDelta();

    virtual int N();
    virtual void compute_diagonal(double* target);
    virtual void compute_row(int row, double* target);
};

class CholeskyLocal : public Cholesky {

protected:
    boost::shared_ptr<Matrix> C_;
public:
    CholeskyLocal(boost::shared_ptr<Matrix> C,
        double delta, unsigned long int memory);
    ~CholeskyLocal();

    virtual int N();
    virtual void compute_diagonal(double* target);
    virtual void compute_row(int row, double* target);
};

} // Namespace psi
#endif
