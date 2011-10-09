#ifndef three_index_denominator_H
#define three_index_denominator_H

namespace psi {

class Matrix;
class Vector;

// Denominator Factorizations (MP2-like for now)
class Denominator {

protected:
    // Denominator (w in rows, ia in column)
    boost::shared_ptr<Matrix> denominator_;

    // Pointer to active occupied orbital eigenvalues
    boost::shared_ptr<Vector> eps_occ_;
    // Pointer to active virtual orbital eigenvalues
    boost::shared_ptr<Vector> eps_vir_;
    // Number of vectors required to obtain given accuracy
    int nvector_;
    // Maximum error norm allowed in denominator
    double delta_;

    virtual void decompose() = 0;
public:
    Denominator(boost::shared_ptr<Vector> eps_occ, boost::shared_ptr<Vector> eps_vir, double delta);
    virtual ~Denominator();

    // Factory method, algorithm should be LAPLACE or CHOLESKY
    static boost::shared_ptr<Denominator> buildDenominator(const std::string& algorithm,
            boost::shared_ptr<Vector> eps_occ, boost::shared_ptr<Vector> eps_vir, double delta);

    double delta() const { return delta_; }
    int nvector() const { return nvector_; }
    virtual void debug();
    boost::shared_ptr<Matrix> denominator() const { return denominator_; }

};

class LaplaceDenominator : public Denominator {

protected:
    // Fully split denominator (w in rows, i in columns)
    boost::shared_ptr<Matrix> denominator_occ_;
    // Fully split denominator (w in rows, a in columns)
    boost::shared_ptr<Matrix> denominator_vir_;

    void decompose();
public:
    LaplaceDenominator(boost::shared_ptr<Vector> eps_occ_, boost::shared_ptr<Vector> eps_vir, double delta);
    ~LaplaceDenominator();
    void debug();
    boost::shared_ptr<Matrix> denominator_occ() const { return denominator_occ_; }
    boost::shared_ptr<Matrix> denominator_vir() const { return denominator_vir_; }

};

class CholeskyDenominator : public Denominator {

protected:
    void decompose();
public:
    CholeskyDenominator(boost::shared_ptr<Vector> eps_occ_, boost::shared_ptr<Vector> eps_vir, double delta);
    ~CholeskyDenominator();
    void debug();

};

class SAPTDenominator {

protected:
    // Denominator (w in rows, ar in column) (monomer A)
    boost::shared_ptr<Matrix> denominatorA_;
    // Denominator (w in rows, bs in column) (monomer B)
    boost::shared_ptr<Matrix> denominatorB_;

    // Pointer to active occupied orbital eigenvalues (monomer A)
    boost::shared_ptr<Vector> eps_occA_;
    // Pointer to active virtual orbital eigenvalues (monomer A)
    boost::shared_ptr<Vector> eps_virA_;
    // Pointer to active occupied orbital eigenvalues (monomer B)
    boost::shared_ptr<Vector> eps_occB_;
    // Pointer to active virtual orbital eigenvalues (monomer B)
    boost::shared_ptr<Vector> eps_virB_;
    // Number of vectors required to obtain given accuracy
    int nvector_;
    // Maximum error norm allowed in denominator
    double delta_;
    // Crap all over the output file?
    bool debug_;

    virtual void decompose() = 0;
    void check_denom(boost::shared_ptr<Vector>, boost::shared_ptr<Vector>,
      boost::shared_ptr<Matrix>);
public:
    SAPTDenominator(boost::shared_ptr<Vector>, boost::shared_ptr<Vector>,
      boost::shared_ptr<Vector>, boost::shared_ptr<Vector>, double, bool);
    virtual ~SAPTDenominator();

    // Factory method, algorithm should be LAPLACE or CHOLESKY
    static boost::shared_ptr<SAPTDenominator> buildDenominator(const std::string& algorithm,
            boost::shared_ptr<Vector> eps_occA, boost::shared_ptr<Vector> eps_virA,
            boost::shared_ptr<Vector> eps_occB, boost::shared_ptr<Vector> eps_virB,
            double delta, bool debug = false);

    double delta() const { return delta_; }
    int nvector() const { return nvector_; }
    virtual void debug();
    boost::shared_ptr<Matrix> denominatorA() const { return denominatorA_; }
    boost::shared_ptr<Matrix> denominatorB() const { return denominatorB_; }

};

class SAPTLaplaceDenominator : public SAPTDenominator {

protected:
    // Fully split denominator (w in rows, a in columns) (monomer A)
    boost::shared_ptr<Matrix> denominator_occA_;
    // Fully split denominator (w in rows, r in columns) (monomer A)
    boost::shared_ptr<Matrix> denominator_virA_;
    // Fully split denominator (w in rows, b in columns) (monomer B)
    boost::shared_ptr<Matrix> denominator_occB_;
    // Fully split denominator (w in rows, s in columns) (monomer B)
    boost::shared_ptr<Matrix> denominator_virB_;

    void decompose();
    void check_split(boost::shared_ptr<Vector>, boost::shared_ptr<Vector>,
      boost::shared_ptr<Matrix>, boost::shared_ptr<Matrix>);
public:
    SAPTLaplaceDenominator(boost::shared_ptr<Vector>, boost::shared_ptr<Vector>,
      boost::shared_ptr<Vector>, boost::shared_ptr<Vector>, double, bool debug = false);
    ~SAPTLaplaceDenominator();

    void debug();
    boost::shared_ptr<Matrix> denominator_occA() const { return denominator_occA_; }
    boost::shared_ptr<Matrix> denominator_virA() const { return denominator_virA_; }
    boost::shared_ptr<Matrix> denominator_occB() const { return denominator_occB_; }
    boost::shared_ptr<Matrix> denominator_virB() const { return denominator_virB_; }

};

class SAPTCholeskyDenominator : public SAPTDenominator {

protected:
    void decompose();
public:
    SAPTCholeskyDenominator(boost::shared_ptr<Vector>, boost::shared_ptr<Vector>,
      boost::shared_ptr<Vector>, boost::shared_ptr<Vector>, double, bool debug = false);
    ~SAPTCholeskyDenominator();
};

} // Namespace psi
#endif
