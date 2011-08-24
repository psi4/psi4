#ifndef THREE_INDEX_QR
#define THREE_INDEX_QR

#include <psi4-dec.h>
#include <psiconfig.h>
#include <vector>

namespace psi {

class Matrix;
class Vector;

/*! 
 * First pass of QR decompostion for rank compression
 * Hermitian only, wastes some memory
 * Uses Householder reflections, as that sounds cooler
 */
class QR {

protected:
    /// Print flag (defaults to 0)
    int print_;
    /// Debug flag (defaults to 0)
    int debug_;

    /// Termination condition 
    double delta_;
    /// Original matrix (untouched) 
    boost::shared_ptr<Matrix> A_; 
    /// Q factor, partial 
    boost::shared_ptr<Matrix> Q_;
    /// R factor, full
    boost::shared_ptr<Matrix> R_;
    /// Pivots
    std::vector<int> pivots_; 
    /// P factor, partial
    boost::shared_ptr<Matrix> P_; 
    /// N factor, partial 
    boost::shared_ptr<Matrix> N_; 

    void form_QR();
    void form_PN();

public:
    QR(boost::shared_ptr<Matrix> A, double delta);
    ~QR();

    void decompose();

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }

    boost::shared_ptr<Matrix> A() const { return A_; }
    boost::shared_ptr<Matrix> Q() const { return Q_; }
    boost::shared_ptr<Matrix> R() const { return R_; }
    boost::shared_ptr<Matrix> P() const { return P_; }
    boost::shared_ptr<Matrix> N() const { return N_; }
    const std::vector<int>& pivots() const { return  pivots_; }

};

} // Namespace psi
#endif
