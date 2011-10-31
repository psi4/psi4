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
    SharedMatrix A_; 
    /// Q factor, partial 
    SharedMatrix Q_;
    /// R factor, full
    SharedMatrix R_;
    /// Pivots
    std::vector<int> pivots_; 
    /// P factor, partial
    SharedMatrix P_; 
    /// N factor, partial 
    SharedMatrix N_; 

    void form_QR();
    void form_PN();

public:
    QR(SharedMatrix A, double delta);
    ~QR();

    void decompose();

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }

    SharedMatrix A() const { return A_; }
    SharedMatrix Q() const { return Q_; }
    SharedMatrix R() const { return R_; }
    SharedMatrix P() const { return P_; }
    SharedMatrix N() const { return N_; }
    const std::vector<int>& pivots() const { return  pivots_; }

};

} // Namespace psi
#endif
