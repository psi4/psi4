#ifndef THREE_INDEX_FITTER
#define THREE_INDEX_FITTER

#include <psi4-dec.h>
#include <psiconfig.h>

namespace psi {

class BasisSet;
class Matrix;
class Vector;

/*! 
 * Small utility class to compute d_A = J_AB^{-1} (Q|mn) D_mn coefficients for QM/MM
 */
class DFChargeFitter {

protected:
    /// Print flag (defaults to 1)
    int print_;
    /// Debug flag (defaults to 0)
    int debug_;

    /// Target coefficients
    SharedVector d_;
    /// Driving density
    SharedMatrix D_;
    
    /// Primary Basis Set
    boost::shared_ptr<BasisSet> primary_;
    /// Auxiliary Basis Set
    boost::shared_ptr<BasisSet> auxiliary_;

public:
    
    DFChargeFitter();
    ~DFChargeFitter();

    SharedVector fit();

    void setD(SharedMatrix D) { D_ = D; }
    void setPrimary(boost::shared_ptr<BasisSet> primary) { primary_ = primary; }
    void setAuxiliary(boost::shared_ptr<BasisSet> auxiliary) { auxiliary_ = auxiliary; }

    SharedVector d() const { return d_; }

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }
};

} // Namespace psi
#endif
