#ifndef _psi_src_lib_libmints_dealias_h_
#define _psi_src_lib_libmints_dealias_h_

#include <cstdio>
#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace psi {

class DealiasBasisSet {

protected:
    /// The primary basis set
    boost::shared_ptr<BasisSet> primary_;
    /// The effective alphas of the primary set [center][am][index]
    std::vector<std::vector<std::vector<double> > > primary_alpha_;
    /// The alphas of the dealias set [center][am][index]
    std::vector<std::vector<std::vector<double> > > dealias_alpha_;
    /// The resultant dealias set
    boost::shared_ptr<BasisSet> dealias_;

    // => Parameters <= //
    /// Base for core and diffuse functions (~2)
    double delta_;
    /// Base for even-tempered cap functions (~3.5)
    double beta_;
    /// Number of intercalaters per window (~1)
    int nintercalater_;
    /// Number of core functions per block (~1)
    int ncore_;
    /// Number of diffuse functions per block (~1)
    int ndiffuse_;
    /// Number of cap functions (~1)
    int ncap_;

    /// Helper functions
    void form_primary_alpha();
    void form_core();
    void form_diffuse();
    void form_intercalater();
    void form_cap();
    void form_basis();

public:
    DealiasBasisSet(boost::shared_ptr<BasisSet> primary_);
    virtual ~DealiasBasisSet();

    /// Parameter entry
    void setDelta(double delta) { delta_ = delta; }
    void setBeta(double beta) { beta_ = beta; }
    void setNCore(double n) { ncore_ = n; }
    void setNCap(double n) { ncap_ = n; }
    void setNIntercalater(double n) { nintercalater_ = n; }
    void setNDiffuse(double n) { ndiffuse_ = n; }
   
    /// Master build routine 
    boost::shared_ptr<BasisSet> buildDealiasBasisSet();

    /// Convenience routine
    boost::shared_ptr<BasisSet> dealiasSet() const { return dealias_; }
};

}

#endif
