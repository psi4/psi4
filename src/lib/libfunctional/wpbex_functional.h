#ifndef wPBEX_FUNCTIONAL_H
#define wPBEX_FUNCTIONAL_H

#include "functional.h"

namespace psi {

void wpbe_F(double rho, double s, double omega, double* F, double* F_rho, double* F_s);


/** 
 * Short-range PBE functional 
 * 
 **/

class wPBEXFunctional : public Functional {

protected:

    // => Specialized Parameters <= //
   
    // Global types (Slater, Thomas-Fermi constants)
    double _K0_;
    double _k0_;

    // Set defaults up internally 
    
    void common_init();

public:

    // => Constructors (Use the factory constructor, or really know what's up) <= //

    wPBEXFunctional();
    virtual ~wPBEXFunctional(); 

    // => Parameters <= //
    
    virtual void set_parameter(const std::string& key, double val);

    // => Computers <= //

    virtual void computeRKSFunctional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha);
    virtual void computeUKSFunctional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha);

    void computeFunctional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha, bool in_alpha, bool out_alpha);
    

};

}

#endif
