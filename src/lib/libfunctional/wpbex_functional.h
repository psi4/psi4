#ifndef wPBEX_FUNCTIONAL_H
#define wPBEX_FUNCTIONAL_H

#include "functional.h"

namespace psi {

void wpbe_F(double rho, double s, double omega, double* F, double* F_rho, double* F_s);


/** 
 * Short-range PBE functional 
 * Following HJS notation of J. Chem. Phys., 128, 194105
 **/

class wPBEXFunctional : public Functional {

protected:

    // => Specialized Parameters <= //
   
    // Global types (Slater, Thomas-Fermi constants)
    double _K0_;
    double _k0_;
    double _s0_;
    double _pi12_;

    double _s_min_tol_;
    double _nu_min_tol_;

    // HJS Parameters 
    double _A_;
    double _B_;
    double _C_;
    double _D_;
    double _E_;
    
    std::vector<double> _Ha_;
    std::vector<double> _Hb_;

    // HJS or HSE?
    bool hjs_;

    // F_HJS^\omega(s,nu) kernel
    void hjs_F(double s, double nu, double* F, double* F_s, double* F_nu);
    // F_HSE^\omega(s,nu) kernel
    void hse_F(double rho, double s, double omega, double* F, double* F_rho, double* F_s);

    // Set defaults up internally 
    
    void common_init();

public:

    // => Constructors (Use the factory constructor, or really know what's up) <= //

    wPBEXFunctional();
    virtual ~wPBEXFunctional(); 

    // => Parameters <= //
    
    virtual void set_parameter(const std::string& key, double val);

    // => Computers <= //

    virtual void compute_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha);
    void compute_sigma_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha, bool spin);

    bool hjs() const { return hjs_; }
    void set_hjs(bool hjs) { hjs_ = hjs; }

};

}

#endif
