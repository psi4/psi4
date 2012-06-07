#ifndef X_FUNCTIONAL_H
#define X_FUNCTIONAL_H

#include "functional.h"

namespace psi {

/** 
 * General exchange-type functional
 * 
 **/

class XFunctional : public Functional {

friend class Functional;

/**
* Fake polymorphic behavior
**/ 
public:
    enum GGA_Type { GGA_None, B88, PBE, RPBE, SOGGA, PW91, B97}; 
    enum Meta_Type { Meta_None, Becke };
    enum SR_Type { SR_None, LSDA, GGA, Meta }; 

protected:

    // => Enhancement factor types <= //
    GGA_Type gga_type_;
    Meta_Type meta_type_;
    SR_Type sr_type_;
    
    // => Specialized Parameters <= //
   
    // Global types (Slater, Thomas-Fermi constants)
    double _K0_;
    double _C0_; 
    double _k0_;
    double _pi12_;

    // > GGA < //

    // B/B3
    double _B88_d_;
    double _B88_a_;

    // PBE
    double _PBE_kp_;
    double _PBE_mu_;

    // PW91
    double _PW91_a1_;
    double _PW91_a2_;
    double _PW91_a3_;
    double _PW91_a4_;
    double _PW91_a5_;
    double _PW91_a6_;

    // B97
    double _B97_gamma_;
    std::vector<double> _B97_a_;

    // > Meta < //
    
    std::vector<double> _Meta_a_;

    // > SR < //

    // Set defaults up internally 
    
    void common_init();

public:

    // => Constructors (Use the factory constructor, or really know what's up) <= //

    XFunctional();
    virtual ~XFunctional(); 

    // => Parameters <= //
    
    virtual void set_parameter(const std::string& key, double val);

    // => Computers <= //

    virtual void compute_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha);

    void compute_sigma_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha, bool spin);
};

}

#endif
