#ifndef C_FUNCTIONAL_H
#define C_FUNCTIONAL_H

#include "functional.h"

namespace psi {

/** 
 * General correlation-type functional
 **/

class CFunctional : public Functional {

friend class Functional;

/**
* Fake polymorphic behavior
**/ 
public:
    enum LSDA_Type { LSDA_None, PW92}; 
    enum GGA_Type  { GGA_None, B97}; 
    enum Meta_Type { Meta_None, B95};

protected:

    // => Enhancement factor types <= //
    LSDA_Type lsda_type_;
    GGA_Type  gga_type_;
    Meta_Type meta_type_;
    
    // => Specialized Parameters <= //

    // PW92 Parameters
    double _c0_;
    double _two13_;
    double _d2fz0_;

    double _c0a_;
    double _a1a_;
    double _b1a_;
    double _b2a_;
    double _b3a_;
    double _b4a_;

    double _c0p_;
    double _a1p_;
    double _b1p_;
    double _b2p_;
    double _b3p_;
    double _b4p_;
   
    double _c0f_;
    double _a1f_;
    double _b1f_;
    double _b2f_;
    double _b3f_;
    double _b4f_;

    void G1(double r, double* G, double* Gr);
    void G2(double r, double* G, double* Gr);
    void G3(double r, double* G, double* Gr);
   
    // B97
    double _B97_ss_gamma_;
    std::vector<double> _B97_ss_a_;
    double _B97_os_gamma_;
    std::vector<double> _B97_os_a_;

    // Set defaults up internally 
    
    void common_init();

public:

    // => Constructors (Use the factory constructor, or really know what's up) <= //

    CFunctional();
    virtual ~CFunctional(); 

    // => Parameters <= //
    
    virtual void set_parameter(const std::string& key, double val);

    // => Computers <= //

    virtual void compute_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha);

    void compute_ss_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha, bool spin);
    void compute_os_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha);
    

};

}

#endif
