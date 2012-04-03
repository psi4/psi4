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
