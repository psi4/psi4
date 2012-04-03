#include <libmints/vector.h>
#include "cfunctional.h"
#include "utility.h"
#include <psi4-dec.h>
#include <cmath>

using namespace psi;

namespace psi {

CFunctional::CFunctional()
{
    common_init();
}
CFunctional::~CFunctional()
{
}
void CFunctional::common_init()
{
    gga_type_ = GGA_None;
    meta_type_ = Meta_None;
    lsda_type_ = LSDA_None;
 
    _B97_ss_gamma_ = 0.0380;
    _B97_os_gamma_ = 0.0031;
}
void CFunctional::set_parameter(const std::string& key, double val) 
{
    parameters_[key] = val;

    if (key == "B97_ss_gamma") {
        _B97_ss_gamma_ = val;
    } else if (key.substr(0,8) == "B97_ss_a") {
        // B97_a0, B97_a1, etc
        int index = atoi(key.substr(8).c_str());
        if (_B97_ss_a_.size() < index + 1) {
            _B97_ss_a_.resize(index + 1);
            _B97_ss_a_[index] = val;
        }
    } else if (key == "B97_os_gamma") {
        _B97_os_gamma_ = val;
    } else if (key.substr(0,8) == "B97_os_a") {
        // B97_a0, B97_a1, etc
        int index = atoi(key.substr(8).c_str());
        if (_B97_os_a_.size() < index + 1) {
            _B97_os_a_.resize(index + 1);
            _B97_os_a_[index] = val;
        }
    } else {
        throw PSIEXCEPTION("Error, unknown generalized correlation functional parameter");    
    }
}
void CFunctional::compute_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha)
{
    compute_ss_functional(in,out,npoints,deriv,alpha,true);
    compute_ss_functional(in,out,npoints,deriv,alpha,false);
    compute_os_functional(in,out,npoints,deriv,alpha);
}
void CFunctional::compute_ss_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha, bool spin)
{
    if (deriv > 1) {
        throw PSIEXCEPTION("CFunctional: 2nd and higher partials not implemented yet.");
    }

    // Overall scale factor
    double A = alpha_ * alpha;

}
void CFunctional::compute_os_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha)
{
    if (deriv > 1) {
        throw PSIEXCEPTION("CFunctional: 2nd and higher partials not implemented yet.");
    }

    // Overall scale factor
    double A = alpha_ * alpha;

}

}
