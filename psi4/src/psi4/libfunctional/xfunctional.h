/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

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
    enum GGA_Type { GGA_None, B88, PBE, RPBE, SOGGA, PW91, B97, B86B, PW86}; 
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

    // B86B
    double _B86B_mu_;
    double _B86B_k_;

    // PW86
    double _PW86_m_;
    double _PW86_b_;
    double _PW86_c_;

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
