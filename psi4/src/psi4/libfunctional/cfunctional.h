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
    enum GGA_Type  { GGA_None, B97, PBE}; 
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

    void PW92_C(double rho, double z, double* PW92, double* PW92_rho, double* PW92_z);

    // PBE
    double _bet_;

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