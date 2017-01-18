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

    // F_HJS^\omega(s,nu) kernel
    void hjs_F(double s, double nu, double* F, double* F_s, double* F_nu);

    // wB88? If so, re-scale s
    bool B88_;

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

    void set_B88(bool B88) { B88_ = B88; }
    bool B88() const { return B88_; }
};

}

#endif