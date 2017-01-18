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

#ifndef wPBEC_FUNCTIONAL_H
#define wPBEC_FUNCTIONAL_H

#include "functional.h"

namespace psi {

/** 
 * Short-range PBE correlation functional 
 **/

class wPBECFunctional : public Functional {

public:
    enum wPBEC_Type { null_wPBEC_type, pw92c_type, pw92c_sr_type, pbec_type, pbec_sr_type };

protected:

    // => Fake Polymorphism <= //

    wPBEC_Type type_;

    // => Utility Computers <= //

    // Usual PW92 epsilon
    void pw92c_eps(
        double rho, 
        double z, 
        double* eps, 
        double* eps_rho, 
        double* eps_z);
    
    // Short-range PW92 epsilon
    void pw92c_sr_eps(
        double omega, 
        double rho, 
        double z, 
        double* eps,
        double* eps_rho, 
        double* eps_z,
        double* eps_sr,
        double* eps_sr_rho, 
        double* eps_sr_z);

    // PW92C functional
    void pw92c_f(
        double rho,
        double z,
        double* f,  
        double* f_rho,
        double* f_z);

    // PBE functional
    void pbec_f(
        double rho,
        double z,
        double s,
        double* f,  
        double* f_rho,
        double* f_z,
        double* f_s);

    // Short-Range PW92C functional
    void pw92c_sr_f(
        double omega,
        double rho,
        double z,
        double* f,  
        double* f_rho,
        double* f_z);
    
    // PBE functional
    void pbec_sr_f(
        double omega,
        double rho,
        double z,
        double s,
        double* f,  
        double* f_rho,
        double* f_z,
        double* f_s);

    // Set defaults up internally 
    void common_init();

public:

    // => Constructors (Use the factory constructor, or really know what's up) <= //

    wPBECFunctional();
    virtual ~wPBECFunctional(); 

    // => Computers <= //

    virtual void compute_functional(const std::map<std::string,SharedVector>& in, const std::map<std::string,SharedVector>& out, int npoints, int deriv, double alpha);

    void set_wPBEC_type(wPBEC_Type type) { type_ = type; common_init(); }
};

}

#endif