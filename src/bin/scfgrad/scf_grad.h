/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef SCF_GRAD_H
#define SCF_GRAD_H

#include <libmints/wavefunction.h>
#include <libmints/typedefs.h>

namespace psi {

namespace scfgrad {

class SCFGrad : public Wavefunction {

protected:

    /// Common initialization
    void common_init();
    
public:
    SCFGrad();
    virtual ~SCFGrad();
    
    double compute_energy() { throw PSIEXCEPTION("SCFGrad needs a rehash, call Rob."); }
   
    SharedMatrix compute_gradient(); 

    SharedMatrix compute_hessian();

    SharedMatrix rhf_hessian_response();
};

}} // Namespaces

#endif
