/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "scf_grad.h"

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libscf_solver/rhf.h"
#include "psi4/libscf_solver/uhf.h"

namespace psi{
namespace scfgrad {

SharedMatrix scfgrad(std::shared_ptr<scf::HF> ref_wfn, Options &options)
{
    tstart();

    SCFDeriv grad(ref_wfn, options);
    SharedMatrix G = grad.compute_gradient();

    tstop();
    return G;
}

SharedMatrix scfhess(std::shared_ptr<scf::HF> ref_wfn, Options &options)
{
    tstart();

    SharedMatrix H;
    if( ref_wfn->same_a_b_orbs() && ref_wfn->same_a_b_dens()) {
        // RHF
        RSCFDeriv hessian_computer(std::dynamic_pointer_cast<scf::RHF>(ref_wfn), options);
        H = hessian_computer.compute_hessian();
    } else {
        USCFDeriv hessian_computer(std::dynamic_pointer_cast<scf::UHF>(ref_wfn), options);
        H = hessian_computer.compute_hessian();
    }
    ref_wfn->set_hessian(H);

    tstop();

    return H;
}

}} // End Namespaces
