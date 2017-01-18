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

#include "psi4/psi4-dec.h"

#include "psi4/liboptions/liboptions.h"
#include "psi4/libciomr/libciomr.h"
#include "scf_grad.h"

namespace psi{
namespace scfgrad {

SharedMatrix scfgrad(SharedWavefunction ref_wfn, Options &options)
{
    tstart();

    SCFGrad grad(ref_wfn, options);
    SharedMatrix G = grad.compute_gradient();

    Process::environment.arrays["SCF TOTAL GRADIENT"] = G;
    Process::environment.arrays["CURRENT GRADIENT"] = G;
    Process::environment.set_gradient(G);

    tstop();
    return G;
}

SharedMatrix scfhess(SharedWavefunction ref_wfn, Options &options)
{
    tstart();

    SCFGrad grad(ref_wfn, options);
    SharedMatrix H = grad.compute_hessian();
    ref_wfn->set_hessian(H);

    tstop();

    return H;
}

}} // End Namespaces
