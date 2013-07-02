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

#include "psi4-dec.h"
#include <libmints/mints.h>
#include <liboptions/liboptions.h>
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include "scf_grad.h"

namespace psi{ 
namespace scfgrad {

PsiReturnType scfgrad(Options &options)
{
    tstart();

    boost::shared_ptr<SCFGrad> grad(new SCFGrad());
    SharedMatrix G = grad->compute_gradient();

    Process::environment.set_gradient(G); 
    Process::environment.wavefunction()->set_gradient(G);

    tstop();

    return Success;
}

PsiReturnType scfhess(Options &options)
{
    tstart();

    boost::shared_ptr<SCFGrad> grad(new SCFGrad());
    SharedMatrix G = grad->compute_hessian();

    tstop();

    return Success;
}

}} // End Namespaces
