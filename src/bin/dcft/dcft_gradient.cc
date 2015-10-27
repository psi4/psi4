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

#include <libtrans/integraltransform.h>
#include <libpsio/psio.hpp>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <libdiis/diismanager.h>
#include "dcft.h"
#include "defines.h"

namespace psi{ namespace dcft{

SharedMatrix DCFTSolver::compute_gradient()
{
    // Print out the header
    outfile->Printf("\n\n\t***********************************************************************************\n");
    outfile->Printf(    "\t*                           DCFT Analytic Gradients Code                          *\n");
    outfile->Printf(    "\t*                by Alexander Sokolov, Andy Simmonett, and Xiao Wang              *\n");
    outfile->Printf(    "\t***********************************************************************************\n");

    // If the system is closed-shell, then ...
    if(options_.get_str("REFERENCE") == "RHF"){
        compute_gradient_RHF();
    }
    // If the system is open-shell, then ...
    else{
        compute_gradient_UHF();
    }

    return SharedMatrix(new Matrix("NULL", 0, 0));
}

}} //End namespaces
