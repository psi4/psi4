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

#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libdiis/diismanager.h"
#include "dcft.h"
#include "defines.h"

namespace psi{ namespace dcft{

SharedMatrix DCFTSolver::compute_gradient()
{
    // Print out the header
    outfile->Printf("\n\n\t***********************************************************************************\n");
    outfile->Printf(    "\t*                           DCFT Analytic Gradients Code                          *\n");
    outfile->Printf(    "\t*                by Alexander Sokolov, Andy Simmonett, and Xiao Wang              *\n");
    outfile->Printf(    "\t***********************************************************************************\n\n");

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
