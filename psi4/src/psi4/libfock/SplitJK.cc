/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include "SplitJK.h"

#include "psi4/libqt/qt.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
//#include "psi4/libqt/qt.h"
//#include "psi4/libmints/basisset.h"
//#include "psi4/libmints/mintshelper.h"

using namespace psi;

namespace psi {

SplitJK::SplitJK(std::shared_ptr<BasisSet> primary, Options& options) : primary_(primary), options_(options) {
    // set defaults
    print_ = 1;
    bench_ = 0;
    debug_ = 0;
    cutoff_ = 1.0E-12;
    
    // change defaults based on options
    if (options["PRINT"].has_changed()) bench_ = options_.get_int("PRINT");
    if (options["BENCH"].has_changed()) bench_ = options_.get_int("BENCH");
    if (options["DEBUG"].has_changed()) debug_ = options_.get_int("DEBUG");
    if (options["INTS_TOLERANCE"].has_changed()) cutoff_ = options_.get_double("INTS_TOLERANCE");
};

SplitJK::~SplitJK() {};

size_t SplitJK::num_computed_shells() {
    outfile->Printf("WARNING: SplitJK::num_computed_shells() was called, but benchmarking is disabled for the chosen SplitJK algorithm.");
    outfile->Printf(" Returning 0 as computed shells count.\n");

    return 0;
}

}
