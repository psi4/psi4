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

//
// irrep.cc
//
// Modifications are
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

/* irrep.cc -- implementation of the point group classes
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      June, 1993
 */

#include "psi4/psi4-dec.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <cstdlib>

using namespace psi;

/////////////////////////////////////////////////////////////////////////

IrreducibleRepresentation::IrreducibleRepresentation()
    : g(0), degen(0), nrot_(0), ntrans_(0), complex_(0), symb(nullptr), csymb(nullptr) {}

IrreducibleRepresentation::IrreducibleRepresentation(int order, int d, const char *lab, const char *clab)
    : g(0), degen(0), nrot_(0), ntrans_(0), complex_(0), symb(nullptr), csymb(nullptr) {
    init(order, d, lab, clab);
}

IrreducibleRepresentation::IrreducibleRepresentation(const IrreducibleRepresentation &ir)
    : g(0), degen(0), nrot_(0), ntrans_(0), complex_(0), symb(nullptr), csymb(nullptr) {
    *this = ir;
}

IrreducibleRepresentation::~IrreducibleRepresentation() { init(); }

IrreducibleRepresentation &IrreducibleRepresentation::operator=(const IrreducibleRepresentation &ir) {
    init(ir.g, ir.degen, ir.symb, ir.csymb);

    nrot_ = ir.nrot_;
    ntrans_ = ir.ntrans_;
    complex_ = ir.complex_;

    for (int i = 0; i < g; i++) rep[i] = ir.rep[i];

    return *this;
}

void IrreducibleRepresentation::init(int order, int d, const char *lab, const char *clab) {
    g = order;
    degen = d;
    ntrans_ = nrot_ = complex_ = 0;

    free(symb);
    if (lab)
        symb = strdup(lab);
    else
        symb = nullptr;

    free(csymb);
    if (clab)
        csymb = strdup(clab);
    else
        csymb = nullptr;

    if (rep.size()) {
        rep.clear();
    }

    if (g) {
        rep = std::vector<SymRep>(g);
        for (int i = 0; i < g; i++) rep[i].set_dim(d);
    }
}

void IrreducibleRepresentation::print(std::string out) const {
    if (!g) return;
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    int i, d;

    printer->Printf("  %-5s", symb);

    for (i = 0; i < g; i++) printer->Printf(" %6.3f", character(i));
    printer->Printf(" | %d t, %d R\n", ntrans_, nrot_);

    for (d = 0; d < nproj(); d++) {
        printer->Printf("       ");
        for (i = 0; i < g; i++) printer->Printf(" %6.3f", p(d, i));
        printer->Printf("\n");
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
