/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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
// chartab.cc
//
// Additional modifications made by Justin Turney <jturney@ccqc.uga.edu>
// for use in PSI4.
//
// Modifictions are
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

/* chartab.cc -- implementation of the point group classes
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

#include "psi4/libmints/pointgrp.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/exception.h"

#include <cctype>

using namespace psi;

////////////////////////////////////////////////////////////////////////

CharacterTable::CharacterTable()
    : nt(0), pg(PointGroups::C1), nirrep_(0), gamma_(nullptr), symop(nullptr), _inv(nullptr), symb(), bits_(0) {}

CharacterTable::CharacterTable(const CharacterTable& ct)
    : nt(0), pg(PointGroups::C1), nirrep_(0), gamma_(nullptr), symop(nullptr), _inv(nullptr), symb(), bits_(0) {
    *this = ct;
}

CharacterTable::~CharacterTable() {
    if (gamma_) delete[] gamma_;
    gamma_ = nullptr;
    if (symop) delete[] symop;
    symop = nullptr;
    if (_inv) delete[] _inv;
    _inv = nullptr;
    nt = nirrep_ = 0;
}

CharacterTable& CharacterTable::operator=(const CharacterTable& ct) {
    nt = ct.nt;
    pg = ct.pg;
    nirrep_ = ct.nirrep_;

    symb = ct.symb;

    if (gamma_) delete[] gamma_;
    gamma_ = nullptr;
    if (ct.gamma_) {
        gamma_ = new IrreducibleRepresentation[nirrep_];
        for (int i = 0; i < nirrep_; i++) {
            gamma_[i].init();
            gamma_[i] = ct.gamma_[i];
        }
    }

    if (symop) delete[] symop;
    symop = nullptr;

    if (ct.symop) {
        symop = new SymmetryOperation[nirrep_];
        for (int i = 0; i < nirrep_; i++) {
            symop[i] = ct.symop[i];
        }
    }

    if (_inv) delete[] _inv;
    _inv = nullptr;

    if (ct._inv) {
        _inv = new int[nirrep_];
        memcpy(_inv, ct._inv, sizeof(int) * nirrep_);
    }

    bits_ = ct.bits_;

    return *this;
}

void CharacterTable::print(std::string out) const {
    if (!nirrep_) return;
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    int i;

    printer->Printf("  point group %s\n\n", symb.c_str());

    for (i = 0; i < nirrep_; i++) gamma_[i].print(out);

    printer->Printf("\n  symmetry operation matrices:\n\n");
    for (i = 0; i < nirrep_; i++) symop[i].print(out);

    printer->Printf("\n  inverse symmetry operation matrices:\n\n");
    for (i = 0; i < nirrep_; i++) symop[inverse(i)].print(out);
}

void CharacterTable::common_init() {
    // first parse the point group symbol, this will give us the order of the
    // point group(g), the type of point group (pg), the order of the principle
    // rotation axis (nt), and the number of irreps (nirrep_)

    if (!symb.length()) {
        // ExEnv::errn() << "CharacterTable::CharacterTable: null point group" << endl;
        // exit(1);
        throw PSIEXCEPTION("CharacterTable::CharacterTable: null point group");
    }

    if (make_table() < 0) {
        throw PSIEXCEPTION("CharacterTable::CharacterTable: could not make table");
    }
}

CharacterTable::CharacterTable(const std::string& cpg)
    : nt(0), pg(PointGroups::C1), nirrep_(0), gamma_(nullptr), symop(nullptr), _inv(nullptr), symb(cpg), bits_(0) {
    // Check the symbol coming in
    if (!PointGroup::full_name_to_bits(cpg, bits_)) {
        outfile->Printf("CharacterTable: Invalid point group name: %s\n", cpg.c_str());
        throw PSIEXCEPTION("CharacterTable: Invalid point group name provided.");
    }
    common_init();
}

CharacterTable::CharacterTable(unsigned char bits) : bits_(bits) {
    symb = PointGroup::bits_to_basic_name(bits);
    common_init();
}

unsigned char CharacterTable::bits() { return bits_; }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
