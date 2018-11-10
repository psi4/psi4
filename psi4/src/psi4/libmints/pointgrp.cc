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
// pointgrp.cc
//
// Additional modifications made by Justin Turney <jturney@ccqc.uga.edu>
// for use in PSI4.
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

/* pointgrp.cc -- implementation of the point group classes
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
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"

#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>

namespace psi {
namespace PointGroups {
void similar(unsigned char bits, unsigned char *sim, char &cnt) {
    static unsigned char cs[3] = {CsX, CsY, CsZ};
    static unsigned char c2v[3] = {C2vZ, C2vY, C2vX};
    static unsigned char c2h[3] = {C2hZ, C2hY, C2hX};
    static unsigned char c2[3] = {C2Z, C2Y, C2X};
    static unsigned char d2h[3] = {D2h, 0, 0};
    static unsigned char d2[3] = {D2, 0, 0};
    static unsigned char ci[3] = {Ci, 0, 0};
    static unsigned char c1[3] = {C1, 0, 0};

    switch (bits) {
        case CsX:
        case CsY:
        case CsZ:
            memcpy(sim, cs, sizeof(char) * 3);
            cnt = 3;
            break;
        case C2vZ:
        case C2vY:
        case C2vX:
            memcpy(sim, c2v, sizeof(char) * 3);
            cnt = 3;
            break;
        case C2hZ:
        case C2hY:
        case C2hX:
            memcpy(sim, c2h, sizeof(char) * 3);
            cnt = 3;
            break;
        case C2Z:
        case C2Y:
        case C2X:
            memcpy(sim, c2, sizeof(char) * 3);
            cnt = 3;
            break;
        case D2h:
            memcpy(sim, d2h, sizeof(char) * 1);
            cnt = 1;
            break;
        case Ci:
            memcpy(sim, ci, sizeof(char) * 1);
            cnt = 1;
            break;
        case C1:
            memcpy(sim, c1, sizeof(char) * 1);
            cnt = 1;
            break;
        case D2:
            memcpy(sim, d2, sizeof(char) * 1);
            cnt = 1;
            break;
        default:
            throw PSIEXCEPTION("Should not have reaced here.");
    }
}
}  // namespace PointGroups

////////////////////////////////////////////////////////////////////////

PointGroup::PointGroup() {
    set_symbol("c1");
    origin_[0] = origin_[1] = origin_[2] = 0;
}

PointGroup::PointGroup(const std::string &s) {
    if (full_name_to_bits(s, bits_) == false) throw PSIEXCEPTION("PointGroup: Unknown point group name provided.");
    set_symbol(bits_to_basic_name(bits_));
    origin_[0] = origin_[1] = origin_[2] = 0;
}

PointGroup::PointGroup(const std::string &s, const Vector3 &origin) {
    if (full_name_to_bits(s, bits_) == false) throw PSIEXCEPTION("PointGroup: Unknown point group name provided.");
    set_symbol(bits_to_basic_name(bits_));
    origin_ = origin;
}

PointGroup::PointGroup(unsigned char bits) : bits_(bits) {
    set_symbol(bits_to_basic_name(bits));
    origin_[0] = origin_[1] = origin_[2] = 0;
}

PointGroup::PointGroup(unsigned char bits, const Vector3 &origin) : bits_(bits) {
    set_symbol(bits_to_basic_name(bits));
    origin_ = origin;
}

PointGroup::PointGroup(const PointGroup &pg) { *this = pg; }

PointGroup::PointGroup(const std::shared_ptr<PointGroup> &pg) { *this = *pg.get(); }

PointGroup::~PointGroup() {}

PointGroup &PointGroup::operator=(const PointGroup &pg) {
    set_symbol(pg.symb);
    origin_ = pg.origin_;
    return *this;
}

void PointGroup::set_symbol(const std::string &sym) {
    if (sym.length()) {
        symb = sym;
    } else {
        set_symbol("c1");
    }
}

CharacterTable PointGroup::char_table() const {
    CharacterTable ret(bits_);
    return ret;
}

int PointGroup::equiv(const std::shared_ptr<PointGroup> &grp, double /*tol*/) const {
    if (symb != grp->symb) return 0;

    return 1;
}

bool PointGroup::full_name_to_bits(const std::string &pg, unsigned char &bits) {
    bool retvalue = true;

    if (iequals(pg, std::string("C1")))
        bits = PointGroups::C1;
    else if (iequals(pg, std::string("Ci")))
        bits = PointGroups::Ci;
    else if (iequals(pg, std::string("C2(x)")) || iequals(pg, std::string("C2x")) || iequals(pg, std::string("C2_x")))
        bits = PointGroups::C2X;
    else if (iequals(pg, std::string("C2(y)")) || iequals(pg, std::string("C2y")) || iequals(pg, std::string("C2_y")))
        bits = PointGroups::C2Y;
    else if (iequals(pg, std::string("C2(z)")) || iequals(pg, std::string("C2z")) || iequals(pg, std::string("C2_z")))
        bits = PointGroups::C2Z;
    else if (iequals(pg, std::string("Cs(x)")) || iequals(pg, std::string("Csx")) || iequals(pg, std::string("Cs_x")))
        bits = PointGroups::CsX;
    else if (iequals(pg, std::string("Cs(y)")) || iequals(pg, std::string("Csy")) || iequals(pg, std::string("Cs_y")))
        bits = PointGroups::CsY;
    else if (iequals(pg, std::string("Cs(z)")) || iequals(pg, std::string("Csz")) || iequals(pg, std::string("Cs_z")))
        bits = PointGroups::CsZ;
    else if (iequals(pg, std::string("D2")))
        bits = PointGroups::D2;
    else if (iequals(pg, std::string("C2v(X)")) || iequals(pg, std::string("C2vx")) ||
             iequals(pg, std::string("C2v_x")))
        bits = PointGroups::C2vX;
    else if (iequals(pg, std::string("C2v(Y)")) || iequals(pg, std::string("C2vy")) ||
             iequals(pg, std::string("C2v_y")))
        bits = PointGroups::C2vY;
    else if (iequals(pg, std::string("C2v(Z)")) || iequals(pg, std::string("C2vz")) ||
             iequals(pg, std::string("C2v_z")))
        bits = PointGroups::C2vZ;
    else if (iequals(pg, std::string("C2h(X)")) || iequals(pg, std::string("C2hx")) ||
             iequals(pg, std::string("C2h_x")))
        bits = PointGroups::C2hX;
    else if (iequals(pg, std::string("C2h(Y)")) || iequals(pg, std::string("C2hy")) ||
             iequals(pg, std::string("C2h_y")))
        bits = PointGroups::C2hY;
    else if (iequals(pg, std::string("C2h(Z)")) || iequals(pg, std::string("C2hz")) ||
             iequals(pg, std::string("C2h_z")))
        bits = PointGroups::C2hZ;
    else if (iequals(pg, std::string("D2h")))
        bits = PointGroups::D2h;

    // Ok, the user gave us Cs, C2v, C2h, C2, but no directionality
    else if (iequals(pg, std::string("Cs")))
        bits = PointGroups::CsX;
    else if (iequals(pg, std::string("C2v")))
        bits = PointGroups::C2vZ;
    else if (iequals(pg, std::string("C2h")))
        bits = PointGroups::C2hZ;
    else if (iequals(pg, std::string("C2")))
        bits = PointGroups::C2Z;

    else
        retvalue = false;

    return retvalue;
}

const char *PointGroup::bits_to_full_name(unsigned char bits) {
    switch (bits) {
        case PointGroups::C1:
            return "C1";
        case PointGroups::Ci:
            return "Ci";
        case PointGroups::C2X:
            return "C2(x)";
        case PointGroups::C2Y:
            return "C2(y)";
        case PointGroups::C2Z:
            return "C2(z)";
        case PointGroups::CsZ:
            return "Cs(Z)";
        case PointGroups::CsY:
            return "Cs(Y)";
        case PointGroups::CsX:
            return "Cs(X)";
        case PointGroups::D2:
            return "D2";
        case PointGroups::C2vX:
            return "C2v(X)";
        case PointGroups::C2vY:
            return "C2v(Y)";
        case PointGroups::C2vZ:
            return "C2v(Z)";
        case PointGroups::C2hX:
            return "C2h(X)";
        case PointGroups::C2hY:
            return "C2h(Y)";
        case PointGroups::C2hZ:
            return "C2h(Z)";
        case PointGroups::D2h:
            return "D2h";
        default:
            outfile->Printf("Unrecognized point group bits: %d\n", bits);
            throw PSIEXCEPTION("Unrecognized point group bits");
    }
}

const char *PointGroup::bits_to_basic_name(unsigned char bits) {
    switch (bits) {
        case PointGroups::C1:
            return "c1";
        case PointGroups::Ci:
            return "ci";
        case PointGroups::C2X:
        case PointGroups::C2Y:
        case PointGroups::C2Z:
            return "c2";
        case PointGroups::CsZ:
        case PointGroups::CsY:
        case PointGroups::CsX:
            return "cs";
        case PointGroups::D2:
            return "d2";
        case PointGroups::C2vX:
        case PointGroups::C2vY:
        case PointGroups::C2vZ:
            return "c2v";
        case PointGroups::C2hX:
        case PointGroups::C2hY:
        case PointGroups::C2hZ:
            return "c2h";
        case PointGroups::D2h:
            return "d2h";
        default:
            outfile->Printf("Unrecognized point group bits: %d\n", bits);
            throw PSIEXCEPTION("Unrecognized point group bits");
    }
}

void PointGroup::print(std::string out) const {
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    printer->Printf("PointGroup: %s\n", symb.c_str());
}

std::string PointGroup::irrep_bits_to_string(int irrep_bits) const {
    std::string irrep_str;
    const CharacterTable c_table = char_table();
    for (int irrep = 0; irrep < c_table.nirrep(); ++irrep) {
        if ((1 << irrep) & irrep_bits) {
            if (!irrep_str.empty()) {
                irrep_str += ", ";
            }
            irrep_str += c_table.gamma(irrep).symbol();
        }
    }
    return irrep_str;
}
}  // namespace psi

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
