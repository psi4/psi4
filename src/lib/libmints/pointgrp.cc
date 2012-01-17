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
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
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

#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <cmath>

#include <libmints/pointgrp.h>

using namespace std;
using namespace psi;

////////////////////////////////////////////////////////////////////////

PointGroup::PointGroup()
{
    set_symbol("c1");
    frame(0,0) = frame(1,1) = frame(2,2) = 1;
    origin_[0] = origin_[1] = origin_[2] =0;
}

PointGroup::PointGroup(const std::string& s)
{
    set_symbol(s);
    frame(0,0) = frame(1,1) = frame(2,2) = 1;
    origin_[0] = origin_[1] = origin_[2] =0;
}

PointGroup::PointGroup(const std::string& s, SymmetryOperation& so)
{
    set_symbol(s);
    frame = so;
    origin_[0] = origin_[1] = origin_[2] =0;
}

PointGroup::PointGroup(const std::string& s, SymmetryOperation& so,
                       const Vector3& origin)
{
    fprintf(outfile, "in PointGroup: %s\n", s.c_str());
    so.print(outfile);
    fprintf(outfile, "origin: %lf %lf %lf\n", origin[0], origin[1], origin[2]);

    set_symbol(s);
    frame = so;
    origin_ = origin;
}

PointGroup::PointGroup(const PointGroup& pg)
{
    *this = pg;
}

PointGroup::PointGroup(const boost::shared_ptr<PointGroup>& pg)
{
    *this = *pg.get();
}

PointGroup::~PointGroup()
{
}

PointGroup&
PointGroup::operator=(const PointGroup& pg)
{
    set_symbol(pg.symb);
    frame = pg.frame;
    origin_ = pg.origin_;
    return *this;
}

void
PointGroup::set_symbol(const std::string& sym)
{
    if (sym.length()) {
        symb = sym;
    } else {
        set_symbol("c1");
    }
}

CharacterTable
PointGroup::char_table() const
{
    CharacterTable ret(symb,frame);
    return ret;
}

int
PointGroup::equiv(const boost::shared_ptr<PointGroup> &grp, double tol) const
{
    if (symb != grp->symb)
        return 0;

    for (int i=0; i < 3; i++) {
        for (int j=0; j < 3; j++) {
            if (fabs(frame(i,j) - grp->frame(i,j)) > tol) return 0;
        }
    }

    return 1;
}

const char* PointGroup::bits_to_full_name(unsigned char bits)
{
    switch(bits) {
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
        throw PSIEXCEPTION("Unrecognized point group bits");
    }
}

const char* PointGroup::bits_to_basic_name(unsigned char bits)
{
    switch(bits) {
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
        throw PSIEXCEPTION("Unrecognized point group bits");
    }
}

void
PointGroup::print(FILE *out) const
{
    fprintf(outfile, "PointGroup: %s\n", symb.c_str());
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
