//
// maketab.cc
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

/* maketab.cc
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

#include <cmath>
#include <cstdio>
#include <cstring>

#include <libmints/pointgrp.h>

using namespace std;
using namespace psi;

/*
 * This function will generate a character table for the point group.
 * This character table is in the order that symmetry operations are
 * generated, not in Cotton order. If this is a problem, tough.
 * Also generate the transformation matrices.
 */

int CharacterTable::make_table()
{
    int i,j,ei,gi;
    char label[4];

    switch (bits_) {
    case PointGroups::C1:
        nirrep_ = 1;
        nt = 1;
        break;
    case PointGroups::CsX:
    case PointGroups::CsY:
    case PointGroups::CsZ:
    case PointGroups::Ci:
        nirrep_ = 2;
        nt = 1;
        break;
    case PointGroups::C2X:
    case PointGroups::C2Y:
    case PointGroups::C2Z:
        nirrep_ = 2;
        nt = 2;
        break;
    case PointGroups::C2hX:
    case PointGroups::C2hY:
    case PointGroups::C2hZ:
    case PointGroups::C2vX:
    case PointGroups::C2vY:
    case PointGroups::C2vZ:
    case PointGroups::D2:
        nirrep_ = 4;
        nt = 2;
        break;
    case PointGroups::D2h:
        nirrep_ = 8;
        nt = 2;
        break;
    default:
        throw PSIEXCEPTION("Should not have receached here!");
    }

    if (!nirrep_) return 0;

    gamma_ = new IrreducibleRepresentation[nirrep_];
    symop = new SymmetryOperation[nirrep_];
    SymmetryOperation so;

    _inv = new int[nirrep_];

    // this array forms a reducible representation for rotations about x,y,z
    double *rot = new double[nirrep_];
    memset(rot,0,sizeof(double)*nirrep_);

    // this array forms a reducible representation for translations along x,y,z
    double *trans = new double[nirrep_];
    memset(trans,0,sizeof(double)*nirrep_);

    // the angle to rotate about the principal axis
    double theta = (nt) ? 2.0*M_PI/nt : 2.0*M_PI;

    // Handle irreducible representations:
    switch (bits_) {
    case PointGroups::C1:
        // no symmetry case
        gamma_[0].init(1,1,"A");
        gamma_[0].nrot_ = 3;
        gamma_[0].ntrans_ = 3;
        gamma_[0].rep[0][0][0] = 1.0;

        break;

    case PointGroups::CsX: // reflection through the yz plane
    case PointGroups::CsY: // reflection through the xz plane
    case PointGroups::CsZ: // reflection through the xy plane
        gamma_[0].init(2,1,"A'","Ap");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].nrot_=1;
        gamma_[0].ntrans_=2;

        gamma_[1].init(2,1,"A\"","App");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] = -1.0;
        gamma_[1].nrot_=2;
        gamma_[1].ntrans_=1;

        break;

    case PointGroups::Ci:
        // equivalent to S2 about the z axis
        gamma_[0].init(2,1,"Ag");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].nrot_=3;

        gamma_[1].init(2,1,"Au");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] = -1.0;
        gamma_[1].ntrans_=3;

        break;

    case PointGroups::C2X:
    case PointGroups::C2Y:
    case PointGroups::C2Z:
        gamma_[0].init(2,1,"A","A");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].nrot_=1;
        gamma_[0].ntrans_=1;

        gamma_[1].init(2,1,"B","B");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] = -1.0;
        gamma_[1].nrot_=2;
        gamma_[1].ntrans_=2;

        break;

    case PointGroups::C2hX:
    case PointGroups::C2hY:
    case PointGroups::C2hZ:
        gamma_[0].init(4,1,"Ag","Ag");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].rep[2][0][0] = 1.0;
        gamma_[0].rep[3][0][0] = 1.0;
        gamma_[0].nrot_=1;
        gamma_[0].ntrans_=0;

        gamma_[1].init(4,1,"Bg","Bg");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] = -1.0;
        gamma_[1].rep[2][0][0] =  1.0;
        gamma_[1].rep[3][0][0] = -1.0;
        gamma_[1].nrot_=2;
        gamma_[1].ntrans_=0;

        gamma_[2].init(4, 1,"Au","Au");
        gamma_[2].rep[0][0][0] =  1.0;
        gamma_[2].rep[1][0][0] =  1.0;
        gamma_[2].rep[2][0][0] = -1.0;
        gamma_[2].rep[3][0][0] = -1.0;
        gamma_[2].nrot_=0;
        gamma_[2].ntrans_=1;

        gamma_[3].init(4, 1,"Bu","Bu");
        gamma_[3].rep[0][0][0] =  1.0;
        gamma_[3].rep[1][0][0] = -1.0;
        gamma_[3].rep[2][0][0] = -1.0;
        gamma_[3].rep[3][0][0] =  1.0;
        gamma_[3].nrot_=0;
        gamma_[3].ntrans_=2;

        break;

    case PointGroups::C2vX:
    case PointGroups::C2vY:
    case PointGroups::C2vZ:
        gamma_[0].init(4,1,"A1","A1");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].rep[2][0][0] = 1.0;
        gamma_[0].rep[3][0][0] = 1.0;
        gamma_[0].nrot_=0;
        gamma_[0].ntrans_=1;

        gamma_[1].init(4,1,"A2","A2");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] =  1.0;
        gamma_[1].rep[2][0][0] = -1.0;
        gamma_[1].rep[3][0][0] = -1.0;
        gamma_[1].nrot_=1;
        gamma_[1].ntrans_=0;

        gamma_[2].init(4, 1,"B1","B1");
        gamma_[2].rep[0][0][0] =  1.0;
        gamma_[2].rep[1][0][0] = -1.0;
        gamma_[2].rep[2][0][0] =  1.0;
        gamma_[2].rep[3][0][0] = -1.0;
        gamma_[2].nrot_=1;
        gamma_[2].ntrans_=1;

        gamma_[3].init(4, 1,"B2","B2");
        gamma_[3].rep[0][0][0] =  1.0;
        gamma_[3].rep[1][0][0] = -1.0;
        gamma_[3].rep[2][0][0] = -1.0;
        gamma_[3].rep[3][0][0] =  1.0;
        gamma_[3].nrot_=1;
        gamma_[3].ntrans_=1;

        break;

    case PointGroups::D2:
        gamma_[0].init(4,1,"A","A");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].rep[2][0][0] = 1.0;
        gamma_[0].rep[3][0][0] = 1.0;
        gamma_[0].nrot_=0;
        gamma_[0].ntrans_=0;

        gamma_[1].init(4,1,"B1","B1");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] =  1.0;
        gamma_[1].rep[2][0][0] = -1.0;
        gamma_[1].rep[3][0][0] = -1.0;
        gamma_[1].nrot_=1;
        gamma_[1].ntrans_=1;

        gamma_[2].init(4, 1,"B2","B2");
        gamma_[2].rep[0][0][0] =  1.0;
        gamma_[2].rep[1][0][0] = -1.0;
        gamma_[2].rep[2][0][0] =  1.0;
        gamma_[2].rep[3][0][0] = -1.0;
        gamma_[2].nrot_=1;
        gamma_[2].ntrans_=1;

        gamma_[3].init(4, 1,"B3","B3");
        gamma_[3].rep[0][0][0] =  1.0;
        gamma_[3].rep[1][0][0] = -1.0;
        gamma_[3].rep[2][0][0] = -1.0;
        gamma_[3].rep[3][0][0] =  1.0;
        gamma_[3].nrot_=1;
        gamma_[3].ntrans_=1;

        break;

    case PointGroups::D2h:

        gamma_[0].init(8,1,"Ag","Ag");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].rep[2][0][0] = 1.0;
        gamma_[0].rep[3][0][0] = 1.0;
        gamma_[0].rep[4][0][0] = 1.0;
        gamma_[0].rep[5][0][0] = 1.0;
        gamma_[0].rep[6][0][0] = 1.0;
        gamma_[0].rep[7][0][0] = 1.0;
        gamma_[0].nrot_=0;
        gamma_[0].ntrans_=0;

        gamma_[1].init(8,1,"B1g","B1g");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] =  1.0;
        gamma_[1].rep[2][0][0] = -1.0;
        gamma_[1].rep[3][0][0] = -1.0;
        gamma_[1].rep[4][0][0] =  1.0;
        gamma_[1].rep[5][0][0] =  1.0;
        gamma_[1].rep[6][0][0] = -1.0;
        gamma_[1].rep[7][0][0] = -1.0;
        gamma_[1].nrot_=1;
        gamma_[1].ntrans_=0;

        gamma_[2].init(8,1,"B2g","B2g");
        gamma_[2].rep[0][0][0] =  1.0;
        gamma_[2].rep[1][0][0] = -1.0;
        gamma_[2].rep[2][0][0] =  1.0;
        gamma_[2].rep[3][0][0] = -1.0;
        gamma_[2].rep[4][0][0] =  1.0;
        gamma_[2].rep[5][0][0] = -1.0;
        gamma_[2].rep[6][0][0] =  1.0;
        gamma_[2].rep[7][0][0] = -1.0;
        gamma_[2].nrot_=1;
        gamma_[2].ntrans_=0;

        gamma_[3].init(8,1,"B3g","B3g");
        gamma_[3].rep[0][0][0] =  1.0;
        gamma_[3].rep[1][0][0] = -1.0;
        gamma_[3].rep[2][0][0] = -1.0;
        gamma_[3].rep[3][0][0] =  1.0;
        gamma_[3].rep[4][0][0] =  1.0;
        gamma_[3].rep[5][0][0] = -1.0;
        gamma_[3].rep[6][0][0] = -1.0;
        gamma_[3].rep[7][0][0] =  1.0;
        gamma_[3].nrot_=1;
        gamma_[3].ntrans_=0;

        gamma_[4].init(8,1,"Au","Au");
        gamma_[4].rep[0][0][0] =  1.0;
        gamma_[4].rep[1][0][0] =  1.0;
        gamma_[4].rep[2][0][0] =  1.0;
        gamma_[4].rep[3][0][0] =  1.0;
        gamma_[4].rep[4][0][0] = -1.0;
        gamma_[4].rep[5][0][0] = -1.0;
        gamma_[4].rep[6][0][0] = -1.0;
        gamma_[4].rep[7][0][0] = -1.0;
        gamma_[4].nrot_=0;
        gamma_[4].ntrans_=0;

        gamma_[5].init(8,1,"B1u","B1u");
        gamma_[5].rep[0][0][0] =  1.0;
        gamma_[5].rep[1][0][0] =  1.0;
        gamma_[5].rep[2][0][0] = -1.0;
        gamma_[5].rep[3][0][0] = -1.0;
        gamma_[5].rep[4][0][0] = -1.0;
        gamma_[5].rep[5][0][0] = -1.0;
        gamma_[5].rep[6][0][0] =  1.0;
        gamma_[5].rep[7][0][0] =  1.0;
        gamma_[5].nrot_=0;
        gamma_[5].ntrans_=1;

        gamma_[6].init(8,1,"B2u","B2u");
        gamma_[6].rep[0][0][0] =  1.0;
        gamma_[6].rep[1][0][0] = -1.0;
        gamma_[6].rep[2][0][0] =  1.0;
        gamma_[6].rep[3][0][0] = -1.0;
        gamma_[6].rep[4][0][0] = -1.0;
        gamma_[6].rep[5][0][0] =  1.0;
        gamma_[6].rep[6][0][0] = -1.0;
        gamma_[6].rep[7][0][0] =  1.0;
        gamma_[6].nrot_=0;
        gamma_[6].ntrans_=1;

        gamma_[7].init(8,1,"B3u","B3u");
        gamma_[7].rep[0][0][0] =  1.0;
        gamma_[7].rep[1][0][0] = -1.0;
        gamma_[7].rep[2][0][0] = -1.0;
        gamma_[7].rep[3][0][0] =  1.0;
        gamma_[7].rep[4][0][0] = -1.0;
        gamma_[7].rep[5][0][0] =  1.0;
        gamma_[7].rep[6][0][0] =  1.0;
        gamma_[7].rep[7][0][0] = -1.0;
        gamma_[7].nrot_=0;
        gamma_[7].ntrans_=1;

        break;
    }

    // Handle symmetry operations
    symop[0].E();

    switch (bits_) {
    case PointGroups::C1:

        // nothing to do.

        break;

    case PointGroups::Ci:

        symop[1].i();

        break;

    case PointGroups::CsX: // reflection through the yz plane

        symop[1].sigma_yz();

        break;

    case PointGroups::CsY: // reflection through the xz plane

        symop[1].sigma_xz();

        break;

    case PointGroups::CsZ: // reflection through the xy plane

        symop[1].sigma_xy();

        break;

    case PointGroups::C2X:

        symop[1].c2_x();

        break;

    case PointGroups::C2Y:

        symop[1].c2_y();

        break;

    case PointGroups::C2Z:

        symop[1].rotation(2);

        break;

    case PointGroups::C2hX:

        symop[1].c2_x();
        symop[2].i();
        symop[3].sigma_yz();

        break;

    case PointGroups::C2hY:

        symop[1].c2_y();
        symop[2].i();
        symop[3].sigma_xz();

        break;

    case PointGroups::C2hZ:

        symop[1].rotation(2);
        symop[2].i();
        symop[3].sigma_xy();

        break;

    case PointGroups::C2vX:

        symop[1].c2_x();
        symop[2].sigma_xy();
        symop[3].sigma_xz();

        break;

    case PointGroups::C2vY:

        symop[1].c2_y();
        symop[2].sigma_xy();
        symop[3].sigma_yz();

        break;

    case PointGroups::C2vZ:

        symop[1].rotation(2);
        symop[2].sigma_xz();
        symop[3].sigma_yz();

        break;

    case PointGroups::D2:

        symop[1].rotation(2);
        symop[2].c2_y();
        symop[3].c2_x();

        break;

    case PointGroups::D2h:

        symop[1].rotation(2);
        symop[2].c2_y();
        symop[3].c2_x();
        symop[4].i();
        symop[5].sigma_xy();
        symop[6].sigma_xz();
        symop[7].sigma_yz();

        break;

    default:
        return -1;

    }

    // now find the inverse of each symop
    for (gi=0; gi < nirrep_; gi++) {
        int gj;
        for (gj=0; gj < nirrep_; gj++) {
            so = symop[gi].operate(symop[gj]);

            // is so a unit matrix?
            if (fabs(1.0-so[0][0]) < 1.0e-8 &&
                    fabs(1.0-so[1][1]) < 1.0e-8 &&
                    fabs(1.0-so[2][2]) < 1.0e-8) break;
        }

        if (gj==nirrep_) {
            // ExEnv::err0() << indent
            //      << "make_table: uh oh, can't find inverse of " << gi << endl;
            // abort();
            throw PSIEXCEPTION("make_table: uh oh, can't find inverse");
        }

        _inv[gi] = gj;
    }

    // Check the bits of the operator make sure they make what
    // we were given.
    unsigned char sym_bits = 0;
    for (int i=0; i<nirrep_; ++i) {
        sym_bits |= symop[i].bit();
    }

    if (sym_bits != bits_)
        throw PSIEXCEPTION("make_table: Symmetry operators did not match the point group given.");

    return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
