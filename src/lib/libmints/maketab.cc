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

    if (!g) return 0;

    gamma_ = new IrreducibleRepresentation[nirrep_];

    symop = new SymmetryOperation[g];
    SymmetryOperation so;

    _inv = new int[g];

    // this array forms a reducible representation for rotations about x,y,z
    double *rot = new double[g];
    memset(rot,0,sizeof(double)*g);

    // this array forms a reducible representation for translations along x,y,z
    double *trans = new double[g];
    memset(trans,0,sizeof(double)*g);

    // the angle to rotate about the principal axis
    double theta = (nt) ? 2.0*M_PI/nt : 2.0*M_PI;

    switch (pg) {

    case C1:
        // no symmetry case
        gamma_[0].init(1,1,"A");
        gamma_[0].nrot_ = 3;
        gamma_[0].ntrans_ = 3;
        gamma_[0].rep[0][0][0] = 1.0;

        symop[0].E();

        break;

    case CI:
        // equivalent to S2 about the z axis
        gamma_[0].init(2,1,"Ag");
        gamma_[0].rep[0][0][0] = 1.0;
        gamma_[0].rep[1][0][0] = 1.0;
        gamma_[0].nrot_=3;

        gamma_[1].init(2,1,"Au");
        gamma_[1].rep[0][0][0] =  1.0;
        gamma_[1].rep[1][0][0] = -1.0;
        gamma_[1].ntrans_=3;

        symop[0].E();
        symop[1].i();

        break;

    case CS: // reflection through the xy plane
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

        symop[0].E();
        symop[1].sigma_h();

        break;

    case CN:
        // clockwise rotation about z axis by theta*i radians
        //
        // for odd n, the irreps are A and E1...E(nir-1)
        // for even n, the irreps are A, B, and E1...E(nir-2)
        //
        gamma_[0].init(g,1,"A");
        for (gi=0; gi < g; gi++)
            gamma_[0].rep[gi][0][0] = 1.0;

        i=1;

        if (!(nt%2)) {
            gamma_[1].init(g,1,"B");
            for (gi=0; gi < g; gi++)
                gamma_[1].rep[gi][0][0] = (gi%2) ? -1.0 : 1.0;

            i++;
        }

        ei=1;
        for (; i < nirrep_; i++, ei++) {
            IrreducibleRepresentation& ir = gamma_[i];

            if (nt==3 || nt==4)
                sprintf(label,"E");
            else
                sprintf(label,"E%d",ei);

            ir.init(g,2,label);
            ir.complex_=1;

            // identity
            ir.rep[0].E();

            // Cn
            ir.rep[1].rotation(ei*theta);

            // the other n-1 Cn's
            for (j=2; j < g; j++)
                ir.rep[j] = ir.rep[j-1].operate(ir.rep[1]);
        }

        // identity
        symop[0].E();

        // Cn
        symop[1].rotation(theta);

        // the other n-2 Cn's
        for (i=2; i < nt; i++)
            symop[i] = symop[i-1].operate(symop[1]);

        for (i=0; i < nt ; i++)
            rot[i] = trans[i] = symop[i].trace();

        break;

    case CNV:
        // clockwise rotation about z axis by theta*i radians, then
        // reflect through the xz plane
        //
        // for odd n, the irreps are A1, A2, and E1...E(nir-2)
        // for even n, the irreps are A1, A2, B1, B2, and E1...E(nir-4)
        //

        gamma_[0].init(g,1,"A1");
        gamma_[1].init(g,1,"A2");

        for (gi=0; gi < nt; gi++) {
            // Cn's
            gamma_[0].rep[gi][0][0] = 1.0;
            gamma_[1].rep[gi][0][0] = 1.0;

            // sigma's
            gamma_[0].rep[gi+nt][0][0] =  1.0;
            gamma_[1].rep[gi+nt][0][0] = -1.0;
        }

        if (!(nt%2)) {
            gamma_[2].init(g,1,"B1");
            gamma_[3].init(g,1,"B2");

            for (gi=0; gi < nt ; gi++) {
                double ci = (gi%2) ? -1.0 : 1.0;

                // Cn's
                gamma_[2].rep[gi][0][0] = ci;
                gamma_[3].rep[gi][0][0] = ci;

                // sigma's
                gamma_[2].rep[gi+nt][0][0] =  ci;
                gamma_[3].rep[gi+nt][0][0] = -ci;
            }
        }

        ei=1;
        for (i = (nt%2) ? 2 : 4; i < nirrep_; i++, ei++) {
            IrreducibleRepresentation& ir = gamma_[i];

            char lab[4];
            if (nt==3 || nt==4)
                sprintf(lab,"E");
            else
                sprintf(lab,"E%d",ei);

            ir.init(g,2,lab);

            // identity
            ir.rep[0].E();

            // Cn
            ir.rep[1].rotation(ei*theta);

            // the other n-2 Cn's
            for (j=2; j < nt; j++)
                ir.rep[j] = ir.rep[j-1].operate(ir.rep[1]);

            // sigma xz
            ir.rep[nt].sigma_xz();

            SymRep sr(2);
            sr.rotation(ei*theta/2.0);

            // the other n-1 sigma's
            for (j=nt+1; j < g; j++)
                ir.rep[j] = ir.rep[j-1].transform(sr);
        }

        // identity
        symop[0].E();

        // Cn
        symop[1].rotation(theta);

        // the other n-2 Cn's
        for (i=2; i < nt; i++)
            symop[i] = symop[i-1].operate(symop[1]);

        // sigma xz
        symop[nt].sigma_xz();

        so.rotation(theta/2.0);

        // the other n-1 sigma's
        for (j=nt+1; j < g; j++)
            symop[j] = symop[j-1].transform(so);

        for (i=0; i < nt ; i++) {
            rot[i] = trans[i] = symop[i].trace();

            rot[i+nt] = -symop[i+nt].trace();
            trans[i+nt] = symop[i+nt].trace();
        }

        break;

    case CNH:
        // lockwise rotation about z axis by theta*i radians,
        // as well as rotation-reflection about same axis

        //
        // for odd n, the irreps are A', A", E1'...E(nir/2-1)', E1"...E(nir/2-1)''
        // for even n, the irreps are Ag, Bg, Au, Bu,
        //                            E1g...E(nir/2-1)g, E1u...E(nir/2-1)u
        //
        gamma_[0].init(g,1, (nt%2) ? "A'" : "Ag", (nt%2) ? "Ap" : 0);
        gamma_[nirrep_/2].init(g,1, (nt%2) ? "A\"" : "Au", (nt%2) ? "Ap" : 0);

        for (gi=0; gi < nt; gi++) {
            // Cn's
            gamma_[0].rep[gi][0][0] = 1.0;
            gamma_[nirrep_/2].rep[gi][0][0] = 1.0;

            // Sn's
            gamma_[0].rep[gi+nt][0][0] = 1.0;
            gamma_[nirrep_/2].rep[gi+nt][0][0] = -1.0;
        }

        if (!(nt%2)) {
            gamma_[1].init(g,1,"Bg");
            gamma_[1+nirrep_/2].init(g,1,"Bu");

            for (gi=0; gi < nt; gi++) {
                double ci = (gi%2) ? -1.0 : 1.0;

                // Cn's
                gamma_[1].rep[gi][0][0] = ci;
                gamma_[1+nirrep_/2].rep[gi][0][0] = ci;

                // Sn's
                gamma_[1].rep[gi+nt][0][0] =  ci;
                gamma_[1+nirrep_/2].rep[gi+nt][0][0] = -ci;
            }
        }

        ei=1;
        for (i = (nt%2) ? 1 : 2; i < nirrep_/2 ; i++, ei++) {
            IrreducibleRepresentation& ir1 = gamma_[i];
            IrreducibleRepresentation& ir2 = gamma_[i+nirrep_/2];

            if (nt==3 || nt==4)
                sprintf(label,(nt%2) ? "E'" : "Eg");
            else
                sprintf(label,"E%d%s", ei, (nt%2) ? "'" : "g");

            ir1.init(g,2,label);

            if (nt==3 || nt==4)
                sprintf(label,(nt%2) ? "E\"" : "Eu");
            else
                sprintf(label,"E%d%s", ei, (nt%2) ? "\"" : "u");

            ir2.init(g,2,label);

            ir1.complex_=1;
            ir2.complex_=1;

            // identity
            ir1.rep[0].E();
            ir2.rep[0].E();

            // Cn
            ir1.rep[1].rotation(ei*theta);
            ir2.rep[1].rotation(ei*theta);

            for (j=2; j < nt; j++) {
                ir1.rep[j] = ir1.rep[j-1].operate(ir1.rep[1]);
                ir2.rep[j] = ir2.rep[j-1].operate(ir2.rep[1]);
            }

            // Sn's
            SymRep sr(2);
            sr.i();

            for (j=nt; j < g; j++) {
                ir1.rep[j] = ir1.rep[j-nt];
                ir2.rep[j] = ir2.rep[j-nt].operate(sr);
            }
        }

        // identity
        symop[0].E();

        // Cn
        symop[1].rotation(theta);

        // the other n-2 Cn's
        for (i=2; i < nt; i++)
            symop[i] = symop[i-1].operate(symop[1]);

        // Sn's, for odd nt, operate on Cn's with sigma_h, for even nt,
        // operate Cn's with i
        if (nt%2)
            so.sigma_h();
        else
            so.i();

        for (i=0; i < nt ; i++) {
            symop[i+nt] = symop[i].operate(so);

            rot[i] = trans[i] = symop[i].trace();
            trans[i+nt] = symop[i+nt].trace();
            rot[i+nt] = -trans[i+nt];
        }

        break;

    case SN:
        // clockwise rotation-reflection by theta*i radians about z axis
        //
        // for odd n/2, the irreps are Ag, Au, E1g...E(nir/2-1)g,E1u...E(nir/2-1)u
        // for even n/2, the irreps are A, B, E1...E(nir-2)
        //
        if ((nt/2)%2) {
            gamma_[0].init(g, 1, "Ag");
            gamma_[nirrep_/2].init(g, 1, "Au");

            for (gi=0; gi < nt/2; gi++) {
                gamma_[0].rep[gi][0][0] = 1.0;
                gamma_[nirrep_/2].rep[gi][0][0] = 1.0;

                gamma_[0].rep[gi+nt/2][0][0] =  1.0;
                gamma_[nirrep_/2].rep[gi+nt/2][0][0] = -1.0;
            }

            ei=1;
            for (i=1; i < nirrep_/2 ; i++, ei++) {
                IrreducibleRepresentation& ir1 = gamma_[i];
                IrreducibleRepresentation& ir2 = gamma_[i+nirrep_/2];

                if (nt==6)
                    sprintf(label,"Eg");
                else
                    sprintf(label,"E%dg",ei);

                ir1.init(g,2,label);
                ir1.complex_=1;

                if (nt==6)
                    sprintf(label,"Eu");
                else
                    sprintf(label,"E%du", ei);

                ir2.init(g,2,label);
                ir2.complex_=1;

                // identity
                ir1.rep[0].E();
                ir2.rep[0].E();

                // C(n/2)
                ir1.rep[1].rotation(ei*theta*2.0);
                ir2.rep[1].rotation(ei*theta*2.0);

                for (j=2; j < nt/2; j++) {
                    ir1.rep[j] = ir1.rep[j-1].operate(ir1.rep[1]);
                    ir2.rep[j] = ir2.rep[j-1].operate(ir2.rep[1]);
                }

                SymRep sr(2);
                sr.i();

                // Sn
                for (j=nt/2; j < nt; j++) {
                    ir1.rep[j] = ir1.rep[j-nt/2];
                    ir2.rep[j] = ir2.rep[j-nt/2].operate(sr);
                }
            }

            // identity
            symop[0].E();

            // Cn
            symop[1].rotation(2.0*theta);

            for (i=2; i < nt/2 ; i++)
                symop[i] = symop[i-1].operate(symop[1]);

            so.i();

            // Sn
            for (i=nt/2; i < nt; i++)
                symop[i] = symop[i-nt/2].operate(so);

            for (i=0; i < nt/2 ; i++) {
                rot[i] = trans[i] = symop[i].trace();

                trans[i+nt/2] = symop[i+nt/2].trace();
                rot[i+nt/2] = -trans[i+nt/2];
            }

        } else {
            gamma_[0].init(g, 1, "A");
            gamma_[1].init(g, 1, "B");

            for (gi=0; gi < nt; gi++) {
                gamma_[0].rep[gi][0][0] = 1.0;
                gamma_[1].rep[gi][0][0] = (gi%2) ? -1.0 : 1.0;
            }

            ei=1;
            for (i=2; i < nirrep_; i++, ei++) {
                IrreducibleRepresentation& ir = gamma_[i];

                if (nt==4)
                    sprintf(label,"E");
                else
                    sprintf(label,"E%d",ei);

                ir.init(g,2,label);
                ir.complex_ = 1;

                // identity
                ir.rep[0].E();

                // Sn
                ir.rep[1].rotation(ei*theta);

                for (j=2; j < nt; j++)
                    ir.rep[j] = ir.rep[j-1].operate(ir.rep[1]);
            }

            // identity
            symop[0].E();

            // Sn
            symop[1].rotation(theta);
            symop[1][2][2] = -1.0;

            for (i=2; i < nt ; i++)
                symop[i] = symop[i-1].operate(symop[1]);

            for (i=0; i < nt ; i++) {
                trans[i] = symop[i].trace();
                rot[i] = (i%2) ? -trans[i] : trans[i];
            }
        }

        break;

    case DN:
        // clockwise rotation about z axis, followed by C2 about x axis

        // D2 is a special case
        if (nt==2) {
            gamma_[0].init(g,1,"A");
            gamma_[1].init(g,1,"B1");
            gamma_[2].init(g,1,"B2");
            gamma_[3].init(g,1,"B3");

            for (i=0; i < g; i++) {
                gamma_[0].rep[i][0][0] = 1.0;
                gamma_[1].rep[i][0][0] = (i < 2) ? 1.0 : -1.0;
                gamma_[2].rep[i][0][0] = (i % 2) ? -1.0 : 1.0;
                gamma_[3].rep[i][0][0] = (i < 2) ?
                            ((i % 2) ? -1.0 : 1.0) : ((i%2) ? 1.0 : -1.0);
            }
        } else {
            // Dn is isomorphic with Cnv
            //
            // for odd n, the irreps are A1, A2, and E1...E(nir-2)
            // for even n, the irreps are A1, A2, B1, B2, and E1...E(nir-4)
            //
            gamma_[0].init(g,1,"A1");
            gamma_[1].init(g,1,"A2");

            for (gi=0; gi < nt; gi++) {
                // Cn's
                gamma_[0].rep[gi][0][0] = 1.0;
                gamma_[1].rep[gi][0][0] = 1.0;

                // C2's
                gamma_[0].rep[gi+nt][0][0] =  1.0;
                gamma_[1].rep[gi+nt][0][0] = -1.0;
            }

            i=2;

            if (!(nt%2)) {
                gamma_[2].init(g,1,"B1");
                gamma_[3].init(g,1,"B2");

                for (gi=0; gi < nt ; gi++) {
                    double ci = (gi%2) ? -1.0 : 1.0;

                    // Cn's
                    gamma_[2].rep[gi][0][0] = ci;
                    gamma_[3].rep[gi][0][0] = ci;

                    // sigma's
                    gamma_[2].rep[gi+nt][0][0] =  ci;
                    gamma_[3].rep[gi+nt][0][0] = -ci;
                }

                i = 4;
            }

            ei=1;
            for (; i < nirrep_; i++, ei++) {
                IrreducibleRepresentation& ir = gamma_[i];

                char lab[4];
                if (nt==3 || nt==4)
                    sprintf(lab,"E");
                else
                    sprintf(lab,"E%d",ei);

                ir.init(g,2,lab);

                // identity
                ir.rep[0].E();

                // Cn
                ir.rep[1].rotation(ei*theta);

                for (j=2; j < nt; j++)
                    ir.rep[j] = ir.rep[j-1].operate(ir.rep[1]);

                // C2(x)
                ir.rep[nt].c2_y();

                SymRep sr(2);
                sr.rotation(ei*theta/2.0);

                for (j=nt+1; j < 2*nt; j++)
                    ir.rep[j] = ir.rep[j-1].transform(sr);
            }
        }

        // identity
        symop[0].E();

        // Cn
        symop[1].rotation(theta);

        for (i=2; i < nt; i++)
            symop[i] = symop[i-1].operate(symop[1]);

        // C2(x)
        symop[nt].c2_y();

        so.rotation(theta/2.0);

        for (i=nt+1; i < 2*nt; i++)
            symop[i] = symop[i-1].transform(so);

        for (i=0; i < 2*nt ; i++)
            rot[i] = trans[i] = symop[i].trace();

        break;

    case DND:
        // rotation reflection about z axis by theta/2 radians, followed
        // by c2 about x axis, then reflection through yz plane
        //
        // for odd n, the irreps are A1g, A2g, A1u, A2u, E1g...E(nir/2-2)g,
        //                                               E1u...E(nir/2-2)u
        // for even n, the irreps are A1, A2, B1, B2, E1...E(nir-4)
        //

        if (nt%2) {
            gamma_[0].init(g,1,"A1g");
            gamma_[1].init(g,1,"A2g");

            for (gi=0; gi < g; gi++) {
                gamma_[0].rep[gi][0][0] = 1.0;
                gamma_[1].rep[gi][0][0] = (gi/nt==0 || gi/nt==2) ? 1.0 : -1.0;
            }

            i=nirrep_/2;
            j=i+1;
            gamma_[i].init(g,1,"A1u");
            gamma_[j].init(g,1,"A2u");

            for (gi=0; gi < g/2; gi++) {
                gamma_[i].rep[gi][0][0] = gamma_[0].rep[gi][0][0];
                gamma_[j].rep[gi][0][0] = gamma_[1].rep[gi][0][0];

                gamma_[i].rep[gi+g/2][0][0] = -gamma_[0].rep[gi][0][0];
                gamma_[j].rep[gi+g/2][0][0] = -gamma_[1].rep[gi][0][0];
            }

            ei=1;

            for (i=2; i < nirrep_/2 ; i++, ei++) {
                IrreducibleRepresentation& ir1 = gamma_[i];
                IrreducibleRepresentation& ir2 = gamma_[i+nirrep_/2];

                if (nt==3) {
                    ir1.init(g,2,"Eg");
                    ir2.init(g,2,"Eu");
                } else {
                    sprintf(label,"E%dg",ei);
                    ir1.init(g,2,label);
                    sprintf(label,"E%du",ei);
                    ir2.init(g,2,label);
                }

                // identity
                ir1.rep[0].E();

                // Cn
                ir1.rep[1].rotation(ei*theta);

                for (j=2; j < nt; j++)
                    ir1.rep[j] = ir1.rep[j-1].operate(ir1.rep[1]);

                // C2(x)
                ir1.rep[nt].c2_y();

                for (j=nt+1; j < 2*nt; j++)
                    ir1.rep[j] = ir1.rep[j-1].transform(ir1.rep[1]);

                for (j=0; j < 2*nt; j++)
                    ir2.rep[j] = ir1.rep[j];

                // Sn and sigma d
                SymRep sr(2);
                sr.i();

                for (j=2*nt; j < g; j++) {
                    ir1.rep[j] = ir1.rep[j-2*nt];
                    ir2.rep[j] = ir2.rep[j-2*nt].operate(sr);
                }
            }

            // identity
            symop[0].E();

            // Cn
            symop[1].rotation(theta);

            for (i=2; i < nt; i++)
                symop[i] = symop[i-1].operate(symop[1]);

            // C2(x)
            symop[nt].c2_y();

            for (i=nt+1; i < 2*nt; i++)
                symop[i] = symop[i-1].transform(symop[1]);

            // i + n-1 S2n + n sigma
            so.i();
            for (i=2*nt; i < g; i++)
                symop[i] = symop[i-2*nt].operate(so);

            for (i=0; i < g; i++) {
                trans[i] = symop[i].trace();
                rot[i] = (i < g/2) ? trans[i] : -trans[i];
            }

        } else { // even nt

            gamma_[0].init(g,1,"A1");
            gamma_[1].init(g,1,"A2");
            gamma_[2].init(g,1,"B1");
            gamma_[3].init(g,1,"B2");

            for (gi=0; gi < 2*nt; gi++) {
                // Sn
                gamma_[0].rep[gi][0][0] = 1.0;
                gamma_[1].rep[gi][0][0] = 1.0;
                gamma_[2].rep[gi][0][0] = (gi%2) ? -1.0 : 1.0;
                gamma_[3].rep[gi][0][0] = (gi%2) ? -1.0 : 1.0;

                // n C2's and n sigma's
                gamma_[0].rep[gi+2*nt][0][0] =  1.0;
                gamma_[1].rep[gi+2*nt][0][0] = -1.0;
                gamma_[2].rep[gi+2*nt][0][0] = (gi%2) ? -1.0 : 1.0;
                gamma_[3].rep[gi+2*nt][0][0] = (gi%2) ? 1.0 : -1.0;
            }

            ei=1;
            for (i=4; i < nirrep_; i++, ei++) {
                IrreducibleRepresentation& ir = gamma_[i];

                if (nt==2)
                    sprintf(label,"E");
                else
                    sprintf(label,"E%d",ei);

                ir.init(g,2,label);

                // identity
                ir.rep[0].E();

                // S2n
                ir.rep[1].rotation(ei*theta/2.0);

                for (j=2; j < 2*nt; j++)
                    ir.rep[j] = ir.rep[j-1].operate(ir.rep[1]);

                // C2(x) + sigma_d
                ir.rep[2*nt].c2_y();

                for (j=2*nt+1; j < g; j++)
                    ir.rep[j] = ir.rep[j-1].operate(ir.rep[1]);
            }

            // identity
            symop[0].E();

            // Sn's
            symop[1].rotation(theta/2.0);
            symop[1][2][2] = -1.0;

            for (i=2; i < 2*nt; i++)
                symop[i] = symop[i-1].operate(symop[1]);

            // C2(x)
            symop[2*nt].c2_y();

            for (i=2*nt+1; i < g; i++)
                symop[i] = symop[i-1].operate(symop[1]);

            for (i=0; i < g; i++) {
                trans[i] = symop[i].trace();
                rot[i] = (i%2) ? -trans[i] : trans[i];
            }
        }

        break;

    case DNH:
        // clockwise rotation and rotation-reflection about z axis,
        // followed by c2 about x axis and then reflection
        // through xz

        i=nirrep_/2; j=i+1;

        if (nt%2) {
            gamma_[0].init(g,1,"A1'");
            gamma_[1].init(g,1,"A2'");
            gamma_[i].init(g,1,"A1\"");
            gamma_[j].init(g,1,"A2\"");
        } else {
            if (nt==2) {
                gamma_[0].init(g,1,"Ag");
                gamma_[1].init(g,1,"B1g");
                gamma_[4].init(g,1,"Au");
                gamma_[5].init(g,1,"B1u");
            } else {
                gamma_[0].init(g,1,"A1g");
                gamma_[1].init(g,1,"A2g");
                gamma_[i].init(g,1,"A1u");
                gamma_[j].init(g,1,"A2u");
            }
        }

        for (gi=0; gi < nt; gi++) {
            // E + n-1 Cn's
            gamma_[0].rep[gi][0][0] = gamma_[1].rep[gi][0][0] =
                    gamma_[i].rep[gi][0][0] = gamma_[j].rep[gi][0][0] = 1.0;

            // n C2's
            gamma_[0].rep[gi+nt][0][0] = gamma_[i].rep[gi+nt][0][0] =  1.0;
            gamma_[1].rep[gi+nt][0][0] = gamma_[j].rep[gi+nt][0][0] = -1.0;

            // i + n-1 S2n's
            gamma_[0].rep[gi+2*nt][0][0] = gamma_[1].rep[gi+2*nt][0][0] =  1.0;
            gamma_[i].rep[gi+2*nt][0][0] = gamma_[j].rep[gi+2*nt][0][0] = -1.0;

            // n sigma's
            gamma_[0].rep[gi+3*nt][0][0] = gamma_[j].rep[gi+3*nt][0][0] =  1.0;
            gamma_[i].rep[gi+3*nt][0][0] = gamma_[1].rep[gi+3*nt][0][0] = -1.0;
        }

        if (!(nt%2)) {
            if (nt==2) {
                gamma_[2].init(g,1,"B2g");
                gamma_[3].init(g,1,"B3g");
                gamma_[6].init(g,1,"B2u");
                gamma_[7].init(g,1,"B3u");
            } else {
                gamma_[2].init(g,1,"B1g");
                gamma_[3].init(g,1,"B2g");
                gamma_[i+2].init(g,1,"B1u");
                gamma_[j+2].init(g,1,"B2u");
            }

            for (gi=0; gi < nt; gi++) {
                // E + n-1 Cn's
                gamma_[2].rep[gi][0][0] = gamma_[3].rep[gi][0][0] =
                        gamma_[i+2].rep[gi][0][0] = gamma_[j+2].rep[gi][0][0] =
                        (gi%2) ? -1.0 : 1.0;

                // n C2's
                gamma_[2].rep[gi+nt][0][0] = gamma_[i+2].rep[gi+nt][0][0] =
                        (gi%2) ? -1.0 : 1.0;
                gamma_[3].rep[gi+nt][0][0] = gamma_[j+2].rep[gi+nt][0][0] =
                        (gi%2) ? 1.0 : -1.0;

                // i + n-1 S2n's
                gamma_[2].rep[gi+2*nt][0][0] = gamma_[3].rep[gi+2*nt][0][0] =
                        (gi%2) ? -1.0 : 1.0;
                gamma_[i+2].rep[gi+2*nt][0][0] = gamma_[j+2].rep[gi+2*nt][0][0] =
                        (gi%2) ? 1.0 : -1.0;

                // n sigma's
                gamma_[2].rep[gi+3*nt][0][0] = gamma_[j+2].rep[gi+3*nt][0][0] =
                        (gi%2) ? -1.0 : 1.0;
                gamma_[i+2].rep[gi+3*nt][0][0] = gamma_[3].rep[gi+3*nt][0][0] =
                        (gi%2) ? 1.0 : -1.0;
            }
        }

        ei=1;
        for (i = (nt%2) ? 2 : 4; i < nirrep_/2 ; i++, ei++) {
            IrreducibleRepresentation& ir1 = gamma_[i];
            IrreducibleRepresentation& ir2 = gamma_[i+nirrep_/2];

            if (nt==3) {
                ir1.init(g,2,"E'");
                ir2.init(g,2,"E\"");
            } else if (nt==4) {
                ir1.init(g,2,"Eg");
                ir2.init(g,2,"Eu");
            } else {
                sprintf(label,"E%d%s", ei, (nt%2) ? "'" : "g");
                ir1.init(g,2,label);

                sprintf(label,"E%d%s", ei, (nt%2) ? "\"" : "u");
                ir2.init(g,2,label);
            }

            // identity
            ir1.rep[0].E();

            // n-1 Cn's
            ir1.rep[1].rotation(ei*theta);

            for (j=2; j < nt; j++)
                ir1.rep[j] = ir1.rep[j-1].operate(ir1.rep[1]);

            // n C2's
            ir1.rep[nt].c2_y();

            SymRep sr(2);
            sr.rotation(ei*theta/2.0);

            for (j=nt+1; j < 2*nt; j++)
                ir1.rep[j] = ir1.rep[j-1].transform(sr);

            sr.i();
            for (j=0; j < 2*nt; j++) {
                ir1.rep[j+2*nt] = ir1.rep[j];
                ir2.rep[j] = ir1.rep[j];
                ir2.rep[j+2*nt] = ir1.rep[j].operate(sr);
            }
        }

        // identity
        symop[0].E();

        // n-1 Cn's
        symop[1].rotation(theta);

        for (i=2; i < nt; i++)
            symop[i] = symop[i-1].operate(symop[1]);

        // n C2's
        symop[nt].c2_y();

        so.rotation(theta/2.0);
        for (i=nt+1; i < 2*nt; i++)
            symop[i] = symop[i-1].transform(so);

        if (nt%2)
            so.sigma_h();
        else
            so.i();

        for (i=2*nt; i < g; i++)
            symop[i] = symop[i-2*nt].operate(so);

        for (i=0,j=2*nt; i < 2*nt ; i++,j++) {
            rot[i] = trans[i] = symop[i].trace();
            trans[j] = symop[j].trace();
            rot[j] = -trans[j];
        }

        break;

    case T:
        t();
        break;

    case TH:
        th();
        break;

    case TD:
        td();
        break;

    case O:
        o();
        break;

    case OH:
        oh();
        break;

    case I:
        this->i();
        break;

    case IH:
        ih();
        break;

    default:
        return -1;

    }

    /* ok, we have the reducible representation of the rotations and
   * translations, now let's use projection operators to find out how many
   * rotations and translations there are in each irrep
   */

    if (pg != C1 && pg != CI && pg != CS && pg != T && pg != TD && pg != TH &&
            pg != O && pg != OH && pg != I && pg != IH) {

        for (i=0; i < nirrep_; i++) {
            double nr=0; double nt=0;

            for (j=0; j < gamma_[i].g; j++) {
                nr += rot[j]*gamma_[i].character(j);
                nt += trans[j]*gamma_[i].character(j);
            }

            gamma_[i].nrot_ = (int) ((nr+0.5)/gamma_[i].g);
            gamma_[i].ntrans_ = (int) ((nt+0.5)/gamma_[i].g);
        }
    }

    delete[] rot;
    delete[] trans;

    // now find the inverse of each symop
    for (gi=0; gi < g; gi++) {
        int gj;
        for (gj=0; gj < g; gj++) {
            so = symop[gi].operate(symop[gj]);

            // is so a unit matrix?
            if (fabs(1.0-so[0][0]) < 1.0e-8 &&
                    fabs(1.0-so[1][1]) < 1.0e-8 &&
                    fabs(1.0-so[2][2]) < 1.0e-8) break;
        }

        if (gj==g) {
            // ExEnv::err0() << indent
            //      << "make_table: uh oh, can't find inverse of " << gi << endl;
            // abort();
            throw PSIEXCEPTION("make_table: uh oh, can't find inverse");
        }

        _inv[gi] = gj;
    }

    return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
