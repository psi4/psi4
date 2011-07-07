/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id: pointgroup.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/
#ifndef MAD_POINTGROUP_H
#define MAD_POINTGROUP_H

/// \file pointgroup.h
/// \brief Implements basic functionality for Abelian point groups (D2h and subgroups)

#include <tensor/tensor.h>
#include <world/array.h>
#include <string>
#include <ostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

class PointGroup;
std::ostream& operator<<(std::ostream& s, const PointGroup& g);

class PointGroup {
    std::string name;           //< group name
    int order;                  //< group order
    std::string irs[8];         //< Names of the irreps
    std::string ops[8];         //< Names of the operators
    int c[8][8];                //< Character table

public:
    typedef madness::Vector<double,3> coordT;

    /// Constructs point group by name (D2h and subgroups only)
    PointGroup(const std::string name) {
        this->name = name;
        if (name == "C1") {
            order = 1;
            irs[0] = "a";
            ops[0] = "e";
            c[0][0] = 1;
        }
        else if (name == "C2") {
            order = 2;
            irs[0] = "a"; irs[1] = "b";
            ops[0] = "e"; ops[1] = "c2z";
            c[0][0] = 1;  c[0][1] = 1;
            c[1][0] = 1;  c[1][1] =-1;
        }
        else if (name == "Ci") {
            order = 2;
            irs[0] = "ag"; irs[1] = "au";
            ops[0] = "e"; ops[1] = "i";
            c[0][0] = 1;  c[0][1] = 1;
            c[1][0] = 1;  c[1][1] =-1;
        }
        else if (name == "Cs") {
            order = 2;
            irs[0] = "a"; irs[1] = "a'";
            ops[0] = "e"; ops[1] = "sxy";
            c[0][0] = 1;  c[0][1] = 1;
            c[1][0] = 1;  c[1][1] =-1;
        }
        else if (name == "C2h") {
            order = 4;
            irs[0] = "ag"; irs[1] = "au"; irs[2] = "bg"; irs[3] = "bu";
            ops[0] = "e";  ops[1] = "c2z";ops[2] = "sxy";ops[3] = "i";
            c[0][0] = 1; c[0][1] = 1; c[0][2] = 1; c[0][3] = 1;
            c[1][0] = 1; c[1][1] = 1; c[1][2] =-1; c[1][3] =-1;
            c[2][0] = 1; c[2][1] =-1; c[2][2] =-1; c[2][3] = 1;
            c[3][0] = 1; c[3][1] =-1; c[3][2] = 1; c[3][3] =-1;
        }
        else if (name == "C2v") {
            order = 4;
            irs[0] = "a1"; irs[1] = "a2"; irs[2] = "b1"; irs[3] = "b2";
            ops[0] = "e";  ops[1] = "c2z";ops[2] = "sxz";ops[3] = "syz";
            c[0][0] = 1; c[0][1] = 1; c[0][2] = 1; c[0][3] = 1;
            c[1][0] = 1; c[1][1] = 1; c[1][2] =-1; c[1][3] =-1;
            c[2][0] = 1; c[2][1] =-1; c[2][2] = 1; c[2][3] =-1;
            c[3][0] = 1; c[3][1] =-1; c[3][2] =-1; c[3][3] = 1;
        }
        else if (name == "D2") {
            order = 4;
            irs[0] = "a1"; irs[1] = "b1"; irs[2] = "b2"; irs[3] = "b3";
            ops[0] = "e";  ops[1] = "c2z";ops[2] = "c2y";ops[3] = "c2x";
            c[0][0] = 1; c[0][1] = 1; c[0][2] = 1; c[0][3] = 1;
            c[1][0] = 1; c[1][1] = 1; c[1][2] =-1; c[1][3] =-1;
            c[2][0] = 1; c[2][1] =-1; c[2][2] = 1; c[2][3] =-1;
            c[3][0] = 1; c[3][1] =-1; c[3][2] =-1; c[3][3] = 1;
        }
        else if (name == "D2h") {
            order = 8;
            irs[0] = "ag"; irs[1] = "au"; irs[2] = "b1g"; irs[3] = "b1u";
            irs[4] = "b2g";irs[5] = "b2u";irs[6] = "b3g"; irs[7] = "b3u";
            ops[0] = "e";  ops[1] = "c2z";ops[2] = "c2y"; ops[3] = "c2x";
            ops[4] = "i";  ops[5] = "sxy";ops[6] = "sxz"; ops[7] = "syz";
            c[0][0] = 1;   c[0][1] = 1;   c[0][2] = 1;   c[0][3] = 1;   c[0][4] = 1;   c[0][5] = 1;   c[0][6] = 1;   c[0][7] = 1;
            c[1][0] = 1;   c[1][1] = 1;   c[1][2] = 1;   c[1][3] = 1;   c[1][4] =-1;   c[1][5] =-1;   c[1][6] =-1;   c[1][7] =-1;
            c[2][0] = 1;   c[2][1] = 1;   c[2][2] =-1;   c[2][3] =-1;   c[2][4] = 1;   c[2][5] = 1;   c[2][6] =-1;   c[2][7] =-1;
            c[3][0] = 1;   c[3][1] = 1;   c[3][2] =-1;   c[3][3] =-1;   c[3][4] =-1;   c[3][5] =-1;   c[3][6] = 1;   c[3][7] = 1;
            c[4][0] = 1;   c[4][1] =-1;   c[4][2] = 1;   c[4][3] =-1;   c[4][4] = 1;   c[4][5] =-1;   c[4][6] = 1;   c[4][7] =-1;
            c[5][0] = 1;   c[5][1] =-1;   c[5][2] = 1;   c[5][3] =-1;   c[5][4] =-1;   c[5][5] = 1;   c[5][6] =-1;   c[5][7] = 1;
            c[6][0] = 1;   c[6][1] =-1;   c[6][2] =-1;   c[6][3] = 1;   c[6][4] = 1;   c[6][5] =-1;   c[6][6] =-1;   c[6][7] = 1;
            c[7][0] = 1;   c[7][1] =-1;   c[7][2] =-1;   c[7][3] = 1;   c[7][4] =-1;   c[7][5] = 1;   c[7][6] = 1;   c[7][7] =-1;
        }
        else {
            throw "PointGroup: unknown group";
        }
    }

    /// Copy constructor
    PointGroup(const PointGroup& other) {
        *this = other;
    }

    /// Assignment
    PointGroup& operator=(const PointGroup& other) {
        if (this != &other) {
            name = other.name;
            order = other.order;
            for (int ir=0; ir<order; ++ir) {
                irs[ir] = other.irs[ir];
                ops[ir] = other.ops[ir];
                for (int op=0; op<order; ++op) {
                    c[ir][op] = other.c[ir][op];
                }
            }
        }
        return *this;
    }

    /// Destructor
    virtual ~PointGroup() {}

    /// Returns irreducible representation corresponding to product
    int irmul(int ir1, int ir2) const {
        return ir1 ^ ir2;
    }

    /// Applies group operator number op (0,1,...,order-1) to point
    coordT apply(int op, const coordT& r) const {
        return apply(ops[op], r);
    }

    /// Applies named operator (e, c2z, c2y, c2x, sxy, sxz, syz, i) to point
    static coordT apply(const std::string& op, const coordT& r)  {
        const double x=r[0], y=r[1], z=r[2];
        coordT q;
        if (op == "e") {
            q[0]=x;  q[1]=y;  q[2]=z;
        }
        else if (op == "c2z") {
            q[0]=-x; q[1]=-y; q[2]=z;
        }
        else if (op == "c2y") {
            q[0]=-x; q[1]=y;  q[2]=-z;
        }
        else if (op == "c2x") {
            q[0]=x;  q[1]=-y; q[2]=-z;
        }
        else if (op == "sxy") {
            q[0]=x;  q[1]=y;  q[2]=-z;
        }
        else if (op == "sxz") {
            q[0]=x;  q[1]=-y; q[2]=z;
        }
        else if (op == "syz") {
            q[0]=-x; q[1]=y;  q[2]=z;
        }
        else if (op == "i") {
            q[0]=-x; q[1]=-y; q[2]=-z;
        }
        else {
            throw "PointGroup: apply_op_by_name: unknown operator name";
        }
        return q;
    }

    /// Returns the irrep of the Cartesian axis (0, 1, 2 = x, y, z)
    int cart_ir(int axis) const {
        coordT r(0.0);
        r[axis] = 1.0;
        for (int ir=0; ir<order; ++ir) {
            double sum = 0.0;
            for (int op=0; op<order; ++op) {
                sum += apply(op,r)[axis] * c[ir][op];
            }
            sum /= order;
            if (sum > 0.9) {
                return ir;
            }
        }
        throw "PointGroup: cart_ir: problem identifying axis";
    }

    /// Returns irreducible cell

    /// The irreducible cell is defined by considering the cube
    /// -1<x<1, -1<y<1, -1<z<1.  Use symmetry operations to map
    /// negative coordinates onto positive coordinates if possible.
    /// This leads to the unique cell for D2h being the cube 0<x<1,
    /// 0<y<1, 0<z<1.
    ///
    /// If a coordinate in the cell is positive, it means that the
    /// irreducible cell has positive values of that coordinate.
    coordT ircell() const {
        double xmin=1.0, ymin=1.0, zmin=1.0;
        // Loop thru corners of the cube
        for (int x=-1; x<=1; x+=2) {
            for (int y=-1; y<=1; y+=2) {
                for (int z=-1; z<=1; z+=2) {
                    // Find the most positive corner it can be mapped into
                    double rx=x, ry=y, rz=z;
                    for (int op=0; op<order; ++op) {
                        coordT r;
                        r[0] = rx; r[1] = ry; r[2] = rz;
                        coordT q = apply(op,r);
                        double xx = q[0], yy=q[1], zz = q[2];
                        if ((xx>rx) || (xx==rx && yy>ry) || (xx==rx && yy==ry && zz>rz)) {
                            rx=xx; ry=yy; rz=zz;
                        }
                    }
                    xmin = std::min(xmin,rx); ymin = std::min(ymin,ry); zmin = std::min(zmin,rz);
                }
            }
        }

        coordT r;
        r[0] = xmin; r[1] = ymin; r[2] = zmin;
        return r;
    }

    static bool test() {
        std::cout << PointGroup("C1");
        std::cout << PointGroup("C2");
        std::cout << PointGroup("Ci");
        std::cout << PointGroup("Cs");
        std::cout << PointGroup("C2h");
        std::cout << PointGroup("C2v");
        std::cout << PointGroup("D2");
        std::cout << PointGroup("D2h");
        return true;
    }

    const std::string& get_name() const {
        return name;
    }

    int get_order() const {
        return order;
    }

    const std::string& get_ir_name(int ir) const {
        return irs[ir];
    }

    const std::string& get_op_name(int op) const {
        return ops[op];
    }

    int ctable(int ir, int op) const {
        return c[ir][op];
    }
};

std::ostream& operator<<(std::ostream& s, const PointGroup& g) {
    int order = g.get_order();
    s << "\n";
    s << "Group " << g.get_name() << " - irreducible cell " << g.ircell() << "\n";
    s << "---------\n";
    s << "        ";
    for (int op=0; op<order; ++op)
        s << " " << std::setw(3) << g.get_op_name(op);
    s << "\n";
    s << "        ";
    for (int op=0; op<order; ++op)
        s << " ---";

    s << "\n";

    int irx = g.cart_ir(0);
    int iry = g.cart_ir(1);
    int irz = g.cart_ir(2);

    for (int ir=0; ir<order; ++ir) {
        s << "  " << std::left << std::setw(3) << g.get_ir_name(ir) << std::right << "  ";
        for (int op=0; op<order; ++op) {
            s << " " << std::setw(3) << g.ctable(ir, op);
        }
        if (ir == irx) s << "   x";
        if (ir == iry) s << "   y";
        if (ir == irz) s << "   z";
        s << "\n";
    }
    s << "\n";
    return s;
}

#endif
