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

  $Id: test_problems.h 1856 2010-04-06 14:03:52Z mgr522 $
*/

/** \file atom.h
    \brief

*/

#ifndef MADNESS_INTERIOR_BC_ATOM_H__INCLUDED
#define MADNESS_INTERIOR_BC_ATOM_H__INCLUDED

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include "basisfunction.h"
#include <string>

using namespace madness;

typedef std::shared_ptr<GaussianBF> BasisFunc;

/** \brief Abstract Atom class. */
class Atom {
    private:
        Atom() {}

    protected:
        std::vector<BasisFunc> basis;
        Vector<double, 3> center;

    public:
        /// \brief Sets up the basis functions for the atom.
        Atom(const Vector<double, 3> &c)
            : basis(0) {

            center[0] = c[0];
            center[1] = c[1];
            center[2] = c[2];
        }

        Atom(const Vector<double, 3> &c, int n)
            : basis(n) {

            center[0] = c[0];
            center[1] = c[1];
            center[2] = c[2];
        }

        virtual ~Atom() {
            for(std::vector<BasisFunc>::iterator iter = basis.begin();
                iter != basis.end();
                ++iter) {

                delete *iter;
            }
        }

        BasisFunc getBasisFunc(unsigned int n) {
            assert(n >= 0 && n < basis.size());
            return basis[n];
        }

        const Vector<double, 3> &getCenter() const {
            return center;
        }

        virtual int dimBasis() const = 0;
};

/** \brief Hydrogen atom */
class Hydrogen : public Atom {
    private:
        Hydrogen() : Atom(0) {}

    public:
        Hydrogen(Vector<double, 3> &c)
            : Atom(c, 2) {

            std::vector<double> c1(3), c2(1);
            std::vector<double> e1(3), e2(1);

            c1[0] = 0.033494604338;
            c1[1] = 0.234726953484;
            c1[2] = 0.813757326146;
            e1[0] = 18.7311370;
            e1[1] = 2.8253944;
            e1[2] = 0.6401217;

            c2[0] = 1.0;
            e2[0] = 0.1612778;

            basis[0] = BasisFunc(new SBF(c1, e1, center));
            basis[1] = BasisFunc(new SBF(c2, e2, center));
        }

        virtual int dimBasis() const { return 2; }
};

/** \brief Carbon atom */
class Carbon : public Atom {
    private:
        Carbon() : Atom(0) {}

    public:
        Carbon(Vector<double, 3> &c)
            : Atom(c, 21) {

            std::vector<double> c1(6), e1(6);
            std::vector<double> c2s(3), c2p(3), e2(3);
            std::vector<double> c3s(1), c3p(1), e3(1);
            std::vector<double> c4(1), e4(1);
            std::vector<double> c5(1), e5(1);

            c1[0] = 0.001834737132;
            c1[1] = 0.014037322813;
            c1[2] = 0.068842622264;
            c1[3] = 0.232184443216;
            c1[4] = 0.467941348435;
            c1[5] = 0.362311985337;
            e1[0] = 3047.5248800;
            e1[1] = 457.3695180;
            e1[2] = 103.9486850;
            e1[3] = 29.2101553;
            e1[4] = 9.2866630;
            e1[5] = 3.1639270;

            c2s[0] = -0.119332419775;
            c2s[1] = -0.160854151696;
            c2s[2] = 1.143456437840;
            c2p[0] = 0.068999066591;
            c2p[1] = 0.316423960957;
            c2p[2] = 0.744308290898;
            e2[0] = 7.8682723;
            e2[1] = 1.8812885;
            e2[2] = 0.5442493;

            c3s[0] = 1.0;
            c3p[0] = 1.0;
            e3[0] = 0.1687145;

            c4[0] = 1.0;
            e4[0] = 1.6000000;

            c5[0] = 1.0;
            e5[0] = 0.4000000;

            basis[0] = BasisFunc(new SBF(c1, e1, center));
            basis[1] = BasisFunc(new SBF(c2s, e2, center));
            basis[2] = BasisFunc(new PBF(c2p, e2, center, PBF::X));
            basis[3] = BasisFunc(new PBF(c2p, e2, center, PBF::Y));
            basis[4] = BasisFunc(new PBF(c2p, e2, center, PBF::Z));
            basis[5] = BasisFunc(new SBF(c3s, e3, center));
            basis[6] = BasisFunc(new PBF(c3p, e3, center, PBF::X));
            basis[7] = BasisFunc(new PBF(c3p, e3, center, PBF::Y));
            basis[8] = BasisFunc(new PBF(c3p, e3, center, PBF::Z));
            basis[9] = BasisFunc(new DBF(c4, e4, center, DBF::XX));
            basis[10] = BasisFunc(new DBF(c4, e4, center, DBF::YY));
            basis[11] = BasisFunc(new DBF(c4, e4, center, DBF::ZZ));
            basis[12] = BasisFunc(new DBF(c4, e4, center, DBF::XY));
            basis[13] = BasisFunc(new DBF(c4, e4, center, DBF::XZ));
            basis[14] = BasisFunc(new DBF(c4, e4, center, DBF::YZ));
            basis[15] = BasisFunc(new DBF(c5, e5, center, DBF::XX));
            basis[16] = BasisFunc(new DBF(c5, e5, center, DBF::YY));
            basis[17] = BasisFunc(new DBF(c5, e5, center, DBF::ZZ));
            basis[18] = BasisFunc(new DBF(c5, e5, center, DBF::XY));
            basis[19] = BasisFunc(new DBF(c5, e5, center, DBF::XZ));
            basis[20] = BasisFunc(new DBF(c5, e5, center, DBF::YZ));
        }

        virtual int dimBasis() const { return 21; }
};

/** \brief Silicon atom */
class Silicon : public Atom {
    private:
        Silicon() : Atom(0) {}

    public:
        Silicon(Vector<double, 3> &c)
            : Atom(c, 25) {

            std::vector<double> c1(6), e1(6);
            std::vector<double> c2s(6), c2p(6), e2(6);
            std::vector<double> c3s(3), c3p(3), e3(3);
            std::vector<double> c4s(1), c4p(1), e4(1);
            std::vector<double> c5(1), e5(1);
            std::vector<double> c6(1), e6(1);

            c1[0] = 0.001949239405;
            c1[1] = 0.014855895466;
            c1[2] = 0.072568877851;
            c1[3] = 0.245654925022;
            c1[4] = 0.486059851646;
            c1[5] = 0.325719900585;
            e1[0] = 16192.1000000;
            e1[1] = 2436.0900000;
            e1[2] = 556.0010000;
            e1[3] = 156.8130000;
            e1[4] = 50.1692000;
            e1[5] = 17.0300000;

            c2s[0] = -0.002829912441;
            c2s[1] = -0.036073731114;
            c2s[2] = -0.116808100749;
            c2s[3] = 0.093576880711;
            c2s[4] = 0.601705518979;
            c2s[5] = 0.422072364043;
            c2p[0] = 0.004433341604;
            c2p[1] = 0.032440211739;
            c2p[2] = 0.133719048387;
            c2p[3] = 0.326780118247;
            c2p[4] = 0.451139163247;
            c2p[5] = 0.264105095568;
            e2[0] = 293.3500000;
            e2[1] = 70.1173000;
            e2[2] = 22.4301000;
            e2[3] = 8.1942500;
            e2[4] = 3.1476800;
            e2[5] = 1.2151500;

            c3s[0] = -0.240599186106;
            c3s[1] = 0.073795050368;
            c3s[2] = 1.040936478740;
            c3p[0] = -0.015177406538;
            c3p[1] = 0.275139118516;
            c3p[2] = 0.783008337279;
            e3[0] = 1.6537000;
            e3[1] = 0.5407600;
            e3[2] = 0.2044060;

            c4s[0] = 1.0;
            c4p[0] = 1.0;
            e4[0] = 0.0723837;

            c5[0] = 1.0;
            e5[0] = 0.7900000;

            c6[0] = 1.0;
            e6[0] = 0.1975000;

            basis[0] = BasisFunc(new SBF(c1, e1, center));
            basis[1] = BasisFunc(new SBF(c2s, e2, center));
            basis[2] = BasisFunc(new PBF(c2p, e2, center, PBF::X));
            basis[3] = BasisFunc(new PBF(c2p, e2, center, PBF::Y));
            basis[4] = BasisFunc(new PBF(c2p, e2, center, PBF::Z));
            basis[5] = BasisFunc(new SBF(c3s, e3, center));
            basis[6] = BasisFunc(new PBF(c3p, e3, center, PBF::X));
            basis[7] = BasisFunc(new PBF(c3p, e3, center, PBF::Y));
            basis[8] = BasisFunc(new PBF(c3p, e3, center, PBF::Z));
            basis[9] = BasisFunc(new SBF(c4s, e4, center));
            basis[10] = BasisFunc(new PBF(c4p, e4, center, PBF::X));
            basis[11] = BasisFunc(new PBF(c4p, e4, center, PBF::Y));
            basis[12] = BasisFunc(new PBF(c4p, e4, center, PBF::Z));
            basis[13] = BasisFunc(new DBF(c5, e5, center, DBF::XX));
            basis[14] = BasisFunc(new DBF(c5, e5, center, DBF::YY));
            basis[15] = BasisFunc(new DBF(c5, e5, center, DBF::ZZ));
            basis[16] = BasisFunc(new DBF(c5, e5, center, DBF::XY));
            basis[17] = BasisFunc(new DBF(c5, e5, center, DBF::XZ));
            basis[18] = BasisFunc(new DBF(c5, e5, center, DBF::YZ));
            basis[19] = BasisFunc(new DBF(c6, e6, center, DBF::XX));
            basis[20] = BasisFunc(new DBF(c6, e6, center, DBF::YY));
            basis[21] = BasisFunc(new DBF(c6, e6, center, DBF::ZZ));
            basis[22] = BasisFunc(new DBF(c6, e6, center, DBF::XY));
            basis[23] = BasisFunc(new DBF(c6, e6, center, DBF::XZ));
            basis[24] = BasisFunc(new DBF(c6, e6, center, DBF::YZ));
        }

        virtual int dimBasis() const { return 25; }
};

#endif // MADNESS_INTERIOR_BC_ATOM_H__INCLUDED
