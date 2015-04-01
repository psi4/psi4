/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_TRIVIALNODES_HH_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_TRIVIALNODES_HH_

#include "Node.h"
namespace psi{
namespace LibMolecule{
class Atom;class H: public Node{
    public:
       unsigned Z()const{return 1;}
       H(int i=0):Node(i,FxnGrpType::H){}
};
class He: public Node{
    public:
       unsigned Z()const{return 2;}
       He(int i=0):Node(i,FxnGrpType::He){}
};
class Li: public Node{
    public:
       unsigned Z()const{return 3;}
       Li(int i=0):Node(i,FxnGrpType::Li){}
};
class Be: public Node{
    public:
       unsigned Z()const{return 4;}
       Be(int i=0):Node(i,FxnGrpType::Be){}
};
class B: public Node{
    public:
       unsigned Z()const{return 5;}
       B(int i=0):Node(i,FxnGrpType::B){}
};
class C: public Node{
    public:
       unsigned Z()const{return 6;}
       C(int i=0):Node(i,FxnGrpType::C){}
};
class N: public Node{
    public:
       unsigned Z()const{return 7;}
       N(int i=0):Node(i,FxnGrpType::N){}
};
class O: public Node{
    public:
       unsigned Z()const{return 8;}
       O(int i=0):Node(i,FxnGrpType::O){}
};
class F: public Node{
    public:
       unsigned Z()const{return 9;}
       F(int i=0):Node(i,FxnGrpType::F){}
};
class Ne: public Node{
    public:
       unsigned Z()const{return 10;}
       Ne(int i=0):Node(i,FxnGrpType::Ne){}
};
class Na: public Node{
    public:
       unsigned Z()const{return 11;}
       Na(int i=0):Node(i,FxnGrpType::Na){}
};
class Mg: public Node{
    public:
       unsigned Z()const{return 12;}
       Mg(int i=0):Node(i,FxnGrpType::Mg){}
};
class Al: public Node{
    public:
       unsigned Z()const{return 13;}
       Al(int i=0):Node(i,FxnGrpType::Al){}
};
class Si: public Node{
    public:
       unsigned Z()const{return 14;}
       Si(int i=0):Node(i,FxnGrpType::Si){}
};
class P: public Node{
    public:
       unsigned Z()const{return 15;}
       P(int i=0):Node(i,FxnGrpType::P){}
};
class S: public Node{
    public:
       unsigned Z()const{return 16;}
       S(int i=0):Node(i,FxnGrpType::S){}
};
class Cl: public Node{
    public:
       unsigned Z()const{return 17;}
       Cl(int i=0):Node(i,FxnGrpType::Cl){}
};
class Ar: public Node{
    public:
       unsigned Z()const{return 18;}
       Ar(int i=0):Node(i,FxnGrpType::Ar){}
};
class K: public Node{
    public:
       unsigned Z()const{return 19;}
       K(int i=0):Node(i,FxnGrpType::K){}
};
class Ca: public Node{
    public:
       unsigned Z()const{return 20;}
       Ca(int i=0):Node(i,FxnGrpType::Ca){}
};
class Sc: public Node{
    public:
       unsigned Z()const{return 21;}
       Sc(int i=0):Node(i,FxnGrpType::Sc){}
};
class Ti: public Node{
    public:
       unsigned Z()const{return 22;}
       Ti(int i=0):Node(i,FxnGrpType::Ti){}
};
class V: public Node{
    public:
       unsigned Z()const{return 23;}
       V(int i=0):Node(i,FxnGrpType::V){}
};
class Cr: public Node{
    public:
       unsigned Z()const{return 24;}
       Cr(int i=0):Node(i,FxnGrpType::Cr){}
};
class Mn: public Node{
    public:
       unsigned Z()const{return 25;}
       Mn(int i=0):Node(i,FxnGrpType::Mn){}
};
class Fe: public Node{
    public:
       unsigned Z()const{return 26;}
       Fe(int i=0):Node(i,FxnGrpType::Fe){}
};
class Co: public Node{
    public:
       unsigned Z()const{return 27;}
       Co(int i=0):Node(i,FxnGrpType::Co){}
};
class Ni: public Node{
    public:
       unsigned Z()const{return 28;}
       Ni(int i=0):Node(i,FxnGrpType::Ni){}
};
class Cu: public Node{
    public:
       unsigned Z()const{return 29;}
       Cu(int i=0):Node(i,FxnGrpType::Cu){}
};
class Zn: public Node{
    public:
       unsigned Z()const{return 30;}
       Zn(int i=0):Node(i,FxnGrpType::Zn){}
};
class Ga: public Node{
    public:
       unsigned Z()const{return 31;}
       Ga(int i=0):Node(i,FxnGrpType::Ga){}
};
class Ge: public Node{
    public:
       unsigned Z()const{return 32;}
       Ge(int i=0):Node(i,FxnGrpType::Ge){}
};
class As: public Node{
    public:
       unsigned Z()const{return 33;}
       As(int i=0):Node(i,FxnGrpType::As){}
};
class Se: public Node{
    public:
       unsigned Z()const{return 34;}
       Se(int i=0):Node(i,FxnGrpType::Se){}
};
class Br: public Node{
    public:
       unsigned Z()const{return 35;}
       Br(int i=0):Node(i,FxnGrpType::Br){}
};
class Kr: public Node{
    public:
       unsigned Z()const{return 36;}
       Kr(int i=0):Node(i,FxnGrpType::Kr){}
};
class Rb: public Node{
    public:
       unsigned Z()const{return 37;}
       Rb(int i=0):Node(i,FxnGrpType::Rb){}
};
class Sr: public Node{
    public:
       unsigned Z()const{return 38;}
       Sr(int i=0):Node(i,FxnGrpType::Sr){}
};
class Y: public Node{
    public:
       unsigned Z()const{return 39;}
       Y(int i=0):Node(i,FxnGrpType::Y){}
};
class Zr: public Node{
    public:
       unsigned Z()const{return 40;}
       Zr(int i=0):Node(i,FxnGrpType::Zr){}
};
class Nb: public Node{
    public:
       unsigned Z()const{return 41;}
       Nb(int i=0):Node(i,FxnGrpType::Nb){}
};
class Mo: public Node{
    public:
       unsigned Z()const{return 42;}
       Mo(int i=0):Node(i,FxnGrpType::Mo){}
};
class Tc: public Node{
    public:
       unsigned Z()const{return 43;}
       Tc(int i=0):Node(i,FxnGrpType::Tc){}
};
class Ru: public Node{
    public:
       unsigned Z()const{return 44;}
       Ru(int i=0):Node(i,FxnGrpType::Ru){}
};
class Rh: public Node{
    public:
       unsigned Z()const{return 45;}
       Rh(int i=0):Node(i,FxnGrpType::Rh){}
};
class Pd: public Node{
    public:
       unsigned Z()const{return 46;}
       Pd(int i=0):Node(i,FxnGrpType::Pd){}
};
class Ag: public Node{
    public:
       unsigned Z()const{return 47;}
       Ag(int i=0):Node(i,FxnGrpType::Ag){}
};
class Cd: public Node{
    public:
       unsigned Z()const{return 48;}
       Cd(int i=0):Node(i,FxnGrpType::Cd){}
};
class In: public Node{
    public:
       unsigned Z()const{return 49;}
       In(int i=0):Node(i,FxnGrpType::In){}
};
class Sn: public Node{
    public:
       unsigned Z()const{return 50;}
       Sn(int i=0):Node(i,FxnGrpType::Sn){}
};
class Sb: public Node{
    public:
       unsigned Z()const{return 51;}
       Sb(int i=0):Node(i,FxnGrpType::Sb){}
};
class Te: public Node{
    public:
       unsigned Z()const{return 52;}
       Te(int i=0):Node(i,FxnGrpType::Te){}
};
class I: public Node{
    public:
       unsigned Z()const{return 53;}
       I(int i=0):Node(i,FxnGrpType::I){}
};
class Xe: public Node{
    public:
       unsigned Z()const{return 54;}
       Xe(int i=0):Node(i,FxnGrpType::Xe){}
};
class Cs: public Node{
    public:
       unsigned Z()const{return 55;}
       Cs(int i=0):Node(i,FxnGrpType::Cs){}
};
class Ba: public Node{
    public:
       unsigned Z()const{return 56;}
       Ba(int i=0):Node(i,FxnGrpType::Ba){}
};
class La: public Node{
    public:
       unsigned Z()const{return 57;}
       La(int i=0):Node(i,FxnGrpType::La){}
};
class Ce: public Node{
    public:
       unsigned Z()const{return 58;}
       Ce(int i=0):Node(i,FxnGrpType::Ce){}
};
class Pr: public Node{
    public:
       unsigned Z()const{return 59;}
       Pr(int i=0):Node(i,FxnGrpType::Pr){}
};
class Nd: public Node{
    public:
       unsigned Z()const{return 60;}
       Nd(int i=0):Node(i,FxnGrpType::Nd){}
};
class Pm: public Node{
    public:
       unsigned Z()const{return 61;}
       Pm(int i=0):Node(i,FxnGrpType::Pm){}
};
class Sm: public Node{
    public:
       unsigned Z()const{return 62;}
       Sm(int i=0):Node(i,FxnGrpType::Sm){}
};
class Eu: public Node{
    public:
       unsigned Z()const{return 63;}
       Eu(int i=0):Node(i,FxnGrpType::Eu){}
};
class Gd: public Node{
    public:
       unsigned Z()const{return 64;}
       Gd(int i=0):Node(i,FxnGrpType::Gd){}
};
class Tb: public Node{
    public:
       unsigned Z()const{return 65;}
       Tb(int i=0):Node(i,FxnGrpType::Tb){}
};
class Dy: public Node{
    public:
       unsigned Z()const{return 66;}
       Dy(int i=0):Node(i,FxnGrpType::Dy){}
};
class Ho: public Node{
    public:
       unsigned Z()const{return 67;}
       Ho(int i=0):Node(i,FxnGrpType::Ho){}
};
class Er: public Node{
    public:
       unsigned Z()const{return 68;}
       Er(int i=0):Node(i,FxnGrpType::Er){}
};
class Tm: public Node{
    public:
       unsigned Z()const{return 69;}
       Tm(int i=0):Node(i,FxnGrpType::Tm){}
};
class Yb: public Node{
    public:
       unsigned Z()const{return 70;}
       Yb(int i=0):Node(i,FxnGrpType::Yb){}
};
class Lu: public Node{
    public:
       unsigned Z()const{return 71;}
       Lu(int i=0):Node(i,FxnGrpType::Lu){}
};
class Hf: public Node{
    public:
       unsigned Z()const{return 72;}
       Hf(int i=0):Node(i,FxnGrpType::Hf){}
};
class Ta: public Node{
    public:
       unsigned Z()const{return 73;}
       Ta(int i=0):Node(i,FxnGrpType::Ta){}
};
class W: public Node{
    public:
       unsigned Z()const{return 74;}
       W(int i=0):Node(i,FxnGrpType::W){}
};
class Re: public Node{
    public:
       unsigned Z()const{return 75;}
       Re(int i=0):Node(i,FxnGrpType::Re){}
};
class Os: public Node{
    public:
       unsigned Z()const{return 76;}
       Os(int i=0):Node(i,FxnGrpType::Os){}
};
class Ir: public Node{
    public:
       unsigned Z()const{return 77;}
       Ir(int i=0):Node(i,FxnGrpType::Ir){}
};
class Pt: public Node{
    public:
       unsigned Z()const{return 78;}
       Pt(int i=0):Node(i,FxnGrpType::Pt){}
};
class Au: public Node{
    public:
       unsigned Z()const{return 79;}
       Au(int i=0):Node(i,FxnGrpType::Au){}
};
class Hg: public Node{
    public:
       unsigned Z()const{return 80;}
       Hg(int i=0):Node(i,FxnGrpType::Hg){}
};
class Tl: public Node{
    public:
       unsigned Z()const{return 81;}
       Tl(int i=0):Node(i,FxnGrpType::Tl){}
};
class Pb: public Node{
    public:
       unsigned Z()const{return 82;}
       Pb(int i=0):Node(i,FxnGrpType::Pb){}
};
class Bi: public Node{
    public:
       unsigned Z()const{return 83;}
       Bi(int i=0):Node(i,FxnGrpType::Bi){}
};
class Po: public Node{
    public:
       unsigned Z()const{return 84;}
       Po(int i=0):Node(i,FxnGrpType::Po){}
};
class At: public Node{
    public:
       unsigned Z()const{return 85;}
       At(int i=0):Node(i,FxnGrpType::At){}
};
class Rn: public Node{
    public:
       unsigned Z()const{return 86;}
       Rn(int i=0):Node(i,FxnGrpType::Rn){}
};
class Fr: public Node{
    public:
       unsigned Z()const{return 87;}
       Fr(int i=0):Node(i,FxnGrpType::Fr){}
};
class Ra: public Node{
    public:
       unsigned Z()const{return 88;}
       Ra(int i=0):Node(i,FxnGrpType::Ra){}
};
class Ac: public Node{
    public:
       unsigned Z()const{return 89;}
       Ac(int i=0):Node(i,FxnGrpType::Ac){}
};
class Th: public Node{
    public:
       unsigned Z()const{return 90;}
       Th(int i=0):Node(i,FxnGrpType::Th){}
};
class Pa: public Node{
    public:
       unsigned Z()const{return 91;}
       Pa(int i=0):Node(i,FxnGrpType::Pa){}
};
class U: public Node{
    public:
       unsigned Z()const{return 92;}
       U(int i=0):Node(i,FxnGrpType::U){}
};
class Np: public Node{
    public:
       unsigned Z()const{return 93;}
       Np(int i=0):Node(i,FxnGrpType::Np){}
};
class Pu: public Node{
    public:
       unsigned Z()const{return 94;}
       Pu(int i=0):Node(i,FxnGrpType::Pu){}
};
class Am: public Node{
    public:
       unsigned Z()const{return 95;}
       Am(int i=0):Node(i,FxnGrpType::Am){}
};
class Cm: public Node{
    public:
       unsigned Z()const{return 96;}
       Cm(int i=0):Node(i,FxnGrpType::Cm){}
};
class Bk: public Node{
    public:
       unsigned Z()const{return 97;}
       Bk(int i=0):Node(i,FxnGrpType::Bk){}
};
class Cf: public Node{
    public:
       unsigned Z()const{return 98;}
       Cf(int i=0):Node(i,FxnGrpType::Cf){}
};
class Es: public Node{
    public:
       unsigned Z()const{return 99;}
       Es(int i=0):Node(i,FxnGrpType::Es){}
};
class Fm: public Node{
    public:
       unsigned Z()const{return 100;}
       Fm(int i=0):Node(i,FxnGrpType::Fm){}
};
class Md: public Node{
    public:
       unsigned Z()const{return 101;}
       Md(int i=0):Node(i,FxnGrpType::Md){}
};
class No: public Node{
    public:
       unsigned Z()const{return 102;}
       No(int i=0):Node(i,FxnGrpType::No){}
};
class Lr: public Node{
    public:
       unsigned Z()const{return 103;}
       Lr(int i=0):Node(i,FxnGrpType::Lr){}
};
class Rf: public Node{
    public:
       unsigned Z()const{return 104;}
       Rf(int i=0):Node(i,FxnGrpType::Rf){}
};
class Db: public Node{
    public:
       unsigned Z()const{return 105;}
       Db(int i=0):Node(i,FxnGrpType::Db){}
};
class Sg: public Node{
    public:
       unsigned Z()const{return 106;}
       Sg(int i=0):Node(i,FxnGrpType::Sg){}
};
class Bh: public Node{
    public:
       unsigned Z()const{return 107;}
       Bh(int i=0):Node(i,FxnGrpType::Bh){}
};
class Hs: public Node{
    public:
       unsigned Z()const{return 108;}
       Hs(int i=0):Node(i,FxnGrpType::Hs){}
};
class Mt: public Node{
    public:
       unsigned Z()const{return 109;}
       Mt(int i=0):Node(i,FxnGrpType::Mt){}
};
class Ds: public Node{
    public:
       unsigned Z()const{return 110;}
       Ds(int i=0):Node(i,FxnGrpType::Ds){}
};
class Rg: public Node{
    public:
       unsigned Z()const{return 111;}
       Rg(int i=0):Node(i,FxnGrpType::Rg){}
};
class Cn: public Node{
    public:
       unsigned Z()const{return 112;}
       Cn(int i=0):Node(i,FxnGrpType::Cn){}
};
class Uut: public Node{
    public:
       unsigned Z()const{return 113;}
       Uut(int i=0):Node(i,FxnGrpType::Uut){}
};
class Fl: public Node{
    public:
       unsigned Z()const{return 114;}
       Fl(int i=0):Node(i,FxnGrpType::Fl){}
};
class Uup: public Node{
    public:
       unsigned Z()const{return 115;}
       Uup(int i=0):Node(i,FxnGrpType::Uup){}
};
class Lv: public Node{
    public:
       unsigned Z()const{return 116;}
       Lv(int i=0):Node(i,FxnGrpType::Lv){}
};
class Uus: public Node{
    public:
       unsigned Z()const{return 117;}
       Uus(int i=0):Node(i,FxnGrpType::Uus){}
};
class Uuo: public Node{
    public:
       unsigned Z()const{return 118;}
       Uuo(int i=0):Node(i,FxnGrpType::Uuo){}
};
}}//End namespaces
#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_TRIVIALNODES_HH_*/
