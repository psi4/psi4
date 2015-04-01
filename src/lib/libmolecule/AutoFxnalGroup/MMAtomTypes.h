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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_MMATOMTYPES_H_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_MMATOMTYPES_H_

namespace psi{
namespace LibMolecule{
/** \brief The atom types my code can find
 *         (look at the *.cc for descriptions)
 *
 *  I have tried to be systematic in naming these things.  The type of
 *  the atoms are just their atomic symbol.
 *
 *  The types of primitives account for their order and their (likely)
 *  bond multiplicity.  Atoms that form only single bonds
 *  follow the pattern \f$XN\f$ where \f$X\f$ is the atomic symbol
 *  of the central atom and \f$N\f$ is the order (a carbon forming four
 *  bonds, to three hydrogens is of order 1, a so-called primary carbon
 *  Put another way the order is the number of heavy atoms something is bonded
 *  to).  For
 *  multiple bonds \f$X\f$ is expanded with a label of DB or TB depending
 *  on whether or not the bond-multiplicity is that of a DB or a TB.  For
 *  a carbon that makes three bonds, two of which are to hydrogens, the
 *  carbon is of double-bond character, and labeled CDB and the order is
 *  1 (a primary double-bonded carbon).  The full symbol would be CDB1.
 *  For each primitive we need a systematic name for the atoms in it, this
 *  is done by appending the atomic symbol of the atom to the group's
 *  sybmol, e.g. the H in a methyl group is HC1 and the C in a primary
 *  double-bond-like carbon is CCDB1.  For nitrogen (and possibly other
 *  atoms) it was necessary to include a charge we did this with the letter
 *  "P" for positive, e.g. the type of the ammonium ion is NP0.
 *
 *  When naming the derived groups originating from multiple bonds we switch
 *  the order of the bond-order and the atom identifications, e.g. a
 *  carbon-carbon double bond would be DBCC\f$N\f$ where \f$N\f$ is again
 *  the order.  A carbon in that double bond would be named CDBCC\f$N\f$.
 *  I originally considered the system: CCCDB1 for a carbon in a primary
 *  carbon-carbon double bond, etc. I thought that such a system was more
 *  natural; however the name of this group, CCDB1, obviously conflicts
 *  with the type of a primary, double-bond like carbon, also CCDB1.  As an
 *  added benefit, the accepted convention avoids having many of the same
 *  letter next to eachother, e.g. the 3 "C"'s in CCCDB1.  For DBCC1 we
 *  have two carbons that are not symmetric we label then numerically
 *  with 1 being the highest priority carbon, 2 the next highest etc.
 *  Priority is first to the atom with highest Z letting an R group
 *  have a Z of infinity, then we break ties by bond order i.e. the carbonyl
 *  oxygen in a carboxyl is higher priority than the hydroxyl oxygen.
 *
 *  Particularly for derived groups containing oxygen I had to get a bit
 *  creative in the names hopefully it's not too hard to figure out.
 *
 */
enum class FxnGrpType{
      NONE,TEMP,
      H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,
      Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,
      Sn,Sb,Te,I,Xe,Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,
      W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,At,Rn,Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,
      Cf,Es,Fm,Md,No,Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,Rg,Cn,Uut,Fl,Uup,Lv,Uus,Uuo,
      C0,
      CC0,HC0,
      C1,
      CC1,HC1,
      C2,
      CC2,HC2,
      C3,
      CC3,HC3,
      C4,
      CC4,
      CDB1,
      CCDB1,HCDB1,
      CDB2,
      CCDB2,HCDB2,
      CDB3,
      CCDB3,
      CTB1,
      CCTB1,HCTB1,
      CTB2,
      CCTB2,
      N0,
      NN0,HN0,
      N1,
      NN1,HN1,
      N2,
      NN2,HN2,
      N3,
      NN3,
      NDB1,
      NNDB1,HNDB1,
      NDB2,
      NNDB2,
      NTB,
      NNTB,
      NP0,
      NNP0,HNP0,
      NP1,
      NNP1,HNP1,
      NP2,
      NNP2,HNP2,
      NP3,
      NNP3,HNP3,
      NP4,
      NNP4,
      O0,
      OO0,HO0,
      O1,
      OO1,HO1,
      O2,
      OO2,
      ODB,
      OODB,
      S0,
      SS0,HS0,
      S1,
      SS1,HS1,
      S2,
      SS2,
      SDB,
      SSDB,
      DBCC0,
      CDBCC0,HDBCC0,
      DBCC1,
      C1DBCC1,H1DBCC1,C2DBCC1,H2DBCC1,
      DBCC2,
      CDBCC2,HDBCC2,
      DBCC2G,
      C1DBCC2G,HDBCC2G,C2DBCC2G,
      DBCC3,
      C1DBCC3,HDBCC3,C2DBCC3,
      DBCC4,
      CDBCC4,
      TBCC0,
      CTBCC0,HTBCC0,
      TBCC1,
      C1TBCC1,HTBCC1,C2TBCC1,
      TBCC2,
      CTBCC2,
      DBCN0,
      CDBCN0,NDBCN0,H1DBCN0,H2DBCN0,
      DBCN1N,
      CDBCN1N,NDBCN1N,HDBCN1N,
      DBCN1C,
      CDBCN1C,NDBCN1C,H1DBCN1C,H2DBCN1C,
      DBCN2,
      CDBCN2,NDBCN2,HDBCN2,
      DBCN2G,
      CDBCN2G,NDBCN2G,HDBCN2G,
      DBCN3,
      CDBCN3,NDBCN3,
      TBCN0,
      CTBCN0,NTBCN0,HTBCN0,
      TBCN1,
      CTBCN1,NTBCN1,
      DBCO0,
      CDBCO0,HDBCO0,ODBCO0,
      DBCO1,
      CDBCO1,HDBCO1,ODBCO1,
      DBCO2,
      CDBCO2,ODBCO2,
      CO2M,
      CCO2M,OCO2M,
      CO2,
      CCO2,O1CO2,O2CO2,
      CO2H,
      CCO2H,O1CO2H,O2CO2H,HCO2H,
      CO3,
      CCO3,O1CO3,O2CO3,
      CONH2,
      CCONH2,OCONH2,NCONH2,HCONH2,
      CHONH,
      CCHONH,OCHONH,NCHONH,H1CHONH,H2CHONH,
      CONH,
      CCONH,OCONH,NCONH,HCONH,
      CHON,
      CCHON,OCHON,NCHON,HCHON,
      CON,
      CCON,OCON,NCON,
      CONHCOH,
      C1CONHCOH,C2CONHCOH,O1CONHCOH,O2CONHCOH,NCONHCOH,H1CONHCOH,H2CONHCOH,
      CHONCOH,
      CCHONCOH,HCHONCOH,OCHONCOH,NCHONCOH,
      CONCOH,
      C1CONCOH,C2CONCOH,O1CONCOH,O2CONCOH,NCONCOH,HCONCOH,
      CONHCO,
      CCONHCO,OCONHCO,NCONHCO,HCONHCO,
      CONCO,
      CCONCO,OCONCO,NCONCO,
      OCH4,
      OOCH4,COCH4,H1OCH4,H2OCH4,
      OCH3,
      OOCH3,COCH3,HOCH3,
      O2H2,
      OO2H2,HO2H2,
      O2H,
      O1O2H,O2O2H,HO2H,
      OO,
      OOO,
      OCHOH,
      O1OCHOH,O2OCHOH,COCHOH,H1OCHOH,H2OCHOH,
      OCOH,
      O1OCOH,O2OCOH,COCOH,HOCOH,
      OCHO,
      OOCHO,COCHO,HOCHO,
      OCO,
      OOCO,COCO,
      COOO,
      CCOOO,OCOOO,
      NNN,
      N1NNN,N2NNN,N3NNN,
      DBNN0,
      NDBNN0,HDBNN0,
      DBNN1,
      N1DBNN1,N2DBNN1,HDBNN1,
      DBNN2,
      NDBNN2,
      OCN,
      COCN,NOCN,OOCN,
      NCO,
      NNCO,CNCO,ONCO,
      NOO,
      ONOO,NNOO,
      NO3,
      O1NO3,O2NO3,NNO3,
      NO,
      NNO,ONO,
      NO2,
      NNO2,O1NO2,O2NO2,
      IPR,
      C1IPR,C2IPR,H1IPR,H2IPR,
      Bz0,
      CBz0,HBz0,
      Bz1,
      C1Bz1,C2Bz1,C3Bz1,C4Bz1,H1Bz1,H2Bz1,H3Bz1,
      Bz2o,
      C1Bz2o,C2Bz2o,C3Bz2o,H1Bz2o,H2Bz2o,
      Bz2m,
      C1Bz2m,C2Bz2m,C3Bz2m,C4Bz2m,H1Bz2m,H2Bz2m,H3Bz2m,
      Bz2p,
      C1Bz2p,C2Bz2p,HBz2p,
      Bz6,
      CBz6,
      PNT,
      NPNT,HNPNT,CAPNT,HCAPNT,CPBPNT,OPBPNT,
      NCT,
      CNCT,ONCT,CANCT,HCANCT,NPBNCT,HPBNCT,
      GLYPNT,
      NGLYPNT,HNGLYPNT,CAGLYPNT,HCAGLYPNT,CPBGLYPNT,OPBGLYPNT,
      GLYNCT,
      NGLYNCT,HNGLYNCT,CAGLYNCT,HCAGLYNCT,CGLYNCT,OGLYNCT,
      AABB,
      CPBAABB,OPBAABB,CCAAABB,HCAAABB,NPBAABB,HPBAABB,
      GLY,
      CPBGLY,OPBGLY,CCAGLY,HCAGLY,NPBGLY,HPBGLY,
      ALAR,
      CALAR,HALAR,
      VALR,
      C1VALR,C2VALR,H1VALR,H2VALR,
      ILER,
      C1ILER,C2ILER,C3ILER,C4ILER,H1ILER,H2ILER,H3ILER,H4ILER,
      LEUR,
      C1LEUR,C2LEUR,C3LEUR,H1LEUR,H2LEUR,H3LEUR,
      METR,
      C1METR,C2METR,C3METR,SMETR,H1METR,H2METR,H3METR,
      PHER,
      C1PHER,C2PHER,C3PHER,C4PHER,C5PHER,H1PHER,H2PHER,H3PHER,H4PHER,
      TYRR,
      C1TYRR,C2TYRR,C3TYRR,C4TYRR,C5TYRR,OTYRR,H1TYRR,H2TYRR,H3TYRR,H4TYRR,
      TRPR,
      C1TRPR,C2TRPR,C3TRPR,C4TRPR,C5TRPR,C6TRPR,C7TRPR,C8TRPR,C9TRPR,
      H1TRPR,H2TRPR,H3TRPR,H4TRPR,H5TRPR,H6TRPR,H7TRPR,NTRPR,
      SERR,
      OSERR,CSERR,H1SERR,H2SERR
      //Need to end in a , to continue seamlessly into IndoleTypes.hh
      ,
      #include "IndoleTypes.hh"

};

///A handy function for printing out the FxnGrpType
std::string FxnGrpType2String(const FxnGrpType& Ty);

}}




#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_MMATOMTYPES_H_ */
