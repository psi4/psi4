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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_FXNALGROUP_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_FXNALGROUP_H_

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
#include "PsiMap.h"
#include "Parameter.h"
namespace psi{
namespace LibMolecule{

/** \brief The recognized functional groups
 *
 * - Primitive Groups: Here is what I mean by each term:
 *   - methane (primary-1 carbon) CH4
 *   - methyl (primary carbon) R-CH3
 *   - methene (secondary carbon) R-CH2-R`
 *   - alkenyl1 (primary alkenyl) R=CH2
 *   - methyne (tertiary carbon) R-(CH-R`)-R```
 *   - alkenyl2 (seconary alkenyl) R=CH-R'
 *   - alkynyl1 (primary alkenyl)  R <triple-bond> CH
 *   - carbon4 (quaternary carbon) R-((C-R`)-R``)-R```
 *   - alkenyl3 (tertiary alkenyl) R=(C-R')-R''
 *   - alkynyl2 (secondary alynyl) R <triple-bond>C-R'
 *   - water   (primary-1 oxygen) OH2
 *   - hydroxyl (priarmy oxygen) R-OH
 *   - oxygen2 (secondary oxygen) R-O-R`
 *   - oxydb   (primary oxygen doubl bond) R=O
 *   - ammonia (primary-1 amine) NH3
 *   - amine1 (primary amine) R-NH2
 *   - amine2 (secondary amine) R-HN-R`
 *   - nitrodb1 (primary N double bond) R=N-H
 *   - amine3 (tertiary amine) R-(N-R`)-R``
 *   - nitrodb2 (secondary double N bond) R-N=R`
 *   - nitrotb (primary triple N bond) R <triple bond> N
 * -Derived Groups.  These are functional groups that can be written as
 *   combinations of primitive groups or in some cases other derived groups.
 *    - hydrogen cyanide (alkynyl1)-(nitrotb)
 *    - nitrile R- (alkynyl2) -(nitrotb)
 *    - formaldehyde (alkenyl1)-(oxydb)
 *    - aldehyde R-(alkenyl2)-(oxydb)
 *    - carbonyl R-[(alkenyl3)-(oxydb)]-R`
 *    - carboxyl R-(carbonyl)-(hydroxyl)
 *    - hydroperoxy R-(secondary oxygen)-(hydroxyl)
 *    - peroxide (hydroxyl)-(hydroxyl)
 *    - methoxy R-(secondary oxygen)(methyl)
 *    - methanol (hydroxyl)-(methyl)
 *    - ccdb4 (quatanary double-bond)
 *            R-((alkenyl3)-R')-((alkenyl3)-R'')-R'''
 *    - ccdb3 (ternary double-bond)
 *            R-(alkenyl2)-((alkenyl3)-R')-R''
 *    - ethenyl2 (secondary ethenyl)(alkenyl1)-((alkenyl3)-R)-R''
 *    - ccdb2 (seconardy double bond)
 *            R-(alkenyl2)-(alkenyl2)-R'
 *    - ethenyl1 (primary double bond)
 *            R-(alkenyl2)-(alkenyl1)
 *    - ethene (alkenyl1)-(alkenyl1)
 *    - cctb (triple-bond) R-(alkynyl2)-(alkynyl2)-R'
 *    - ethynyl R-(alkynyl2)-(alkynyl1)
 *    - ethyne (alkynyl1)-(alkynyl1)
 *    - ketimine1 (primary ketimine) R-((alkenyl3)-R')-(nitrodb1)
 *    - ketimine2 (secondary ketimine) R-(alkenyl3)-R')-(nitrodb2)-R''
 *    - aldimine1 (primary aldimine) R-(alkenyl2)-(nitrodb1)
 *    - aldimine2 (secondary aldimine) R-(alkenyl2)-(nitrodb2)-R'
 *    - methanimine (alkenyl1)-(nitrodb1)
 *    - aromaticring
 *
 */
enum FxnGroup_t{NO_GROUP,//The "NULL" Group, useful as a default
                //The atoms lying by themselves
                HYDROGEN,CARBON,NITROGEN,OXYGEN,
                FLUORINE,CHLORINE,BROMINE,IODINE,
                //Primitive Groups by atom:
                METHANE,
                METHYL,
                METHENE,ALKENYL1,
                METHYNE,ALKENYL2,ALKYNYL1,
                CARBON4,ALKENYL3,ALKYNYL2,//C
                WATER,
                HYDROXYL,
                OXYGEN2,OXYDB,//O
                AMMONIA,
                AMINE1,
                AMINE2,NITRODB1,
                AMINE3,NITRODB2,NITROTB,//N
                HYDROGENFLUORIDE,FLUORINE1,//F
                HYDROGENCHLORIDE,CHLORINE1,//Cl
                HYDROGENBROMIDE,BROMINE1,//Br
                HYDROGENIODIDE,IODINE1,//I
                //Derived Groups in a relatively random order:
                HYDROGENCYANIDE,NITRILE,//C w/ amine3
                FORMALDEHYDE,ALDEHYDE,CARBONYL, //C w/ oxygen2
                CARBOXYL,PEROXIDE,HYDROPEROXY,METHOXY,
                METHANOL,//Other things w/oxygen2
                CCDB4,CCDB3,CCDB2,ETHENYL1,ETHENYL2,ETHENE,
                CCTB,ETHYNYL,ETHYNE,
                KETIMINE1,KETIMINE2,ALDIMINE1,ALDIMINE2,METHANIMINE,
                AROMATICRING
                /*Still need added
                 *
                 * Haloformyl R(C=O)X
                 * Carbonate ester ROC=OOR'
                 * Carboxylate R(C-O)O
                 * Ester R(C=O)OR'
                 * Peroxy ROOR`
                 * ether ROR`
                 *
                 */
};


/* \brief The base class for all functional groups
 *
 *  Fundamentally a molecule is nothing more than a mathematical graph. The
 *  atoms are nodes, and the bonds are edges.  When certain node and edge
 *  patterns occur we attribute chemical significance to them and call them
 *  functional groups.  This class is meant to serve as the base class
 *  for all such functional groups.
 *
 *  Here's how this class works.  For each functional group we define it
 *  in terms of the smallest functional groups that uniquely define it.
 *  For example we can think of a carboxyl as a carbonyl bonded to a
 *  hydroxyl.  In this example, we term the hydroxyl a "primitive"
 *  functional group b/c it's just a heavy atom, and its hydrogens.  The
 *  carbonyl, on the other hand, is made up of two primitives:
 *  a secondary oxygen (I don't know if that's really the term for an oxygen
 *  with no hydrogens attached, but tough sh*t), and a quatanary carbon
 *  (again not sure if that's the right term, but what I mean by that is
 *  defined above).  Anyways, this class is really meant to be the guts
 *  of the primitive functional groups; the derived class group below
 *  is for, you guessed it, derived groups.
 *
 * For primitive functional groups other functional groups always attach
 * to the first atom inserted into the object (i.e. Members_[0]).  Derived
 * groups are free to override this as they please.
 *
 *
 */
class FxnalGroup:public GeometricQuantity{
   private:
      ///The type of this group
      FxnGroup_t Type_;
   protected:
      ///The number of R groups this group can support
      int Order_;
      ///An array of the atoms in this group
      std::vector<int> Members_;
      ///The approved names of each functional group (just a look-up table)
      static PsiMap<FxnGroup_t,std::string> Names_;
      ///Constructor for derived groups that want to set the members
      FxnalGroup(const FxnGroup_t Type,const int Order):
         Type_(Type),Order_(Order){}
   public:
      ///Returns the order (primary, secondary, etc.)
      int Order()const{return Order_;}
      ///Returns the identity of atom i
      int operator[](const int i)const{return Members_[i];}
      ///Returns the identity of the atom R groups attach to
      virtual int AttachPoint(const int i=0)const{return (*this)[0];}
      ///Returns the number of attachment points
      virtual int NAttachPoint()const{return 1;}
      ///Returns the number of members in this group
      int size()const{return Members_.size();}
      ///Returns the type of this group
      FxnGroup_t Type()const{return Type_;}
      ///Makes a group of type Type, with NMems Members
      FxnalGroup(const FxnGroup_t Type,const int Order,
                 const int NMems,const int* Members);
      ///No memory to free
      virtual ~FxnalGroup(){}
      ///Prints out the final group
      virtual std::string PrintOut()const{return PrintOut("");}
      ///Allows for nesting of groups
      virtual std::string PrintOut(const std::string& spaces)const;
};

/* \brief A base class for groups that can be thought of as the union of
 *        several primitive functional groups
 *
 *
 */
class DerivedFxnalGrp:public FxnalGroup{
   private:
      std::vector<boost::shared_ptr<const FxnalGroup> > Groups_;
      int NPrims_;
   protected:
      typedef std::vector<boost::shared_ptr<const FxnalGroup> > Group_t;
   public:
      int GetNPrims()const{return NPrims_;}
      ///Returns -1 if the group isn't attached to i points
      int AttachPoint(const int i=0)const{
         return (i<Groups_[i]->size()?Groups_[i]->AttachPoint():-1);
      }
      virtual FxnGroup_t GetPrimI(int i)const{return Groups_[i]->Type();}
      DerivedFxnalGrp(const FxnGroup_t Type,int NPrims,int Order):
         FxnalGroup(Type,Order),NPrims_(NPrims){}
      void SetGroups(const Group_t& other);
      std::string PrintOut(const std::string& Spaces="")const;
};

class AromaticRing:public DerivedFxnalGrp{
   private:
      std::vector<int> Attachment_;
   public:
      AromaticRing(const std::vector<int>& NAttach,const int Order,
                   const Group_t& Groups);
      int NAttachPoint()const{return Attachment_.size();}
      int AttachPoint(const int i)const{return Attachment_[i];}
      std::string PrintOut(const std::string& Spaces="")const;
};

}}//End namespaces
#include "PrimitiveGroups.hh"
#include "DerivedGroups.hh"




#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_FXNALGROUP_H_ */
