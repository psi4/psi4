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
#include "FxnalGroupTypes.hh"
namespace psi{
namespace LibMolecule{

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
//class FxnalGroup:public GeometricQuantity{
class FxnalGroup{
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
      virtual int AttachPoint(const int i=0)const{return (*this)[i];}
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
      ///Allows for nesting of groups
      virtual std::string PrintOut(const std::string& spaces="")const;
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


template<typename MyType, typename T, typename...Args>
class FindDerGroup;

template<int Valency, typename GroupName,typename...Args>
class FindPrimGroup;

/** \brief Helper class to reduce copy/paste code
 *
 *  \param[in] MyType the type of the resulting group
 *  \param[in] i The number of groups comprising the resulting group
 *
 *  The rest of the parameters are the groups that comprise the new
 *  group, in the order they should be looked for
 *
 */
template<FxnGroup_t MyName, typename MyType,int order,typename...Args>
class GenDerGrp:public DerivedFxnalGrp{
private:
	typedef DerivedFxnalGrp Base_t;
public:
      GenDerGrp<MyName,MyType,order,Args...>():
    		  Base_t(MyName,(sizeof...(Args)),order){}
      typedef FindDerGroup<MyType,Args...> FindMe;
};

template<int NBonds,FxnGroup_t MyName, typename MyType,
          int order,typename...Args>
class GenPrimGrp:public GenDerGrp<MyName,MyType,order,Args...>{
private:
	typedef GenDerGrp<MyName,MyType,order,Args...> Base_t;
public:
	typedef FindPrimGroup<NBonds,MyType,Args...> FindMe;
};

}}//End namespaces
#include "AtomicTypes.hh"
#include "PrimitiveGroups.hh"
#include "DerivedGroups.hh"




#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_FXNALGROUP_H_ */
