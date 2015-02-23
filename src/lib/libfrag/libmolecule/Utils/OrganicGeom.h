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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ORGANICGEOM_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ORGANICGEOM_H_
#include "Geometry.h"
#include "PsiMap.h"
namespace psi{
namespace LibMolecule{
class FxnalGroup;
class ConnGroups;
class GroupIt{
   private:
      typedef PsiMap<int,boost::shared_ptr<const FxnalGroup> > PsiMap_t;
      PsiMap_t::const_iterator It_;
      PsiMap_t::const_iterator ItEnd_;
   public:
      boost::shared_ptr<const FxnalGroup> operator*()const{return (It_->second);}
      GroupIt& operator++();
      bool operator!=(const GroupIt& other)const{return It_!=other.It_;}
      bool operator==(const GroupIt& other)const{return !((*this)!=other);}
      GroupIt(const PsiMap<int,boost::shared_ptr<const FxnalGroup> >::
                     const_iterator& It,
              const PsiMap<int,boost::shared_ptr<const FxnalGroup> >::
                     const_iterator& End):It_(It),ItEnd_(End){}
};

class ConnGroups:public PsiMap<int,boost::shared_ptr<const FxnalGroup> >{
   private:
      typedef PsiMap<int,boost::shared_ptr<const FxnalGroup> > Base_t;
   public:
      ///Returns the group that has atom i as an attachment point
      boost::shared_ptr<const FxnalGroup> AttachedGroup(const int i)const{
         return Base_t::operator[](i);
      }
      GroupIt begin()const{return GroupIt(Base_t::begin(),Base_t::end());}
      GroupIt end()const{return GroupIt(Base_t::end(),Base_t::end());}
};

class OrganicGeom:public Geometry{
   private:
      ConnGroups FxnalGroups_;
   protected:
      virtual ConnGroups MakeFxnGroups()const;
   public:
      /** \brief Returns a reference to the Fxnal Groups this class found
       *
       *  The returned groups are a little wrapper class that also know
       *  their connectivity.  This is handy b/c you will often want to
       *  know what group is bonded to what group
       */
      const ConnGroups& GetGroups()const{return FxnalGroups_;}
      ///Nothing to be done clean-up wise
      virtual ~OrganicGeom(){}
      ///Given a presumably organic molecule, figures out functional groups
      OrganicGeom(const Molecule* Mol,const double BondDef=1.8,
            const int MaxBonds=4);
      std::string PrintOut()const;
};

}}//End namespaces




#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ORGANICGEOM_H_ */
