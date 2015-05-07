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
#include "AutoFxnalGroup/Graph.h"
namespace psi{
namespace LibMolecule{

class OrganicGeom:public Geometry{
   protected:
      ///The fxn that assigns the functional groups
      virtual Graph MakeFxnGroups(bool FindAAs=false)const;
      ///The actual groups we found
      Graph FxnalGroups_;
   public:
      /** \brief Returns a reference to the Fxnal Groups this class found
       *
       *  The returned groups are a little wrapper class that also know
       *  their connectivity.  This is handy b/c you will often want to
       *  know what group is bonded to what group
       */
      const Graph& GetGroups()const{return FxnalGroups_;}
      ///Nothing to be done clean-up wise
      virtual ~OrganicGeom(){}
      ///Given a presumably organic molecule, figures out functional groups
      OrganicGeom(const Molecule& Mol,bool FindAAs=false);
      ///Groups to a pretty (fyi:beauty is in the eye of the beholder) string
      std::string PrintOut()const;
};

}}//End namespaces




#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ORGANICGEOM_H_ */
