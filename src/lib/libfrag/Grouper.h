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
#ifndef GROUPER_H_
#define GROUPER_H_

#include "AtomSet.h"
#include "Connections.h"
#include "LibFragTypes.h"
#include <vector>
#include <boost/shared_ptr.hpp>

namespace psi{
namespace LibFrag{

class Grouper{
   private:
      /** \brief Determines if Path contains a ring
       *
       * This function figures out if CurrGroup is in Path assuming it's not
       * Path[Path.size()-2].  If it is, and say it is Path[2], then the Groups
       * array returned will have Group[Path[2]] set to the union of Group[Path[2]]
       * through Group[Path[Path.size()-1]] and Group[Path[Path[3]] through
       * Group[Path[Path.size()-1]] will be removed.
       */
      bool CheckRing(const std::vector<int>&Path, GroupType& Groups, const
            int CurrGroup)const;
      void RingRecursion(const int MaxSize,const Connections& GroupTable,
            GroupType& Groups, std::vector<int>& Path,bool& rfound)const;
   protected:
      ///Makes H and its heavy atom a group
      void CondenseHydrogens(GroupType& Groups,std::vector<bool>& Assigned,
            const Connections& CTable,const SharedMol& Mol2Frag)const;
      void FindRings(const int MaxSize, GroupType& Groups,
            const Connections& CTable)const;
   public:
      GroupType MakeGroups(const Connections& CTable,
            const SharedMol& Mol2Frag)const;

};


}}//End namespaces




#endif /* GROUPER_H_ */
