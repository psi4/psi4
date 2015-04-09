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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RINGFINDER_H_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RINGFINDER_H_
#include "LinearSearch.h"

namespace psi{
namespace LibMolecule{

/** \brief The class that finds rings
 *
 *   This class is basically a linear search except that it ensures that
 *   the groups on the ends of the ring are attached.
 *
 */
template<typename...Groups>
class RingFinder: public LinearSearch<Groups...>{
   private:
      typedef LinearSearch<Groups...> Base_t;
   public:
      template<typename...Types>
      RingFinder(Types...MMTypes):Base_t(MMTypes...){}

      bool FindMe(boost::shared_ptr<Node> Nodes){
         typedef boost::shared_ptr<Node> SharedNode;
         //Can't possibly make a ring
         if(Nodes->NEdges()<2)return false;
         if(!Base_t::FindMe(Nodes))return false;
         //Check if NodeI is connected to NodeJ
         SharedNode NodeI=this->SubNodes_[0],NodeJ=this->SubNodes_.back();
         std::vector<SharedNode>::iterator
           It=NodeI->ConnNodes_.begin(),ItEnd=NodeI->ConnNodes_.end();
         for(;It!=ItEnd;++It)if(It->get()==NodeJ.get())return true;
         return false;
      }
};

}}//end namespaces




#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RINGFINDER_H_ */
