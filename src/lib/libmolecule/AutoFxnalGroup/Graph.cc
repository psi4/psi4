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
#include<sstream>
#include "Graph.h"

namespace psi{
namespace LibMolecule{
typedef boost::shared_ptr<Node> SharedNode;

std::string Graph::PrintOut()const{
   std::stringstream Message;
   Graph::const_iterator It=this->begin(),ItEnd=this->end();
   for(;It!=ItEnd;++It)Message<<(*It)->PrintOut();
   return Message.str();
}

void Graph::AddNode(SharedNode NewNode){
   std::vector<SharedNode>::iterator
      It=NewNode->SubNodes_.begin(),ItEnd=NewNode->SubNodes_.end();
   this->push_back(NewNode);
   if(It==ItEnd)return;
   for(;It!=ItEnd;++It)this->remove(*It);
}
}}//End namespaces
