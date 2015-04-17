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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RADIALSEARCH_H_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RADIALSEARCH_H_

#include "LinearSearch.h"

namespace psi{
namespace LibMolecule{
/** \brief The class responsible for creating a node via a Radial search
 *
 *  \param[in] NBonds The number of bonds the central node should make
 *  \param[in] CenterType The type of the central node
 *  \param[in] RadialTypes The types of the radial nodes, in the order you
 *                          want them found
 *
 *  A radial search is defined as a search that examines all nodes connected
 *  to a particular node, e.g. for a methane molecule a radial search
 *  about the carbon would look at the four hydrogens attached to it.
 *  It is worth noting that a radial search for two atoms, a central
 *  atom and single atom bonded to it, is equivalent to a linear search
 *  for two atoms.  Either algorithm should work, but for aesthetic reasons
 *  we opt to use radial searches for two atoms when they are trivial nodes
 *  and linear searches when at least one is a primitive or derived node.
 *
 *  This class can handle a radial search for one atom, i.e. determining
 *  if a certain node makes a set number of bonds, e.g. the search
 *  for C4, which is of type RadialSearch<4,C>.
 *
 *  The magic of this class happens entirely within the FindMe set of
 *  functions.  Initially the public version is called with a shared
 *  pointer to a node that may be the central node of the group we
 *  are after.  That node is checked for a) to ensure it's the right type
 *  and b) that it forms the correct number of bonds.  If it passes these
 *  checks we enter recursion, else we return false.
 *
 *  If we have entered recursion we build a deque that will house our
 *  nodes as we find them.  The first node is always the central node
 *  and the remainder of the queue is the nodes attached to the central
 *  node as we find them.  At each stage of recursion we loop over the
 *  nodes attached to the central node and check if they are of the
 *  correct type and ensure they are not already in the queue.  If both
 *  of those steps pass recursion continues.  If at any point recursion
 *  fails we remove the last element in the queue and continue iterating
 *  (that's why we use a deque and not a queue).
 *
 *  Finally, once we hit the RadialSearch base class we check for the
 *  final group.  If we find it we need to update the graph.  Updating
 *  the graph entails several things:
 *  - We must consolidate all of the nodes in the queue into a single
 *    node.  We choose to do that in the "this" node of the object
 *    that called FindMe.
 *  - Inside this node we now have to figure out how many edges we have
 *    left.  This is done automatically for us by the AddSubNode command
 *  - Finally, the nodes connected to "this" node need to know that "this"
 *    node changed.  That's the user's job and can be done automatically
 *    by calling Graph::AddNode
 *
 * Implementation note: A radial search is just a linear search that needs
 * to first check if the center atom makes a certain number of bonds and
 * then needs to always check the connections of the central atom vs.
 * the back of the queue.
 */
template<unsigned NBonds,typename... RadialTypes>
class RadialSearch : public LinearSearch<RadialTypes...>{
   protected:
      ///Convenient typedef of the base class
      typedef LinearSearch<RadialTypes...> Base_t;
      ///The internal call that continues the recursion
      bool FindMe(std::deque<boost::shared_ptr<Node> >& FoundNodes);
      bool IsLinear()const{return false;}
      virtual std::vector<size_t> DeteriminePrior(
            std::deque<boost::shared_ptr<Node> >& FoundNodes,
            bool& Symm)const{
         Symm=true;
         return std::vector<size_t>(FoundNodes.size(),0);
      }
   public:
      ///Same arguments as Node's constructor
      RadialSearch<NBonds,RadialTypes...>
         (const std::string& BaseAbbrv,const std::string& BaseName):
         Base_t(BaseAbbrv,BaseName){}

      ///Given a starting node, tries to find the desired group
      bool FindMe(boost::shared_ptr<Node> Nodes);
};

/************* Implementations ******************/
template<unsigned NBonds,typename... RadialTypes>
bool RadialSearch<NBonds,RadialTypes...>::
   FindMe(boost::shared_ptr<Node> NodeI){
   if(NodeI->NEdges()!=NBonds)return false;
   return Base_t::FindMe(NodeI);
}

}}//End namespaces

#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RADIALSEARCH_H_ */
