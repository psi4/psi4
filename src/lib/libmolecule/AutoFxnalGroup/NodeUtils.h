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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_NODEUTILS_H_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_NODEUTILS_H_
#include <vector>
#include <boost/shared_ptr.hpp>
#include "Node.h"
namespace psi{
namespace LibMolecule{

/** \brief A useful utility function that when given iterators to a
 *         container of objects of type: boost::shared_ptr<Node> (henceforth
 *         known as SharedNode_t) returns an array of the edges such that
 *         he first pointer is the node that was in the container and the
 *         second is its connected node
 *
 *
 *         Important note, two SharedNode_t objects are the same if their
 *         underlying pointers point to the same object.
 *
 *         Given a list, \f$S\f$, of SharedNode_t objects, the \f$i\f$-th
 *         of which is labeled \f$N_i\f$, this function loops over nodes
 *         connected to \f$N_i\f$ and returns the set of ordered pairs:
 *         \f$ \lbrace \left( N_i,N_i^j\right)\ : N_i^j \notin S\rbrace\f$
 *
 *         Guarantees: The relative ordering of the returned pairs is
 *         consistent with that of the input container.  If the returned
 *         result is empty, then the input set of nodes is a closed graph.
 *
 *         A typical usage of this function is something like:
 *         \code
 *         typedef boost::shared_ptr<Node> SharedNode_t;
 *         typedef std::pair<SharedNode_t,SharedNode_t> Pair_t;
 *         std::vector<SharedNode_t> FoundNodes=FxnThatMakesListOfNodes();
 *         std::vector<Pair_t> Connections=
 *         DetermineConns(FoundNodes.begin(),FoundNodes.end());
 *         \endcode
 *
 *  \param[in] T The type of the container holding the
 *               SharedNode_t list. e.g. std::vector<SharedNode_t> or
 *               std::set<SharedNode_t>, etc.
 *  \param[in] It An iterator to the first element of a container of Nodes
 *  \param[in] ItEnd An iterator just past the last element in the container
 *                   corresponding to It.
 *  \return An std::vector of ordered pairs of nodes, such that the first
 *          element of the pair is in the input container and the second
 *          element is connected to it, but not in the container.
 */
template<typename T>
std::vector<std::pair<boost::shared_ptr<Node>,boost::shared_ptr<Node> > >
DetermineConns(T It,T ItEnd);

/** \brief Determines which node is of higher priority
 *
 *  Given two nodes, \f$Ni\f$ and \f$Nj\f$ this function determines which
 *  has a higher priority.
 *
 *  \param[in] Ni One of the nodes to compare
 *  \param[in] Nj The other node to compare
 *  \return null if the nodes are of the same priority, otherwise
 *          the higher priority node
 */
const Node* Priority(const Node* Ni,const Node* Nj);

/********************Implementations********************/
template<typename T>
std::vector<std::pair<boost::shared_ptr<Node>,boost::shared_ptr<Node> > >
DetermineConns(T It,T ItEnd){
   //Our node types
   typedef boost::shared_ptr<Node> SharedNode_t;
   //The type of a pair of Nodes
   typedef std::pair<SharedNode_t,SharedNode_t> Pair_t;
   std::vector<Pair_t> Result;
   T ItStart=It,Itj;
   for(;It!=ItEnd;++It){
      Node::const_ConnItr Iti=(*It)->ConnBegin(),ItiEnd=(*It)->ConnEnd();
      for(;Iti!=ItiEnd;++Iti){
         for(Itj=ItStart;Itj!=ItEnd;++Itj)
            if(Itj->get()==Iti->second.get())break;
         if(Itj!=ItEnd)continue;
         Result.push_back(Pair_t(*It,Iti->second));
      }
  }
  return Result;
}






}}



#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_NODEUTILS_H_ */
