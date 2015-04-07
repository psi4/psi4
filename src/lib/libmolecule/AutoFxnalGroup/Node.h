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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_NODE_H_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_NODE_H_
#include<queue>
#include<vector>
#include<set>
#include<boost/shared_ptr.hpp>
#include "ParameterType.h"
namespace psi{
namespace LibMolecule{

/** \brief The class that serves as the basis for determining a molecule's
 *         identity
 *
 *  Every molecule can be thought of in a graph-theory sense as a series of
 *  nodes and edges.  The atoms in
 *  the molecule comprise nodes, and the bonds (regardless of bond-order)
 *  are the edges.  Chemistry is about assigning significance to certain node
 *  and edge patterns, the pattern being termed a functional
 *  group.  For example, an oxygen double-bonded to a carbon is a carbonyl.
 *
 *  Once we have found a functional group, we replace the nodes that comprise
 *  it with a single node.  For concreteness assume our graph looks like:
 *  \verbatim
 *     \   /            \ /
 *      B-C-     ---->   A-
 *     /   \            / \
 *  \endverabtim
 *  In this ASCII drawing we have combined two connected nodes, labeled B and
 *  C to form a new node A.  Note that the edges leading to B or C, that
 *  did not become part of A, are now edges that lead to A.  Going back to
 *  our carbonyl example:
 *  \verbatim
 *      O
 *     | |
 *      C     -->        (CO)
 *     / \              /    \
 *    R   R`           R      R'
 *  \endverbatim
 *  where the quantity (CO) is interpreted as the new node.
 *
 *  In general as we progress through this algorithm we will amass a large
 *  number of composite nodes (a single node that has taken the place of
 *  many nodes).  With each new composite node we make, the search space
 *  becomes smaller.  The number of edges our new composite node possesses is
 *  equal to the number of edges of each node in the new node, less the number
 *  they have in common.
 *
 *  From a chemistry perspective, each node may only possess a certain number
 *  of edges (equal to its valency).
 *
 *  Ultimately, the entire algorithm hinges on the connectivity and determining
 *  if two or more nodes are connected.  Within a node we thus keep pointers
 *  to the nodes connected to it.  In this way we can update the local
 *  environment of our graph from any node.
 *
 *  The algorithm starts with each atom being assigned to a trivial node
 *  type that is just it's atomic symbol and ensures each atom is on the
 *  graph.  We then go back over the atoms and create primitive nodes, which
 *  are the heavy atom and its hydrogens.  It is at this stage we account
 *  for bond order, e.g. a carbon connected to three nodes, has double-bond
 *  character, whereas a carbon connected to two nodes has triple-bond
 *  character.  I specify "character" because it is conceivable that the
 *  carbon does not actually make a double or triple bond, but is rather a
 *  radical.  The presence of what I will call an alkenyl carbon (carbon w/ 3
 *  bonds), but no double bond indicates that either your structure is f'ed
 *  up or it contains a radical/anionic carbon (and I would still call
 *  it f'ed up).
 *
 *  As a systematic way of naming things, a node that makes \f$n\f$ bonds
 *  maximally, and is bonded to \f$n-1\$ hydrogens is said to be a primary
 *  node.  If it is bonded to \f$n-2\f$ then it is a secondary node, etc.
 *  For example:
 *  \verbatim
 *  H     R
 *   \   /
 *    C=C      ---->  (Alkenyl1)-(Alkenyl2)-R ----> (CCDB1)-R
 *   /   \
 *  H     H
 *  \endverbatim
 *  , this functional group is decomposed into a primary and secondary
 *  alkenyl carbon, which are then joined to form a primary carbon-carbon
 *  double bond.
 *
 *  Every time a series of nodes are combined to form a new node each atom
 *  must be assigned a new name.  In the above example we have:
 *  \verbatim
 *  H--->HAlkenyl1--->HCCDB1
 *  C--->CAlkenyl1--->CCCDB1  C--->CAlkenyl2--->CCCDB1
 *  H--->HAlkenyl1--->HCCDB1  H--->HAlkenyl2--->HCCDB1
 *  \endverbatim
 *  where the placement of the atom in the diagram is meant to correspond
 *  to which atom it is.
 *
 *  In order to facilitate naming atoms, the atom order must be fixed for
 *  each fxnal group.  This is trivial for the lone atoms and is enforced
 *  for all other groups by making the atoms appear in the search order
 *  requested.  For example, when we look for the Alkenyl1 group above we
 *  will look for Carbon,Hydrogen,Hydrogen, and when we look for Alkenyl2 we
 *  will look for Carbon,Hydrogen.  Thus for Alkenyl1 we would specify types
 *  of CAlkenyl1, HAlkenyl1, HAlkenyl1 and for Alkenyl2 we would specify
 *  types of CAlkenyl2, HAlkenyl2.  When we look for the double bond we
 *  will search for Alkenyl2, Alkenyl1, and our types would thus then be:
 *  CCCDB1,HCCDB1,CCCDB1,HCCDB1,HCCDB1.  In general the search string should
 *  always start with a node that will contribute an edge to resulting
 *  node, this facilitates the next part.
 *
 *  Some nodes will have multiple edges that originate from different
 *  subnodes, e.g. double bonds, or disubstitued benzenes.  This leads to
 *  two problems: an ambiguity in the order of the subnodes and an
 *  ambiguity in how they connect.  Because all nodes that have edges start
 *  with one of the subnodes that have edges we can enforce some order
 *  by reordering the subnodes so that the first one is the one that
 *  attaches to our group.  For something like a double bond this resolves
 *  the ambiguity because if we ask for
 *
 *  All of the above examples suggest that there are two ways to search
 *  for nodes
 *  in order to make a new node: 1) Grab nodes around the current node,
 *  we term this a radial search or 2) Grab nodes along a consecutive path,
 *  which we term a linear search.
 *
 *  Operationally, the node class has two constructors.  They both take a
 *  series of types, the first of which is always the type of the current
 *  node, the remainder are the types of the atoms in the node in the
 *  predefined order (see above).  One of the constructors also takes
 *  an integer argument.  This constructor is the one that should be
 *  used initially when atoms are being turned into nodes and the integer
 *  corresponds to the index of that atom in the molecule.
 */
class Node{
   private:
      ///An array of the atom indices comprising this Node
      std::vector<int> MyMember_;
      ///An array of the primitive Nodes in this Node
      std::vector<boost::shared_ptr<Node> > PrimNodes_;
      ///These are the sub-nodes of the current node
      std::vector<boost::shared_ptr<Node> > SubNodes_;
      ///These are the nodes connected to the current node
      std::vector<boost::shared_ptr<Node> > ConnNodes_;
   protected:
      ///As a node picks up names, they are deposited here
      std::vector<ParamT> MyTypes_;
      ///Assigns the types to a subnode
      void AddTypes(const ParamT& Parent,
                    const size_t Order,
                    const size_t Priority);
      /** \brief Adds a subnode to the current node
       *
       *  In addition to adding the sub node, this fxn
       *  also updates the current node's connections so that they
       *  are the union of the connected nodes less the subnodes.
       */
      void AddSubNode(boost::shared_ptr<Node> NewNode,
                      const size_t Order,
                      const size_t Prior);
      ///This is the set of active subnodes (ones still with edges)
      std::set<size_t> ActiveSubNodes_;
   public:
      typedef std::vector<boost::shared_ptr<Node> >::iterator iterator;
      typedef std::vector<boost::shared_ptr<Node> >::const_iterator
            const_iterator;
      /** \defgroup Iterators @{*/
      iterator PrimBegin(){return PrimNodes_.begin();}
      iterator SubBegin(){return SubNodes_.begin();}
      iterator ConnBegin(){return ConnNodes_.begin();}
      iterator PrimEnd(){return PrimNodes_.end();}
      iterator SubEnd(){return SubNodes_.end();}
      iterator ConnEnd(){return ConnNodes_.end();}
      const_iterator PrimBegin()const{return PrimNodes_.begin();}
      const_iterator SubBegin()const{return SubNodes_.begin();}
      const_iterator ConnBegin()const{return ConnNodes_.begin();}
      const_iterator PrimEnd()const{return PrimNodes_.end();}
      const_iterator SubEnd()const{return SubNodes_.end();}
      const_iterator ConnEnd()const{return ConnNodes_.end();}
      /** @}*/


      /** \brief Returns the subnode of this connected to other
       *
       *   This function is intended for use when trying to figure out
       *   if a connection has a certain quality.  Like is this connection
       *   between Alkenyl3 carbon of an CCDB3 and the Alkenyl2 carbon
       *   of another CCDB3?  It simply returns the first subnode of this
       *   node that is connected to the Other node.  In particular this
       *   will only work if both subnodes are still be updated, i.e. they
       *   are active edges of the nodes they are part of.
       */
      boost::shared_ptr<Node> GetConnSubNode(boost::shared_ptr<Node> Other);

      ///Prim Node constructor takes index of the atom, its Sym and Name
      Node(int i,const std::string& BaseAbbrv,const std::string& BaseName):
         MyMember_(i),MyTypes_(1,ParamT(BaseAbbrv,BaseName)){}
      ///The sym of this group, and its Name
      Node(const std::string&BaseAbbrv,const std::string& BaseName):
         MyTypes_(1,ParamT(BaseAbbrv,BaseName)){}

      ///The number of atoms in this node
      unsigned size()const;
      ///The maximum atomic number of an atom in this node
      virtual unsigned Z()const;
      ///Returns the i-th atom in the node
      int operator[](const unsigned i)const;

      ///Given a boost shared pointer to the current object modifies conns
      void UpdateConns(boost::shared_ptr<Node> NewMe);
      ///Returns the number of types this node has
      size_t NTypes()const{return MyTypes_.size();}
      ///Returns the i-th type of the node (i<0 for main type)
      ParamT Type(const int i=-1)const{
         return (i>=0?MyTypes_[i]:MyTypes_.back());
      }
      ///Returns the number of edges for the node
      unsigned NEdges()const{return ConnNodes_.size();}
      ///Adds NewNode, optionally in place of OldNode (if OldNode!=Null)
      void AddConn(boost::shared_ptr<Node> NewNode,
                   boost::shared_ptr<Node> OldNode=
                         boost::shared_ptr<Node>());
      ///Prints out the node
      std::string PrintOut(const std::string spaces="")const;
      ///STFU compiler
      virtual ~Node(){}
};



}}//End namespaces

#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_NODE_H_ */
