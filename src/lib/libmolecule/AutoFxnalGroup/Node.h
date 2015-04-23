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
 *  number of composite nodes (a node that has taken the place of
 *  other nodes).  In general, with each new composite node we make,
 *  the search space becomes smaller.  Ultimately, the entire algorithm
 *  hinges on the connectivity, which is assumed to have been given to us.
 *
 *  Our algorithm is designed to automatically determine as much as it
 *  can.  Inside the header files: TrivialNodes.hh, PrimitiveNodes.hh,
 *  DerivedNodes.hh,etc. we take care of all of the definitions of what
 *  each chemically meaningful node is.  If you want a new node to be
 *  found those are the places to add it.
 *
 *  Ultimately our class structure looks something like:
 *  \verbatim
 *          Node
 *            |
 *            |
 *            |
 *           \ /
 *       LinearSearch
 *         |        |
 *         |        |
 *         |        |
 *        \ /      \ /
 * RadialSearch    RingSearch
 * \endverbatim
 *
 * off of each of these classes we derive our search patterns.  First up
 * are the groups of nodes that can be found by a linear search.  These
 * are node patterns that look like:
 * \verbatim
 * A---B---C---D
 * \endverbatim
 * i.e. the nodes lie along a continuous path that doesn't double-back
 * or loop on itself.
 *
 * Radial searches are stuff like:
 * \verbatim
 *      A
 *      |
 *   D--B--C
 *      |
 *      E
 * \endverbatim
 * where we are explicitly looking for a B node surrounded by exactly four
 * nodes, and those nodes are of types A, C, E, and D, the relative
 * ordering of which is irrelevant.  This is a linear search where instead
 * of checking if the node is connected to the last found node we always
 * check if it's connected to the central node (i.e. the above is four,
 * two-node linear searches (it may be worth considering in future
 * developments implementing such a search as two linear searches: A-B-E,
 * then D-(A-B-E)-C)).
 *
 * Finally we have ring search which is a linear search that additionally
 * requires that the last node found shares an edge with the first node
 * found.  There exists a degeneracy between certain search patterns, i.e.
 * a radial search for two nodes can be done with a linear search for three
 * nodes, but this should not prove to be any complication.
 *
 * After finding a node many of the details are directly calculable from
 * it's connectivity and identity, e.g. we know it's order, can deduce
 * the priorities of the atoms, etc.  We use this as the basis for
 * automatically defining our parameter set (more details of which
 * can be found in ParameterType.h).  This attempt at auto defining
 * parameters brings us to an important point: the distinction between
 * a search pattern and an instance of the pattern.  Consider
 * the search pattern for a methene unit:
 * \code
 * class methene: public
 * RadialSearch<4,Carbon,Hydrogen,Hydrogen>{};
 * \endcode
 *
 * To find a methene we look for an instance of a carbon,
 * that makes four bonds, two and only two of which are to hydrogen instances.
 * This is the search pattern.  The class methene
 * defines the search pattern.  When we instantiate a methene class from
 * three particular atoms, we have an instance.  This instance has more
 * information available to it then the search pattern because it's actually
 * in the molecule.  For example our instance knows where the non-hydrogen
 * edges go, and can establish a direction (which way is a higher priority
 * group).  Our instance also knows it's order.  Although all methenes
 * have an order of two, we won't be so blessed for more complicated
 * groups like carbon-carbon double bonds which can have an order ranging
 * from 0 (i.e. ethene) to 4 (a quadruply substituted double bond).
 *
 * Conceptually, the compiler writes out a series of logical if statements
 * that compare instances to a search patterns.  This means that
 * the last point in the previous paragraph is very important, when we
 * specify a more complicated search pattern, such as that of a carbonyl
 * \verbatim
 * class carbonyl: LinearSearch<Alkenyl3,ODB>{};
 * \endverbatim
 * The two classes on the right, (Alkenyl3 and ODB) are also search
 * patterns, but they are going to be compared to actual instances.
 * This means that the equality will fail because the equality check
 * of a ParameterType uses more information than is available to
 * just the search pattern.  For the vast majority of search patterns
 * it is possible to infer the missing information (a carbonyl is always
 * a secondary group, a methyne group is always tertiary etc.), but in
 * general it is not always possible.
 *
 * To circumvent this problem we define our search patterns as variadic
 * templates.  The parameter pack is to be of size 2n, where n is the
 * number of subnodes that contribute edges.  If n==0 than we have the
 * search pattern, for all other values of n we have an instance to be
 * compared against.  The
 * parameters in the pack are pairs of integers such that the first
 * value is the subnode number and the second is its order, for example
 * a tertiary carbon-carbon double bond would be:
 * \code
 * class CCDB<0,2,1,1>:LinearSearch<Alkenyl,Alkenyl>{};
 * \endcode
 *
 * The above example illustrates two important ideas: wildcards and
 * priorities.  The former is easier to explain and we do so first.
 * In general there are 6 types of carbon-carbon double bonds (we don't
 * distinguish between cis or trans):
 * \verbatim
 *
 * H        H      H        R
 *  \      /        \      /
 *   C====C          C====C
 *  /      \        /      \
 * H        H      H        H
 *
 * R        H      R        R
 *  \      /        \      /
 *   C====C          C====C
 *  /      \        /      \
 * R        H      H        H
 *
 * R        H      R        R
 *  \      /        \      /
 *   C====C          C====C
 *  /      \        /      \
 * R        R      R        R
 * \endverbatim
 * It's possible to code up the search pattern for each of these six
 * types, but rather redundant.  Instead we realize that any time we have
 * two Alkenyl like carbons next to each other, we have a carbon-carbon
 * double bond.  We thus create a wildcard called Alkenyl that is a stand
 * in for Alkenyl1, Alkenyl2, or Alkenyl3 (internally the search loops
 * over all the possible types, but that's irrelevant).  Now all six
 * double bonds have the same search pattern (for aesthic reasons we
 * siphon ethene off, but in principle it could be found this way as
 * well):
 * \code
 * LinearSearch<Alkenyl,Alkenyl>
 * \endcode
 * Although it wouldn't be too bad to code up all of the double bond
 * search patterns, coding up all of the substitutions of say indole,
 * would require much more work (for indole there are 256 possible
 * substitutions (in general there are \f$2^{N}\f$ possible substitutions
 * where \f$N\f$ is the number of edges leading from the fully
 * "R-substituted" functional group (some of them may be equivalent through
 * symmetry))).  Wildcards allow us to specify a set of functional groups
 * quite concisely.  It's also worth mentioning that wildcards automatically
 * substitute the correct instances.
 *
 * When a search pattern is being used as an instance, we need to be
 * able to specify, at compile time, what the order of each subnode is
 * in order to ensure we have the instance we truly want.  Above we wanted
 * the tertiary instance of a carbon-carbon double bond, which we requested
 * by:
 * \code
 * CCDB<0,2,1,1>::LinearSearch<Alkenyl,Alkenyl>{};
 * \endcode
 * Because of the use of wildcards, both:
 * \verbatim
 * Alkenyl3-Alkenyl2
 *
 * Alkenyl2-Alkenyl3
 * \endverabtim
 * match the search pattern.  This in theory means we could have two
 * CCDB types: CCDB<0,2,1,1> and CCDB<0,1,1,2>; however, this doesn't
 * actually happen.  The reason is because of priorities.  Basically,
 * similar to the same concept in organic chemistry we assign different
 * groups different priorities.  We use these priorities to always store
 * the subnodes in the same order.  For a radial search no priorities
 * are needed and we store the subnodes in the order requested by the
 * search pattern.  For linear searches, the subnodes are either stored
 * in the requested order, or the reverse requested order, whichever places
 * the first unique, highest priority group at a lower index.  Because
 * the first criteria for priority is order, that is how we knew the geminal
 * carbon was stored first.
 *
 * Rings are a bit more complicated.  Let the maximum number of edges coming
 * from any one node, and leading outside the ring, be \f$N\f$.  First
 * the ring is rotated such that one of the subnodes with \f$N\f$ edges
 * leading outside the ring is in the first spot.  The priorities are
 * then calculated using the linear criteria.  Next the substitution
 * pattern of the result returned by the linear criteria is evaluated,
 * i.e. is it 1,2-benzene, 3,4,5-benzene, etc.
 * Next, this process is repeated for each subnode with \f$N\f$ such
 * edges.  This ensures that our substitution pattern is always 1,... and
 * the reliance on the linear priority criteria ensures we always traverse
 * the ring in a direction that leads to the next substituted node fastest.
 * If one of the resulting patterns is the lowest, e.g. 1,2,4 is lower
 * than 1,2,5, we use it, otherwise we traverse the ring until one of the
 * patterns hits a group with higher priority than the other patterns. If
 * the latter doesn't occur the patterns are assumed equal and one is
 * chosen arbitrarily.
 *
 *
 * In general priority works such that for two groups, A and B, A is
 * higher priority than group B if:
 * - A has a higher order
 * - A has an atom of higher atomic number than those in B
 * - A contains more subnodes than B
 * The highest criteria in this list, that is not a tie, is used to
 * establish priority.
 *
 */
class Node{
   private:
      ///Convenient internal typedef of a shared_ptr to this class
      typedef boost::shared_ptr<Node> SharedNode_t;
      ///Convenient internal typedef of a pair of SharedNodes
      typedef std::pair<SharedNode_t,SharedNode_t> Pair_t;
      ///An array of the atom indices comprising this Node
      std::vector<int> MyMember_;
      ///An array of the primitive Nodes in this Node
      std::vector<SharedNode_t> PrimNodes_;
      ///These are the sub-nodes of the current node
      std::vector<SharedNode_t> SubNodes_;
      ///These are the nodes connected to the current node
      std::vector<Pair_t> ConnNodes_;
      ///A numeric mapping from ConnNodes_ to SubNodes_
      std::vector<size_t> ActiveSubNodes_;
      ///The priority pattern of node
      std::vector<size_t> Prior_;
      ///Whether the pattern is symmetric
      bool IsSymm_;
   protected:
      virtual bool IsPrim()const{return false;}
      ///As a node picks up names, they are deposited here
      std::vector<ParamT> MyTypes_;
      /** \brief Function that updates the types of a node
       *
       *  In order to update the information we need several things,
       *  first we need to know the full name and abbreviation of
       *  the new type.  This is obtained from:
       *
       *  \param[in] Parent The parameterization of the super node
       *
       *  Next are the orders, which are a map of the subnodes that
       *  still have edges, and the number of such edges.
       *
       *  \param[in] Order The nodes that still have edges and how many
       *                   in the supernode
       *
       *  Finally, the hard part, determining the priority of the new node.
       *  For the supernode we know its priority pattern and whether or
       *  not it is symmetric.  If the supernode is not symmetric than
       *  the subnodes may loose their symmetry were before, e.g.
       *  a 1,2 carbon-carbon double bond is symmetric
       *  before the R's are known, once the R's are resolved to two different
       *  groups this symmetry is gone.
       *
       *  The case when the symmetry breaks is the hardest.  Going back
       *  to the C-C double bond. Initially, both carbons have a priority
       *  of 0, because they are symmetric.  Once the R's are known the
       *  symmetry is lost, hence when numbering the C's we need to establish
       *  a direction.  To do this, we pass in the highest priority node
       *  attached to the subnode (i.e. the node left of this one in the
       *  array).  The side of the current subnode closest to it is the
       *  highest priority and numbering continues from there.
       *
       *  An additional complication comes from
       *
       *  If the symmetry is maintained, then we assign priority starting
       *  from 0, and following the symmetry pattern.  Once we are done
       *  with the first subnode we pass the final priority into the next
       *  subnode and continue until our symmetry reveals itself, at which
       *  point we undo the previous steps sequentially..
       *
       *  \param[in,out] Priority In, the starting priority of the current
       *                          subnode.  Out, the last priority used in the
       *                          current subnode.
       *  \param[in] Symm Whether the current subnode maintains its current
       *                  symmetry or not
       *  \param[in] Left The group I'm attached to with higher priority
       *
       */
      void AddTypes(const ParamT& Parent,
                    const PsiMap<size_t,size_t>& Order,
                    size_t& Priority,
                    bool Symm,
                    boost::shared_ptr<Node> Left);
      /** \brief Adds a subnode to the current node
       *
       *  In addition to adding the sub node, this fxn
       *  also updates the current node's connections so that they
       *  are the union of the connected nodes less the subnodes.
       */
      void AddSubNode(boost::shared_ptr<Node> NewNode);
      /** \brief Returns a vector of the priorities of the current node.
       *
       * This function is overridden by the various searches
       * (and if need be by the actual functional groups) so that it
       * returns a vector of length \f$N\f$ such that element \f$i\f$
       * is the priority of SubNode \f$i\f$.  Note 0 is the highest
       * priority and if two nodes possess the same priority they
       * should be symmetric in the actual group.
       *
       * The returned priorities should be in an ascending order and in
       * the order of the nodes of the deque.   This means this function
       * may change the order of the deque if need be.
       *
       * \param[in] FoundNodes The nodes we want the priority of
       * \param[out] Symm True if the nodes possess symmetry
       *
       * \return An array of the priorities
       *
       *
       */
      virtual std::vector<size_t> DeterminePrior(
            std::deque<boost::shared_ptr<Node> >& FoundNodes,
            bool& Symm)const;
   public:
      typedef std::vector<SharedNode_t>::iterator iterator;
      typedef std::vector<SharedNode_t>::const_iterator const_iterator;
      typedef std::vector<Pair_t>::iterator ConnItr;
      typedef std::vector<Pair_t>::const_iterator const_ConnItr;
      /** \defgroup Iterators @{*/
      iterator PrimBegin(){return PrimNodes_.begin();}
      iterator SubBegin(){return SubNodes_.begin();}
      ConnItr ConnBegin(){return ConnNodes_.begin();}
      iterator PrimEnd(){return PrimNodes_.end();}
      iterator SubEnd(){return SubNodes_.end();}
      ConnItr ConnEnd(){return ConnNodes_.end();}
      const_iterator PrimBegin()const{return PrimNodes_.begin();}
      const_iterator SubBegin()const{return SubNodes_.begin();}
      const_ConnItr ConnBegin()const{return ConnNodes_.begin();}
      const_iterator PrimEnd()const{return PrimNodes_.end();}
      const_iterator SubEnd()const{return SubNodes_.end();}
      const_ConnItr ConnEnd()const{return ConnNodes_.end();}
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
      Node(const int i,const std::string& BaseAbbrv,const std::string& BaseName):
         MyMember_(1,i),MyTypes_(1,ParamT(BaseAbbrv,BaseName)),
         Prior_(1,0),IsSymm_(true){}
      ///Takes a list of the subnodes
      void FillNode(std::deque<boost::shared_ptr<Node> >& NewNode);
      ///Takes the abbreviated and full names of this node
      Node(const std::string&BaseAbbrv,
           const std::string& BaseName);

      ///The number of atoms in this node
      unsigned size()const;
      ///The maximum atomic number of an atom in this node
      virtual unsigned Z()const;
      ///Returns the i-th atom in the node
      int operator[](const unsigned i)const;

      ///Tells nodes connected to this about this's new shared pointer
      void UpdateConns(boost::shared_ptr<Node> NewMe);
      ///Returns the number of types this node has
      size_t NTypes()const{return MyTypes_.size();}
      ///Returns the i-th type of the node (i<0 for main type)
      ParamT Type(const int i=-1)const{
         return (i>=0?MyTypes_[i]:MyTypes_.back());
      }
      ///Returns the number of edges for the node
      unsigned NEdges()const{return ConnNodes_.size();}
      /** \brief Adds NewNode to this node
       *
       *   This function performs two seemingly different functions
       *   depending on the identity of OldNode.  First when initializing
       *   the trivial set of nodes, if OldNode.get()==this AddConn adds
       *   NewNode as a connected node to the current node.  Otherwise,
       *   NewNode replaces oldnode.  This latter function is what is used
       *   the majority of the time as we update our graph.
       */
      void AddConn(boost::shared_ptr<Node> NewNode,
                   boost::shared_ptr<Node> OldNode);
      ///Prints out the node
      std::string PrintOut(const std::string spaces="",const size_t Level=1)const;
      ///STFU compiler
      virtual ~Node(){}
};



}}//End namespaces

#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_NODE_H_ */
