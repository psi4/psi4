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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_GRAPH_H_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_GRAPH_H_

#include<list>
#include "Node.h"
namespace psi{
namespace LibMolecule{

class GraphItr{
   private:
      typedef std::list<boost::shared_ptr<Node> > List_t;
      List_t& List_;
      List_t::iterator It_;
      List_t::iterator ItEnd_;
      Node::iterator NodeIt_;
      Node::iterator NodeItEnd_;
   public:
      GraphItr(std::list<boost::shared_ptr<Node> >& List):List_(List){}
      GraphItr& begin();
      GraphItr& end();
      boost::shared_ptr<Node> operator*(){return (*NodeIt_);}
      const GraphItr& operator++();
      bool operator==(const GraphItr& other)const{
         return (NodeIt_==other.NodeIt_&&It_==other.It_);
      }
      bool operator!=(const GraphItr& other)const{return !((*this)==other);}
};


class Graph: public std::list<boost::shared_ptr<Node> >{
   public:
      ///Prints the graph out
      std::string PrintOut(const size_t Level=1)const;
      GraphItr PrimBegin(){return GraphItr(*this).begin();}
      GraphItr PrimEnd(){return GraphItr(*this).end();}
      ///Convenience function for adding a node and removing its subnodes
      void AddNode(boost::shared_ptr<Node> NewNode);

};





}}




#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_GRAPH_H_ */
