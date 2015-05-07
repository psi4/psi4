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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_PRIMRUNNER_H_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_PRIMRUNNER_H_
#include "Graph.h"

namespace psi{
namespace LibMolecule{

template <typename T1, typename ...Ts>
class PrimRunner:public PrimRunner<Ts...> {
   private:
      typedef PrimRunner<Ts...> Base_t;
   public:
      static bool Run(Graph::iterator& It, Graph& Nodes) {
         boost::shared_ptr<T1> Temp(new T1);
         if (Temp->FindMe(*It)) {
            Temp->UpdateConns(Temp);
            Nodes.AddNode(Temp);
            It=Nodes.begin();
            return true;
         }
         return Base_t::Run(It, Nodes);
      }
};

template <typename T1>
class PrimRunner<T1> {
   public:
      static bool Run(Graph::iterator& It, Graph& Nodes) {
         boost::shared_ptr<T1> Temp(new T1);
         if (Temp->FindMe(*It)) {
            Temp->UpdateConns(Temp);
            Nodes.AddNode(Temp);
            It=Nodes.begin();
            return true;
         }
         return false;
      }
};

template<typename T1,typename...Rest>
class SetRunner: public SetRunner<Rest...>{
   private:
      typedef SetRunner<Rest...> Base_t;
   protected:
      static bool Run(Graph::iterator& It, Graph& Nodes){
         if(T1::Run(It,Nodes))return true;
         else return Base_t::Run(It,Nodes);
      }
   public:
      static void Run(Graph& Nodes){
         Graph::iterator It=Nodes.begin();
         while(It!=Nodes.end()){
            if(T1::Run(It,Nodes))continue;
            else if(!Base_t::Run(It,Nodes))++It;
         }
      }
};

template<typename T1>
class SetRunner<T1>{
   protected:
      static bool Run(Graph::iterator& It, Graph& Nodes){
         return T1::Run(It,Nodes);
      }
   public:
      static void Run(Graph& Nodes){
         Graph::iterator It=Nodes.begin();
         while(It!=Nodes.end()){
            if(T1::Run(It,Nodes))continue;
            else ++It;
         }
      }
};

}}//end namespaces

#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_PRIMRUNNER_H_ */
