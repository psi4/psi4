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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_PARAMETERTYPE_H_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_PARAMETERTYPE_H_
#include <string>
#include "PsiMap.h"
namespace psi{
namespace LibMolecule{

/** \brief class to automate the creation of parameter types
 *
 *   The idea of this class is to automate the way the parameter names
 *   are made.  This is best done by example.  Consider a hydroxyl.
 *   We give the group a base descriptor, in this case just "O".
 *   Then we give the order of the group (number of R's it connects
 *   to, such that R=/=H).  Hence for a hydroxyl we come up with the
 *   name "O1" (water would be "O0", and an ether like oxygen (R-O-R)
 *   would be "O2".  This is the (abbreviated) type of the node.  For
 *   the purposes of MM we also need (abbreviated) types of each atom.
 *   Again we attempt a systematic scheme by prepending the atomic symbol
 *   to the type of the node.  For our hydroxyl the  atomic types look like:
 *   \verbatim
 *   O_O1
 *   H_O1
 *   \endverbatim
 *
 *   For non-radial groups things get a bit more complicated.  Consider
 *   finding a secondary carbon-carbon double bond of the pattern:
 *   \verbatim
 *   alkenyl2=alkenyl2
 *   \endverbatim
 *   there is no ambiguity at this point as the carbons and hydrogens are
 *   equivalent (our general philosophy is to never consider the identity
 *   of the R groups). Now consider the other possible secondary
 *   carbon-carbon double bond:
 *   \verbatim
 *   alkenyl3-alkenyl1
 *   \verbatim
 *   now the two carbons are no longer symmetric and it's possible we may
 *   want to be able to tell them apart.  This is the job of priority.
 *
 *
 *
 */
class ParamT{
   private:
      ///Base[0]=symbolic base Base[1]=full base
      std::string Base_[2];
      ///Atom[0]=atomic symbol, Atom[1]=full name
      std::string Atom_[2];
      ///The number of edges each subnode has outside this node
      PsiMap<size_t,size_t> Order_;
      ///The importance, in an O-Chem sense of this atom
      size_t Priority_;
      ///Saves results for quicker access
      std::string Full_[2];
      ///Function that generates the symbol or pretty printed name
      std::string Name(const size_t SymOrName)const;
      bool IsBase_;
   public:
      const PsiMap<size_t,size_t>& Order()const{return Order_;}
      size_t Priority()const{return Priority_;}
      const std::string* Atom()const{return Atom_;}
      const std::string* Base()const{return Base_;}
      ///Constructor for the base atom
      ParamT(const std::string& BaseAbbrv,
             const std::string& BaseName,
             const PsiMap<size_t,size_t>& Order=PsiMap<size_t,size_t>(),
             const size_t Priority=0);
      ///Used when "upping" an existing parameter
      ParamT(const ParamT& Other,
             const std::string& BaseAbbrv,
             const std::string& BaseName,
             const PsiMap<size_t,size_t>& Order,
             const size_t Priority);
      /** \brief Returns true if the two parameters are equal
       *
       *   Two parameters are equal if Full_[0]==other.Full_[0].
       *   In particular this means that not only do they need to
       *   have the same abbreviation they need to have the same
       *   order (if it's the type of the entire group) and if this
       *   is the type of an atomic subnode then they also need to have
       *   the same atomic symbols and priorities; however, such a
       *   comparision is never needed.
       *
       */
      bool operator==(const ParamT& other)const;
      ///Returns true if two parameters are not equal
      bool operator!=(const ParamT& other)const{
         return !((*this)==other);
      }
      ///Returns the pretty print of this parameter
      std::string PrintOut()const{return Full_[1];}
      ///Returns the MM type of this parameter
      std::string GenMMType()const{return Full_[0];}
};


}}
#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_PARAMETERTYPE_H_ */
