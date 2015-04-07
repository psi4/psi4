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
namespace psi{
namespace LibMolecule{

/** \brief class to automate the creation of parameter types
 *
 *   The idea of this class is to automate the way the parameter names
 *   are made.  This is best done by example.  Consider a hydroxyl.
 *   We give the group a base descriptor, in this case just "O".
 *   Then we give the order of the group (number of R's it connects
 *   to).  If we were naming either the O or the H in it that would
 *   be the next argument, and finally we specify the priority in
 *   an O-Chem sense.
 *
 *   For hydroxyl our groups look like:
 *   O1
 *   O_O1
 *   H_O1
 *
 */
class ParamT{
   private:
      ///Base[0]=symbolic base Base[1]=full base
      std::string Base_[2];
      ///Atom[0]=atomic symbol, Atom[1]=full name
      std::string Atom_[2];
      ///The number of R-groups attached to this
      int Order_;
      ///The importance, in an O-Chem sense of this atom
      int Priority_;
      ///Saves results for quicker access
      std::string Full_[2];
      ///Function that generates the symbol or pretty printed name
      std::string Name(const size_t SymOrName)const;
   public:
      const std::string* Atom()const{return Atom_;}
      const std::string* Base()const{return Base_;}
      ///Constructor for the base atom
      ParamT(const std::string& BaseAbbrv,
             const std::string& BaseName,
             const int Order=-1,
             const int Priority=-1);
      ///Used when "upping" an existing parameter
      ParamT(const ParamT& Other,
             const std::string& BaseAbbrv,
             const std::string& BaseName,
             const int Order,
             const int Priority);
      /** \brief Returns true if the two parameters are equal
       *
       *   A huge portion of the process of locating functional groups is
       *   checking for equality of two nodes, which in turn involves
       *   comparing parameters.  In the most rigorous sense when we
       *   want to compare parameters we want to know if the abbreviation,
       *   order and priority are the same.  However, sometimes we don't
       *   particularly care if its say a secondary alkenyl or a tertiary
       *   one. To allow customization to the level of the comparison we
       *   propose the following: abb
       */
      bool operator==(const ParamT& other)const{
         return (Atom_[0]==other.Atom_[0]&&Base_[0]==other.Base_[0]);
      }
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
