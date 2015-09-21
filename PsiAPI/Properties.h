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
#ifndef PSIAPI_PROPERTIES_H_
#define PSIAPI_PROPERTIES_H_
#include <string>
#include <boost/move/unique_ptr.hpp>
#include "Utils/PsiMap.h"
namespace PsiAPI{

/** \brief The class to hold the results of a derivative calculation
 *
 *   Operationally this class behaves nearly identical to the Reference
 *   class except that it allows you to name your returned properties.
 *   This design decision stems from the fact that there are likely to
 *   be a plethora of properties a method may want to return versus
 *   a small set that are essential for restart/initialization. It also
 *   is more slimed down (it doesn't contain a Molecule, BasisSet, or
 *   Options object).  The logic behind the latter decision is that
 *   the Properties object is associated with the Reference object and the
 *   two are assumed synched.  Another difference is we do not support
 *   multiple copies of the same tag (i.e. this is not a multimap).
 *
 *   As with the Reference object, the Properties object owns its memory,
 *   that is you should allocate the memory for the each array inside the
 *   Properties object, but not delete it.
 *
 *   For consistency we recommend the following naming conventions:
 *   -All energies are <method-name>_ENERGY
 *   -All gradients are <method-name>_GRAD
 *   -All Hessians are <method-name>_HESS
 *
 *   Typical usage:
 *   \code
 *   //Create a Properties instance
 *   Properties MyProps;
 *   //Allocate memory for SCF Energy
 *   MyProps.AddProp("SCF_ENERGY",new double[1]);
 *   //Get a stored SCF gradient (dimension is 3*NAtoms, which comes from
 *   // Reference class)
 *   const double* Grad=MyRef.GetProp("SCF_GRADIENT");
 *   \endcode
 *
 */
class Properties: private PsiMap<std::string,boost::unit::unique_ptr<double*> >{
   private:
      typedef boost::unit::unique_ptr<double *> Data_t;
      typedef PsiMap<std::string,Data_t> Base_t;
   public:
      ///Returns the property with the specified name (must already exist)
      const double* GetProp(const std::string& Str)const{
         return (*this)[Str].get();
      }
      ///Returns the property in a modifiable state (must already exist)
      double*& GetProp(const std::string& Str)const{
         return (*this)[Str].get();
      }
      ///Returns true if this instance contains a given property
      bool HasProp(const std::string& Str)const{return count(Str)==1;}
      ///Allocates memory for a property
      void AddProp(const std::string& Str,const double*& Address){
         Base_t::push_back(Data_t(Address));
      }
};

}



#endif /* PSIAPI_PROPERTIES_H_ */
