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
#ifndef PSIAPI_PSIAPIBASE_H_
#define PSIAPI_PSIAPIBASE_H_
#include<string>

/** \brief The namespace for the PsiAPI project*/
namespace PsiAPI{
/** \brief The base class of the PsiAPI library
 *
 *  I am stemming all objects in PsiAPI from a common base class for
 *  the following reasons:
 *
 *  - I want a uniform printing method, which in particular allows all
 *    objects to be passed to std::cout or a similar object.  I implement
 *    such a printing scheme here.
 *  - The C++ way of controlling how much memory is allocated is to
 *    override:
      \code
      T::operator new (std::size_t)
      \endcode
 *    for an object of type T.  Because all classes inherit from this class
 *    we need only overwrite it in this class (not actually done yet).
 *
 */
class PsiAPIBase{
   public:
      ///Should be overridden by your class to return a string
      virtual operator std::string() const=0;
      ///Declaration of virtual destructor to remove warning
      virtual ~PsiAPIBase(){}
};

static std::ostream& operator<<(std::ostream& os, const PsiAPIBase& Base){
   return os<<std::string(Base);
}


}



#endif /* PSIAPI_PSIAPIBASE_H_ */
