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
#ifndef SRC_LIB_LIBPSIUTIL_IMPLEMENTATIONS_UNITSGUTS_HH_
#define SRC_LIB_LIBPSIUTIL_IMPLEMENTATIONS_UNITSGUTS_HH_
#include "../PsiMap.h"
#include <boost/shared_ptr.hpp>
namespace psi{
///The base class for the converters
template<typename T>
class Converter{
   protected:
      static boost::shared_ptr<
                PsiMap<T,PsiMap<T,double> >
             > Conversions_;
   public:
      double operator()(const T& From,const T& To)const{
         return (*Conversions_)[From][To];
      }
};

template<typename T>
boost::shared_ptr<PsiMap<T,PsiMap<T,double> > > Converter<T>::Conversions_;


}//End namespace psi




#endif /* SRC_LIB_LIBPSIUTIL_IMPLEMENTATIONS_UNITSGUTS_HH_ */
