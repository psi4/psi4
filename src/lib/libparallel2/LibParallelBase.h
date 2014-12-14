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
#ifndef SRC_LIB_LIBPARALLEL2_LIBPARALLELBASE_H_
#define SRC_LIB_LIBPARALLEL2_LIBPARALLELBASE_H_
#include <iostream>
#include <limits>
#include "exception.h"
#include "../libPsiUtil/PsiRdmNumGen.h"
namespace psi{
namespace LibParallel{
class ParallelEnvironment;
/** \brief Basically this is my interface to however the program wants
 *         to do stuff.
 *
 *  I will want to print, throw errors, and generate random numbers.
 *  Your implementation needs to fill in the details in this class.
 */
class LibParallelBase{
   protected:
      boost::shared_ptr<ParallelEnvironment> Env_;
   public:
      LibParallelBase();
      ///How I'm printing
      template<typename T>
      const LibParallelBase& operator<<(const T& in)const{
         std::cout<<in;
         return *this;
      }
      ///How I'm throwing errors
      void Error(const std::string& Error)const{
         throw PSIEXCEPTION(Error);
      }
      ///How I'm generating RdmNumbers
      int RdmNum(const int Max=std::numeric_limits<int>::max(),const int Min=0)const{
         PsiRdmNumGen<> Gen(Max,Min);
         return Gen();
      }
};

}}//End namespaces


#endif /* SRC_LIB_LIBPARALLEL2_LIBPARALLELBASE_H_ */
