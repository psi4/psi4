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
#ifndef PSIAPI_PSIMETHOD_H_
#define PSIAPI_PSIMETHOD_H_
#include "Reference.h"

namespace PsiAPI{
class Options;
class Reference;

/** \brief The class serving as the base class for any electronic structure
 *         method capable of calculating an energy.
 *
 *  As far as your method is concerned the relevant sequence of events unfolds
 *  like this:
 *  - InitializeOptions will be called, which your method will implement so
 *    that the Options object passed in is set to your default values.
 *  - Based on the Options you set, and the user's initial input, you will
 *    be given an initial reference via the PsiMethod::SetReferencer
 *
 *
 */

class PsiMethod{
   protected:
      Reference Ref_;
   public:
      PsiMethod(){}
      virtual ~PsiMethod(){}
      ///Called by driver to set your method's initial reference
      void SetReference(const Reference& Ref){Ref_=Ref;}
      ///Allows access to your Reference, primarily for checkpointing
      const Reference& GetRef()const{return Ref_;}
      /** \brief Function that sets your method's default options
       *
       *  This function is to be implemented by your method when it
       *  derives from this base class.  In particular it should modify
       *  the Options object so that it is consistent with your method's
       *  default options.  See the documentation of the Options class
       *  for more information how to use it.
       */
      virtual void InitializeOptions(const Options&)=0;
      /** \brief Function that returns the "Order"-th derivative of your
       *         method.
       *
       *  Let us assume that your method has a function: MyMethodEnergy
       *  which returns a dynamically allocated double of sizeMinimally this function should be modified something akin to:
       *  \code
       *  MyMethod::Deriv(size_t Order=0){
       *    double * ReturnValue=
       *     (Order==0?ReturnValue=EnergyMyMethod():PsiMethod::Deriv(Order));
       *
       *  }
       *  \endcode
       *
       */
      double* Deriv(size_t Order=0);
};


}//End namespace PSIAPI

#endif /* PSIAPI_PSIMETHOD_H_ */
