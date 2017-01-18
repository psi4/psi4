/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#include "StreamBase.h"

namespace psi{

/** \brief Specialization of PsiStream to input streams
 *  (I doubt this class will ever be used by itself, but I wanted to mirror
 *  the OutStream hierarchy)
 *
 *  Theoretically, this class would allow us to read from cin; however, its
 *  real purpose is to serve the same role as PsiOutStream, except for
 *  InputFiles.
 */
class PsiInStream:public PsiStreamBase<std::istream>{
   protected:
      typedef PsiStreamBase<std::istream> BaseType;
      void ReadFromStream();
   public:
      PsiInStream(const PsiInStream& other):BaseType(other){}
      const PsiInStream& operator=(const PsiInStream& other){
         BaseType::operator=(other);
         return *this;
      }
      virtual ~PsiInStream(){}

      PsiInStream(SharedInStream Stream=SharedInStream());

      template<typename T>
      std::istream& operator>>(T& thing){
         Buffer_>>thing;
         return Buffer_;
      }

      ///Buffer is false when we have everything from it
      operator bool(){return (bool)Buffer_;}

};
}