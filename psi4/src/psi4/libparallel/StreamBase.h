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
#ifndef STREAMBASE_H_
#define STREAMBASE_H_
#include<iostream>
#include<sstream>
#include<map>
#include<memory>
#include "BasesBase.h"
//This is the signature of std::endl and other sorts of iomanip things
typedef std::ostream& (*StreamManips)(std::ostream&);

///A shared output stream
typedef std::shared_ptr<std::ostream> SharedOutStream;

///A shared input stream
typedef std::shared_ptr<std::istream> SharedInStream;



namespace psi{
/** \brief The base class for Psi4's new parallel safe printing system
 *
 *  At this level we take care of the "who gets to print" sort of details.
 *  The call PsiStreamBase::ImSpecial is the function that determines this; if
 *  down the road people want to change this, this is where to do it at.  Other
 *  than that, this class takes care of copying the Buffer_ and synching it
 *  across MPI processes.
 *
 *  template parameter T Either std::istream or std::ostream for input/output streams
 *                       respectively
 */
template <typename T>
class PsiStreamBase:public BasesBase{
   private:
      typedef PsiStreamBase<T> MyType;

      ///Makes this a copy of other, copies data in Buffer_, not &Buffer
      void Clone(const MyType& other){
         this->Buffer_.str(other.Buffer_.str());
         this->Stream_=other.Stream_;
      }


   protected:
      ///This is where each MPI task ultimately writes from or to
      std::stringstream Buffer_;

      ///The actual stream
      std::shared_ptr<T> Stream_;

      ///I always forget how to empty a stringstream
      void EmptyBuffer(){
         Buffer_.str("");
         Buffer_.clear();
      }

   public:

      ///Calls Clone for actual copy
      PsiStreamBase<T>(const MyType& other){this->Clone(other);}

      ///Calls Clone for assignment iff this!=&other, returns *this
      const MyType& operator=(const MyType& other){
         if(this!=&other)this->Clone(other);
         return *this;
      }

      ///Default constructor of stringstream is called
      PsiStreamBase<T>(){}

      ///Memory not worried about until we get down the class tree to files
      virtual ~PsiStreamBase<T>(){}
};
}//End psi namespace
#endif /* STREAMBASE_H_ */
