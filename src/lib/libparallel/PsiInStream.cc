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

#include "PsiInStream.h"
#include "psi4-dec.h"

namespace psi{
void Destructor(std::istream* Stream_){}
PsiInStream::PsiInStream(SharedInStream Stream){
   if(this->ImSpecial()){
      Stream_=(Stream?Stream:SharedInStream(&std::cin,Destructor));
   }
}

void PsiInStream::ReadFromStream(){
   if(this->ImSpecial()){
      std::string line;
      while(std::getline(*Stream_,line))
         Buffer_<<line;
   }
   int BufferLength=Buffer_.str().size();
   WorldComm->bcast(&BufferLength,1,WhoIsSpecial());
   char* tempbuffer=new char[BufferLength];
   WorldComm->bcast(tempbuffer,BufferLength,WhoIsSpecial());
   Buffer_<<tempbuffer;
   delete [] tempbuffer;
}


}

