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
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "PsiInStream.h"
#include "psi4/psi4-dec.h"
//#include "../libparallel2/Communicator.h"
//#include "../libparallel2/ParallelEnvironment.h"

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
   /*int BufferLength=Buffer_.str().size();
   std::shared_ptr<const LibParallel::Communicator> Comm=WorldComm->GetComm();
   Comm->Bcast(&BufferLength,1,WhoIsSpecial());
   char* tempbuffer=new char[BufferLength];
   Comm->Bcast(tempbuffer,BufferLength,WhoIsSpecial());
   Buffer_<<tempbuffer;
   delete [] tempbuffer;*/
}


}
