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
#include <stdarg.h>
#include "psi4-dec.h"
#include "StreamBase.h"

namespace psi{
/** \brief Destructor for Stream_
 *
 * We don't want to free the cout/cerr/etc. so we use this fxn
 * as boost::shared_ptr<std::ostream>'s destructor
 */
void Destructor(std::ostream* Stream_){}


void PsiStreamBase::Clone(const PsiStreamBase& other){
   this->Buffer_.str(other.Buffer_.str());
}



const PsiStreamBase& PsiStreamBase::operator=(const PsiStreamBase& other){
   if(this!=&other)this->Clone(other);
   return *this;
}

void PsiOutStream::DumpBuffer(){
   if(this->ImSpecial()){
        (*Stream_)<<Buffer_.rdbuf();
        this->Flush();
     }
     Buffer_.str("");
     Buffer_.clear();
}

void PsiOutStream::Clone(const PsiOutStream& other){
   this->Stream_=other.Stream_;
}

PsiOutStream::PsiOutStream(const PsiOutStream& other):PsiStreamBase(other){
   this->Clone(other);
}

const PsiOutStream& PsiOutStream::operator=(const PsiOutStream& other){
   PsiStreamBase::operator=(other);
   if(this!=&other)Clone(other);
   return *this;
}


void PsiOutStream::MakeBanner(const std::string& message,const char delimiter,
      const int width){
      std::string symbols(width,delimiter);
      (*this)<<symbols<<std::endl;
      int size=message.size();
      int nspaces=2;//Number of spaces between message, on each side
      if(size<width-2*(nspaces+1)){
         //Divide message in half, extra char goes left
         int lsize=(size-(size%2))/2+(size%2);
         int rsize=size-lsize;
         //Number of times to print character
         int lchars=width/2-nspaces-lsize;
         int rchars=width/2-nspaces-rsize;
         std::string lsym(lchars,delimiter);
         std::string rsym(rchars,delimiter);
         std::string spaces(nspaces,' ');
         (*this)<<lsym<<spaces<<message<<spaces<<rsym<<std::endl;
      }
      (*this)<<symbols<<std::endl;
}


PsiOutStream::PsiOutStream(SharedStream Stream){
   if(this->ImSpecial()){
      Stream_=(Stream?Stream:SharedStream(&std::cout,Destructor));
   }
}

bool PsiStreamBase::ImSpecial()const{
   return(WorldComm->me("COMM_WORLD")==0);
}

void PsiOutStream::Printf(const char* format,...){
   char buffer[1000];
   va_list args;
   va_start (args, format);
   int left=vsnprintf(buffer,1000,format,args);
   //if(left>1000)outfile->Printf("WARNING::Your entire message was not printed");
   va_end(args);
   Write2Buffer(buffer);
}

void PsiOutStream::Flush(){
   if(this->ImSpecial()){
      Stream_->flush();
   }
}

std::ostream& PsiOutStream::Write2Buffer(StreamManips fp){
   Buffer_<<fp;
   this->DumpBuffer();
   return Buffer_;
}
}//End psi namespace

