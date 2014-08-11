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

#include "LibFragBase.h"
#include "psi4-dec.h"
#include <cstdarg>

namespace LibFrag{

void LibFragBase::PrintBanner(const char symbol,const std::string& message)const{
   std::string sym(Width_,symbol);
   sym+='\n';
   this->Print(sym.c_str());
   int size=message.size();
   int nspaces=2;//Number of spaces between message, on each side
   if(size<Width_-2*(nspaces+1)){//minus two from mandating symbol prints once
      //Divide message in half, extra char to left
      int lsize=(size-(size%2))/2+(size%2);
      int rsize=size-lsize;
      //Number of times to print character
      int lchars=Width_/2-nspaces-lsize;
      int rchars=Width_/2-nspaces-rsize;
      std::string lsym(lchars,symbol);
      std::string rsym(rchars,symbol);
      std::string spaces(nspaces,' ');
      rsym+='\n';
      this->Print(lsym.c_str());
      this->Print(spaces.c_str());
      this->Print(message.c_str());
      this->Print(spaces.c_str());
      this->Print(rsym.c_str());
   }
   this->Print(sym.c_str());
}

void LibFragBase::Print(const char* format,...)const{
   va_list args;
   va_start(args,format);
   //psi::vfprintf(psi::outfile,format,args);
   va_end(args);
}

}//End namespace LibFrag Base


