/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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
#include<vector>
#include<execinfo.h>
#include<cxxabi.h>
#include<stdlib.h>
#include "../Exception2.h"
namespace psi{

PsiException2::PsiException2(const std::string& arg,
              const char* file,
              const int line,
              const size_t NCalls){
   std::stringstream Error;
   Error<<std::endl<<"Fatal Error: "<<arg<<std::endl
        <<"Error occurred in file: "<<file<<" on line: "<<line<<std::endl;
   std::vector<void *> Stack(NCalls);
   char **strings;
   int size=backtrace(&Stack[0],NCalls),status=-1;
   Error<<"The most recent "<<(size<NCalls?size:NCalls)
        <<" function calls were:"<<std::endl<<std::endl;
   strings=backtrace_symbols(&Stack[0],size);
   for(int i=0;i<size;i++){
      //This part from https://panthema.net/2008/0901-stacktrace-demangled/
      char *begin_name = NULL, *begin_offset = NULL, *end_offset =NULL;
      for (char *p = strings[i]; *p; ++p){
         if (*p == '(') begin_name = p;
         else if (*p == '+')begin_offset = p;
         else if (*p == ')' && begin_offset) {
            end_offset = p;
            break;
         }
      }
      if (begin_name&&begin_offset&&end_offset
            &&begin_name<begin_offset){
             *begin_name++ = '\0';
             *begin_offset++ = '\0';
             *end_offset = '\0';
             char* demangledname=
                   abi::__cxa_demangle(begin_name,0,0,&status);
             if(status==0)
                Error<<demangledname<<std::endl;
             free(demangledname);
      }
   }
   Error_=Error.str();
}

}//End namespace
