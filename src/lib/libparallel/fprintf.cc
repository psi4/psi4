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

#include <stdio.h>
#include <stdarg.h>
#include "psi4-dec.h"

namespace psi {
//extern "C"{
/*void p_fprintf(FILE * __restrict __stream, const char * __restrict __format,
      ...) {
   va_list args;
   va_start(args, __format);
   int status=0;

   status=::vfprintf(__stream, __format, args);

   va_end(args);
}

int vfprintf(FILE * __restrict __stream, const char * __restrict __format,va_list& args){
   int status=0;
   if ((WorldComm.get()!=NULL)&&(WorldComm->me("COMM_WORLD")==0)) {
         status= ::vfprintf(__stream, __format, args);
         fflush(__stream);
      }
      else if (WorldComm.get()==NULL) {
         p_fprintf(stderr, "Communicator object does not exist.\n");
         status= ::vfprintf(__stream, __format, args);
         fflush(__stream);
      }
   return status;
}

int fprintf(FILE * __restrict __stream, const char * __restrict __format, ...) {
   va_list args;
   va_start(args, __format);
   int status=psi::vfprintf(__stream,__format,args);
   va_end(args);
   return status;
}
//}*/
}
