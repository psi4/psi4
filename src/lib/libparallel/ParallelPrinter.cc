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
#include "psi4-dec.h"
#include "ParallelPrinter.h"
#include <stdio.h>
#include <stdarg.h>
namespace psi {

void ParallelStream::CloseImpl() {
   bool Good2Write=(IsOpen() && Mode==WRITE);
   std::cout<<Buffer.str()<<std::endl;
   if (WorldComm->me()==0 && Good2Write)
      File<<Buffer.rdbuf();
   if(IsOpen())File.close();
   //Empty buffer
   Buffer.clear();
   Buffer.str(std::string());
}

stdFMode ParallelStream::GetFileOptions(const FileOptions& Options){
   stdFMode mode;
   switch(Options){
      case(NOOPTIONS):{break;}
      case(END):{mode=std::fstream::ate;break;}
      case(APPEND):{mode=std::fstream::app;break;}
      case(TRUNCATE):{mode=std::fstream::trunc;break;}
      case(BINARY):{mode=std::fstream::binary;break;}
   }
   return mode;
}

stdFMode ParallelStream::ReadorWrite(const FileModes& Mode){
   stdFMode mode;
   switch(Mode){
      case(NOMODE):{break;}
      case(READ):{mode=std::fstream::in;break;}
      case(WRITE):{mode=std::fstream::out;break;}
   }
   return mode;
}

void ParallelStream::OpenImpl(const std::string& name,const FileModes InMode,
      const FileOptions Options) {
   Mode=InMode;
   if(WorldComm->me()==0) {
      std::cout<<"Opening "<<name<<std::endl;
      File.open(name.c_str(),ReadorWrite(Mode)|GetFileOptions(Options));
      if (Mode==READ){
         if (WorldComm->me()==0) Buffer<<File.rdbuf();
         std::string temp=Buffer.str();
         WorldComm->bcast(temp, 0);
         Buffer.clear();
         Buffer.str(temp);
      }
   }

}

OutFile::OutFile(const std::string& name,const FileOptions& options) {
   if (name!="NONE") Open(name, options);
}

void OutFile::Open(const std::string& name,const FileOptions& options) {
   ParallelStream::OpenImpl(name,WRITE,options);
}

void OutFile::Printf(const char* format,...){
   char buffer[1000];
   va_list args;
   va_start (args, format);
   int left=vsnprintf(buffer,1000,format,args);
   if(left>1000)fprintf(outfile,"WARNING::Your entire message was not printed");
   va_end(args);
   Read(buffer);
}

InFile::InFile(const std::string& name, const FileOptions& options):
      RealEoF(false){
   if (name!="NONE") Open(name, options);
}

void InFile::Open(const std::string& name,
      const FileOptions& options) {
   ParallelStream::OpenImpl(name,READ,options);
}

} //End namespace psi
