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
#ifndef PSIFILEIMPL_H_
#define PSIFILEIMPL_H_

#include "psi4/libpsi4util/exception.h"
#include <fstream>

///Typedef of std::file modes
typedef std::fstream::openmode stdFMode;
namespace psi{

///For each value the first name is preferred, but let's face it typing stinks
enum FileMode{NOFILEMODE=0,END=1,ATE=1,APPEND=2,APP=2,A=2,TRUNCATE=3,TRUNC=3,
              T=3,BINARY=4,BIN=4};

///The machinery common to input and output files
template<typename T>
class PsiFileImpl{
   private:
      ///Maps enumerated FileModes to those in c++ std library
      std::map<FileMode,stdFMode> FOptions_;

      ///Function that loads FOptions_ up
      void LoadFOptions();

      ///No copying of file streams
      PsiFileImpl<T>(const PsiFileImpl<T>& /*other*/){}

      ///No assignment of file streams
      const PsiFileImpl<T>& operator=(const PsiFileImpl<T>& /*other*/){
         return *this;
      }

   protected:
      template<typename T2>
      void Open(const std::string& filename, const FileMode& Mode,
            std::shared_ptr<T2>& FileStream,const bool ImSpecial){
         if (ImSpecial&&filename!="NULL") {
               this->Close(FileStream);
               FileStream=std::shared_ptr<T>(
                     (Mode==NOFILEMODE? new T(filename.c_str()):
                     new T(filename.c_str(), FOptions_[Mode])));
               if (!FileStream) {
                  std::string error="Could not open file: "+filename;
                  throw PSIEXCEPTION(error.c_str());
               }
            }
      }
      ///We cheat to avoid the upcast
      template<typename T2>
      void Close(std::shared_ptr<T2>& FileStream){
         if(FileStream)FileStream.reset();
      }
   public:
      PsiFileImpl<T>(){
         LoadFOptions();
      }

};

template<typename T>
void PsiFileImpl<T>::LoadFOptions() {
   FOptions_[NOFILEMODE]=std::fstream::out;
   FOptions_[END]=std::fstream::ate;
   FOptions_[APPEND]=std::fstream::app;
   FOptions_[TRUNCATE]=std::fstream::trunc;
   FOptions_[BINARY]=std::fstream::binary;
}
}



#endif /* PSIFILEIMPL_H_ */
