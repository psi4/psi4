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
#include "exception.h"

namespace psi {

void OutFile::LoadFOptions() {
   FOptions_[NOFILEMODE]=std::fstream::out;
   FOptions_[END]=std::fstream::ate;
   FOptions_[APPEND]=std::fstream::app;
   FOptions_[TRUNCATE]=std::fstream::trunc;
   FOptions_[BINARY]=std::fstream::binary;
}

OutFile::OutFile(const std::string& filename, const FileMode& mode) :
      PsiOutStream(SharedStream()) {
   this->LoadFOptions();
   if (filename!="") Open(filename, mode);
}

void OutFile::Close() {
   if (Stream_) Stream_.reset();
   Buffer_.str("");
   Buffer_.clear();
}

void OutFile::Open(const std::string& filename, const FileMode& mode) {
   if (ImSpecial()) {
      this->Close();
      Stream_=SharedStream(
            new std::ofstream(filename.c_str(), FOptions_[mode]));
      if (!Stream_) {
         std::string error="Could not open file: "+filename;
         throw PSIEXCEPTION(error.c_str());
      }
   }
}

} //End namespace psi
