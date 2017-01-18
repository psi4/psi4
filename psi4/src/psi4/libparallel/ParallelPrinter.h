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
#ifndef PARALLELPRINTER_H_
#define PARALLELPRINTER_H_

#include "PsiOutStream.h"
#include "PsiFileImpl.h"
#include <fstream>

namespace psi{


/** \brief The interface to the new Psi4 ASCII outfile writer
 *
 * First thing to note is that copying of this class is disallowed and is
 * enforced by making the assignment operator and copy constructor private.
 * This is because the C++ standard does not allow copying of streams, which
 * is surprising to most people when they first hear it, but in the long run
 * makes a lot of sense.  Streams are, as the name suggests, the mechanism
 * through which data flows; they are not the data itself.  Most likely when
 * someone wants to copy a stream, they want to copy the data itself (this is
 * what we do in the copy constructor of PsiStreamBase) and not the mechanism
 * by which the data flows.  In order to avoid this confusion the C++ standard
 * makes copying prohibited.  Now in reality we could allow copying because we
 * are dealing with pointers and not the actual objects; copying of a pointer
 * is of course not prohibited, but I'll stick to the example laid forth by
 * the C++ standard.
 *
 * The other thing to note is that this class keeps the file open until it
 * goes out of scope or you tell it to close.  With the exception of Psi4's main
 * output file, this means you probably want to print what you gotta' print and
 * then close it, don't keep it open any longer then you have to.
 *
 */
class OutFile: public PsiOutStream,private PsiFileImpl<std::ofstream>{
   private:
      typedef PsiFileImpl<std::ofstream> FileBase;

   public:
      /** \brief Constructor that opens a file with name "filename"
       *
       *  Upon construction passes NULL to base, so that the pointer is
       *  not set to cout.  This avoids free-ing cout in the Open call.
       *  By default no file is opened, but supplying a filename will open
       *  that file in APPEND mode.  Supplying a mode as well as a filename
       *  opens the file in that mode instead.  The mechanism of the
       *  constructor is to call Open so behavior is identical between
       *  initializing an object with a file and calling Open after calling
       *  the default constructor.
       *
       *  \param[in] filename Name of file to open (path relative to directory
       *                      Psi was called from)
       *
       *  \param[in] mode     The mode the file will be opened in
       *
       *
       */
      OutFile(const std::string& filename="",const FileMode& mode=APPEND);

      ///Calls Close to release memory
      ~OutFile(){Close();}

      ///Opens "filename" in write mode "mode"
      void Open(const std::string& filename,const FileMode& mode);

      ///Closes file and frees Stream_ pointer
      void Close();
};

}//End namespace psi
#endif /* PARALLELPRINTER_H_ */